module mod_hmm
use,intrinsic :: ISO_FORTRAN_ENV, only : REAL64,STDOUT=>OUTPUT_UNIT
use mod_rps
implicit none
private
public :: rps_hmm
public ::  assignment (=)
!
type rps_hmm
  integer                         :: D = 1
  integer                         :: K = 9
  real(REAL64)                    :: L = -HUGE(0D0)
  real(REAL64)                    :: strategy = 100D0
  real(REAL64),allocatable        :: c(:)
  real(REAL64),allocatable        :: A(:,:)
  real(REAL64),allocatable        :: traceback(:)
  integer                         :: wld(3) = 0
  type(rps_model),allocatable     :: phi(:)
contains
  procedure    :: init           => rps_hmm_init
  procedure    :: fit            => rps_hmm_fit
  procedure    :: pred           => rps_hmm_pred
  procedure    :: test           => rps_hmm_test
  procedure    :: show_status    => rps_hmm_show_status
  procedure    :: log_likelihood => rps_hmm_log_likelihood
  procedure    :: winning_ratio  => rps_hmm_winning_ratio
  procedure    :: n_episode      => rps_hmm_n_episode
  procedure    :: n_dimension    => rps_hmm_n_dimension
  procedure    :: clear          => rps_hmm_clear
  final        :: rps_hmm_destroy
end type rps_hmm
!
  interface assignment (=)
    module procedure rps_hmm_assign
  end interface assignment (=)
!
contains
!
  subroutine rps_hmm_init(this,D)
  class(rps_hmm),intent(inout) :: this
  integer,intent(in),optional  :: D
!
    call this%clear()
!
    if(present(D))   this%D   = D
!
    ALLOCATE( this%c(this%D), this%A(this%D,this%D), &
   &          this%phi(this%D),                      &
   &          this%traceback(0) )
!
    this%c(:)   = 0D0
    this%A(:,:) = 1D0 / REAL(this%D,REAL64)
!
  end subroutine rps_hmm_init
!
  subroutine rps_hmm_assign(this,rhs)
  class(rps_hmm),intent(inout) :: this
  class(rps_hmm),intent(in)    :: rhs
!
    call this%init(rhs%D)
!
    this%K      = rhs%K
    this%L      = rhs%L
    this%wld    = rhs%wld
    this%c(:)   = rhs%c(:)
    this%A(:,:) = rhs%A(:,:)
    this%phi(:) = rhs%phi(:)
!
  end subroutine rps_hmm_assign
!
  subroutine rps_hmm_fit(this,x,maxiter,threshold,n_episode,n_event,lambda)
  class(rps_hmm),intent(inout)     :: this
  integer,intent(in)               :: x(:,:)
  integer,intent(in),optional      :: maxiter,n_episode,n_event
  real(REAL64),intent(in),optional :: threshold,lambda
  real(REAL64),allocatable         :: gam(:,:)
  real(REAL64),allocatable         :: xi(:,:,:)
  real(REAL64),allocatable         :: p(:)
  real(REAL64),allocatable         :: a(:,:)
  real(REAL64),allocatable         :: b(:,:)
  real(REAL64),allocatable         :: c(:)
  real(REAL64),allocatable         :: w(:,:)
  real(REAL64),allocatable         :: pred(:)
  real(REAL64)                     :: prev,new,ths
  integer                          :: i,j,k
  integer                          :: n,d
  integer                          :: n_em,nmax
!
    d    = this%d
    n    = SIZE(x,2)
!
    if(PRESENT(n_event)) this%K = n_event
!
    ALLOCATE( gam(d,n),    &
   &          xi(d,d,n-1), &
   &          p(d),        &
   &          a(d,n),      &
   &          b(d,n),      &
   &          c(n),        &
   &          w(d,n),      &
   &          pred(this%K) &
   &         )
!
    call RANDOM_SEED()
!
    do i=1,this%D
!
      call this%phi(i)%init( N=n_episode,K=n_event,lam=lambda)
!
      call RANDOM_NUMBER( new )
      j = INT( new * ( n - this%phi(i)%n_episode() - 1 ) ) + 1
      k = j + this%phi(i)%n_episode()
!
      call this%phi(i)%maximize(x(:,j:k))
!
    enddo
!
    gam(:,:)  = 1D0 / REAL(d,REAL64)
!    xi(:,:,:) = 0D0 ; xi(:,:,1) = this%A(:,:)
    xi(:,:,:) = 0D0 ; xi(:,:,1) = 1D0 / REAL(d,REAL64)
    a(:,:)    = 1D0 / REAL(d,REAL64)
!
    nmax = 1000 ; if(PRESENT(maxiter))   nmax = maxiter
    ths  = 1D-8 ; if(PRESENT(threshold)) ths  = threshold
!
    new  = - HUGE(0D0) / 2D0
    prev = - HUGE(0D0)
!
    this%wld = 0
!
    do n_em=1,nmax
!
!     get_max_p ... PRML (13.18) Skiped
!     get_max_A ... PRML (13.19)
!
      this%A(:,:) = 0D0
      do i=1,n-1
        this%A(:,:) = this%A(:,:) + xi(:,:,i)
      enddo
      do i=1,d
        this%A(:,i) = this%A(:,i) / SUM( this%A(:,i) )
      enddo
!
      if( new - prev < ths ) EXIT
      prev  = new
!
!     get_weight
!
      w(:,:) = 0D0
      do i=1,d
        do j=1,n
          w(i,j) = this%phi(i)%pdf( x(:,j) )
        enddo
      enddo
!
!     get_alpha
!     PRML (13.37)
!
      a(:,1) = w(:,1) * gam(:,1) / SUM( gam(:,1) )
      c(1)   = SUM( a(:,1) )

!     PRML (13.36) -> (13.59)
      do j=1,n-1
        do i=1,d
          a(i,j+1) = SUM( this%A(i,:) * a(:,j) ) * w(i,j+1)
        enddo
        c(j+1)   = SUM( a(:,j+1) )
        a(:,j+1) = a(:,j+1) / c(j+1)
      enddo
!
!     get_beta
!     PRML (13.39)
!
      b(:,n) = 1D0
!
!     PRML (13.38) -> (13.62)
!
      do j=n-1,1,-1
        do i=1,d
          b(i,j) = SUM( this%A(:,i) * b(:,j+1) * w(:,j+1) )
        enddo
        b(:,j) = b(:,j) / c(j+1)
      enddo
!     PRML (13.63)
!
      this%L = SUM( LOG( c ) )
!
      new = this%log_likelihood()
      if(new/=new) RETURN        ! terminate when llh is NaN
!
!     get_max_gamma ... PRML (13.33) -> (13.64)
!
      gam(:,:) = a(:,:) * b(:,:)
!
!     get_max_xi    ... PRML (13.43) -> (13.65)
!
      do j=1,n-1
        do i=1,d
          xi(:,i,j) = a(i,j) * w(:,j+1) * b(:,i+1) * this%A(:,i) / c(j+1)
        enddo
      enddo
!
!     maximize_theta
!
      do i=1,d
        call this%phi(i)%maximize( x, gam(i,:) )
        this%c(i) = this%phi(i)%n_data()
      enddo
      this%c(:) = this%c(:) / SUM( this%c(:) )
!
      this%traceback = [this%traceback,new]
!
    enddo
!
    call linear_search(this,x)
    this%wld = this%test(x)
!
  contains
!
    subroutine linear_search(this,x)
    type(rps_hmm),intent(inout) :: this
    integer,intent(in)          :: x(:,:)
    real(REAL64)                :: s(0:2),f(2)
    real(REAL64)                :: sr,fr ! reflect
    real(REAL64)                :: se,fe ! expand
    real(REAL64)                :: si,fi ! inner contract
    real(REAL64)                :: so,fo ! outer contract
    real(REAL64)                :: conv,thre
    integer                     :: best,worst
    integer                     :: i,maxiter
!
! nelder_mead parameters
!
!    alpha = 1d0
!    beta  = 2d0
!    gam   = 5d-1
!
      s(1)   =  1D0 ; f(1) = objective(this,x,s(1))
      s(2)   = -1D0 ; f(2) = objective(this,x,s(2))
      s(0)   = s(MINLOC(f,1))
!
      conv    = 0D0
      thre    = 1D-4
      maxiter = 100
!
      do i=1,maxiter
!
        best  = MINLOC(f,1) ; worst = MAXLOC(f,1)
!
! X_c = X_best in 1D simplex
! sr is refrected points : S = s(best) + alpha * (s(best) - s(worst))
!
        sr = s(best) + s(best) - s(worst)
        fr = objective(this,x,sr)
!
        if(f(best) <= fr .and. fr < f(worst))then
!
        s(worst) = sr ; f(worst) = fr
!
        elseif(fr < f(best))then
!
!! se is expanded point : S = s(best) + beta * (sr - s(best))
!
          se     = sr + sr - s(best)
          fe = objective(this,x,se)
!
          if(fe < fr)then
            s(worst) = se ; f(worst) = fe
          else
            s(worst) = sr ; f(worst) = fr
          endif
!
        else
!
!! se is inner or outer contraction point : S = s(best) + gam * (sr - s(best))
!
          si   = 0.5D0 * (s(best) + sr)
          fi = objective(this,x,si)
!
!!  as S = s(best) - gam * (sr - s(best))
!
          so   = 1.5D0 * s(best) - 0.5D0 * sr
          fo = objective(this,x,so)
!
          if(fi<fo)then
            s(worst) = si ; f(worst) = fi
          else
            s(worst) = so ; f(worst) = fo
          endif
        endif
!
        s(0) = s(1)*s(1) + s(2)*s(2)
        conv = abs(conv - s(0))/s(0)
        if(conv < thre) EXIT
        conv = s(0)
!
      enddo
!
      this%strategy = s(best)
!
    end subroutine linear_search
!
    pure function objective(this,x,s) result(res)
    class(rps_hmm),intent(in)        :: this
    integer,intent(in)               :: x(:,:)
    real(REAL64),intent(in)          :: s
    integer                          :: wl(3)
    real(REAL64)                     :: res
      wl  = this%test(x,s)
      res = - REAL( wl(1)-wl(2), REAL64 ) &
     &    / REAL( wl(1)+wl(2)+wl(3), REAL64 )
    end function objective
!
  end subroutine rps_hmm_fit
!
  pure function rps_hmm_test(this,x,strategy) result(res)
  class(rps_hmm),intent(in)        :: this
  integer,intent(in)               :: x(:,:)
  real(REAL64),intent(in),optional :: strategy
  integer                          :: res(3)
  real(REAL64)                     :: pred(this%K)
  real(REAL64)                     :: a(this%d,2)
  real(REAL64)                     :: s
  integer                          :: n,d,wl
  integer                          :: i,j,k
!
    res   = 0
!
    d     = this%d
    n     = SIZE(x,2)
    s     = this%strategy
    if(PRESENT(strategy)) s = strategy
!
    pred = 0D0
    do j=1,this%K
      do i=1,d
        pred(j) = pred(j) + this%phi(i)%n_data() * this%phi(i)%pdf( [j,x(2:,1)] )
      enddo
    enddo
!
    wl     = rps_versus( &
   &           rps_decision( pred(:), s ), &
   &           rps_gethand( x(1,1) ) )
    if( wl==rpsHAND_WIN )  res(1) = res(1) + 1
    if( wl==rpsHAND_LOSE ) res(2) = res(2) + 1
    if( wl==rpsHAND_DRAW ) res(3) = res(3) + 1
!
    do i=1,this%d
      a(i,2) = this%phi(i)%n_data() * this%phi(i)%pdf( x(:,1) )
    enddo
    a(:,1) = a(:,2) / SUM( a(:,2) )
!
    do k=1,n-1
!
      pred = 0D0
      do j=1,this%K
        do i=1,d
          pred(j) = pred(j) + SUM( this%A(i,:) * a(:,1) ) * this%phi(i)%pdf( [j,x(:,k)] )
        enddo
      enddo
!
      wl     = rps_versus( &
     &           rps_decision( pred(:), s ), &
     &           rps_gethand( x(1,k+1) ) )
      if( wl==rpsHAND_WIN )  res(1) = res(1) + 1
      if( wl==rpsHAND_LOSE ) res(2) = res(2) + 1
      if( wl==rpsHAND_DRAW ) res(3) = res(3) + 1
!
      do i=1,this%d
        a(i,2) = SUM( this%A(i,:) * a(:,1) ) * this%phi(i)%pdf( x(:,k+1) )
      enddo
      a(:,1) = a(:,2) / SUM( a(:,2) )
!
    enddo
!
  end function rps_hmm_test
!
  function rps_hmm_pred(this,x) result(res)
  class(rps_hmm),intent(in) :: this
  integer,intent(in)        :: x(:,:)
  real(REAL64)              :: res(this%K)
  real(REAL64)              :: a(this%d,2)
  integer                   :: n
  integer                   :: i,j
!
    res(:) = 0D0
    n      = SIZE(x,2)
!
!   get_alpha
!   PRML (13.37)
!   PRML (13.36) -> (13.59)
!
    do i=1,this%d
      a(i,2) = this%phi(i)%n_data() * this%phi(i)%pdf( x(:,1) )
    enddo
    a(:,1) = a(:,2) / SUM( a(:,2) )
!
    do j=1,n-1
      do i=1,this%d
        a(i,2) = SUM( this%A(i,:) * a(:,1) ) * this%phi(i)%pdf( x(:,j+1) )
      enddo
      a(:,1) = a(:,2) / SUM( a(:,2) )
    enddo
!
    do i=1,this%d
      a(i,2) = SUM( this%A(i,:) * a(:,1) )
    enddo
!
    do j=1,this%K
!
      do i = 1,this%d
        res(j) = res(j) + a(i,2) * this%phi(i)%pdf( [j,x(:,n)] )
      enddo
!
    enddo
!
    res(:) = res(:) / SUM( res(:) )
!
  end function rps_hmm_pred
!
  pure elemental function rps_hmm_log_likelihood(this) result(res)
  class(rps_hmm),intent(in) :: this
  real(REAL64)                    :: res
!
    res = this%L
!
  end function rps_hmm_log_likelihood
!
  pure elemental function rps_hmm_winning_ratio(this) result(res)
  class(rps_hmm),intent(in) :: this
  real(REAL64)              :: res
!
    res = REAL( this%wld(1), REAL64 ) / REAL( this%wld(1)+this%wld(2), REAL64 )
!
  end function rps_hmm_winning_ratio
!
  pure elemental function rps_hmm_n_episode(this) result(res)
  class(rps_hmm),intent(in) :: this
  integer                   :: res
!
    if(.not.ALLOCATED(this%phi))then
      res = 0
    else
      res = this%phi(1)%n_episode()
    endif
!
  end function rps_hmm_n_episode
!
  pure elemental function rps_hmm_n_dimension(this) result(res)
  class(rps_hmm),intent(in) :: this
  integer                   :: res
!
    res = this%D
!
  end function rps_hmm_n_dimension
!
  subroutine rps_hmm_show_status(this,dev)
  class(rps_hmm),intent(in)   :: this
  integer,intent(in),optional :: dev
  integer                     :: i,ld
!
    ld = STDOUT ; if( PRESENT(dev) ) ld = dev
!
    write(ld,'(2(A,f16.9))',ERR=100) '  log_likelihood = ',this%log_likelihood(), &
   &                                 '  winning_ratio = ',this%winning_ratio()
    do i=1,this%d
      write(ld,'(*(f16.9))',ERR=100) this%A(:,i)
    enddo
    write(ld,'(a)',ERR=100) ''
!
    do i=1,this%d
      call this%phi(i)%show_status(dev)
    enddo
!
100 RETURN
!
  end subroutine rps_hmm_show_status
!
  subroutine rps_hmm_clear(this)
  class(rps_hmm),intent(inout)   :: this
    this%D      = 0
    this%K      = 9
    this%wld    = 0
    this%strategy = 0D0
    this%L        =  -HUGE(0D0)
    if( ALLOCATED(this%A) )         DEALLOCATE( this%A )
    if( ALLOCATED(this%c) )         DEALLOCATE( this%c )
    if( ALLOCATED(this%phi) )       DEALLOCATE( this%phi )
    if( ALLOCATED(this%traceback) ) DEALLOCATE( this%traceback )
  end subroutine rps_hmm_clear
!
  subroutine rps_hmm_destroy(this)
  type(rps_hmm),intent(inout)   :: this
    call this%clear()
  end subroutine rps_hmm_destroy
!
end module mod_hmm
