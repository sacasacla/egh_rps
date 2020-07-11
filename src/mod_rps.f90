module mod_rps
use,intrinsic :: ISO_FORTRAN_ENV, only : REAL64,STDOUT=>OUTPUT_UNIT
implicit none
private
public :: rps_data
public :: rps_model
public :: rps_versus
public :: rps_hand_index
public :: rps_gethand
public :: rps_transition
public :: rps_decision
public :: rpsHAND_DRAW
public :: rpsHAND_LOSE
public :: rpsHAND_WIN
public :: assignment (=)
!
  real(REAL64),parameter   :: DEF_LAM = 0D0
!
  integer,parameter        :: rpsHAND_ROCK = 1
  integer,parameter        :: rpsHAND_SSRS = 2
  integer,parameter        :: rpsHAND_PPER = 3
!
  integer,parameter        :: rpsHAND_DRAW = 0
  integer,parameter        :: rpsHAND_LOSE = 1
  integer,parameter        :: rpsHAND_WIN  = 2
!
type rps_data
  integer      :: i,y,m,d,hand
  character(1) :: h
end type rps_data
!
type rps_node
  private
  integer,allocatable :: path(:)
  real(REAL64)        :: C = 0D0
contains
  procedure    :: init      => rps_node_init
  final        :: rps_node_destroy
end type rps_node
!
type rps_model
  private
  integer                    :: N = 0             ! Maximum of episode
  integer                    :: K = 9             ! Size of event
  type(rps_node),allocatable :: node(:)
  real(REAL64)               :: lam    = DEF_LAM  ! noise
  real(REAL64)               :: ndata  = 0D0
  real(REAL64)               :: nconst = 0D0
contains
  procedure    :: init        => rps_model_init
  procedure    :: maximize    => rps_model_maximize
  procedure    :: show_status => rps_model_show_status
  procedure    :: pdf         => rps_model_pdf
  procedure    :: logpdf      => rps_model_logpdf
  procedure    :: n_episode   => rps_model_n_episode
  procedure    :: n_data      => rps_model_n_data
  procedure    :: clear       => rps_model_clear
  final        :: rps_model_destroy
end type rps_model
!
  interface assignment (=)
    module procedure rps_model_assign
  end interface assignment (=)
contains
!
  pure elemental function rps_hand_index(d) result(res)
  class(rps_data),intent(in) :: d
  integer                    :: res
    select case(d%h)
    case('r')    ; res = rpsHAND_ROCK
    case('s')    ; res = rpsHAND_SSRS
    case('p')    ; res = rpsHAND_PPER
    case default ; res = 0
    end select
  end function rps_hand_index
!
  pure elemental function rps_transition(d0,d1) result(res)
  class(rps_data),intent(in) :: d0,d1
  integer                    :: res
!
    res = rps_hand_index(d1)
    res = res * 3 - 3
    res = res + rps_hand_index(d0)
!
  end function rps_transition
!
  pure elemental function rps_gethand(hand) result(res)
  integer,intent(in) :: hand
  integer            :: res
!
    res = ( hand - 1 ) / 3 + 1
!
  end function rps_gethand
!
  pure elemental function rps_versus(l,r) result(res)
  integer,intent(in) :: l,r
  integer            :: res
!
    res = MODULO( l - r, 3 )
!
  end function rps_versus
!
  pure function rps_decision(pred,strategy) result(res)
  real(REAL64),intent(in)     :: pred(9),strategy
  real(REAL64)                :: plty(3),prob(3),s1,s2,s
  integer                     :: res
!
    if( strategy<-3D1 )then
      s1    = 1D0
      s2    = 0D0
    elseif( strategy>3D1 )then
      s1    = 0D0
      s2    = 1D0
    else
      s1    = EXP(-strategy)
      s2    = EXP(strategy)
    endif
    s       = 1D0 / ( s1 + s2 )
    s1      = s1 * s
    s2      = s2 * s
!
    plty(1) = pred(7) + pred(8) + pred(9)
    plty(2) = pred(1) + pred(2) + pred(3)
    plty(3) = pred(4) + pred(5) + pred(6)
    prob(1) = plty(3)
    prob(2) = plty(1)
    prob(3) = plty(2)
!
    prob(:) = prob(:) * s1 - plty(:) * s2
!
    res     = MAXLOC( prob, 1 )
!
  end function rps_decision
!
  pure elemental subroutine rps_node_init(this,N)
  class(rps_node),intent(inout)   :: this
  integer,intent(in)               :: N
    call rps_node_destroy(this)
    ALLOCATE( this%path(N) )
  end subroutine rps_node_init
!
  pure subroutine rps_node_destroy(this)
  type(rps_node),intent(inout)   :: this
    if(ALLOCATED(this%path)) DEALLOCATE(this%path)
    this%C = 0D0
  end subroutine rps_node_destroy
!
  pure subroutine rps_model_init(this,N,K,lam)
  class(rps_model),intent(inout)   :: this
  integer,intent(in),optional      :: N,K
  real(REAL64),intent(in),optional :: lam
!
    call this%clear()
!
    if(present(N))   this%N   = N
    if(present(K))   this%K   = K
    if(present(lam)) this%lam = lam
    ALLOCATE( this%node(1) )
    call this%node(1)%init( this%N )
!
  end subroutine rps_model_init
!
  pure elemental subroutine rps_model_assign(this,rhs)
  class(rps_model),intent(inout) :: this
  class(rps_model),intent(in)    :: rhs
!
    call this%clear()
!
    this%N      = rhs%N
    this%K      = rhs%K
    this%lam    = rhs%lam
    this%ndata  = rhs%ndata
    this%nconst = rhs%nconst
    ALLOCATE( this%node(SIZE(rhs%node)) )
    this%node(:) = rhs%node(:)
  end subroutine rps_model_assign
!
  subroutine rps_model_maximize(this,x,gam)
  class(rps_model),intent(inout)   :: this
  integer,intent(in)               :: x(:,:)
  real(REAL64),intent(in),optional :: gam(:)
  integer                          :: i
!
    this%node(:)%C    = -1D0
!
    if(present(gam))then
!
      do i=1,SIZE(x,2)
        call countup_path( this%node, x(:this%N,i), gam(i) )
      enddo
!
    else
!
      do i=1,SIZE(x,2)
        call countup_path( this%node, x(:this%N,i), 1D0 )
      enddo
!
    endif
!
    this%ndata  = SUM( this%node(:)%C, this%node(:)%C>=0D0 )
    call linear_search(this,x)
    this%nconst = rps_model_n_const(this,this%lam)
!
  contains
!
    pure subroutine countup_path(node,ep,c)
    type(rps_node),intent(inout),allocatable :: node(:)
    integer,intent(in)                       :: ep(:)
    real(REAL64),intent(in)                  :: c
    integer                                  :: i,n
!
      do i=1,SIZE(node)
        if( ALL( node(i)%path(:)==ep(:) ) )then
          if(node(i)%C<=0D0) node(i)%C = 0D0
          node(i)%C = node(i)%C + c
          RETURN
        endif
      enddo
!
      do i=1,SIZE(node)
        if(node(i)%C<=0D0)then
          node(i)%path(:)  = ep(:)
          node(i)%C        = c
          RETURN
        endif
      enddo
!
      n    = SIZE(node) + 1
      node = [node(:),node(:)]
      do i=n,SIZE(node)
        node(i)%path(:) = -1
        node(i)%C       = -1D0
      enddo
!
      node(n)%path(:)  = ep(:)
      node(n)%C        = c
!
    end subroutine countup_path
!
    pure function objective(this,x,s) result(res)
    type(rps_model),intent(in) :: this
    integer,intent(in)         :: x(:,:)
    real(REAL64),intent(in)    :: s
    real(REAL64)               :: res,lambda
    integer                    :: i
      res     = 0D0
      lambda  = EXP( s )
      do i=1,SIZE(x,2)
        res = res - this%logpdf(x(:,i),lambda)
      enddo
      if(res/=res) res = HUGE(0D0)
    end function objective
!
    subroutine linear_search(this,x)
    type(rps_model),intent(inout) :: this
    integer,intent(in)            :: x(:,:)
    real(REAL64)                  :: s(0:2),f(2)
    real(REAL64)                  :: sr,fr ! reflect
    real(REAL64)                  :: se,fe ! expand
    real(REAL64)                  :: si,fi ! inner contract
    real(REAL64)                  :: so,fo ! outer contract
    real(REAL64)                  :: conv,thre
    integer                       :: best,worst
    integer                       :: i,maxiter
!
! nelder_mead parameters
!
!    alpha = 1d0
!    beta  = 2d0
!    gam   = 5d-1
!
      s(1)   = 0D0  ; f(1) = objective(this,x,s(1))
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
      if(s(best)>0D0)then
        this%lam = 1D0
      else
        this%lam = EXP( s(best) )
      endif
!
    end subroutine linear_search
!
  end subroutine rps_model_maximize
!
  subroutine rps_model_show_status(this,dev)
  class(rps_model),intent(in) :: this
  integer,intent(in),optional :: dev
  logical                     :: m(SIZE(this%node))
  integer                     :: i,p,ld
!
    ld = STDOUT ; if( PRESENT(dev) ) ld = dev
!
    write(ld,'(2(A,f16.3))',ERR=100) '  n_data  = ',this%ndata,   &
   &                                 '  n_const = ',this%nconst
    write(ld,'(A,G16.3)',ERR=100)    '  lambda  = ',this%lam
!
    m = .TRUE.
    do i=1,SIZE(this%node)
!
      p = MAXLOC(this%node(:)%C,1,m)
      if(p<1) EXIT
      if(this%node(p)%C/this%ndata<=1D-8) EXIT
!
      write(ld,'(2f16.9,*(i4))',ERR=100) this%node(p)%C,this%node(p)%C/this%ndata,this%node(p)%path
      m(p) = .FALSE.
!
    enddo
!
    write(ld,'(a)',ERR=100) ''
!
100 RETURN
!
  end subroutine rps_model_show_status
!
  pure function rps_model_pdf(this,ep,lambda) result(res)
  class(rps_model),intent(in)      :: this
  integer,intent(in)               :: ep(:)
  real(REAL64),intent(in),optional :: lambda
  real(REAL64)                     :: res
    if(this%lam>1D-12)then
      res = EXP( this%logpdf(ep,lambda) )
    else
      if(PRESENT(lambda))then
        res = rps_model_countup(this,ep,lambda) / rps_model_n_const(this,lambda)
      else
        res = rps_model_countup(this,ep,this%lam) / this%nconst
      endif
    endif
  end function rps_model_pdf
!
  pure function rps_model_logpdf(this,ep,lambda) result(res)
  class(rps_model),intent(in)      :: this
  integer,intent(in)               :: ep(:)
  real(REAL64),intent(in),optional :: lambda
  real(REAL64)                     :: res
!
    if(PRESENT(lambda))then
      res = LOG( rps_model_countup(this,ep,lambda) )   - LOG( rps_model_n_const(this,lambda) )
    else
      res = LOG( rps_model_countup(this,ep,this%lam) ) - LOG( this%nconst )
    endif
!
  end function rps_model_logpdf
!
  pure function rps_model_countup(this,ep,lambda) result(res)
  class(rps_model),intent(in) :: this
  integer,intent(in)          :: ep(:)
  real(REAL64),intent(in)     :: lambda
  real(REAL64)                :: tmp(this%N)
  real(REAL64)                :: res
  integer                     :: i
!
    res = 0D0
!
    do i=1,SIZE(this%node(:))
      if(this%node(i)%C<=0D0) CYCLE
      where(this%node(i)%path(:)==ep(:this%N))
        tmp(:) = ( 1D0 - lambda )
      elsewhere
        tmp(:) = lambda
      endwhere
      res = res + this%node(i)%C * PRODUCT( tmp(:) )
    enddo
!
  end function rps_model_countup
!
  pure elemental function rps_model_n_const(this,lambda) result(res)
  class(rps_model),intent(in) :: this
  real(REAL64),intent(in)     :: lambda
  real(REAL64)                :: res
  integer                     :: i,j,fn
!
    res = 0D0
!
    fn = factorial(this%N)
!
    do i=0,this%N
      j = this%N - i
      res = res                                  &
     &    + fn / ( factorial(i) * factorial(j) ) &
     &    * ( 1D0 - lambda )**i                  &
     &    * ( ( this%K - 1 ) * lambda )**j
    enddo
!
    res = res * this%ndata
!
  contains
!
    pure recursive function factorial(n) result(res)
    integer,intent(in) :: n
    integer            :: res
      if(n<=1)then
        res = 1
      else
        res = n * factorial(n-1)
      endif
    end function factorial
!
  end function rps_model_n_const
!
  pure elemental function rps_model_n_episode(this) result(res)
  class(rps_model),intent(in) :: this
  integer                     :: res
!
    res = this%n
!
  end function rps_model_n_episode
!
  pure elemental function rps_model_n_data(this) result(res)
  class(rps_model),intent(in) :: this
  integer                     :: res
!
    res = NINT(this%ndata)
!
  end function rps_model_n_data
!
  pure subroutine rps_model_clear(this)
  class(rps_model),intent(inout)   :: this
    this%N      = 0
    this%K      = 9
    this%lam    = DEF_LAM
    this%ndata  = 0D0
    this%nconst = 0D0
    if( ALLOCATED(this%node) ) DEALLOCATE( this%node )
  end subroutine rps_model_clear
!
  pure subroutine rps_model_destroy(this)
  type(rps_model),intent(inout)   :: this
    call this%clear()
  end subroutine rps_model_destroy
!
end module mod_rps
