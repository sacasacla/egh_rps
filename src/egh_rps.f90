program main
  use mod_rps
  use mod_hmm
implicit none
!
type(rps_hmm)              :: h
type(rps_data),allocatable :: dat(:)
character(256)             :: arg
integer                    :: ndm,nep,nmc
integer                    :: iyear
!
  call GET_COMMAND_ARGUMENT(1,arg)
  call load(arg,dat)
!
  call GET_COMMAND_ARGUMENT(2,arg)
  read(arg,*,iostat=iyear) ndm
  if(iyear/=0) ndm = 2
!
  call GET_COMMAND_ARGUMENT(3,arg)
  read(arg,*,iostat=iyear) nep
  if(iyear/=0) nep = 2
!
  call GET_COMMAND_ARGUMENT(4,arg)
  read(arg,*,iostat=iyear) nmc
  if(iyear/=0) nmc = 1
!
  call get_map(dat,50)
!
  call hmm(h,pack(dat,dat%y<2011),pack(dat,dat%y>2010),d=ndm,n=nep,nmc=nmc)
!
  do iyear=1992,2020
!
    call hmm_test(iyear,h,dat)
!
  enddo
!
  call hmm_test(-1,h,pack(dat,dat%y>2010))
!
contains
!
  subroutine get_map(dat,tau)
  type(rps_data),intent(in) :: dat(:)
  integer,intent(in)        :: tau
  real(8)                   :: x(tau,0:2),z(2,tau,0:2)
  integer                   :: i,j,p
!
    x(:,:) = 0D0
!
    do j=1,tau
!
      do i=j+1,SIZE(dat)
        p      = rps_versus( dat(i)%hand, dat(i-j)%hand )
        x(j,p) = x(j,p) + 1D0
      enddo
!
      x(j,:) = x(j,:) / SUM( x(j,:) )
      x(j,:) = ( x(j,:) - 1D0/3D0 )
!
    enddo
!
    call dct(x(:,0),z(:,:,0))
    call dct(x(:,1),z(:,:,1))
    call dct(x(:,2),z(:,:,2))
!
    do i=1,tau
      print'(*(f9.3))',x(i,:),z(:,i,:)
    enddo
    print*
    FLUSH(6)
!
  end subroutine get_map
!
  subroutine dct(x,z)
  real(8),intent(in)  :: x(:)
  real(8),intent(out) :: z(:,:)
  real(8),parameter   :: pi=ACOS(-1D0)*0.25D0
  real(8)             :: phase,pp
  integer             :: i,j
!
    z = 0D0
    z(1,1) = SUM( x(:) )
!
    do j=2,SIZE(x)
      phase = pi / REAL(j-1,8)
      pp    = 0D0
      do i=1,SIZE(x)
        z(1,j) = z(1,j) + x(i) * COS( pp )
        z(2,j) = z(2,j) + x(i) * SIN( pp )
        pp   = pp + phase
      enddo
    enddo
!
  end subroutine dct
!
  pure subroutine get_fv(dat,x)
  type(rps_data),intent(in) :: dat(:)
  integer,intent(inout)     :: x(:,:)
  integer                   :: i,m
    m = SIZE(x,1)
    do i=1,m
      x(i,1) = rps_transition( dat(m-i+2), dat(m-i+1) )
    enddo
    do i=2,size(dat)-m
      x(:,i)  = [ rps_transition( dat(i+m-1), dat(i+m) ), x(:m-1,i-1) ]
    enddo
  end subroutine get_fv
!
  subroutine hmm(h,dat,test,d,n,nmc)
  type(rps_hmm),intent(inout) :: h
  type(rps_data),intent(in)   :: dat(:),test(:)
  integer,intent(in)          :: d,n,nmc
  type(rps_hmm)               :: htest
  integer                     :: x(n,size(dat)-n)
  integer                     :: x_test(n,size(test)-n)
  integer                     :: i,j
  integer                     :: wtest(3),wl,wb
!
    call get_fv(dat,x)
    call get_fv(test,x_test)
!
    wb   = 0
!
    call htest%init( d )
!
    do j=2,n
!
      do i=1,nmc
!
        call htest%fit(x(:j,:),n_episode=j)
        wtest = htest%test(x_test(:j,:))
        wl    = wtest(1)-wtest(2)
        if(wl>wb)then
          h  = htest
          wb = wl
          print'(7i6,2f24.9,A)',d,j,i,wtest,wl,h%log_likelihood(),h%winning_ratio(),'  *'
        else
          print'(7i6,2f24.9)', d,j,i,wtest,wl,htest%log_likelihood(),h%winning_ratio()
        endif
!
        FLUSH(6)
!
      enddo
!
    enddo
!
    call h%show_status()
!
  end subroutine hmm
!
  subroutine hmm_test(year,h,dat)
  integer,intent(in)          :: year
  type(rps_hmm),intent(inout) :: h
  type(rps_data),intent(in)   :: dat(:)
  integer,allocatable         :: x(:,:)
  logical                     :: mask(SIZE(dat))
  integer                     :: wtest(3),ip,i,n
!
    n = COUNT(dat%y==year) ; if(year<0) n = SIZE(dat) - h%n_episode()
!
    if(n<1) RETURN
!
    if(year<0)then
      mask = .TRUE.
    else
      ip   = MINVAL([(i,i=1,SIZE(dat))],1,dat%y==year) - h%n_episode()
      mask = [(i,i=1,SIZE(dat))]>ip .and. dat%y==year-1 .or. dat%y==year
    endif
!
    allocate( x(h%n_episode(),n) )
!
    call get_fv(PACK(dat,mask),x)
    wtest = h%test(x)
!
    if(year<0)then
      print'(8X,4i8,f16.9)',wtest,SIZE(dat),                       &
   &                     REAL(wtest(1),8)/REAL(wtest(1)+wtest(2),8)
    else
      print'(5i8,f16.9)',year,wtest,COUNT(dat%y==year), &
   &                     REAL(wtest(1),8)/REAL(wtest(1)+wtest(2),8)
    endif
!
  end subroutine hmm_test
!
  subroutine load(path,dat)
  character(*),intent(in)                  :: path
  type(rps_data),intent(inout),allocatable :: dat(:)
  type(rps_data)                           :: buff(10000)
  integer                                  :: i,n
!
    open(10,file=path)
    read(10,'(a)')
    read(10,'(a)')
    do i=1,SIZE(buff)
      read(10,*,end=100) buff(i)%i,buff(i)%y,buff(i)%m,buff(i)%d,buff(i)%h
      buff(i)%hand = rps_hand_index(buff(i))
    enddo
!
100 close(10)
    n = i - 1
!
    if(allocated(dat)) deallocate(dat)
    allocate(dat(n))
!
    dat(:n) = buff(:n)
!
  end subroutine load
!
end program main
