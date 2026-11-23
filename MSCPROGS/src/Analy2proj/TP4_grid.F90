module TP4_grid

  public :: TP4_vertical_inter

contains
  
  subroutine TP4_vertical_inter(Tprof,Tlevs,Kdm,undef2,Olevs,Outv,Klen)
    use mod_spline_calc
    implicit none
    integer, intent(in) :: Kdm
    real,    intent(in) :: undef2
    real,dimension(:), intent(in) :: Tprof,Tlevs
    real,dimension(:), intent(in) :: Olevs
    real,dimension(:),intent(out) :: Outv
    integer,          intent(out) :: Klen
    ! temperary variables
    real,allocatable, dimension(:) :: depthint,Invar
    integer :: k
    integer :: k0      ! used by model stratification
    integer :: kk      ! used by observation level
    logical :: VT0,VT1,VT2

    Klen=0
    allocate(depthint(Kdm),Invar(Kdm))
    depthint=Tlevs(1:Kdm)
    Invar=Tprof(1:Kdm)
    ! first to check the model stratification
    k=1; k0=0
    do while (k<Kdm)
      ! valid profile 
      VT0=Invar(k)>50.or.Invar(k)<-4
      if(VT0) then
        k=Kdm
      else
        k0=k
        ! reach to bottom
        VT1=depthint(k+1)-depthint(k)<1.0
        if(VT1) then
          k=Kdm
        endif 
      endif 
      k=k+1
    end do 
    if (k0<1) then
      return
    endif
!    do k=1,k0
!      print '(i5,2F10.3)', k,depthint(k),Invar(k)
!    end do
 

    ! second to check the observed profile 
    kk=0
    if (maxval(depthint(1:k0))>=maxval(Olevs)) then
      kk=size(Olevs)
    else 
      k=1
      do while (k<=size(Olevs))
        if (Olevs(k)>=maxval(depthint(1:k0))) then
          k=size(Olevs)
        else
          kk=k
        endif
        k=k+1
      end do
    endif
    if (kk<1) then
      return
    endif
!    do k=1,kk
!      print '(i5,F10.3)', k,Olevs(k)
!    end do

    Klen=kk
    call spline_calc_ini_frominput('spline',Olevs(1:kk),kk)
    call spline_calc_1d(Invar(1:k0),depthint(1:k0),Outv,kk,k0)
    do k=1,kk
      print '(2F10.3)', Olevs(k),Outv(k)
    end do
    deallocate(depthint,Invar)
  end subroutine


end module
