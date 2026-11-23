module woa13_grid
  use netcdf
  use mod_read_netcdf
 
  private
  real,parameter :: Mingrid=.25  ! minimal grid for degree 
  real,parameter :: Mindis=30    ! minimal grid for kilometer  
  character*200, parameter :: WOADIR='/work/xiejp/work_2015/Data/WOA13/' 

  public :: read_woa13gridLen, read_woa13grid,read_woa13annual,woa13_nearpoint
  public :: woa_vertical_inter

contains
   
  subroutine read_woa13gridLen(Nlon,Nlat,Nlev)
    implicit none
    integer, intent(out) :: Nlat,Nlon,Nlev
    !-------
    integer :: ncid
    !-------
    character*200 :: ncfile
    ncfile=trim(WOADIR)//'woa13_decav_t00_04v2.nc'
    call openNCfile(trim(ncfile),ncid)
    call readNCDim(ncid,'lat',Nlat)
    call readNCDim(ncid,'lon',Nlon)
    call readNCDim(ncid,'depth',Nlev)
!    print *, Nlat, Nlon,Nlev
    call closeNCfile(ncid)
  end subroutine

  subroutine read_woa13grid(Wlat,Wlon,Wdep,Nlon,Nlat,Nlev)
    implicit none
    integer, intent(in) :: Nlat,Nlon,Nlev
    real,dimension(:,:),intent(out) :: Wlat,Wlon
    real,dimension(:),  intent(out) :: Wdep
    !-------
    integer :: ncid
    integer :: varid
    integer,dimension(nf90_max_var_dims) :: vardims
    !-------
    integer :: N1,N2,N3
    real,allocatable,dimension(:,:) :: var0,var1,var2
    integer i,j,k
    character*200 :: ncfile
    ncfile=trim(WOADIR)//'woa13_decav_t00_04v2.nc'
    call openNCfile(trim(ncfile),ncid)
    call handle_err(nf90_inq_varid(ncid,'lat_bnds',varid),.true.)
    call handle_err(nf90_inquire_variable(ncid,varid,dimids=vardims),.true.)
    call handle_err(nf90_inquire_dimension(ncid,vardims(1),len=N1))
    call handle_err(nf90_inquire_dimension(ncid,vardims(2),len=N2))
!    print *, 'lat_bnds:',N1,N2
    allocate(var0(N1,N2))
    call handle_err(nf90_get_var(ncid,varid,var0),.true.)

    call handle_err(nf90_inq_varid(ncid,'lon_bnds',varid),.true.)
    call handle_err(nf90_inquire_variable(ncid,varid,dimids=vardims),.true.)
    call handle_err(nf90_inquire_dimension(ncid,vardims(1),len=N1))
    call handle_err(nf90_inquire_dimension(ncid,vardims(2),len=N2))
!    print *, 'lon_bnds:',N1,N2
    allocate(var1(N1,N2))
    call handle_err(nf90_get_var(ncid,varid,var1),.true.)

    call handle_err(nf90_inq_varid(ncid,'depth_bnds',varid),.true.)
    call handle_err(nf90_inquire_variable(ncid,varid,dimids=vardims),.true.)
    call handle_err(nf90_inquire_dimension(ncid,vardims(1),len=N1))
    call handle_err(nf90_inquire_dimension(ncid,vardims(2),len=N2))
!    print *, 'depth_bnds:',N1,N2
    allocate(var2(N1,N2))
    call handle_err(nf90_get_var(ncid,varid,var2),.true.)
    do i=1,Nlon
      do j=1,Nlat
        Wlat(i,j)=sum(var0(:,j))/2.0
        Wlon(i,j)=sum(var1(:,i))/2.0
      end do
    end do
    do k=1,N2
       Wdep(k)=var2(2,k)
    end do
    deallocate(var0,var1,var2)
    call closeNCfile(ncid)
  end subroutine

  subroutine read_woa13annual(WOAfile,var3d,Nlon,Nlat,Nlev,Ivar)
    implicit none
    integer, intent(in) :: Nlat,Nlon,Nlev
    integer, intent(in) :: Ivar   ! =0 for temperature; =1 for salinity;
    character(*), intent(in) :: WOAfile
    real,dimension(:,:,:),intent(out) :: var3d
    !-------
    integer :: ncid
    integer :: varid
    integer,dimension(nf90_max_var_dims) :: vardims
    !-------
    integer :: N1,N2,N3
    real,allocatable,dimension(:,:,:) :: var0
    integer i,j,k
    character*20 :: varnam 
    character*200 :: ncfile 
    ncfile=trim(WOADIR)//trim(WOAfile)
    call openNCfile(trim(ncfile),ncid)
    select case (Ivar)
      case (0)
        varnam='temp'
        call handle_err(nf90_inq_varid(ncid,'t_an',varid),.true.)
      case (1)
        varnam='saln'
        call handle_err(nf90_inq_varid(ncid,'s_an',varid),.true.)
    end select
    call handle_err(nf90_inquire_variable(ncid,varid,dimids=vardims),.true.)
    call handle_err(nf90_inquire_dimension(ncid,vardims(1),len=N1))
    call handle_err(nf90_inquire_dimension(ncid,vardims(2),len=N2))
    call handle_err(nf90_inquire_dimension(ncid,vardims(3),len=N3))
    print *, trim(varnam)//':',N1,N2,N3
    allocate(var0(N1,N2,N3))
    call handle_err(nf90_get_var(ncid,varid,var0),.true.)
    var3d=var0
    deallocate(var0)
    call closeNCfile(ncid)
  end subroutine



! works for the formal grids
  subroutine woa13_nearpoint(Wlat,Wlon,Wmsk,Nlon,Nlat,olon,olat,ipiv,jpiv)
    implicit none
    integer,intent(in) :: Nlon,Nlat
    logical,dimension(:,:), intent(in) :: Wmsk
    real,   dimension(:,:), intent(in) :: Wlat,Wlon
    real,   intent(in) :: olon,olat
    integer,intent(out) :: ipiv,jpiv
    real, external :: spherdist
    !-------------------
    integer :: i,j,k
    integer :: i1,i2,j1,j2
    integer :: Vali(Nlon,Nlat)
    real    :: mgrd,tmpvar
    ! default values:
    ipiv=0;  jpiv=0
    j=1
    do while (j<Nlat)
      if (olat<=Wlat(1,1)) then
        j1=1; j2=1; j=Nlat
      elseif (olat>=Wlat(1,j).and.olat<=Wlat(1,j+1)) then
        j1=j; j2=j+1;  j=Nlat
      elseif (olat>Wlat(1,Nlat)) then
        j1=Nlat; j2=Nlat
      endif
      j=j+1
    end do
    i=1
    do while (i<Nlon)
      if (olon<Wlon(1,1)) then
        i1=Nlon;  i2=1;  i=Nlon
      elseif (olon>=Wlon(i,1).and.olon<=Wlon(i+1,1)) then
        i1=i;  i2=i+1;  i=Nlon
      elseif (olon>Wlon(Nlon,1)) then
        i1=Nlon;  i2=1;  i=Nlon
      endif
      i=i+1
    end do
    Vali=0; 
    if (i2>i1) then
      where (Wlat>=Wlat(1,j1).and.Wlat<=Wlat(1,j2).and.Wmsk==.true..and.Wlon>=Wlon(i1,1).and.Wlon<=Wlon(i2,1))
        Vali=1
      endwhere
    else
      where (Wlat>=Wlat(1,j1).and.Wlat<=Wlat(1,j2).and.Wmsk==.true..and.(Wlon>=Wlon(Nlon,1).or.Wlon<=Wlon(1,1)))
        Vali=1
      endwhere
    endif
    if (sum(Vali)>0) then
      mgrd=Mindis
      do k=1,4
        select case (k)
          case (1)
            i=i1; j=j1
          case (2)
            i=i2; j=j1
          case (3)
            i=i2; j=j2
          case (4)
            i=i1; j=j2
        end select
        if (Vali(i,j)==1) then
          tmpvar=spherdist(Wlon(i,j),Wlat(i,j),olon,olat)/1000.
          if (tmpvar<mgrd) then
             ipiv=i;  jpiv=j;
             mgrd=tmpvar 
          endif
        endif
      end do
    endif 
!    print *,ipiv,jpiv
!    print '(2F10.3)',olon,olat
!    print '(2F10.3)',Wlon(ipiv,jpiv),Wlat(ipiv,jpiv)
  end subroutine

! extending the first layer to surface
  subroutine woa_vertical_inter(Wprof,Wlevs,Nlev,undef,Olevs,Outv,Klen)
    use mod_spline_calc
    implicit none
    integer, intent(in) :: Nlev
    real,    intent(in) :: undef
    real,dimension(:), intent(in) :: Wprof,Wlevs
    real,dimension(:), intent(in) :: Olevs
    real,dimension(:),intent(out) :: Outv
    integer,          intent(out) :: Klen
    
    ! keep the valid values in the Wprof
    real,allocatable,dimension(:) :: Nprof,Ndeps
    integer :: k0,k,kk
    Klen=0
    allocate(Nprof(Nlev),Ndeps(Nlev));
    k=1; k0=0
    do while (k<=Nlev)
      if (abs(Wprof(k))>undef) then
        k=Nlev
      else
        k0=k0+1
        Nprof(k0)=Wprof(k);  
        Ndeps(k0)=Wlevs(k); 
      endif
      k=k+1 
    end do
    if (k0<1) then
       return
    endif
!    do k=1,k0
!      print '(i5,2F10.3)', k,Ndeps(k),Nprof(k)
!    end do
    ! saving the depths level into file
    kk=0
    if (maxval(Ndeps(1:k0))>=maxval(Olevs)) then
      kk=size(Olevs)
    else
      k=1
      do while(k<=size(Olevs)) 
        if (Olevs(k)>=maxval(Ndeps(1:k0))) then
          k=size(Olevs)
        else
          kk=k
        endif
         k=k+1
      end do
    endif
    ! begin to interpolation
!    call spline_calc_ini_fromfile('spline')
    if (kk<1) then
      return
    endif
    Klen=kk
    call spline_calc_ini_frominput('spline',Olevs(1:kk),kk)
    call spline_calc_1d(Nprof(1:k0),Ndeps(1:k0),Outv,kk,k0)
    do k=1,kk
      print '(2F10.3)', Olevs(k),Outv(k)
    end do
    deallocate(Nprof,Ndeps);
  end subroutine



end module

