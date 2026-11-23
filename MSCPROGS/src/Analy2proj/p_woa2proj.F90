!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine creates NetCDF files from (various) hycom files. The netcdf
! files hold the data on a projection grid at specified depthlevels. Routine 
! works on several HYCOM files; restart, archive, nersc weekly average and nersc
! daily average.
!
! -------------------------------------------------------------------------
! Preparation for data processing:
!
!
! LB  05052015 -   variables
!---------------------------------------------------------------------------------
program p_woa2proj
   use netcdf
   use mod_spline_calc
   use mod_read_netcdf
   use woa13_grid
   use m_read_mean_ssh
   use m_nearestpoint
   implicit none
   integer,parameter :: itest=600,jtest=555
   integer,parameter :: Lmax=200000
   character*80  :: Edir,Fdir,Finf
   character*100 :: Fwoa,Fobs
   character*5   :: Jname,cvar

   ! for reading the raw field
   real,allocatable,    dimension(:)  :: woapre
   real,allocatable,dimension(:,:,:)  :: woatem,woasal
   real,allocatable,dimension(:,:)  :: plon,plat,plev 
    
   integer      :: idm,jdm,kdm,Jdy
   integer      :: k,nc0,nc1

   ! for reading the observation file
   integer      :: ivar,Olen
   real,allocatable,   dimension(:)  :: Olon,Olat,Odep,Ovar
   real,allocatable,dimension(:)  :: Oip,Ojp
   integer,allocatable,dimension(:)  :: Omask 
   ! find the S0 profiles by searching the observations
   integer       :: S0 
   real,allocatable,   dimension(:) :: Pflon,Pflat
   integer,allocatable,dimension(:) :: Pfi,Pfj  ! corresponding the nearest location(ipiv,jpiv)
   integer,allocatable,dimension(:) :: Pflev  ! corresponding 
   logical,allocatable,dimension(:,:) :: woamsk    ! corresponding the nearest location(ipiv,jpiv)
   integer       :: ii,jj,kk
   real          :: a1,a2,a3,a4
   logical       :: lvali
   integer       :: K1,K2                         ! record the depth levels for one profile
   integer       :: Oz
   ! begin to interpolation
   real,allocatable, dimension(:,:,:) :: hy3d   
   real,allocatable,     dimension(:) :: Outdeep,Outvar
   logical :: VT0,VT1,VT2  
   integer :: K0     ! deepest layer
   real    :: w0

   ! output the location information and the interpolated result
   character*100 :: Foutinf,Foutvar    ! I/O: 20, 21, ...

   logical :: lev,ex1

   ! convert reading the WOA simulation
   integer :: kz
   real,allocatable,dimension(:) :: woadepth

   character(512) :: options
   integer, external :: iargc




   if (iargc() /= 1) then
     print *, 'Usage: clim2proj <Julianday>'
     print *, '       clim2proj -h'
     print *, 'Options:'
     print *, '  -h -- describe parameter fie format'
   else
     call getarg(1, options)
     if (trim(options) == "-h") then
       print *,'Interpolation WOA13 into the observed profile'
       stop
     else 
       read(options,*) Jdy
       print *,'Julianday: ',Jdy
     end if
   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  required the analy.in for defining the raw file and profile location
   Finf='./clim.in';
   inquire(exist=ex1,file=trim(Finf))
   if(ex1==.true.) then
     open(10,file=trim(Finf))
       read(10,'(a)') Edir 
       read(10,'(a)') Fdir
     close(10)
   else
     print *,'File '//trim(Finf)//' does not exist. '
     stop('Broken!')
   endif
   print *, kdm,Jdy
   inquire(exist=ex1,file=trim(Edir))
   if(ex1==.false.) then
     call system('mkdir '//trim(Edir))
   endif
   write(Jname,'(i5)') Jdy
   call read_woa13gridLen(idm,jdm,kdm)
   allocate(plon(idm,jdm),plat(idm,jdm),woadepth(kdm));
   call read_woa13grid(plat,plon,woadepth,idm,jdm,kdm);
   allocate(woatem(idm,jdm,kdm),woasal(idm,jdm,kdm))
   call read_woa13annual('woa13_decav_t00_04v2.nc',woatem,idm,jdm,kdm,0)
   call read_woa13annual('woa13_decav_s00_04v2.nc',woasal,idm,jdm,kdm,1)
 !  do kk=1,kdm
 !    ii=itest;  jj=jtest
 !    print '(2F10.3,i5)',woatem(ii,jj,kk),woasal(ii,jj,kk),kk
 !  end do
   ! obtain the mask for WOA13
   woamsk=.false.
   allocate(woamsk(idm,jdm),plev(idm,jdm))
   plev=woatem(:,:,1)
   where (plev<100)
     woamsk=.true.
   endwhere
    
  
   allocate(woapre(kdm+1))
   allocate(hy3d(idm,jdm,kdm))

   Foutinf=trim(Edir)//'clim-proflocat-'//trim(Jname)//'.txt'
   open(20,file=trim(Foutinf))

!!!!!!
! begin to read the related observations file
   do ivar=1,2
     if(ivar==1) then
        cvar='TEM'
       do kk=1,kdm
         hy3d(:,:,kk)=woatem(:,:,kk)
       end do
     else
       cvar='SAL'
       print *,trim(cvar)
       do kk=1,kdm
         hy3d(:,:,kk)=woasal(:,:,kk)
       end do
     endif
     Foutvar=trim(Edir)//'clim-prof-'//trim(cvar)//trim(Jname)//'.txt'
     open(21,file=Foutvar)

    Fobs=trim(Fdir)//'observations-'//trim(cvar)//'-'//trim(Jname)//'.nc'
    print *,trim(Fobs)
    allocate(Olon(Lmax),Olat(Lmax),Oip(Lmax),Ojp(Lmax),Odep(Lmax),Ovar(Lmax))
    call openNCfile(Fobs,nc1)
      call readNC1vec(nc1,'d',2,Ovar,Olen)
      call readNC1vec(nc1,'lon',2,Olon,Olen)
      call readNC1vec(nc1,'lat',2,Olat,Olen)
      call readNC1vec(nc1,'depth',2,Odep,Olen)
      call readNC1vec(nc1,'ipiv',1,Oip,Olen)
      call readNC1vec(nc1,'jpiv',1,Ojp,Olen)
    call closeNCfile(nc1)
    allocate(Omask(Olen),Pflon(Olen),Pflat(Olen),Pfi(Olen),Pfj(Olen),Pflev(Olen))
    Omask(1:Olen)=0;  Pflev=0; K1=0;
    do k=1,Olen
      if(k==1) then
        print *,'near: ...'
        call woa13_nearpoint(plat,plon,woamsk,idm,jdm,Olon(k),Olat(k),ii,jj)
        if(ii>0.and.jj>0) then
          S0=1;  Omask(k)=S0
          Pflon(S0)=Olon(k); Pflat(S0)=Olat(k) 
          Pfi(S0)=ii; Pfj(S0)=jj;  !Pflev(S0)=Odep(k)
          K1=k;  K2=k;  Pflev(S0)=1
          !print '(a,i5,a,f7.2)', 'saving ', S0, ' profile for depth of ',Odep(k)
        endif
      else
        if(abs(Olon(k)-Pflon(S0))<0.001.and.abs(Olat(k)-Pflat(S0))<0.001) then
!          print '(a,i5,a,f7.2)', 'saving ', S0, ' profile for depth of ',Odep(k)
          K2=k;  Pflev(S0)=K2-K1+1
          if(K1>0.and.K2>=K1.and.k==Olen) then
            Oz=K2-K1+1
            print *,'deal with ',S0,' profile of '//trim(cvar),Oz
            allocate(Outdeep(Oz),Outvar(Oz))
            Outdeep=Odep(K1:K2) 
            call woa_vertical_inter(hy3d(Pfi(S0),Pfj(S0),:),woadepth,kdm,100.0,Outdeep,Outvar,Oz)
            do kk=1,Oz
               write(21,'(a5,i6,2F9.3,2i5,F7.1,2F9.3)') trim(cvar),S0,Pflon(S0),Pflat(S0),Pfi(S0),Pfj(S0),Outdeep(kk),Outvar(kk),Ovar(K1+kk-1) 
            end do
            deallocate(Outdeep,Outvar)
          endif
        else
          if(K1>0.and.K2>=K1) then
            Oz=K2-K1+1
            allocate(Outdeep(Oz),Outvar(Oz))
            print *,'deal with ',S0,' profile of '//trim(cvar),Oz
            Outdeep=Odep(K1:K2) 
            call woa_vertical_inter(hy3d(Pfi(S0),Pfj(S0),:),woadepth,kdm,100.0,Outdeep,Outvar,Oz)
            if (Oz>0) then 
              do kk=1,Oz
                 write(21,'(a5,i6,2F9.3,2i5,F7.1,2F9.3)') trim(cvar),S0,Pflon(S0),Pflat(S0),Pfi(S0),Pfj(S0),Outdeep(kk),Outvar(kk),Ovar(K1+kk-1) 
              end do
            endif
            deallocate(Outdeep,Outvar)
          endif
          print *,'New profile: ',Olen,'(',k,')'
          ! change to deal with another profile
          call woa13_nearpoint(plat,plon,woamsk,idm,jdm,Olon(k),Olat(k),ii,jj)
          if(ii>0.and.jj>0) then
            S0=S0+1
            Pflon(S0)=Olon(k); Pflat(S0)=Olat(k) 
            Pfi(S0)=ii; Pfj(S0)=jj;   
            Omask(k)=S0
            K1=k; K2=k;  Pflev(S0)=K2-K1+1
          endif 
        endif  
      endif 
   enddo

   if(S0>0) then
     do kk=1,S0
       write(20,'(a8,i5,2F9.3,3i5)') trim(cvar),kk,Pflon(kk),Pflat(kk),Pfi(kk),Pfj(kk),Pflev(kk) 
     end do
   endif

   deallocate(Omask,Pflon,Pflat,Pfi,Pfj,Pflev)
   deallocate(Olon,Olat,Oip,Ojp,Odep,Ovar)

   close(21)
  end do   ! observation type cycle

  deallocate(hy3d);
  deallocate(woapre,woatem,woasal,woamsk,woadepth);
  close(20)
end
