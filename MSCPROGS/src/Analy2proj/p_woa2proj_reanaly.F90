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
! 2) Specify depth levels to which you want to interpolate in the file
! "depthlevels.in". Here you specify the number of depth levels on the first
! line, then the actual depthlevels on the following lines. Depth levels must be
! increasing (down is positive direction).
!
! LB  05052015 -   variables
!---------------------------------------------------------------------------------
program p_woa2proj
   use mod_parameters
   use mod_xc, only: idm,jdm, xcspmd
   use mod_grid, only : get_grid, depths, plon, plat, qlat, qlon
   use mod_spline_calc
   use netcdf
   use mod_rotate
   use mod_za , only: zaiost, zaioempty
   use mod_read_netcdf
   use m_read_mean_ssh
   use m_nearestpoint
   implicit none
   integer,parameter :: itest=467,jtest=701
   integer,parameter :: Lmax=200000
   character*80  :: Edir,Fdir,Finf,Fdep
   character*100 :: Fmod,Fobs
   character*5   :: Jname,cvar

   ! for reading the raw field
   real,allocatable,    dimension(:)  :: modpre
   real,allocatable,dimension(:,:,:)  :: modtem,modsal,moduu,modvv
    
   integer      :: kdm,Jdy
   integer      :: k,nc0,nc1

   ! for reading the observation file
   integer      :: ivar,Olen
   real,allocatable,   dimension(:)  :: Olon,Olat,Odep,Ovar
   real,allocatable,dimension(:)  :: Oip,Ojp
   integer,allocatable,dimension(:)  :: Omask 
   ! find the S0 profiles by searching the observations
   integer       :: S0 
   real,allocatable,   dimension(:) :: Pflon,Pflat
   integer,allocatable,dimension(:) :: Pfi,Pfj,Pfk  ! corresponding the nearest location(ipiv,jpiv)
   integer,allocatable,dimension(:) :: Pflev  ! corresponding 
   logical,allocatable,dimension(:,:) :: modmsk    ! corresponding the nearest location(ipiv,jpiv)
   integer       :: ii,jj,kk
   real          :: a1,a2,a3,a4
   logical       :: lvali
   integer       :: K1,K2                         ! record the depth levels for one profile
   integer       :: Oz
   ! begin to interpolation
   real,allocatable, dimension(:,:,:) :: hy3d   
   real,allocatable,     dimension(:) :: Outvar
   real,allocatable,     dimension(:) :: Invar,depthint   
   logical :: VT0,VT1,VT2  
   integer :: K0     ! deepest layer
   real    :: w0

   ! output the location information and the interpolated result
   character*100 :: Foutinf,Foutvar    ! I/O: 20, 21, ...

   logical :: lev,ex1

   ! convert reading the WOA simulation
   integer :: kz
   real,dimension(:) :: woadepth


   !Turn on/off the output of levitus climatology
   lev=.false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  required the analy.in for defining the raw file and profile location
   Finf='./clim.in';
   Fdep='./depthlevels.in';
   inquire(exist=ex1,file=trim(Finf))
   if(ex1==.true.) then
     open(10,file=trim(Finf))
       read(10,'(a)') Edir 
       read(10,'(a)') Fdir
       read(10,'(2i6)') kdm,Jdy
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
   Fmod='./analysis-'//trim(Jname)//'.nc'
   print *,trim(Fmod)
   inquire(exist=ex1,file=trim(Fmod))
   if(ex1==.false.) then
     stop
   endif

   ! Init file io
   call xcspmd()
   call zaiost()
   ! Get bathymetry - lon/lat fields and sets up bigrid
   call get_grid ()
   allocate(modmsk(idm,jdm))
   modmsk(:,:)=.false.
   do ii=1,idm
     do jj=1,jdm
       if(depths(ii,jj)>0.or.depths(ii,jj)<20000) then
         modmsk(ii,jj)=.true.
       endif
     end do
   end do

   call openNCfile(Fmod,nc0)
     !print *, 'begin to read kdm=', kdm
     call readNCDim(nc0,'kz',kdm)
     allocate(modpre(kdm+1))
     allocate(modtem(idm,jdm,kdm),modsal(idm,jdm,kdm))
     allocate(moduu(idm,jdm,kdm),modvv(idm,jdm,kdm))
     allocate(hy3d(idm,jdm,kdm))
     select case (kdm)
       case (17)
         modpre(2:kdm+1)=(/0, 5, 10, 15, 20, 30, 50, 75, 100, 125, 150, &
               200, 250, 300, 400, 500, 700/)
       case (23)
         modpre(2:kdm+1)=(/5, 10, 20, 30, 50, 75, 100,125, 150, 200, 250,&
               300, 400, 500, 600, 700, 800,&
               900, 1000, 1250, 1500, 2000, 3000/)
       case (24)
         modpre(2:kdm+1)=(/5, 10, 15, 20, 30, 50, 75, 100, 125, 150, 200,&
               250, 300, 400, 500, 600, 700, 800,&
               900, 1000, 1250, 1500, 2000, 3000/)
       case (26)
         modpre(2:kdm+1)=(/5, 10, 20, 40, 60, 80, 100, 120, 140, 160,180, &
               200, 250, 300, 350, 400, 450, 500, 600, 700, 800,&
               900, 1000, 1500, 2000, 3000/)
     end select
     !print '(30F10.1)',(modpre(kk),kk=1,kdm)
     call readNC3d(nc0,'u_z',idm,jdm,kdm,moduu)
     call readNC3d(nc0,'v_z',idm,jdm,kdm,modvv)
     call readNC3d(nc0,'temp_z',idm,jdm,kdm,modtem)
     call readNC3d(nc0,'salt_z',idm,jdm,kdm,modsal)
     modpre(1)=0
   call closeNCfile(nc0)

   Foutinf=trim(Edir)//'clim-proflocat-'//trim(Jname)//'.txt'
   open(20,file=trim(Foutinf))

!!!!!!
! begin to read the related observations file
   do ivar=1,2
     if(ivar==1) then
        cvar='TEM'
       do kk=1,kdm
         hy3d(:,:,kk)=modtem(:,:,kk)
       end do
     else
       cvar='SAL'
       print *,trim(cvar)
       do kk=1,kdm
         hy3d(:,:,kk)=modsal(:,:,kk)
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
    allocate(Omask(Olen),Pflon(Olen),Pflat(Olen),Pfi(Olen),Pfj(Olen),Pfk(Olen),Pflev(Olen))
    Omask(1:Olen)=0;  Pflev=0; K1=0;
    do k=1,Olen
      if(k==1) then
        ii=Oip(k);   jj=Ojp(k);
        call nearestpoint(plon,plat,idm,jdm,Olon(k),Olat(k),ii,jj,a1,a2,a3,a4,modmsk,lvali)
        if(lvali) then
          S0=1 
          Pflon(S0)=Olon(k); Pflat(S0)=Olat(k) 
          Pfi(S0)=ii; Pfj(S0)=jj;  Pfk(S0)=kdm
          Omask(k)=S0
          K1=k;  K2=k; Pflev(S0)=K2-K1+1
        endif
      else
        if(abs(Olon(k)-Pflon(S0))<0.001.and.abs(Olat(k)-Pflat(S0))<0.001.and.k<Olen) then
          if(K1>0) then
            K2=k;  Pflev(S0)=K2-K1+1
            !print '(a,i5,a,f7.2)', 'saving ', S0, ' profile for depth of ',Odep(k)
          endif
        else
          if(K1>0.and.K2>=K1) then
          ! step one: output the depth levels into file
            Oz=K2-K1+1
            allocate(Outvar(Oz),depthint(kdm),Invar(kdm))
            open(10,file=trim(Fdep))
            write(10,'(i3,a)') Oz,'                 # Number of z levels'
            do kk=K1,K2
              write(10,'(F6.1)') Odep(kk)
            end do
            close(10)
          ! step two: spline the simulation into the observation profile
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call spline_calc_ini_fromfile('spline')
            ii=Pfi(S0);  jj=Pfj(S0)
            depthint(:)=modpre(2:kdm+1)
            Invar(:)=hy3d(ii,jj,1:kdm)
            ! controling the values to be smoothly and removeing the anomaly value
            do kk=2,kdm
              if(kk<Pfk(S0)) then
                VT0=abs(Invar(kk))>40
                VT1=Invar(kk-1)>-2.and.Invar(kk-1)<40
                VT2=Invar(kk+1)>-2.and.Invar(kk+1)<40
                if(VT0.and. VT1.and.VT2) then
                  w0=(depthint(kk)-depthint(kk-1))/(depthint(kk+1)-depthint(kk-1))
                  Invar(kk)=w0*Invar(kk+1)+(1-w0)*Invar(kk-1)
                endif         
              else
                VT0=abs(Invar(kk))>40
                if(VT0) then
                  Invar(kk)=undef
                  if(Pfk(S0)>=kk-1) Pfk(S0)=kk-1
                endif
              endif       
            end do
            kk=1;
            do while(kk<kdm)
              VT0=abs(Invar(kk))>40
              VT1=depthint(kk+1)-depthint(kk)<0.5
              if(VT0) then
                kk=kdm
              else
                K0=kk
                if(VT1) then
                  kk=kdm
                endif
              endif 
              kk=kk+1
            end do
            call spline_calc_1d(Invar(1:K0),depthint(1:K0),Outvar,ndeep,K0)
          ! step three: write the simulation profile into txt-file
            print *,'deal with ',S0,' profile of '//trim(cvar),K0,Pfk(S0)
         !   do kk=1,kdm
         !     if(kk<=K0) then
         !       print '(a8,4i5,F9.3,F8.1,F9.3)',trim(cvar),S0,ii,jj,kk,Invar(kk),depthint(kk),modtem(ii,jj,kk)
         !     endif
         !   end do
            if(K0<Pfk(S0)) Pfk(S0)=K0
            do kk=1,Oz
              if(deeps(kk)<=depthint(Pfk(S0))) then
                 write(21,'(a5,i6,2F9.3,2i5,F7.1,2F9.3)') trim(cvar),S0,Pflon(S0),Pflat(S0),Pfi(S0),Pfj(S0),deeps(kk),Outvar(kk),Ovar(K1+kk-1) 
              endif
            end do

            deallocate(Outvar,depthint,Invar)

          endif
          ! change to deal with another profile
          ii=Oip(k);   jj=Ojp(k);
          call nearestpoint(plon,plat,idm,jdm,Olon(k),Olat(k),ii,jj,a1,a2,a3,a4,modmsk,lvali)
          if(lvali) then
            S0=S0+1
            Pflon(S0)=Olon(k); Pflat(S0)=Olat(k) 
            Pfi(S0)=ii; Pfj(S0)=jj;  Pfk(S0)=kdm 
            Omask(k)=S0
            K1=k; K2=k;  Pflev(S0)=K2-K1+1
          endif 
        endif  
      endif 
!     print '(4F10.3,3i6)', Ovar(k),Oip(k),Ojp(k),Odep(k),S0,ii,jj
   enddo

   if(S0>0) then
     do kk=1,S0
       write(20,'(a8,i5,2F9.3,3i5)') trim(cvar),kk,Pflon(kk),Pflat(kk),Pfi(kk),Pfj(kk),Pflev(kk) 
     end do
   endif

   deallocate(Omask,Pflon,Pflat,Pfi,Pfj,Pfk,Pflev)
   deallocate(Olon,Olat,Oip,Ojp,Odep,Ovar)

   close(21)
  end do   ! observation type cycle

  deallocate(hy3d);
  deallocate(modpre,modtem,modsal,moduu,modvv,modmsk);
  close(20)
end
