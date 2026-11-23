!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine creates NetCDF files from the daily average files. The netcdf
! files hold the data on a regular lon/lat grid at specified depthlevels. The
! grid and depth levels are user specified.
!
! To use this you should do the following:
!
! -------------------------------------------------------------------------
! HYCOM  Model preparation:
!
! 1) Make the hycom model dump daily averages. This is done by defining the
! cpp flag "DAILY_AVERAGE" in  MODEL.CPP for Hycom. Daily averages are then
! dumped in unformatted files named RRRDAILY_yyyy_ddd.uf.
!
! -------------------------------------------------------------------------
!
! 2) specify daily average files you want to create netcdf files for. This is
! done in the file "infiles.daily". One daily average file per line.
!
! 3) Specify a regular grid in longitude latitude coordinates in the file named 
! "regugrid.in". Here you ! specify the eastern and western longitudes, along with 
! longitude grid spacing, then the southern and northern latitudes, along with 
! latitude grid spacing. Each value on a separate line in this infile.
!
! 4) Specify depth levels to which you want to interpolate in the file
! "depthlevels.in". Here you specify the number of depth levels on the first
! line, then the actual depthlevels on the following lines. Depth levels must be
! increasing (down is positive direction).
!
! 5) Several other datafiles will be needed as well. These are files used when
! HYCOM runs, so if you run "daily2regunc" in the same catalogue as you run
! hycom, everything should be ok
!
! 6) Run "daily2regunc" to generate netcdf files. For each daily average file, a
! netcdf file will be created. The file name is the same as for the daily
! average files, except that ".uf" is replaced with ".nc". To browse the data
! file you can use the program "ncview" installed on tre.
!
!
!
!NB: Sample infiles are given in the catalogue SAMPLE_INFILES
!
!NB2: To use this program you need








! If you make changes to this program, update the version number!!!


program p_daily2pak
   use mod_common
   use mod_year_info
   use m_year_day
   use m_get_daily_infiles
   use m_datetojulian
   use m_get_grid
   use m_prepak26_hack
   use m_pakmsk
   use m_mkfldh
   use mod_xc, only: mod_xc_idm => idm, mod_xc_jdm => jdm
   use mod_za , only: zaiost
   use mod_read_dailyab
   implicit none

   integer*4, external :: iargc

   ! Version info of this program
   character(len=*), parameter :: cver = 'V0.9'

   ! Input file with days
   character(len=80), dimension(:), pointer :: daily_files
   integer  :: maxfiles,numfiles

   type(year_info)    :: rt,rti,rtd,rtb       ! Time info type
   integer            :: year,day,hour,find,k,ios,imem,nrmem,ifile, &
                         counter, ilon, jlat, refyear,totmem
   real               :: tinfo,undef
   logical            :: ex,readok,l_process_mean_ssh,lok

   ! 2D vars -- on model grid
   real   , dimension(:,:), allocatable :: &
      varout2d,mostrf, ssh, strmf, mld1, mld2, tmp,  emnp, &
      ficem,hicem,hsnwm, uice,vice, surflx,sswflx,salflx,meanssh, &
      ubavg,vbavg,taux,tauy,wice,util1,work2d
   integer, dimension(:,:), allocatable :: varout2dint

   ! 3D vars -- on model grid
   real, dimension(:,:,:), allocatable :: &
      p,  temp, saln, utot, vtot, varout3d, dp

   ! 3D vars -- on regular grid -- on specified z levels
   real, dimension(:,:,:), allocatable :: fldlevel

   ! 3D vars -- on regular grid -- on model depth levels
   real, dimension(:,:,:), allocatable ::  fldint,depthint,depthintm


   character(len=8)   :: dateinfo
   character(len=3)   :: rungen     ! Run ID
   character(len=2),allocatable   :: util(:)
   character(len=80)  :: pakfile    ! Name of new netcdf file
   character(len=80)  :: fbase
   character(len=80)  :: fieldhead
   character(len=40)  :: tmpchar,vartitle,fname
   character(len=160) :: tmpcharlong


   integer :: diffday
   integer :: fname_iday,fname_iyear,fname_dday,fname_dyear
   integer :: est_fsize,finduf,findab
   integer :: i,j, irec
   integer :: n2d, n3d, maxfld, length, idat(3), nfields
   integer , parameter :: maxflddim=1000
   character(len=5) varlist(maxflddim)
   real, allocatable :: sigma(:)

   real :: pi, thref, radian


   include 'stmt_funcs.h'


   ! This sends hycom routines output to "fort.66"
   lp=66

   pi = acos(0.)*2 
   radian= 180./pi
   thref=1e-3


   


   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization stuff. 3 infiles are required:
! infiles.daily   -- Which daily average files to process
! depths.daily    -- Depth levels used 
! regugrid.daily  -- Specification of regular grid


   call get_daily_infiles(daily_files,numfiles,maxfiles)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization stuff for model  -- allocation galore
   print *

   ! Get bathymetry - lon/lat fields and sets up bigrid
   call get_grid 

   ! Get file size from a c-routine call. the string sent is null-terminated
   ! in order for the c-routine to know its end.
   tmpcharlong=daily_files(1);
   finduf=index(daily_files(1),'.uf')
   findab=max(index(daily_files(1),'.a'),index(daily_files(1),'.b'))
   !print *,findab,finduf
   if (finduf>0) then ! Assume all files are ".uf" files

      print *,'Routine does not support .uf-files'
      stop '(daily2pak)'

   elseif (findab>0) then ! Assume all files are ".[ab]" files
      findab=max(index(daily_files(1),'.a'),index(daily_files(1),'.b'))

      !call get_new_filetype_kdm(daily_files(ifile),idm,jdm,nz)
      fbase=daily_files(1)(1:findab-1)
      call daily_average_read_header(trim(fbase),rti,rtd,nrmem,idm,jdm,kdm)
      nz=kdm
      kk=kdm
      mod_xc_idm=idm
      mod_xc_jdm=jdm
      do ifile=1,numfiles
         findab=max(index(daily_files(ifile),'.b'),index(daily_files(ifile),'.a'))
         finduf=index(daily_files(ifile),'.uf')
         if(finduf>0) then
            print *,'NB! Do NOT mix old (.uf) and new (.[ab]) files'
            stop
         end if
      end do
   else
      print *,'File contains no daily files of old or new type...'
      stop '(daily2pak)'
   end if

   print '(a,i3)', ' Number of layers: nz = ',nz


   call zaiost

   ! mod_common
   allocate(p    (idm,jdm,kdm+1),temp(idm,jdm,kdm))
   allocate(saln (idm,jdm,kdm),taux (idm,jdm), dp(idm,jdm,kdm))
   allocate(tauy (idm,jdm),ubavg(idm,jdm),vbavg(idm,jdm))
   allocate(sswflx(idm,jdm), salflx(idm,jdm), meanssh(idm,jdm),surflx(idm,jdm))

   ! mod_common_ice
   allocate(ficem (idm,jdm),hicem (idm,jdm),vice (idm,jdm), wice(idm,jdm))
   allocate(uice  (idm,jdm))


   ! Variables local to this routine
   allocate(utot(idm,jdm,kdm), vtot(idm,jdm,kdm))
   allocate(ssh(idm,jdm), strmf(idm,jdm), mld1(idm,jdm), mld2(idm,jdm))
   allocate(tmp(idm,jdm), emnp(idm,jdm), work2d(idm,jdm))
   !print *,'got grid vars'






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! More Initialization stuff. Allocate grids and initialize
! Splines

   allocate(sigma(kdm))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Phew... Ready to step through files to process

   ! Were ready...
   print *
   do ifile=1,numfiles
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read raw data


      write (6,'(a)',advance='yes') 'Reading from '//trim(daily_files(ifile))

      ! Check for type ...
      finduf=index(daily_files(ifile),'.uf')
      findab=max(index(daily_files(ifile),'.a'),index(daily_files(ifile),'.b'))
      if (finduf>0) then
         print *,'Routine does not support .uf-files'
         stop '(daily2pak)'
      elseif (findab>0) then
         print *,'ab-files in implementation'


         ! Read different fields - returns undef if not defined
         fbase=daily_files(ifile)(1:findab-1)
         call daily_average_read_header(trim(fbase),rti,rtd,nrmem,idm,jdm,kdm)

         call read_field2d(trim(fbase),'ubavg   ',ubavg ,idm,jdm,0,undef)
         where(ubavg/=undef) ubavg=ubavg/nrmem

         call read_field2d(trim(fbase),'vbavg   ',vbavg ,idm,jdm,0,undef)
         where(vbavg/=undef) vbavg=vbavg/nrmem

         call read_field2d(trim(fbase),'ssh     ',ssh   ,idm,jdm,0,undef)
         where(ssh/=undef) ssh=ssh/nrmem

         call read_field2d(trim(fbase),'surflx  ',surflx,idm,jdm,0,undef)
         where(surflx/=undef) surflx=surflx/nrmem

         call read_field2d(trim(fbase),'swflx   ',sswflx,idm,jdm,0,undef)
         where(sswflx/=undef) sswflx=sswflx/nrmem

         call read_field2d(trim(fbase),'salflx  ',salflx,idm,jdm,0,undef)
         where(salflx/=undef) salflx=salflx/nrmem

         call read_field2d(trim(fbase),'emnp    ',emnp  ,idm,jdm,0,undef)
         where(emnp/=undef) emnp=emnp/nrmem

         call read_field2d(trim(fbase),'taux    ',taux  ,idm,jdm,0,undef)
         where(taux/=undef) taux=taux/nrmem

         call read_field2d(trim(fbase),'tauy    ',tauy  ,idm,jdm,0,undef)
         where(tauy/=undef) tauy=tauy/nrmem

         call read_field2d(trim(fbase),'fice    ',ficem ,idm,jdm,0,undef)
         where(ficem/=undef) ficem=ficem/nrmem

         call read_field2d(trim(fbase),'hice    ',hicem ,idm,jdm,0,undef)
         where(hicem/=undef) hicem=hicem/nrmem

         call read_field2d(trim(fbase),'uice    ',uice  ,idm,jdm,0,undef)
         where(uice/=undef) uice=uice/nrmem

         call read_field2d(trim(fbase),'vice    ',vice  ,idm,jdm,0,undef)
         where(vice/=undef) vice=vice/nrmem
         !print *,maxval(vice),minval(vice)

         call read_field2d(trim(fbase),'wice    ',wice  ,idm,jdm,0,undef)
         where(wice/=undef) wice=sqrt(wice/nrmem)

         p=0.
         call read_field3d(trim(fbase),'utot    ',utot  ,idm,jdm,kdm,undef)
         where(utot/=undef) utot=utot/nrmem

         call read_field3d(trim(fbase),'vtot    ',vtot  ,idm,jdm,kdm,undef)
         where(vtot/=undef) vtot=vtot/nrmem

         call read_field3d(trim(fbase),'saln    ',saln  ,idm,jdm,kdm,undef)
         where(saln/=undef) saln=saln/nrmem

         call read_field3d(trim(fbase),'temp    ',temp  ,idm,jdm,kdm,undef)
         where(temp/=undef) temp=temp/nrmem

         call read_field3d(trim(fbase),'pres    ',p(:,:,2:kdm+1) ,idm,jdm,kdm,undef)
         where(p/=undef) p=p/nrmem

      else
         print *,'File contains no daily files of old or new type...'
         stop '(daily2pak)'


      end if



      ! We now have the averaged estimate
      write (6,'(a,i4,a)') 'Fields based on ',nrmem,' members'









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Start processing the raw data


      ! Extend to old/new timesteps here!!
      !mld1=undef; mld2=undef
      !print '(a)','Calculating mld'
      !call mixlayer_depths(temp,saln,p(:,:,1:kdm+1)/onem,mld1,mld2,idm,jdm,kdm)

      !p to dp
      do k=1,kdm
         dp(:,:,k)=p(:,:,k+1)-p(:,:,k)
      end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Start Punching out fields
      rungen=daily_files(ifile)(1:3)
      pakfile=rungen//'DAILY_b'// &
            rti%cyy//rti%cmm//rti%cdm//'_f'//     &
            rtd%cyy//rtd%cmm//rtd%cdm

      ! Prepare to dump pak-file
      n2d=0
      n3d=0
      nfields=0
      varlist(nfields+1:nfields+11)=(/'SSH  ','UBAVG','VBAVG','SURFL', &
         'SWFLX','SALFL','EMNP ','UTAU ','VTAU ','HICEM','FICEM'/)
      n2d=n2d+11

      nfields=n2d
      varlist(nfields+1:nfields+6)=(/'DP   ','TEM  ','SAL  ','TH3D ','UT   ','VT   '/)
      n3d=n3d+6
      maxfld=n2d+n3d
      nfields=maxfld

      write(6,'(a,a)')'I prepak, outname: ',pakfile
      call prepak26_hack(rungen,pakfile,n2d,n3d,varlist,maxfld,idat,rt,ii,jj,kk,sigma)

      irec=1
      allocate(util(ii*jj+14))
      allocate(util1(ii,jj))
      inquire(iolength=j)fieldhead,util(1:ii*jj+14)
      open(10,file=pakfile,status='replace',form='unformatted',access='direct',recl=j)

! --- ---------------------------------------------------------
! --- Punch out fields
! --- ---------------------------------------------------------

      

      call pakmsk(ip,real(ssh),util1,ii,jj,util,length)
      call mkfldh('SSH',1,1,1,ii,ii,jj,length,fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1
     
      call pakmsk(ip,real(ubavg),util1,ii,jj,util,length)
      call mkfldh('UBAVG',1,1,1,ii,ii,jj,length,fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1

      call pakmsk(ip,real(vbavg),util1,ii,jj,util,length)
      call mkfldh('VBAVG',1,1,1,ii,ii,jj,length,fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1

      call pakmsk(ip,surflx,util1,ii,jj,util,length)
      call mkfldh('SURFL',1,1,1,ii,ii,jj,length,fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1

      call pakmsk(ip,sswflx,util1,ii,jj,util,length)
      call mkfldh('SWFLX',1,1,1,ii,ii,jj,length,fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1

      call pakmsk(ip,salflx,util1,ii,jj,util,length)
      call mkfldh('SALFL',1,1,1,ii,ii,jj,length,fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1

      call pakmsk(ip,emnp,util1,ii,jj,util,length)
      call mkfldh('EMNP',1,1,1,ii,ii,jj,length,fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1

      call pakmsk(ip,taux,util1,ii,jj,util,length)
      call mkfldh('UTAU',1,1,1,ii,ii,jj,length,fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1

      call pakmsk(ip,tauy,util1,ii,jj,util,length)
      call mkfldh('VTAU',1,1,1,ii,ii,jj,length,fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1

      call pakmsk(ip,hicem,util1,ii,jj,util,length)
      call mkfldh('HICEM',1,1,1,ii,ii,jj,length,fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1

      call pakmsk(ip,ficem,util1,ii,jj,util,length)
      call mkfldh('FICEM',1,1,1,ii,ii,jj,length,fieldhead)
      write(10,rec=irec)fieldhead,util ; irec=irec+1


      do k=1,kk
         call pakmsk(ip,real(dp(:,:,k)/onem),util1,ii,jj,util,length)
         call mkfldh('DP',1,1,k,ii,ii,jj,length,fieldhead)
         write(10,rec=irec)fieldhead,util ; irec=irec+1

         call pakmsk(ip,temp(:,:,k),util1,ii,jj,util,length)
         call mkfldh('TEM',1,1,k,ii,ii,jj,length,fieldhead)
         write(10,rec=irec)fieldhead,util ; irec=irec+1
         
         call pakmsk(ip,saln(:,:,k),util1,ii,jj,util,length)
         call mkfldh('SAL',1,1,k,ii,ii,jj,length,fieldhead)
         write(10,rec=irec)fieldhead,util ; irec=irec+1

         do j=1,jj
         do i=1,ii
            work2d(i,j) = ip(i,j)*sig(real(saln(i,j,k)),real(temp(i,j,k)))
         end do
         end do
         call pakmsk(ip,work2d,util1,ii,jj,util,length)
         call mkfldh('TH3D',1,1,k,ii,ii,jj,length,fieldhead)
         write(10,rec=irec)fieldhead,util ; irec=irec+1

         call pakmsk(ip,real(utot(:,:,k)),util1,ii,jj,util,length)
         call mkfldh('UT',1,1,k,ii,ii,jj,length,fieldhead)
         write(10,rec=irec)fieldhead,util ; irec=irec+1

         call pakmsk(ip,real(vtot(:,:,k)),util1,ii,jj,util,length)
         call mkfldh('VT',1,1,k,ii,ii,jj,length,fieldhead)
         write(10,rec=irec)fieldhead,util ; irec=irec+1
      end do
      close (10)
      deallocate(util,util1)
   end do


end program p_daily2pak
