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


program p_daily2regunc
   use mod_common
   use mod_year_info
   use mod_toregugrid
   use mod_spline_calc
   use netcdf
   use m_handle_err
   use m_rotate
   use m_strmf_eval
   use m_mixlayer_depths
   use m_pack_short
   use m_year_day
   use m_mosf
   use m_mersea_prepare
   use m_get_daily_infiles
   use m_get_depth_levels
   use m_get_regugrid
   use m_datetojulian
   use m_read_mean_ssh
   use m_get_grid
   use mod_xc, only: mod_xc_idm => idm, mod_xc_jdm => jdm
   use mod_za , only: zaiost
   use mod_read_dailyab
   !use mod_za
   implicit none

   ! These are variable names in NetCDF files
   character(len=*), parameter :: temname    ='temperature'
   character(len=*), parameter :: sshname    ='ssh'
   character(len=*), parameter :: meansshname='meanssh'
   character(len=*), parameter :: salname    ='salinity'
   character(len=*), parameter :: vname      ='v'
   character(len=*), parameter :: uname      ='u'
   character(len=*), parameter :: tauyname   ='tauy'
   character(len=*), parameter :: tauxname   ='taux'
   character(len=*), parameter :: mldname    ='mlp'
   character(len=*), parameter :: qtname     ='qtot'
   character(len=*), parameter :: qswname    ='qsw'
   character(len=*), parameter :: fwname     ='emp'
   character(len=*), parameter :: ubname     ='ub'
   character(len=*), parameter :: vbname     ='vb'
   character(len=*), parameter :: streamfname='bsfd'
   character(len=*), parameter :: mltemname  ='mld'
   character(len=*), parameter :: ficename   ='fice'
   character(len=*), parameter :: hicename   ='hice'

   integer*4, external :: iargc

   ! Version info of this program
   character(len=*), parameter :: cver = 'V0.9'

   ! Input file with days
   character(len=80), dimension(:), pointer :: daily_files
   integer  :: maxfiles,numfiles

   ! Input files for depth levels 
   character(len=*), parameter  :: infile_depths   = 'depthlevels.in'
   integer :: nzlevel
   real, pointer :: zlevels(:)

   ! Input file for regular grid specification
   character(len=*), parameter  :: infile_regugrid = 'regugrid.in'
   real :: lon_sw,lat_sw,dellon,dellat,lon_ne,lat_ne


   type(year_info)    :: rt,rti,rtd,rtb       ! Time info type
   integer            :: year,day,hour,find,k,ios,imem,nrmem,ifile, &
                         counter, ilon, jlat, refyear,totmem
   real               :: tinfo,undef
   logical            :: ex,readok,l_process_mean_ssh,lok

   ! 2D vars -- on model grid
   real   , dimension(:,:), allocatable :: &
      varout2d,mostrf, ssh, strmf, mld1, mld2, tmp,  emnp, &
      ficem,hicem,hsnwm, uice,vice, surflx,sswflx,salflx,meanssh, &
      ubavg,vbavg,taux,tauy,wice
   integer, dimension(:,:), allocatable :: varout2dint

   ! 3D vars -- on model grid
   real, dimension(:,:,:), allocatable :: &
      p,  temp, saln, utot, vtot, & 
      varout3d

   ! 3D vars -- on regular grid -- on specified z levels
   real, dimension(:,:,:), allocatable :: fldlevel

   ! 3D vars -- on regular grid -- on model depth levels
   real, dimension(:,:,:), allocatable ::  fldint,depthint,depthintm


   character(len=8)   :: dateinfo
   character(len=3)   :: rungen     ! Run ID
   character(len=80)  :: ncfil      ! Name of new netcdf file
   character(len=80)  :: ncfil2     ! Name of new netcdf file
   character(len=80)  :: fbase,filename_ice
   character(len=40)  :: tmpchar,vartitle,fname
   character(len=160) :: tmpcharlong


   ! Netcdf variable ids
   integer :: lon_dimension_id,lat_dimension_id,k_dimension_id 
   integer :: lon_dimension_id2, lat_dimension_id2, vars_2d2(2), lon_id2, lat_id2
   integer :: lon_id,lat_id,int_id,temp_id,sal_id,u_id,v_id,pb_id, &
              pbu_id,pbv_id,taux_id,tauy_id,hflx_id,z_id,ssh_id,   &
              strf_id, salflx_id, swflx_id, mld1_id, mld2_id,      &
              emnp_id, meanssh_id, ub_id,vb_id, fice_id, hice_id
   integer :: ncid,ncid2 ! netcdf id of netcdf file
   integer :: vars_2d(2), vars_3d(3)  ! holder for dimension ids

   integer :: NF90_MyFloat, fillval
   logical :: lprocess,lmersea, l_analysis,lemofor,lusefilename
   real    :: scale_factor, add_offset
   !integer :: i,j
   integer :: diffday,iceind
   real    :: basinssh

   integer :: fname_iday,fname_iyear,fname_dday,fname_dyear
   integer :: est_fsize,finduf,findab


   ! This sends hycom routines output to "fort.66"
   lp=66
   l_process_mean_ssh=.false.

   ! Set the fortran 90 datatype here
   NF90_MyFloat=NF90_FLOAT


   ! For new type files

   ! Is this generated for MERSEA ?
   lmersea = .false.
   lemofor = .false.
   lusefilename=.false.
   if (iargc()>=1 .and. iargc() <=2 ) then
      call getarg(1,tmpchar) 
      if ( trim(tmpchar) == "MERSEA" ) then
         print *,'MERSEA file names will be used!'
         lmersea = .true.
      elseif ( trim(tmpchar) == "EMOFOR" ) then
         print *,'EMOFOR file names will be used!'
         lemofor = .true.
      else
         print *,'Unknown naming convention '//trim(tmpchar)
      end if
      if (iargc() == 2) then
         call getarg(2,tmpchar) 
         if ( trim(tmpchar) == "USEFILENAME" ) then
            lusefilename=.true.
         end if
      end if
   elseif (iargc()>2) then
      print *,'Only one or two arguments !!'
      stop '(daily2regunc)'
   end if





   


   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization stuff. 3 infiles are required:
! infiles.daily   -- Which daily average files to process
! depths.daily    -- Depth levels used 
! regugrid.daily  -- Specification of regular grid


   call get_daily_infiles(daily_files,numfiles,maxfiles)

   call get_depth_levels(zlevels,nzlevel)

   call get_regugrid(lon_sw,lon_ne,dellon, &
                     lat_sw,lat_ne,dellat)

   !print *,lon_sw,lon_ne,lat_sw,lat_ne

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

         print *,'uf-files no longer supported by this routine ...'
         stop

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
      stop '(daily2regunc)'
   end if

   print '(a,i3)', ' Number of layers: nz = ',nz


   call zaiost

   ! mod_common
   allocate(p    (idm,jdm,kdm+1),temp(idm,jdm,kdm))
   allocate(saln (idm,jdm,kdm),taux (idm,jdm))
   allocate(tauy (idm,jdm),ubavg(idm,jdm),vbavg(idm,jdm))
   allocate(sswflx(idm,jdm), salflx(idm,jdm), meanssh(idm,jdm),surflx(idm,jdm))

   ! mod_common_ice
   allocate(ficem (idm,jdm),hicem (idm,jdm),vice (idm,jdm), wice(idm,jdm))
   allocate(uice  (idm,jdm))


   ! Variables local to this routine
   allocate(utot(idm,jdm,kdm), vtot(idm,jdm,kdm))
   allocate(ssh(idm,jdm), strmf(idm,jdm), mld1(idm,jdm), mld2(idm,jdm))
   allocate(tmp(idm,jdm), emnp(idm,jdm))
   !print *,'got grid vars'






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! More Initialization stuff. Allocate grids and initialize
! Splines

   ! Initialize for regular grid and for splines
   undef = -1e14
   call regugrid_ini( lon_sw, lon_ne,dellon,lat_sw,lat_ne,dellat,undef)
   call spline_calc_ini(zlevels,nzlevel,undef)
   !print *,nlons,nlats
   allocate(varout2d   (nlons,nlats))
   allocate(varout2dint(nlons,nlats))
   allocate(varout3d   (nlons,nlats,ndeep))


   ! These will hold 3d vars on a regular lon-lat grid
   allocate(fldint   (nlons,nlats,kdm))
   allocate(depthint (nlons,nlats,kdm))
   allocate(depthintm(nlons,nlats,kdm))

   ! These will hold 3d vars on a regular lon-lat grid
   ! after spline- operation
   allocate(fldlevel(nlons,nlats,ndeep))

   ! Merid overturning streamfunction
   allocate(mostrf(nlats,nzlevel))
                  


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
      if (finduf>0) then ! old type

         print *,'uf-files no longer supported by this routine ...'
         stop

      else if (findab>0) then


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

         call read_field3d(trim(fbase),'utot    ',utot  ,idm,jdm,kdm,undef)
         where(utot/=undef) utot=utot/nrmem

         call read_field3d(trim(fbase),'vtot    ',vtot  ,idm,jdm,kdm,undef)
         where(vtot/=undef) vtot=vtot/nrmem

         call read_field3d(trim(fbase),'saln    ',saln  ,idm,jdm,kdm,undef)
         where(saln/=undef) saln=saln/nrmem

         call read_field3d(trim(fbase),'temp    ',temp  ,idm,jdm,kdm,undef)
         where(temp/=undef) temp=temp/nrmem

         p=0.
         call read_field3d(trim(fbase),'pres    ',p(:,:,2:kdm+1) ,idm,jdm,kdm,undef)
         where(p/=undef) p=p/nrmem

      end if



      ! We now have the averaged estimate
      write (6,'(a,i4,a)') 'Fields based on ',nrmem,' members'





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! lprocess tells us wether this file should be processed.
      ! It is modified by MERSEA, processin, but is set to "on" 
      ! as default for all other files
      lprocess=.true.
      l_analysis=.false.
      rungen=daily_files(ifile)(1:3)


      ! For MERSEA
      if (lmersea) then


         ! Use this in case of problems -- overrides rtb and rti read from file
         ! and uses filenames
         if (lusefilename) then
            print *,'Using file name for date calculation'
            tmpchar=daily_files(ifile) 
            read(tmpchar(10:13),'(i4)') fname_iyear
            read(tmpchar(15:17),'(i3)') fname_iday 
            read(tmpchar(19:22),'(i4)') fname_dyear
            read(tmpchar(24:26),'(i3)') fname_dday 
            print *,fname_iday,fname_iyear
            print *,fname_dday,fname_dyear
            call year_day(float(fname_iday),fname_iyear,rti,'ecmwf')
            call year_day(float(fname_dday),fname_dyear,rtd,'ecmwf')
         end if

         call mersea_prepare(rungen,rtb,rti,rtd,lprocess,l_analysis)

         ! For mersea class 1, check hat we are using Levitus levels:
         if (lmersea) then
            if (nzlevel/=12) then
               print *,'NB: You must use NAT depth levels for MERSEA class 1!'
               stop '(daily2average)'
            end if
         end if


         ncfil='topaz_V1_mersea_nat_grid1to8_da_class1_b'// &
               rtb%cyy//rtb%cmm//rtb%cdm//'_f'//           &
               rtd%cyy//rtd%cmm//rtd%cdm//'9999.nc'

         ncfil2='topaz_V1_mersea_nat_grid1to8_da_class1_b'// &
                rti%cyy//rti%cmm//rti%cdm//'_f'//           &
                rtd%cyy//rtd%cmm//rtd%cdm//'9999.nc'
      elseif (lemofor) then

         if (nzlevel/=1) then
            print *,'NB: You must use EMOFOR depth levels for EMOFOR!'
            stop '(daily2regunc)'
         end if

         ! Days start from "0". Add 1 to get real 
         ! dates
         rtb=rti
         write(rtb%cdm,'(i2.2)') rtb%idm+1
         write(rtb%cdd,'(i3.3)') rtb%idd+1
         write(rtd%cdm,'(i2.2)') rtd%idm+1
         write(rtd%cdd,'(i3.3)') rtd%idd+1
         write(rti%cdm,'(i2.2)') rti%idm+1
         write(rti%cdd,'(i3.3)') rti%idd+1 


         ncfil='WafricaDaily_start'// &
               rti%cyy//rti%cmm//rti%cdm//'_dump'//           &
               rtd%cyy//rtd%cmm//rtd%cdm//'.nc'

         lprocess=.true.
      else

         ! Days start from "0". Add 1 to get real 
         ! dates
         rtb=rti
         write(rtb%cdm,'(i2.2)') rtb%idm+1
         write(rtb%cdd,'(i3.3)') rtb%idd+1
         write(rtd%cdm,'(i2.2)') rtd%idm+1
         write(rtd%cdd,'(i3.3)') rtd%idd+1
         write(rti%cdm,'(i2.2)') rti%idm+1
         write(rti%cdd,'(i3.3)') rti%idd+1



         find=index(daily_files(ifile),'DAILY_')
         find=find+5
         ncfil=daily_files(ifile)(1:find)//                           &
               'start'//rti%cyy//rti%cmm//rti%cdm// &
               '_dump'//rtd%cyy//rtd%cmm//rtd%cdm// &
               '.nc'

      end if
      !print *,'rtb (bulletin) : ',rtb
      !print *,'rti (init date): ',rti
      !print *,'rtd (fielddate): ',rtd



      if (lprocess) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Start processing the raw data


      ! Rotate velocities
      print '(a)','Rotating velocities'
      do k=1,kdm
         call rotate(utot(:,:,k),vtot(:,:,k),modlat,modlon,nx,ny,'m2l')
      end do

      ! Rotate surface drag
      call rotate(taux(:,:),tauy(:,:),   modlat,modlon,nx,ny,'m2l')

      ! Barotropic streamfunction
      print '(a)','Calculating streamfunction'
      call strmf_eval(idm,jdm,strmf,ubavg,vbavg,depths,mqlat,mqlon)

      ! Rotate barotropic velocities (after streamf calculation)
      print '(a)','Rotating ubavg '
      call rotate(ubavg(:,:),vbavg(:,:),modlat,modlon,nx,ny,'m2l')

      ! Extend to old/new timesteps here!!
      !mld1=undef; mld2=undef
      print '(a)','Calculating mld'
      call mixlayer_depths(temp,saln,p(:,:,1:kdm+1)/onem,mld1,mld2,idm,jdm,kdm)
      !print '(a)','Done Calculating mld'
      !print *,p(20,25,:)/onem





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- Start NC definition section 

      !Create netcdf file and fill it with data
      if (NF90_CREATE(ncfil,NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening the netcdf file'
         stop '(restart2netcdf)'
      end if


      ! Define dimensions -- idm,jdm,kdm
      ! Following MERCATOR name definitions for netcdf var name
      call handle_err(NF90_DEF_DIM(ncid,'longitude',nlons,lon_dimension_id))
      call handle_err(NF90_DEF_DIM(ncid,'latitude',nlats,lat_dimension_id))
      call handle_err(NF90_DEF_DIM(ncid,'depth'  ,ndeep,k_dimension_id))

      ! Define dimension "holders"
      vars_3d=(/lon_dimension_id,lat_dimension_id,k_dimension_id/)
      vars_2d=(/lon_dimension_id,lat_dimension_id/)

! Define variables -- put attributes

      !Global attributes
      rt=rtd
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'title',        &
         'TOPAZ daily average'))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'institution',        &
         'NERSC, Edvard Griegs Vei 3A, Bergen, Norway'))
      call date_and_time(date=dateinfo)
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Dataset_generation_date', &
         dateinfo(1:8)))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'history', &
         'Created '//dateinfo(1:8)//' by program regudaily2nc '//cver))
      if (lmersea) then
         call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'references', &
            'http://www.mersea.eu.org ; http://topaz.nersc.no'))
      else
         call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'references', &
            'http://topaz.nersc.no'))
      endif
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'field_type', &
         'Daily average fields'))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Conventions', &
                      'COARDS'))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'field_date', &
         rtd%cyy//'-'//rtd%cmm//'-'//rtd%cdm))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'bulletin_date', &
         rtb%cyy//'-'//rtb%cmm//'-'//rtb%cdm))

      if (lmersea.or.lmersea) then
         if (rungen=="FOR") then
            call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'bulletin_type', &
                            'Forecast'))
            call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Forecast_range', &
                            '7 days'))
         elseif (rungen=="RTH") then
            call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'bulletin_type', &
                            'Hindcast'))
         end if
      end if


      

      ! Longitude
      vartitle='longitude'
      call handle_err(NF90_DEF_VAR(ncid,'longitude',NF90_MyFloat,lon_dimension_id,lon_id))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'long_name',trim(vartitle)))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'units','degrees_east'))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'valid_range',(/-180.0,180.0/)))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'missing_value',undef))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'standard_name','longitude'))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'axis','X'))

      ! Latitude
      vartitle= 'latitude'
      call handle_err(NF90_DEF_VAR(ncid,'latitude',NF90_MyFloat,lat_dimension_id,lat_id))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'long_name',trim(vartitle)))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'units','degrees_north'))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'valid_range',(/-90.0, 90.0/)))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'missing_value',undef))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'standard_name','latitude'))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'axis','Y'))

      ! Depth levels
      vartitle= 'depth'
      call handle_err(NF90_DEF_VAR(ncid,'depth',NF90_MyFloat,k_dimension_id,z_id))
      call handle_err(NF90_PUT_ATT(ncid,z_id,'long_name',trim(vartitle)))
      call handle_err(NF90_PUT_ATT(ncid,z_id,'units','m'))
      call handle_err(NF90_PUT_ATT(ncid,z_id,'valid_range',(/0., maxval(zlevels)+1/)))
      call handle_err(NF90_PUT_ATT(ncid,z_id,'missing_value',undef))
      call handle_err(NF90_PUT_ATT(ncid,z_id,'standard_name','depth'))
      call handle_err(NF90_PUT_ATT(ncid,z_id,'positive','down'))
      call handle_err(NF90_PUT_ATT(ncid,z_id,'axis','Z'))

      ! Temperature
      call handle_err(NF90_DEF_VAR(ncid,temname,NF90_SHORT,vars_3d,temp_id))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'long_name','potential temperature'))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'units','degree_Celsius'))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'standard_name', &
                      'sea_water_potential_temperature'))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'scale_factor',1.))

      ! Salinity
      call handle_err(NF90_DEF_VAR(ncid,salname,NF90_SHORT,vars_3d,sal_id))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'long_name',salname))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'units','psu'))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'standard_name', &
                      'sea_water_salinity'))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'scale_factor',1.))

      ! U-velocity
      call handle_err(NF90_DEF_VAR(ncid,uname,NF90_SHORT,vars_3d,u_id))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'long_name','Eastward velocity'))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'units','meter second-1'))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'standard_name', &
                      'eastward_sea_water_velocity'))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'scale_factor',1.))

      ! V-velocity
      call handle_err(NF90_DEF_VAR(ncid,vname,NF90_SHORT,vars_3d,v_id))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'long_name','Northward velocity'))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'units','meter second-1'))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'standard_name', &
                      'northward_sea_water_velocity'))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'scale_factor',1.))

      ! ssh
      call handle_err(NF90_DEF_VAR(ncid,sshname,NF90_SHORT,vars_2d,ssh_id))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'long_name','Sea Surface Height'))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'units','meter'))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'standard_name', &
                      'sea_surface_elevation'))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'_FillValue',0_2)) ! Short

      ! Barotropic Streamfunction
      call handle_err(NF90_DEF_VAR(ncid,streamfname,NF90_SHORT,vars_2d,strf_id))
      call handle_err(NF90_PUT_ATT(ncid,strf_id,'long_name','Barotropic Streamfunction'))
      call handle_err(NF90_PUT_ATT(ncid,strf_id,'units','meter3 second-1'))
      call handle_err(NF90_PUT_ATT(ncid,strf_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,strf_id,'standard_name', &
                      'ocean_barotropic_streamfunction '))
      call handle_err(NF90_PUT_ATT(ncid,strf_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,strf_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,strf_id,'_FillValue',0_2)) ! Short

      ! MLD1
      call handle_err(NF90_DEF_VAR(ncid,mltemname,NF90_SHORT,vars_2d,mld1_id))
      call handle_err(NF90_PUT_ATT(ncid,mld1_id,'long_name','Temperature defined mixlayer'))
      call handle_err(NF90_PUT_ATT(ncid,mld1_id,'units','meter'))
      call handle_err(NF90_PUT_ATT(ncid,mld1_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,mld1_id,'standard_name', &
                      'ocean_mixed_layer_thickness'))
      call handle_err(NF90_PUT_ATT(ncid,mld1_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,mld1_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,mld1_id,'_FillValue',0_2)) ! Short

      ! MLD2
      call handle_err(NF90_DEF_VAR(ncid,mldname,NF90_SHORT,vars_2d,mld2_id))
      call handle_err(NF90_PUT_ATT(ncid,mld2_id,'long_name','Density defined mixlayer'))
      call handle_err(NF90_PUT_ATT(ncid,mld2_id,'units','meter'))
      call handle_err(NF90_PUT_ATT(ncid,mld2_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,mld2_id,'standard_name', &
                      'ocean_mixed_layer_thickness'))
      call handle_err(NF90_PUT_ATT(ncid,mld2_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,mld2_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,mld2_id,'_FillValue',0_2)) ! Short

      ! Atmospheric stress on surface
      call handle_err(NF90_DEF_VAR(ncid,tauxname,NF90_SHORT,vars_2d,taux_id))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'long_name','Atm Stress Eastwards'))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'units','pascal'))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'standard_name', &
                      'surface_downward_eastward_stress'))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'_FillValue',0_2)) ! Short

      ! Atmospheric stress on surface
      call handle_err(NF90_DEF_VAR(ncid,tauyname,NF90_SHORT,vars_2d,tauy_id))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'long_name','Atm Stress Northwards'))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'units','pascal'))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'standard_name', &
                      'surface_downward_northward_stress'))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'_FillValue',0_2)) ! Short

      ! Net heatflux at surface
      call handle_err(NF90_DEF_VAR(ncid,qtname,NF90_SHORT,vars_2d,hflx_id))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'long_name','net surface heat flux into ocean'))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'units','watt meter-2'))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'standard_name',  &
                      'surface_downward_heat_flux_in_air'))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'_FillValue',0_2)) ! Short

      ! Shrtwave heatflux at surface
      call handle_err(NF90_DEF_VAR(ncid,qswname,NF90_SHORT,vars_2d,swflx_id))
      call handle_err(NF90_PUT_ATT(ncid,swflx_id,'long_name','Shortwave heat flux into ocean'))
      call handle_err(NF90_PUT_ATT(ncid,swflx_id,'units','watt meter-2'))
      call handle_err(NF90_PUT_ATT(ncid,swflx_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,swflx_id,'standard_name',  &
                      'net_downward_shortwave_flux_in_air '))
      call handle_err(NF90_PUT_ATT(ncid,swflx_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,swflx_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,swflx_id,'_FillValue',0_2)) ! Short

      ! (evap - precip - rivers)
      call handle_err(NF90_DEF_VAR(ncid,fwname,NF90_SHORT,vars_2d,emnp_id))
      call handle_err(NF90_PUT_ATT(ncid,emnp_id,'long_name','Water Flux'))
      call handle_err(NF90_PUT_ATT(ncid,emnp_id,'units','kilogram meter-2 second-1'))
      call handle_err(NF90_PUT_ATT(ncid,emnp_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,emnp_id,'standard_name',  &
                      'water_flux_into_ocean'))
      call handle_err(NF90_PUT_ATT(ncid,emnp_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,emnp_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,emnp_id,'_FillValue',0_2)) ! Short

      if (.not.lmersea) then
         call handle_err(NF90_DEF_VAR(ncid,ubname,NF90_SHORT,vars_2d,ub_id))
         call handle_err(NF90_PUT_ATT(ncid,ub_id,'long_name','Eastward Barotropic velocity'))
         call handle_err(NF90_PUT_ATT(ncid,ub_id,'units','meter second-1'))
         call handle_err(NF90_PUT_ATT(ncid,ub_id,'missing_value',0_2))
         call handle_err(NF90_PUT_ATT(ncid,ub_id,'standard_name',  &
                         'eastward_barotropic_sea_water_velocity'))
         call handle_err(NF90_PUT_ATT(ncid,ub_id,'add_offset',0.))
         call handle_err(NF90_PUT_ATT(ncid,ub_id,'scale_factor',1.))
         call handle_err(NF90_PUT_ATT(ncid,ub_id,'_FillValue',0_2)) ! Short

         call handle_err(NF90_DEF_VAR(ncid,vbname,NF90_SHORT,vars_2d,vb_id))
         call handle_err(NF90_PUT_ATT(ncid,vb_id,'long_name','Northward Barotropic velocity'))
         call handle_err(NF90_PUT_ATT(ncid,vb_id,'units','meter second-1'))
         call handle_err(NF90_PUT_ATT(ncid,vb_id,'missing_value',0_2))
         call handle_err(NF90_PUT_ATT(ncid,vb_id,'standard_name',  &
                         'northward_barotropic_sea_water_velocity'))
         call handle_err(NF90_PUT_ATT(ncid,vb_id,'add_offset',0.))
         call handle_err(NF90_PUT_ATT(ncid,vb_id,'scale_factor',1.))
         call handle_err(NF90_PUT_ATT(ncid,vb_id,'_FillValue',0_2)) ! Short
      end if


      if (.not.lmersea) then
         call handle_err(NF90_DEF_VAR(ncid,'hice',NF90_SHORT,vars_2d,hice_id))
         call handle_err(NF90_PUT_ATT(ncid,hice_id,'long_name','Ice Thickness'))
         call handle_err(NF90_PUT_ATT(ncid,hice_id,'units','meter'))
         call handle_err(NF90_PUT_ATT(ncid,hice_id,'missing_value',0_2))
         call handle_err(NF90_PUT_ATT(ncid,hice_id,'standard_name',  &
                         'sea_ice_thickness'))
         call handle_err(NF90_PUT_ATT(ncid,hice_id,'add_offset',0.))
         call handle_err(NF90_PUT_ATT(ncid,hice_id,'scale_factor',1.))
         call handle_err(NF90_PUT_ATT(ncid,hice_id,'_FillValue',0_2)) ! Short

         call handle_err(NF90_DEF_VAR(ncid,'fice',NF90_SHORT,vars_2d,fice_id))
         call handle_err(NF90_PUT_ATT(ncid,fice_id,'long_name','Ice Concentration'))
         call handle_err(NF90_PUT_ATT(ncid,fice_id,'units','meter'))
         call handle_err(NF90_PUT_ATT(ncid,fice_id,'missing_value',0_2))
         call handle_err(NF90_PUT_ATT(ncid,fice_id,'standard_name',  &
                         'sea_ice_thickness'))
         call handle_err(NF90_PUT_ATT(ncid,fice_id,'add_offset',0.))
         call handle_err(NF90_PUT_ATT(ncid,fice_id,'scale_factor',1.))
         call handle_err(NF90_PUT_ATT(ncid,fice_id,'_FillValue',0_2)) ! Short
      end if

      call handle_err(NF90_ENDDEF(ncid))
      !print *,'end def'

      ! For processing meanssh -- only needed once
      if (l_process_mean_ssh)  then
         print '(a)','  Generating meanssh file'
         if (NF90_CREATE("topaz_nat_grid1to8_meanssh.nc",NF90_CLOBBER,ncid2) /= NF90_NOERR) then
            print *,'An error occured when opening the netcdf file'
            stop '(restart2netcdf)'
         end if

         !Global attributes
         rt=rtd
         call handle_err(NF90_PUT_ATT(ncid2,NF90_GLOBAL,'title',        &
            'TOPAZ Mean SSH'))
         call handle_err(NF90_PUT_ATT(ncid2,NF90_GLOBAL,'institution',        &
            'NERSC, Edvard Griegs Vei 3A, Bergen, Norway'))
         call date_and_time(date=dateinfo)
         call handle_err(NF90_PUT_ATT(ncid2,NF90_GLOBAL,'Dataset_generation_date', &
            dateinfo(1:8)))
         call handle_err(NF90_PUT_ATT(ncid2,NF90_GLOBAL,'history', &
            'Created '//dateinfo(1:8)//' by program regudaily2nc '//cver))
         if (lmersea) then
            call handle_err(NF90_PUT_ATT(ncid2,NF90_GLOBAL,'references', &
               'http://www.mersea.eu.org ; http://topaz.nersc.no'))
         else
            call handle_err(NF90_PUT_ATT(ncid2,NF90_GLOBAL,'references', &
               'http://topaz.nersc.no'))
         endif
         call handle_err(NF90_PUT_ATT(ncid2,NF90_GLOBAL,'field_type', &
            'Constant in time'))
         call handle_err(NF90_PUT_ATT(ncid2,NF90_GLOBAL,'Conventions', &
                         'COARDS'))

         ! Define dimensions -- idm,jdm,kdm
         ! Following MERCATOR name definitions for netcdf var name
         call handle_err(NF90_DEF_DIM(ncid2,'longitude',nlons,lon_dimension_id2))
         call handle_err(NF90_DEF_DIM(ncid2,'latitude',nlats,lat_dimension_id2))

         vars_2d2=(/lon_dimension_id2,lat_dimension_id2/)

         ! Longitude
         vartitle='longitude'
         call handle_err(NF90_DEF_VAR(ncid2,'longitude',NF90_MyFloat,lon_dimension_id2,lon_id2))
         call handle_err(NF90_PUT_ATT(ncid2,lon_id2,'long_name',trim(vartitle)))
         call handle_err(NF90_PUT_ATT(ncid2,lon_id2,'units','degrees_east'))
         call handle_err(NF90_PUT_ATT(ncid2,lon_id2,'valid_range',(/-180.0,180.0/)))
         call handle_err(NF90_PUT_ATT(ncid2,lon_id2,'missing_value',undef))
         call handle_err(NF90_PUT_ATT(ncid2,lon_id2,'standard_name','longitude'))
         call handle_err(NF90_PUT_ATT(ncid2,lon_id2,'axis','X'))

         ! Latitude
         vartitle= 'latitude'
         call handle_err(NF90_DEF_VAR(ncid2,'latitude',NF90_MyFloat,lat_dimension_id2,lat_id2))
         call handle_err(NF90_PUT_ATT(ncid2,lat_id2,'long_name',trim(vartitle)))
         call handle_err(NF90_PUT_ATT(ncid2,lat_id2,'units','degrees_north'))
         call handle_err(NF90_PUT_ATT(ncid2,lat_id2,'valid_range',(/-90.0, 90.0/)))
         call handle_err(NF90_PUT_ATT(ncid2,lat_id2,'missing_value',undef))
         call handle_err(NF90_PUT_ATT(ncid2,lat_id2,'standard_name','latitude'))
         call handle_err(NF90_PUT_ATT(ncid2,lat_id2,'axis','Y'))


         ! Mean ssh
         call handle_err(NF90_DEF_VAR(ncid2,meansshname,NF90_MyFloat,vars_2d2,meanssh_id))
         call handle_err(NF90_PUT_ATT(ncid2,meanssh_id,'long_name','Mean Sea Surface Height'))
         call handle_err(NF90_PUT_ATT(ncid2,meanssh_id,'units','meter'))
         call handle_err(NF90_PUT_ATT(ncid2,meanssh_id,'missing_value',undef))
         call handle_err(NF90_PUT_ATT(ncid2,meanssh_id,'standard_name', &
                         'mean_sea_surface_elevation'))
         call handle_err(NF90_PUT_ATT(ncid2,meanssh_id,'_FillValue',real(undef,4))) ! Short
         call handle_err(NF90_ENDDEF(ncid2))

         call handle_err(NF90_PUT_VAR(ncid2,lon_id,lons))
         call handle_err(NF90_PUT_VAR(ncid2,lat_id,lats))

         call read_mean_ssh(meanssh)
         call to_regugrid(meanssh,varout2d)

         ! Remove som troublesome points
         where (varout2d>10) varout2d=undef;

         basinssh=sum(varout2d,mask=(abs(varout2d-undef)>1e-4))
         basinssh=basinssh/max(1,count(abs(varout2d-undef)>1e-4))
         where ((abs(varout2d-undef)>1e-4)) varout2d=varout2d-basinssh

         call handle_err(NF90_PUT_VAR(ncid2,meanssh_id,varout2d))
         call handle_err(NF90_CLOSE(ncid2))
         print '(a)',' Done Generating meanssh file'

         ! Only needed once. If the flag is on (which it is) unset it
         l_process_mean_ssh=.false.
      end if


! --- End definition section 





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- Start variable put section 

      call handle_err(NF90_PUT_VAR(ncid,lon_id,lons))
      call handle_err(NF90_PUT_VAR(ncid,lat_id,lats))
      call handle_err(NF90_PUT_VAR(ncid,z_id,deeps))

      ! Regu grid and spline calculations for 3D fields
      call to_regugrid(p(:,:,2:kdm+1)/onem,depthint)
      call to_regugrid((p(:,:,1:kdm)+p(:,:,2:kdm+1))/(2*onem),depthintm)

      write(6,'(a)',advance='no') '3D temp interpolation:'
      call to_regugrid(temp(:,:,1:kdm),fldint)
      call spline_calc(fldint,depthint,nlons,nlats,ongrid,fldlevel,nzlevel)
      call pack_short(fldlevel,nlons,nlats,ndeep,undef,ncid,temp_id,temname)
      print*

      write(6,'(a)',advance='no') '3D salt interpolation:'
      call to_regugrid(saln(:,:,1:kdm),fldint)
      call spline_calc(fldint,depthint,nlons,nlats,ongrid,fldlevel,nzlevel)
      call pack_short(fldlevel,nlons,nlats,ndeep,undef,ncid,sal_id ,salname)
      print*

      write(6,'(a)',advance='no') '3D u    interpolation:'
      call to_regugrid(utot(:,:,1:kdm),fldint)
      call spline_calc(fldint,depthint,nlons,nlats,ongrid,fldlevel,nzlevel)
      call pack_short(fldlevel,nlons,nlats,ndeep,undef,ncid,u_id   ,uname)
      print*

      write(6,'(a)',advance='no') '3D v    interpolation:'
      call to_regugrid(vtot(:,:,1:kdm),fldint)
      call spline_calc(fldint,depthint,nlons,nlats,ongrid,fldlevel,nzlevel)
      call pack_short(fldlevel,nlons,nlats,ndeep,undef,ncid,v_id   ,vname)
      !call spline_calc(uint,vint,dint,tint,sint,nlons,nlats,ongrid, &
      !                 ulevel,vlevel,dlevel,tlevel,slevel,nzlevel)
      print*




      ! Experiment with add_offset and scale_factor
      varout2d=0.
      call to_regugrid(ssh,varout2d)
      basinssh=sum(ssh,mask=(abs(ssh-undef)>1e-4))
      basinssh=basinssh/max(1,count(abs(ssh-undef)>1e-4))
      where ((abs(ssh-undef)>1e-4)) ssh=ssh-basinssh
      call pack_short(varout2d,nlons,nlats,1,undef,ncid,ssh_id,sshname)

      !print *,'max/min strf:',maxval(strmf),minval(strmf)
      call to_regugrid(strmf(:,:),varout2d)
      call pack_short(varout2d,nlons,nlats,1,undef,ncid,strf_id,streamfname)

      call to_regugrid(mld1(:,:),varout2d)
      call pack_short(varout2d,nlons,nlats,1,undef,ncid,mld1_id,mltemname)

      call to_regugrid(mld2(:,:),varout2d)
      call pack_short(varout2d,nlons,nlats,1,undef,ncid,mld2_id,mldname)

      call to_regugrid(taux(:,:),varout2d)
      call pack_short(varout2d,nlons,nlats,1,undef,ncid,taux_id,tauxname)

      call to_regugrid(tauy(:,:),varout2d)
      call pack_short(varout2d,nlons,nlats,1,undef,ncid,tauy_id,tauyname)

      call to_regugrid(surflx(:,:),varout2d)
      call pack_short(varout2d,nlons,nlats,1,undef,ncid,hflx_id,qtname)

      call to_regugrid(sswflx(:,:),varout2d)
      call pack_short(varout2d,nlons,nlats,1,undef,ncid,swflx_id,qswname)

      call to_regugrid(-emnp(:,:)*1000.,varout2d)
      call pack_short(varout2d,nlons,nlats,1,undef,ncid,emnp_id,fwname)

      if (.not.lmersea) then
         call to_regugrid(ubavg,varout2d)
         call pack_short(varout2d,nlons,nlats,1,undef,ncid,ub_id,ubname)

         call to_regugrid(vbavg,varout2d)
         call pack_short(varout2d,nlons,nlats,1,undef,ncid,vb_id,vbname)

         !print *,'max/min hice:',maxval(hicem),minval(hicem)
         call to_regugrid(hicem(:,:),varout2d)
         call pack_short(varout2d,nlons,nlats,1,undef,ncid,hice_id,hicename)

         call to_regugrid(ficem(:,:),varout2d)
         call pack_short(varout2d,nlons,nlats,1,undef,ncid,fice_id,ficename)
      end if

      !print *,'2dvars'
! --- End variable put section 

      !Close up
      call handle_err(NF90_CLOSE(ncid))
      print '(a)','Generated '//trim(ncfil)

      ! Special case
      if (l_analysis.and.lmersea) then
         print *,'NB: Special case -- processing "Analysis"'
         print '(a)','Generating '//trim(ncfil2)
         tmpcharlong="cp "//trim(ncfil)//" "//trim(ncfil2)
         print *,trim(tmpcharlong)
         call system( trim(tmpcharlong) )
      end if

      end if !lprocess
      print *

   enddo ! ifile-loop

      

end program p_daily2regunc
