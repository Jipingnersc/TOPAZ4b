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
! dumped in unformatted files named RRRDAILY_yyyy_ddd_YYYY_DDDD.[ab]
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
! 6) New - option to dump only specified fields. Create a file
! "fields.daily2regunc" where you specify which fields to dump. The name of the
! fields to specify is the same as those in the ".b" header files.
!
! 7) Run "daily2regunc" to generate netcdf files. For each daily average file, a
! netcdf file will be created. The file name is the same as for the daily
! average files, except that ".uf" is replaced with ".nc". To browse the data
! file you can use the program "ncview" installed on tre.
!
!
!
!NB: Sample infiles are given in the catalogue SAMPLE_INFILES
!
! To add new fields, create a new var_info entry in the file "mod_netcdf_pars",
! and an entry for the variable in the routine below. See how Temperature is
! processed for an example on how 3D vars are processed. Look at "hice" for how
! 2D vars are processed (same, but without vertical interpolation.








! If you make changes to this program, update the version number!!!


program p_daily2regunc
   use mod_common
   use mod_year_info
   use mod_toregugrid
   use mod_spline_calc
   use mod_netcdf_pars
   use netcdf
   use m_handle_err
   use m_rotate
   use m_strmf_eval
   use m_mixlayer_depths
   use m_pack_short
   use m_year_day
   use m_mersea_prepare
   use m_get_daily_infiles
   use m_get_depth_levels
   use m_get_regugrid
   use m_datetojulian
   use m_get_grid
   use mod_xc, only: mod_xc_idm => idm, mod_xc_jdm => jdm
   use mod_za , only: zaiost, zaioempty
   use mod_read_dailyab
   use m_def_netcdf_var
   implicit none

   !integer*4, external :: iargc

   ! Version info of this program
   character(len=*), parameter :: cver = 'V1.1'

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
   integer            :: year,day,hour,find,k,ios,nrmem,ifile, &
                         counter, refyear
   real               :: undef
   logical            :: ex

   real, allocatable, dimension(:,:,:) :: hy3d, regu3d, regusp3d, depthint, &
      hy3d2,hy3d3
   real, allocatable, dimension(:,:)   :: hy2d, hy2d2, regu2d, strmf, &
      mld1, mld2, hy2d3,hy2d4

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
   integer :: var_id
   integer :: ncid,ncid2 ! netcdf id of netcdf file
   integer :: vars_2d(2), vars_3d(3)  ! holder for dimension ids

   integer :: NF90_MyFloat
   logical :: lmersea, lemofor,lusefilename
   real    :: basinssh

   integer :: fname_iday,fname_iyear,fname_dday,fname_dyear
   integer :: finduf,findab
   integer :: i, i2, i3


   ! This sends hycom routines output to "fort.66"
   lp=66

   ! Set the fortran 90 datatype here
   NF90_MyFloat=NF90_FLOAT


   ! For new type files

!   ! Is this generated for MERSEA ?
!   lmersea = .false.
!   lemofor = .false.
!   lusefilename=.false.
!   if (iargc()>=1 .and. iargc() <=2 ) then
!      call getarg(1,tmpchar) 
!      if ( trim(tmpchar) == "MERSEA" ) then
!         print *,'MERSEA file names will be used!'
!         lmersea = .true.
!      elseif ( trim(tmpchar) == "EMOFOR" ) then
!         print *,'EMOFOR file names will be used!'
!         lemofor = .true.
!      else
!         print *,'Unknown naming convention '//trim(tmpchar)
!      end if
!      if (iargc() == 2) then
!         call getarg(2,tmpchar) 
!         if ( trim(tmpchar) == "USEFILENAME" ) then
!            lusefilename=.true.
!         end if
!      end if
!   elseif (iargc()>2) then
!      print *,'Only one or two arguments !!'
!      stop '(daily2regunc)'
!   end if





   


   
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

!   print '(a,i3)', ' Number of layers: nz = ',nz


   call zaiost

   ! Variables local to this routine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! More Initialization stuff. Allocate grids and initialize
! Splines

   ! Initialize for regular grid and for splines
   undef = -1e14
   call regugrid_ini( lon_sw, lon_ne,dellon,lat_sw,lat_ne,dellat,undef)
   call spline_calc_ini(zlevels,nzlevel,undef)
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize netcdf parameter info
   call init_var_info()

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

      else if (.not. findab>0) then

         print *,'No .ab files'
         stop

      end if


      ! Read different fields - returns undef if not defined
      fbase=daily_files(ifile)(1:findab-1)
      call daily_average_read_header(trim(fbase),rti,rtd,nrmem,idm,jdm,kdm)

      !mod_xc_idm=idm
      !mod_xc_jdm=jdm
      !call zaiost

      ! We now have the averaged estimate
      write (6,'(a,i4,a)') 'Fields based on ',nrmem,' members'


      ! Start allocating fileds now..
      allocate(hy3d     (idm,jdm,kdm))        ! Holds 3D vars on original grid
      allocate(hy3d2    (idm,jdm,kdm))        ! Holds 3D vars on original grid
      allocate(hy3d3    (idm,jdm,kdm+1))        ! Holds 3D vars on original grid
      allocate(regu3d   (nlons,nlats,kdm))    ! 3D vars on regular grid, original vertical grid
      allocate(regusp3d (nlons,nlats,ndeep))  ! 3D vars on regular grid, z-level vertical grid
      allocate(strmf    (idm,jdm))            ! 2D stream func
      allocate(mld1     (idm,jdm))            ! MLD1
      allocate(mld2     (idm,jdm))            ! MLD2
      allocate(hy2d     (idm,jdm))            ! Holds 2D vars on original grid
      allocate(hy2d2    (idm,jdm))            ! Holds 2D vars on original grid
      allocate(hy2d3    (idm,jdm))            ! Holds 2D vars on original grid
      allocate(hy2d4    (idm,jdm))            ! Holds 2D vars on original grid
      allocate(regu2d   (nlons,nlats))        ! 2D vars on regular gri
      allocate(depthint (nlons,nlats,kdm))    ! Keeps bottom interface on regu grid





!         call read_field2d(trim(fbase),'salflx  ',salflx,idm,jdm,0,undef)
!         where(salflx/=undef) salflx=salflx/nrmem












!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rungen=daily_files(ifile)(1:3)


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

      ! NB: For MERSEA stuff - look at old versions - not impl yet



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- Start NC definition section 

      !Create netcdf file and fill it with data
      if (NF90_CREATE(ncfil,NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening the netcdf file'
         stop '(restart2netcdf)'
      end if



      ! Regu grid and spline calculations for 3D fields
      !p=0.
      call read_field3d(trim(fbase),'pres    ',hy3d ,idm,jdm,kdm,undef)
      where(hy3d/=undef) hy3d=hy3d/nrmem
      call to_regugrid(hy3d/onem,depthint)
      !call to_regugrid((p(:,:,1:kdm)+p(:,:,2:kdm+1))/(2*onem),depthintm)

      ! used later on
      hy3d3=0.
      hy3d3(:,:,2:kdm+1)=hy3d




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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- Start defining vars - put them as we move along
! Longitude
      vartitle='longitude'
      call handle_err(NF90_DEF_VAR(ncid,'longitude',NF90_MyFloat,lon_dimension_id,var_id))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name',trim(vartitle)))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'units','degrees_east'))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'valid_range',(/-180.0,180.0/)))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','longitude'))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'axis','X'))
      call handle_err(nf90_enddef(ncid))
      call handle_err(NF90_PUT_VAR(ncid,var_id,lons))

! Latitude
      call handle_err(nf90_redef(ncid))
      vartitle= 'latitude'
      call handle_err(NF90_DEF_VAR(ncid,'latitude',NF90_MyFloat,lat_dimension_id,var_id))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name',trim(vartitle)))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'units','degrees_north'))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'valid_range',(/-90.0, 90.0/)))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','latitude'))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'axis','Y'))
      call handle_err(nf90_enddef(ncid))
      call handle_err(NF90_PUT_VAR(ncid,var_id,lats))

      

! Depth levels
      call handle_err(nf90_redef(ncid))
      vartitle= 'depth'
      call handle_err(NF90_DEF_VAR(ncid,'depth',NF90_MyFloat,k_dimension_id,var_id))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name',trim(vartitle)))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'units','m'))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'valid_range',(/0., maxval(zlevels)+1/)))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','depth'))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'positive','down'))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'axis','Z'))
      call handle_err(nf90_enddef(ncid))
      call handle_err(NF90_PUT_VAR(ncid,var_id  ,deeps))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Temperature
      ! Gets "database" info for what we want to process
      i=var_info_getindex('temp') 

      if (i/=-1.and.var_info(i)%process) then  ! I we found it ...

         ! Define netcdf variable - using what we found in the database
         call def_netcdf_var(ncid,var_id,vars_3d,3,var_info(i))   

         !  Read field from data file
         call read_field3d(trim(fbase),'temp    ',hy3d  ,idm,jdm,kdm,undef) 
         where(hy3d/=undef) hy3d=hy3d/nrmem


         ! Interpolate to regular grid (horizontal)
         call to_regugrid(hy3d(:,:,1:kdm),regu3d)

         ! Interpolate to fixed vertical levels 
         write(6,'(a)',advance='no') '3D temp interpolation:'
         call spline_calc(regu3d,depthint,nlons,nlats,ongrid,regusp3d,nzlevel)
         print *

         ! Pack data as short type, and put it into netcdf file
         call pack_short(regusp3d,nlons,nlats,ndeep,undef,ncid,var_id,var_info(i))
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Salinity 
      i=var_info_getindex('saln') 
      if (i/=-1.and.var_info(i)%process) then
         call def_netcdf_var(ncid,var_id,vars_3d,3,var_info(i))  
         call read_field3d(trim(fbase),'saln    ',hy3d2 ,idm,jdm,kdm,undef)
         where(hy3d2/=undef) hy3d2=hy3d2/nrmem
         call to_regugrid(hy3d2(:,:,1:kdm),regu3d)
         write(6,'(a)',advance='no') '3D salt interpolation:'
         call spline_calc(regu3d,depthint,nlons,nlats,ongrid,regusp3d,nzlevel)
         print *
         call pack_short(regusp3d,nlons,nlats,ndeep,undef,ncid,var_id,var_info(i))
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate mlds from temp/saln
      ! hy3d =temp
      ! hy3d2=saln
      ! hy3d3=pres
      i=var_info_getindex('temp')
      i2=var_info_getindex('saln')
      if (i /=-1.and.var_info(i )%process .and.  i2/=-1.and.var_info(i2)%process ) then
         print '(a)','Calculating mld'
         call mixlayer_depths(hy3d,hy3d2,hy3d3(:,:,1:kdm+1)/onem,mld1,mld2,idm,jdm,kdm)
      else
         i=var_info_getindex('mld1') 
         if (i/=-1) var_info(i)%process=.false.
         i=var_info_getindex('mld2')
         if (i/=-1) var_info(i)%process=.false.
      end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vlocities - must be rotated first
      i =var_info_getindex('utot')
      i2=var_info_getindex('vtot')
      if (i/=-1.and.var_info(i)%process .and.  i2/=-1.and.var_info(i2)%process) then
         call read_field3d(trim(fbase),'utot    ',hy3d  ,idm,jdm,kdm,undef)
         call read_field3d(trim(fbase),'vtot    ',hy3d2 ,idm,jdm,kdm,undef)
         where(hy3d /=undef) hy3d =hy3d /nrmem
         where(hy3d2/=undef) hy3d2=hy3d2/nrmem

         do k=1,kdm
            call rotate(hy3d(:,:,k),hy3d2(:,:,k),   modlat,modlon,nx,ny,'m2l')
         end do

! U-velocity
         call def_netcdf_var(ncid,var_id,vars_3d,3,var_info(i))
         call to_regugrid(hy3d(:,:,1:kdm),regu3d)
         write(6,'(a)',advance='no') '3D U interpolation:'
         call spline_calc(regu3d,depthint,nlons,nlats,ongrid,regusp3d,nzlevel)
         print *
         call pack_short(regusp3d,nlons,nlats,ndeep,undef,ncid,var_id,var_info(i))

! V-velocity

         call def_netcdf_var(ncid,var_id,vars_3d,3,var_info(i2))
         call to_regugrid(hy3d2(:,:,1:kdm),regu3d)
         write(6,'(a)',advance='no') '3D V interpolation:'
         call spline_calc(regu3d,depthint,nlons,nlats,ongrid,regusp3d,nzlevel)
         print *
         call pack_short(regusp3d,nlons,nlats,ndeep,undef,ncid,var_id,var_info(i2))
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ssh
      i =var_info_getindex('ssh')
      if (i/=-1 .and. var_info(i)%process) then
         call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i))
         ! Experiment with add_offset and scale_factor
         call read_field2d(trim(fbase),'ssh     ',hy2d   ,idm,jdm,0,undef)
         where(hy2d/=undef) hy2d=hy2d/nrmem
         call to_regugrid(hy2d,regu2d)
         !basinssh=sum(ssh,mask=(abs(ssh-undef)>1e-4))
         !basinssh=basinssh/max(1,count(abs(ssh-undef)>1e-4))
         !where ((abs(ssh-undef)>1e-4)) ssh=ssh-basinssh
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dump mld1
      i =var_info_getindex('mld1')
      if (i/=-1 .and. var_info(i)%process) then
         ! MLD1
         call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i))
         call to_regugrid(mld1(:,:),regu2d) !mld1 processed above....
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dump mld1
      i =var_info_getindex('mld2')
      if (i/=-1 .and. var_info(i)%process) then
         ! MLD2
         call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i))
         call to_regugrid(mld2(:,:),regu2d)
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))
      end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Surface drag - must be rotated
      i =var_info_getindex('taux')
      i2=var_info_getindex('tauy')
      if (i/=-1.and.var_info(i)%process .and.  i2/=-1.and.var_info(i2)%process) then

         call read_field2d(trim(fbase),'taux    ',hy2d  ,idm,jdm,0,undef)
         call read_field2d(trim(fbase),'tauy    ',hy2d2  ,idm,jdm,0,undef)
         where(hy2d2/=undef) hy2d2=hy2d2/nrmem
         where(hy2d /=undef) hy2d =hy2d /nrmem

         call rotate(hy2d(:,:),hy2d2(:,:),   modlat,modlon,nx,ny,'m2l')

!Eastward atm stress
         call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i))
         call to_regugrid(hy2d(:,:),regu2d)
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))

! Northward Atmospheric stress on surface
         call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i2))
         call to_regugrid(hy2d2(:,:),regu2d)
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i2))
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Surface net heat flux
      i =var_info_getindex('surflx')
      if (i/=-1.and.var_info(i)%process) then
         call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i))
         call read_field2d(trim(fbase),'surflx  ',hy2d,idm,jdm,0,undef)
         where(hy2d/=undef) hy2d=hy2d/nrmem
         call to_regugrid(hy2d(:,:),regu2d)
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Surface sw heat flux
      i =var_info_getindex('swflx')
      if (i/=-1.and.var_info(i)%process) then
         call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i))
         call read_field2d(trim(fbase),'swflx   ',hy2d,idm,jdm,0,undef)
         where(hy2d/=undef) hy2d=hy2d/nrmem
         call to_regugrid(hy2d(:,:),regu2d)
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Evaporation minus precipitation
      i =var_info_getindex('emnp')
      if (i/=-1.and.var_info(i)%process) then
         call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i))
         call read_field2d(trim(fbase),'emnp    ',hy2d  ,idm,jdm,0,undef)
         where(hy2d/=undef) hy2d=hy2d/nrmem
         call to_regugrid(-hy2d(:,:)*1000.,regu2d)
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Barotropic streamfunction - must be rotated - dump ub/vb if necessary
      i =var_info_getindex('ubavg')
      i2=var_info_getindex('vbavg')
      i3=var_info_getindex('bsfd')
      if ((i /=-1.and.var_info(i )%process .and.  i2/=-1.and.var_info(i2)%process) .or. &
          (i3/=-1.and.var_info(i3)%process)) then

         call read_field2d(trim(fbase),'ubavg   ',hy2d  ,idm,jdm,0,undef)
         call read_field2d(trim(fbase),'vbavg   ',hy2d2  ,idm,jdm,0,undef)
         where(hy2d /=undef) hy2d =hy2d/nrmem
         where(hy2d2/=undef) hy2d2=hy2d2/nrmem

         print '(a)','Calculating streamfunction'
         if (i3/=-1.and.var_info(i3)%process) then
            call strmf_eval(idm,jdm,strmf,hy2d,hy2d2,depths,mqlat,mqlon)
            call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i3))
            call to_regugrid(strmf,regu2d)
            call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i3))
         end if
         call rotate(hy2d(:,:),hy2d2(:,:),   modlat,modlon,nx,ny,'m2l')

         if (i/=-1.and.var_info(i)%process) then
            call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i))
            call to_regugrid(hy2d,regu2d)
            call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))
         end if

         if (i2/=-1.and.var_info(i2)%process) then
            call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i2))
            call to_regugrid(hy2d2,regu2d)
            call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i2))
         end if
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ice thickness
      i =var_info_getindex('hice')
      if (i/=-1.and.var_info(i)%process) then
         call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i))
         call read_field2d(trim(fbase),'hice    ',hy2d  ,idm,jdm,0,undef)
         where(hy2d/=undef) hy2d=hy2d/nrmem
         call to_regugrid(hy2d(:,:),regu2d)
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ice concentration
      i =var_info_getindex('fice')
      if (i/=-1.and.var_info(i)%process) then
         call def_netcdf_var(ncid,var_id,vars_2d,2,var_info(i))
         call read_field2d(trim(fbase),'fice    ',hy2d  ,idm,jdm,0,undef)
         where(hy2d/=undef) hy2d=hy2d/nrmem
         call to_regugrid(hy2d(:,:),regu2d)
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ice velocity - must be rotated
      i =var_info_getindex('uice')
      i2=var_info_getindex('vice')
      if ((i /=-1.and.var_info(i )%process .and.  i2/=-1.and.var_info(i2)%process))then

         call read_field2d(trim(fbase),'uice    ',hy2d ,idm,jdm,0,undef)
         call read_field2d(trim(fbase),'vice    ',hy2d2,idm,jdm,0,undef)
         where(hy2d /=undef) hy2d =hy2d /nrmem
         where(hy2d2/=undef) hy2d2=hy2d2/nrmem
         hy2d3=hy2d
         hy2d4=hy2d2

         call rotate(hy2d3(:,:),hy2d4(:,:),   modlat,modlon,nx,ny,'m2l')
         where(hy2d ==undef) hy2d3=undef
         where(hy2d2==undef) hy2d4=undef
         hy2d =hy2d3
         hy2d2=hy2d4

         call to_regugrid(hy2d(:,:),regu2d)
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))

         call to_regugrid(hy2d2(:,:),regu2d)
         call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i2))
      end if

      !if  (.false.) then
      !   call handle_err(nf90_redef(ncid))
      !   call handle_err(NF90_DEF_VAR(ncid,'wice',NF90_SHORT,vars_2d,var_id))
      !   call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name','Ice Speed'))
      !   call handle_err(NF90_PUT_ATT(ncid,var_id,'units','m s-1'))
      !   call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',0_2))
      !   call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name',  &
      !                   'sea_ice_speed'))
      !   call handle_err(NF90_PUT_ATT(ncid,var_id,'add_offset',0.))
      !   call handle_err(NF90_PUT_ATT(ncid,var_id,'scale_factor',1.))
      !   call handle_err(NF90_PUT_ATT(ncid,var_id,'_FillValue',0_2)) ! Short
      !   call handle_err(nf90_enddef(ncid))
      !   call read_field2d(trim(fbase),'wice    ',hy2d ,idm,jdm,0,undef)
      !   where(hy2d/=undef) hy2d=hy2d/nrmem
      !   call to_regugrid(hy2d(:,:),regu2d)
      !   call pack_short(regu2d,nlons,nlats,1,undef,ncid,var_id,var_info(i))
      !end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !Close up
     call handle_err(NF90_CLOSE(ncid))
     print '(a)','Generated '//trim(ncfil)
     print *
     print *


     ! Start deallocating fileds now..
     deallocate(hy3d    )
     deallocate(hy3d2   )
     deallocate(hy3d3   )
     deallocate(regu3d  )
     deallocate(regusp3d)
     deallocate(strmf   )
     deallocate(mld1    )
     deallocate(mld2    )
     deallocate(hy2d    )
     deallocate(hy2d2   )
     deallocate(hy2d3   )
     deallocate(hy2d4   )
     deallocate(regu2d  )
     deallocate(depthint)

     !call zaioempty ! deallocated za - fields - ready for next file

  enddo ! ifile-loop

      

end program p_daily2regunc
