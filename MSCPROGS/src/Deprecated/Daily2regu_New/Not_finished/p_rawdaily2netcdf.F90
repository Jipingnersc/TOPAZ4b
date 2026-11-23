program p_rawdaily2netcdf
   use mod_xc, only: mod_xc_idm => idm, mod_xc_jdm => jdm
   use mod_za , only: zaiost
   use mod_common
   use mod_year_info
   use netcdf
   use m_handle_err
   use m_get_daily_infiles
   use m_pack_short
   use m_get_grid
   use mod_read_dailyab
   implicit none

   integer :: est_fsize

   ! These are variable names in NetCDF files
   character(len=*), parameter :: temname    ='temperature'
   character(len=*), parameter :: sshname    ='ssh'
   character(len=*), parameter :: salname    ='salinity'
   character(len=*), parameter :: ifacename  ='interface'
   character(len=*), parameter :: vname      ='v'
   character(len=*), parameter :: uname      ='u'
   character(len=*), parameter :: tauyname   ='tauy'
   character(len=*), parameter :: tauxname   ='taux'
   character(len=*), parameter :: mldname    ='mlp'
   character(len=*), parameter :: qtname     ='qtot'
   character(len=*), parameter :: qswname    ='qsw'
   character(len=*), parameter :: fwname     ='emp'
   character(len=*), parameter :: streamfname='bsfd'
   character(len=*), parameter :: mltemname  ='mld'
   character(len=*), parameter :: vbname     ='vb'
   character(len=*), parameter :: ubname     ='ub'
   character(len=*), parameter :: ficename   ='fice'
   character(len=*), parameter :: hicename   ='hice'

   real, parameter :: undef=-1e20

   ! 2D vars -- on model grid
   real   , dimension(:,:), allocatable :: &
      varout2d,mostrf, ssh, strmf, mld1, mld2, tmp,  emnp, &
      ficem,hicem,hsnwm, uice,vice, surflx,sswflx,salflx,meanssh, &
      ubavg,vbavg,taux,tauy,wice
   integer, dimension(:,:), allocatable :: varout2dint

   ! 3D vars -- on model grid
   real, dimension(:,:,:), allocatable :: &
      p,  temp, saln, utot, vtot, varout3d,dp, intface


   !integer         , parameter  :: maxfiles=100
   !character(len=*), parameter  :: infile_daily = 'infile.daily'

   character(len=60), pointer, dimension(:) :: daily_files
   integer :: maxfiles, numfiles

   type(year_info)    :: rt,rtd   ! Time info type
   integer            :: year,day,hour,find,k,ios,imem,nrmem,ifile, &
                         counter, finda, findb, finduf,findab
   real               :: tinfo
   logical            :: ex,readok
   type (year_info) :: rti,rtb


   character(len=7)   :: tag7
   character(len=8)   :: dateinfo
   character(len=3)   :: rungen     ! Run ID
   character(len=4)   :: cyy        ! Year in chars
   character(len=3)   :: cdd        ! Day in chars
   character(len=2)   :: chh        ! Hour in chars
   character(len=20)  :: depthfile  ! Name of bathymetry file
   character(len=40)  :: ncfil      ! Name of new netcdf file
   character(len=80)  :: tmpchar,vartitle,fname,fbase
   character(len=160) :: tmpcharlong
   logical :: lok


   ! Netcdf variable ids
   integer :: i_dimension_id,j_dimension_id,k_dimension_id 
   integer :: lon_id,lat_id,int_id,temp_id,sal_id,u_id,v_id,pb_id, &
              pbu_id,pbv_id,taux_id,tauy_id,hflx_id, emnp_id,ssh_id
   integer :: ncid ! netcdf id of netcdf file
   integer :: vars_2d(2), vars_3d(3)  ! holder for dimension ids


   ! Get infiles
   call get_daily_infiles(daily_files,numfiles,maxfiles)

   ! Get bathymetry - lon/lat fields and sets up bigrid
   call get_grid 


   ! Get file size from a c-routine call. the string sent is null-terminated
   ! in order for the c-routine to know its end.
   tmpcharlong=daily_files(1);
   finduf=index(daily_files(1),'.uf')
   findab=max(index(daily_files(1),'.a'),index(daily_files(1),'.b'))
   !print *,findab,finduf
   if (finduf>0) then ! Assume all files are ".uf" files
      print *,'This file no longer works on .uf (old) type files'
      stop '(daily2regunc)'
   elseif (findab>0) then ! Assume all files are ".[ab]" files
      findab=max(index(daily_files(1),'.a'),index(daily_files(1),'.b'))

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
   allocate(p    (idm,jdm,kdm+1),temp(idm,jdm,kdm),dp(idm,jdm,kdm))
   allocate(saln (idm,jdm,kdm),taux (idm,jdm))
   allocate(tauy (idm,jdm),ubavg(idm,jdm),vbavg(idm,jdm))
   allocate(sswflx(idm,jdm), salflx(idm,jdm), meanssh(idm,jdm),surflx(idm,jdm))

   ! mod_common_ice
   allocate(ficem (idm,jdm),hicem (idm,jdm),vice (idm,jdm), wice(idm,jdm))
   allocate(uice  (idm,jdm))
   allocate(intface  (idm,jdm,1:kdm+1))


   ! Variables local to this routine
   allocate(utot(idm,jdm,kdm), vtot(idm,jdm,kdm))
   allocate(ssh(idm,jdm), strmf(idm,jdm), mld1(idm,jdm), mld2(idm,jdm))
   allocate(tmp(idm,jdm), emnp(idm,jdm))
   !print *,'got grid vars'






   ! Were ready...
   print *
   do ifile=1,numfiles
      
      write (6,'(a)',advance='yes') 'Reading from '//trim(daily_files(ifile))

      ! Check for type ...
      finduf=index(daily_files(ifile),'.uf')
      findab=max(index(daily_files(ifile),'.a'),index(daily_files(ifile),'.b'))
      if (finduf>0) then ! old type
         print *,'This file no longer works on .uf (old) type files'
         stop '(daily2regunc)'
      else if (findab>0) then


         ! Read different fields - returns undef if not defined
         fbase=daily_files(ifile)(1:findab-1)
         call daily_average_read_header(trim(fbase),rti,rtd,nrmem,idm,jdm,kdm)

         call read_field2d(trim(fbase),'ubavg   ',ubavg ,idm,jdm,0,undef)
         where(ubavg/=undef.or.depths<.5) ubavg=ubavg/nrmem

         call read_field2d(trim(fbase),'vbavg   ',vbavg ,idm,jdm,0,undef)
         where(vbavg/=undef.or.depths<.5) vbavg=vbavg/nrmem

         call read_field2d(trim(fbase),'ssh     ',ssh   ,idm,jdm,0,undef)
         where(ssh/=undef.or.depths<.5) ssh=ssh/nrmem

         call read_field2d(trim(fbase),'surflx  ',surflx,idm,jdm,0,undef)
         where(surflx/=undef.or.depths<.5) surflx=surflx/nrmem

         call read_field2d(trim(fbase),'swflx   ',sswflx,idm,jdm,0,undef)
         where(sswflx/=undef.or.depths<.5) sswflx=sswflx/nrmem

         call read_field2d(trim(fbase),'salflx  ',salflx,idm,jdm,0,undef)
         where(salflx/=undef.or.depths<.5) salflx=salflx/nrmem

         call read_field2d(trim(fbase),'emnp    ',emnp  ,idm,jdm,0,undef)
         where(emnp/=undef.or.depths<.5) emnp=emnp/nrmem

         call read_field2d(trim(fbase),'taux    ',taux  ,idm,jdm,0,undef)
         where(taux/=undef.or.depths<.5) taux=taux/nrmem

         call read_field2d(trim(fbase),'tauy    ',tauy  ,idm,jdm,0,undef)
         where(tauy/=undef.or.depths<.5) tauy=tauy/nrmem

         call read_field2d(trim(fbase),'fice    ',ficem ,idm,jdm,0,undef)
         where(ficem/=undef.or.depths<.5) ficem=ficem/nrmem

         call read_field2d(trim(fbase),'hice    ',hicem ,idm,jdm,0,undef)
         where(hicem/=undef.or.depths<.5) hicem=hicem/nrmem

         call read_field2d(trim(fbase),'uice    ',uice  ,idm,jdm,0,undef)
         where(uice/=undef.or.depths<.5) uice=uice/nrmem

         call read_field2d(trim(fbase),'vice    ',vice  ,idm,jdm,0,undef)
         where(vice/=undef.or.depths<.5) vice=vice/nrmem

         call read_field2d(trim(fbase),'wice    ',wice  ,idm,jdm,0,undef)
         where(wice/=undef.or.depths<.5) wice=sqrt(wice/nrmem)

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

         do k=1,nz
            where(depths<.5)
               temp(:,:,k)=undef
               saln(:,:,k)=undef
               utot(:,:,k)=undef
               vtot(:,:,k)=undef
               p   (:,:,k)=undef
            end where
         end do

         

      end if



      ! We now have the averaged estimate
      write (6,'(a,i4,a)') 'Fields based on ',nrmem,' members'

         



      !Get Layer thickness 
      do k=1,kdm
         dp(:,:,k) = p(:,:,k+1)-p(:,:,k)
      enddo


      !Set new filename (chaneg ".uf" to ".nc"
      finda=index(daily_files(ifile),'.a')
      findb=index(daily_files(ifile),'.b')
      find = max(finda,findb)
      ncfil=daily_files(ifile)
      ncfil(find:find+2)='.nc'


! --- Start NC definition section 

      !Create netcdf file and fill it with data
      if (NF90_CREATE(ncfil,NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening the netcdf file'
         stop '(restart2netcdf)'
      else
         write(lp,'(a)') 'Creating '//trim(ncfil)
      end if


      ! Define dimensions -- idm,jdm,kdm
      call handle_err(NF90_DEF_DIM(ncid,'idm',idm,i_dimension_id))
      call handle_err(NF90_DEF_DIM(ncid,'jdm',jdm,j_dimension_id))
      call handle_err(NF90_DEF_DIM(ncid,'kdm',kdm,k_dimension_id))

      ! Define dimension "holders"
      vars_3d=(/i_dimension_id,j_dimension_id,k_dimension_id/)
      vars_2d=(/i_dimension_id,j_dimension_id/)

! Define variables -- put attributes

      !Global attributes
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'title',        &
         'NetCDF dataset generated from NERSC HYCOM restart file'// &
         trim(fname)))
      call date_and_time(date=dateinfo)
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Dataset_generation_date', &
         dateinfo(1:8)))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Model_Year',rt%iyy))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Model_DayOfYear',rt%idd))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Model_DayOfYear_Note', &
         '1. January is DayOfYear 0'))
      call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Model_Hour',rt%ihh))
      

      ! Longitude
      vartitle='Longitude of pressure gridcell center'
      call handle_err(NF90_DEF_VAR(ncid,'lon',NF90_FLOAT,vars_2d,lon_id))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'title',trim(vartitle)))
      call handle_err(NF90_PUT_ATT(ncid,lon_id,'units','Degrees East'))

      ! Latitude
      vartitle= 'Latitude of pressure gridcell center'
      call handle_err(NF90_DEF_VAR(ncid,'lat',NF90_FLOAT,vars_2d,lat_id))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'title',trim(vartitle)))
      call handle_err(NF90_PUT_ATT(ncid,lat_id,'units','Degrees North'))

      ! Layer interfaces
      call handle_err(NF90_DEF_VAR(ncid,ifacename,NF90_SHORT,vars_3d,int_id))
      call handle_err(NF90_PUT_ATT(ncid,int_id,'title','Bottom Interface of layer'))
      call handle_err(NF90_PUT_ATT(ncid,int_id,'units','m'))
      call handle_err(NF90_PUT_ATT(ncid,int_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,int_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,int_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,int_id,'add_offset',0.))

      ! Temperature
      call handle_err(NF90_DEF_VAR(ncid,temname,NF90_SHORT,vars_3d,temp_id))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'title','Temperature of layer'))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'units','Celsius'))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,temp_id,'add_offset',0.))


      ! Salinity
      call handle_err(NF90_DEF_VAR(ncid,salname,NF90_SHORT,vars_3d,sal_id))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'title','Salinity of layer'))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'units','ppm'))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,sal_id,'add_offset',0.))


      ! U-velocity
      call handle_err(NF90_DEF_VAR(ncid,uname,NF90_SHORT,vars_3d,u_id))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'title','Baroclinic U-velocity of layer'))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'units','m/s'))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,u_id,'add_offset',0.))

      ! V-velocity
      call handle_err(NF90_DEF_VAR(ncid,vname,NF90_SHORT,vars_3d,v_id))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'title','Baroclinic V-velocity of layer'))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'units','m/s'))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,v_id,'add_offset',0.))

      ! Barotropic U-velocity
      call handle_err(NF90_DEF_VAR(ncid,ubname,NF90_SHORT,vars_2d,pbu_id))
      call handle_err(NF90_PUT_ATT(ncid,pbu_id,'title','Barotropic U-velocity'))
      call handle_err(NF90_PUT_ATT(ncid,pbu_id,'units','m/s'))
      call handle_err(NF90_PUT_ATT(ncid,pbu_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,pbu_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,pbu_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,pbu_id,'add_offset',0.))

      ! Barotropic V-velocity
      call handle_err(NF90_DEF_VAR(ncid,vbname,NF90_SHORT,vars_2d,pbv_id))
      call handle_err(NF90_PUT_ATT(ncid,pbv_id,'title','Barotropic V-velocity'))
      call handle_err(NF90_PUT_ATT(ncid,pbv_id,'units','m/s'))
      call handle_err(NF90_PUT_ATT(ncid,pbv_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,pbv_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,pbv_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,pbv_id,'add_offset',0.))

      ! Atmospheric stress on surface
      call handle_err(NF90_DEF_VAR(ncid,tauxname,NF90_SHORT,vars_2d,taux_id))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'title','Atmospheric surface drag'))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'units','xxx'))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,taux_id,'add_offset',0.))

      ! Atmospheric stress on surface
      call handle_err(NF90_DEF_VAR(ncid,tauyname,NF90_SHORT,vars_2d,tauy_id))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'title','Atmospheric surface drag'))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'units','xxx'))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,tauy_id,'add_offset',0.))

      ! Atmospheric stress on surface
      call handle_err(NF90_DEF_VAR(ncid,sshname,NF90_SHORT,vars_2d,ssh_id))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'title','Sea Surface Height'))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'units','m'))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,ssh_id,'add_offset',0.))


      ! Heatflux at surface
      call handle_err(NF90_DEF_VAR(ncid,qtname,NF90_SHORT,vars_2d,hflx_id))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'title','Net heat flux'))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'units','xxx'))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,hflx_id,'add_offset',0.))

      call handle_err(NF90_ENDDEF(ncid))
      !print *,'end def'
! --- End definition section 

!stop '(must check for undefined points !)'

      intface=p(:,:,:)/onem
      do k=1,kdm
         where (abs(depths - p(:,:,k)/onem)<1.) 
            temp   (:,:,k) = undef
            saln   (:,:,k) = undef
            utot   (:,:,k) = undef
            vtot   (:,:,k) = undef
         end where
         where (depths<.5) intface(:,:,k)=undef
      end do

      
      where (depths<.5) 
         ubavg(:,:)  =undef
         vbavg(:,:)  =undef
         ssh  (:,:)  =undef
         emnp (:,:)  =-undef
         taux (:,:)  =undef
         tauy (:,:)  =undef
         surflx=undef
         intface(:,:,kdm+1)  = undef
      end where



! --- Start variable put section 
      call handle_err(NF90_PUT_VAR(ncid,lon_id,modlon))
      call handle_err(NF90_PUT_VAR(ncid,lat_id,modlat))
      call pack_short(intface(:,:,2:kdm+1),idm,jdm,kdm,undef,ncid,int_id,ifacename)
      call pack_short(temp(:,:,1:kdm),idm,jdm,kdm,undef,ncid,temp_id,temname)
      call pack_short(saln(:,:,1:kdm),idm,jdm,kdm,undef,ncid,sal_id,salname)
      call pack_short(utot(:,:,1:kdm),idm,jdm,kdm,undef,ncid,u_id,uname)
      call pack_short(vtot(:,:,1:kdm),idm,jdm,kdm,undef,ncid,v_id,vname)
      call pack_short2D(vbavg(:,:),idm,jdm,undef,ncid,pbv_id,vbname)
      call pack_short2D(ubavg(:,:),idm,jdm,undef,ncid,pbu_id,ubname)
      call pack_short2D(taux(:,:),idm,jdm,undef,ncid,taux_id,tauxname)
      call pack_short2D(tauy(:,:),idm,jdm,undef,ncid,tauy_id,tauyname)
      call pack_short2D(surflx(:,:),idm,jdm,undef,ncid,hflx_id,qtname)
      call pack_short2D(ssh(:,:),idm,jdm,undef,ncid,ssh_id,sshname)
! --- End variable put section 

      !Close up
      call handle_err(NF90_CLOSE(ncid))

      print *
      print *

   enddo ! ifile-loop

end program p_rawdaily2netcdf
