program p_icedrift
! Calculates ice drift over a given number of days. input
! is julian days for start and end of integration. Also needed
! is start time of baseline run. This is needed since the drift is
! based on daily averages.
!
! Ice drift is originating in "ifremer" drift points, and model velocities
! are interpolated to model lagrangian point at each time step. Runge-Kutta
! 2nd order method
   use netcdf
   use nfw_mod
   use mod_xc
   use mod_za
   use mod_grid , only : plon,plat,depths, scpx, scpy, get_grid
   use mod_year_info
   use mod_confmap
   use mod_hycomfile_io , only: hycomfile, initHF, HFReadField,  &
      getfiletype, undef
   use m_ncvar_dims
   use m_ncvar_read
   use m_rk2
   implicit none

   character(len=80) :: tmparg,driftfile,flnmdep,flnmgrd,cline,hycfile, &
                        drfilemodel1,drfilemodel2, ftype
   character(len=80), allocatable :: drfile_list(:)
   character(len= 3) :: rungen,cmem
   character(len= 8) :: cbtime
   integer :: byear,bjuld

   real :: time,rkdt

   real :: hmina,hminb,hmaxa,hmaxb,hmin,hmax
   real :: a1,a2,a3,a4
   real :: x1,y1,x2,y2
   real :: wx,wy,wt

   integer ::i, lgth, j, i2, j2, mini2, minj2
   integer ::ip,jp
#if defined(IARGC)
   integer*4, external :: iargc
#endif
   type(year_info) :: rtdump,rtinit

   ! Drift input fields
   integer :: drnx, drny, drdt
   real, dimension(:,:), allocatable ::  &
      drlon,drlat,x,y,drzon,drmer,qflag, &
      x0,y0,   &
      drmer_mgrid, &
      drzon_mgrid, tmp, x_x0_mgrid,y_y0_mgrid, &
      modmer, modzon
   real*4, dimension(:,:), allocatable ::  modzonio,modmerio
   integer, dimension(:,:), allocatable ::  qmask
   real, dimension(1) :: scalefac, fillval, addoffset

   real, dimension(:), allocatable :: tmplon, tmplat 
   integer                         :: tmptime(2)
   integer                         :: ncid, lon_id, lat_id,mer_id,zon_id,time_id
!   real                            :: fillval,scalefac,addoffset

   real, dimension(:,:,:), allocatable :: blcoeff
   integer, dimension(:,:), allocatable :: pivoti, pivotj
   logical :: storepivots,ex

   ! RK velocity estimates
   real, dimension(:,:,:), allocatable :: u,v
   real, dimension(:,:,:), allocatable :: qvar1,qvar2
   real                                :: Vzon,Vmer 

   real :: delt
   integer ljulday,fjulday
   integer :: istep,irk,istep2
   integer :: irecl
   logical :: first

   integer :: fyear, fmonth, fdm
   integer :: lyear, lmonth, ldm

   integer , dimension(NF90_MAX_VAR_DIMS) :: dimsizes
   integer :: recdim, numdim
   integer :: nsteps, imem,findnc

   integer :: summask

   real :: dist, mindist
   real :: lat_o_1,lon_o_1,lat_o_2,lon_o_2

   real :: flon, flat,tlon,tlat
   integer :: nfiles
   integer :: fcount
   character(len=8) :: mode
   character(len=8) :: filetype
   integer :: indxa, indxb

   integer :: iyear, iday, imonth, iyear2, idayom
   character(len=20) :: cdate

   integer :: kdm,nrmem

   real*4, allocatable, dimension(:,:) :: iou,iov
   real, external :: spherdist
   type(hycomfile) :: hfile

#if defined(MATLAB)
#include "fintrf.h"
   MWPOINTER :: mxCreateNumericMatrix, mxGetPr, mxClassIDFromClassName, matopen,  &
      mxCreateDoubleMatrix, matPutVariableAsGlobal, mp, pa1
   integer matputvariable, matclose
   real*8, dimension(:,:), allocatable :: iofld,iofld2
   integer :: status,iret
   integer,parameter :: IKIND=8
#endif




   storepivots=.true.

!-------------------------------------------------------------------
!  Process input arguments
!-------------------------------------------------------------------
   driftfile=''
   if (iargc()==5) then
       
      ! Rungen for model id
      call getarg(1,tmparg) ; read(tmparg,*) rungen

      ! bulletin refyear from model
      call getarg(2,tmparg) ; read(tmparg,*) byear

      ! bulletin time from model
      call getarg(3,tmparg) ; read(tmparg,*) bjuld

      ! Input ice drift file. Needed for initial grid positions,
      ! start and end times of model integrations
      call getarg(4,tmparg) ; read(tmparg,*) driftfile

      call getarg(5,tmparg) ; read(tmparg,*) imem

      filetype='DRIFT'

   elseif (iargc()==0) then

      ! This mode reads position from input file 
      inquire(exist=ex,file='icedrift.in')
      drnx=1
      drny=1
      if (ex) then
         open(10,file='icedrift.in',status='old') 
         read(10,*) mode
         if (trim(mode)=='POSITION') then
            read(10,*) flon,flat
         else if (trim(mode)=='CERSAT') then
            read(10,*) driftfile
            print *, '(icedrift_cersat:CERSAT mode not yet implemented for infile)'
            call exit(1)
         else 
            print *, '(icedrift_cersat:CERSAT or POSITION mode )'
            call exit(1)
         end if
         print *,'bef read imem'
         read(10,*) imem
         read(10,*) filetype
         read(10,*) nfiles
         allocate(drfile_list(nfiles))
         do i=1,nfiles
            read(10,*) drfile_list(i)
         end do
      else
         print *,'icedrift_cersat:zero args - icedrift.in must be present' 
         print *,'---Template----'
         print *,'POSITION     # Position or CERSAT'
         print *,'EAST NORTH   # lon lat pos or name of driftfile'
         print *,'1            # always set to 1'
         print *,'DAILY        # file type , DAILY or DRIFT'
         print *,'NFILES       # total number of files'
         print *,'file_1       '
         print *,'.......      '
         print *,'file_NFILES  '
         print *,'---End Template----'
         call exit(1)
      end if
   else
      print *,'icedrift_cersat calculates ice drift based on a specified input file' 
      print *,'Containing daily average files, or using a CERSAT ice drift file' 
      print *
      print *,'Called with five or zero args' 
      print *,'Usage - 5 args:'
      print *,'icedrift_cersat rungen bull_year bull_day cersat_file member'
      print *
      print *,'Or, call with no args to calculate individual ice particle'
      print *,'drift. You must then have a icedrift.in file present (see '
      print *,'Infiles in src directory for an example '
      call exit(1)
   end if
   write(cmem,'(i3.3)') imem


!-------------------------------------------------------------------
!  Process date - find correct ice drift file
!-------------------------------------------------------------------

   if (trim(driftfile)/='') then
      ! btime is in format yyyymmdd
      call year_day(real(bjuld),byear,rtinit,'ecmwf')
      write(rtinit%cdm,'(i2.2)') rtinit%idm+1

      ! check the datat time contained in the observation NC file
      ! open the nc file:
      call nfw_open(driftfile, nf_nowrite, ncid)
      call nfw_inq_varid(driftfile,ncid,'time',time_id)
      
      call nfw_get_var_int(driftfile, ncid, time_id,tmptime)
      ljulday= maxval(tmptime)/24.      ! included this day for 2-day drifting 
      fjulday= ljulday

      call juliantodate(ljulday,lyear,lmonth,ldm,1900,1,1)
      call juliantodate(fjulday,fyear,fmonth,fdm,1900,1,1)

!      ! Convert to julian days ref start year
      fjulday=datetojulian(fyear,fmonth,fdm,fyear,1,1)
      ljulday=datetojulian(lyear,lmonth,ldm,fyear,1,1)

      Print *,'First year, month, day - last year, month, day'
      print *,fyear,fmonth,fdm
      print *,lyear,lmonth,ldm

      ! Get dimensions of drift file
      !call ncvar_dims(driftfile,'zonal_motion',dimsizes,numdim,recdim)
      call ncvar_dims(driftfile,'eastward_sea_ice_velocity',dimsizes,numdim,recdim)

      ! Name for resulting model drift file
      findnc=index(driftfile,'.nc')

      ! Drift file - on data grid - using zonal/meridional motions
      drfilemodel1=rungen//driftfile(1:findnc-1)//'.uf'

      ! Drift file - on model grid . usin model u/v motions
      ! KAL -- This is only used by ensstat
      drfilemodel2=rungen//driftfile(1:findnc-1)//'_diag.uf' ! Used by ensstat

   else

      
      rungen=drfile_list(1)(1:3)
      fjulday=0
      ljulday=nfiles-1

      drfilemodel1=rungen//'drift.uf'
      drfilemodel2=rungen//'drift_diag.uf'
   end if


!-------------------------------------------------------------------
!  read observations file
!-------------------------------------------------------------------

   if (trim(driftfile)/='') then

      drnx=dimsizes(1)
      drny=dimsizes(2)
      drdt=dimsizes(3)
      allocate(tmplon(drnx))
      allocate(tmplat(drny))

      allocate(drlon(drnx,drny))
      allocate(drlat(drnx,drny))
      allocate(drmer(drnx,drny))
      allocate(drzon(drnx,drny))
      allocate(qflag(drnx,drny))
      allocate(qmask(drnx,drny))
      allocate(qvar1(drnx,drny,2))
      allocate(qvar2(drnx,drny,2))

      call nfw_inq_varid(driftfile,ncid,'longitude',lon_id)
      call nfw_inq_varid(driftfile,ncid,'latitude',lat_id)
      call cersat_readfield(driftfile, ncid, lon_id, tmplon, drnx,0)
      call cersat_readfield(driftfile, ncid, lat_id, tmplat, drny,0)

      do i=1,drnx
        drlat(i,:)=tmplat(:)
      end do
      do j=1,drny
        drlon(:,j)=tmplon(:)
      end do
      do i=1,drnx
         do j=1,drny
            if (drlon(i,j)<0) then
              drlon(i,j)=drlon(i,j)+360 
            endif
         end do
      end do
  
      qmask=1
      call nfw_inq_varid(driftfile,ncid,'eastward_sea_ice_velocity',zon_id)
      call nfw_inq_varid(driftfile,ncid,'northward_sea_ice_velocity',mer_id)

      if (drdt==2) then
         call cersat_readfield(driftfile, ncid, zon_id, qvar1,drnx*drny*2,0)
      else
         call cersat_readfield(driftfile, ncid, zon_id, qvar1(:,:,1),drnx*drny,0)
         qvar1(:,:,2)=qvar1(:,:,1)
      endif
      call nfw_get_att_double(driftfile, ncid, zon_id, '_FillValue', fillval)
      call nfw_get_att_double(driftfile, ncid, zon_id, 'scale_factor', scalefac)
      call nfw_get_att_double(driftfile, ncid, zon_id, 'add_offset', addoffset)
      Vzon=fillval(1) * scalefac(1) + addoffset(1)
      qvar1=qvar1*scalefac(1)+addoffset(1)

      if (drdt==2) then
         call cersat_readfield(driftfile, ncid, mer_id, qvar2,drnx*drny*2,0)
      else
         call cersat_readfield(driftfile, ncid, mer_id, qvar2(:,:,1),drnx*drny,0)
         qvar2(:,:,2)=qvar2(:,:,1)
      endif
      call nfw_get_att_double(driftfile, ncid, mer_id, '_FillValue', fillval)
      call nfw_get_att_double(driftfile, ncid, mer_id, 'scale_factor', scalefac)
      call nfw_get_att_double(driftfile, ncid, mer_id, 'add_offset', addoffset)
      Vmer=fillval(1) * scalefac(1) + addoffset(1)
      qvar2=qvar2*scalefac(1)+addoffset(1)

      qmask=1
      drzon=qvar1(:,:,1); drmer=qvar2(:,:,1)
      where (abs(drzon+0.65536)<0.001.or.abs(drmer+0.65536)<0.001) 
         qmask=0
      end where
      drzon=qvar1(:,:,2); drmer=qvar2(:,:,2)
      where (abs(drzon+0.65536)<0.001.or.abs(drmer+0.65536)<0.001) 
         qmask=0
      end where
      deallocate(tmplon,tmplat)
      deallocate(qvar1,qvar2)

   else

      ! Read data from drift file
      allocate(drlon(drnx,drny))
      allocate(drlat(drnx,drny))
      allocate(drmer(drnx,drny))
      allocate(drzon(drnx,drny))
      allocate(qflag(drnx,drny))
      allocate(qmask(drnx,drny))

      drlon=flon
      drlat=flat
      drmer=0.
      drzon=0.
      qflag=0.
      qmask=0.
   end if


!-------------------------------------------------------------------
!  read model grid info
!-------------------------------------------------------------------
   ! Read grid dimensions
   call xcspmd()
   call zaiost()
   call get_grid
   allocate(iou (idm,jdm))
   allocate(iov (idm,jdm))
   ! Get longitude and latidude of cell midpoints - also depths

   call initconfmap(idm,jdm)


   print *,'Drift details read',drnx,drny

! ------------------------------------------------------------------
! Allocate remaining arrays
! ------------------------------------------------------------------
   allocate(u   (idm,jdm,2))
   allocate(v   (idm,jdm,2))

!-------------------------------------------------------------------
!  Calculate initial positions
!-------------------------------------------------------------------
   
   allocate(x0   (drnx,drny))
   allocate(y0   (drnx,drny))
   allocate(x    (drnx,drny))
   allocate(y    (drnx,drny))

   x0=undef
   y0=undef
   do j=1,drny
   do i=1,drnx

      if (qmask(i,j)==1) then
         call oldtonew(drlat(i,j),drlon(i,j),x0(i,j),y0(i,j))
         call pivotp(y0(i,j),x0(i,j),ip,jp)
         if (ip>0.and.ip<idm.and.jp>0 .and. jp<jdm-1) then
            x0(i,j)=ip 
            y0(i,j)=jp
      !      call ll2gind(drlon(i,j),drlat(i,j),x0(i,j),y0(i,j))

         end if
      end if

   end do
   end do

   print *,'Initial Grid positions ',fjulday,'-',ljulday,':'
!   print *,minval(x0,x0/=undef),maxval(x0,x0/=undef)
!!   print *,minval(y0,y0/=undef),maxval(y0,y0/=undef)


!-------------------------------------------------------------------
!  RK forward integration - model drift on model grid
!-------------------------------------------------------------------

   ! Start the time loop from start to end time. We use 6 hours
   ! as the RK(2) time step
   delt=86400/4.
!  nsteps=ljulday-fjulday+1    ! for 2 days drifting
   nsteps=3       ! for 3 days drift
   time=fjulday
   first=.true.
   x=x0;
   y=y0;
   fcount=1

   if (drnx*drny==1) then

      ! Deduce date from file name
      hycfile=drfile_list(fcount)
      read(hycfile(19:22),'(i4)') iyear
      read(hycfile(24:26),'(i3)') iday

      ! Actual date
      call juliantodate(iday,iyear2,imonth,idayom,iyear,1,1)

      ! Encode to string
      write(cdate,'(i2.2,"/",i2.2,"-",i4.4)') idayom,imonth,iyear2
      !print *,cdate

      open(999,file='tsdrift.asc',status='replace')
      write(999,'(4f14.4,x,a)') x(1,1),y(1,1),drlon(1,1),drlat(1,1),trim(cdate)
   end if

   ! This cycles days - requires day/nstep to be integer ....
   u=0.
   v=0.
   do istep=1, nsteps
      !print *,'day step',istep
      ! Initially ; Read new field into time slot 0 and 1

      if (trim(driftfile)/='') then
         ! Get date from julian day
!         call year_day(time+1+(istep-1),fyear,rtdump,'ecmwf')
!         write(rtdump%cdm,'(i2.2)') rtdump%idm+1

         call year_day(time+(istep-1),fyear,rtdump,'ecmwf')
         write(rtdump%cdm,'(i2.2)') rtdump%idm

         ! HYCOM drift file
         hycfile=rungen//'DAILY_'//rtinit%cyy//'_'//rtinit%cdd// &
                 '_'//rtdump%cyy//'_'//rtdump%cdd// &
                 '_ICEDRIFT.uf'
         !print *,rtdump%cyy,rtdump%cmm,rtdump%cdm
         print *,hycfile
      else
         hycfile=drfile_list(fcount)
         fcount=fcount+1
      end if


      ! print *,'reading '//trim(hycfile)
      if (trim(filetype)=='DRIFT') then
         inquire(iolength=irecl)iou,iov
         open(10,file=trim(hycfile),access='direct',status='old', &
                 recl=irecl)
         read(10,rec=imem) iou,iov
         close(10)
         if (istep==1) then
            u(:,:,1)=iou
            v(:,:,1)=iov
            u(:,:,2)=iou
            v(:,:,2)=iov
         else
            u(:,:,1)=u(:,:,2)
            v(:,:,1)=v(:,:,2)
            u(:,:,2)=iou
            v(:,:,2)=iov
         endif
      else if (trim(filetype)=='DAILY') then
         hycfile=drfile_list(fcount)
         indxa=index(hycfile,'.a')
         indxb=index(hycfile,'.b')
         hycfile=hycfile(1:max(indxa,indxb)-1)

         !call daily_average_read_header(trim(hycfile),rtinit,rtdump,nrmem,idm,jdm,kdm)
         !call read_field2d(trim(hycfile),'uice    ',u(:,:,2),idm,jdm,0,undef)
         !call read_field2d(trim(hycfile),'vice    ',v(:,:,2),idm,jdm,0,undef)
         !u(:,:,2)=u(:,:,2)/nrmem
         !v(:,:,2)=v(:,:,2)/nrmem
         ftype=getfiletype(trim(hycfile)//'.a')
         call initHF(hfile,trim(hycfile)//'.a',trim(ftype))
         call HFReadField(hfile,u(:,:,1),idm,jdm,'uice    ',0,1)
         call HFReadField(hfile,v(:,:,1),idm,jdm,'vice    ',0,1)
         call HFReadField(hfile,u(:,:,2),idm,jdm,'uice    ',0,1)
         call HFReadField(hfile,v(:,:,2),idm,jdm,'vice    ',0,1)

      else
         print *,'icedrift_cersat:Unknown filetype '//trim(filetype)
         call exit(1)
      end if


      ! Deduce date from file name
      read(hycfile(19:22),'(i4)') iyear
      read(hycfile(24:26),'(i3)') iday

      ! Actual date
      call juliantodate(iday,iyear2,imonth,idayom,iyear,1,1)

      ! Encode to string
      write(cdate,'(i2.2,"/",i2.2,"-",i4.4)') idayom,imonth,iyear2
      !print *,cdate
      !print *, scpx(430,660),u(430,660,1), x(430,50),v(430,660,1), y(430,50),drnx

      ! One days worth of rk2 integration here
      call rk2(u,v,scpx,scpy,idm,jdm,x,y,drnx,drny,delt,undef)
      
      ! Dumps daily values of position
      if (drnx*drny==1) then
         call gind2ll(x(1,1),y(1,1),flon,flat)
         write(999,'(4f14.4,x,a)') x(1,1),y(1,1),flon,flat,trim(cdate)
      end if

   enddo ! All days loop

   if (drnx*drny==1) then
      close(999)
   end if

   !print *
   !print *,'Displacement max in grid coordinates:'
   !print *,'x:',maxval(abs(x-x0)),maxloc(y-y0)
   !print *,'y:',maxval(abs(y-y0))


   ! Find meridional / zonal motion
   allocate(modmer(drnx,drny))
   allocate(modzon(drnx,drny))
   allocate(modmerio(drnx,drny))
   allocate(modzonio(drnx,drny))
   modzon=0.
   modmer=0.
   do j=2,drny-1
   do i=2,drnx-1
      if (qmask(i,j)==1) then
      call gind2ll(x0(i,j),y0(i,j),lon_o_1,lat_o_1)
      call gind2ll(x (i,j),y (i,j),lon_o_2,lat_o_2)
      !print *,x0(i,j),y0(i,j),lon_o_1,lat_o_1
      !print *,x(i,j),y(i,j),lon_o_2,lat_o_2

      ! Put in same lon lat interval
      if (lon_o_2 >  lon_o_1+180.) lon_o_2=lon_o_2-360.
      if (lon_o_2 <= lon_o_1-180.) lon_o_2=lon_o_2+360.

      ! Meridional motion
      !modmer(i,j)=lat_o_2-lat_o_1
      modmer(i,j)=spherdist(lon_o_1,lat_o_1,lon_o_1,lat_o_2)* &
                  sign(1.,lat_o_2-lat_o_1)*1e-3

      ! Zonal motion 
      !modzon(i,j)=lon_o_2-lon_o_1
      modzon(i,j)=spherdist(lon_o_1,lat_o_1,lon_o_2,lat_o_1)* &
                  sign(1.,lon_o_2-lon_o_1)*1e-3

      endif

      !if (j==28.and.i==683) then  
      !    print *,x0(i,j),y0(i,j),lon_o_1,lat_o_1
      !    print *,x(i,j),y(i,j),lon_o_2,lat_o_2
      !    print *,lon_o_1-lon_o_2,lat_o_1-lat_o_2, lat_o_2
      !    print *, modzon(i,j)
      !    print *, modmer(i,j)
      !    print *, spherdist(1.0,1.0,1.0,2.0)
      !endif
   end do
   end do
   print *,'Displacement max in km:'
   print *,maxval(sqrt( ((modmer))**2  + ((modzon))**2 ))
   !print *,minval(modmer),maxval(modmer)
   !print *,minval(modzon),maxval(modzon)


#if defined (MATLAB)
if (trim(driftfile)/='') then
   ! Matlaby diags
   ! matlab file
   print *,'Dumping to matlab file'
   allocate(iofld (drnx,drny))
   allocate(iofld2(idm , jdm))
   mp=matopen('drift'//cmem//'.mat','w')

   iofld=drlon;
   pa1=mxCreateNumericMatrix(int(drnx,kind=IKIND),int(drny,kind=IKIND),mxClassIDFromClassName('double'),0)
   call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),int(drnx*drny,kind=IKIND))
   status = matPutVariable(mp, 'o_longitude', pa1)

   iofld=drlat;
   pa1=mxCreateNumericMatrix(int(drnx,kind=IKIND),int(drny,kind=IKIND),mxClassIDFromClassName('double'),0)
   call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),int(drnx*drny,kind=IKIND))
   status = matPutVariable(mp, 'o_latitude', pa1)

   iofld=drzon;
   pa1=mxCreateNumericMatrix(int(drnx,kind=IKIND),int(drny,kind=IKIND),mxClassIDFromClassName('double'),0)
   call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),int(drnx*drny,kind=IKIND))
   status = matPutVariable(mp, 'o_zonal', pa1)

   iofld=drmer;
   pa1=mxCreateNumericMatrix(int(drnx,kind=IKIND),int(drny,kind=IKIND),mxClassIDFromClassName('double'),0)
   call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),int(drnx*drny,kind=IKIND))
   status = matPutVariable(mp, 'o_merid', pa1)

   iofld=modzon;
   pa1=mxCreateNumericMatrix(int(drnx,kind=IKIND),int(drny,kind=IKIND),mxClassIDFromClassName('double'),0)
   call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),int(drnx*drny,kind=IKIND))
   status = matPutVariable(mp, 'm_zonal', pa1)

   iofld=modmer;
   pa1=mxCreateNumericMatrix(int(drnx,kind=IKIND),int(drny,kind=IKIND),mxClassIDFromClassName('double'),0)
   call mxCopyReal8ToPtr(iofld,mxGetPr(pa1),int(drnx*drny,kind=IKIND))
   status = matPutVariable(mp, 'm_merid', pa1)

   iret=matClose(mp)
end if
#endif


   ! Produce output file - ice drift over data file drift period
   ! -- on model grid -- This is our model observation, basically
   !inquire(iolength=irecl)x,y
   inquire(iolength=irecl)modzonio,modmerio
   open(10,file=trim(drfilemodel1),access='direct',status='unknown', &
           recl=irecl)
   !write(10,rec=imem) x-x0,y-y0
   modzonio=modzon
   modmerio=modmer
   write(10,rec=imem) modzonio,modmerio
   close(10)




! Only do this in nc mode ...
if (trim(driftfile)/='' .and. .false.) then
!-------------------------------------------------------------------
!  Needed by ensstat routine ..
!
!  Interpolate observered and modeled drift on model grid
!  Note that this is mostly for diag purposes. The observation grid
!  should be used for assimilation.
! 
! The pivot point calculations could be vastly improved 
!-------------------------------------------------------------------

   ! Interpolate observed drift to model grid 
   allocate(drzon_mgrid(idm,jdm))
   allocate(drmer_mgrid(idm,jdm))
   allocate(x_x0_mgrid (idm,jdm))
   allocate(y_y0_mgrid (idm,jdm))
   allocate(tmp        (idm,jdm))
   allocate(pivoti     (idm,jdm))
   allocate(pivotj     (idm,jdm))
   allocate(blcoeff    (idm,jdm,4))

   pivoti=1
   pivotj=1
   blcoeff=0.


   if (.not.storepivots .or. imem==1) then
      do j=1,jdm
      print *,j
      do i=1,idm

        ! Find closest x0 and y0 larger than i and j
         mindist=1e8
         mini2=1
         minj2=1
         do j2=1,drny-1
         do i2=1,drnx-1

            dist=sqrt( (x0(i2,j2)-i)**2+(y0(i2,j2)-j)**2)
            if (dist<mindist) then
              mindist=dist
              mini2=i2
              minj2=j2
            end if
         end do
         end do


         ! If point is valid, interpolate from obs to model grid
         i2=mini2
         j2=minj2
         if (mindist<sqrt(2.)) then
            summask = qmask(i2  ,j2  ) + qmask(i2+1,j2  )  &
                    + qmask(i2+1,j2+1) + qmask(i2  ,j2+1) 
         end if

         
         if (mindist<sqrt(2.)) then

            pivoti(i,j)=i2
            pivotj(i,j)=j2

!            ! bilinear coefficients - spatial weight
!            if (x0(i2+1,j2  )-x0(i2,j2)>0.) then 
!               wx=(i-x0(i2,j2))/x0(i2+1,j2  )-x0(i2,j2)
!            else
!               wx=.5
!               print *,'x error ',i,j,x0(i2,j2),x0(i2+1,j2)
!            end if
!
!            if (y0(i2  ,j2+1)-y0(i2,j2)>0.) then
!               wy=(j-y0(i2,j2))/abs(y0(i2  ,j2+1)-y0(i2,j2))
!            else
!               wy=.5
!               print *,'y error ',i,j,x0(i2,j2),y0(i2,j2+1)
!            end if
!
!            a1=(1-wx)*(1-wy);
!            a2=wx*(1-wy);
!            a3=wx*wy;
!            a4=(1-wx)*wy;
!

            ! Nearestpoint for now
            a1=1.
            a2=0.
            a3=0.
            a4=0.

            blcoeff(i,j,1)=a1
            blcoeff(i,j,2)=a2
            blcoeff(i,j,3)=a3
            blcoeff(i,j,4)=a4

            if (summask==4) then
               drzon_mgrid(i,j)=a1*drzon(i2  ,j2  ) +&
                                a2*drzon(i2+1,j2  ) +&
                                a3*drzon(i2+1,j2+1) +&
                                a4*drzon(i2  ,j2+1)

               drmer_mgrid(i,j)=a1*drmer(i2  ,j2  ) +&
                                a2*drmer(i2+1,j2  ) +&
                                a3*drmer(i2+1,j2+1) +&
                                a4*drmer(i2  ,j2+1)

               else
                  drzon_mgrid(i,j)=undef
                  drmer_mgrid(i,j)=undef
            end if

            ! Model drift on model grid
            x_x0_mgrid (i,j)=a1*(x(i2  ,j2  )-x0(i2  ,j2  )) +&
                             a2*(x(i2+1,j2  )-x0(i2+1,j2  )) +&
                             a3*(x(i2+1,j2+1)-x0(i2+1,j2+1)) +&
                             a4*(x(i2  ,j2+1)-x0(i2  ,j2+1))

            y_y0_mgrid (i,j)=a1*(y(i2  ,j2  )-y0(i2  ,j2  )) +&
                             a2*(y(i2+1,j2  )-y0(i2+1,j2  )) +&
                             a3*(y(i2+1,j2+1)-y0(i2+1,j2+1)) +&
                             a4*(y(i2  ,j2+1)-y0(i2  ,j2+1))
         else
            drzon_mgrid(i,j)=undef
            drmer_mgrid(i,j)=undef
            x_x0_mgrid(i,j) =0.
            y_y0_mgrid(i,j) =0.
         end if


      end do
      end do

      ! Store pivots
      if (imem==1 .and. storepivots) then
         open(10,file='icedrift_pivots.uf',form='unformatted',status='replace')
         write(10) pivoti,pivotj,blcoeff
         close(10)
      end if

      ! Rotate to model grid
      tmp = drzon_mgrid
      where (drzon_mgrid==undef)
         drzon_mgrid=0
         drmer_mgrid=0
      end where

      ! Produce drift from data file - rotate meridional and zonal motions
      ! to ocean model grid.
      call rotate(drzon_mgrid,drmer_mgrid,plat,plon,idm,jdm,'l2m')
      where (tmp==undef)
         drzon_mgrid=undef
         drmer_mgrid=undef
      end where

      ! Produce test file - obs ice drift over data file drift period  - 
      ! interpolated and rotated onto model grid
      print *,'Dump obs drift'
      open(10,file='obsdrift_mgrid.tst',status='replace')
      do j=1,jdm
      do i=1,idm
         if (drzon_mgrid(i,j)/=undef) then
            !write(10,'(2i5,4e15.4)'), i,j,drzon_mgrid(i,j),drmer_mgrid(i,j),  &
            !x_x0_mgrid(i,j), y_y0_mgrid(i,j)
            write(10,'(2i5,12e15.4)') i,j,drzon_mgrid(i,j),drmer_mgrid(i,j),  &
            x_x0_mgrid(i,j), y_y0_mgrid(i,j), &
            real(pivoti(i,j)), real(pivotj(i,j)), &
            blcoeff(i,j,1), blcoeff(i,j,2),  &
            blcoeff(i,j,3), blcoeff(i,j,4)
         end if
      end do
      end do
      close(10)
    

   else

      inquire(exist=ex,file='icedrift_pivots.uf')

      if (ex) then

         open(10,file='icedrift_pivots.uf',form='unformatted',status='old')
         read(10) pivoti,pivotj,blcoeff
         close(10)

         do j=1,jdm
         do i=1,idm

            a1=blcoeff(i,j,1)
            a2=blcoeff(i,j,2)
            a3=blcoeff(i,j,3)
            a4=blcoeff(i,j,4)

            i2=pivoti(i,j)
            j2=pivotj(i,j)

            ! Model drift on model grid
            x_x0_mgrid (i,j)=a1*(x(i2  ,j2  )-x0(i2  ,j2  )) +&
                             a2*(x(i2+1,j2  )-x0(i2+1,j2  )) +&
                             a3*(x(i2+1,j2+1)-x0(i2+1,j2+1)) +&
                             a4*(x(i2  ,j2+1)-x0(i2  ,j2+1))

            y_y0_mgrid (i,j)=a1*(y(i2  ,j2  )-y0(i2  ,j2  )) +&
                             a2*(y(i2+1,j2  )-y0(i2+1,j2  )) +&
                             a3*(y(i2+1,j2+1)-y0(i2+1,j2+1)) +&
                             a4*(y(i2  ,j2+1)-y0(i2  ,j2+1))
         end do
         end do
      end if
   end if


   ! Produce output file - ice drift over data file drift period
   inquire(iolength=irecl)x_x0_mgrid,y_y0_mgrid
   open(10,file=trim(drfilemodel2),access='direct',status='unknown', &
           recl=irecl)
   write(10,rec=imem) x_x0_mgrid,y_y0_mgrid
   close(10)
else
#if defined (AIX)
   call exit_(0)
#else
   call exit (0)
#endif
end  if

   end program

  subroutine cersat_readfield(fname, ncid, varid, v, vlen,Tscale)
    use nfw_mod
    implicit none

    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid,Tscale
    integer :: vlen
    real(8), intent(out) :: v(vlen)

    real(8) :: scale_factor(1)
    real(8) :: offset(1)
    call nfw_get_var_double(fname, ncid, varid, v)
    if (Tscale==1) then
      call nfw_get_att_double(fname, ncid, varid, 'scale_factor', scale_factor)
      call nfw_get_att_double(fname, ncid, varid, 'add_offset', offset)
      v = v * scale_factor(1) + offset(1)
    endif

  end subroutine cersat_readfield

