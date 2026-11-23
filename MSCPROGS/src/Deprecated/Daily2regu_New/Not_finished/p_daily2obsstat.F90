!
!If you make changes to this program, update the version number!!!

! This routine estimates the rms values over a set region. RMS is estimated
! wrt observations


program p_daily2obsstat
   use mod_common
   use mod_daily_average
   use mod_year_info
   use mod_sphere_tools
   use m_rotate
   use m_strmf_eval
   use m_year_day
   use m_mersea_prepare
   use m_get_depth_interval
   use m_get_regions
   use m_get_daily_infiles
   use m_tecfld
   use m_datetojulian

   use mod_measurement
   use mod_meanssh
   use m_read_mean_ssh
   use m_get_obs_d
   use m_get_nrobs_d
   implicit none

   integer*4, external :: iargc

   ! Version info of this program
   character(len=*), parameter :: cver = 'V0.1'
   real            , parameter :: undef=1e-9

   real, parameter :: epsiloon=1e-7

   ! Input file with days
   integer          :: maxfiles
   character(len=50), dimension(:), pointer :: daily_files
   character(len=20), dimension(:), pointer :: name_region

   ! 
   integer :: point_counter,numpoints,iregion,j,i,iint, region_counter
   integer :: nregion, npoints, maxpoints

   integer d2j,j2d
   ! 
   real, pointer :: lon_region(:,:), lat_region(:,:)
   real, allocatable, dimension(:,:) :: RMS_diff, Autocorr, AVE_MOD, AVE_diff
   real :: modval,omdiff,tmpval,obsval,pntlat, pntlon
   integer, pointer :: numvertex_region(:)
   real :: tmpmask(idm,jdm)
   real, allocatable :: int_dp(:,:,:)
   real, allocatable :: tmplon(:),tmplat(:)
   integer :: tmpnvert



   ! Z intervals..
   integer :: nzint
   real, pointer :: zinterval(:,:)

   !
   real :: upper , lower 



   type(year_info)    :: rt,rti,rtd,rtb       ! Time info type
   integer            :: year,day,hour,find,k,ios,imem,nrmem,ifile, &
                         numfiles, counter, ilon, jlat, refyear,month
   logical            :: ex, lok

   real, dimension(:,:,:), allocatable :: uinterval,vinterval,sinterval,tinterval,dinterval
   real, dimension(idm,jdm,kdm) :: utot, vtot
   real, dimension(idm,jdm) :: ssh, strmf, mld1, mld2, tmp, dsaln, emnp,totmask
   real, dimension(idm,jdm,kdm+1) :: pnew


   character(len=7)   :: tag7
   character(len=5)   :: cjj,obsid
   character(len=8)   :: dateinfo
   character(len=3)   :: rungen     ! Run ID
   character(len=4)   :: cyy        ! Year in chars
   character(len=3)   :: cdd,creg,cint,cbox        ! Day in chars
   character(len=2)   :: chh        ! Hour in chars
   character(len=20)  :: depthfile  ! Name of bathymetry file
   character(len=40)  :: tmpchar,vartitle,fname,fileobs
   character(len=160) :: tmpcharlong
   character(len=70)  :: fileout,persfile
   character(len=20)  :: rmstype


   integer :: year1, year2, day1, day2
   real    :: time1, time2
   logical :: lprocess,lmersea, l_analysis

   type(daily_averages) :: da1


   ! Observation state
   type(measurement) :: obs
   integer, allocatable :: nrobs_region(:,:)
   integer :: obsinregion,obsinint,nrobs
   integer :: reclO

   real :: a1, a2, a3, a4, ftime
   integer :: ipiv, jpiv, numint, iobs, retval, ind




   ! This sends hycom routines output to "fort.99"
   print *,'HYCOM routines output sent to "fort.66"'
   lp=66


   ! Is this generated for MERSEA ?
   lmersea = .false.


   if (iargc()>=1 .and. iargc()<=2 ) then

      ! Get RMS "type"
      call getarg(1,rmstype) 

      if (iargc()==2) then
         call getarg(2,tmpchar) 
         if ( trim(tmpchar) == "MERSEA" ) then
            print *,'MERSEA file names will be used!'
            lmersea = .true.
         else
            print *,'Unknown naming convention '//trim(tmpchar)
         end if
      end if
   elseif (iargc()>2 .or. iargc()==0) then
      print *,'Only one or two  arguments !!'
      stop '(daily2regunc)'
   end if

   ! So far -- two allowed rms-types
   if ( trim(rmstype)/="INNOV1" .and. trim(RMSTYPE)/="PERSISTENCE") then
      print *,'Allowed rms types are INNOV1 and PERSISTENCE'
      stop '(daily2obsstat)'
   end if






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization stuff for model

   ! Get bathymetry - lon/lat fields and sets up bigrid
   call get_grid 

   ! Get file size from a c-routine call. the string sent is null-terminated
   ! in order for the c-routine to know its end.
   tmpcharlong=daily_files(1);
   finduf=index(daily_files(1),'.uf')
   findab=max(index(daily_files(1),'.a'),index(daily_files(1),'.b'))
   !print *,findab,finduf
   if (finduf>0) then ! Assume all files are ".uf" files
      tmpcharlong(len_trim(tmpcharlong)+1:len_trim(tmpcharlong)+1)=char(0)
      call fsize(est_fsize,tmpcharlong)
      if (est_fsize<0) then
         print *,'The file size cant be estimated .. '
         print *,'You probably gave a non-existent file in infiles.daily'
         stop '(daily2regunc)'
      end if
      print *,'File size of '//trim(tmpcharlong)//' :',est_fsize
      nz=guess_nz(est_fsize)
      kk=nz
      kdm=nz

      ! Make sure .a files and .uf files are not mixed...
      do ifile=1,numfiles
         findab=index(daily_files(ifile),'.a')
         findab=max(findab,index(daily_files(ifile),'.b'))
         if(findab>0) then
            print *,'NB! Do NOT mix old (.uf) and new (.[ab]) files'
            stop
         end if
      end do

      print *,'Estimated z-dimension is ',nz
      if (nz>40 .or. nz<1) then
         print *,'Estimated z-dimension is > 40 or < 1'
         print *,'This probably means the file size estimate failed'
         print *,'Check the size of the first file in infiles.daily'
         stop '(daily2regunc)'
      end if

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


   ! Allocation ....
   allocate(temp(idm,jdm,kdm))
   allocate(saln(idm,jdm,kdm))
   allocate(utot(idm,jdm,kdm))
   allocate(vtot(idm,jdm,kdm))
   allocate(p   (idm,jdm,kdm+1))




   print *
   print *
   print *


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization stuff. 3 infiles are required:
! infiles.daily   -- Which daily average files to process
! depths.daily    -- Depth INTERVALS used 
! depths.daily    -- Region specification


!!!!!!!!!!!!!!!!!! Infile specification !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

   call get_daily_infiles(daily_files,numfiles,maxfiles)

!!!!!!!!!!!!!!!!!! Depths specification !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

   call get_depth_interval(zinterval,nzint)

!!!!!!!!!!!!!!!!!! Region specification !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            
   call  get_regions(lon_region,lat_region,numvertex_region,name_region,nregion,maxpoints)
   allocate(tmplon(1:maxval(numvertex_region)))
   allocate(tmplat(1:maxval(numvertex_region)))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! More Initialization stuff. Allocate grids and initialize




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ready to step through files to process

   ! Were ready...
   print *
   do ifile=1,numfiles
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read raw data -- Check for rms type (persistence or ordinary ...)

      ! In case RMS type is "PERSISTENCE", we get the data from the
      ! analysis corresponding to the current bulletin date... 
      ! The daily infiles only serve as a guideline for the times where 
      ! we calculate the rms values

      write (6,'(a)',advance='no') 'Reading from '//trim(daily_files(ifile))

      ! Check for type ...
      finduf=index(daily_files(ifile),'.uf')
      findab=max(index(daily_files(ifile),'.a'),index(daily_files(ifile),'.b'))
      if (finduf>0) then ! old type

         ! daily_average state
         lok=da_allocate(da)

         readok=daily_average_read(da,trim(daily_files(ifile)),1)

         rti=da%rtinit
         rtd=da%rtdump
         nrmem=da%nrmem

         ! Average -- best guess ...
         if (nrmem==0) then
            print '(a)', 'No members read from '//trim(daily_files(ifile))
            stop '(rawdaily2netcdf)'
         end if

         ! Convert from daily_averages to type used here
         p=0.; p(:,:,2:kdm+1)=da%p(:,:,:)
         utot=da%u ; vtot=da%v ; 
         saln(:,:,1:kdm)=da%s ; 
         temp(:,:,1:kdm)=da%t
         taux(:,:)=da%taux ;
         tauy(:,:)=da%tauy; 
         surflx=da%heatflux
         sswflx=da%swflx ; 
         salflx=da%salflx ; 
         emnp=da%emnp ; 
         ssh = da%ssh
         ubavg(:,:)=da%ub ; 
         vbavg(:,:)=da%vb

         ! Free up memory in the daily average type -- we dont need it anymore
         lok=da_deallocate(da)

      else if (findab>0) then
         print *,'ab-files in implementation'

        call read_field3d(trim(fbase),'utot    ',utot  ,idm,jdm,kdm,undef)
         where(utot/=undef) utot=utot/nrmem

         call read_field3d(trim(fbase),'vtot    ',vtot  ,idm,jdm,kdm,undef)
         where(vtot/=undef) vtot=vtot/nrmem

         call read_field3d(trim(fbase),'saln    ',saln  ,idm,jdm,kdm,undef)
         where(saln/=undef) saln=saln/nrmem

         call read_field3d(trim(fbase),'temp    ',temp  ,idm,jdm,kdm,undef)
         where(temp/=undef) temp=temp/nrmem

         call read_field3d(trim(fbase),'pres    ',p(:,:,2:kdm+1) ,idm,jdm,kdm,undef)
         where(p/=undef)
            p=p/nrmem
         elsewhere
            p=p/onem
         end if

      end if



      rti=da1%rtinit
      rtb=da1%rtdump
      rtd=da1%rtdump




      ! We now have the averaged estimate
      write (6,'(i4,a)') nrmem,' members'
      da_rtinit=rti
      da_rtdump=rtd







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! --- Set filename for output file. This is done already here, because
! --- some MERSEA files do not really need to be processed.

      !Set new filename (change ".uf" to ".nc")
      ! This is the DEFAULT filename. MERSEA might change it below
      print '(a)','Start File name generation'
      find=index(daily_files(ifile),'.uf')

      ! lprocess tells us wether this file should be processed.
      ! It is modified by MERSEA, processin, but is set to "on" 
      ! as default for all other files
      lprocess=.true.
      l_analysis=.false.

      ! For MERSEA
      if (lmersea) then

         call mersea_prepare(daily_files(ifile)(1:3),rtb,rti,rtd,lprocess,l_analysis)

      end if







      if (lprocess) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Start processing the raw  daily data

         ! For zero layers !

         ! MERSEA bug... files processed earlier had midpoint of first layer
         ! in stead of lower interface. This is a quick hack to circumvent this
         call interface_bug(p,pnew)
         p=pnew





         !Get Layer thickness 
         do k=1,kdm
            dp(:,:,k) = p(:,:,k+1)-p(:,:,k)
         enddo


         ! Baroclinic velocities
         do k=1,kdm
            u (:,:,k)   = utot (:,:,k) - ubavg(:,:,1)
            v (:,:,k)   = vtot (:,:,k) - vbavg(:,:,1)
         end do

         ! Rotate velocities
         print '(a)','Rotating velocities'
         do k=1,kdm
            call rotate(utot(:,:,k),vtot(:,:,k),modlat,modlon,nx,ny,'m2l')
            call rotate(u   (:,:,k),v   (:,:,k),modlat,modlon,nx,ny,'m2l')
         end do






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get observations here -- all observations must be of the same type
! in observations.uf


         !print *,rtd
         d2j=datetojulian(rtd%iyy,rtd%imm,rtd%idm,1950,1,1)
         print *,'datetojulian:',d2j

         !
         call juliantodate(d2j,year,month,day,1950,1,1)
         print *,'juliantodate:',year,month,day
         write(cjj,'(i5.5)') d2j

         ! time in fraction of year
         ftime = rtd%iyy+float(rtd%idd)/rtd%daysinyear



         ! Get number of observations
         nrobs = get_nrobs_d(rtd%cyy//rtd%cmm//rtd%cdm)
         if (nrobs==0) then
            print *,'No observations -- Skipping to next file!'
            print *
            cycle
         else
            print *,' Nrobs= ',nrobs
         end if

         ! Read first observation
         fileobs='observations'//rtd%cyy//rtd%cmm//rtd%cdm//'.uf'
         inquire(exist=ex,file=trim(fileobs))
         if (.not.ex) then
            print *,'File '//trim(fileobs)//' does not exist'
            stop
         end if
         inquire(iolength=reclO)obs
         open(11,file=trim(fileobs),form='unformatted',access='direct',recl=reclO)
         read(11,rec=1,iostat=ios)obs
         if (ios/=0) then
            print *,'Read error on '//trim(fileobs)
            print *,'iostat is ',ios
            stop
         end if
         obsid=obs%id
         print *,obsid
         print *,obs
         

         ! Set number of levels we actually use (numint)
         if (trim(obsid) == 'SST' .or. trim(obsid) == 'TIDE') then
            numint=1
         elseif (trim(obsid) == 'SLA') then
            call  read_mean_ssh(meanssh)
            numint=1
         elseif (trim(obsid) == 'TEMP' .or. trim(obsid) == 'SALN') then
            numint=nzint 
         else
            print *,'Unknown observation '//trim(obsid)
            stop
         end if
         !rewind(11)


         allocate(RMS_diff(nregion,numint))
         allocate(Autocorr(nregion,numint))
         allocate(AVE_diff(nregion,numint))
         allocate(AVE_MOD (nregion,numint))
         allocate(nrobs_region(nregion,numint))
         RMS_diff=0.
         autocorr=0.
         ave_diff=0.
         ave_mod=0.
         nrobs_region=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate model  to observation points -- Temporary calculation
! of rms values, bias average and model average in different boxes and depth
! intervals
!
! NB: After prep_obs, all observations are on model grid
         print *,'Interpolating observations in regions and intervals'
         print *,'Observations are:',obsid,' '
         open (12,file='regions.'//trim(obsid)//'.dat',status='replace')
         do iobs=1,nrobs
            
            read(11,rec=iobs,iostat=ios)obs

            if (ios/=0) then
               print *,'Read error on '//trim(fileobs)
               print *,'iostat is ',ios
               stop
            end if

            if (mod(iobs,nrobs/5)==0) print *,'Progress...',iobs,nrobs

            if (trim(obs%id) /= trim(obsid)) then
               print *,'All observations must be the same !'
               print *,obsid, obs%id
               stop
            end if


            ! Check if observation is in a region -- this also works 
            ! If there are overlapping regions
            do iregion=1,nregion

               obsinregion=0
               obsinint   =0


               ! NB -- lon_region and lat_region are pointers... 
               ! Convert to reals to pass on to "inbox" routine.
               ! (Failure to do that WILL inflate memory ...
               tmpnvert=numvertex_region(iregion)
               tmplon(1:tmpnvert)  =lon_region(iregion,1:tmpnvert)
               tmplat(1:tmpnvert)  =lat_region(iregion,1:tmpnvert)
               if (inbox(tmplon(1:tmpnvert),tmplat(1:tmpnvert), &
                         tmpnvert, obs%lon, obs%lat )) then
                  obsinregion=iregion
               end if

               if     (trim(obsid) == 'SST' .or. trim(obsid) == 'SLA' .or.  trim(obsid)=='TIDE') then
                  obsinint=1
               elseif (trim(obsid) == 'TEMP'.or. trim(obsid) == 'SALN') then
                  do iint=1,nzint
                     if (obs%depths>zinterval(iint,1) .and. obs%depths<zinterval(iint,2)) then
                        obsinint=iint
                     end if
                  end do

               else
                  print *,'Unknown observation '//trim(obsid)
                  stop
               end if

               if (obsinregion/=0 .and. obsinint/=0) then


                  ! obs attributes
                  ipiv = obs%ipiv
                  jpiv = obs%jpiv
                  a1   = obs%a1
                  a2   = obs%a2
                  a3   = obs%a3
                  a4   = obs%a4

                  ! Find model layer...
                  ka1=0; ka2=0; ka3=0; ka4=0;
                  do klevel=1,kdm
                     if (obs%depths>p(ipiv  ,jpiv  ,k) .and. obs%depths<=p(ipiv  ,jpiv  ,k) ) &
                        ka1=k
                     if (obs%depths>p(ipiv+1,jpiv  ,k) .and. obs%depths<=p(ipiv+1,jpiv  ,k) ) &
                        ka2=k
                     if (obs%depths>p(ipiv+1,jpiv+1,k) .and. obs%depths<=p(ipiv+1,jpiv+1,k) ) &
                        ka3=k
                     if (obs%depths>p(ipiv  ,jpiv+1,k) .and. obs%depths<=p(ipiv  ,jpiv+1,k) ) &
                        ka4=k
                  end do

                  ! Adjust coeffs depending on wether surrounding points are
                  ! over depth interval
                  if (ka1==0) then
                     ka1=1
                     a1=0
                  endif 
                  if (ka2==0) then
                     ka2=1
                     a2=0
                  endif 
                  if (ka3==0) then
                     ka3=1
                     a3=0
                  endif 
                  if (ka4==0) then
                     ka4=1
                     a4=0
                  endif 
                  asum=a1+a2+a3+a4
                  a1=a1/asum
                  a2=a2/asum
                  a3=a3/asum
                  a4=a4/asum

                  if (asum<.5) then
                     print *,'daily2obsstat'
                     print *,'ka error ..'
                     stop
                  end if

                  


                  ! Calculate model values for observation point
                  if (trim(obsid) == 'SST') then
                     modval= &
                        temp(ipiv  ,jpiv  ,1)*a1 + &
                        temp(ipiv+1,jpiv  ,1)*a2 + &
                        temp(ipiv+1,jpiv+1,1)*a3 + &
                        temp(ipiv  ,jpiv+1,1)*a4 
                  elseif (trim(obsid) == 'SLA' .or. trim(obsid)=='TIDE') then
                     modval =&
                        ssh     (ipiv  ,jpiv  )*a1 + &
                        ssh     (ipiv+1,jpiv  )*a2 + &
                        ssh     (ipiv+1,jpiv+1)*a3 + &
                        ssh     (ipiv  ,jpiv+1)*a4 
                     tmpval= &
                        meanssh(ipiv  ,jpiv  )*a1 + &
                        meanssh(ipiv+1,jpiv  )*a2 + &
                        meanssh(ipiv+1,jpiv+1)*a3 + &
                        meanssh(ipiv  ,jpiv+1)*a4 
                     modval=modval -tmpval
                  elseif (trim(obsid) == 'TEMP') then

                     modval= &
                        temp(ipiv  ,jpiv  ,ka1)*a1 + &
                        temp(ipiv+1,jpiv  ,ka2)*a2 + &
                        temp(ipiv+1,jpiv+1,ka3)*a3 + &
                        temp(ipiv  ,jpiv+1,ka4)*a4 
                  elseif (trim(obsid) == 'SALN') then
                     modval= &
                        saln(ipiv  ,jpiv  ,ka1)*a1 + &
                        saln(ipiv+1,jpiv  ,ka2)*a2 + &
                        saln(ipiv+1,jpiv+1,ka3)*a3 + &
                        saln(ipiv  ,jpiv+1,ka4)*a4 
                  else
                     print *,'Unknown observation '//trim(obsid)
                     stop
                  end if
                  omdiff= modval-obs%d
                  obsval= obs%d


                  nrobs_region(iregion,obsinint)=nrobs_region(iregion,obsinint) +1

                  ! RMS... 
                  RMS_diff(iregion,obsinint) =  RMS_diff(iregion,obsinint) + omdiff**2

                  ! Average difference (obs points)
                  AVE_diff(iregion,obsinint) =  AVE_diff(iregion,obsinint) + omdiff

                  ! Average model values in observation points
                  AVE_MOD (iregion,obsinint) =  AVE_MOD (iregion,obsinint) + modval

                  ! For testing Region masks
                  write (12,'(4f10.2,2i5)') obs%lon,obs%lat, modval,obsval,iregion,obsinint
                  !print *,obs%id,obsinregion,obsinint,obs%lon,obs%lat,modval, &
                  !   obs%d,nrobs_region(iregion,obsinint)
               end if
            end do
         end do
         close(11)
         close(12)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final calculation of RMS values for different regions and depth intervals.
! Dump the result to formatted ascii files. One file for each time. Not
! practical per se, but through use of post-processing routines, it is...
         do iregion = 1,nregion
         do iint=1,numint
            npoints = nrobs_region(iregion,iint)
            if (npoints>0) then

               ! RMS... 
               RMS_diff(iregion,iint) = sqrt(  RMS_diff(iregion,iint) / npoints)

               ! Average difference (obs points)
               AVE_diff(iregion,iint) =  AVE_diff(iregion,iint) / npoints

               ! Average model values in observation points
               AVE_MOD (iregion,iint) =  AVE_MOD (iregion,iint) / npoints


               write (creg,'(i3.3)') iregion
               write (cint,'(i3.3)') iint
               if (obsid == 'SLA' .or. obsid == 'SST' .or. obsid == 'TIDE' ) then
               fileout='RMS.'//trim(rmstype)//'.b'//rtb%cyy//rtb%cmm//rtb%cdm//'_f'//   &
                       rtd%cyy//rtd%cmm//rtd%cdm//'.'//         &
                       trim(obsid)//'.r'//creg//'.dat'
               else
               fileout='RMS.'//trim(rmstype)//'.b'//rtb%cyy//rtb%cmm//rtb%cdm//'_f'//   &
                       rtd%cyy//rtd%cmm//rtd%cdm//'.'//         &
                       trim(obsid)//'.r'//creg//'.i'//cint//'.dat'
               end if


               ! Create file header
               write (cbox,'(i3.3)') numvertex_region(iregion)
               open(10,file=trim(fileout),status='replace')
               write(10,'(2a)') 'Region :', trim(name_region(iregion))
               write(10,'(2a)') 'Data Type:',trim(obsid)
               write(10,'(a,'//cbox//'f10.3)') 'Region - Longitude :', &
                  lon_region(iregion,1:numvertex_region(iregion))
               write(10,'(a,'//cbox//'f10.3)') 'Region - Latitude  :', &
                  lon_region(iregion,1:numvertex_region(iregion))
               !if (trim(obsid) == 'SST' .or. trim(obsid)=='SLA') then
               if (obsid == 'SLA' .or. obsid == 'SST' .or. obsid == 'TIDE' ) then
                  write(10,'(2a)') 'Depth interval     :','-'
               elseif (trim(obsid) == 'SALN' .or. trim(obsid)=='TEMP') then
                  write(10,'(a,2f9.2)') 'Depth interval :',zinterval(iint,1),zinterval(iint,2)
               else
                  print *,'Unknown observation '//trim(obsid)
                  stop
               end if
               write(10,'(2a)') 'Column Legend :', '1)bulletin 2)forecast 3) Float-time' //&
                  '4) Model average 5)Bias 6)RMS difference 7) nrobs'
               write(10,'(a,f10.4,3f10.2,i6)')  &
                  rtb%cyy//rtb%cmm//rtb%cdm//'   '//rtd%cyy//rtd%cmm//rtd%cdm, &
                  ftime, AVE_MOD (iregion,iint), AVE_diff(iregion,iint),       &
                  RMS_diff(iregion,iint), npoints
               close(10)

            end if 
         end do 
         end do 






         deallocate(RMS_diff)
         deallocate(Autocorr)
         deallocate(AVE_diff)
         deallocate(AVE_MOD )
         deallocate(nrobs_region)

      end if !lprocess
      print *
      print *


   enddo ! ifile-loop
   print *,'done ..'

end program p_daily2obsstat
