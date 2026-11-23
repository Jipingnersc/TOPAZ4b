program p_argocmp
   use mod_common
   use mod_year_info
   use mod_levitus
   use netcdf
   use m_year_day
   use m_datetojulian
   use m_get_grid
   use mod_xc, only: mod_xc_idm => idm, mod_xc_jdm => jdm
   use mod_za , only: zaiost, zaioempty
   use mod_read_dailyab
   use mod_read_rstab
   use mod_read_weekab
   use mod_spline_calc
   use m_read_field_wrapper
   use m_read_coriolis
   use m_dump_coriolis
   use m_initconfmap
   use m_pivotp
   use m_oldtonew
   implicit none

   !integer*4, external :: iargc

   ! Version info of this program
   character(len=*), parameter :: cver = 'V0.1'

   ! Input file with days
   character(len=80) :: hyc_file
   integer  :: maxfiles,numfiles

   type(year_info)    :: rt,rti,rtd,rtb       ! Time info type
   integer            :: year,day,hour,find,k,ios,nrmem,ifile, &
                         counter, refyear
   real               :: undef
   logical            :: ex

   real, allocatable, dimension(:,:  ) :: dp, pres, presold, saln, temp
   logical, allocatable, dimension(:) :: valid_profile

   character(len=8), dimension(:), pointer   :: floatid
   character(len=8)   :: dateinfo,cdate
   character(len=11)  :: tag7
   character(len=3)   :: rungen     ! Run ID
   character(len=80)  :: modcorfile,corfile    ! Name of new netcdf file
   character(len=80)  :: fbase,filename_ice,cprof
   character(len=40)  :: tmpchar,vartitle,fname
   character(len=40)  :: grid_mapping,coords_mapping
   character(len=160) :: tmpcharlong
   character(len=  4) :: char_fac
   character(len=3) :: cgg,cii
   character(len=12) :: clon, clat, cjuld
   character(len=4)  :: cyear
   character(len=2)  :: cmonth
   character(len=40) :: fulldir



   integer :: fname_iday,fname_iyear,fname_dday,fname_dyear
   integer :: finduf,findab,findrst,finddaily,findweek
   integer :: i, i2, i3,ilon,jlat
   integer :: hycfile
   integer :: tmpdim(2)

   integer :: igroup, istat,nstat

   real, pointer, dimension(:,:) :: csaln, cspres, ctemp, ctpres
   real, pointer, dimension(:,:) :: msaln , mtemp, msaln1, mtemp1, mpres1
   real, pointer, dimension(:,:) :: mpresd, msalnd, mtempd
   real, pointer, dimension(:,:) :: levsaln, levtemp

   ! Error flags - used locally by nersc
   logical, pointer, dimension(:,:) :: errorflag_sal,errorflag_temp

   ! invalid profile - per mersea definitions
   logical, pointer, dimension(:)   :: invalid_sal,invalid_temp

   real, pointer, dimension(:)   :: cslon,ctlon,cslat,ctlat,csjuld,ctjuld
   real :: lat_n,lon_n

   integer :: iprof, izlev
   integer :: nprof, nzlev
   integer :: ipiv,jpiv
   integer :: validk
   integer :: k2
   integer :: iyear, imonth, iday
   logical :: l_ncdump=.false.

   integer :: fid_num
   character(len=10) :: cfid_num

   undef = -999


   
   if (iargc()/=2 .and. iargc()/=3 ) then
      print *,'Usage:'
      print *,'  argocmp argo-file hycom-file [ncdump]'
      stop
   else
      call getarg(1,corfile)
      call getarg(2,hyc_file)
      if (iargc()==3) then
         call getarg(3,tmpchar)
         if (trim(tmpchar)=='-ncdump') then
            l_ncdump=.true.
         else
            print *,'Usage:'
            print *,'  argocmp argo-file hycom-file [-ncdump]'
            stop
         end if
      end if
   end if



   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization stuff for model  -- allocation galore
   print *

   ! Get bathymetry - lon/lat fields and sets up bigrid
   call get_grid 



! Get the coriolis schtuff - for now we can deal with one input file only
   !corfile='20060131_prof.nc'
   call read_coriolis(trim(corfile),csaln,cspres,cslon,cslat,csjuld,floatid, &
                      errorflag_sal, invalid_sal, 'PSAL',undef)
   call read_coriolis(trim(corfile),ctemp,ctpres,ctlon,ctlat,ctjuld,floatid, &
                      errorflag_temp,invalid_temp,'TEMP',undef)
   if (any(cspres/=ctpres)) then
      print *,'Temperature and salinity is defined at different levels. Fix me.'
      stop
   end if
   nzlev = size(csaln,1)
   nprof = size(csaln,2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start working with hycom data file

   !hyc_file='FORDAILY_2006_108_2006_108.a'
   write (6,'(a)',advance='yes') 'Reading from '//trim(hyc_file)
   rungen=hyc_file(1:3)

   ! Check for type ...
   finduf=index(hyc_file,'.uf')
   findab=max(index(hyc_file,'.a'),index(hyc_file,'.b'))
   if (finduf>0) then ! old type
      print *,'uf-files no longer supported by this routine ...'
      stop
   else if (.not. findab>0) then
      print *,'No .ab files'
      stop
   else
      ! We have a .ab-file. Now figure out what type.
      findrst  =index(hyc_file,'restart')
      finddaily=index(hyc_file,'DAILY')
      findweek =index(hyc_file,'AVE')
      if (findrst==4) then
         hycfile=1 ! 1 for restart files
      elseif (finddaily==4) then
         hycfile=2 ! 1 for daily files
      elseif (findweek==4) then
         hycfile=3 ! 1 for weekly files
      else
         print *,'Can not deduce file type from  file name'
         stop
      end if
   end if


   ! Read different fields - returns undef if not defined
   fbase=hyc_file(1:findab-1)
   if (hycfile==2) then
      call daily_average_read_header(trim(fbase),rti,rtd,nrmem,idm,jdm,kdm)
   else if (hycfile==1) then
      call rst_read_header(trim(fbase),rti,nrmem,idm,jdm,kdm)
      rtd=rti
      nrmem=1
   else if (hycfile==3) then
      call week_read_header(trim(fbase),rti,idm,jdm,kdm)
      rtd=rti
      nrmem=1
   end if

   ! Init file io
   mod_xc_idm=idm
   mod_xc_jdm=jdm
   nx=idm ! Get rid of these
   ny=jdm ! Get rid of these
   nz=kdm ! Get rid of these
   call zaiost
   call initconfmap


   ! We now have the averaged estimate
   if (hycfile==2) then ! daily
      write (6,'(a,i4,a)') 'Fields based on ',nrmem,' members'
   end if


   ! Start allocating fields now..
   allocate(dp       (idm,jdm))            ! Holds 2D vars on original grid
   allocate(pres    (idm,jdm))            ! Holds 2D vars on original grid
   allocate(presold (idm,jdm))            ! Holds 2D vars on original grid
   allocate(temp     (idm,jdm))            ! Holds 2D vars on original grid
   allocate(saln     (idm,jdm))            ! Holds 2D vars on original grid
   allocate(msaln1(kdm,nprof))
   allocate(mtemp1(kdm,nprof))
   allocate(mpres1(kdm,nprof))
   allocate(msalnd(nzlev,nprof))
   allocate(mtempd(nzlev,nprof))
   allocate(mpresd(nzlev,nprof))
   allocate(msaln(nzlev,nprof))
   allocate(mtemp(nzlev,nprof))




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   ! Days start from "0". Add 1 to get real 
   ! dates
   rtb=rti
   write(rtb%cdm,'(i2.2)') rtb%idm+1
   write(rtb%cdd,'(i3.3)') rtb%idd+1
   write(rtd%cdm,'(i2.2)') rtd%idm+1
   write(rtd%cdd,'(i3.3)') rtd%idd+1
   write(rti%cdm,'(i2.2)') rti%idm+1
   write(rti%cdd,'(i3.3)') rti%idd+1

   ! Start working through the model layers
   mpres1=0.
   pres=0.
   presold=0.
   allocate(valid_profile(nprof))
   valid_profile=.true.
   do k=1,kdm


      ! Cater to different needs
      if (hycfile==2) then
         call read_field2d(trim(fbase),'pres    ',pres ,idm,jdm,k,undef)
         pres=pres / ( onem * nrmem)
         dp=pres-presold
      elseif (hycfile==3) then
         ! mod_average typo, pres is actually dp
         call read_weekfield2d(trim(fbase),'pres    ',dp ,idm,jdm,k,undef)
         pres=presold+dp
      elseif (hycfile==1) then
         call read_rstfield2d(trim(fbase),'dp      ',dp ,idm,jdm,k,undef)
         dp=dp/onem
         pres=presold+dp
      end if

      ! Read temp & saln
      call read_field_wrapper_2d(saln ,'saln    ',trim(fbase),idm,jdm,k,hycfile,nrmem,undef,dp)
      call read_field_wrapper_2d(temp ,'temp    ',trim(fbase),idm,jdm,k,hycfile,nrmem,undef,dp)

      !print *,hycfile
      !print *,minval(dp,dp>0.),maxval(dp,dp<1e27)
      !print *,minval(temp),maxval((temp))

      ! Put into model profile array 1. Same vertical levels as model
      ! Go through coriolis profiles -- find corresponding points in model
      do iprof=1,nprof

         


         ! Get closest point in model
         call oldtonew(cslat(iprof),cslon(iprof),lat_n,lon_n)
         call pivotp(lon_n,lat_n,ipiv,jpiv)

         if ( ipiv>1 .and. ipiv < idm .and. jpiv > 1 .and. jpiv < jdm .and. &
             depths(ipiv,jpiv)>1.) then

            if (k==1) then 
               print *,'model / obs point :',modlon(ipiv,jpiv),modlat(ipiv,jpiv), &
               cslon(iprof),cslat(iprof)
            end if

            mtemp1(k,iprof)=temp(ipiv,jpiv)
            msaln1(k,iprof)=saln(ipiv,jpiv)
            mpres1(k,iprof)=pres(ipiv,jpiv)
         else if (.not.(ipiv>1.and.ipiv<idm.and.jpiv>1.and.jpiv<jdm)) then
            if (k==1) print *,' obs point out of domain :', cslon(iprof),cslat(iprof)
            valid_profile(iprof)=.false.
         elseif (depths(ipiv,jpiv)<1.) then
            if (k==1) print *,' obs point on model land :', cslon(iprof),cslat(iprof)
            valid_profile(iprof)=.false.
         end if

      end do

      presold=pres
   end do




   ! We have the model on its own vertical levels. Interpolate to observation
   ! points
   mtemp=undef
   msaln=undef
   do iprof=1,nprof

      ! 1) Get valid depth points from argo profile
      ! TODO - Use error flag
      validk=0
      do k=1,nzlev
         !print *,k,iprof
         !print *,cspres(k,iprof)
         if (cspres(k,iprof)/=undef) then
            validk=validk+1
            !print *,validk
            cspres(validk,iprof)=cspres(k,iprof)
            csaln (validk,iprof)=csaln (k,iprof)
            ctemp (validk,iprof)=ctemp (k,iprof)
         end if
      end do

      !print *,'prof, validk ',iprof,validk

      if (validk==0) then
         valid_profile(iprof)=.false.

      else if (valid_profile(iprof)) then

         if (validk<nzlev) then
            cspres(validk+1:nzlev,iprof)=undef
         end if

         call spline_calc_ini(cspres(1:validk,iprof),validk,undef,.true.)

         !print *
         !print *,'cspres ',cspres(1:validk,iprof)
         !print *,'mpres  ',mpres1(:,iprof)

         ! Interpolate model to obs levels
         call spline_calc_1d(mtemp1(:,iprof),mpres1(:,iprof),mtemp(1:validk,iprof),validk)
         call spline_calc_1d(msaln1(:,iprof),mpres1(:,iprof),msaln(1:validk,iprof),validk)

         call spline_calc_release()

         ! Also get discrete model layer values
         mtempd(:,iprof)=undef
         msalnd(:,iprof)=undef
         do k=1,validk
         do k2=1,kdm
            
            if (k2==1) then
               if ( mpres1(k2,iprof)>cspres(k,iprof)) then
                  mtempd(k,iprof)=mtemp1(k2,iprof)
                  msalnd(k,iprof)=msaln1(k2,iprof)
                  mpresd(k,iprof)=mpres1(k2,iprof)*.5
               end if
            else
               if ( mpres1(k2,iprof)>cspres(k,iprof) .and.  mpres1(k2-1,iprof)<cspres(k,iprof)) then
                  mtempd(k,iprof)=mtemp1(k2,iprof)
                  msalnd(k,iprof)=msaln1(k2,iprof)
                  mpresd(k,iprof)=(mpres1(k2,iprof) + mpres1(k2-1,iprof))*.5
               end if
            end if
         end do
         end do

      end if

   end do





   ! Now copy the original file to a new one, then dump model values in there
   if (l_ncdump) then

      call dump_coriolis(trim(corfile),'HYCOM_','PSAL',msaln,nzlev,nprof,undef)
      call dump_coriolis(trim(corfile),'HYCOM_','TEMP',mtemp,nzlev,nprof,undef)

      ! Read levitus data on its own levels 
      call levitus_setup('WOA2001')
      call levitus_interp('temperature',csjuld,cslon,cslat,cspres,levtemp,nprof,nzlev,undef)
      call levitus_interp('salinity',   csjuld,cslon,cslat,cspres,levsaln,nprof,nzlev,undef)
      !modcorfile='TOPAZ_'//trim(corfile)
      !call system('cp '//corfile//' '//modcorfile)
      print *,'WOA2001 salinity '
      call dump_coriolis(trim(corfile),'WOA2001_','PSAL',levsaln,nzlev,nprof,undef)
      print *,'WOA2001 temperature '
      call dump_coriolis(trim(corfile),'WOA2001_','TEMP',levtemp,nzlev,nprof,undef)

      ! Read levitus data on its own levels 
      call levitus_setup('WOA2005')
      call levitus_interp('temperature',csjuld,cslon,cslat,cspres,levtemp,nprof,nzlev,undef)
      call levitus_interp('salinity',   csjuld,cslon,cslat,cspres,levsaln,nprof,nzlev,undef)
      !modcorfile='TOPAZ_'//trim(corfile)
      !call system('cp '//corfile//' '//modcorfile)
      print *,'WOA2005 salinity '
      call dump_coriolis(trim(corfile),'WOA2005_','PSAL',levsaln,nzlev,nprof,undef)
      print *,'WOA2005 temperature '
      call dump_coriolis(trim(corfile),'WOA2005_','TEMP',levtemp,nzlev,nprof,undef)
   end if



   !print *,cslon(2),cslat(2)
   !print *
   !print *,csaln(:,2)
   !print *
   !print *,msaln(:,2)
   !print *
   !print *,msaln1(:,2)

   call system("mkdir -p ArgoCmp")

   where (errorflag_sal ) csaln=undef
   where (errorflag_temp) ctemp=undef
   do iprof=1,nprof
     if (valid_profile(iprof)) then

        write(cii,'(i3.3)') iprof
        write(clon,'(f7.2)') cslon(iprof)     ; clon=adjustl(clon)
        write(clat,'(f7.2)') cslat(iprof)     ; clat=adjustl(clat)
        write(cjuld,'(i10.10)') floor(csjuld(iprof))
        call juliantodate(floor(csjuld(iprof)) ,iyear,imonth,iday,1950,1,1)
        write(cdate,'(i4.4,i2.2,i2.2)') iyear,imonth,iday
        write(cyear,'(i4.4)') iyear
        write(cmonth,'(i2.2)') imonth

        ! Read float id number
        print '(a,a,3i5)','floatid , date is :',floatid(iprof)(:),iyear, imonth, iday
        read(floatid(iprof)(:),*,iostat=ios)   fid_num


        ! This neglects floats which are not all-digits....
        if (ios==0) then 
           write(cfid_num,'(i10.10)') fid_num

           ! KAL -- create new directory for files
           fulldir="ArgoCmp/"//cyear//"/"//cmonth
           call system("[ ! -d "//trim(fulldir)//" ] && mkdir -p "//trim(fulldir))

           !cprof='profile_date'//trim(cdate)//'_id'//trim(floatid(iprof)(:))//'_pos'//trim(clon)//'Ex'//trim(clat)//'N'
           cprof='profile_date'//trim(cdate)//'_id'//cfid_num//'_pos'//trim(clon)//'Ex'//trim(clat)//'N'
           !cprof='ArgoCmp/'//cprof
           cprof=trim(fulldir)//'/'//cprof


           ! KAL - TODO - use error flag here
           open(10,file=trim(cprof))
           do izlev=1,nzlev
           if (ctpres(izlev,iprof)/=undef) then
              write(10,'(f8.1,7f12.5)') ctpres(izlev,iprof),csaln(izlev,iprof),msaln(izlev,iprof), &
                                   ctemp(izlev,iprof),mtemp(izlev,iprof), &
                                   mtempd(izlev,iprof),msalnd(izlev,iprof),mpresd(izlev,iprof)
              
           end if
           end do
           close(10)

           ! For now, intended cleanup of code, later on the above will dump
           ! directly here
           !call system("[ ! -d "//trim(fulldir)//" ] && mkdir -p "//trim(fulldir))
           !call system("cp "//trim(cprof)//" "//trim(fulldir))
           

           

        else
           print *,'No ArgoCmp data produced for this profile -- argocmp does not handle non-digit profile ids (yet)'
        end if
     end if
  end do

   ! Dump a profile legend file as well
   open(10,file='ArgoCmp/profile_legend')
   write(10,*) 'Files contain the following data columnes for each profile'
   write(10,*) 'Column ' , ' Quantity'
   write(10,*) '1)     ' , ' Argo  pressure'
   write(10,*) '2)     ' , ' Argo  Salinity'
   write(10,*) '3)     ' , ' Model Salinity (in Argo point)'
   write(10,*) '4)     ' , ' Argo  Temperature'
   write(10,*) '5)     ' , ' Model Temperature (in Argo point)'
   write(10,*) '6)     ' , ' Model Temperature (in model layer center)'
   write(10,*) '7)     ' , ' Model Salinity    (in model layer center)'
   write(10,*) '8)     ' , ' Model Pressure    (in model layer center)'
   write(10,*)
   write(10,*) 'Filename Legend:'
   write(10,*) 'profile_date[date]_id[float id]_pos[position]'
   write(10,*) '[date] is  in format year month day of month (yyyymmdd)'
   write(10,*) '[float id] is float identifier (WMO float id)'
   write(10,*) 'position is in decimal degrees east x decimal degrees North'
   close(10)
   print *,'Profiles dumped  in directory "ArgoCmp"'




end program p_argocmp
