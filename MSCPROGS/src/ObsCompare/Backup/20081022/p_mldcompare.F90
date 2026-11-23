      program p_mldcompare
      ! Routine to compare modeled sla with observations 
      ! Program takes three arguments  - first is model
      ! field (only DODS link for now), second is SLA
      ! field - also DODS. Third is some form of time 
      ! unit (to be decided)
      !
      ! This program must be linked against netcdf
      use mod_xc, only: idm,jdm, xcspmd
      use mod_za, only: zaiost
      use m_dbm_mld_read
      use mod_common
      use mod_sphere_tools
      use m_get_grid 
      use m_getfiletype
      use m_getfileinfo
      use m_read_field_wrapper
      use m_oldtonew
      use m_pivotp
      use m_initconfmap
      use mod_year_info
      use m_datetojulian
      use netcdf
      use m_parse_blkdat
      use mod_regions
      implicit none

      real, parameter :: undef=-1e14

      character(len=80) :: hycfile, fbase
      integer :: filetype,nrmem

      integer       :: nblon, nblat
      real, pointer :: obsmld(:,:),mldlon(:),mldlat(:)


      character(len=3) :: creg
      integer, allocatable :: regflag(:,:)
      real, allocatable :: modmld(:,:), temp(:,:), dp(:,:), dpsum(:,:), nmodmld(:,:) , temp1(:,:)
      character(len=11) tag7
      logical :: ex,regionflag,havedata
      integer :: i,j,ipiv,jpiv,jday,ireg,diy
      real :: lat_n,lon_n
      type(year_info) :: rti, rtd


      integer :: stat,k
      integer :: ncid, dimid_idm, dimid_jdm, obs_varid, mod_varid, maskvarid
      logical, allocatable :: mask(:,:)
      real :: rdummy


      ! Regional stats
      integer :: ios, iosa, iosb, iosc, ierr
      character(len=80) :: regname
      integer :: numregion
      integer :: npoints, npoints2
      real :: diffsq, meanobs, meanmod, RMS, bias, yfrac




      ! Process input args
      if (iargc()/=1) then
          print *,'Usage:'
          print *,'  slacmp hycom-file'
          stop
       else
          call getarg(1,hycfile)
          !print *,hycfile
       end if

      ! Deduce hycom filetype
      call getfiletype(hycfile,fbase,filetype)

      ! Deduce data time
      call getfileinfo(fbase,filetype,rti,rtd,nrmem)

      ! get julian day rel this year-1-1
      jday=datetojulian(rtd%iyy,rtd%imm,rtd%idm,rtd%iyy,1,1)
      diy =datetojulian(rtd%iyy,12,31,rtd%iyy,1,1)

      ! "Float year"
      yfrac=rtd%iyy + float(jday)/float(diy)



      ! Read duacs data from opendap
      call dBM_MLD_read(obsmld,mldlon,mldlat,nblon,nblat,undef,jday)
      print *

      ! Get model grid
      call xcspmd()
      call zaiost()
      call get_grid() ; where(depths>1e20) depths=0.
      call parse_blkdat('kdm   ','integer',rdummy,kdm)
      call initconfmap()

      ! Read model temp field
      allocate(nmodmld(nblon,nblat))
      allocate(modmld(idm,jdm))
      allocate(temp  (idm,jdm))
      allocate(temp1 (idm,jdm))
      allocate(dp    (idm,jdm))
      allocate(dpsum (idm,jdm))

      ! A bit complicated because we must cycle through layers to get MLD
      k=1
      modmld=0.
      temp=0.
      dpsum=0.
      do while (k<=kdm)
         call read_field_wrapper_2d(dp    ,'pres    ',fbase,idm,jdm,k, &
                  filetype, nrmem,undef)
         call read_field_wrapper_2d(temp  ,'temp    ',fbase,idm,jdm,k, &
                  filetype, nrmem,undef,dplayer=dp)
         if (k==1) temp1=temp
         if (filetype/=4) then
            dpsum=dp
         else
            dpsum=dpsum+dp
         end if
         where(depths>.1 .and. abs(temp-temp1)>0.2 .and. modmld==0.)
            modmld=dpsum
         end where
         k=k+1
      end do

      where(depths>.1 .and. modmld==0.)
         modmld=dpsum
      end where


         
        




      ! Find model points corr to mld points
      nmodmld=undef
      do i=1,nblon
      do j=1,nblat
         call oldtonew(mldlat(j),mldlon(i),lat_n,lon_n)
         call pivotp(lon_n,lat_n,ipiv,jpiv)
         if ( ipiv>1 .and. ipiv<idm .and.  jpiv>1 .and. jpiv<jdm ) then
            if (depths(ipiv,jpiv)>1. ) then
               nmodmld(i,j)=modmld(ipiv,jpiv)   
            end if
         end if
      end do
      end do

      ! Correction factor for pressure coords
      if (maxval(nmodmld,mask=nmodmld>undef) > 10.*9806) then
         where(nmodmld>undef) nmodmld=nmodmld/9806.
      end if


      ! Quick dump of data points to netcdf file
      stat = NF90_CREATE('mldcmp.nc', NF90_CLOBBER , ncid)
      if (stat /= nf90_noerr) then
         print *,nf90_strerror(stat)
         call exit(3)
      end if
      stat = NF90_DEF_DIM (NCID, 'nblon', nblon, dimid_idm)
      stat = NF90_DEF_DIM (NCID, 'nblat', nblat, dimid_jdm)
      stat = NF90_DEF_VAR(NCID, 'mldobs', NF90_Float, (/dimid_idm,dimid_jdm/), obs_varid)
      stat = NF90_DEF_VAR(NCID, 'mldmod', NF90_Float, (/dimid_idm,dimid_jdm/), mod_varid)
      stat = NF90_PUT_ATT(NCID, obs_varid, '_FillValue', undef )
      stat = NF90_PUT_ATT(NCID, mod_varid, '_FillValue', undef )
      STAT = NF90_ENDDEF(NCID)
      stat = NF90_PUT_VAR(NCID,mod_VARID,nmodmld,(/1,1/),(/nblon,nblat/))
      if (stat /= nf90_noerr) then
         print *,nf90_strerror(stat)
         call exit(3)
      end if
      stat = NF90_PUT_VAR(NCID,obs_varid,obsmld,(/1,1/),(/nblon,nblat/))
      if (stat /= nf90_noerr) then
         print *,nf90_strerror(stat)
         call exit(3)
      end if


      ! Initialize regions for RMS calc
      call read_regions()
      numregion=getnregions()
      print *,numregion
      allocate(mask   (nblon,nblat))
      allocate(regflag(nblon,nblat))
      do ireg=1,numregion

         print *,'iregion: ',ireg

         ! Get mask for this region
         call getmask(ireg,mask,mldlon,mldlat,nblon,nblat)

         where (mask .and. obsmld/=undef .and. nmodmld/=undef)
            mask=.true.
         elsewhere
            mask=.false.
         end where
         print *,count(mask)

         meanmod=0.
         meanobs=0.
         diffsq=0.
         npoints=0
         do j=1,nblat
         do i=1,nblon
         if (mask(i,j)) then
            npoints=npoints+1
            meanobs=meanobs+obsmld(i,j)
            meanmod=meanmod+nmodmld(i,j)
            diffsq=diffsq+(obsmld(i,j)-nmodmld(i,j))**2
            regflag(i,j)=regflag(i,j)+2**ireg ! Can be binary decoded
         end if
         end do
         end do

         bias=(meanmod-meanobs)/npoints
         RMS = sqrt( diffsq/npoints - bias**2)
         meanmod=meanmod/npoints
         meanobs=meanobs/npoints
         call dumprms_ascii(ireg,yfrac,meanmod,meanobs,bias,rms,npoints,'mldRMS_')
         print *,'meanmod:',meanmod
         print *,'meanobs:',meanobs
         print *,'bias   :',bias   
         print *,'RMS    :',RMS    
      end do


      if (numregion>0) then

         ! put mask in netcdf file
         STAT = NF90_REDEF(NCID)
         stat = NF90_DEF_VAR(NCID, 'regionmask',NF90_INT,(/dimid_idm,dimid_jdm/), maskvarid)
         if (stat /= nf90_noerr) then
            print *,nf90_strerror(stat)
            call exit(3)
         end if
         STAT = NF90_ENDDEF(NCID)

         stat = NF90_PUT_VAR(NCID,maskvarid,regflag,(/1,1/),(/nblon,nblat/))
         if (stat /= nf90_noerr) then
            print *,nf90_strerror(stat)
            call exit(3)
         end if
      end if


      stat = NF90_CLOSE(ncid)





      end program

