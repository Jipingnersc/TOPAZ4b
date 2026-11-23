      program p_slacompare
      ! Routine to compare modeled sla with observations 
      ! Program takes three arguments  - first is model
      ! field (only DODS link for now), second is SLA
      ! field - also DODS. Third is some form of time 
      ! unit (to be decided)
      !
      ! This program must be linked against nc-dap
      use mod_xc, only: idm,jdm, xcspmd
      use mod_za, only: zaiost
      use mod_common
      use mod_sphere_tools
      use m_get_grid 
      use m_duacsread
      use m_getfiletype
      use m_getfileinfo
      use m_read_field_wrapper
      use m_oldtonew
      use m_pivotp
      use m_initconfmap
      use mod_year_info
      use m_datetojulian
      implicit none

      include 'nc-dap.inc'

      integer, parameter :: nblon=1080 ,     nblat=915 
      real, parameter :: undef=-1e14 , noasciidata=1e4
      real :: lon(nblon), lat(nblat), sla(nblat, nblon)

      character(len=80) :: hycfile, fbase
      integer :: filetype,nrmem

      character(len=3) :: creg
      integer, allocatable :: regflag(:,:)
      real, allocatable :: modssh(:,:), modsla(:,:),newsla(:,:)
      real*8, allocatable :: meanssh(:,:)
      character(len=11) tag7
      logical :: ex,regionflag,havedata
      integer :: i,j,ipiv,jpiv,jday,ireg
      real :: lat_n,lon_n
      type(year_info) :: rti, rtd


      integer :: stat
      integer :: ncid, dimid_idm, dimid_jdm, obs_varid, mod_varid, 
     &   mod_msshvarid, mod_sshvarid, maskvarid


      ! Regional stats
      integer, parameter :: maxreg=20, maxcorner=20
      integer :: ios, iosa, iosb, iosc, ierr
      integer :: ncorners
      real, dimension(:), allocatable :: clons, clats
      character(len=80) :: regname
      integer :: npoints, npoints2
      real :: diffsq, meanobs, meanmod, meanmodonly




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

      ! get julian day rel 1950-1-1
      jday=datetojulian(rtd%iyy,rtd%imm,rtd%idm,1950,1,1)



      ! Read duacs data from opendap
      !print *,rtd
      call duacsread(sla,lon,lat,nblon,nblat,undef,jday)
      print *

      ! Get model grid
      call xcspmd()
      call zaiost()
      call get_grid() ; where(depths>1e20) depths=0.
      call initconfmap()

      ! Read model ssh field
      allocate(modssh(idm,jdm))
      call read_field_wrapper_2d(modssh,'ssh     ',fbase,idm,jdm,0,
     &         filetype, nrmem,undef)
      where (depths<.1) modssh=undef

      ! Read mean model ssh
      allocate(meanssh(idm,jdm))
      if (idm>999 .or. jdm>999) then
         write(tag7,'(i5.5,a,i5.5)')idm,'x',jdm
      else
         write(tag7,'(i3.3,a,i3.3)')idm,'x',jdm
      end if
      inquire(exist=ex, file='meanssh'//trim(tag7)//'.uf')
      if (ex) then
         open(11,file  ='meanssh'//trim(tag7)//'.uf', 
     &           status='old',form='unformatted')
         read (11)meanssh
         close(11)
      else
         print *,'Could not find mean ssh file '
         call exit(1)
      end if
      !meanssh=meanssh/100.
      !print *,meanssh
      where (depths<.1) meanssh=undef


      print *,minval(modssh), maxval(modssh)
      print *,minval(meanssh), maxval(meanssh)
      allocate(modsla(idm,jdm))
      modsla=0.
      modsla=modssh-real(meanssh,kind=4)
      where (depths<.1) modsla=undef


      ! Find model points corr  to sla points
      allocate(newsla(idm,jdm))
      newsla=undef
      do j=1,nblon
      do i=1,nblat
      !if (slaisok?)

         call oldtonew(lat(i),lon(j),lat_n,lon_n)
         call pivotp(lon_n,lat_n,ipiv,jpiv)

         if ( ipiv>1 .and. ipiv<idm .and. 
     &        jpiv>1 .and. jpiv<jdm ) then

            if (depths(ipiv,jpiv)>1. ) then
               newsla(ipiv,jpiv)=sla(i,j)   
            end if

         end if

      end do
      end do


      ! Quick dump of data points to netcdf file
      stat = NF_CREATE('slacmp.nc', NF_CLOBBER , ncid)
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if
      stat = NF_DEF_DIM (NCID, 'idm', idm, dimid_idm)
      stat = NF_DEF_DIM (NCID, 'jdm', jdm, dimid_jdm)

      stat = NF_DEF_VAR(NCID, 'slaobs', NF_REAL, 2,
     &   (/dimid_idm,dimid_jdm/), obs_varid)

      stat = NF_DEF_VAR(NCID, 'slamod', NF_REAL, 2,
     &   (/dimid_idm,dimid_jdm/), mod_varid)

      stat = NF_DEF_VAR(NCID, 'modmeanssh', NF_REAL, 2,
     &   (/dimid_idm,dimid_jdm/), mod_msshvarid)

      stat = NF_DEF_VAR(NCID, 'modssh', NF_REAL, 2,
     &   (/dimid_idm,dimid_jdm/), mod_sshvarid)

      stat = NF_PUT_ATT_REAL  (NCID, obs_VARID, '_FillValue',
     &   NF_REAL,1, undef )
      stat = NF_PUT_ATT_REAL  (NCID, mod_VARID, '_FillValue',
     &   NF_REAL, 1, undef )
      stat = NF_PUT_ATT_REAL  (NCID, mod_msshVARID, '_FillValue',
     &   NF_REAL, 1, undef )
      stat = NF_PUT_ATT_REAL  (NCID, mod_sshVARID, '_FillValue',
     &   NF_REAL, 1, undef )

      STAT = NF_ENDDEF(NCID)


      stat = NF_PUT_VARA_REAL(NCID,mod_VARID,(/1,1/),(/idm,jdm/),
     &   modsla)
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if

      stat = NF_PUT_VARA_REAL(NCID,mod_msshVARID,(/1,1/),(/idm,jdm/),
     &   real(meanssh,kind=4))
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if

      stat = NF_PUT_VARA_REAL(NCID,obs_VARID,(/1,1/),(/idm,jdm/),
     &   newsla)
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if


      stat = NF_PUT_VARA_REAL(NCID,mod_sshVARID,(/1,1/),(/idm,jdm/),
     &   modssh)
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if






!!!!!!! Format for regional stats file:
!regname
!ncorners
!clons
!clats
!!!!!!!
! Example: 
!NORTHATLANTIC
!4
!-80 20 20 -80
!40  40 60  60
!!!!!!!
! You can have several region defs in one file, if they follow eachother
      inquire(exist=ex,file='slareduce.in')
      havedata=.not.all(modssh==undef)
      if (ex) then
         allocate(regflag(idm,jdm))
         open(10,file='slareduce.in',action='read')
         ios=0
         ierr=0
         ireg=0
         regflag=0
         do while(ios==0 .and. ierr==0)
            read(10,*,iostat=ios)    regname
            if (ios==0) then
               read(10,*,iostat=iosa)ncorners                     !; print *,ncorners
               allocate(clons(ncorners))
               allocate(clats(ncorners))
               read(10,*,iostat=iosb)clons
               read(10,*,iostat=iosc)clats

               ierr=maxval(abs( (/iosa,iosb,iosc/) ))
               if (ierr==0) then
                  print *,'-------'//trim(regname)//'---------'
                  print *,ncorners
                  print *,clons
                  print *,clats
                  npoints=0
                  npoints2=0
                  diffsq=0.
                  meanobs=0.
                  meanmod=0.
                  meanmodonly=0.
                  ireg=ireg+1
                  do j=1,jdm
                  do i=1,idm
                     regionflag=inbox(clons,clats,ncorners,
     &                          modlon(i,j),modlat(i,j))
                     if (regionflag .and. depths(i,j)>1 .and.
     &                   newsla(i,j)/=undef) then
                        npoints=npoints+1
                        meanobs=meanobs+newsla(i,j)
                        meanmod=meanmod+modsla(i,j)*100.
                        diffsq=diffsq+(newsla(i,j)-modsla(i,j)*100.)**2
                     end if

                     if (regionflag .and. depths(i,j)>1 ) then
                        npoints2=npoints2+1
                        meanmodonly=meanmodonly+modsla(i,j)*100.
                     end if

                     if (regionflag) then
                        regflag(i,j)=regflag(i,j)+2**ireg
                     end if
                  end do
                  end do

                  if (npoints2>0) then
                     meanmodonly=meanmodonly/npoints2
                  end if

                  if (.not.havedata) then
                     open(11,file=trim(regname)//'.asc',
     &                       position='append',
     &                       status='unknown')
                     write(11,'(5f14.5,i8)')
     &                           rtd%iyy+float(rtd%idd)/
     &                                  float(rtd%daysinyear),
     &                           rti%iyy+float(rti%idd)/
     &                                  float(rti%daysinyear),
     &                           noasciidata,
     &                           noasciidata,
     &                           noasciidata,
     &                           0
                     close(11)
                  elseif (npoints>0) then
                     print *,'# points ',npoints
                     print *,'mean obs ',meanobs/npoints
                     print *,'mean mod ',meanmod/npoints
                     print *,'bias     ',(meanmod-meanobs)/npoints
                     print *,'RMS      ',sqrt(diffsq/npoints)

                     open(11,file=trim(regname)//'.asc',
     &                       position='append',
     &                       status='unknown')
                     write(11,'(5f14.5,i8)')
     &                           rtd%iyy+float(rtd%idd)/
     &                                  float(rtd%daysinyear),
     &                           rti%iyy+float(rti%idd)/
     &                                  float(rti%daysinyear),
     &                           meanmodonly,              
     &                           (meanmod-meanobs)/npoints,
     &                           sqrt(diffsq/npoints), 
     &                           npoints
                     close(11)
                  end if

               end if
               deallocate (clons)
               deallocate (clats)
            end if
         end do

         if (ireg>0) then

            ! put mask in netcdf file
            STAT = NF_REDEF(NCID)
            stat = NF_DEF_VAR(NCID, 'regionmask',NF_INT,
     &         2,(/dimid_idm,dimid_jdm/), maskvarid)
            if (stat /= nf_noerr) then
               print *,nf_strerror(stat)
               call exit(3)
            end if
            STAT = NF_ENDDEF(NCID)

            where(depths<.1) regflag=-2
            stat = NF_PUT_VARA_INT(NCID,maskvarid,
     &                (/1,1/),(/idm,jdm/),regflag)
            if (stat /= nf_noerr) then
               print *,nf_strerror(stat)
               call exit(3)
            end if
         end if

      end if

      stat = NF_CLOSE(ncid)





      end program

