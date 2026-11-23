      module m_nestcondread
      contains

C --- Reads one nesting record, and shifts the others

#if defined (NEST_OFFLINE)
      subroutine nestcondread(rtin,iens,nestposition)
#else
      subroutine nestcondread(rtin,iens)
#endif
      use mod_xc
#if defined (ICE_NEST) && defined (ICE)
      use mod_common_ice
#endif
      use mod_nesting
      use mod_year_info
      use m_year_day
      use m_filenesting
      implicit none

      ! Input variables
      type(year_info), intent(in) :: rtin
      integer        , intent(in) :: iens
#if defined (NEST_OFFLINE)
      character(len=*), intent(in) :: nestposition
      logical :: warnincons=.true.
      real*8 :: rtime,rntime
#endif

      type(year_info) :: rt
      integer irec,irec_offset,j,k,i,klevel
      character(len=40) nestname
      character(len=100) nestingfilei1,nestingfileii,
     &   nestingfilej1,nestingfilejj,nestingfilehdr,tmpfile
      character(len=5) char5
      character(len=3) char3
      logical exi1,exii,exj1,exjj,exhdr,exall
      real  , dimension(nx,ny) :: globalfld,gdepths
      real*4, dimension(nx,ny) :: iofld
      integer :: ntecsave,irec_dummy, irec_offset_dummy,
     &           iyy, idd, ihh, kdm_dummy, irec_offset_dummy_old
      integer ::ios10,ios11,ios12,ios13,ios16
      integer ::iitdm,jjtdm
      integer :: nvar_tot_file
      real    :: hmin, hmax
      logical :: first, gotfld(12,0:kdm),inputerr
      character(len=80) :: c80
      integer :: indxnest
      logical, parameter :: diagtec=.false.


      include 'common_blocks.h'

      gotfld=.false.

      ! Initalize variable reading time
      ! This may change depending on nesting option. Old
      ! scheme uses unchanged rt
      rt=rtin

#if defined(NEST_OFFLINE)
      ! Check if we should do anything or return
      if (nest_option==2 .or. nest_option==1) then
         if (mnproc==1) write(lp,'(a,i1,a)') 'nest_option=', 
     &      nest_option , ' in action'

         if (nest_option==1) then
            rtime = rt%iyy + (rt%idd+rt%ihh/24d0+rt%iss/86400d0)/
     &         rt%daysinyear
         else if (nest_option==2) then
            rtime = rt%idd+rt%ihh/24d0+rt%iss/86400d0
         end if

         if (nest_slot_index(1) < 0 .or. nest_slot_index(2)< 0. ) then ! restart - always read
            if (mnproc==1) print *,'Reading after  restart ...'
         ! Check if julian day > next item in nesting time slot
         elseif (nest_slot_time(2) < rtime) then

            ! modulus for nesting option 2 - if index is last item
            ! we allow the anove condition - otherwise ... 
            if (nest_option==2) then
               if (nest_slot_index(1) == nrnestfiles .and. 
     &            rtime > nest_slot_time(1)) then
                  return
               end if
            end if

            if (mnproc==1) print *,'Regular nesting read'
         else
            if (mnproc==1) write(lp,*) 'No read required'
            return
         end if
      end if




      ! option=1 use nesting filelis2
      if (nest_option==2 .or. nest_option==1) then

         if (mnproc==1) write(lp,'(a,i1,a)') 'nest_option=', 
     &      nest_option , ' in action'

         ! Find closest match  - climatology so its modulus-like
         if (nestposition=='before')  then ! fill nesting index 2 with fld before "rt"

            if (nest_option==2) then
               indxnest=nrnestfiles
            elseif (nest_option==1) then
               indxnest=-1
            end if

            do i=1,nrnestfiles
               if (nest_option==2) then
                  rntime=nestday(i)+.5
               else if (nest_option==1) then
                  call year_day(real(nestday(i),kind=8),
     &               nestyear(indxnest),rt,'ecmwf')
                  rntime =(nestday(i)+.5)/rt%daysinyear
     &                   +nestyear(i)
               end if

               if (rntime<=rtime) indxnest=i
               !print *,i,nrnestfiles,nestday(i),nestyear(i)
            end do

            if (indxnest==-1) then
               write(lp,*) 'Found no nesting files before date !!'
               call xcstop('nestcondread')
               stop '(nestcondread)'
            end if

            call year_day(real(nestday(indxnest),kind=8),
     &         nestyear(indxnest),rt,'ecmwf')



         else if (nestposition=='after')  then ! fill nesting index 2 with fld after "rt"

            if (nest_option==2) then
               indxnest=1
            elseif (nest_option==1) then
               indxnest=-1
            end if

            do i=nrnestfiles,1,-1

               if (nest_option==2) then
                  rntime=nestday(i)+.5
               else if (nest_option==1) then
                  call year_day(real(nestday(i),kind=8),
     &               nestyear(i),rt,'ecmwf')
                  rntime =(nestday(i)+.5)/rt%daysinyear
     &                   +nestyear(i)
               end if
               print *,i,rtime,rntime,indxnest
               if (rntime> rtime) indxnest=i
            end do

            if (indxnest==-1) then
               write(lp,*) 'Found no nesting files after date !!'
               call xcstop('nestcondread')
               stop '(nestcondread)'
            end if

            call year_day(real(nestday(indxnest),kind=8),
     &         nestyear(indxnest),rt,'ecmwf')

         else
            if (mnproc==1) then
               write(lp,*) 'Unknown nestposition '//trim(nestposition)
            end if
            call xcstop('nestcondread')
            stop '(nestcondread)'
         end if
         if (mnproc==1) then
            print *,'READING NEST OFFLINE OPTION 2'
            print *,'*****************************'
            print *,'rt  :',rt%iyy,rt%idd
            print *,'slot position "'//trim(nestposition)//'"'
            print *,'rtin :',rtin%iyy,rtin%idd
            print *
            print *,indxnest,nestyear(indxnest),nestday(indxnest)
            print *,'NB! offline nesting - check rforce'
            print *,'NB! offline nesting - check rforce'
         end if
      else if (nest_option/=0) then
         if (mnproc==1) then
            write(lp,*) 'Unknown nesting option ',nest_option
         end if
         call xcstop('nestcondread')
         stop '(nestcondread)'
      end if
#endif




      ! Open header file
      nestname=''
      nestingfilehdr=
     &   'Nest/'//trim(filenesting(rt,trim(nestname)))//'.hdr'
      if (mnproc==1) then
         write(lp,'(a)') 'Reading from files .. header is '//
     &      trim(nestingfilehdr)
      end if
      inquire(file=trim(nestingfilehdr),exist=exhdr)
      if (.not.exhdr) then
         if (mnproc==1)  then
            write(lp,*) 'nestcondread: Nest file header '//
     &      ' do not exist'
            print *,trim(nestingfilehdr),exhdr
         end if
         call xcstop('nestcondread')
         stop '(nestcondread)'
      endif
      open(16,file=trim(nestingfilehdr),form='formatted',
     &        access='direct',recl=100,status='old',action='read')


      ! Get number of records per "dump" from header
      irec_offset_dummy    =1
      irec_offset_dummy_old=1
      ios16=0
      first=.true.
      irec=1
      do while (irec_offset_dummy==irec_offset_dummy_old.and.ios16==0)
         irec_offset_dummy_old=irec_offset_dummy

         read(16,'(a80)',rec=irec,iostat=ios16) c80

         !print *,c80 
         if (ios16==0) then
            read(c80,103)
     &         char3,klevel,irec_dummy,irec_offset_dummy,iyy,idd,ihh, 
     &         kdm_dummy,hmin,hmax
            if (first) then
               irec_offset_dummy_old=irec_offset_dummy
               first=.false.
            end if
            !print *,char3,irec,irec_dummy,irec_offset_dummy,ios16
            irec=irec+1
         end if
      end do
      irec=irec-1
      close(16)
      if (ios16/=0.and.irec_offset_dummy/=1) then ! <0 is eof 
         if (mnproc==1) then
            write(lp,*) 'Inconsistency when reading header file'
            write(lp,*) trim(nestingfilehdr)
            print *,irec_dummy
            print *,ios16
         end if
         call xcstop('(nestcondread)')
         stop '(nestcondread)'
      else if (ios16==0) then
         irec_dummy=irec_dummy-1 ! File contains more data
      end if


      nvar_tot_file=irec_dummy
      if (mnproc==1) then
         write(lp,'(a,i4)') 
     &      'Number of records pr dump in input nesting file is :',
     &      nvar_tot_file
      end if
C     call xcstop('(nestcondread testing)')


      ! Offset between dumps - irec is start record for this data dump
#ifdef NEST_OFFLINE
      if (mnproc == 1) print *,'Adjust offset for nesting options !'
#endif
      if (nestdti < 24) then
         !irec_offset=(iens-1)*24/nestdti+rt%ihh/nestdti+1 !iens has never been used
         irec_offset=rt%ihh/nestdti+1
      else
         irec_offset=1
      endif
      irec=(irec_offset-1)*nvar_tot_file+1

      ! Error check on nest dimensions
      if (kknest/=kdm) then
         if (mnproc==1) then
            write(lp,'(a)') 'kknest is not=kdm in dimensions.H'
            write(lp,'(a)') 'Fix this in order to use nesting '
            call flush(lp)
         end if
         call xcstop('(nestcondread)')
         stop '(nestcondread)'
      end if

      ! For testing
      call xcaget(gdepths,depths,0)

      ! Shift nesting - Now we use the HYCOM nest fields...
      pnest(:,:,:,1)=pnest(:,:,:,2)
      unest(:,:,:,1)=unest(:,:,:,2)
      vnest(:,:,:,1)=vnest(:,:,:,2)
      snest(:,:,:,1)=snest(:,:,:,2)
      tnest(:,:,:,1)=tnest(:,:,:,2)

      ubnest(:,:,1)=ubnest(:,:,2)
      vbnest(:,:,1)=vbnest(:,:,2)
      pbnest(:,:,1)=pbnest(:,:,2)

#if defined (ICE_NEST) && defined (ICE)
      uicenest(:,:,1)=uicenest(:,:,2)
      vicenest(:,:,1)=vicenest(:,:,2)
      ficenest(:,:,1)=ficenest(:,:,2)
      hicenest(:,:,1)=hicenest(:,:,2)
#endif

#if defined(NEST_OFFLINE)
      if (nest_option>0) then
         ! Update nesting time slot
         nest_slot_time (1)=nest_slot_time (2)
         nest_slot_index(1)=nest_slot_index(2)
      end if
#endif
      


      ntecsave=1

      !
      iitdm=itdm-inest+1
      jjtdm=jtdm-inest+1

      ! Inquire for nesting files
      exall=.true.
      nestname=''
      !KAL print *,rt%iyy,rt%idd

      nestingfilei1='Nest/'//trim(filenesting(rt,trim(nestname)))//'_i1'
      inquire(file=trim(nestingfilei1),exist=exi1)

      nestingfileii='Nest/'//trim(filenesting(rt,trim(nestname)))//'_ii'
      inquire(file=trim(nestingfileii),exist=exii)

      nestingfilej1='Nest/'//trim(filenesting(rt,trim(nestname)))//'_j1'
      inquire(file=trim(nestingfilej1),exist=exj1)

      nestingfilejj='Nest/'//trim(filenesting(rt,trim(nestname)))//'_jj'
      inquire(file=trim(nestingfilejj),exist=exjj)

      nestingfilehdr=
     &   'Nest/'//trim(filenesting(rt,trim(nestname)))//'.hdr'
      inquire(file=trim(nestingfilehdr),exist=exhdr)

      tmpfile='Nest/'//trim(filenesting(rt,trim(nestname)))
      if (mnproc==1) 
     &  write(lp,'(a,a,a,i4,a,i4)') 'file=',trim(tmpfile)//'_xx',
     &                         ' offset= '     ,irec_offset,
     &                         ' start record=',irec


      exall=exi1.and.exii.and.exj1.and.exjj.and.exhdr
      if (.not.exall) then
         if (mnproc==1)  then
            write(lp,*) 'nestcondread: One or more nest files '//
     &      ' do not exist'
            print *,trim(nestingfilei1 ),exi1
            print *,trim(nestingfileii ),exii
            print *,trim(nestingfilej1 ),exj1
            print *,trim(nestingfilejj ),exjj
            print *,trim(nestingfilehdr),exhdr
         end if
         call xcstop('nestcondread')
         stop '(nestcondread)'
      endif

      ! Open files for io - all nodes.. TODO - one read & distribute
      inquire(iolength=j) iofld(1:inest,1:jtdm)
      open(10,file=trim(nestingfilei1),form='unformatted',
     &   access='direct',recl=j,status='old',action='read')
      !KAL print *,j

      inquire(iolength=j) iofld(1:inest,1:jtdm)
      open(11,file=trim(nestingfileii),form='unformatted',
     &   access='direct',recl=j,status='old',action='read')
      !KAL print *,j

      inquire(iolength=j) iofld(1:itdm,1:inest)
      open(12,file=trim(nestingfilej1),form='unformatted',
     &   access='direct',recl=j,status='old',action='read')

      inquire(iolength=j) iofld(1:itdm,1:inest)
      open(13,file=trim(nestingfilejj),form='unformatted',
     &   access='direct', recl=j,status='old',action='read')

      open(16,file=trim(nestingfilehdr),form='formatted',
     &        access='direct',recl=100,status='old',action='read')

      if (mnproc==1 .and. diagtec) 
     &   open(15,file='nestread.tec',form='formatted', status='replace')



      !Loop through variables
      iofld=0.
      do irec=(irec_offset-1)*nvar_tot_file+1,irec_offset*nvar_tot_file

         read(10,rec=irec,iostat=ios10) 
     &      iofld(1:inest,1:jtdm)

         read(11,rec=irec,iostat=ios11) 
     &      iofld(iitdm:itdm,1:jtdm)

         read(12,rec=irec,iostat=ios12) 
     &      iofld(1:itdm,1:inest)

         read(13,rec=irec,iostat=ios13) 
     &      iofld(1:itdm,jjtdm:jtdm)

         read(16,103,rec=irec,iostat=ios16)
     &      char3,klevel,irec_dummy,irec_offset_dummy,iyy,idd,ihh, 
     &      kdm_dummy,hmin,hmax

         globalfld=iofld

         !KAL print *,irec,char3,klevel
         if (ios10+ios11+ios12+ios13+ios16.ne.0) then
            if (mnproc==1) then
               write(lp,'(a)')'Error occured while reading nest files..'
               print *,'IOSTAT values:'
               print *,trim(nestingfilei1 ),ios10
               print *,trim(nestingfileii ),ios11
               print *,trim(nestingfilej1 ),ios12
               print *,trim(nestingfilejj ),ios13
               print *,trim(nestingfilehdr),ios16
            end if
            call xcstop('nestcondread')
            stop '(nestcondread)'
         end if

#if defined (NEST_OFFLINE)
         if (mnproc==1) then
            print *,'**char3=',char3,' klevel=',klevel,' irec=',irec
         end if
#endif

         !Sanity check
         !if (rt%iyy/=iyy .or. rt%idd/=idd .or. rt%ihh/=ihh .or.
         if (rt%iyy/=iyy .or. rt%ihh/=ihh .or.
     &       irec_dummy/=irec .or. irec_offset/=irec_offset_dummy .or.
     &       kdm_dummy/=kdm) then
#if defined (NEST_OFFLINE)
            if (mnproc==1.and.warnincons) then
               print *,'**Overriding incons check for NESTOFFLINE**'
               warnincons=.false.
            end if
#else
            if (mnproc==1) then
               write(lp,'(a)') 'Nest file header mismatch ...'
               write(lp,'(a,i5,i5)') 'Year :',iyy,rt%iyy
               write(lp,'(a,i5,i5)') 'Day  :',idd,rt%idd
               write(lp,'(a,i5,i5)') 'Hour :',ihh,rt%ihh
               write(lp,'(a,i5,i5)') 'irec :',irec_dummy,irec
               write(lp,'(a,i5,i5)') 'irec_offset :',irec_offset_dummy,
     &            irec_offset
               write(lp,'(a,i5,i5)') 'kdm :',kdm_dummy,kdm
            end if
            call xcstop('nestcondread')
            stop '(nestcondread)'
#endif
         elseif ( rt%idd/=idd ) then
            if (mnproc==1.and.irec==(irec_offset-1)*nvar_tot_file+1)then
               write(lp,'(a)') 'WARNING:'
               write(lp,'(a)') 'Nest file header mismatch for day'
               write(lp,'(a)') 'If you used a link ignore this message'
            end if
         elseif (klevel>kdm .or.klevel<0) then ! 2D vars are 0
            if (mnproc==1) then
               write(lp,'(a)') 'vertical level is out of range for '//
     &                         'nesting data:'
               write(lp,'(a,i3)') 'Variable is '//char3//', level is ',
     &            klevel
            end if
            call xcstop('nestcondread')
            stop '(nestcondread)'
         end if
            



         ! Use header file to place data in nest arrays
         if     (trim(char3)=='SAL') then
            call xcaput(globalfld,snest(:,:,klevel,2),0)
            gotfld(1,klevel)=.true.
         elseif (trim(char3)=='TEM') then ! Temperature
            call xcaput(globalfld,tnest(:,:,klevel,2),0)
            gotfld(2,klevel)=.true.
         elseif (trim(char3)=='INT') then ! Interface values
            call xcaput(globalfld,pnest(:,:,klevel,2),0)
            gotfld(3,klevel)=.true.
         elseif (trim(char3)=='UT') then ! Velocities
            call xcaput(globalfld,unest(:,:,klevel,2),0)
            gotfld(4,klevel)=.true.
         elseif (trim(char3)=='VT') then ! Velocities
            call xcaput(globalfld,vnest(:,:,klevel,2),0)
            gotfld(5,klevel)=.true.
         elseif (trim(char3)=='PB') then ! ! Baro Pressure
            call xcaput(globalfld,pbnest(:,:,2),0)
            gotfld(6,klevel)=.true.
         elseif (trim(char3)=='UB') then ! ! Baro Velocity
            call xcaput(globalfld,ubnest(:,:,2),0)
            gotfld(7,klevel)=.true.
         elseif (trim(char3)=='VB') then ! ! Baro Velocity
            call xcaput(globalfld,vbnest(:,:,2),0)
            gotfld(8,klevel)=.true.
#if defined (ICE_NEST) && defined (ICE) 
         elseif (trim(char3)=='UI') then ! ! Ice Velocity
            call xcaput(globalfld,uicenest(:,:,2),0)
            gotfld(9,klevel)=.true.
         elseif (trim(char3)=='VI') then ! ! Ice Velocity
            call xcaput(globalfld,vicenest(:,:,2),0)
            gotfld(10,klevel)=.true.
         elseif (trim(char3)=='HI') then ! ! Ice Thickness
            call xcaput(globalfld,hicenest(:,:,2),0)
            gotfld(11,klevel)=.true.
         elseif (trim(char3)=='FI') then ! ! Ice Concentration
            call xcaput(globalfld,ficenest(:,:,2),0)
            gotfld(12,klevel)=.true.
#else
         elseif ( trim(char3)=='UI' .or.  trim(char3)=='VI' .or.
     &            trim(char3)=='FI' .or.  trim(char3)=='HI') then ! ! Ice Velocity
           if (mnproc==1) write(lp,*) 'Nest ice variable '//
     &        trim(char3)//' skipped'
#endif
         else
            if (mnproc==1) then
               write(lp,'(a)') 'Unknown Nesting field ID '//char3
            end if
            call xcstop('nestcondread')
            stop '(nestcondread)'
         end if

         !! For testing
         !call xcaget(globalfld,tnest(:,:,klevel,2),0)

         ! Store for later
         if (mnproc==1.and.diagtec) then
         if (trim(char3)=='TEM') then ! Temperature
            write(char5,'(a3,i2.2)') 'TEM',klevel
            if (ntecsave==1) then
               WRITE(15,'(''TITLE= "Read Nested fields"'')')
               write(15,'(a)')'VARIABLES="i" "j" "depths" "variable"'
               WRITE(15,'(a,a5,a,i4,a,i4,a)') 'ZONE T=',trim(char5),
     &               ' I=',itdm,' J=', jtdm, ' F=BLOCK'
               WRITE(15,101)((i,i=1,itdm),j=1,jtdm)
               WRITE(15,101)((j,i=1,itdm),j=1,jtdm)
               WRITE(15,102)((gdepths(i,j),i=1,itdm),j=1,jtdm)
            else
               WRITE(15,'(a,a5,a,i4,a,i4,a)') 'ZONE T=',trim(char5),
     &               ' I=',itdm,' J=', jtdm, ' F=BLOCK'
               write(15,'(a)')'D=(1,2,3)'
            end if
            WRITE(15,102)((globalfld(i,j),i=1,itdm),j=1,jtdm)
            ntecsave=ntecsave+1
         end if
         end if

      end do

      if (mnproc==1 .and. diagtec) close(15)
      close(10)
      close(11)
      close(12)
      close(13)
      close(16)


      ! Test that everything is read
      inputerr=.false.
      if (.not.all(gotfld(1,1:kdm))) then
         if (mnproc==1) write(lp,*) 'Not all salinity    values read'
         inputerr=.true.
      end if

      if (.not.all(gotfld(2,1:kdm))) then
         if (mnproc==1) write(lp,*) 'Not all temperature values read'
         inputerr=.true.
      end if

      if (.not.all(gotfld(3,1:kdm))) then
         if (mnproc==1) write(lp,*) 'Not all interface   values read'
         inputerr=.true.
      end if

      if (.not.all(gotfld(4,1:kdm))) then
         if (mnproc==1) write(lp,*) 'Not all u-velocity  values read'
         inputerr=.true.
      end if

      if (.not.all(gotfld(5,1:kdm)))  then
         if (mnproc==1) write(lp,*) 'Not all v-velocity  values read'
         inputerr=.true.
      end if

      if (.not.gotfld(6,0)) then
         if (mnproc==1) write(lp,*) 'Barotropic pressure   not  read'
         inputerr=.true.
      end if

      if (.not.gotfld(7,0)) then
         if (mnproc==1) write(lp,*) 'Barotropic u-velocity not  read'
         inputerr=.true.
      end if

      if (.not.gotfld(8,0)) then
         if (mnproc==1) write(lp,*) 'Barotropic v-velocity not  read'
         inputerr=.true.
      end if

#if defined (ICE_NEST) && defined (ICE) 
      if (.not.gotfld(9,0)) then
         if (mnproc==1) write(lp,*) 'Ice u-velocity not  read'
         inputerr=.true.
      end if

      if (.not.gotfld(10,0)) then
         if (mnproc==1) write(lp,*) 'Ice v-velocity not  read'
         inputerr=.true.
      end if

      if (.not.gotfld(11,0)) then
         if (mnproc==1) write(lp,*) 'Ice thickness  not  read'
         inputerr=.true.
      end if

      if (.not.gotfld(12,0)) then
         if (mnproc==1) write(lp,*) 'Ice fraction   not  read'
         inputerr=.true.
      end if
#endif

      if (inputerr) then
         if (mnproc==1) write(lp,*) 'Error when reading nesting files.'
         call xcstop('(nestcondread)')
      !else
      !   if (mnproc==1) write(lp,'(a)') 'Nesting files read - all ok'
      end if
C     call xcstop('(nestcondread testing)')

#if defined(NEST_OFFLINE)
      ! Check if we should do anything or return
      if (nest_option==1) then
         if (mnproc==1) write(lp,*) 'nest_option=1 in action'
         call year_day(real(nestday(indxnest),kind=8),
     &      nestyear(indxnest),rt,'ecmwf')
         nest_slot_index(2)=indxnest
         nest_slot_time (2)=(nestday(indxnest)+.5)/rt%daysinyear
     &                     +nestyear(indxnest)
      elseif (nest_option==2) then
         if (mnproc==1) write(lp,*) 'nest_option=2 in action'
         nest_slot_index(2)=indxnest
         nest_slot_time (2)=nestday(indxnest)+.5
      end if
#endif

  101 format(30i5)      
  102 format(10e14.3)
C 103 format(a3," level=",i4," record=",i5," offset=",i5," date=",i4,
C    &       " ",i3," ",i2," kdm=",i3," min/max:",2e12.2)
  103 format(a3,7x,i4,8x,i5,8x,i5,6x,i4,
     &       1x,i3,1x,i2,5x,i3,9x,2e12.2)

      end subroutine nestcondread
      end module m_nestcondread
