      module m_NOR05_sed_restart

      use mod_necessary_ecovars
          
      implicit none
c --- Shield everything not explicitly declared public
      private     

      public ::  save_restart_sed, read_restart_sed
                
      contains
      

C --- -----------------------------------------------------------
C --- Restart logic moved from hycom and mod_restart in here
C --- -----------------------------------------------------------
      subroutine read_restart_sed(rt)
      use mod_year_info, only : year_info
      use mod_xc
      implicit none
      type(year_info),  intent(in) :: rt
      character(len=80) :: tmpchar
      integer :: find
      logical :: exa, exb
      if (mnproc==1) print*, 'Sedflag: ', sedflg 
      if (sedflg==1) then
         tmpchar=filesed(rt)
         inquire(exist=exa,file=trim(tmpchar)//'.a')
         inquire(exist=exb,file=trim(tmpchar)//'.b')
         if (.not.exa.and..not.exb) then
           if (mnproc==1) write(lp,'(a)') 'no sediment restart'
           call xcstop('read_restart_sed')
         end if
c
         if (mnproc==1) write(lp,'(a)') 'reading:'//trim(tmpchar)//
     &                                      '.[ab]'
         call read_sed(trim(tmpchar))
      end if
      end subroutine

C --- -----------------------------------------------------------
C --- Restart logic moved from hycom in here, the sediment is 
C --- saved in separate file
C --- -----------------------------------------------------------
      subroutine save_restart_sed(rt)
      use mod_year_info, only : year_info
      implicit none
      type(year_info),  intent(in) :: rt
      character(len=80) :: tmpchar
      integer :: find
      if (sedflg==1) then
         if (mnproc==1) print *,'Dumping sed fields'
         tmpchar=filesed(rt)
         if (mnproc==1) 
     &      write(lp,'(a)',advance='no') 'save_restart_ave: saving:'//
     &      trim(tmpchar)//'.[ab]'
         call save_sed(trim(tmpchar))
         if (mnproc==1)
     &    write(lp,'(a,i5,a)') '     SED-file done'
      endif
      end subroutine

C --- -----------------------------------------------------------
C --- -----------------------------------------------------------
C --- read_sed -- Read sediments and distribute it to nodes
C --- Result is put in ave_week, private to this module
C --- -----------------------------------------------------------
C --- -----------------------------------------------------------
      subroutine read_sed(filebase)
         use mod_za
         implicit none
         character(len=*), intent(in) :: filebase
c
         real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmpfld
         real :: xmin,xmax,xmin2,xmax2,coord
         integer :: iyy, imm, nop,dummyint,tmpcount
         integer :: l_kdm, l_itdm, l_jtdm,ifld,ios,klevel
         character(len=8) a8
         character(len=80) a80
         logical :: exa,exb

         include 'common_blocks.h'

         ! Open files .. 
         inquire(exist=exa,file=filebase//'.a')
         inquire(exist=exb,file=filebase//'.b')
         if (.not.(exa.and.exb)) then
            if (mnproc==1) 
     &      write(lp,*) 'Can not find '//trim(filebase)//'.b'
            call xcstop('(read_sed)')
            stop '(read_sed)'
         end if


         nop=99
         open (unit=nop,file=filebase//'.b',
     &         status='old')
         call zaiopf(filebase//'.a','old',nop)
         
         ! Header
!         read(nop,116) ctitle,iversn,iexpt,yrflag,l_itdm,l_jtdm, 
!     &      l_kdm,imm,iyy,tmpcount
         read(nop,'(a80)') a80;! print *,a80 !Header
         read(nop,'(a80)') a80;! print *,a80 !Header
         read(nop,'(a80)') a80;! print *,a80!Header
         read(nop,'(a80)') a80;! print *,a80!Header
         read(nop,'(a80)') a80;! print *,a80!iversn
         read(nop,'(a80)') a80;! print *,a80!iexpt
         read(nop,'(a80)') a80;! print *,a80!yrflag
         call sedini(l_itdm   ,'idm   ')
         call sedini(l_jtdm   ,'jdm   ')
         call sedini(l_kdm    ,'kdm   ')
         read(nop,*) a80!header

         ! Check dims....
         if (l_kdm/=kdm .or. l_itdm/=itdm .or. l_jtdm /= jtdm) then
            if (mnproc==1) then
               write(lp,'(a)') 'read_sed:dimension mismatch....'
               print *,itdm,l_itdm
               print *,jtdm,l_jtdm
               print *,kdm ,l_kdm
            end if
            call xcstop('(read_sed)')
            stop '(read_sed)'
         end if

         ! Start reading -- for each field check layer and var name
         ios=0
         do while (ios==0)
            read (nop,117,iostat=ios) a8,nstep,time,klevel,coord,
     &                                xmin2,xmax2
Cdiag       print *,a8,time,klevel

            if (ios/=0) exit

            call zaiord(tmpfld(1-nbdy,1-nbdy),ip,.false.,
     &          xmin,xmax, nop)

            ! Check xmin/xmax
            if (abs(xmin-xmin2)>abs(xmin)*1e-4 .or.
     &          abs(xmax-xmax2)>abs(xmax)*1e-4) then
               if     (mnproc.eq.1) then
               write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &           'error - .a and .b files not consistent:',
     &           '.a,.b min = ',xmin,xmin2 ,xmin2-xmin ,
     &           '.a,.b max = ',xmax,xmax2 ,xmax2-xmax
                  print *,'Variable :',a8
               endif
               call xcstop('(read_sed)')
               stop '(read_sed)'
            end if

            if ( trim(a8) == 'sednit' .and. klevel/=0) then
               bot_layer(:,:,1,init)=tmpfld
               bot_layer(:,:,2,init)=tmpfld
            elseif ( trim(a8) == 'sedpho' .and. klevel/=0) then
               bot_layer(:,:,1,ipho)=tmpfld
               bot_layer(:,:,2,ipho)=tmpfld
            elseif ( trim(a8) == 'sedsil' .and. klevel/=0) then
               bot_layer(:,:,1,isil)=tmpfld
               bot_layer(:,:,2,isil)=tmpfld
            elseif ( trim(a8) == 'seddet' .and. klevel/=0) then
               bot_layer(:,:,1,idet)=tmpfld
               bot_layer(:,:,2,idet)=tmpfld
            elseif ( trim(a8) == 'sedsis' .and. klevel/=0) then
               bot_layer(:,:,1,isis)=tmpfld
               bot_layer(:,:,2,isis)=tmpfld
            elseif ( trim(a8) == 'sedfla' .and. klevel/=0) then
               bot_layer(:,:,1,ifla)=tmpfld
               bot_layer(:,:,2,ifla)=tmpfld
            elseif ( trim(a8) == 'seddia' .and. klevel/=0) then
               bot_layer(:,:,1,idia)=tmpfld
               bot_layer(:,:,2,idia)=tmpfld
            elseif ( trim(a8) == 'sedoxy' .and. klevel/=0) then
               bot_layer(:,:,1,ioxy)=tmpfld
               bot_layer(:,:,2,ioxy)=tmpfld
#if defined (YELLOWSED)
            elseif ( trim(a8) == 'sedspm' .and. klevel/=0) then
               bot_layer(:,:,1,ised)=tmpfld
               bot_layer(:,:,2,ised)=tmpfld
            elseif ( trim(a8) == 'sedyel' .and. klevel/=0) then
               bot_layer(:,:,1,iyel)=tmpfld
               bot_layer(:,:,2,iyel)=tmpfld
#endif /*YELLOWSED*/
#ifdef CHATONELLA
            elseif ( trim(a8) == 'sedcha' .and. klevel/=0) then
               bot_layer(:,:,1,icha)=tmpfld
               bot_layer(:,:,2,icha)=tmpfld
#endif
#if defined (DETPHO)
            elseif ( trim(a8) == 'seddetp' .and. klevel/=0) then
               bot_layer(:,:,1,idetp)=tmpfld
               bot_layer(:,:,2,idetp)=tmpfld
#endif
            else
               if (mnproc==1) 
     &            write(lp,*) 'read_sed: Unknown field '//trim(a8)
            end if

         end do

         ! Close unit 
         close(unit=nop)
         call zaiocl(nop)




 116  format (a80/a80/a80/a80/
     & i5,4x/
     & i5,4x/
     & i5,4x/
     & i5,4x/
     & i5,4x/
     & i5,4x/
     & 'field       time step  model day',
     & '  k  dens        min              max')
 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)
      end subroutine

C --- -----------------------------------------------------------
C --- -----------------------------------------------------------
C --- save_sed -- Save sediments and collect it from nodes
C --- -----------------------------------------------------------
C --- -----------------------------------------------------------
      ! Gather sediments and save it to file
      ! NB -- file saving mechanism is now similar
      ! to the one used in standard hycom...
      subroutine save_sed(filebase)
         use mod_za
         implicit none
         character(len=*), intent(in) :: filebase
c
         real :: coord,xmin,xmax
         integer :: nop,kbio
         character(len=8) :: char8

         include 'common_blocks.h'

         ! Open files .. 
         nop=710
         if (mnproc==1) 
     &      open (unit=nop,file=filebase//'.b',
     &            status='replace')
         call zaiopf(filebase//'.a','replace',nop)
         
         ! Header
         if (mnproc==1) 
     &      write(nop,116) ctitle,iversn,iexpt,yrflag,itdm,jtdm, 
     &         kdm

! --     Sediment fields
         coord=0.0
#if defined (ZOOPL)
#if defined (DETPHO)
         do kbio=1,nbio-2
#else
         do kbio=1,nbio-3
#endif
#else
         do kbio=1,nbio
#endif
           ! sediment nitrate
           if (kbio==init) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='sednit  '
           ! sediment phophate
           elseif (kbio==ipho) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='sedpho  '
           ! sediment silicate
           elseif (kbio==isil) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='sedsil  '
           ! sediment detritus
           elseif (kbio==idet) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='seddet  '
           ! sediment silicious detritus
           elseif (kbio==isis) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='sedsis  '
           ! sediment flagellates
           elseif (kbio==ifla) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='sedfla  '
           ! sediment diatoms
           elseif (kbio==idia) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='seddia  '
           ! sediment oxygen
           elseif (kbio==ioxy) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='sedoxy  '
           ! sediment spm
           elseif (kbio==ised) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='sedspm  '
           ! sediment yellow matter
           elseif (kbio==iyel) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='sedyel  '
           elseif (kbio==icha) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='sedcha  '
#if defined (DETPHO)
           elseif (kbio==idetp) then
             call zaiowr(bot_layer(1-nbdy,1-nbdy,1,kbio),ip,.false.,
     &           xmin,xmax, nop, .false.)
             char8='seddetp '
#endif
           endif 
           if     (mnproc.eq.1) then
             write (nop,117) char8,nstep,time,kdm,coord,xmin,xmax
             call flush(nop)
           endif !1st tile
         end do

         ! Close unit 
         close(unit=nop)
         call zaiocl(nop)

 116  format (a80/a80/a80/a80/
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''iexpt '' = experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & i5,4x,'''kdm   '' = Vertical     array size'/
     & 'field       time step  model day',
     & '  k  dens        min              max')
 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)

      end subroutine

c*************************************************
      function filesed(rt)
         use mod_year_info, only : year_info
         use mod_hycom_nersc
         implicit none
         character(len=24) filesed
         type(year_info) rt
         filesed(1:3)=rungen(1:3)
         filesed(4:10)='restart'
         filesed(11:14)=rt%cyy(1:4)
         filesed(15:15)='_'
         filesed(16:18)=rt%cdd(1:3)
         filesed(19:19)='_'
         filesed(20:21)=rt%chh(1:2)
         filesed(22:24)='SED'
      end function filesed

      end module m_NOR05_sed_restart 

c*************************************************
      subroutine sedini(ivar,cvar)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value
c
      character*6 cvarin
c
      read(99,*) ivar,cvarin
      if (mnproc.eq.1) then
      write(lp,6000) cvarin,ivar
      call flush(lp)
      endif !1st tile
c
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in sedini - input ',cvarin,
     +                       ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(sedini)')
               stop
      endif
      return
 6000 format(a6,' =',i10)
      end subroutine sedini
