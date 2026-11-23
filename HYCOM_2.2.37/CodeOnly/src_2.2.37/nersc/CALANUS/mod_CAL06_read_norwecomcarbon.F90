module mod_CAL06_read_norwecomcarbon
   use mod_year_info
   use mod_xc
   type daily_averages
      integer         :: nrmem
      integer         :: lmem(1000) ! To decide which members are saved
      type(year_info) :: rtinit
      type(year_info) :: rtdump
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
           p, dia, fla
   end type daily_averages

   contains

   subroutine read_daily_carbon(rt)
      use mod_za
      use mod_year_info
      use mod_CAL06_calnusparticles,only : mod_carb, pm_rid
      implicit none
      include 'comcal.h'
      type(daily_averages) :: da
      type(year_info), intent(in) :: rt

      real :: xmin,xmax,xmin2,xmax2,coord,loctime
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmpfld
      integer :: idummy, l_itdm,l_jtdm,l_kdm,k,nop, ios,n, &
                locnstep, lociversn, lociexpt, locyrflag
      character(len=80) :: locctitle(4)
      character(len=60) :: filenamea,filenameb
      character(len=8) :: char8
      logical :: ex

      real :: mod_carb_ave(kdm)
      integer :: mod_carb_num(kdm)
      integer :: i,j

!AS   check if the small carbon files exists and read this if it does

      inquire(exist=ex,file='NORWECOM_carbon/Carbon_'//rt%cyy//'_'//rt%cdd//'.b')
!      print*, 'NORWECOM_carbon/Carbon_'//rt%cyy//'_'//rt%cdd//'.b', ex
      if (ex) then
         filenamea=trim('NORWECOM_carbon/Carbon_'//rt%cyy//'_'//rt%cdd//'.a')
         filenameb=trim('NORWECOM_carbon/Carbon_'//rt%cyy//'_'//rt%cdd//'.b')

         if (mnproc==1) then
            write(lp,'(a)') 'Reading '// &
            trim(da_filebase(pm_rid,rt%cyy,rt%cdd))//'.[ab]'
            call flush(lp)
         end if
         ! Open files .. 
         nop=710
         open (unit=nop,file=trim(filenameb),&
            status='old')
!         print*, filenamea
!         print*, filenameb

         call zaiopf(trim(filenamea),'old',nop)

         ! Read .b Header
         read(nop,116) locctitle,lociversn,lociexpt,locyrflag, &
           l_itdm,l_jtdm, &
           l_kdm,da%rtinit%iyy,da%rtinit%idd, &
           da%rtdump%iyy,da%rtdump%idd,da%nrmem

         if (l_itdm/=itdm.and.l_jtdm/=jtdm.and.l_kdm/=kdm) then
            Print *,'Grid dimensions are wrong!!!'
            call xcstop ('(daily_read_ave)')
            stop         '(daily_read_ave)'
         end if

        ! Start reading fields
         ios=0
         do while (ios==0)
            read (nop,117,iostat=ios) char8,locnstep,loctime, &
                                   k,coord,xmin2,xmax2
!            print *,char8,locnstep,loctime,k,ios

            if (ios==0) then

            call zaiord(tmpfld(1-nbdy,1-nbdy),ip,.false., &
               xmin,xmax, nop)

            ! Check xmin/xmax
            if (abs(xmin-xmin2)>abs(xmin)*1e-4 .or.  &
                abs(xmax-xmax2)>abs(xmax)*1e-4) then
               if     (mnproc.eq.1) then
               write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')  &
                'error - .a and .b files not consistent:',   &
                '.a,.b min = ',xmin,xmin2 ,xmin2-xmin ,      &
                '.a,.b max = ',xmax,xmax2 ,xmax2-xmax
                  print *,'Variable :',char8
               endif
               call xcstop('(mod_CAL06_read_norwecomcarbon)')
               stop '(mod_CAL06_read_norwecomcarbon)'
            end if

            ! Put fields into its place:
            if ( trim(char8)=='davp' .and. k /= 0 ) then
               da%p(:,:,k)=tmpfld
            elseif ( trim(char8)=='davc' .and. k /= 0 ) then
               mod_carb(:,:,k)=tmpfld
            endif

           end if  !ios==0
         enddo

         call zaiocl(nop)
         if (mnproc==1) close(nop)

      else
!!! Fix her nr du faar filnavn
         filenamea='NORWECOM_carbon/'//trim(da_filebase(pm_rid,rt%cyy,rt%cdd))//'.a'
         filenameb='NORWECOM_carbon/'//trim(da_filebase(pm_rid,rt%cyy,rt%cdd))//'.b'

         if (mnproc==1) then
            write(lp,'(a)') 'Reading '// &
            trim(da_filebase(pm_rid,rt%cyy,rt%cdd))//'.[ab]'
            call flush(lp)
         end if


         ! Open files .. 
         nop=710
         open (unit=nop,file=trim(filenameb),&
            status='old')
!         print*, filenamea
!         print*, filenameb

         call zaiopf(trim(filenamea),'old',nop)

         ! Read .b Header
         read(nop,116) locctitle,lociversn,lociexpt,locyrflag, &
           l_itdm,l_jtdm, &
           l_kdm,da%rtinit%iyy,da%rtinit%idd, &
           da%rtdump%iyy,da%rtdump%idd,da%nrmem

         if (l_itdm/=itdm.and.l_jtdm/=jtdm.and.l_kdm/=kdm) then
            Print *,'Grid dimensions are wrong!!!'
            call xcstop ('(daily_read_ave)')
            stop         '(daily_read_ave)'
         end if

        ! Read the local ensemble members into ".b" header
         do n=1,1000
            read(nop,118) idummy,da%lmem(n)
         end do
         read(nop,119)

         ! Start reading fields
         ios=0
         do while (ios==0)
            read (nop,117,iostat=ios) char8,locnstep,loctime, &
                                   k,coord,xmin2,xmax2
!            print *,char8,locnstep,loctime,k,ios

            if (ios==0) then

            call zaiord(tmpfld(1-nbdy,1-nbdy),ip,.false., &
               xmin,xmax, nop)

            ! Check xmin/xmax
            if (abs(xmin-xmin2)>abs(xmin)*1e-4 .or.  &
                abs(xmax-xmax2)>abs(xmax)*1e-4) then
               if     (mnproc.eq.1) then
               write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')  &
                'error - .a and .b files not consistent:',   &
                '.a,.b min = ',xmin,xmin2 ,xmin2-xmin ,      &
                '.a,.b max = ',xmax,xmax2 ,xmax2-xmax
                  print *,'Variable :',char8
               endif
               call xcstop('(mod_CAL06_read_norwecomcarbon)')
               stop '(mod_CAL06_read_norwecomcarbon)'
            end if

            ! Put fields into its place:
            if ( trim(char8)=='pres' .and. k /= 0 ) then
               da%p(:,:,k)=tmpfld
            elseif ( trim(char8)=='DIA' .and. k /= 0 ) then
               da%dia(:,:,k)=tmpfld
            elseif ( trim(char8)=='FLA' .and. k /= 0 ) then
               da%fla(:,:,k)=tmpfld
            endif

           end if  !ios==0
         enddo

         call zaiocl(nop)
         if (mnproc==1) close(nop)

!AS         mod_carb=(da%fla(:,:,:)+da%dia(:,:,:))*red*12.0/(14.0)
!prøver kun med diatomer:
         mod_carb=(da%dia(:,:,:))*red*12.0/(14.0)
         where(mod_carb*14.0/(red*12.0) < 0.2) mod_carb = 0.2*red*12.0/14.0

!save carbon content and dp to file:
         if (mnproc==1) &
           open (unit=nop,file='NORWECOM_carbon/Carbon_'//rt%cyy//'_'//rt%cdd//'.b', &
                 status='replace')
         call zaiopf('NORWECOM_carbon/Carbon_'//rt%cyy//'_'//rt%cdd//'.a','replace',nop)

         ! Write .b Header
         write(nop,116) locctitle,lociversn,lociexpt,locyrflag, &
           l_itdm,l_jtdm, &
           l_kdm,da%rtinit%iyy,da%rtinit%idd, &
           da%rtdump%iyy,da%rtdump%idd,da%nrmem

         do k=1,kdm
           call zaiowr(da%p(1-nbdy,1-nbdy,k),ip,.false., &
                xmin,xmax, nop, .false.)
           char8='davp    '
           if     (mnproc.eq.1) then
             write (nop,117) char8,locnstep,loctime,k,coord,xmin,xmax
             call flush(nop)
           endif !1st tile  '
         end do

         do k=1,kdm
           call zaiowr(mod_carb(1-nbdy,1-nbdy,k),ip,.false., &
              xmin,xmax, nop, .false.)
           char8='davc    '
           if     (mnproc.eq.1) then
             write (nop,117) char8,locnstep,loctime,k,coord,xmin,xmax
             call flush(nop)
           endif !1st tile  '
         end do
!
         call zaiocl(nop)
         if (mnproc==1) close(nop)
       end if

!       print*, 'MOD_CARB1 ', mod_carb(50,50,1:2)
!       print*, 'MOD_CARB2 ', mod_carb(100,100,1:2)
!       print*, 'MOD_CARB3 ', mod_carb(150,150,1:2)

!do k=1,kdm
!print*, 'Carbon',mod_carb(30,30,k),mod_carb(100,100,k),mod_carb(80,80,k),mod_carb(50,50,k)
!enddo
!print*, 'Fla', da%fla(30,30,:,6)
!print*, 'Dia', da%dia(30,30,:,7)
!pause

 116  format (a80/a80/a80/a80/ &
      i5,4x,'''iversn'' = hycom version number x10'/ &
      i5,4x,'''iexpt '' = experiment number x10'/ &
      i5,4x,'''yrflag'' = days in year flag'/ &
      i5,4x,'''idm   '' = longitudinal array size'/ &
      i5,4x,'''jdm   '' = latitudinal  array size'/ &
      i5,4x,'''kdm   '' = Vertical     array size'/ &
      i5,4x,'''syear '' = Year of integration start '/ &
      i5,4x,'''sday  '' = Day of integration start'/ &
      i5,4x,'''dyear '' = Year of this dump      '/ &
      i5,4x,'''dday  '' = Day of this dump     '/ &
      i5,4x,'''count '' = Ensemble counter       ')
!     'field       time step  model day',&
!     '  k  dens        min              max')
 117  format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)
 118  format ('member ',i5.5,' = ',i1,' Ensemble member flag')
 119  format('field         time step  model day',&
      '  k  dens        min              max')

   end subroutine read_daily_carbon

   function da_filebase(rungen,year,day)
      character(len=3)    , intent(in) :: rungen
      character(len=4)    , intent(in) :: year 
      character(len=3)    , intent(in) :: day
      character(len=60) :: da_filebase

      da_filebase=rungen//'DAILY_'// &
        year//'_'//day

   end function da_filebase
     
end module mod_CAL06_read_norwecomcarbon