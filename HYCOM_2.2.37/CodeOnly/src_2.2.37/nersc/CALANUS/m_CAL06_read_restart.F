      module m_CAL06_read_restart
      use mod_year_info

      contains

      subroutine read_restart_calanus(rungen,rt,nflt)
      use mod_xc
!      use mod_hycomfuncs
!      use mod_syntheticfloats
      use mod_CAL06_calnusparticles
      implicit none
      
!AS      logical ::  read_restart_calanus
      type(year_info) , intent(in) :: rt
      character(len=3), intent(in) :: rungen
      integer, intent(in) :: nflt

      character(len=70) :: calname
      logical ::ex
      integer :: i,j,iostat
      type(calfin) :: tmpcal
      integer :: calanusdd, calanusyy
      real :: otherdata(15)
      character*12 :: tmpstr1 
      character*8  :: tmpstr2
      character*45 :: tmpstr3 
      integer :: itmp, nfl
      integer :: mortnn, mortdd,mortyy
      
      real :: rtmp
       
      do i=1,mpop
         do j=1,3
            cfs(i)%pos(j)=-999.0
         end do
      end do 

!AS      read_restart_calanus=.true.
      calname=
     & ''//rungen//'CAL_y'//rt%cyy//'_d'//rt%cdd//'_h'//rt%chh//'.uf'
      inquire(file=trim(calname),exist=ex)
      if (ex) then
         open(10,file=trim(calname),form='unformatted',action='read')
         read(10) nflt
         !print*, 'nflt', nflt
         do nfl=1,nflt
            read(10) cfs(nfl), itmp, itmp
!            print*, cfs(nfl)%sv(1)
!            stop
         end do
         close(10)
!         print*, cfs(1)
      else 
         if (mnproc==1) then
            write(lp,*) 'Restart file for CALANUS does not exist.'
            write(lp,*) 'I continue with fields from input file'
         end if
!AS         read_restart_calanus=.false.
      endif

CAS if wake-up day is one day restart day, the upated vertical migation flag will not be
CAS read because output is only once a day  
      do i=1,nflt
         if (int(cfs(i)%sv(1)).eq.rt%idd-1) cfs(i)%av(11)=3
      end do

CAS      modify the mortality output files so the data is written in the same place
      print*, 'Starting mortality read/write'
 101  format(a12,i8,a8,i8,a45)
      iostat=0
      open(810,file='tmpmortality_out', 
     &         form='formatted', status='replace')
      open(811,file='mortality_out', form='formatted', status='old',
     &         iostat=iostat)
      do while (mortdd.le.rt%idd.and.iostat.eq.0)   
         read(811,101,iostat=iostat) tmpstr1, mortnn, tmpstr2, 
     &                                        mortdd, tmpstr3
         if ((mortdd.le.rt%idd.or.calanusyy.lt.rt%iyy)
     &        .and.iostat.eq.0) then
            write(810,101) tmpstr1, mortnn, tmpstr2, mortdd, tmpstr3
         end if
      end do
      close(810)
      close(811)

      iostat=0
      open(810,file='tmpmortality_out', form='formatted', status='old',
     &               iostat=iostat)
      open(811,file='mortality_out', form='formatted', status='replace')
      do while (iostat.eq.0) 
         read(810,101,iostat=iostat)tmpstr1, mortnn, tmpstr2, 
     &                                       mortdd, tmpstr3 
         if (iostat.eq.0) then 
            write(811,101) tmpstr1, mortnn, tmpstr2, mortdd, tmpstr3
         end if
      end do
      close(810)
      close(811)
      print*, 'Ending mortality read/write'

CAS      modify the biological mortality output files so the data is written in the same place
      print*, 'Starting biological mortality read/write'
      iostat=0
      open(810,file='tmpfile',form='unformatted', status='replace')
      open(811,file='Biol_mortality.dat',form='unformatted',
     &    status='old',iostat=iostat)
      mortdd=0
      do while (mortdd.le.rt%idd.and.iostat.eq.0)   
         read(811,iostat=iostat) mortyy,mortdd
         read(811,iostat=iostat) zm_cal
         read(811,iostat=iostat) zn_cal
         read(811,iostat=iostat) zb_cal
         if ((mortdd.le.rt%idd.or.calanusyy.lt.rt%iyy)
     &      .and.iostat.eq.0) then
           write(810) mortyy,mortdd,zm_cal,zn_cal,zb_cal
         endif
      end do
      close(810)
      close(811)

      iostat=0
      open(810,file='tmpfile', form='unformatted', status='old',
     &               iostat=iostat)
      open(811,file='Biol_mortality.dat', form='unformatted', 
     &               status='replace')
      do while (iostat.eq.0) 
        read(810,iostat=iostat) mortyy,mortdd,zm_cal,zn_cal,zb_cal
          if (iostat.eq.0) then 
            write(811) mortyy,mortdd
            write(811) zm_cal
            write(811) zn_cal
            write(811) zb_cal
          end if
        end do
      close(810)
      close(811)
      print*, 'Ending biological mortality read/write'

      end subroutine read_restart_calanus

      end module m_CAL06_read_restart


