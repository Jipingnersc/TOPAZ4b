      module m_NOR05_relaxation
      use mod_xc
      use mod_necessary_ecovars, only: init, ipho, isil, idet, isis,
     &                                 ifla, idia, ioxy, ised, iyel,
     &                                 icha, idetp,imeso,imicro
      contains

!**********************************************************************

      subroutine NOR05_readnutrelax()
      use mod_year_info, only : year_info, rt, refyear, year_day

      implicit none
      include 'common_blocks.h'
c
      integer month
      call year_day(time,refyear,rt,yrflag)

      month=rt%imm

      
      call NOR05_nitrrelax(mod(month+10,12)+1,1)
      call NOR05_nitrrelax(month             ,2)
      call NOR05_nitrrelax(mod(month,   12)+1,3)
      call NOR05_nitrrelax(mod(month+1, 12)+1,4)

      call NOR05_phosrelax(mod(month+10,12)+1,1)
      call NOR05_phosrelax(month             ,2)
      call NOR05_phosrelax(mod(month,   12)+1,3)
      call NOR05_phosrelax(mod(month+1, 12)+1,4)

      call NOR05_silirelax(mod(month+10,12)+1,1)
      call NOR05_silirelax(month             ,2)
      call NOR05_silirelax(mod(month,   12)+1,3)
      call NOR05_silirelax(mod(month+1, 12)+1,4)

      call NOR05_oxygrelax(mod(month+10,12)+1,1)
      call NOR05_oxygrelax(month             ,2)
      call NOR05_oxygrelax(mod(month,   12)+1,3)
      call NOR05_oxygrelax(mod(month+1, 12)+1,4)

      rmutr(:,:,idet)=0.0
      rmutr(:,:,isis)=0.0
      rmutr(:,:,ifla)=0.0
      rmutr(:,:,idia)=0.0
#if defined (YELLOWSED)
      rmutr(:,:,ised)=0.0
      rmutr(:,:,iyel)=0.0
#endif /*YELLOWSED*/
      rmutr(:,:,idetp)=0.0
      rmutr(:,:,imeso)=0.0
      rmutr(:,:,imicro)=0.0

      if (mnproc==3) then
      
      print*, 'month', month, 'oxy trwall', trwall(100,1:5,1,1,ioxy), 
     &    'mnproc', mnproc      
      print*, 'month', month, 'oxy trwall', trwall(100,1:5,1,2,ioxy), 
     &    'mnproc', mnproc     
      print*, 'month', month, 'oxy trwall', trwall(100,1:5,1,3,ioxy), 
     &    'mnproc', mnproc      
      print*, 'month', month, 'oxy trwall', trwall(100,1:5,1,4,ioxy), 
     &    'mnproc', mnproc
      
      endif
      end subroutine NOR05_readnutrelax

!**********************************************************************

      subroutine NOR05_nitrrelax(month,lslot)

      use mod_za
c
      implicit none
      include 'common_blocks.h'
c
      integer lslot,month
c
      integer   i,k,ios,lgth 
      integer   mnth
      character preambl(5)*79, cline*80
      real      xmin, xmax

      lgth = len_trim(flnmfor)
      mnth=mod(month-1,12)+1

      call zaiopf(flnmforw(1:lgth)//'relax.nitr.a', 'old', 931)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
      open (unit=uoff+931,file=flnmforw(1:lgth)//'relax.nitr.b',
     &        status='old', action='read')
      read (uoff+931,'(a79)') preambl
      endif !1st tile
      call preambl_print(preambl)

      do i=1,month
        do k= 1,kk
          call zaiord(trwall(1-nbdy,1-nbdy,k,lslot,init),ip,.false.,
     &          xmin,xmax,931)
        enddo
      enddo

      trwall(:,:,:,lslot,init)=trwall(:,:,:,lslot,init)*14.0
      
      call zaiocl(931)
      close(uoff+931)

       
      end subroutine NOR05_nitrrelax

!**********************************************************************

      subroutine NOR05_phosrelax(month,lslot)

      use mod_za
c
      implicit none
      include 'common_blocks.h'
c
      integer lslot,month
c
      integer   i,k,ios,lgth
      integer   mnth
      character preambl(5)*79, cline*80
      real      xmin, xmax

      lgth = len_trim(flnmfor)
      mnth=mod(month-1,12)+1

       call zaiopf(flnmforw(1:lgth)//'relax.phos.a', 'old', 932)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
      open (unit=uoff+932,file=flnmforw(1:lgth)//'relax.phos.b',
     &        status='old', action='read')
      read (uoff+932,'(a79)') preambl
      endif !1st tile
      call preambl_print(preambl)

       do i=1,month
        do k= 1,kk
          call zaiord(trwall(1-nbdy,1-nbdy,k,lslot,ipho),ip,.false.,
     &          xmin,xmax,932)
        enddo
      enddo     

      trwall(:,:,:,lslot,ipho)=trwall(:,:,:,lslot,ipho)*31.0

      call zaiocl(932)
      close(uoff+932)

      end subroutine NOR05_phosrelax

!**********************************************************************
c      
      subroutine NOR05_silirelax(month,lslot)

      use mod_za
c
      implicit none
      include 'common_blocks.h'
c
      integer lslot,month
c
      integer   i,k,ios,lgth
      integer   mnth
      character preambl(5)*79, cline*80
      real      xmin, xmax

      lgth = len_trim(flnmfor)
      mnth=mod(month-1,12)+1

       call zaiopf(flnmforw(1:lgth)//'relax.sili.a', 'old', 933)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
      open (unit=uoff+933,file=flnmforw(1:lgth)//'relax.sili.b',
     &        status='old', action='read')
      read (uoff+933,'(a79)') preambl
      endif !1st tile
      call preambl_print(preambl)

       do i=1,month
        do k= 1,kk
          call zaiord(trwall(1-nbdy,1-nbdy,k,lslot,isil),ip,.false.,
     &          xmin,xmax,933)
        enddo
      enddo     
 
      trwall(:,:,:,lslot,isil)=trwall(:,:,:,lslot,isil)*28.0

      call zaiocl(933)
      close(uoff+933)
     
      end subroutine NOR05_silirelax

!**********************************************************************
c      
      subroutine NOR05_oxygrelax(month,lslot)

      use mod_za
c
      implicit none
      include 'common_blocks.h'
c
      integer lslot,month
c
      integer   i,k,ios,lgth
      integer   mnth
      character preambl(5)*79, cline*80
      real      xmin, xmax

      lgth = len_trim(flnmfor)
      mnth=mod(month-1,12)+1

       call zaiopf(flnmforw(1:lgth)//'relax.oxyg.a', 'old', 934)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
      open (unit=uoff+934,file=flnmforw(1:lgth)//'relax.oxyg.b',
     &        status='old', action='read')
      read (uoff+934,'(a79)') preambl
      endif !1st tile
      call preambl_print(preambl)

      do i=1,month
        do k= 1,kk
          call zaiord(trwall(1-nbdy,1-nbdy,k,lslot,ioxy),ip,.false.,
     &          xmin,xmax,934)
        enddo
      enddo     
     
      trwall(:,:,:,lslot,ioxy)=trwall(:,:,:,lslot,ioxy)*1330.0


      call zaiocl(934)
      close(uoff+934)
     
      end subroutine NOR05_oxygrelax

!**********************************************************************

      end module m_NOR05_relaxation
