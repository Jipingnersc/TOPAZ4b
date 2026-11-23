      module m_NOR05_riverloads
      use mod_xc
      contains

!**********************************************************************

      subroutine NOR05_readnutrloads()
      implicit none

      call NOR05_nitrloads()
      call NOR05_phosloads()
      call NOR05_sililoads()
!AS      call NOR05_donoloads()

      end subroutine NOR05_readnutrloads

!**********************************************************************

      subroutine NOR05_rivloads(n)
      use mod_necessary_ecovars
      implicit none
      include 'common_blocks.h'

      integer,intent(in) :: n
      integer :: i,j,k

!AS      print*, 'maxnit', maxval(rivnit)
!AS      print*, 'maxpho', maxval(rivpho)
!AS      print*, 'maxsil', maxval(rivsil)
!AS      print*, 'maxdon', maxval(rivdon)



!AS divide by dp-layer depth when using TRIP-type nutrients
! distribute over the upper 5 layers
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,5
      if (sum(dp(i,j,1:5,n))>0) then
         bio(i,j,k,n,init)=bio(i,j,k,n,init)+
     &               rivnit(i,j)*biodt_sou*onem/sum(dp(i,j,1:5,n))

         bio(i,j,k,n,ipho)=bio(i,j,k,n,ipho)+
     &               rivpho(i,j)*biodt_sou*onem/sum(dp(i,j,1:5,n))

         bio(i,j,k,n,isil)=bio(i,j,k,n,isil)+
     &               rivsil(i,j)*biodt_sou*onem/sum(dp(i,j,1:5,n))

!AS  take out because it causes detp to be negative
!AS         bio(i,j,k,n,idet)=bio(i,j,k,n,idet)+
!AS     &               rivdon(i,j)*biodt_sou*onem/sum(dp(i,j,1:5,n))
      end if
      end do
      end do
      end do

      end subroutine NOR05_rivloads

!**********************************************************************

      subroutine NOR05_nitrloads()
       
      use mod_necessary_ecovars
      use mod_xc
      use mod_za
      use mod_year_info, only : year_info, rt
      use mod_checknan

      implicit none
      include 'common_blocks.h'
      
      logical lexist
      character(len=20)::filename 
      character(len=4) :: char8
      character(len=24) :: txt
      character preambl(5)*79
      real :: xmin,xmax
      real :: xmin2,xmax2
      integer :: month
      integer :: nop1,nop2,nop3
      integer :: nfields,ifld
      integer :: k
      
!--- Reading nitrate forcing field
      nop1=710
      nop2=2710

      filename='forcing.rivnitr'
      inquire(file=trim(filename)//'.a',exist=lexist) 
      if(.not.lexist) then
         write(lp,*)
         write(lp,*) '***** no nitrate riverload *****'
         write(lp,*)
      else
        if(mnproc==1) then
           write(lp,'(a,i4)')'read nitrate riverloads'
        endif
      endif
      call zaiopf(trim(filename)//'.a','old',nop1)
      open(nop2,file=trim(filename)//'.b',status='old',action='read')
      read(nop2,'(a79)') preambl
      call preambl_print(preambl)
c
c --- read desired month
c
      month=rt%imm
      if (mnproc.eq.1) print*, 'month', rt
      if (mnproc.eq.1) print*, 'month', month
      do ifld=1,month
        read(nop2,'(a24,i2.2,2e16.8)')  
     &                  txt,k,xmin2, xmax2
        call zaiord(rivnit(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax,nop1)
          if     (mnproc.eq.1) then
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &           '.a and .b files max  min:',
     &           '.a,.b min = ',xmin,xmin2 ,xmin2-xmin ,
     &           '.a,.b max = ',xmax,xmax2 ,xmax2-xmax
            print *,'Month :',ifld
          endif
        ! Check xmin/xmax
        if (abs(xmin-xmin2)>abs(xmin)*1e-4 .or.
     &      abs(xmax-xmax2)>abs(xmax)*1e-4) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &           'error - .a and .b files not consistent:',
     &           '.a,.b min = ',xmin,xmin2 ,xmin2-xmin ,
     &           '.a,.b max = ',xmax,xmax2 ,xmax2-xmax
            print *,'Month :',ifld
          endif
          call xcstop('(NOR05_nitrloads)')
          stop '(NOR05_nitrloads)'
        end if
      enddo
      call zaiocl(nop1)
      close(nop2)

      end subroutine NOR05_nitrloads

!**********************************************************************

      subroutine NOR05_phosloads()

      use mod_necessary_ecovars
      use mod_xc
      use mod_za
      use mod_year_info, only : year_info, rt
      use mod_checknan

      implicit none
      include 'common_blocks.h'

      logical lexist
      character(len=20)::filename 
      character(len=4) :: char8
      character(len=24) :: txt
      character preambl(5)*79
      real :: xmin,xmax
      real :: xmin2,xmax2
      integer :: month
      integer :: nop1,nop2,nop3
      integer :: nfields,ifld
      integer :: k
      integer :: i,j
      
!--- Reading phosphate forcing field
      nop1=710
      nop2=2710

      filename='forcing.rivphos'
      inquire(file=trim(filename)//'.a',exist=lexist) 
      if(.not.lexist) then
         write(lp,*)
         write(lp,*) '***** no phosphate riverload *****'
         write(lp,*)
      else
         if(mnproc==1) then
            write(lp,'(a,i4)')'read phosphate riverloads'
         endif
      endif
      call zaiopf(trim(filename)//'.a','old',nop1)
      open(nop2,file=trim(filename)//'.b',status='old',action='read')
      read(nop2,'(a79)') preambl
      call preambl_print(preambl)
c
c --- read desired month
c
      month=rt%imm
      if (mnproc.eq.1) print*, 'month', month
      do ifld=1,month
        read(nop2,'(a24,i2.2,2e16.8)')  
     &                  txt,k,xmin2, xmax2
        call zaiord(rivpho(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax,nop1)
        ! Check xmin/xmax
        if (abs(xmin-xmin2)>abs(xmin)*1e-4 .or.
     &      abs(xmax-xmax2)>abs(xmax)*1e-4) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &           'error - .a and .b files not consistent:',
     &           '.a,.b min = ',xmin,xmin2 ,xmin2-xmin ,
     &           '.a,.b max = ',xmax,xmax2 ,xmax2-xmax
            print *,'Month :',ifld
          endif
          call xcstop('(NOR05_phosloads)')
          stop '(NOR05_phosloads)'
        end if
      enddo
      call zaiocl(nop1)
      close(nop2)

      end subroutine NOR05_phosloads

!**********************************************************************
c      
      subroutine NOR05_sililoads()
     
      use mod_necessary_ecovars
      use mod_xc
      use mod_za
      use mod_year_info, only : year_info, rt
      use mod_checknan

      implicit none
      include 'common_blocks.h'
      
      logical lexist
      character(len=20)::filename 
      character(len=4) :: char8
      character(len=24) :: txt
      character preambl(5)*79
      real :: xmin,xmax
      real :: xmin2,xmax2
      integer :: month
      integer :: nop1,nop2,nop3
      integer :: nfields,ifld
      integer :: k
      integer :: i,j
      
!--- Reading silicate forcing field
      nop1=710
      nop2=2710

      filename='forcing.rivsili'
      inquire(file=trim(filename)//'.a',exist=lexist) 
      if(.not.lexist) then
         write(lp,*)
         write(lp,*) '***** no silicate riverload *****'
         write(lp,*)
      else
         if(mnproc==1) then
            write(lp,'(a,i4)')'read silicate riverloads'
         endif
      endif
      call zaiopf(trim(filename)//'.a','old',nop1)
      open(nop2,file=trim(filename)//'.b',status='old',action='read')
      read(nop2,'(a79)') preambl
      call preambl_print(preambl)
c
c --- read desired month
c
      month=rt%imm
      if (mnproc.eq.1) print*, 'month', month
      do ifld=1,month
        read(nop2,'(a24,i2.2,2e16.8)')  
     &                  txt,k,xmin2, xmax2
        call zaiord(rivsil(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax,nop1)
        ! Check xmin/xmax
        if (abs(xmin-xmin2)>abs(xmin)*1e-4 .or.
     &      abs(xmax-xmax2)>abs(xmax)*1e-4) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &           'error - .a and .b files not consistent:',
     &           '.a,.b min = ',xmin,xmin2 ,xmin2-xmin ,
     &           '.a,.b max = ',xmax,xmax2 ,xmax2-xmax
            print *,'Month :',ifld
          endif
          call xcstop('(NOR05_sililoads)')
          stop '(NOR05_sililoads)'
        end if
      enddo
      call zaiocl(nop1)
      close(nop2)
 
      end subroutine NOR05_sililoads

!**********************************************************************
c      
      subroutine NOR05_donoloads()
     
      use mod_necessary_ecovars
      use mod_xc
      use mod_za
      use mod_year_info, only : year_info, rt
      use mod_checknan

      implicit none
      include 'common_blocks.h'
      
      logical lexist
      character(len=20)::filename 
      character(len=4) :: char8
      character(len=24) :: txt
      character preambl(5)*79
      real :: xmin,xmax
      real :: xmin2,xmax2
      integer :: month
      integer :: nop1,nop2,nop3
      integer :: nfields,ifld
      integer :: k
      integer :: i,j
      
!--- Reading DON forcing field
      nop1=710
      nop2=2710

      filename='forcing.rivdono'
      inquire(file=trim(filename)//'.a',exist=lexist) 
      if(.not.lexist) then
         write(lp,*)
         write(lp,*) '***** no DON riverload *****'
         write(lp,*)
      else
         if(mnproc==1) then
            write(lp,'(a,i4)')'read DON riverloads'
         endif
      endif
      call zaiopf(trim(filename)//'.a','old',nop1)
      open(nop2,file=trim(filename)//'.b',status='old',action='read')
      read(nop2,'(a79)') preambl
      call preambl_print(preambl)
c
c --- read desired month
c
      month=rt%imm
      if (mnproc.eq.1) print*, 'month', month
      do ifld=1,month
        read(nop2,'(a24,i2.2,2e16.8)')  
     &                  txt,k,xmin2, xmax2
        call zaiord(rivdon(1-nbdy,1-nbdy),ip,.false.,
     &                  xmin,xmax,nop1)
        ! Check xmin/xmax
        if (abs(xmin-xmin2)>abs(xmin)*1e-4 .or.
     &      abs(xmax-xmax2)>abs(xmax)*1e-4) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &           'error - .a and .b files not consistent:',
     &           '.a,.b min = ',xmin,xmin2 ,xmin2-xmin ,
     &           '.a,.b max = ',xmax,xmax2 ,xmax2-xmax
            print *,'Month :',ifld
          endif
          call xcstop('(NOR05_donoloads)')
          stop '(NOR05_donoloads)'
        end if
      enddo
      call zaiocl(nop1)
      close(nop2)
 
      end subroutine NOR05_donoloads

!**********************************************************************

      end module m_NOR05_riverloads
