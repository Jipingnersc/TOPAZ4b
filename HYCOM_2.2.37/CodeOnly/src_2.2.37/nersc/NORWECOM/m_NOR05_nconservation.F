      module m_NOR05_nconservation
      
      contains
      
      subroutine nconservation(string,n)
      

      use mod_xc
      use mod_necessary_ecovars

      implicit none

      character(len=*), intent(in) :: string
      integer, intent(in) :: n
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   tmp
#ifdef NOR05      
      real*8 :: totsum_nit,totsum_pho,totsum_sil 
#endif      
      real*8 :: totsum
      real :: pho_nit, sil_nit, sis_nit, oxy_nit
      integer i,j,l,k,o

      include 'common_blocks.h'
      include 'sedcom.h'
      include 'biocom.h'
      
#ifdef NOR05      
      totsum_nit=0.
      totsum_pho=0.
      totsum_sil=0.
#endif      
      tmp=0.
      totsum=0.

! compute the total amount of nitrogen in the model
      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point 
      !margin=0
#ifdef NOR05
!--- CH:need to convert the consentrations of phosphate, oxygen, silica
!--- and biogenic silica to be in N and not P,Si, and O. 
!KAL      do j=1-margin,jj+margin
!KAL      do l=1,isp(j)
!KAL      do i=ifp(j,l),ilp(j,l)
!KAL      do k=1,kdm
!KAL!         totsum_nit=totsum_nit+bio(i,j,k,n,1)+bio(i,j,k,n,4)
!KAL!     &             +bio(i,j,k,n,6)+bio(i,j,k,n,7)
!KAL!         totsum_sil=totsum_sil+(1/cc(2))*bio(i,j,k,n,3)+bio(i,j,k,n,7)
!KAL!     &             +(1/cc(2))*bio(i,j,k,n,5)    
!KAL!         totsum_pho=totsum_pho+(1/cc(1))*bio(i,j,k,n,2)+bio(i,j,k,n,6)
!KAL!     &              +bio(i,j,k,n,7)+bio(i,j,k,n,4)    
!KAL      totsum_nit=totsum_nit+bio(i,j,k,n,1)
!KAL      totsum_pho=totsum_pho+bio(i,j,k,n,2)
!KAL      totsum_sil=totsum_sil+bio(i,j,k,n,3)
!KAL      enddo
!KAL      enddo
!KAL      enddo
!KAL      enddo
!KAL -- built in routines for global sum
      call xcsum(totsum_nit,bio(i,j,k,n,init),ip)
      call xcsum(totsum_pho,bio(i,j,k,n,ipho),ip)
      call xcsum(totsum_sil,bio(i,j,k,n,isil),ip)
      if (mnproc==1) 
     &   print*,'Total amount of N in the model; ',
     &      string, ': ', totsum_nit,totsum_pho,totsum_sil
#else
      do o=1,nbio
      do k=1,kdm
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=ifp(j,l),ilp(j,l)
         tmp(i,j)=bio(i,j,k,n,o)*scp2(i,j)*dp(i,j,k,n)/onem
      enddo
      enddo
      enddo
      enddo
      enddo
!      
      call xcsum(totsum,tmp,ip)
      if (mnproc==1) 
     &   print*,'Total amount of N in the model; ',
     &      string, ': ', totsum

!      write(*,'(a33,a20,a2,f25.4)')'Total amount of N in the model; ',
!     &      string, ': ', totsum
#endif      


      end subroutine nconservation

      subroutine nconservation_2_file(thistime,n)

      use mod_xc
      use mod_necessary_ecovars

      implicit none

      integer, intent(in) :: thistime
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   tmp1,tmp2,tmp3
#ifdef NOR05      
      real*8 :: totsum_nit,totsum_pho,totsum_sil
#endif      
      real*8 :: totsum
      real :: pho_nit, sil_nit, sis_nit, oxy_nit
      integer, intent(in) :: n
      integer i,j,l,k,kbio
      
      include 'common_blocks.h'
      include 'biocom.h'
      include 'sedcom.h'

#ifdef NOR05
      totsum_nit=0.
      totsum_sil=0.
      totsum_pho=0.
#endif      
      totsum=0.
      tmp1=0.
      tmp2=0.
      tmp3=0.

! compute the total amount of nitrogen in the model
      if (mnproc==1) print*,
     &   'RB: NOR05_nconservation: mulieg noko muffens her'
      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point 
      !margin=0
#ifdef NOR05
!--- CH:need to convert the consentrations of phosphate, oxygen, silica
!--- and biogenic silica to be in N and not P,Si, and O. 
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=ifp(j,l),ilp(j,l)
      do k=1,kdm
         tmp1(i,j)=bio(i,j,k,n,init)+bio(i,j,k,n,idet)
     &             +bio(i,j,k,n,ifla)+bio(i,j,k,n,idia)
         tmp2(i,j)=(1/cc(2))*bio(i,j,k,n,isil)+bio(i,j,k,n,idia)
     &             +(1/cc(2))*bio(i,j,k,n,isis)    
         tmp3(i,j)=(1/cc(1))*bio(i,j,k,n,ipho)+bio(i,j,k,n,ifla)
     &              +bio(i,j,k,n,idia)+bio(i,j,k,n,idet)    
      enddo
      enddo
      enddo
      enddo
      !KAL - built-in routines
      call xcsum(totsum_nit,tmp1,ip)
      call xcsum(totsum_sil,tmp2,ip)
      call xcsum(totsum_pho,tmp3,ip)

      if (mnproc==1) then
         open(10,file='total_nitrogen.dat',position='append',
     &        status='unknown')
         write(10,'(i5.5,x,3g14.4)')
     &      thistime, totsum_nit,totsum_pho,totsum_sil
         close(10)
      end if
#else
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=ifp(j,l),ilp(j,l)
      do kbio=1,nbio
      do k=1,kdm

        tmp1(i,j)=bio(i,j,k,n,kbio)*scp2(i,j)*dp(i,j,k,n)/onem

      enddo
      enddo
      enddo
      enddo
      enddo

      !KAL - built-in routines
      call xcsum(totsum,tmp1,ip)

      ! Only one process writes
      if (mnproc==1) then
         open(10,file='total_nitrogen.dat',position='append',
     &        status='unknown')
         write(10,*) thistime, totsum
         close(10)
      end if
#endif      

      

      end subroutine nconservation_2_file
      end module m_NOR05_nconservation
