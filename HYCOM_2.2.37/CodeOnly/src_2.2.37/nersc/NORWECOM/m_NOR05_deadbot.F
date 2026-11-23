      module m_NOR05_deadbot

      contains

      subroutine NOR05_deadbot(n,k_list)
!
!***BEGIN PROLOGUE DEADBOT
!***DATE WRITTEN   29/1-96
!***AUTHOR
!     HENRIK SOILAND. henrik@imr.no
!       Institute of Marine Research,
!       Postboks 1870, N-5024 Bergen-Nordnes, Norway
!
!***PURPOSE DEADBOT.
!
!       
!       Tranfer all diatoms and flagellates and chatonella at the bottom to 
!       detritus and diatoms also form silicate  shell.
!
!***END PROLOGUE DEADBOT
      use mod_xc
      use mod_necessary_ecovars
!
      implicit none
      include 'sedcom.h'
      include 'biocom.h'
      include 'common_blocks.h'
      
      integer,intent(in),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & k_list
      integer,intent(in):: n
      integer i,j,kmax

      margin=0
      

      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         kmax=k_list(i,j)
         if(dp(i,j,kmax,n).lt.10.) then
            bio(i,j,kmax,n,idet)=bio(i,j,kmax,n,idet) + 
     &               bio(i,j,kmax,n,idia) + bio(i,j,kmax,n,ifla)
#if defined (DETPHO)
            bio(i,j,kmax,n,idetp)=bio(i,j,kmax,n,idetp) + 
     &               cc(1)*(bio(i,j,kmax,n,idia) + bio(i,j,kmax,n,ifla))
#endif
            bio(i,j,kmax,n,isis)=bio(i,j,kmax,n,isis) + 
     &               cc(2)*bio(i,j,kmax,n,7)
            bio(i,j,kmax,n,idia)=0.0
            bio(i,j,kmax,n,ifla)=0.0
#ifdef CHATONELLA
         bio(i,j,kmax,n,idet)=bio(i,j,kmax,n,idet) + 
     &               bio(i,j,kmax,n,icha)
#if defined (DETPHO)
         bio(i,j,kmax,n,idetp)=bio(i,j,kmax,n,idetp) + 
     &               cc(1)*bio(i,j,kmax,n,icha)
#endif
         bio(i,j,kmax,n,icha)=0.0
#endif
         endif
      enddo
      enddo
!
      return
      end subroutine NOR05_deadbot
      end module m_NOR05_deadbot
