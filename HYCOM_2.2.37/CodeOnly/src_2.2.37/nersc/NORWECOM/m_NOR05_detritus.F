      module m_NOR05_detritus

      contains

      subroutine NOR05_detritus(k_list,n)
 
C***BEGIN PROLOGUE DETRITUS
C***DATE WRITTEN   920402 (YYMMDD)
C***AUTHOR
C     K{\aa}re B. Ulvestad, Department of Fisheries and
C     Marine Biology, University of Bergen,
C     Bergen High Technology Center,
C     N-5020 BERGEN, Norway
C
C***REVISION DATE 25/8-92
C     Morten D. Skogen, Institute of Marine Research
C     Postbox 1870, N-5024 Bergen-Nordnes
C     Email: morten@imr.no
C
C     Added include and prologue
C
C***REVISION DATE 23/3-93
C     Morten D. Skogen, Institute of Marine Research
C     Changed death to equals CC(3) if the concentration somewhere
C     (not only at the surface) is greater than the minimum value
C     This enables the algae to extinct at the surface, and therefore
C     possible better chances for subsurface production.
C
C***PURPOSE DETRITUS  puts in a death rate where the concentration
C     in the surface cell exceeds a minimum value.  The
C     biological model is light limited in winter. In order to
C     prevent the algae in the model to extinct, we have to
C     prevent further death of the algae when their concentration
C     becomes to small.
C
C***ROUTINES CALLED : NONE
C
C***END PROLOGUE DETRITUS
C
C   Global variables
C
C   On entry :
C
C     See description in globals.doc
C
C   On exit :
C
C     DEATH_DIA : field for the death-rate of diatoms. Equals 0 if
C       the concentration of diatoms in the surface cell < DIAMIN
C     DEATH_FLA : same purpose for flagellates
C     DEATH_CHA : same purpose for chatonella
C      
      use mod_xc
      use mod_necessary_ecovars
C
      implicit none
C
      include 'common_blocks.h'
      include 'biocom.h'
C
      integer,intent(in) :: n
      integer,intent(in),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)::
     &          k_list
C
C Local variables
C
      integer i,j,k
C
      death_dia(:,:,:)=0.
      death_fla(:,:,:)=0.
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      if(depths(i,j).gt.0.5)then
C
          do k=1,k_list(i,j)
          if(bio(i,j,k,n,idia).gt.diamin)then
#if defined ZOOPL
      death_dia(i,j,k)=cc(3)
#else
      death_dia(i,j,k)=cc(3)*bio(i,j,k,n,idia)/15.0
#endif
          endif
C
          if(bio(i,j,k,n,ifla).gt.flamin)then
#if defined ZOOPL
      death_fla(i,j,k)=cc(3)
#else
      death_fla(i,j,k)=cc(3)*bio(i,j,k,n,ifla)/15.0
#endif
          endif
          enddo
      endif
      enddo
      enddo
C
!#ifdef CHATONELLA
!          death_cha(i,j,:)=0.
!          do 130 k=1,k_list(i,j)
!          if(cha(i,j,k).gt.chamin)then
!               death_cha(i,j,k)=0.1*cc(3)
!               goto 135
!          endif
! 130      continue
! 135      continue
!#endif
C
!       endif
 100  continue
 
      return
      end subroutine NOR05_detritus
      end module m_NOR05_detritus
