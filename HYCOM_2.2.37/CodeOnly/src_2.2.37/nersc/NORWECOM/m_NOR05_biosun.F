      module m_NOR05_biosun
      contains
      subroutine NOR05_biosun(daynu,ttime_nor,theta_nor,sunti,sunhei)
C
C***
C***BEGIN PROLOGUE BIOSUN
C***DATE WRITTEN  23/3-93
C***AUTHOR
C      Morten D. Skogen, Institute of Marine Research
C      Postbox 1870, N-5024 Bergen-Nordnes
C      Email: morten@imr.no
C
C***PURPOSE BIOSUN computes the suntime and sunheight in all
C     horizontal points.
C
C***DESCRIPTION :
C     BIOSUN does the work of SUNTI and SUNELEV. The routine computes
C     the sunelevation and suntime in one call instead of calling
C     the routines IM*JM times, and are possible more effective.
C
C***END PROLOGUE
      use mod_xc
      implicit none

      include 'common_blocks.h'
C
      integer,intent(in) ::  daynu
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: sunhei,sunti
      real,intent(in) :: theta_nor,ttime_nor
C
      integer i,j
      real nnorti,a,b,c,d
      real :: daynr
c --- CH:081105:adding 1 to daynr
      daynr=real(daynu)+1.
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C First find the suntime, 2D-version of subroutine SUNTIME
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      nnorti= ttime_nor * 360.0 / 24.0
 
      b=(0.9856*(daynr-80.0)) * pi / 180.0
      b=0.4083 * sin(b)
 
      c=(0.9856*(daynr-80.0)) * pi / 180.0
      c=1.7958 * cos(c)
 
      d=(1.9712*(daynr-80.0)) * pi / 180.0
      d=2.4875 * sin(d)
 
      do  j=1-margin,jj+margin
      do  i=1-margin,ii+margin
         a = plon(i,j) -15.
         sunti(i,j) = nnorti + a - b - c + d 
      enddo   
      enddo   
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Then find the sunelevation, 2D-version of subroutine SUNELEV
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      a= theta_nor * pi / 180.0
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         b = plat(i,j) * pi/180.
         c = sunti(i,j) * pi/180.
         sunhei(i,j) = (sin(a) * sin(b)) - (cos(a) * cos(b) * cos(c))
         sunhei(i,j) = asin(sunhei(i,j))
      enddo
      enddo
C
C***END subroutine BIOSUN
C
      return
      end subroutine NOR05_biosun
      end module m_NOR05_biosun
