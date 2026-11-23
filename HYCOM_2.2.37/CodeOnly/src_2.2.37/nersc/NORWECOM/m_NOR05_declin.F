      module m_NOR05_declin

      contains
 
      subroutine NOR05_declin(daynu,theta_nor)
C***
C***BEGIN PROLOGUE DECLIN
C***DATE WRITTEN 3/2-88
C***AUTHOR
C      Sam Wakili
C
C***Revised date 1/9-92 Morten D. Skogen, morten@imr.no
C     Added prologue
C
C***PURPOSE DECLIN computes the declination of the sun as a
C     function of the daynr. The declination is the elevation of
C     sun measured over the North pole given in degrees.
C
C***REFERENCE : Meteorological report series, University of Bergen,
C     1-1987. OLSETH, J.A. & A. SKARTVEIT, equation 3, p.4.
C
C***ROUTINES CALLED : NONE
C
C***END PROLOGUE
C
C   Global variables:
C
C   On entry :
C
C     DAYNR : day number (1-365)
C
C   On exit :
C
C     THETA : the declination
C
      implicit none

      real theta_nor
      integer,intent(in)::  daynu
      real              :: daynr
C
C   Local variables
C
      real pi,a,b
c --- CH:081105: adding one day to daynr, as hycom starts from 0 and pom
c                from 1
      daynr=real(daynu)+1.
      
C
C***FIRST executable statement DECLIN
C
      pi=3.141593
C
      a = (0.9856*daynr) * pi / 180.0
      a = sin(a)
      b = (0.9856 * (daynr-80.0)) + (1.9171 * (a - 0.98112))
      b = b * pi / 180.0
C
      theta_nor = asin (0.3979 * sin(b))
      theta_nor = theta_nor * 180.0 / pi
C
      return
      end subroutine NOR05_declin
      end module m_NOR05_declin
