      module m_NOR05_fx

      contains

      real function NOR05_fx(d,e,sunhei)
C
C     COMPUTE THE SOLAR ALTITUDE FUNCTION
C
C     VARAIBLES GOING INTO THE FUNCTION
C     D,E   : COEFFISIENTS, SEE OLSETH ANS SKAARTVEIT
C     SUNHEI: SOLAR ALTITUDE
      implicit none
C
      real :: d,e,sunhei
C
      NOR05_fx=1-d-e+d*sin(sunhei)+e*sqrt(sin(sunhei))
C
      return
      end function NOR05_fx
      end module m_NOR05_fx
