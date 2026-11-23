       module m_NOR05_snell

       contains

       real function NOR05_snell(SUNHEI)
!***
!***BEGIN PROLOGUE SNELL
!***DATE WRITTEN 22/9-92
!***AUTHOR
!      Morten D. Skogen, Institute of Marine Research
!      Postbox 1870, N-5024 Bergen-Nordnes
!      Email: morten@imr.no
!
!***PURPOSE SNELL computes zenith angel of direct light in the water
!      column using Snells formula
!
!***ROUTINES CALLED : NONE
!
!***END PROLOGUE
!
!   Global variables:
!
!   On entry :
!
!     SUNHEI : The sunheight
!
!   On exit :
!
!     SNELL : The computed angel due to Snell's law
       use mod_xc
!
      implicit none 
!
      REAL SUNHEI
      include 'biocom.h'
!
!***First executable statement SLIM
!
      NOR05_SNELL = ASIN(COS(SUNHEI)*RC1/RC2)
!
      RETURN
      end function NOR05_snell
      end module m_NOR05_snell
