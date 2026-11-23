      module m_NOR05_slim

      contains
  
      real function NOR05_slim(alpha,pmax,s)
C***
C***BEGIN PROLOGUE SLIM
C***AUTHOR
C      K{\aa}re B. Ulvestad, Department of Fisheries and
C      Marine Biology, Bergen High Technology Centre,
C      5014 Bergen, Norway
C
C***Revised date 26/8-92 Morten D. Skogen, morten@imr.no
C     Added prologue
C
C***PURPOSE SLIM calculates the Michaels-Menten limitation [0..1)
C
C***ROUTINES CALLED : NONE
C
C***END PROLOGUE
C
C   Global variables:
C
C   On entry :
C
C      ALPHA : Affinity
C      PMAX  : Production maximum
C      S     : Substrate concentration
C
C   On exit :
C
C     SLIM : The  Michaels-Menten limitation
C
      implicit none 
      REAL alpha,pmax,s
C
C Local variables
C
C      K : Half-saturation
C
      REAL k
C
C First executable statement SLIM
C
      k=pmax/alpha
 
      if(s.lt.0)then
        NOR05_slim=0.
      else
        NOR05_slim=s/(k+s)
      endif
C
      return
      end function NOR05_slim
      end module m_NOR05_slim
