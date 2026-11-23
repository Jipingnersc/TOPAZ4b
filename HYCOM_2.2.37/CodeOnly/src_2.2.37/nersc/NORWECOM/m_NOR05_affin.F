      module m_NOR05_affin

      contains

      subroutine NOR05_affin
!
!***BEGIN PROLOGUE AFFIN
!***AUTHOR
!     K{\aa}re B. Ulvestad,
!     Institute of Fishery and Marine Biology
!     University of Bergen, Norway
!
!***REVISED DATE 25/8-92
!     Morten D. Skogen, Institute of Marine Research
!     Postbox 1870, N-5024 Bergen-Nordnes
!     Email: morten@imr.no
!
!     Added include and prologue
!
!***PURPOSE AFFIN is a subroutine for the computation of the affinity
!     at TMEAN degrees Celsius, when the half saturation is given
!     in mikroMolar
!
!***ROUTINES CALLED : NONE
!
!***END PROLOGUE AFFIN
!
!   Global variables
!
!   On entry :
!
!
!   On exit :
!
!     AA(1)    : Diatom production maximum at 0 degree celsius (s^{-1})
!     AA(2)    : Temperature dependent PMAX for the diatoms (C^{-1})
!     AA(3)    : Flagellate production maximum at 0 degree celsius (s^{-1})
!     AA(4)    : Temperature dependent PMAX for the flagellates (C^{-1})
!     AA(5)    : Respiration rate at 0 degrees celcius (s^{-1})
!     AA(6)    : Respiration rate temperature dependence (C^{-1})
!     KN_FLA   : Half saturation nitrogen flagellates (uMN)
!     KP_FLA   : Half saturation phosphate flagellates (uMP)
!     KR_FLA   : Half saturation radiation flagellates (uE/m^2/s)
!     KN_DIA   : Half saturation nitrogen diatoms (uMN)
!     KP_DIA   : Half saturation phosphate diatoms (uMP)
!     KS_DIA   : Half saturation silicate diatoms (uMSi)
!     KR_DIA   : Half saturation radiation diatoms (uE/m^2/s)
!     KN_CHA   : Half saturation nitrogen chatonella (uMN)
!     KP_CHA   : Half saturation phosphate chatonella (uMP)
!     KR_CHA   : Half saturation radiation chatonella (uE/m^2/s)
      use mod_xc
!
      implicit none
      include 'biocom.h'
!
!   Local variables :
!
!     pmax : production maximum at tmean degrees celsius
!     TMEAN: mean temperature
!
      real  pmax,tmean
!  
!
      AA(1)=1.15e-5!1.53E-5
      AA(2)=0.063
      AA(3)=0.76e-5!1.02E-5
      AA(4)=0.063
      AA(5)=8.05E-7
      AA(6)=0.070 
! 
      KN_DIA=2.0
      KP_DIA=0.125
      KS_DIA=1.4
      KR_DIA=96.
      KN_FLA=1.5
      KP_FLA=0.094
      KR_FLA=209.
!
      TMEAN=13.0
!
! Affinity for diatoms
!
      pmax=aa(1)*exp(tmean*aa(2))
 
      kn_dia=pmax/(cnit*kn_dia)
      kp_dia=pmax/(cpho*kp_dia)
      ks_dia=pmax/(csil*ks_dia)
      kr_dia=pmax/kr_dia
!
! Affinity for flagellates
!
      pmax=aa(3)*exp(tmean*aa(4))
 
      kn_fla=pmax/(cnit*kn_fla)
      kp_fla=pmax/(cpho*kp_fla)
      kr_fla=pmax/kr_fla
!
!#ifdef CHATONELLA
!
! Half saturation for Chatonella
! Nutrients from flagellates, light from Naustvoll unpubl.data
!
!      kn_cha = 1.500*cnit
!      kp_cha = 0.094*cpho
!      kr_cha =30.
!
!      kp_cha = 0.04*cpho
!      kn_cha = 1.000*cnit
!#endif
!
      return

      end subroutine NOR05_affin
      end module m_NOR05_affin
