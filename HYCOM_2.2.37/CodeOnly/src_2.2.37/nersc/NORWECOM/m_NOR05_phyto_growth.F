      module m_NOR05_phyto_growth

      contains
 
      subroutine NOR05_phyto_growth(i,j,n)
C      SUBROUTINE PHYTO_GROWTH(I,J,T,S)
C***
C***BEGIN PROLOGUE PHYTO_GROWTH
C***DATE WRITTEN   920406 (YYMMDD)
C***AUTHOR
C      K{\aa}re B. Ulvestad, Department of Fisheries and
C      Marine Biology, Bergen High Technology Centre,
C      5014 Bergen, Norway
C
C***Revised date 26/8-92 Morten D. Skogen, morten@imr.no
C      Added prologue
C
C***PURPOSE PHYTO_GROWTH computes the growth and respiration
C      of diatoms and flagellates and chatonella.
C
C***DESCRIPTION
C      The computation is temperature dependent, and a function
C      of nutrient and light regime. The function is called once
C      for each pair of indexes (i,j). At the first call of
C      PHYTO_GROWTH the parameters AA() and the affinities is set.
C      The Michaelis-Menten limitation for nutrients and light
C      is first found, together with the maximum temparature
C      dependent growth. Finally the phytoplankton production
C      and respiration is computed.
C
C***ROUTINES CALLED :  NOR05_SLIM
C
C***END PROLOGUE
C
C   Global variables:
C
C   On entry :
C
C      I,J    : The indexes of the column where to compute phyto_growth
C      RAD    : Array for the irradiation in the different layers
C      T      : 3D-array for temperature
C      T      : 3D-array for salinity
C      others : see the files globals.doc for description
C
C   On exit :
C
C      PROD_DIA : 3D-array for the production of diatoms
C      PROD_FLA : 3D-array for the production of flagellates
C      PROD_CHA : 3D-array for the production of chatonella
C      RED_DIA  : 3D-array for the respiration of diatoms
C      RED_FLA  : 3D-array for the respiration of flagellates
C      RED_CHA  : 3D-array for the respiration of chatonella
C
      use mod_xc
      use mod_necessary_ecovars
      use m_NOR05_biosun
      use m_NOR05_slim
      use mod_checknan
      implicit none
C
      integer,intent(in) :: i,j,n
      integer k
C
C   Local variables:
C
C     AA(1)    : Diatom production maximum at 0 degree celsius (s^{-1})
C     AA(2)    : Temperature dependent PMAX for the diatoms (C^{-1})
C     AA(3)    : Flagellate production maximum at 0 degree celsius (s^{-1})
C     AA(4)    : Temperature dependent PMAX for the flagellates (C^{-1})
C     AA(5)    : Respiration rate at 0 degrees celcius (s^{-1})
C     AA(6)    : Respiration rate temperature dependence (C^{-1})
C     KN_FLA   : Affinity nitrogen flagellates (uMN)
C     KP_FLA   : Affinity phosphate flagellates (uMP)
C     KR_FLA   : Affinity radiation flagellates (uE/m^2/s)
C     KN_DIA   : Affinity nitrogen diatoms (uMN)
C     KP_DIA   : Affinity phosphate diatoms (uMP)
C     KS_DIA   : Affinity silicate diatoms (uMSi)
C     KR_DIA   : Affinity radiation diatoms (uE/m^2/s)
C     KN_CHA   : Half saturation nitrogen chatonella (uMN)
C     KP_CHA   : Half saturation phosphate chatonella (uMP)
C     KR_CHA   : Half saturation radiation chatonella (uE/m^2/s)
C     RAD_LIM  : Michaelis-Menten limitation on LIGHT
C     PHO_LIM  : Michaelis-Menten limitation on PHOSPHATE
C     SIL_LIM  : Michaelis-Menten limitation on SILICATE
C     NIT_LIM  : Michaelis-Menten limitation on NITRATE
C     NUT_LIM  : Minimum of pho_lim,sil_lim and nit_lim
C
      logical   first
C
      real nut_lim,nit_lim,sil_lim,pho_lim,rad_lim,prod_max
      real oxymin

      include 'biocom.h' 
      include 'common_blocks.h'
 
      oxymin=10.0
C
C First executable statement PHYTO_GROWTH
C
C Compute photosynthetic production and respiration for
C the DIATOMS. There is no respiration when the concentration
C of cellular nitrate is less than some minimum value.
C
c
      do 100 k=1,kdm-1
C
        prod_max=aa(1)*exp(aa(2)*temp(i,j,k,n))
        rad_lim=NOR05_slim(kr_dia,prod_max,rad(i,j,k))
!AS try formulation that is not so limiting (Approx. Christian, 2002)
!AS (1-e(alpha*rad/prod_max): alpha=0.01, prod_max=2.0
!AS        rad_lim=(1-exp(-0.01*rad(i,j,k)))
!AS end changes
        nit_lim=NOR05_slim(kn_dia,prod_max,bio(i,j,k,n,init))
        pho_lim=NOR05_slim(kp_dia,prod_max,bio(i,j,k,n,ipho))
        sil_lim=NOR05_slim(ks_dia,prod_max,bio(i,j,k,n,isil))
C NB!!!!!!!!! RADIATION INCORPORATED IN NUT_LIM NB!!!!!!!!!!!
        nut_lim=min(nit_lim,pho_lim,sil_lim,rad_lim)
CAS        nut_lim=min(nit_lim,pho_lim,sil_lim)
CHS        NUT_LIM=min(NIT_LIM,PHO_LIM,SIL_LIM)
C
        if(bio(i,j,k,n,idia).le.diamin)then
          red_dia(i,j,k) =0.
C NB!!!!!!!!! RADIATION INCORPORATED IN NUT_LIM NB!!!!!!!!!!!
          prod_dia(i,j,k)=prod_max*nut_lim*max(0.,bio(i,j,k,n,idia))
CHS          PROD_DIA(I,J,K)=PROD_MAX*RAD_LIM*NUT_LIM*max(0.,bio(I,J,K))
        elseif(bio(i,j,k,n,ioxy).lt.oxymin)then
           red_dia(i,j,k)=0.
       else
C NB!!!!!!!!! RADIATION INCORPORATED IN NUT_LIM NB!!!!!!!!!!!
          prod_dia(i,j,k)=prod_max*nut_lim*bio(i,j,k,n,idia)
CHS          PROD_DIA(I,J,K)=PROD_MAX*RAD_LIM*NUT_LIM*bio(I,J,K)
          red_dia(i,j,k) =aa(5)*bio(i,j,k,n,idia)*
     &                    exp(aa(6)*temp(i,j,k,n))
       endif
C-----------------------------------------
C Do the same for FLAGELLATES
C-----------------------------------------
        prod_max=aa(3)*exp(aa(4)*temp(i,j,k,n))

        rad_lim=NOR05_slim(kr_fla,prod_max,rad(i,j,k))
        nit_lim=NOR05_slim(kn_fla,prod_max,bio(i,j,k,n,init))
        pho_lim=NOR05_slim(kp_fla,prod_max,bio(i,j,k,n,ipho))
C NB!!!!!!!!! RADIATION INCORPORATED IN NUT_LIM NB!!!!!!!!!!!
        nut_lim=min(nit_lim,pho_lim,rad_lim)
C        NUT_LIM=min(NIT_LIM,PHO_LIM)
C
        if(bio(i,j,k,n,ifla).le.flamin)then
          red_fla(i,j,k) =0.
C NB!!!!!!!!! RADIATION INCORPORATED IN NUT_LIM NB!!!!!!!!!!!
          prod_fla(i,j,k)=prod_max*nut_lim*max(0.,bio(i,j,k,n,ifla))
CHS          PROD_FLA(I,J,K)=PROD_MAX*RAD_LIM*NUT_LIM*max(0.,FLA(I,J,K))
        elseif(bio(i,j,k,n,ioxy).lt.oxymin) then
           red_fla(i,j,k)=0.
       else
C NB!!!!!!!!! RADIATION INCORPORATED IN NUT_LIM NB!!!!!!!!!!!
          prod_fla(i,j,k)=prod_max*nut_lim*bio(i,j,k,n,ifla)
CHS          PROD_FLA(I,J,K)=PROD_MAX*RAD_LIM*NUT_LIM*FLA(I,J,K)
          red_fla(i,j,k) =aa(5)*bio(i,j,k,n,ifla)*
     &                    exp(aa(6)*temp(i,j,k,n))
       endif
C
c      #ifdef CHATONELLA
C-----------------------------------------
C Do the same for CHATONELLA
C-----------------------------------------
c        PROD_MAX = 0.27562160924845 
c     +            -0.00329561997731*temp(i,j,k,n) 
c     +            -0.07123134308130*saln(i,j,k,n) 
c     +            +0.00413658308488*temp(i,j,k,n)*temp(i,j,k,n) 
c     +            +0.00102524245676*temp(i,j,k,n)*saln(i,j,k,n) 
c     +            +0.00483423417638*saln(i,j,k,n)*saln(i,j,k,n) 
c     +            -0.00021627146528*temp(i,j,k,n)*temp(i,j,k,n)
c     &            *temp(i,j,k,n) 
c     +            -0.00010754985072*temp(i,j,k,n)*temp(i,j,k,n)
c     &            *saln(i,j,k,n) 
c     +            +0.00001546340522*temp(i,j,k,n)*saln(i,j,k,n)
c     &            *saln(i,j,k,n) 
c     +            -0.00008742600567*saln(i,j,k,n)*saln(i,j,k,n)
c     &            *saln(i,j,k,n)
c        PROD_MAX = max(0.,PROD_MAX)/86400.
cC
c        RAD_LIM=max(0.,RAD(I,J,K)/(RAD(I,J,K)+KR_CHA))
c        NIT_LIM=max(0.,bio(i,j,k,n,init)/(bio(i,j,k,n,init)+KN_CHA))
c        PHO_LIM=max(0.,bio(i,j,k,n,ipho)/(bio(i,j,k,n,ipho)+KP_CHA))
c        NUT_LIM=min(NIT_LIM,PHO_LIM,RAD_LIM)
C
c        IF(CHA(I,J,K).LE.CHAMIN)THEN
c          RED_CHA(I,J,K) =0.
c          PROD_CHA(I,J,K)=PROD_MAX*NUT_LIM*max(0.,CHA(I,J,K))
c        ELSE
c          PROD_CHA(I,J,K)=PROD_MAX*NUT_LIM*CHA(I,J,K)
c          RED_CHA(I,J,K) =AA(5)*CHA(I,J,K)*EXP(AA(6)*T(I,J,K))
c        ENDIF
C
c      #endif
! save primary production for the record
        grosspp(i,j,k)=prod_dia(i,j,k) + prod_fla(i,j,k);
        netpp(i,j,k)= prod_dia(i,j,k) + prod_fla(i,j,k) 
     &            - red_dia(i,j,k) - red_fla(i,j,k);
C
 100  continue
C
      return
      end subroutine NOR05_phyto_growth
      end module m_NOR05_phyto_growth
