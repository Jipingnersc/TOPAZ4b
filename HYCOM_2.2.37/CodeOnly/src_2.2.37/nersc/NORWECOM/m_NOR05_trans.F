      module m_NOR05_trans

      contains
 
      subroutine NOR05_trans(daynu,beamtr,difftr,
     &                       global,low,maxrad,meanrad)
C
C***
C***BEGIN PROLOGUE TRANS
C***DATE WRITTEN  5/9-92
C***AUTHOR
C      K{\aa}re B. Ulvestad, Department of Fisheries and
C      Marine Biology, Bergen High Technology Centre,
C      5014 Bergen, Norway
C
C***Revised date 16/12-92 Morten D. Skogen, morten@imr.no
C     Added comments and prologue and made some small changes.
C***Revised date 28/8-96 Morten D. Skogen, morten@imr.no
C     Spooling in data file when reading GLOBAL
C***Revised date 22/8-97 Morten D. Skogen, morten@imr.no
C     Reading of data from file (18) moved to main program
C***Revised date 1/11-01 Morten D. Skogen, morten@imr.no
C     Radiation data is made 2D field, thus several of the parameters
C     to TRANS are also 2D from this version (GLOBAL,MAXRAD,MEANRAD,
C     BEAMTR,DIFFTR,LOW). The routine is changed acording to this.
C
C***PURPOSE TRANS computes the beam transmittance and the diffuse
C     transmittance given the day number (1..365).
C
C***REFERENCE : Modelling slope irradiance at high latitudes,
C                Olseth, J.A. & Skartveit A., pp. 335
C
C***ROUTINES CALLED : NONE
C
C***END PROLOGUE
C
C   Global variables:
C
C   On entry :
C
C     DAYNR  : Integer, day number
C     GLOBAL : 2D-Real, global irradiance from data
C     MAXRAD : 2D-Real, global irradiance using clear sky constants
C     MEANRAD: 2D-Real, global irradiance using all average constants
C     
C   On exit :
C
C     BEAMTR : 2D-Real, beam transmittance
C     DIFFTR : 2D-Real, diffuse transmittance
C     LOW    : 2D-Logical, true if GLOBAL < MEANGLO
C
      use mod_xc

      implicit none
      integer,intent(in) :: daynu
      logical,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: low
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: maxrad,meanrad
     &               ,beamtr,difftr,global
C
C   Local variables:
C
C     AA,BB,CC : Reals, coefficients used in the trans computation
C     A/B/C-CLEAR/CLOUD : Reals, clear and cloud sky conditions
C          used in the computation of A, B and C. The clear sky
C          is values for Bergen, the cloud is average between 
C          Bergen and Montreal.
C
      integer   i,j
      real      frac1,frac2,aa,aclear,acloud,bb,bclear,bcloud
      real      cc,cclear,ccloud
      real   :: daynr
C
      include 'common_blocks.h'     
C      
      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point 
      !margin = 0 
c --- CH:081105:adding 1 to daynr
      daynr=real(daynu)+1.
C
C***FIRST executable statement TRANS
C
C
C Compute the transmitance
C
C  If GLOBAL < MEANRAD we don't need the transmitance. The irradiance
C  is computed in a special way in this case, see LOWLIGHT and IRRSUR
C
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
C
      if(global(i,j).lt.meanrad(i,j))then
        low(i,j)=.TRUE.
        goto 101
      else
        low(i,j)=.FALSE.
C
C  If GLOBAL > MAXRAD we compute both diffuse and beam transmitance
C  using clear sky constants.
C
        if(global(i,j).gt.maxrad(i,j))then
          beamtr(i,j) = 1.+0.044*cos((daynr-21)/365.*2.*pi)
          beamtr(i,j) = 0.722*beamtr(i,j)
          difftr(i,j) = 1.+0.052*cos((daynr-64)/365.*2.*pi)
          difftr(i,j) = 0.094*difftr(i,j)
          goto 101
        endif 
C
C  If neither of this MEANRAD < GLOBAL < MAXRAD, we interpolate to find
C  the right values for A, B and C
C
C  First in the BEAM-case
C
        aclear = 0.722
        acloud = 0.2805
        bclear = 0.044
        bcloud = 0.0665
        cclear = 21
        ccloud = 121
C
        frac1 = (global(i,j)-meanrad(i,j))/(maxrad(i,j)-meanrad(i,j))
        frac2 = (maxrad(i,j)-global(i,j))/(maxrad(i,j)-meanrad(i,j))
        aa = frac1*aclear + frac2*acloud
        bb = frac1*bclear + frac2*bcloud
        cc = frac1*cclear + frac2*ccloud
C
        beamtr(i,j) = 1.+bb*cos((daynr-cc)/365.*2.*pi)
        beamtr(i,j) = aa*beamtr(i,j)
C
C Same procedure in the diffuse case
C
        aclear = 0.094
        acloud = 0.232
        bclear = 0.052
        bcloud = 0.03
        cclear = 64
        ccloud = 50
C
        aa = frac1*aclear + frac2*acloud
        bb = frac1*bclear + frac2*bcloud
        cc = frac1*cclear + frac2*ccloud
C
        difftr(i,j) = 1.+bb*cos((daynr-cc)/365.*2.*pi)
        difftr(i,j) = aa*difftr(i,j)
C
      endif
 101  continue
C
      enddo
      enddo 
C
C***END subroutine TRANS
C
      return
      end subroutine NOR05_trans
      end module m_NOR05_trans
