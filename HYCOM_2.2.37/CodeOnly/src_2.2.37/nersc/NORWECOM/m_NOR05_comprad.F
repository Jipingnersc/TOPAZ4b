      module m_NOR05_comprad

      contains
  
CKAL  subroutine NOR05_comprad(daynu,theta_nor,maxglo,meanglo,
CKAL &                         sumdf,dti,modlat,modlon)
      subroutine NOR05_comprad(daynu,theta_nor,maxglo,meanglo,
     &                         sumdf,dti)
C***
C***BEGIN PROLOGUE COMPRAD
C***DATE WRITTEN  21/9-92
C***AUTHOR
C      K{\aa}re B. Ulvestad, Department of Fisheries and
C      Marine Biology, Bergen High Technology Centre,
C      5014 Bergen, Norway
C
C***Revised date 4/12-92 Morten D. Skogen, morten@imr.no
C     Added comments and prologue and made some small changes.
C***Revised date 27/8-96 Morten D. Skogen, morten@imr.no
C     MEANRAD and MAXRAD combined to COMPRAD
C***Revised date 1/11-01 Morten D. Skogen, morten@imr.no
C     Radiation data is made 2D field, thus the surface irradiance is also
C     extended to 2D (GLORAD). This means that instead of using parameters
C     LONGIT and LATIT within COMPRAD, these are now 2D variables in the
C     call to the routine. In addition instead of 2 calls (FLAG=1,2)
C     in the previous version, one call is now enough, and GLORAD is 
C     substituted with MAXGLO and MEANGLO, with two separate parameter sets
C     xxMAX and xxMEAN. To save computing time the calls to SUNTIME and 
C     SUNELEV is changed to a 2D version, BIOSUN. Finally, the routine
C     LOWLIGHT has been logically included by computing SUMDF in COMPRAD. 
C
C***PURPOSE COMPRAD computes the totally integrated global
C     surface irradiance at daynymber DAYNR using a set
C     of all constants.
C***PURPOSE LOWLIGHT computes the irradiance when the total diurnal
C     integrated irradiance is lower than the middle value in the
C     respective month. When this is true, we assume that the beam
C     irradiance is equal to zero and that the global irradiance is
C     equal to the diffuse irradiance.
C
C***DESCRIPTION :
C     H(h,n) = I0(n) \times Tr(n) \times F(h),
C        where
C     H(h,n) is the global irradiance at surface,
C     I0(n) is the solar irradiance at normal incident just 
C        outside the atmosphere
C     Tr(n) the transmittance
C     F(h) the solar elevation function,
C     h = sunheight, and n = daynr.
C
C***REFERENCE : Modelling slope irradiance at high latitudes,
C     Olseth, J.A. & Skartveit A., pp 335-337
C
C***ROUTINES CALLED : FX, BIOSUN
C
C***END PROLOGUE
C
C   Global variables:
C
C   On entry :
C     DAYNR : Integer, day number (1..365)
C     THETA_NOR : Real, declination
C     DTI   : Real, internal timestep
C     LATIT : 2D-Real, latitude at points to integrate at
C     LONGIT: 2D-Real, longitude at points to integrate at
C
C   On exit :
C
C     MAXGLO : 2D-Real, total global irradiance using MAX values
C     MEANGLO: 2D-Real, total global irradiance using MEAN values
C     SUMDF  : 2D real, sum of the diffuse irradiance functions
      use mod_xc
      use m_NOR05_biosun
      use m_NOR05_fx
C
      implicit none
C
      integer,intent(in) ::  daynu
      real dti,theta_nor
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)::maxglo,meanglo,
     &                                                 sumdf
CKAL  real,intent(in),dimension(1:itdm,1:jtdm) ::
CKAL &               modlon,modlat 
C
C   Local variables:
C
C     AA,BB,CC,DD,EE,FF : Reals, helping variables for modeling
C         the light. Here : values for Bergen. See reference.
C     DD,EE  : Reals, constants used in the diffuse light comp.
C     GLOTR  : Real, global transmittance
C     I0     : Real, solar constant
C     IRREX  : Real, solar irradiance at normal incident just
C         outside the atmosphere
C     SUNHEI : 2D-Real, sun elevation above horizon
C     SUNTI  : 2D-Real, true solar time
C     TTIME  : Real, time (1..24)
C
      integer i,j
      real pi
      parameter(pi=3.14159265)
      real aamax,bbmax,ccmax,ddmax,eemax,glotrmax
      data aamax/0.821/,bbmax/0.048/,ccmax/24./,ddmax/1.409/,
     +     eemax/-0.470/
      real aamean,bbmean,ccmean,ddmean,eemean,glotrmean
      data aamean/0.5135/,bbmean/0.041/,ccmean/104./,ddmean/1.503/,
     +     eemean/-0.59/
      real ff,rad_top,irrex,t_comp
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)::sunhei,sunti
      real :: daynr

C      
      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point 
      !margin = 0
c --- CH:081105:adding 1 to daynr
      daynr=real(daynu)+1.
C
C***FIRST executable statement COMPRAD
C
C Compute the light transmission considering the light
C formulation between Bergen and Montreal as valid in whole domain
C (Note longitude/latitude input to BIOSUN)
C
      glotrmax=1.+bbmax*cos((daynr-ccmax)/365.*2.*pi)
      glotrmax=aamax*glotrmax
C
      glotrmean=1.+bbmean*cos((daynr-ccmean)/365.*2.*pi)
      glotrmean=aamean*glotrmean
C
C Compute the extraterrestrial global irradiance at normal incidence
C
      rad_top=1367.
      ff=(0.9856*(daynr-3.0))*pi/180.
      ff=cos(ff)
      ff=1.+0.03346*ff
      irrex=rad_top*ff
C
C Find the solar elevation and integrate it, timestep DTI
C
      do 80 j=1-margin,jj+margin
      do 80 i=1-margin,ii+margin
         maxglo(i,j) = 0.
         meanglo(i,j)= 0.
         sumdf(i,j)  = 0.
 80   continue
C
         do 100 t_comp=1.,24.,dti/3600.
CKAL - In MPI mode, plon/plat should be used (modlon, modlat are global arrays),
CKAL - plon/plate are available in common_blocks, so not passed as args
CKAL       call NOR05_biosun(daynu,t_comp,theta_nor,modlon,modlat
CKAL &     ,sunti,sunhei)
           call NOR05_biosun(daynu,t_comp,theta_nor,sunti,sunhei)
         do 90 j=1-margin,jj+margin
         do 90 i=1-margin,ii+margin
         if(sunhei(i,j)*180./pi.ge.5)then
             maxglo(i,j) = maxglo(i,j) + NOR05_fx(ddmax,eemax,
     &        sunhei(i,j))
             meanglo(i,j)= meanglo(i,j) + NOR05_fx(ddmean,eemean,
     &        sunhei(i,j))
         elseif(sunhei(i,j)*180./pi.gt.0)then
             maxglo(i,j) = maxglo(i,j) + 
     &            sunhei(i,j)*180./(5.*pi)*NOR05_fx(ddmax,eemax,
     &             5.*pi/180.)
             meanglo(i,j) = meanglo(i,j) + 
     &            sunhei(i,j)*180./(5.*pi)*NOR05_fx(ddmean,eemean,
     &             5.*pi/180.)
         endif
  90  continue
 100  continue
C
C Returns GLORAD (W/m^2)
C
      do 110 j=1-margin,jj+margin
      do 110 i=1-margin,ii+margin
         sumdf(i,j)  = meanglo(i,j)
         maxglo(i,j) = maxglo(i,j)*irrex*glotrmax*dti
         meanglo(i,j)= meanglo(i,j)*irrex*glotrmean*dti
 110  continue
C
C***END subroutine COMPRAD
C
      return
      end subroutine NOR05_comprad
      end module m_NOR05_comprad
