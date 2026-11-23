      module m_NOR05_irrsur

      contains
 
      subroutine NOR05_irrsur(beamtr,difftr,daynu,diffir,beamir,
     &                        global,low,
     &                        maxglo,meanglo,sumdf,sunhei,dti)
C
C***
C***BEGIN PROLOGUE IRRSUR
C***DATE WRITTEN  18/9-92
C***AUTHOR
C      K{\aa}re B. Ulvestad, Department of Fisheries and
C      Marine Biology, Bergen High Technology Centre,
C      5014 Bergen, Norway
C
C***Revised date 7/12-92 Morten D. Skogen, morten@imr.no
C     Added comments and prologue and made some small changes.
C***Revised date 24/3-93 Morten D. Skogen, morten@imr.no
C     Extended to 2D, i.e. the irradiance is computed in all
C     surface points with one call.
C***Revised date 1/11-01 Morten D. Skogen, morten@imr.no
C     Radiation data is made 2D field, thus several of the input parameters
C     to IRRSUR are also 2D from this version (BEAMTR,DIFFTR,GLOBAL,
C     MAXGLO,MEANGLO,LOW). The routine is changed acording to this.
C
C***PURPOSE IRRSUR computes the surface irradiance given the
C     number of day (1..365) and solar altitude above the horizon.
C
C***DESCRIPTION :
C     H(h,n) = sunR(n) \times Tr(n) \times F(h),
C        where
C     H(h,n) is the global irradiance at surface,
C     sunR(n) is the solar irradiance at normal incident just
C        outside the atmosphere
C     Tr(n) the transmittance
C     F(h) the solar elevation function,
C        further
C     h = sunheight, and n = daynr.
C
C***REFERENCE : Modelling slope irradiance at high latitudes,
C     Olseth, J.A. & Skartveit A., pp 335-337
C
C***ROUTINES CALLED : NOR05_fx 
C
C***END PROLOGUE
C
C   Global variables:
C
C   On entry :
C
C     DTI    : Real, internal time step
C     DAYNR  : Integer, Julian day (1..365)
C     BEAMTR : 2D-Real, beam transmitance
C     DIFFTR : 2D-Real, difuse transmittance
C     GLOBAL : 2D-Real, global irradiance at DAYNR
C     LOW    : 2D-Logical, true if GLOBAL is less than MEANGLO
C     MAXGLO : 2D-Real, global irradiance using clear sky constants
C     MEANGLO: 2D-Real, global irradiance using all average constants
C     SUMDF  : 2D-real, sum of solar elevation func. if LOW
C     SUNHEI : 2D-real, sunheight
C
C   On exit :
C
C     BEAMIR : 2D-real, direct surface irradiance
C     DIFFIR : 2D-real, diffuse surface irradiance
      use mod_xc
      use m_NOR05_fx
C
      implicit none
C
      logical,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: low
      integer,intent(in) ::   daynu
      real      dti
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: beamir,diffir,
     &               beamtr,difftr,
     &          sunhei,global,maxglo,
     &          meanglo,sumdf
C
C   Local variables:
C
C     AA_isur,DD,EE : Reals, variabels in the irradiance computation
C     DCLEAR,DCLOUD,ECLEAR,ECLOUD : Reals, parameters used to find DD and EE
C     sunR    : Real, solar constant
C     IRREX : Real, solar irradiance at normal incident just outside
C               the atmosphere
C
      integer   i,j
      real :: daynr
      real      frac1,frac2,aa_isur,sunr,irrex
      real      dd,dclear,dcloud,ee,eclear,ecloud
      parameter(sunr=1367.)
      
      include 'common_blocks.h'
C
      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point 
      !margin=0
c --- CH:081105:adding 1 to daynr
      daynr=real(daynu)+1.

C***FIRST executable statement IRRSUR
C
C  If the irradiance is lower than the mean irradiance for the
C  current day, we compute the irradiance by (8) in Skartveit & Olseth.
C  SUMDF is computed in subroutine LOWLIGHT and equals the sum of Fx.
C  We multiply by DTI to get the total integrated one.
C  In this case the direct irradiance equals 0.
C
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
C
      if(low(i,j))then
        dd = 1.503
        ee = -0.59
        diffir(i,j) = global(i,j)/(dti*sumdf(i,j))*
     +        NOR05_fx(dd,ee,max(sunhei(i,j),5.*pi/180.))
C
         beamir(i,j) = 0.
         if(sunhei(i,j).le.0.)then
            diffir(i,j)=0.
         else if(sunhei(i,j)*180./pi.lt.5.)then
            diffir(i,j) = diffir(i,j)*sunhei(i,j)*180./(5.*pi)
         endif 
C
         goto 100
       endif 
C
C  Compute solar irradiance at normal incidence just outside
C  the atmosphere (solar constant = 1367W/m**2).
C  Reference: Meterological report series, University of Bergen 7-1988
C
      aa_isur=(0.9856*(daynr-3.0))*pi/180.
      aa_isur=cos(aa_isur)
      aa_isur=1.+0.03346*aa_isur
      irrex=sunR*aa_isur
C
C  First we compute irradiance in the BEAM case,
C
C  Initiate the constants DD and EE, the clear sky constants are
C  values for Bergen, the cloud sky are average for Bergen/Montreal
C
      dclear=1.643
      dcloud=1.7515
      eclear=-0.748
      ecloud=-0.8755
C
C  If the global irradiance is greater than modelled light using
C  clear sky constants, we use the clear sky constants.
C  If not, we interpolate to get new constants, DD and EE
C
      if(global(i,j).gt.maxglo(i,j))then
        dd = dclear
        ee = eclear
      else
        frac1 = (global(i,j)-meanglo(i,j))/(maxglo(i,j)-meanglo(i,j))
        frac2 = (maxglo(i,j)-global(i,j))/(maxglo(i,j)-meanglo(i,j))
        dd = frac1*dclear + frac2*dcloud
        ee = frac1*eclear + frac2*ecloud
      endif
C
C  Compute the direct irradiance
C
      beamir(i,j) = irrex*beamtr(i,j)*
     +       NOR05_fx(dd,ee,max(sunhei(i,j),5.*pi/180.))
C
      if(sunhei(i,j).le.0.)then
         beamir(i,j) = 0.
      else if(sunhei(i,j)*180./pi.lt.5)then
         beamir(i,j) = beamir(i,j)*sunhei(i,j)*180./(5.*pi)
      endif 
C
C  Same procedure for the DIFFUSE light
C
      dclear=-0.432
      dcloud=1.1895
      eclear=1.718
      ecloud=-0.2255
C
      if(global(i,j).gt.maxglo(i,j))then
        dd = dclear
        ee = eclear
      else
        frac1 = (global(i,j)-meanglo(i,j))/(maxglo(i,j)-meanglo(i,j))
        frac2 = (maxglo(i,j)-global(i,j))/(maxglo(i,j)-meanglo(i,j))
        dd = frac1*dclear + frac2*dcloud
        ee = frac1*eclear + frac2*ecloud
      endif
C
      diffir(i,j) = irrex*difftr(i,j)*
     +       NOR05_fx(dd,ee,max(sunhei(i,j),5.*pi/180.))
C
      if(sunhei(i,j).le.0.)then
         diffir(i,j) = 0.
      else if(sunhei(i,j)*180./pi.lt.5)then
         diffir(i,j) = diffir(i,j)*sunhei(i,j)*180./(5.*pi)
      endif 
C
C Transform from J/m2 to muEinstein/m2/s
C
 100  continue
C
      beamir(i,j)=beamir(i,j)/0.217
      diffir(i,j)=diffir(i,j)/0.217
C
      enddo 
      enddo 
C
C***END subroutine IRRSUR
C
      return
      end subroutine NOR05_irrsur
      end module m_NOR05_irrsur
