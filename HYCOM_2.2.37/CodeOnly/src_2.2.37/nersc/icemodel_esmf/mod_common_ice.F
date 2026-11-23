      module mod_common_ice
      use mod_xc
      implicit none


! --- common blocks for the ice and snow part of the model
      common /icesnw1/                                                  &
     & albi_m                                                           &! --                 max albedo over melting ice
     &,albi_d                                                           &! --                 max albedo over dry ice
     &,albi                                                             &! --                 albedo over ice
     &,albs                                                             &! --                 albedo over snow
     &,albs_m                                                           &! --                 albedo over melting snow
     &,albw_d                                                           &! --                 albedo for diff light over water
     &,albw                                                             &! --                 albedo for direct light over water
     &,rhoice                                                           &! kg / m^3           density of ice
     &,rhosnw                                                           &! kg / m^3           density of snow
     &,rhowat                                                           &! kg / m^3           density of pure water
     &,emiss                                                            &! --                 emissivity of water
     &,gasconst                                                         &! pa m^3 / (k kg)    gas constant
     &,rkice                                                            &! w / (m k)          ice conductivity
     &,hocond                                                           &! j / kg             heat of condensation/vap
     &,fusi                                                             &! j / kg             heat of fusion of ice
     &,fuss                                                             &! j / kg             heat of fusion of snow
     &,fice_max                                                         &! --                 maximum fractional ice cover
     &,tice_m                                                           &! k                  melting point of ice
     &,tsnw_m                                                           &! k                  melting point of snow
     &,hice_min                                                         &! m                  minimum ice thickness
     &,epsmol                                                           &! --                 molecular weight of h2o/dry air
     &,albsa                                                             ! --                 array containing snow albedo       
      save/icesnw1/

      common /icesnw2/                                                  &
     & sice                                                             &! per mil            salinity of seaice
     &,rksnw                                                            &! w / (m k)          snow conductivity
     &,cpair                                                            &! j / (k kg)         specific heat of dry air
     &,cpsw                                                             &! j / (k kg)         specific heat of seawater
     &,stefanb                                                          &! w / (m^2 k^4)      stefan-boltzman constant
     &,aice                                                             &! --                 vapor pressure parameters
     &,bice                                                             &! k                    ..
     &,awater                                                           &! --                   ..
     &,bwater                                                           &! k                    ..
     &,t0deg                                                             ! k                  zero deg celcius in k
      save/icesnw2/

      real albi_m,albi,albs,albs_m,albw,albw_d, albi_d,                 &
     &     rhoice,rhosnw,rhowat,emiss,gasconst,rkice,hocond,            &
     &     fusi,fuss,fice_max,tice_m,tsnw_m,hice_min,epsmol,sice,       &
     &     rksnw,cpair,cpsw,stefanb,aice,bice,awater,bwater,t0deg       &
     &    ,albsa(12)

      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
     &  ficem,                                                          &
     &  hicem,                                                          &
     &  tsrfm,                                                          &
     &  hsnwm,                                                          &
     &  ticem,                                                          &
     &  qbrine,                                                         &
     &  cawdir,                                                         &
     &  pemnp,                                                          &
     &  clat,                                                           &
     &  radfl0,                                                         &
     &  tauxice,                                                        &
     &  tauyice,                                                        &
     &  qfrz,                                                           &
     &  iceU,                                                           &
     &  iceV,                                                           &
     &  delta_icevol,                                                   &
     &  delta_snwvol,                                                   &
     &  surf_albedo_sum,                                                &
     &  surf_qsw_sum                                                    
        

      common/iceparam/  ficem,hicem,tsrfm,hsnwm,ticem,                  &
     & qbrine,cawdir,pemnp,clat,radfl0,tauxice,tauyice,                 &
     & qfrz,iceU,iceV,delta_icevol,delta_snwvol,                        &
     & surf_albedo_sum,surf_qsw_sum
      save/iceparam/
#if defined (ICE_NEST) && defined (ICE)
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::             &
     &  uicenest,                                                       &
     &  vicenest,                                                       &
     &  hicenest,                                                       &
     &  ficenest
#endif
#if defined (TEST_ICE_AGE)
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
     &  fy_frac ,                                                       &
     &  fy_age,                                                         &
     &  rdg_frac
#endif

      contains

         
      subroutine icedat
      implicit none
      albi_m   =.60      ! --                 max albedo over melting ice
      albi_d   =.73      ! --                 max albedo over dry ice
      albs     =.80      ! --                 albedo over snow
      albs_m   =.71      ! --                 albedo over melting snow
      albw_d   =.065     ! --                 albedo over water; diff light
      rhoice   =900.     ! kg / m^3           density of ice
      rhosnw   =330.     ! kg / m^3           density of snow
      rhowat   =1000.    ! kg / m^3           density of pure water 
      emiss    =.97      ! --                 emissivity of water
      gasconst =.287e3   ! pa m^3 / (k kg)    gas constant
      rkice    =2.04     ! w / (m k)          ice conductivity
      hocond   =2.5e6    ! j / kg             heat of condensation/vap
      fusi     =3.02e8   ! j / m^3            heat of fusion of ice
      fuss     =1.10e8   ! j / m^3            heat of fusion of snow
      fice_max =.995     ! --                 maximum fractional ice cover
      tice_m   =273.05   ! k                  melting point of ice
      tsnw_m   =273.15   ! k                  melting point of snow
      hice_min =.1       ! m                  minimum ice thickness
      epsmol   =.622     ! --                 molecular weight of h2o/dry air
      !albsa    =12*.75 
      albsa    =albs     ! Seasonal albedo -- set to fixed value albs above 

      sice     = 6.      ! per mil            salinity of seaice
      rksnw    =.31      ! w / (m k)          snow conductivity
      cpair    =1004.    ! j / (k kg)         specific heat of dry air
      cpsw     =3987.    ! j / (k kg)         specific heat of seawater
      stefanb  =5.67e-8  ! w / (m^2 k^4)      stefan-boltzman constant
      aice     =9.5      ! --                 vapor pressure parameters
      bice     =7.66     ! k                    ..
      awater   =7.5      ! --                   ..
      bwater   =35.86    ! k                    ..
      t0deg    =273.15   ! k                  zero deg celcius in k
      end subroutine icedat
C
C
C
      subroutine iceinit
      use mod_forcing_nersc
      implicit none
      integer :: i,j
      real    :: flagi1,flagi2,tfrz, flagi, flaghemi
      include 'common_blocks.h'

      if (mnproc==1) then
         write(lp,*) 'Calling iceinit '
      end if
!$OMP PARALLEL DO PRIVATE(j,i,flagi,flaghemi,tfrz) 
!$OMP&      SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         tfrz=273.216-.057*saln(i,j,1,1)
         tfrz=tfrz-t0deg
c
         ! 1 if temperature < freezing point + .3
         ! 0 otherwise
         flagi=(1.+sign(1.,tfrz+.3 -temp(i,j,1,1)))/2.
c
         ! 0  in -50 < lat < 50
         flagi=flagi*(1.+sign(1.,abs(plat(i,j))-50.))/2.
c
         ! Hemisphere flag, 1 if plat > 0
         flaghemi=(1.+sign(1.,plat(i,j)))/2.
c
         ! Initial ice thickness - different for N/S hemispheres
         hicem(i,j)=flagi*(flaghemi*2.5 + (1.-flaghemi)*1.0)
c
         ficem(i,j)=flagi*.95
         hsnwm(i,j)=0.
         ticem(i,j)=temp(i,j,1,1)-flagi*4.+273.15
         tsrfm(i,j)=ticem(i,j)
         qfrz(i,j)=0.0
         pemnp(i,j)=0.
#if defined (TEST_ICE_AGE)
         fy_age  (i,j)=0.
         fy_frac (i,j)=0.
         rdg_frac(i,j)=0.
#endif
      enddo
      enddo
!$OMP END PARALLEL DO
c
      tauxice=0.
      tauyice=0.
      relsurflx=0.
      relsalflx=0.
      end subroutine iceinit
c
      end module mod_common_ice
