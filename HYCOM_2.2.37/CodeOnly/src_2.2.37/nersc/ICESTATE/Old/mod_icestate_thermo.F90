! ======================================================================
!  This module contains all definitions and parameters required to run 
!  the sea ice redistribution process (through freezing, opening,
!  ridging and rafting).
! ======================================================================
   !Written by: David Salas y Melia  and  Ragnhild Hansen
   !Last changed 10.03.2000 by: Knut Arild Liseter

MODULE mod_icestate_thermo
   implicit none

   ! Albedo values
   REAL, save      :: albi_max = .80      ! Max albedo over dry ice []
   real, save      :: albi_mlt = .60      ! Albedo over melting ice []
   REAL, save      :: albsmin  = .75      ! Minimum albedo over snow []
   REAL, save      :: albsmax  = .85      ! Maximum albedo over snow []
   REAL, PARAMETER :: albw     = .065     ! Albedo over water []
   REAL, PARAMETER :: coalbw= 1.-albw     ! Water SW absorbtion


   !  Specific heats [J/(kg K)]
   REAL, PARAMETER :: cpice    = 2080.    ! ice 
   REAL, PARAMETER :: cpsw     = 3987.    ! water 
   REAL, PARAMETER :: cpw      = 4180.    ! fresh water 


   ! Phase transition heat parameters
   REAL, PARAMETER :: hofusn0  = 3.02E+08  ! Heat of fusion of ice [J/m^3]
   REAL, PARAMETER :: hofusni0= 1./hofusn0 ! Inverse of heat of fusion [m^3 / J]


   ! Ice and snow density [kg/m^3]
   real, parameter :: rhoice   = 910      ! Ice
   REAL, PARAMETER :: rhosnwmax= 330.     ! Maximum density for old snow Maykut (1971)
   !REAL, PARAMETER :: rhosnwmax_sum= 450.! Maximum density for old snow (summer) Maykut (1971)
   REAL, PARAMETER :: rhosnwmin= 150.     ! fresh snow 
   REAL, PARAMETER :: rhow     = 1000.    ! fresh water
   REAL, save      :: rhoref              ! Ocean reference specific density -- from OGCM
   REAL, save      :: thref               ! 1/rhoref [m^3/kg]


   ! Phase transition temperatures
   REAL, PARAMETER :: tice_m   = 273.05   ! Melting point of ice [K]
   REAL, PARAMETER :: t0deg    = 273.15   ! 0 C in [K]


   ! Decay coeffs for snow albedo
   REAL, PARAMETER :: taua     = .008     ! Linear decay coef. for snow albedo []
   REAL, PARAMETER :: tauf     = .24      ! Exponential decay coef. for snow albedo []



   ! Various constants
   real, parameter :: g=9.82              ! g ...
   REAL, PARAMETER :: rkice    = 2.04     ! Ice conductivity [W/(mK)]
   REAL, PARAMETER :: sice     = 6.       ! Salt concentrtation per mil for sea ice []


   ! "Tunable" parameters
   REAL, PARAMETER :: wnew     = 2.E-3    ! Minimum snowfall that allows snow albedo to be refreshed to albsmax [m]
   real, save      :: snwlim  = .3        ! Limiting snow thickness
   real, save      :: qst_frac = 0.3      ! Max brine storage relative to total latent heat of ice
   real, save      :: fice_max = 0.999    ! Max allowed ice conc.
   real, parameter :: fice_min = 0.02     ! Min allowed ice conc.


   !Some useful parameters :
   real,    parameter :: pi2 = 6.283185307178 ! 2 times pi
   real,    parameter :: pi = 3.141592653589  !  pi
   REAL,    PARAMETER :: rmnthi   = 3.86E-07  ! Precipitations: 1m/month in m/s [m/s]
   real,    parameter :: radian=57.2957795


   ! Time specifics
   real, save :: dtt       ! Timestep in seconds for thermo routine
   !real, save :: srelax    ! Salinity surface relaxation timescale
   !real, save :: trelax    ! Temperature surface relaxation timescale

   ! Horizontal melt parameterization flags
   logical, save :: hmelt_hak,hmelt_solar

end module mod_icestate_thermo
