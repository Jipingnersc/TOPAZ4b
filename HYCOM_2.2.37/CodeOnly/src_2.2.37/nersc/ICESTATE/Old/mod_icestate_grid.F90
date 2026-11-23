! ======================================================================
! Grid definitions for the icestate model. The size of the domain is taken
! from the ocean model here. 
! ======================================================================
!Written by: David Salas y Melia  and  Ragnhild Hansen
!Module layout changed 03.01.2002 by Knut Arild Liseter
!Changed for MPI - commenced 24.09.2004 - Knut Liseter



MODULE mod_icestate_grid
   use mod_xc 
   implicit none

   ! Test parameters
   integer, parameter :: itst=13, jtst=30


   ! Latitude. (Used in solar flux and to deduce coriolis force). 
   REAL,    save, DIMENSION(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rlatm      ! Latitude in radians


end module mod_icestate_grid
