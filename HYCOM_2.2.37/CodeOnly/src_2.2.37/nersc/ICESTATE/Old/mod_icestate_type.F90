! ======================================================================
! This module contains all of the types needed by other modules
! 
! ======================================================================
!Written by: David Salas y Melia  and  Ragnhild Hansen
!Last changed 10.03.2000 by: Knut Arild Liseter

MODULE mod_icestate_type
   use mod_icestate_classes
   implicit none


   !  - Type describing the state of the ice class
   type t_ice
      integer :: nlay                   ! Number of vtp layers
      real    :: qstore                 ! Energy stored in brine cells
      real    :: albs                   ! Albedo of snow
      real    :: fice                   ! Cell fraction with ice 
      real    :: hice                   ! Average thickness of ice 
      REAL    :: hsnw                   ! Average thickness of snow
      REAL    :: rhosnw                 ! Density of snow
      REAL    :: tsrf                   ! Temperature of surface 
      real, dimension(1:nlaymax) :: vtp ! Vert. temp profile of ice
   end type t_ice


   !  - Type describing the gridcell (ice and ml)
   type t_istate_cell
      type(t_ice),    dimension(nthick) :: ice
      real                              :: hml
      real                              :: sml
      real                              :: tml
   end type t_istate_cell


end module mod_icestate_type
