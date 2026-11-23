! ======================================================================
!  This module contains the definition of the ice thickness classes
! ======================================================================
   !Written by: David Salas y Melia  and  Ragnhild Hansen
   !Last changed 02.01.2002 by Knut Arild Liseter

MODULE mod_icestate_classes
   use mod_icestate_thermo , ONLY : epsil1
   implicit none

   ! ====================================================================
   ! ====================definition of ice classes ======================
   ! ====================================================================
   ! Vector thickl represents the boundaries of each ice class. For 
   ! example, ice class number i consists of ice whose thickness is 
   ! between thickl(i) and thickl(i+1). The thickest ice has boundaries
   ! thickl(nlay) and +infinity.
   ! Vector nlay defines vertical layers of each ice class

#ifdef ICESTATE_1CLASS
   ! - 1 class
   integer, parameter              :: nthick=1
   INTEGER, PARAMETER              :: hklim =1   ! Number of "thin" ice categories
   REAL, save, DIMENSION(nthick)   :: thickl    
   DATA thickl/ epsil1  /
   integer, parameter ::nlaymax = 1
   integer, parameter, dimension(nthick)   :: nlay=(/ 1  /)       
#endif




#ifdef ICESTATE_2CLASS
   ! - 2 classes
   integer, parameter              :: nthick=2
   INTEGER, PARAMETER              :: hklim =2   ! Number "thin" ice categories
   REAL, save, DIMENSION(nthick)   :: thickl     ! thickl describe the the lower vert.boundaries 
   DATA thickl/ epsil1 ,  1.0 /
   integer, parameter ::nlaymax = 3
   integer, parameter, dimension(nthick)   :: nlay=(/ 1 , 3  /)       
#endif




#ifdef ICESTATE_5CLASS
   ! - 5 classes
   integer, parameter              :: nthick = 5
   INTEGER, PARAMETER              :: hklim  = 2   ! Number of "thin" ice categories
   REAL, save, DIMENSION(nthick)   :: thickl      
   DATA thickl/ epsil1 , .5 , 1.0 , 2.0 , 5.0 /
   integer, parameter ::nlaymax = 5
   integer, parameter, dimension(nthick)   :: nlay=(/ 1 , 1 , 3 , 4 , 5 /)       
#endif




   ! The definition of "thin" and "thick" ice is given in module mod_icestate_init
   logical, save, dimension(nthick)                :: thin       ! True if ice is 'thin'



end module mod_icestate_classes
