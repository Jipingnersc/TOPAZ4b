      module m_NOR05_wqual

      contains

      subroutine NOR05_wqual(KAPPA)
!***
!***BEGIN PROLOGUE WQUAL
!***AUTHOR
!      K{\aa}re B. Ulvestad, Department of Fisheries and
!      Marine Biology, Bergen High Technology Centre,
!      5014 Bergen, Norway
!
!***Revised date 26/8-92 Morten D. Skogen, morten@imr.no
!     Added prologue
!
!***PURPOSE WQUAL computes the extinction coefficient in the water
!     column as a function of water type.
!
!***DESCRIPTION the water type is defined in terms of the salinity
!     at the surface of the water column.
!
!***ROUTINES CALLED : NONE
!
!***END PROLOGUE
!
!   Global variables:
!
!     S      : 3D-field, salinity
!     others : see the files globals.doc 
!
!   On exit :
!
!     KAPPA : 2D-field, the extinction coefficient in each watercolumn
      use mod_xc
!
      implicit none
      INCLUDE 'biocom.h'
!
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: KAPPA
!
!   Local variables
!
      INTEGER I,J
!
! First executable statement WQUAL
!
      DO 100 j=1-margin,jj+margin
      DO 100 i=1-margin,ii+margin
!       IF(depths(I,J).GT.0.5)THEN
!           IF(saln(I,J,1,n).LT.31)THEN
              KAPPA(I,J)=b2
!           ELSEIF(saln(I,J,1,n).GE.31.AND.saln(I,J,1,n).LT.32)THEN
!             KAPPA(I,J)=b2
!           ELSEIF(saln(I,J,1,n).GE.32.AND.saln(I,J,1,n).LT.33)THEN
!             KAPPA(I,J)=b2
!           ELSEIF(saln(I,J,1,n).GE.33.AND.saln(I,J,1,n).LT.34)THEN
!             KAPPA(I,J)=b2
!           ELSEIF(saln(I,J,1,n).GE.34.AND.saln(I,J,1,n).LT.35)THEN
!             KAPPA(I,J)=b2
!           ELSEIF(saln(I,J,1,n).GE.35)THEN
!             KAPPA(I,J)=b2
!           ENDIF
!       ENDIF
 100  CONTINUE
!
      RETURN
      end subroutine NOR05_wqual
      end module m_NOR05_wqual
