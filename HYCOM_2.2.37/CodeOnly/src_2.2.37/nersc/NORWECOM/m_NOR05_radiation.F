      module m_NOR05_radiation

      contains

      subroutine NOR05_radiation(beamir,diffir,i,j,kappa,phi,kmax,n)
C***
C***BEGIN PROLOGUE RADIATION
C***DATE WRITTEN 1/9-92
C***AUTHOR
C      Morten D. Skogen, Institute of Marine Research
C      Postbox 1870, N-5024 Bergen-Nordnes
C      Email: morten@imr.no
C
C***PURPOSE RADIATION computes the radiation (Einstein/m^2/s) in a
C      column determined by (i,j).
C
C***DESCIPTION the radiation is computed as the sum of the direct
C      and the diffuse light.
C
C      I_dif = PAR x R_dif x exp(-kappa/mu)
C
C      where R_dif is the diffuse component of the surface irradiance,
C      PAR converts from incident diffuse irradiation to photosynthetic
C      active irradiation, Z is depth, MU is mean cosine of the diffuse
C      light and KAPPA is the attenuation coefficient which is a function
C      of the water quality (at input) and the chlorophyll concentration
C                            z
C      kappa = b_2 z + nu x int( DIA + FLA) dz
C                            0
C      NU is the light extinction coefficient.
C      A similar description is given for the direct light substituting
C      MU with COS(PHI), where phi is the zenith angle of the direct light,
C      and R_dif with R_dir, the surface direct irradiance.
C
C***ROUTINES CALLED : NONE
C
C***END PROLOGUE
C
C   Global variables:
C
C   On entry :
C
C      I,J    : The indexes of the water column where to compute the radiation
C      BEAMIR : The surface direct irradiance
C      DIFFIR : The surface diffuse irradiance
C      KAPPA  : 2D-field for the coefficient of light extinction depending
C         on water quality (b2)
C      PHI    : Zenit angle in water of the direct light
C      others : see the files globals.doc for description
C
C   On exit :
C
C      RAD    : 3D-field for the irradiation in the water column (i,j)
      use mod_xc
      use mod_necessary_ecovars
C
      implicit none
C
      integer,intent(in) :: i,j,kmax,n
      real beamir,diffir,phi
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: kappa
C
C   Local variables:
C
C      DMP   : Light attenuation (from the integral) at (S,T)-points
C      DDMP  : Light attenuation temporary variable
C      DELTAZ: Real physical thickness of a layer
C      COSPHI: Real, cos(phi)
C
      real dmp,ddmp,cosphi
      real,dimension(kdm) :: thickn_box 
      integer   k
C
      include 'common_blocks.h'
      include 'biocom.h'
C First executable statement RADIATION
C
C Compute the depth of each cell
C
      do k=1,kdm
        thickn_box(k) = dp(i,j,k,n)/onem
      enddo
C
C Compute damping on the boundaries and in the middle of each cell.
C Then compute the irradiance in each cell by adding the diffuse
C and the direct contributions
C
      cosphi = cos(phi)
C
      dmp = 0.
      do 101 k=1,kmax
        ddmp = (kappa(i,j)+ny_nor/n2chla*(bio(i,j,k,n,idia)
     &             +bio(i,j,k,n,ifla)))
     +             *0.5*thickn_box(K)
        dmp = dmp + ddmp
        rad(i,j,k)=beamir*exp(-dmp/cosphi)
        rad(i,j,k)=rad(i,j,k)+diffir*exp(-dmp/my_nor)
        rad(i,j,k)=par*rad(i,j,k)
        dmp = dmp + ddmp
C
 101  continue
      rad(i,j,kmax:kdm)=0.0
C
      return
      end subroutine NOR05_radiation
      end module m_NOR05_radiation 
