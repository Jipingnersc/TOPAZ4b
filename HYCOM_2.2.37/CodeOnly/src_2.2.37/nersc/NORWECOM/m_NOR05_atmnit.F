      module m_NOR05_atmnit

      contains

      subroutine NOR05_atmnit(dti,n)
C***BEGIN PROLOGUE ATMNIT
C***DATE WRITTEN 27/3-92
C***AUTHOR
C            Jarle Berntsen, Institute of Marine Research,
C            Postboks 1870,
C            N-5024 Bergen-Nordnes, Norway
C            Email..  jarleb@imr.no
C
C***PURPOSE ATMNIT places atmospheric deposition of
C           nitrogen in the upper cell of sea cells.
C           The deposition is greatest in coastal cells.
C
C   The input is specified in the file atmos.w51 in
C   Guidance Document for the NSTF
C   modelling workshop 6-8 May, 1992, The Hague.
C
C***DESCRIPTION
C
C   ON ENTRY
C
C   NIT    Real array of dimension (IM,JM,KB).
C          Nitrogen before adding deposition.
C
C   ON EXIT
C   
C   NIT    The atmospheric deposition is added in the
C          surface cell.
C
C***ROUTINES CALLED-INTPOLYG
C***END PROLOGUE ATMNIT
C
      use mod_xc
      use mod_necessary_ecovars
C     
      implicit none
C
      include 'common_blocks.h'
C   Global variables.
C
      integer,intent(in):: n
      real,intent(in)   :: dti
C
C   Local variables.
C
C     CSEA   Atmospheric deposition in off coastal cells.
C     CCOAST  Atmospheric deposition in coastal cells.
C            The unit is mg Nitrogen /m**2/s.
C
      integer i,j
      real csea,ccoast,cons,antro
      data ccoast/6.12E-5/
      data csea/3.14E-5/
C
      antro = 1.0
      !margin=0
      !KAL -- Should stick with existing margin
C
C   We loop over all internal grid cells.
C
      do 200 j = 1-margin,jj+margin
      do 200 i = 1-margin,ii+margin
      if (depths(i,j).gt.0.) then
C
C   Check whether we have a coastal sea cell or not.
C
      if(ip(I-1,J-1)*ip(I-1,J)*ip(I-1,J+1)*
     +           ip(I,J-1)*ip(I,J+1)*ip(I+1,J-1)*
     +           ip(I+1,J)*ip(I+1,J+1).lt.0.5)then
             cons = antro*ccoast
      else
             cons = antro*csea
      endif
C
C   Update the Nitrogen content of the surface layer.
C
           bio(i,j,1,n,init) = bio(i,j,1,n,init)+DTI*CONS/(dp(i,j,1,n)
     &      *depths(i,j))
         END IF
200   CONTINUE
C
      RETURN

      end subroutine NOR05_atmnit
      end module m_NOR05_atmnit
