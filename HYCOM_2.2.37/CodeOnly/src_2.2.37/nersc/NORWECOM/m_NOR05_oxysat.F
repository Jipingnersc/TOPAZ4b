      module m_NOR05_oxysat 

      contains
 
      subroutine NOR05_oxysat(n)!AS,rr,salt,tem)
!
!***BEGIN PROLOGUE OXYSAT
!***DATE WRITTEN   29/1-96
!***REVISON DATE 31/10-01, morten@imr.no
!      Corrected error in density computation
!***AUTHOR
!     HENRIK SOILAND. henrik@imr.no
!       Institute of Marine Research,
!       Postboks 1870, N-5024 Bergen-Nordnes, Norway
!***PURPOSE OXYSAT.
!     Calculate the concentration of oxygen (mg/l) in one layer
!     assuming saturation in contact with one normal atosphere.
!***REFERENCE: Kester, D.R.,1975. Dissolved gases other than CO2.
!              In: Riley, J.P. and G. Skirrow (eds), Chemical Oceanography
!              2.nd ed., Vol. 1, Academic Press, New York.
!              Formula for micro mol/kg is used,and it is converted to mg/l.
!***END PROLOGUE OXYSAT 
! 
! Global variables
! 
!     Input : RR,SALT,TEMP is the Density, Salinity and Temperature in a layer
!     Exit  : OX       is the saturation at this layer
      use mod_xc
      use mod_necessary_ecovars
      use mod_common_ice, only: hicem

      implicit none

      include 'sedcom.h'
      include 'biocom.h'
      include 'common_blocks.h'


!AS      real,intent(in) ,dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy)  :: 
!AS     &               rr,salt,tem
!AS      real,intent(out),dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) :: ox
      integer,intent(in) :: n
!
!   LOCAL VARIABLES
!
      real :: a1,a2,a3,a4,b1_oxy,b3,b4,tkel,ts,ss,c1,oxya,oxyb,dens
      integer :: i,j
      data a1/-173.9894/,a2/255.5907/,a3/146.4813/,a4/-22.2040/
      data b1_oxy/-0.037376/,b4/0.016504/,b3/-0.0020564/,tkel/273.15/

      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point 
      !margin=0
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
        if (hicem(i,j).eq.0.0) then
         ts=temp(i,j,1,n)+tkel
         ss=saln(i,j,1,n)
         c1=a1+a2*100./ts+a3*alog(ts/100.)+a4*ts/100. 
     &        + ss*(b1_oxy+b4*ts/100.+b3*ts*ts/10000.)
! OXY1 micro mol/kg
         oxya=exp(c1)
! OXY2 micro mol/m3
!AS         dens = 1000. + 1000*(rr(i,j)+.025)
         dens = 1025. + th3d(i,j,1,n)
         oxyb=oxya*dens*ip(i,j)
! mg/m3
!AS         ox(i,j)=coxy2*oxyb*1.0E-3
         bio(i,j,1,n,ioxy)=coxy2*oxyb*1.0E-3
        end if
      enddo
      enddo
!
      return
      end subroutine NOR05_oxysat
      end module m_NOR05_oxysat
