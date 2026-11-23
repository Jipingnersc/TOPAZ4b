      module m_NOR05_zoo_growth
!
      contains
! 
      subroutine NOR05_zoo_growth(k_list,temp,n)

      use mod_xc
      use mod_necessary_ecovars
!      use m_NOR05_biosun
!      use m_NOR05_slim
!      use mod_checknan
      implicit none
      include 'sedcom.h'
      include 'biocom.h'
!      include 'common_blocks.h'

!
      integer,intent(in) :: n
      integer,intent(in),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)::
     &      k_list 
      real,intent(in),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)::
     &      temp
      integer         :: i,j,k
! Run 15:
      real, parameter :: pi11 = 0.450, pi12 = 0.275, pi13 = 0.275
      real, parameter :: pi21 = 0.334, pi23 = 0.333, pi24 = 0.333
!AS      real, parameter :: pi21 = 0.633, pi23 = 0.0, pi24 = 0.367
      real, parameter :: k3 = 1.0, k6 = 0.2
      real            :: g11, g12, g13, g21, g23, g24, zin, zen
      real            :: denum,p11,p12,p13,p21,p23, p24,tfac,t0,tmp
      real            :: beta, mju1,mju2, delta, eps, g
!-------------------
!
! Computing G11, G12, G13 and estimating the terms for mesozooplankton:
!
! p1n_zen : grazing from mesozooplankton on dia
! zin_zen : grazing from mesozooplankton on micro 
! d1n_zen : grazing from mesozooplankton on det
! zen_d1n : fecal pellets prod by mesozoo mortality into det (actually sum of d1n and d2n)
! zen_don : excretion of DON and NH4 by mesozooplankton
!
!-------------------
!initilize arrays:
      p1n_zen=0.0
      zin_zen=0.0
      d1n_zen=0.0
      zen_d1n=0.0
      zen_d1n=0.0
      zen_don=0.0

      t0    = 10.0
      beta  = 0.75
      mju1  = 0.2
      mju2  = 0.2
      delta = 0.4
      delta = 0.5
      eps   = 0.4
      g     = 0.4  !  grazing param mesozooplankton
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j)
         zen = bio(i,j,k,n,imeso)
         tfac = 1.5**((temp(i,j,k)-t0)/t0)
         denum = pi11*bio(i,j,k,n,idia) + pi12*bio(i,j,k,n,imicro) + 
     &                   pi13*bio(i,j,k,n,idet) 
         p11 = pi11*bio(i,j,k,n,idia)/denum
         p12 = pi12*bio(i,j,k,n,imicro)/denum
         p13 = pi13*bio(i,j,k,n,idet)/denum
         tmp =  tfac*g / 
     &       (cnit*k3 + p11*bio(i,j,k,n,idia) + 
     &               p12*bio(i,j,k,n,imicro) +  p13*bio(i,j,k,n,idet) )
!
         G11 = tmp*p11*bio(i,j,k,n,idia)   * zen
         G12 = tmp*p12*bio(i,j,k,n,imicro) * zen
         G13 = tmp*p13*bio(i,j,k,n,idet)   * zen
!
         p1n_zen(i,j,k) = g11 * dayinv
         zin_zen(i,j,k) = g12 * dayinv
         d1n_zen(i,j,k) = g13 * dayinv
         zen_d1n(i,j,k) = ( (delta+eps)
     &                     *mju1*(zen/(zen+cnit*k6))*zen 
     &                     + (1.-beta)*(g11+g12+g13) ) * dayinv
         zen_don(i,j,k) = ( (1.-delta-eps)
     &                     *mju1*(zen/(zen+cnit*k6))*zen ) 
     &                     * dayinv
      end do
      end do
      end do
!-------------------
!
! Computing g21, g23 and estimating the terms for microzooplankton:
!
! p2n_zin : grazing from microzooplankton on fla
! p1n_zin : grazing from microzooplankton on dia
! d1n_zin : grazing from microzooplankton on det
! zin_d1n : fecal pellets prod by microzoo mortality into det (actually sum of d1n and d2n)
! zin_don : excretion of DON and NH4 by mesozooplankton
!
!-------------------
      g     = 1.0  !  grazing param microzooplankton
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j) 
         zin = bio(i,j,k,n,imicro)
         tfac = 1.5**((temp(i,j,k)-t0)/t0)
         denum = pi21*bio(i,j,k,n,ifla) + pi23*bio(i,j,k,n,idia) + 
     &           pi24*bio(i,j,k,n,idet) 
         p21 = pi21*bio(i,j,k,n,ifla)/denum
         p23 = pi23*bio(i,j,k,n,idia)/denum
         p24 = pi24*bio(i,j,k,n,idet)/denum

         tmp = tfac*g / 
     &      (cnit*k3 + p21*bio(i,j,k,n,ifla) + p23*bio(i,j,k,n,ifla) + 
     &      p24*bio(i,j,k,n,idet) )
!
         g21 = tmp*p21*bio(i,j,k,n,ifla) * zin
         g23 = tmp*p23*bio(i,j,k,n,idet) * zin
         g24 = tmp*p24*bio(i,j,k,n,idet) * zin
!
         p2n_zin(i,j,k) = g21 * dayinv
         p1n_zin(i,j,k) = g23 * dayinv
         d1n_zin(i,j,k) = g24 * dayinv
         zin_d1n(i,j,k) = ( (delta+eps)
     &                     *mju2*(zin/(zin+cnit*k6))*zin 
     &                     + (1.-beta)*(g21+g23+g24) ) * dayinv
         zin_don(i,j,k) = ( (1.-delta-eps)
     &                     *mju2*(zin/(zin+cnit*k6))*zin ) 
     &                     * dayinv
      end do
      end do
      end do
!
! Adjust grazing due to low phytoplankton biomass
!
      where(bio(:,:,:,n,idia).lt.diamin)
         p1n_zen = 0.
      end where
!
      where(bio(:,:,:,n,ifla).lt.flamin)
         p1n_zin = 0.
         p2n_zin = 0.
      end where
!
      return
      end subroutine NOR05_zoo_growth
!
      end module m_NOR05_zoo_growth
