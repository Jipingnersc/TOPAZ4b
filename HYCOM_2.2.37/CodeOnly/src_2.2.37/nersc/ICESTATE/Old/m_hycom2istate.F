      module mod_icestate_transfer
      contains

      subroutine hycom2istate(m,n)
      use mod_icestate , only : icestate
      use mod_icestate_fluxes 
      use mod_common_ice
      use mod_forcing_nersc
      implicit none

      integer, intent(in) :: m,n

      integer i,j, patch,l
      real :: tmlfac,smlfac

      real, dimension(itdm,jtdm) :: gdp1, gsaln1, gtemp1
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: work2d

      include 'common_blocks.h'


c --- slp -- mbar
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Islp  (i,j)=slp(i,j,l0)*w0+slp(i,j,l1)*w1
     &              +slp(i,j,l2)*w2+slp(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- radfl= net       radiative thermal flux (w/m**2) positive into ocean
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iradflx(i,j)=radflx(i,j,l0)*w0+radflx(i,j,l1)*w1
     &               +radflx(i,j,l2)*w2+radflx(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- swfl = shortwave radiative thermal flux (w/m**2) positive into ocean
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iswflx(i,j)=swflx (i,j,l0)*w0+swflx (i,j,l1)*w1
     &              +swflx (i,j,l2)*w2+swflx (i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- airt = air temperature (C)
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iairtmp(i,j)=airtmp(i,j,l0)*w0+airtmp(i,j,l1)*w1
     &               +airtmp(i,j,l2)*w2+airtmp(i,j,l3)*w3
         Iairtmp(i,j)=Iairtmp(i,j)+t0deg ! Must be Kelvin in icestate
      end do
      end do
!$OMP END PARALLEL DO
c ---   vpmx = water vapor mixing ratio (kg/kg)
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Ivapmix(i,j)=vapmix(i,j,l0)*w0+vapmix(i,j,l1)*w1
     &               +vapmix(i,j,l2)*w2+vapmix(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- prcp = precipitation (m/sec; positive into ocean)
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iprecip(i,j)=precip(i,j,l0)*w0+precip(i,j,l1)*w1
     &               +precip(i,j,l2)*w2+precip(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- relhum = relhum ()
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Irelhum(i,j)=relhum(i,j,l0)*w0+relhum(i,j,l1)*w1
     &               +relhum(i,j,l2)*w2+relhum(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- wind = windspeed ()
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iwndspd(i,j)=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
     &               +wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- clouds 
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iclouds(i,j)=clouds(i,j,l0)*w0+clouds(i,j,l1)*w1
     &               +clouds(i,j,l2)*w2+clouds(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- taux 
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Itaux(i,j)=taux(i,j,l0)*w0+taux(i,j,l1)*w1
     &             +taux(i,j,l2)*w2+taux(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- tauy 
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Itauy(i,j)=tauy(i,j,l0)*w0+tauy(i,j,l1)*w1
     &              +tauy(i,j,l2)*w2+tauy(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- tauxice 
      Itauxice=tauxice
c --- tauyice 
      Itauyice=tauyice

C$OMP PARALLEL DO PRIVATE(i,j,l,patch) 
C$OMP&SCHEDULE(STATIC,1)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         !icestate(i,j)%hml = dp  (i,j,1,n)  / onem
         !icestate(i,j)%sml = saln(i,j,1,n)
         !icestate(i,j)%tml = temp(i,j,1,n) + t0deg
         ! Avoids decoupling issues
         icestate(i,j)%hml = .5*(dp  (i,j,1,n) + dp  (i,j,1,m))  / onem
         icestate(i,j)%sml = .5*(saln(i,j,1,n) + saln(i,j,1,m))
         icestate(i,j)%tml = .5*(temp(i,j,1,n) + temp(i,j,1,m)) + t0deg
      enddo
      enddo
COMP END PARALLEL DO

      end subroutine hycom2istate


      subroutine istate2hycom
      use mod_icestate        , ONLY :  icestate
      use mod_icestate_fluxes
      use mod_forcing_nersc
      implicit none

      !TODO: cekman not really needed anymore
      !real, parameter :: cekman=0.7
      integer i,j,l


      include 'common_blocks.h'


!$OMP PARALLEL DO PRIVATE(j,i) 
!$OMP&SCHEDULE (STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         ! Salinity flux
         salflx(i,j) = 0.

         ! Heat fluxes
         sswflx(i,j) = 0.
         surflx(i,j) = 0.

         ! Corio from HYCOM
         ustar (i,j) = 0.
         hekman(i,j) = 0.

      end do
      end do
!OMP END PARALLEL DO



!$OMP PARALLEL DO PRIVATE(j,l,i) 
!$OMP&SCHEDULE (STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      
            ! Salinity flux
            salflx(i,j) = Isalflx(i,j) 

            ! Heat fluxes
            sswflx(i,j) = Isswflx(i,j)
            surflx(i,j) = Isurflx(i,j)

            ! Corio from HYCOM
            ustar (i,j) = Iustar (i,j)
            hekman(i,j) = cekman*ustar(i,j)/ abs(corio(i,j))

      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
!      print *,'KAL:Istate2hycom sswflx :',mnproc,
!     &   maxval(sswflx),minval(sswflx)
!      print *,'KAL:Istate2hycom sswflx :',mnproc,
!     &   maxval(Isswflx),minval(Isswflx)
!      print *,'KAL:Istate2hycom surflx :',mnproc,
!     &   maxval(surflx),minval(surflx)

!      if (mnproc==2) print *,'point ssw',i0+52,j0+82,52,82,
!     &   sswflx(52,82),Isswflx(52,82),radfl_day(52,82)

      end subroutine istate2hycom


      module m_evp2istate
   
      subroutine evp2istate
      use mod_xc
      use mod_icestate_fluxes , ONLY : ISTATE_divu=>divu
      use mod_evp          ,    ONLY : EVP_divu=>divu
      implicit none
      ! No conversion factor needed. Both ICESTATE and EVP operate with SI units
      
      call xcaput(EVP_divu,ISTATE_divu,0)

      end subroutine evp2istate



      end module mod_icestate_transfer
