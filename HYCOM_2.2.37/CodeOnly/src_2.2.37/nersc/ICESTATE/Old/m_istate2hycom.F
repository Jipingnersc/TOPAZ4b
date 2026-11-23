      module m_istate2hycom
      contains
   
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
      end module m_istate2hycom
