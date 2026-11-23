module m_istate2evp
contains

   subroutine istate2evp
      use mod_icestate_classes 
      use mod_icestate 
      use mod_evp          , ONLY : ilo,ihi,jlo,jhi, ncat, hice, hsnow, compact, &
                                    rhosnw, albs, vtp, nlaymax,tsrf
      use mod_evp_tools    , ONLY : boundary
      implicit none

      integer :: si, sj  !  ICESTATE indices
      integer :: ei, ej  !  EVP indices
      integer :: k, l


      do k=1,nthick
         call xcaget(hsnow  (ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%hsnw  ,0)
         call xcaget(compact(ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%fice  ,0)
         call xcaget(hice   (ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%hice  ,0)
         call xcaget(rhosnw (ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%rhosnw,0)
         call xcaget(albs   (ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%albs  ,0)
         call xcaget(tsrf   (ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%tsrf  ,0)
         do l=1,nlaymax
            call xcaget(vtp    (ilo:ihi,jlo:jhi,l,k),icestate(:,:)%ice(k)%vtp(l),0)
         end do
      end do



      ! Set boundary values
      do k=1,nthick
         call boundary( hsnow  (:,:,k)   )
         call boundary( hice   (:,:,k)   )
         call boundary( compact(:,:,k)   )
         call boundary( rhosnw (:,:,k)   )
         call boundary( albs   (:,:,k)   )
         do l=1,nlaymax
            call boundary( vtp (:,:,l,k) )
         end do
      end do

   end subroutine istate2evp

end module m_istate2evp
