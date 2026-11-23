module m_evp2istate
contains
   
   subroutine evp2istate
      use mod_icestate         
      use mod_icestate_fluxes , ONLY : ISTATE_divu=>divu
      use mod_evp          ,    ONLY : ilo, ihi, jlo, jhi, ncat, hice, compact, hsnow, &
                                       EVP_divu=>divu, rhosnw,albs,vtp,tsrf
      implicit none
      ! No conversion factor needed. Both ICESTATE and EVP operate with SI units
      
      integer k,l


      do k=1,nthick
         call xcaput(hsnow  (ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%hsnw,0)
         call xcaput(compact(ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%fice,0)
         call xcaput(hice   (ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%hice,0)
         call xcaput(rhosnw (ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%rhosnw,0)
         call xcaput(tsrf   (ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%tsrf,0)
         call xcaput(albs   (ilo:ihi,jlo:jhi,k),icestate(:,:)%ice(k)%albs,0)
         do l=1,nlaymax
            call xcaput(vtp   (ilo:ihi,jlo:jhi,l,k),icestate(:,:)%ice(k)%vtp(l),0)
         enddo
      enddo

      call xcaput(EVP_divu,ISTATE_divu,0)

   end subroutine evp2istate
end module m_evp2istate
