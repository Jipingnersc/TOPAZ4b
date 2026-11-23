module m_icestate_calcahcum


contains

   function icestate_calcahcum(icem)
      use mod_icestate_redist
      use mod_icestate_type
      use mod_icestate_classes
      use mod_icestate_thermo  , ONLY : fice_max
      implicit none

      type(t_ice),intent(in) :: icem(nthick)
      real :: icestate_calcahcum(0:nthick)

      real :: sigf,ficetot,gh 
      integer :: cl_lim,hk


      ficetot = SUM(icem%fice)

      ! note that after advection ficetot might be > fice_max  here
      ! sigf denotes this fraction (which must be redistributed in
      ! some way)
      sigf = max(0.,fice_max - ficetot)


      ! Calculate classes included in redistribution  (cl_lim)
      DO hk = 1,nthick
         sigf = sigf + icem(hk)%fice
         IF (sigf > gstar) EXIT 
      END DO
      cl_lim = hk - 1  


      ! Calculation of a(h), probability that ice of thickness h should
      ! participate to rafting or ridging. Notice that the following was
      ! adapted from Thorndike, 1975 (JGR), to a discrete case. Density 
      ! probabilities were actually integrated on ice classes.
      ! Function int_b is the integral of function b(h) as defined by 
      ! Thorndike (relation (7)). (See further down in code)
      icestate_calcahcum = 0.
      gh             = max(0.,fice_max - ficetot)


      ! No redistribution of open water area ....
      icestate_calcahcum(0)          = b(min(gh,gstar))


      ! icestate_calcahcum(hk) for every class. All of these classes are entirely involved in ridging.
      DO hk = 1, cl_lim     
         icestate_calcahcum(hk)    = b(gh) * icem(hk)%fice
         gh                     = gh + icem(hk)%fice
      END DO


      ! icestate_calcahcum(hk) for last class (which could be partially involved in ridging)
      icestate_calcahcum(cl_lim+1) = b(gstar)*min(icem(cl_lim+1)%fice,sigf-gstar)

   end function icestate_calcahcum




   ! =======================================================================
   ! ==============================   b   ==================================
   ! =======================================================================
   !
   ! Function b is the function b(x) as defined by 
   ! Thorndike (relation (7)). b(x) is a weighting distribution which
   ! is used to determine how much of g(h) is involved in redistrinution
   ! a(h) = b(h) * g(h)
   !
   real function  b(x)
      use mod_icestate_redist
      implicit none
      real, intent(in) :: x
      b= (1. - x/ gstar)*2./gstar
   end function b


end module m_icestate_calcahcum
