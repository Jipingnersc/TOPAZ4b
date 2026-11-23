      module m_icemodels_save_restart
      contains
c
      subroutine icemodels_save_restart(rt,nstep,dtime,final)
      use mod_xc
      use mod_year_info
      use mod_common_ice
      use mod_hycom_nersc
#if defined(ICE) || defined(ICESTATE)
      use mod_restart
#endif
#ifdef ICESTATE
      use mod_icestate_io
#endif
#if defined (EVP)
      use m_evp_restart
#endif
      implicit none
      type(year_info),  intent(in) :: rt
      integer,          intent(in) :: nstep
      real*8,           intent(in) :: dtime
      logical,          intent(in) :: final
c
      logical :: ltmp
      integer :: itmp
c
#ifdef ICE
#if defined (NEW_ICE_RESTART)
      itmp=save_restart_ice_new(rt,nstep,dtime)
#else
      ltmp=save_restart_ice(rt)
#endif
#endif
c
#ifdef ICESTATE
      ltmp = save_restart_icestate(rt)
#endif
c
#if defined (EVP)
      itmp=evp_save_restart(rt,nstep,dtime)
#endif
      end subroutine icemodels_save_restart
      end module m_icemodels_save_restart
