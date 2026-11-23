      module m_icemodels_read_restart
      contains

      subroutine icemodels_read_restart(rt,nstep,dtime)
      use mod_xc
      use mod_year_info
      use mod_common_ice
      use mod_hycom_nersc
#if defined(ICE) || defined(ICESTATE)
      use mod_restart
#endif
#if defined (ICESTATE)
      !use mod_icestate_init
      !use mod_icestate_io
#endif
#if defined (EVP)
      use m_evp_restart
#endif
      implicit none
      type (year_info) , intent(in) :: rt
c
      real, dimension(itdm,jtdm) :: gdp1, gsaln1, gtemp1, 
     &   gficem, ghicem, ghsnwm, gtsrfm
      integer, intent(in) :: nstep
      real*8, intent(in) :: dtime

      logical :: ltmp
      integer :: reslt

c --- First, initialize from model fields
      call iceinit

#if defined (ICE) || defined (ICESTATE)
      if (.not.read_restart_ice(rt)) then
         if (mnproc==1) then 
            write(lp,*)'Error reading ice restart file ens number ',imem
            call flush(lp)
         end if
      end if
#endif

#ifdef ICESTATE
c      ! First attempt (always succeeds if ICE is read)
c      call icestate_MICOM_init(hicem,ficem,hsnwm,tsrfm,
c     &   dp(:,:,1,1)/onem,saln(:,:,1,1),temp(:,:,1,1),rhosnw)
c
c      ! Second attempt -- Only succeeds if ICESTATE restart file is present
c      ltmp = read_restart_icestate(rt)
#error ICESTATE not included yet
#endif
c
#if defined (EVP)
      reslt=evp_read_restart(rt,nstep,dtime)
      if (reslt/=0 .and. mnproc==1) then
         print *,'Can not read EVP restart file. I continue however...'
         write(lp,*) 'Result of evp_read_restart is ',reslt
      end if
#endif
      end subroutine icemodels_read_restart
      end module m_icemodels_read_restart


