      module m_icemodels_init
C --- This routine takes care of all ice model initialization 
C --- giving a CLEANER HYCOM code.
      contains
      subroutine icemodels_init()
      use mod_xc
      use mod_year_info, only: year_info
#ifdef ICESTATE
C     use mod_icestate_init
#error  ICESTATE not set up yet
#endif
      use mod_common_ice
      implicit none
      include 'common_blocks.h'
c
c --- Initialize variables in mod_common_ice
      call icedat
c
c --- Initialize Hunke&Dukowicz EVP dynamics model
#if defined (EVP)
      if (mnproc==1) print *,' Initializing evp from hycom'
      call evp_init_from_hycom()
      if(mnproc==1) print *,' Initializing evp'
      call evp_init()
#endif
c
c --- Initialize ICESTATE thermodynamic/kinematic model
#ifdef ICESTATE
C     if (mnproc==1) print '(a)','ICESTATE init:'
C     call icestate_read_infile
C     ISTstep = baclin
C     call icestate_init(ISTstep,depths,plat/radian,
C    &                   thref)
#error  ICESTATE not set up yet
#endif
      end subroutine icemodels_init
      end module m_icemodels_init
