c#if defined MP_TEST_DUMP
c#define MP_TEST_DUMP_CALL call mptestdump
c#else
c#define MP_TEST_DUMP_CALL !!! call mptestdump(
c#endif
      module m_NOR05_biochm

      contains 
      subroutine NOR05_biochm(m,n,rt)!,mp_test_step)
!
!***BEGIN PROLOGUE BIOLOGY
!***DATE WRITTEN 
!***AUTHOR Morten D. Skogen, morten@imr.no
!***PURPOSE BIOLOGY is the main routine for the chemical/biological and
!      the sedimentation/resuspension part of NORWECOM
!***REVISION DATE 9/1-97, MDS
!      Variable RAD is made 3D, and added concurrent comment to SGI
!***REVISION DATE 16/6-98, MDS
!      Extended to version 1.53, i.e. *all* biology is moved from MAIN
!      to BIOLOGY to improve portability. 
!***REVISION DATE 2/11-01, MDS
!      Incoming radiation data GLOBAL is made 2D (earlier one point).
!      Due to this several other arrays are extended to 2D :
!      LOW, MAXGLO, MEANGLO, BEAMTR, DIFFTR. This has resulted in changes in
!      the following routines in the light formulation : COMPRAD, IRRSUR
!      and TRANS. In addition LOWLIGHT has been included in COMPRAD and is no
!      longer necessary. Calling sequence to COMPRAD is changed.
!***REVISION DATE 26/2-02, MDS
!      Corrected bug, SUMDF, DAYNR and THETA new SAVE variables
!***REVISION DATE 26/5-04, MDS
!      Added CHATONELLA module
!
!***ROUTINES CALLED : MANY
!***END PROLOGUE BIOLOGY
!
      use mod_xc
      use mod_year_info, only : year_info
      use mod_necessary_ecovars
      use mod_common_ice, only : radfl0,ficem
      use m_NOR05_botstress
      use m_NOR05_affin
      use m_NOR05_atmnit
      use m_NOR05_biosun
      use m_NOR05_botsource
      use mod_checknan
      use m_NOR05_comprad  
      use m_NOR05_declin   
      use m_NOR05_detritus 
      use m_NOR05_irrsur   
      use m_NOR05_nconservation
      use m_NOR05_oxysat   
      use m_NOR05_phyto_growth
      use m_NOR05_zoo_growth
      use m_NOR05_phyto_sed   
!      use m_NOR05_phytosink   
      use m_NOR05_micomsink   
      use m_NOR05_radiation   
      use m_NOR05_resusp      
      use m_NOR05_slim        
      use m_NOR05_snell       
      use m_NOR05_source      
      use m_NOR05_trans       
      use m_NOR05_wqual    
!      use m_mptestdump

      implicit none

      include 'sedcom.h'
      include 'biocom.h'
      include 'common_blocks.h'

      integer,intent(in) :: m,n
      !integer,intent(in) :: mp_test_step !KAL for testing
      type(year_info), intent(in) :: rt
      
      real,dimension(kdm) :: dpmm
!      real,dimension(kdm) :: sink,w
      integer :: daynr
!
!
! Biology/light
!
      real dtib,umol2,phi
!
! Sediments
!
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: sissurf,oxysurf,
     &               yelsurf,sedsurf
      common/surfsed/sissurf,oxysurf,yelsurf,sedsurf
!
! Local variables
!
      integer :: i,j,k,icode,dag
      integer :: kbio,k_k
      real light,ttime_nor,theta_nor 
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)::
     &     kappa,sunti,sunhei,diffus, 
     &     direct,sumdf
     &     ,bsed,byel
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm)::source, 
     &     nmean 
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm,2) ::
     &     tmp 
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm):: bsource
CKAL  real,dimension(1:itdm,1:jtdm) :: modlon, modlat
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: tbs
!#ifdef CHATONELLA
!      real bcha(im,jm)
!#endif
! New 2D variables 
      logical,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: low
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: global,
     &               maxglo,meanglo,beamtr,difftr    
!
      logical first
      data first/.true./
      save first,nmean,low,global,maxglo,meanglo,beamtr, 
     &     difftr,sumdf,theta_nor
     
      integer :: kmax
      integer :: ndivis
      integer :: siad
      real    :: srdia
      real    :: max_nit,max_pho,max_sil,max_det,max_sis,max_fla,
     &           max_dia,max_oxy,max_sed,max_yel 
      real    :: min_nit,min_pho,min_sil,min_det,min_sis,min_fla,
     &           min_dia,min_oxy,min_sed,min_yel 
      real    :: max_bdet,max_bnit 
      real    :: min_bdet,min_bnit
      integer :: loc_m_nit(4)     
      integer :: loc_m_sil(4)     

!AS      sediment = .true.

      !KAL - Should use existing margin, 
      !margin=0 
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,kdm
         source(i,j,k)=0.
         bsource(i,j,k)=0.
      enddo
      enddo
      enddo
!--- Number of internal timesteps per hour      
      ndivis=3600./biodt_sou
      
CKAL  call xcaget(modlon,plon,0)
CKAL  call xcaget(modlat,plat,0)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,kdm
         dpmm(k)=max(onemm,dp(i,j,k,n))
      enddo
      do k=kdm,1,-1
      if(dpmm(k) .gt. tencm) then
         exit
      endif   
      enddo
      klist(i,j)=max(k,2)
      enddo
      enddo
      

! --- CH: calling botstress routine to get the total bottom stress for
! ---     sediment routines.
      if(sedflg==1) then
         call NOR05_botstress(n,tbs)
      endif

      
      !***FIRST EXECUTABLE STATEMENT BIOLOGY
!
!
      if(first)then
         call NOR05_affin
         first=.FALSE.
      end if
 

      phi=0.
      do 401 j=1-margin,jj+margin
      do 401 i=1-margin,ii+margin
          kappa(i,j)     = 0.
          diffus(i,j)    = 0.
          direct(i,j)    = 0.
          death_dia(i,j,:) = 0.
          death_fla(i,j,:) = 0.
!#ifdef CHATONELLA
!          DEATH_CHA(I,J) = 0.
!#endif
 401  continue
!
!   Read new climatologies in the frs-zones every 2'nd day.
!   Note that if BIVIS.NE.1, this must be syncronized with main program...
!
! --- CH:311005: removes this at the moment...
!      IF (MOD(nstep/ndivis,NDIVIS*48).EQ.0 .OR. nstep/ndivis.EQ.1) THEN
!         CALL clibio(rt%idd,H,ZZ)
!         IF(SEDIMENT)THEN
!            CALL clised(rt%idd,ZZ)
!         END IF
!      ENDIF
!
!***START light computation
!
! Find the declination of the sun and the transmitivity once every day
! Note that if BIVIS.NE.1, this must be syncronized with main program...
! --- CH:081105: Something wrong with ndivis/nstep:using bioday instead
          siad=int(60.0*60.0*24.0/biodt_sou)
          daynr=rt%idd
!
#if defined (DIURNAL)
!AS28092009 do this every 3 hours
          if (mod(nstep-1,siad/8)==0) then
#else
          if (mod(nstep-1,siad)==0) then
#endif
             global = radfl0*86400 
             call NOR05_declin(daynr,theta_nor)
CKAL         call NOR05_comprad(daynr,theta_nor,maxglo,meanglo,
CKAL &               sumdf,biodt_sou,modlon,modlat)
             call NOR05_comprad(daynr,theta_nor,maxglo,meanglo,
     &               sumdf,biodt_sou)
             call NOR05_trans(daynr,beamtr,difftr,global,low,
     &                        maxglo,meanglo)
!             call nconservation_2_file(daynr,n)
          endif
c        MP_TEST_DUMP_CALL ('ZNbc01' , n , mp_test_step) 
      
!
! Phytoplankton death rate
!
      call NOR05_detritus(klist,n)
!
      call NOR05_wqual(kappa)
      ttime_nor=real(rt%ihh)+mod(nstep/ndivis,ndivis)/(1.*ndivis)
!
CKAL  call NOR05_biosun(daynr,ttime_nor,theta_nor,modlon,modlat,
CKAL &                  sunti,sunhei)
CKAL - In MPI mode, plon/plat should be used (modlon, modlat are global arrays),
CKAL - plon/plate are available in common_blocks, so not passed as args
      call NOR05_biosun(daynr,ttime_nor,theta_nor,
     &                  sunti,sunhei)
      call NOR05_irrsur(beamtr,difftr,daynr, 
     &               diffus,direct,global,low,maxglo,meanglo, 
     &               sumdf,sunhei,biodt_sou)
!
c        MP_TEST_DUMP_CALL ('ZNbc02' , n , mp_test_step) 
      
!***END light computation
!
!***START primary-production
!
#if defined (ZOOPL)
      call NOR05_zoo_growth(klist,temp(:,:,:,n),n)
#endif
!
      do 691 j=1-margin,jj+margin
      do 691 i=1-margin,ii+margin
      kmax=klist(i,j)
      if(depths(i,j) .gt. 10.)then
         if(sunhei(i,j).lt.0)then
            do k=1,kdm
              rad(i,j,k)=0.
            enddo
         elseif(ficem(i,j) > 1e-6) then
            do k=1,kdm
              rad(i,j,k)=0.
            enddo
         elseif(low(i,j))then
            call NOR05_radiation(direct(i,j),diffus(i,j),i,j,
     &        kappa,0.,kmax,n)
         elseif(global(i,j).eq.0.)then
CKAL  ?? k is not set
            do k=1,kdm !AS
               rad(i,j,k)=0.
            enddo   !AS
         else
            phi=NOR05_snell(sunhei(i,j))
            call NOR05_radiation(direct(i,j),diffus(i,j),i,j,
     &       kappa,phi,kmax,n)
         endif
          call NOR05_phyto_growth(i,j,n)
       endif
 691  continue
!      !call nconservation('after phytogrowth',n)
c        MP_TEST_DUMP_CALL ('ZNbcpp' , n , mp_test_step) 

!
!***END primary-production
!
!  Find the total mass of the accumulated sediment in each gridpoint
!
       if (sedflg==1) then
          do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            kmax=klist(i,j)
!AS            totsed(i,j)=factw*bio(i,j,kmax,n,ised)+detw/niw*
!AS     &                      bio(i,j,kmax,n,idet) 
!AS     &                    + sisw/siw*bio(i,j,kmax,n,isis)+
!AS     &                      factw*bio(i,j,kmax,n,iyel)
            totsed(i,j)=detw/niw*bio(i,j,kmax,n,idet) 
     &                    + sisw/siw*bio(i,j,kmax,n,isis)
          enddo
          enddo
       endif
!      call nconservation('before source',n)
!
!-----------------------------------------------------------------------
!***START Adveksjon/diffusjon/synking av plankton og naering
!-----------------------------------------------------------------------
!
! First recompute W for fulfilment of the continuity equation
! from U and V used by the biology/sediments
!
! CH:271005      CALL VERTVL(2.*DTI,IM,JM,kmax,W,U,V,DT,ETF,ETB,scpx,scpy,A,C)
! CH:271005      CALL BOUND(5,DTE2,scpx,scpy,ELF,ELB,FLUXVA,D, 
!     &           UAF,VAF,FSM,DUM,DVM,UF,VF,S,T,Q2,Q2L,W,rt%idd)
! Update some variables used by superbeef advection
! CH:271005      CALL UPDATEDD(DT,DUM,DVM,DDZ,DUDZ,DVDZ)
!----------------------------------------
! Silicate  (SIL) and oxygen (OXY)
!----------------------------------------
      icode = isil
! Source terms
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
        bio(i,j,k,n,isil)=bio(i,j,k,n,isil)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)

         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,isil)=bio(i,j,k,n,isil)
     &                       +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
!      !call nconservation('after icode=3',n)
! Diffusion
! Update boundary-zone
!      CALL BIOBOUND(ICODE,rt%ihh)
!
! Oxygen
      
      icode = ioxy
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,ioxy)=bio(i,j,k,n,ioxy)+source(i,j,k)*biodt_sou
         bio(i,j,k,n,ioxy)=max(1.0e-4,bio(i,j,k,n,ioxy))
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,ioxy)=bio(i,j,k,n,ioxy)
     &                       +bsource(i,j,k)*biodt_sou
            bio(i,j,k,n,ioxy)=max(1.0e-4,bio(i,j,k,n,ioxy))
         enddo
         enddo
         enddo
      endif
!
! Oxygen saturation at the surface every 24hrs
!
!AS - april 2009: modified routine to work MPI
      if (mod(nstep-1,siad).eq.0 ) then
         call NOR05_oxysat(n) 

      endif
c        MP_TEST_DUMP_CALL ('ZNbcso' , n , mp_test_step) 
!      !call nconservation('after icode=ioxy',n)
!      CALL BIOBOUND(ICODE,rt%ihh)
!----------------------------------------
! Nitrogen og fosfor
!----------------------------------------
      icode = init
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,init)=bio(i,j,k,n,init)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,init)=bio(i,j,k,n,init)
     &                       +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
      call NOR05_atmnit(biodt_sou,n)
!      CALL BIOBOUND(ICODE,rt%ihh)
!      !call nconservation('after icode=1',n)
!
      icode = ipho
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,ipho)=bio(i,j,k,n,ipho)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,ipho)=bio(i,j,k,n,ipho)
     &                       +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
!      !call nconservation('after icode=ipho',n)
!      CALL BIOBOUND(ICODE,rt%ihh)
c        MP_TEST_DUMP_CALL ('ZNbcnf' , n , mp_test_step) 
!----------------------------------------
! Silikate shells (SIS) and Detritus (DET)
!----------------------------------------
      icode = isis
      call NOR05_source(icode,source,klist,n)
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,isis)=bio(i,j,k,n,isis)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,isis)=bio(i,j,k,n,isis)
     &                       +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
!
      if(sedflg==1)then
!
! Sedimentation of silicate shells
!
         call NOR05_phyto_sed(icode,biodt_sou,klist,n,tbs)
!
!  Resuspension of silicate shells
!
        call NOR05_resusp(icode,biodt_sou,n,klist,tbs)
      endif
!      CALL BIOBOUND(ICODE,rt%ihh)
!
      icode = idet
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,idet)=bio(i,j,k,n,idet)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,idet)=bio(i,j,k,n,idet)
     &                       +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
      if(sedflg==1)then
         call NOR05_phyto_sed(icode,biodt_sou,klist,n,tbs)
         call NOR05_resusp(icode,biodt_sou,n,klist,tbs)
      endif
!      !call nconservation('after icode=idet',n)
!      CALL BIOBOUND(ICODE,rt%ihh)
c        MP_TEST_DUMP_CALL ('ZNbcsd' , n , mp_test_step) 
!----------------------------------------
! Diatomeer og flagellater and chatonella
!----------------------------------------
      icode = idia
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,idia)=bio(i,j,k,n,idia)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,idia)=bio(i,j,k,n,idia)
     &                       +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
      if(sedflg==1)then
         call NOR05_phyto_sed(icode,biodt_sou,klist,n,tbs)
!        CALL NOR05_resusp(icode,biodt_sou,n,klist,tbs)
      endif
!      !call nconservation('after icode=idia',n)
!      CALL BIOBOUND(ICODE,rt%ihh)
!-START FLAgellates----------------------------------------------------
      icode = ifla
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,ifla)=bio(i,j,k,n,ifla)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,ifla)=bio(i,j,k,n,ifla)
     &                       +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
      if(sedflg==1)then
         call NOR05_phyto_sed(icode,biodt_sou,klist,n,tbs)
!        CALL NOR05_resusp(icode,biodt_sou,n,klist,tbs)
      endif
!      !call nconservation('after icode=ifla',n)
!      CALL BIOBOUND(ICODE,rt%ihh)
!
!#ifdef CHATONELLA
!!-START CHAtonella----------------------------------------------------
!      ICODE = icha
!      CALL NOR05_source(icode,source,klist,n)
!      IF (SEDIMENT) THEN
!         CALL NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
!         bio(:,:,:,n,?)=bio(:,:,:,n,?)+bsource(:,:)*biodt_sou
!      ENDIF
!      tmp = CHA
!      CALL NOR05_phytosink(icode,tmp,biodt_sou,n)
!      IF(SEDIMENT)THEN
!         tmp = CHA
!         CALL NOR05_phyto_sed(icode,tmp,biodt_sou,klist,n,tbs)
!!!!        CALL NOR05_resusp(icode,biodt_sou,n,klist,tbs)
!      ENDIF
!      CALL BIOBOUND(ICODE,rt%ihh)
!#endif
!----------------------------------------
! Yellow substance (YEL) and mineral sediment (SED)
!----------------------------------------
#if defined (YELLOWSED)
      icode = iyel
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,iyel)=bio(i,j,k,n,iyel)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,iyel)=bio(i,j,k,n,iyel)
     &                        +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
      if(sedflg==1)then
         call NOR05_phyto_sed(icode,biodt_sou,klist,n,tbs)
         call NOR05_resusp(icode,biodt_sou,n,klist,tbs)
      endif
!     ! call nconservation('after icode=iyel',n)
!      CALL BIOBOUND(ICODE,rt%ihh)
!
      icode = ised
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,ised)=bio(i,j,k,n,ised)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,ised)=bio(i,j,k,n,ised)
     &                       +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
      if(sedflg==1)then
         call NOR05_phyto_sed(icode,biodt_sou,klist,n,tbs)
         call NOR05_resusp(icode,biodt_sou,n,klist,tbs)
      endif
#endif /*YELLOWSED*/
!      call nconservation('after source',n)
!      CALL BIOBOUND(ICODE,rt%ihh)
#if defined DETPHO
      icode = idetp
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,idetp)=bio(i,j,k,n,idetp)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,idetp)=bio(i,j,k,n,idetp)
     &                       +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
      if(sedflg==1)then
         call NOR05_phyto_sed(icode,biodt_sou,klist,n,tbs)
         call NOR05_resusp(icode,biodt_sou,n,klist,tbs)
      endif
#endif
#if defined ZOOPL
      icode = imeso
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,imeso)=bio(i,j,k,n,imeso)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,imeso)=bio(i,j,k,n,imeso)
     &                       +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
!
      icode = imicro
      call NOR05_source(icode,source,klist,n)
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
      do k=1,klist(i,j)
         bio(i,j,k,n,imicro)=bio(i,j,k,n,imicro)+source(i,j,k)*biodt_sou
      enddo
      enddo
      enddo
      if (sedflg==1) then
         call NOR05_botsource(icode,bsource,biodt_sou,klist,n,tbs)
         do i=1-margin,ii+margin
         do j=1-margin,jj+margin
         do k=1,klist(i,j)
            bio(i,j,k,n,imicro)=bio(i,j,k,n,imicro)
     &                       +bsource(i,j,k)*biodt_sou
         enddo
         enddo
         enddo
      endif
#endif!
!***END Adveksjon/diffusjon/synking av plankton og naering
c        MP_TEST_DUMP_CALL ('ZNbcys' , n , mp_test_step) 
      
!---CH:sinking
      if(sinkfl==1) then
         call NOR05_micomsink(n,biodt_sou)
      endif
!*** DIA and FLA becomes DET and SIS as they hit the bottom.
!
      
      max_nit=-1000.
      max_pho=-1000.
      max_sil=-1000.
      max_det=-1000.
      max_sis=-1000.
      max_fla=-1000.
      max_dia=-1000.
      max_oxy=-1000.
      max_sed=-1000.
      max_yel=-1000.
      min_nit=10000.
      min_pho=10000.
      min_sil=10000.
      min_det=10000.
      min_sis=10000.
      min_fla=10000.
      min_dia=10000.
      min_oxy=10000.
      min_sed=10000.
      min_yel=10000.
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      if(ip(i,j)==1)then
         min_bdet=min(min_bdet,bot_layer(i,j,n,idet))
         min_bnit=min(min_bnit,bot_layer(i,j,n,init))
         max_bdet=max(max_bdet,bot_layer(i,j,n,idet))
         max_bnit=max(max_bnit,bot_layer(i,j,n,init))
      do k=1,klist(i,j)
      if(ip(i,j)==1.and.depths(i,j)>1.0) then
         max_nit=max(max_nit,bio(i,j,k,n,init))
         min_nit=min(min_nit,bio(i,j,k,n,init))
         max_pho=max(max_pho,bio(i,j,k,n,ipho))
         min_pho=min(min_pho,bio(i,j,k,n,ipho))
         max_sil=max(max_sil,bio(i,j,k,n,isil))
         min_sil=min(min_sil,bio(i,j,k,n,isil))
         max_det=max(max_det,bio(i,j,k,n,idet))
         min_det=min(min_det,bio(i,j,k,n,idet))
         max_sis=max(max_sis,bio(i,j,k,n,isis))
         min_sis=min(min_sis,bio(i,j,k,n,isis))
         max_fla=max(max_fla,bio(i,j,k,n,ifla))
         min_fla=min(min_fla,bio(i,j,k,n,ifla))
         max_dia=max(max_dia,bio(i,j,k,n,idia))
         min_dia=min(min_dia,bio(i,j,k,n,idia))
         max_oxy=max(max_oxy,bio(i,j,k,n,ioxy))
         min_oxy=min(min_oxy,bio(i,j,k,n,ioxy))
#if defined (YELLOWSED)
         max_sed=max(max_sed,bio(i,j,k,n,ised))
         min_sed=min(min_sed,bio(i,j,k,n,ised))
         max_yel=max(max_yel,bio(i,j,k,n,iyel))
         min_yel=min(min_yel,bio(i,j,k,n,iyel))
#endif /*YELLOWSED*/
      endif
      enddo
      endif
      enddo
      enddo

      return
      end subroutine NOR05_biochm
      end module m_NOR05_biochm  
