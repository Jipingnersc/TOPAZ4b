      module m_NOR05_phyto_sed

      contains 
 
      subroutine NOR05_phyto_sed(icode,dti,k_list,n,tbs)
C***
C***BEGIN PROLOGUE PHYTO_SED
C***DATE WRITTEN 13/9-95
C***AUTHOR
C            Henrik Soiland, Institute of Marine Research,
C            P.O. Box 1870,
C            N-5024 Bergen-Nordnes, Norway
C            Email:  henrik@imr.no
C***PURPOSE PHYTO_SED performs sedimentation of dead and living algae
C     to the bottom.
C
C***DESCRIPTION
C     First the sinking velocity is
C     found and then the fields are updated according to what
C     fraction of the field will sink to the next cell during
C     timestep DTI.
C
C***ROUTINES CALLED : NONE
C
C***END PROLOGUE
C
C   Global variables:
C
C   On entry :
C
C      ICODE :  = 4 (DIA)
C               = 5 (FLA)
C               = 6 (DET)
C               = 7 (SIS)
C               = 9 (SED)
C               =10 (YEL)
C               =11 (CHA)
C      TMP   : 3D-real, temporary storage, should have a copy of
C         the field to be sinked on entry
      use mod_xc
      use mod_necessary_ecovars
C
      implicit none
      include 'sedcom.h'
      include 'biocom.h'
      include 'common_blocks.h'
C
      integer icode
      integer,intent(in)::n
      integer,intent(in),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)::
     &  k_list
      real,intent(in):: dti
!      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm) ::tmp
      real,intent(in),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tbs
C
C Local variables
C
      integer :: i,j,k
      integer :: kmax
      real :: wdia,wsis,wdet,wfla,wcha,wsed,wyel,srdia
#if defined DETPHO
      real :: wdetp
#endif
      real :: sink,rho0,ltresh
      real :: layer_thickn
C
      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point 
      !margin=0
C
C***FIRST executable statement SEDIMENT
      rho0=1000.
C
C     LTRESH is a maximum bottom stress value for sedimentation to occur
c  Values below are from Pohlmann and Puls 1995 p 373.
c  Below experimental values from Eberhardt Mayer (See resusp routine also)
       ltresh=0.064 ! this value used in experiment 980327b
c       ltresh=0.1   ! Boon et al 1997
C
C     Since the stress is in m2/s2 it has to be dividede by rho
C
      ltresh=ltresh/rho0
C
C     In the main prog. a thickness for the bottom layer, DZ(KB) is defined.
C     A concentration is assosiated with this and thus a measure
C     for the amount of accumulated sediment at the bottom.
C
C     
C     Sedimentation only if the total bottom stress is
C     below a treshold value
C
      if(icode.eq.idet)then
c --- CH: Fra mortens filer: kb-1=kmax 
C
C Sinking of DET
C
        do j=1-margin,jj+margin
        do i=1-margin,ii+margin
        if(ip(i,j)==1)then
           kmax=k_list(i,j)
           layer_thickn=dp(i,j,kmax,n)/onem
           if(tbs(i,j) .lt. ltresh) then
              wdet = srdet*dayinv
              sink = wdet*bio(i,j,kmax,n,idet)
              bio(i,j,kmax,n,idet) =(bio(i,j,kmax,n,idet)*layer_thickn
     &                           -sink*dti)
     &                           /(layer_thickn+wdet*dti) 
              bot_layer(i,j,n,idet) =bot_layer(i,j,n,idet)
     &                           +sink*dti
           endif
        endif
        enddo
        enddo
C
        elseif(icode.eq.ifla)then
C
C Sinking of FLA
C
        do j=1-margin,jj+margin
        do i=1-margin,ii+margin
        if(ip(i,j)==1)then   
           kmax=k_list(i,j)
           layer_thickn=dp(i,j,kmax,n)/onem
           if(tbs(i,j) .lt. ltresh) then
              wfla = srfla*dayinv
              sink = wfla*bio(i,j,kmax,n,ifla)
              bio(i,j,kmax ,n,ifla) =(bio(i,j,kmax ,n,ifla)*layer_thickn
     &                           -sink*dti)
     &                           /(layer_thickn+wfla*dti)
              bot_layer(i,j,n,idet) = bot_layer(i,j,n,idet)
     &                           +sink*dti
           endif
        endif
        enddo
        enddo
C
        elseif(icode.eq.idia)then
C
C Sinking of DIA
C
        do j=1-margin,jj+margin
        do i=1-margin,ii+margin
        if(ip(i,j)==1)then
           kmax=k_list(i,j)
           layer_thickn=dp(i,j,kmax,n)/onem
           if(tbs(i,j) .lt. ltresh) then
              if(bio(i,j,kmax,n,isil)/csil.lt.sib)then
                 srdia = srdiamax
              else
                 srdia = srdiamin+(srdiamax-srdiamin)/(bio(i,j,kmax,n,3)
     &            /csil)
              srdia = (srdiamin+a3l/(bio(i,j,kmax,n,isil)/csil))
              endif
!CH: operating with the smallest dia-sinking at the moment              
!AS030709 - try the old one            srdia=srdiamin
              wdia = srdia*dayinv
              sink = wdia*bio(i,j,kmax,n,idia)
              bio(i,j,kmax ,n,idia) =(bio(i,j,kmax ,n,idia)*layer_thickn
     &                           -sink*dti)
     &                           /(layer_thickn+wdia*dti)         
              bot_layer(i,j,n,idet) = bot_layer(i,j,n,idet)
     &                            +sink*dti
!CH:--- The diatoms contributes also to the concentration of SIS at the
!CH:--- bottom              
              bot_layer(i,j,n,isis) = bot_layer(i,j,n,isis)
     &                            +cc(2)*sink*dti         
           endif
        endif
        enddo
        enddo
C
        elseif(icode.eq.isis)then
C
C Sinking of SIS
C
        do j=1-margin,jj+margin
        do i=1-margin,ii+margin
        if(ip(i,j)==1)then
           kmax=k_list(i,j)
           layer_thickn=dp(i,j,kmax,n)/onem
           if(tbs(i,j) .lt. ltresh) then
              wsis = srsis*dayinv
              sink = wsis*bio(i,j,kmax,n,isis)
              bio(i,j,kmax ,n,isis) =(bio(i,j,kmax ,n,isis)*layer_thickn
     &                           -sink*dti)
     &                           /(layer_thickn+wsis*dti)       
              bot_layer(i,j,n,isis) = bot_layer(i,j,n,isis)
     &                            +sink*dti
           endif
        endif
        enddo
        enddo

#if defined (YELLOWSED)
C
        elseif(icode.eq.ised)then
C
C Sinking of SED
C
        do j=1-margin,jj+margin
        do i=1-margin,ii+margin
        if(ip(i,j)==1)then
           kmax=k_list(i,j)
           layer_thickn=dp(i,j,kmax,n)/onem
           if(tbs(i,j).lt.ltresh.and.
     &       bio(i,j,kmax,n,ised).gt.sedmin)then
              wsed = srsed*dayinv
              sink = wsed*bio(i,j,kmax,n,ised)
              bio(i,j,kmax ,n,ised) =(bio(i,j,kmax ,n,ised)*layer_thickn
     &                           -sink*dti)
     &                           /(layer_thickn+wsed*dti)         
              bot_layer(i,j,n,ised) = bot_layer(i,j,n,ised)
     &                            +sink*dti
           endif
        endif
        enddo
        enddo
C
        elseif(icode.eq.iyel)then
C
C Sinking of YEL
C
        do j=1-margin,jj+margin
        do i=1-margin,ii+margin
        if(ip(i,j)==1)then 
           kmax=k_list(i,j)
           layer_thickn=dp(i,j,kmax,n)/onem
           if(tbs(i,j) .lt. ltresh .and.
     &                 bio(i,j,kmax,n,iyel) .gt. yelmin) then
              wyel = sryel*dayinv
              sink = wyel*bio(i,j,kmax,n,iyel)
              bio(i,j,kmax ,n,iyel) =(bio(i,j,kmax,n,iyel)*layer_thickn
     &                            -sink*dti)
     &                            /(layer_thickn+wyel*dti)          
              bot_layer(i,j,n,iyel) = bot_layer(i,j,n,iyel)
     &                             +sink*dti
           endif
        endif
        enddo
        enddo
#endif /*YELLOWSED*/

C
!      #ifdef CHATONELLA
!      elseif(icode.eq.icha)then
C
C Sinking of CHA
C
!      do 638 j=1-margin,jj+margin
!      do 638 i=1-margin,ii+margin
!        if(ip(i,j)==1) then
!           kmax=k_list(i,j)
!            layer_thickn=dp(i,j,kmax,n)/onem
!           if(tbs(i,j) .lt. ltresh) then
!              wcha = min(srcha*dayinv*dti/(dp(i,j,kmax,n)/onem),.8)
!              sink = wcha*tmp(i,j,kmax)
!              cha(i,j,kmax ) = cha(i,j,kmax )-sink*dti
!              cha(i,j,kdm+1) = cha(i,j,kdm+1)
!     &                        +sink*dti*(dp(i,j,kmax,n)/onem)
!           endif
!        endif
! 638    continue
!      #endif
C
#if defined DETPHO
      elseif(icode.eq.idetp)then
C
C Sinking of DETP
C
        do j=1-margin,jj+margin
        do i=1-margin,ii+margin
        if(ip(i,j)==1)then
           kmax=k_list(i,j)
           layer_thickn=dp(i,j,kmax,n)/onem
           if(tbs(i,j) .lt. ltresh) then
              wdetp = srdetp*dayinv
              sink = wdetp*bio(i,j,kmax,n,idetp)
              bio(i,j,kmax,n,idetp) =(bio(i,j,kmax,n,idetp)*layer_thickn
     &                           -sink*dti)
     &                           /(layer_thickn+wdetp*dti) 
              bot_layer(i,j,n,idetp) =bot_layer(i,j,n,idetp)
     &                           +sink*dti
           endif
        endif
        enddo
        enddo
#endif
C
      else
C
         write(*,*)'*** FATAL ERROR in PHYTOSED ***'
         write(*,*)'*** UNKNOWN VALUE of ICODE =',ICODE
         stop
      endif
C
C***END subroutine SEDIMENT
C
      return
      end subroutine NOR05_phyto_sed
      end module m_NOR05_phyto_sed
