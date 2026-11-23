      module m_NOR05_resusp

      contains

      subroutine NOR05_resusp(icode,dti,n,k_list,tbs)
C***
C***BEGIN PROLOGUE RESUSP
C***DATE WRITTEN 13/9-95
C***AUTHOR
C            Henrik Soiland, Institute of Marine Research,
C            P.O. Box 1870,
C            N-5024 Bergen-Nordnes, Norway
C            Email:  henrik@imr.no
C***PURPOSE RESUSP performs resuspension  of particles from the bottom.
C
C***DESCRIPTION
C   
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
C         the field to be sunk on entry
      use mod_xc
      use mod_necessary_ecovars
C
      implicit none
      include 'sedcom.h'
      include 'common_blocks.h'
C
      integer icode
      real dti
      integer,intent(in) :: n
      integer,intent(in),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)::
     & k_list
      real,intent(in),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: 
     & tbs 
C
C Local variables
C
      integer :: i,j,k
      integer :: kmax
      real :: bflux,bfluxt,wdia,wfla,wcha,sink,rho0,utresh
      real :: ltresh,fluxfact,tau1,c2,nidet,sisis,sesed
      real :: p_h,p_l,p_b,s_l,d_sl
      real :: factor
      logical :: first 
      save first
      data first/.true./
C
C***FIRST executable statement SEDIMENT
      rho0=1000.0
C
      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point 
      !margin=0
      d_sl=5. !thickness of sediment-contribution
C     UTRESH is the minimum bottom stress value for resuspension to occur
C
c  Values below are from Pohlmann and Puls 1995 p 378.
c  Below experimental values from Eberhardt Mayer (See sediment routine also)
      utresh=0.78/rho0
      tau1=0.78/rho0
C
C     C2 is the slope of the linear increase in bottom flux.
C        [s/m]
C
      c2=(1.0E-4)*rho0
C
C      The formulas give the bottom flux in kg/m2, since mg/m3 is used
C      the flux has to be multiplied by a factor 1000000.
c
      c2=(1.0E+6)*c2
C
C     FLUXFACT is a factor that CAN  beused to tune to give resonable
C     bottom flux
C
      fluxfact=1.0
C
C NIDET is the inverse of the fraction of nitrogen in detritus
C SISIS is the inverse of the fraction of silisium in SIS (biogenic opal)
      nidet=detw/niw
      sisis=sisw/siw
C
C     
C     Resuspension only if the total bottom stress is
C     above a treshold value
C
      if(icode.eq.idet)then
C
C Resuspension of DET
C
        do j=1-margin,jj+margin
        do i=1-margin,ii+margin
        kmax=k_list(i,j)
        if(bot_layer(i,j,n,idet) .gt. 0.0) then
           if(tbs(i,j) .gt. utresh) then
              bfluxt=dti*c2*(tbs(i,j)-tau1)
     &              *bot_layer(i,j,n,idet)/totsed(i,j)*fluxfact
C  The actual formula looks like this, but NIDET drops out
C      &              *nidet/totsed(i,j)*bio(i,j,kdm+1,n,idet)/nidet*fluxfact
              bflux=min(bfluxt,bot_layer(i,j,n,idet))
c --- resuspension divided on the lowest 1 m.              
              do k=1,kdm
              if(dp(i,j,k,n)/onem.gt.0.01)then
                 p_h=p(i,j,k)
                 p_l=p(i,j,k+1)
                 p_b=p(i,j,kk+1)
                 s_l=p_b-d_sl*onem
                 if(p_h-s_l.gt.0)then
                    factor=1.
                 elseif(p_h-s_l.lt.0..and.p_l-s_l.gt.0.)then
                    factor=(p_l-s_l)/(p_l-p_h)
                 else
                    factor=0.
                 endif
                 sink=(bflux/d_sl)*factor*dp(i,j,k,n)/(onem)
                 sink=sink/(dp(i,j,k,n)/onem)
              else
                 sink=0.0
              endif
                 bio(i,j,k,n,idet) = bio(i,j,k,n,idet)+sink
              enddo
              bot_layer(i,j,n,idet) = bot_layer(i,j,n,idet)-bflux
           endif
        endif
        enddo
        enddo
C
      elseif(icode.eq.ifla)then
C
C Resuspension of FLA
C
         write(*,*)' NO RESUSPENSION OF FLAGELLATES IN THIS VERSION'
         stop
C
      elseif(icode.eq.idia)then
C
C Resuspension of DIA
C
         write(*,*)' NO RESUSPENSION OF DIATOMS IN THIS VERSION'
         stop
C
      elseif(icode.eq.isis)then
C
C Resuspension of SIS
C
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
        kmax=k_list(i,j)
        if(bot_layer(i,j,n,isis) .gt. 0.0) then
           if(tbs(i,j) .gt. utresh) then
              bfluxt=dti*c2*(tbs(i,j)-tau1)
     &              *bot_layer(i,j,n,isis)/totsed(i,j)*fluxfact
C  The actual formula looks like this, but SISIS drops out
C     &              *sisis/totsed(i,j)*bio(i,j,kdm+1)/sisis*fluxfact
              bflux=min(bfluxt,bot_layer(i,j,n,isis))
              do k=1,kdm
              if(dp(i,j,k,n)/onem.gt.0.01)then
                 p_h=p(i,j,k)
                 p_l=p(i,j,k+1)
                 p_b=p(i,j,kk+1)
                 s_l=p_b-d_sl*onem
                 if(p_h-s_l.gt.0.)then
                    factor=1.
                 elseif(p_h-s_l.lt.0 .and.p_l-s_l.gt.0)then
                    factor=(p_l-s_l)/(p_l-p_h) 
                 else
                    factor=0.
                 endif
                 sink=(bflux/d_sl)*factor*dp(i,j,k,n)/(onem)
                 sink=sink/(dp(i,j,k,n)/onem)
              else
                 sink=0.0
              endif
                 bio(i,j,k,n,isis)=bio(i,j,k,n,isis)+sink
              enddo
              bot_layer(i,j,n,isis) = bot_layer(i,j,n,isis)-bflux
           endif
        endif
      enddo
      enddo
C
C
#if defined (YELLOWSED)
      elseif(icode.eq.ised)then
C
C Resuspension of SED
C
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
        kmax=k_list(i,j)
        if(bot_layer(i,j,n,ised) .gt. 0.0) then
           if(tbs(i,j) .gt. utresh) then
C    Note that there is multiplication and division by FACTW, i.e.
C    it can be omitted, as is done below.(totsed is in mg sed in C2)
              bfluxt=dti*c2*(tbs(i,j)-tau1)
     &      *bot_layer(i,j,n,ised)/totsed(i,j)*fluxfact
              bflux=min(bfluxt,bot_layer(i,j,n,ised))
              do k=1,kdm
              if(dp(i,j,k,n)/onem.gt.0.01)then
                 p_h=p(i,j,k)
                 p_l=p(i,j,k+1)
                 p_b=p(i,j,kk+1)
                 s_l=p_b-d_sl*onem
                 if(p_h-s_l.gt.0.)then
                    factor=1.
                 elseif(p_h-s_l.lt.0 .and.p_l-s_l.gt.0)then
                    factor=(p_l-s_l)/(p_l-p_h) 
                 else
                    factor=0.
                 endif
                 sink=(bflux/d_sl)*factor*dp(i,j,k,n)/(onem)
                 sink=sink/(dp(i,j,k,n)/onem)
              else
                 sink=0.0
              endif
                 bio(i,j,k,n,ised) = bio(i,j,k,n,ised)+sink
              enddo
              bot_layer(i,j,n,ised) = bot_layer(i,j,n,ised)-bflux
           endif
        endif
        enddo
        enddo
C
      elseif(icode.eq.iyel)then
C
C Resuspension of YEL
C
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
        kmax=k_list(i,j)
        if(bot_layer(i,j,n,iyel) .gt. 0.0) then
           if(tbs(i,j) .gt. utresh) then
C    Note that there is multiplication and division by FACTW, i.e.
C    it can be omitted, as is done below.(totsed is in mg and C2 in mg
              bfluxt=dti*c2*(tbs(i,j)-tau1)
     &        *bot_layer(i,j,n,iyel)/totsed(i,j)*fluxfact
              bflux=min(bfluxt,bot_layer(i,j,n,iyel))
              do k=1,kdm
              if(dp(i,j,k,n)/onem.gt.0.01)then
                 p_h=p(i,j,k)
                 p_l=p(i,j,k+1)
                 p_b=p(i,j,kk+1)
                 s_l=p_b-d_sl*onem
                 if(p_h-s_l.gt.0.)then
                    factor=1.
                 elseif(p_h-s_l.lt.0 .and.p_l-s_l.gt.0)then
                    factor=(p_l-s_l)/(p_l-p_h) 
                 else
                    factor=0.
                 endif
                 sink=(bflux/d_sl)*factor*dp(i,j,k,n)/(onem)
                 sink=sink/(dp(i,j,k,n)/onem)
              else
                 sink=0.0
              endif
                 bio(i,j,k,n,iyel)=bio(i,j,k,n,iyel)+sink
              enddo
              bot_layer(i,j,n,iyel) = bot_layer(i,j,n,iyel)-bflux
           endif
        endif
        enddo
        enddo
#endif /*YELLOWSED*/
C
c      #ifdef CHATONELLA
c      elseif(icode.eq.icha)then
C
C Resuspension of CHA
C
c         write(*,*)' NO RESUSPENSION OF CHATONELLA IN THIS VERSION'
c         stop
c      #endif
C

#ifdef DETPHO
      elseif(icode.eq.idetp)then
C
C Resuspension of DETPHO
C
        do j=1-margin,jj+margin
        do i=1-margin,ii+margin
        kmax=k_list(i,j)
        if(bot_layer(i,j,n,idetp) .gt. 0.0) then
           if(tbs(i,j) .gt. utresh) then
              bfluxt=dti*c2*(tbs(i,j)-tau1)
     &              *bot_layer(i,j,n,idetp)/totsed(i,j)*fluxfact
              bflux=min(bfluxt,bot_layer(i,j,n,idetp))
c --- resuspension divided on the lowest 1 m.              
              do k=1,kdm
              if(dp(i,j,k,n)/onem.gt.0.01)then
                 p_h=p(i,j,k)
                 p_l=p(i,j,k+1)
                 p_b=p(i,j,kk+1)
                 s_l=p_b-d_sl*onem
                 if(p_h-s_l.gt.0)then
                    factor=1.
                 elseif(p_h-s_l.lt.0..and.p_l-s_l.gt.0.)then
                    factor=(p_l-s_l)/(p_l-p_h)
                 else
                    factor=0.
                 endif
                 sink=(bflux/d_sl)*factor*dp(i,j,k,n)/(onem)
                 sink=sink/(dp(i,j,k,n)/onem)
              else
                 sink=0.0
              endif
                 bio(i,j,k,n,idetp) = bio(i,j,k,n,idetp)+sink
              enddo
              bot_layer(i,j,n,idetp) = bot_layer(i,j,n,idetp)-bflux
           endif
        endif
        enddo
        enddo
#endif
C
      else
C
         write(*,*)'*** FATAL ERROR in RESUSP ***'
         write(*,*)'*** UNKNOWN VALUE of ICODE =',ICODE
         stop
      endif
C
      return
      end subroutine NOR05_resusp
      end module m_NOR05_resusp
