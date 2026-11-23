      module m_NOR05_botsource

      contains
  
      subroutine NOR05_botsource(icode,source,dti,k_list,n,tbs)
C***
C***BEGIN PROLOGUE BOTSOURCE
C***AUTHOR
C        HENRIK SOILAND. henrik@imr.no
C       Institute of Marine Research,
C       Postboks 1870, N-5024 Bergen-Nordnes, Norway
C      Based on BIOSOURCE routine by
C      K{\aa}re B. Ulvestad, Department of Fisheries and
C      Marine Biology, Bergen High Technology Centre,
C      5014 Bergen, Norway
C
C
C***PURPOSE BOTSOURCE computes the source terms for all chemical
C     and biological equations in the sediment layer (kmax) and an additional
C     source term in the lowest sigma layer due to leakage from the
C     bottom (NIT, PHO and SIL).
C     NB! the source term due to algea growth and death in layer kmax-1 is
C     already accounted for in the biosource routine,
C     i.e. the leakage from the bottom is added to this term for k=kmax-1
C
C***DESCRIPTION all source-terms are given relative to the
C     source-term for NITrate
C
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
C     ICODE    : Integer code defining the source function
C          =1  : Nitrate
C          =2  : Phosphate
C          =3  : Silicate
C          =4  : Diatoms
C          =5  : Flagellates
C          =6  : Detritus
C          =7  : Silicate particles
C          =8  : Oxygen
C          =11 : Chatonella
C          =12 : DetritusP
C     others   : see the files globals.doc and sedcom.h 
C                for description
C
C   On exit :
C
C     SOURCE : 3D-field, source term for the actual chemical or
C              biological constitutent
C              For k=1,kmax-2 no change,
C              for k=kmax-1 the leakage from the bottom added,
C              for k=kmax the source term is computed.
C              NB! source term in the sediment is different from in the water.
C
      use mod_xc
      use mod_necessary_ecovars
C
      implicit none

      include 'sedcom.h'
      include 'biocom.h'
      include 'common_blocks.h'
C
      integer icode
      real,intent(out),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm)
     &                            :: source
      real :: dti
      integer,intent(in) :: n
      integer,intent(in),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &                 k_list 
      real,intent(in),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tbs
C
C   Local variables
C
      integer ::  i,j,k
      real :: bflux,rho0,c1,c2
      real :: maxnit,minnit,maxpho,minpho,maxsil,minsil,maxtau
      real :: p_h,p_l,p_b,s_l,d_sl
      real :: factor
      data  minnit/0./,maxnit/0. /,minpho/0. /,maxpho/0. /
      data  minsil/0.1 /,maxsil/2.5 /,maxtau/2.0E-3 /
      logical first
      save first
      data first/.true./
C
      real,dimension(kdm) :: dpmm     
      integer :: kmax
C
C First executable statement BOTSOURCE
C
      rho0=1000.0
      d_sl=5.0 !thickness of layer where sediments are resuspended(m)
C
      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point 
      !margin=0
      source(:,:,:)=0.
C
c
      if(icode.eq.init)then
C
C     NIT
C
C         C2= (MAXNIT-MINNIT)/MAXTAU
         c1=(scc(6)-scc(5))/maxtau
C
C     BOTTOM
C
         do  j=1-margin,jj+margin
         do  i=1-margin,ii+margin
            kmax=k_list(i,j)
            bflux=0.0
            if(depths(i,j).gt.1.0) then
            if(bot_layer(i,j,n,init) .gt. 0.0) then
c leakage a function of bottom stress
               c2=min(scc(5)+c1*tbs(i,j),scc(6))
               bflux=c2*bot_layer(i,j,n,init)
               do k=1,kdm
               if(dp(i,j,k,n)/onem.gt.0.01) then
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
                  source(i,j,k)=(bflux/d_sl)*factor*dp(i,j,k,n)
     &                         /(onem)
                  source(i,j,k)=source(i,j,k)/(dp(i,j,k,n)/onem)
               else
                  source(i,j,k)=0.0  
               endif
               enddo
            endif
c Denitfrication is a constant fraction of detritus mineralized at the bottom
            bot_layer(i,j,n,init)=bot_layer(i,j,n,init)+
     &                        dti*((1.0-denitfr)*cc(4)*
     &                        min(bot_layer(i,j,n,idet),detmax)
     &                        -bflux)
c Corrected 1/3-99
c Denit is accumulated amount nitr. that has been denitr.
             denit(i,j)=denit(i,j)+dti*(denitfr)*
     &                  cc(4)*min(bot_layer(i,j,n,idet),detmax)
c   denitification is timedependent
c                 bio(i,j,kdm+1,n,init)=bio(i,j,kdm+1,n,init)+DTI*(CC(4)*bio(i,j,kdm+1,n,idet)-
c     &                         BFLUX-SCC(7)*bio(i,j,kdm+1,n,init))

             else
                source(i,j,:)=0.
             endif
          enddo
          enddo
C
          elseif(icode.eq.ipho)then
C
C     PHO
C
C         C2=(MAXPHO-MINPHO)/MAXTAU
         c1=(scc(6)-scc(5))/maxtau
C     BOTTOM
C
         do j=1-margin,jj+margin 
         do i=1-margin,ii+margin
            kmax=k_list(i,j)
            bflux=0.0
            if(depths(i,j).gt.1.0) then
            if(bot_layer(i,j,n,ipho) .gt. 0.0) then
c leakage a function of bottom stress
              c2=min(scc(5)+c1*tbs(i,j),scc(6))
              bflux=c2*(bot_layer(i,j,n,ipho))
C              BFLUX=min(MINPHO+C2*TTAU(I,J),
C     &                    MAXPHO)/(dp(i,j,kmax,n)/onem)
              do k=1,kdm
              if(dp(i,j,k,n)/onem.gt.0.01) then
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
c---CH:to avoid large values                 
                 source(i,j,k)=(bflux/d_sl)*factor*dp(i,j,k,n)
     &                        /(onem)          
                 source(i,j,k)=source(i,j,k)/(dp(i,j,k,n)/onem)
              else
                 source(i,j,k)=0.0
              endif
              enddo
            endif
            bot_layer(i,j,n,ipho)=bot_layer(i,j,n,ipho)+
#if defined (DETPHO)
     &               dti*(1.3*cc(4)*min(bot_layer(i,j,n,idetp),detpmax)
#else
     &               dti*(cc(1)*cc(4)*min(bot_layer(i,j,n,idet),detmax)
#endif
     &                      -bflux) 
            else
               source(i,j,:)=0.
            endif
         enddo
         enddo
C
         elseif(icode.eq.isil)then
C
C     SIL
C
C         C2= (MAXSIL-MINSIL)/MAXTAU
         c1=(scc(6)-scc(5))/maxtau
         source(:,:,:)=0.0
C
C     BOTTOM
C
         do j=1-margin,jj+margin 
         do i=1-margin,ii+margin 
            kmax=k_list(i,j)
            bflux=0.0
            if(depths(i,j).gt.1.0) then
            if(bot_layer(i,j,n,isil) .gt. 0.0) then
c leakage a function of bottom stress
               c2=min(scc(5)+c1*tbs(i,j),scc(6))
               bflux=c2*(bot_layer(i,j,n,isil))

C                  BFLUX=min(MINSIL+C2*TTAU(I,J),
C     &                       MAXSIL)/(dp(i,j,kmax,n)/onem)
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
c---CH:to avoid large values at bottom
                  source(i,j,k)=(bflux/d_sl)*factor*dp(i,j,kmax,n)
     &                         /(onem)
                  source(i,j,k)=source(i,j,k)/(dp(i,j,k,n)/onem)
               else
                  source(i,j,k)=0.0
               endif
               enddo
            endif
               bot_layer(i,j,n,isil)=bot_layer(i,j,n,isil)+dti*(scc(4)*
     &                           (min(bot_layer(i,j,n,isis),sismax))-
     &                           bflux)
            else
               source(i,j,:)=0.
            endif
         enddo
         enddo
C
         elseif(icode.eq.idia)then
         source(:,:,:)=0.
C
C     DIA (IN THE SIMPLE MODEL, ALL DIATOMS ARE CONVERTED TO 
C     DETRITUS AS THEY HIT THE BOTTOM, THIS IS DONE IN DEADBOT)
C
         elseif(icode.eq.ifla)then
         source(:,:,:)=0.
C
C     FLA (IN THE SIMPLE MODEL, ALL FLAGELLATES ARE CONVERTED TO 
C     DETRITUS AS THEY HIT THE BOTTOM, THIS IS DONE IN DEADBOT)
C
         elseif(icode.eq.idet)then
         source(:,:,:)=0.
C
C     BOTTOM
C
         do j=1-margin,jj+margin 
         do i=1-margin,ii+margin 
         if(depths(i,j).gt.1.0)then
            kmax=k_list(i,j)
c
c Amount burried in the sediemnt
c
            detbu(i,j)=detbu(i,j)+
     &                 dti*scc(8)*max(bot_layer(i,j,n,idet)-detbul,0.)
c
            bot_layer(i,j,n,idet)=bot_layer(i,j,n,idet)-
     &                     dti*(cc(4)*min(bot_layer(i,j,n,idet),detmax)+
     &                     scc(8)*max(bot_layer(i,j,n,idet)-detbul,0.))
c            
         endif
         enddo
         enddo
C
         elseif(icode.eq.isis)then
C
C     SIS
C
C
C     BOTTOM
C
         do j=1-margin,jj+margin 
         do i=1-margin,ii+margin 
         if(depths(i,j).gt.1.0)then
            kmax=k_list(i,j)
c Amount burried in the sediemnt
            sisbu(i,j)=sisbu(i,j)+
     &                 dti*scc(8)*max(bot_layer(i,j,n,isis)-sisbul,0.)

            bot_layer(i,j,n,isis)=bot_layer(i,j,n,isis)-
     &                     dti*(scc(4)*min(bot_layer(i,j,n,isis),
     &                     sismax)+
     &                     scc(8)*max(bot_layer(i,j,n,isis)-sisbul,0.))
         endif
         enddo
         enddo
C
         elseif(icode.eq.ioxy)then
C
C     OXY
C
C
C     BOTTOM (Oxidation at the bottom use oxygen in the water)
C     No oxygen budget in the sediment!
C
         do j=1-margin,jj+margin 
         do i=1-margin,ii+margin 
         if(depths(i,j).gt.1.0) then
            kmax=k_list(i,j)
c               SOURCE(I,J,kmax)=SOURCE(I,J,kmax)+
c     &          (SCC(3)*SCC(7)*bio(i,j,kdm+1,n,init)-CC(4)*SCC(1)*bio(i,j,kdm+1,idet))
c     &          /(dp(i,j,kmax,n))
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
               source(i,j,k)=
     &                     (((1.0-denitfr)*scc(1)+denitfr*
     &                     (scc(1)-scc(3)))*cc(4)*
     &                     min(bio(i,j,kmax,n,idet),detmax))
     &                     *(dp(i,j,kmax,n)*factor/(onem*d_sl))
               source(i,j,k)=source(i,j,k)/(dp(i,j,k,n)/onem)
               source(i,j,k)=-source(i,j,k)
            else
               source(i,j,k)=0.0   
            endif
            enddo
         else
            source(i,j,:)=0.
         endif
         enddo
         enddo
C
         elseif(icode.eq.ised)then
         source(:,:,:)=0.
         return
         elseif(icode.eq.iyel)then
         source(:,:,:)=0.
         return
#ifdef DETPHO
         elseif(icode.eq.idetp)then
C
C     DETP
C
C
C     BOTTOM
C
         do j=1-margin,jj+margin 
         do i=1-margin,ii+margin 
         if(depths(i,j).gt.1.0)then
            kmax=k_list(i,j)
c Amount burried in the sediemnt
            detpbu(i,j)=detpbu(i,j)+
     &                 dti*scc(8)*max(bot_layer(i,j,n,idetp)-detpbul,0.)

            bot_layer(i,j,n,idetp)=bot_layer(i,j,n,idetp)-
     &                   dti*(1.3*cc(4)*min(bot_layer(i,j,n,idetp),
     &                   detpmax)+
     &                   scc(8)*max(bot_layer(i,j,n,idetp)-detpbul,0.))
         endif
         enddo
         enddo
#endif      
c#ifdef CHATONELLA
c         elseif(icode.eq.icha)then
C
C     CHA (IN THE SIMPLE MODEL, ALL CHATONELLA ARE CONVERTED TO 
C     DETRITUS AS THEY HIT THE BOTTOM, THIS IS DONE IN DEADBOT)
C
c         return
c#endif
#if defined (ZOOPL)
      elseif (icode.eq.imeso) then
         source(:,:,:)=0.
         return
      elseif (icode.eq.imicro) then
         source(:,:,:)=0.
         return
#endif
C
      else
         write(*,*)'*** FATAL ERROR in BOTSOURCE ***'
         write(*,*)'*** UNKNOWN VALUE of ICODE =',ICODE
         stop
      endif
C
      return
      end subroutine NOR05_botsource
      end module m_NOR05_botsource
