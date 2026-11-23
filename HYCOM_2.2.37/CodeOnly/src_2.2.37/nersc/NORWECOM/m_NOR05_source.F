      module m_NOR05_source 

      contains

      subroutine NOR05_source(icode,source,k_list,n)
!***
!***BEGIN PROLOGUE BIOSOURCE
!***AUTHOR
!      K{\aa}re B. Ulvestad, Department of Fisheries and
!      Marine Biology, Bergen High Technology Centre,
!      5014 Bergen, Norway
!
!***Revised date 26/8-92 Morten D. Skogen, morten@imr.no
!     Added prologue
!***Revised date 27/10-05 Morten D. Skogen, morten@imr.no
!     Added DETP. This has 30% faster regeneration than DET,
!     besides only 75 of the death goes there, while the rest 25%
!     goes directly back to PHO
!
!***PURPOSE BIOSOURCE computes the source terms for all chemical
!     and biological equations
!
!***DESCRIPTION all source-terms is given relative to the
!     source-term for NITrate
!
!     NIT = RED_DIA+RED_FLA+RED_CHA + CC(4)*DET - P_DIA-P_FLA-P_CHA
!
!     Phosphate and silicate (NB! no silicate in the flagellates, and
!     no regeneration of the silicate either),
!     has a source term equals a constant times NIT's source-term
!     (the constant equals the Redfield ratio P/N, Si/N).
!     DET source-term is the change in the internal-pool of nitrate,
!     i.e. DEATH - REGENERATION, and the source term for the
!     phytoplankton is PROD - EXCRETION - DEATH
!
!***ROUTINES CALLED : NONE
!
!***END PROLOGUE
!
!   Global variables:
!
!   On entry :
!
!     ICODE    : Integer code defining the source function
!          =1  : Nitrate
!          =2  : Phosphate
!          =3  : Silicate
!          =4  : Diatoms
!          =5  : Flagellates
!          =6  : Detritus
!          =7  : Silicate particles
!          =8  : Oxygen
!          =9  : Sediment
!          =10 : Yellow substances
!          =11 : Chatonella
!          =12 : DetritusP
!     PROD_DIA : 3D-field, photosynthetic production of diatoms
!     PROD_FLA : 3D-field, photosynthetic production of flagellates
!     PROD_CHA : 3D-field, photosynthetic production of chatonella
!     RED_DIA  : 3D-field, excretion of diatoms
!     RED_FLA  : 3D-field, excretion of flagellates
!     RED_CHA  : 3D-field, excretion of chatonella
!     DEATH_DIA: 2D-field, death rate diatoms
!     DEATH_FLA: 2D-field, death rate flagellates
!     DEATH_CHA: 2D-field, death rate chatonella
!     others   : see the files globals.doc for description
!
!   On exit :
!
!     SOURCE : 3D-field, source term for the actual chemical or
!              biological constitutent
      use mod_xc
      use mod_necessary_ecovars
!
      implicit none
      include 'sedcom.h'
      include 'biocom.h'
!
      integer icode
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm)::
     &      source
      integer,intent(in):: n
      integer,intent(in),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)::
     &      k_list 
!
!   Local variables
!
      integer   i,j,k,o
!
! First executable statement BIOSOURCE
      
!
      source = 0.
      if(icode.eq.init)then
!
!     NIT
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j)
#ifdef NEWPRODUCTION
         source(i,j,k)=-1.*(prod_dia(i,j,k)+prod_fla(i,j,k))
#else
         source(i,j,k)=red_dia(i,j,k)+red_fla(i,j,k)+ 
!#ifdef CHATONELLA
!     &                    RED_CHA(I,J,K)-PROD_CHA(I,J,K)+ 
!#endif
     &                    cc(4)*bio(i,j,k,n,idet) 
     &                    -prod_dia(i,j,k)-prod_fla(i,j,k)
#ifdef ZOOPL
         source(i,j,k) = source(i,j,k) + zen_don(i,j,k) + zin_don(i,j,k)
#endif
#endif
!
!AS#ifdef NH4_N
!AS         SOURCE(I,J,K) = SOURCE(I,J,K) + &
!AS     &          0.10*( DEATH_DIA(I,J,K)*TBIO(I,J,K,IDIA)+ &
!AS#ifdef CHATONELLA
!AS     &          DEATH_CHA(I,J)*TBIO(I,J,K,ICHA)+ &
!AS#endif
!AS     &          DEATH_FLA(I,J,K)*TBIO(I,J,K,IFLA) )
!AS#endif
!AS!
!AS#endif
      enddo
      enddo
      enddo
!
      elseif(icode.eq.ipho)then
!
!     PHO
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j)
#ifdef NEWPRODUCTION
         source(i,j,k)=-cc(1)*(prod_dia(i,j,k)+prod_fla(i,j,k))
#else
         source(i,j,k)=cc(1)*(red_dia(i,j,k)+red_fla(i,j,k)+ 
!#ifdef CHATONELLA
!     &                    red_cha(i,j,k)-prod_cha(i,j,k)+ 
!#endif
#ifdef DETPHO
     &                    1.3*cc(4)/cc(1)*bio(i,j,k,n,idetp)- 
#else
     &                    cc(4)*bio(i,j,k,n,idet)- 
#endif /* DETPHO */
     &                    prod_dia(i,j,k)-prod_fla(i,j,k))
!
#ifdef DETPHO
          source(i,j,k) = source(i,j,k) + 
     &          cc(1)*0.25*(death_dia(i,j,k)*bio(i,j,k,n,idia)+ 
!#ifdef CHATONELLA
!     &          DEATH_CHA(I,J)*TBIO(I,J,K,ICHA)+ 
!#endif
     &          death_fla(i,j,k)*bio(i,j,k,n,ifla) )
#endif /* DETPHO */
!
#endif
#ifdef ZOOPL
         source(i,j,k)=source(i,j,k) +
     &                 cc(1)*(zen_don(i,j,k) + zin_don(i,j,k))
#endif
      enddo
      enddo
      enddo
!
      elseif(icode.eq.isil)then
!
!     SIL
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do  k=1,k_list(i,j)
         source(i,j,k)=-cc(2)*prod_dia(i,j,k)+scc(4)*bio(i,j,k,n,isis)
      enddo
      enddo
      enddo
!
      elseif(icode.eq.idia)then
!
!     DIA
!
      do  j=1-margin,jj+margin
      do  i=1-margin,ii+margin
      do  k=1,k_list(i,j)
         source(i,j,k)=prod_dia(i,j,k)-red_dia(i,j,k)- 
     &         death_dia(i,j,k)*bio(i,j,k,n,idia)
#ifdef ZOOPL
         source(i,j,k)=source(i,j,k) - p1n_zen(i,j,k) - p1n_zin(i,j,k)
#endif
      enddo
      enddo
      enddo
!
      elseif(icode.eq.ifla)then
!
!     FLA
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j)
         source(i,j,k)=prod_fla(i,j,k)-red_fla(i,j,k)- 
     &         death_fla(i,j,k)*bio(i,j,k,n,ifla)
#ifdef ZOOPL
         source(i,j,k)=source(i,j,k) - p2n_zin(i,j,k)
#endif!
      enddo
      enddo
      enddo

      elseif(icode.eq.idet)then
!
!     DET
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j)
         source(i,j,k)=
     &        death_dia(i,j,k)*bio(i,j,k,n,idia)+ 
!#ifdef CHATONELLA
!     &           DEATH_CHA(I,J)*CHA(I,J,K)+ 
!#endif
     &        death_fla(i,j,k)*bio(i,j,k,n,ifla)- 
     &        cc(4)*bio(i,j,k,n,idet)
#ifdef ZOOPL
         source(i,j,k) = source(i,j,k) + zen_d1n(i,j,k) + 
     &               zin_d1n(i,j,k) - d1n_zen(i,j,k) - d1n_zin(i,j,k)
#endif
      enddo
      enddo
      enddo
!
      elseif(icode.eq.isis)then
!
!     SIS
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do  k=1,k_list(i,j)
        source(i,j,k)=
     &  cc(2)*(death_dia(i,j,k)*bio(i,j,k,n,idia)
     &  +red_dia(i,j,k))-scc(4)*bio(i,j,k,n,isis)
#ifdef ZOOPL
        source(i,j,k) = source(i,j,k) + cc(2)*p1n_zen(i,j,k)
#endif
      enddo
      enddo
      enddo
!
      elseif(icode.eq.ioxy)then
!
!     OXY
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j)
         source(i,j,k)= scc(1)*(prod_dia(i,j,k)+prod_fla(i,j,k) 
!#ifdef CHATONELLA
!     &                   +PROD_CHA(I,J,K)-RED_CHA(I,J,K) 
!#endif
#ifdef ZOOPL
     &                 -zen_don(i,j,k) - zin_don(i,j,k)  ! litt usikker her
#endif
     &                  -red_dia(i,j,k)-red_fla(i,j,k) 
     &                  -cc(4)*bio(i,j,k,n,idet))
      enddo
      enddo
      enddo
!
      elseif(icode.eq.ised)then
!
!     SED
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j)
         source(i,j,k)=0.0
      enddo
      enddo
      enddo
!
      elseif(icode.eq.iyel)then
!
!     YELLOW
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j)
         source(i,j,k)=0.0
      enddo
      enddo
      enddo
!
!#ifdef CHATONELLA
!      ELSEIF(ICODE.EQ.icha)THEN
!!
!     CHAtonella
!
!      DO 121 k=1,kdm-1
!      DO 121 j=1-margin,jj+margin
!      DO 121 i=1-margin,ii+margin
!         SOURCE(I,J,K)=PROD_CHA(I,J,K)-RED_CHA(I,J,K)- 
!     &                    DEATH_CHA(I,J)*CHA(I,J,K)
!121   CONTINUE
!#endif
#ifdef DETPHO
      elseif(icode.eq.idetp) then
!
!     DETP
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j)
         source(i,j,k)=cc(1)*( 
     &         0.75*( death_dia(i,j,k)*bio(i,j,k,n,idia)+ 
!#ifdef CHATONELLA
!     &          DEATH_CHA(I,J)*TBIO(I,J,K,ICHA)+ &
!#endif
     &    death_fla(i,j,k)*bio(i,j,k,n,ifla) ) ) - 
     &    1.3*cc(4)*bio(i,j,k,n,idetp)
#ifdef ZOOPL
         source(i,j,k) = source(i,j,k) + 
     &            cc(1)*(zen_d1n(i,j,k) + zin_d1n(i,j,k) - 
     &                   d1n_zen(i,j,k) - d1n_zin(i,j,k))
#endif
      enddo
      enddo
      enddo

#endif
!
#ifdef ZOOPL
      elseif(icode.eq.imeso) then
!
!     MESOZOOPLANKTON
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j)
         source(i,j,k) = p1n_zen(i,j,k) + d1n_zen(i,j,k) + 
     &          zin_zen(i,j,k) - zen_d1n(i,j,k) - zen_don(i,j,k)
      enddo
      enddo
      enddo
!
      elseif(icode.eq.imicro) then
!
!     MICROZOOPLANKTON
!
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      do k=1,k_list(i,j)
         source(i,j,k) = p2n_zin(i,j,k) + p1n_zin(i,j,k) + 
     &                   d1n_zin(i,j,k) - 
     &          zin_zen(i,j,k) - zin_d1n(i,j,k) - zin_don(i,j,k)
      enddo
      enddo
      enddo
#endif
!
      else
!
         write(*,*)'*** FATAL ERROR in SOURCE ***'
         write(*,*)'*** UNKNOWN VALUE of ICODE =',icode
         stop
      endif
      
!
      return
      end subroutine NOR05_source
      end module m_NOR05_source  
