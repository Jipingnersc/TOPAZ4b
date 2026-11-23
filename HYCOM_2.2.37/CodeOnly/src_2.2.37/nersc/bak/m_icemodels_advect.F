!!****f* HYCOM_2.1.03/icemodels_advect
!!
!! NAME
!!  icemodels_advect - Step ice model one time step forward
!!
!! SYNOPSIS
!!  subroutine icemodels_advect()
!!
!!
!! DESCRIPTION
!!  This is the ice advection routine in NERSC-HYCOM. It is the top-
!!  level logic for advectong scalar fields (hicem, ficem etc).
!!
!!
!! INPUTS   
!!
!! SIDE EFFECTS
!!  called.
!!
!! WARNINGS
!!  Different CPP options will change this routine. CPP options currently
!!  implemented are ICE, EVP_MPI and ICESTATE.
!!
!!
!! PARAMETERS
!!  fice_max        Maximum allowed ice concentration. From module
!!                  mod_common_ice, set in subroutine icedat.
!!
!!  fice_min        Set locally.
!!
!!
!! AUTHOR
!!  Knut Arild Liseter
!!
!! CREATION DATE
!!  25th Feb 2008
!!
!! MODIFICATIONS
!! 9th Jun  2008 - Bug fix in diffusion
!! 9th Jun  2008 - multiplied advection variables with tenm to get
!!                 advected variables on order of 1 
!! SOURCE
!!
      module m_icemodels_advect
      use mod_xc 
      contains


      subroutine icemodels_advect()
#ifdef ICE
      use mod_common_ice
#endif
      use mod_advem
      implicit none 
      integer :: i,j,l,k,oldmargin
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: 
     &   dummy1, snwvol, icevol, iutil3, iutil4, dummy2,
     &   iuflx, ivflx,newicevol, farea, fyarea, fyage, ridgearea,
     &   fyageold,ifld, ficem_old
#if defined (TEST_ICE_AGE)
     &  ,crdg, c_l
      real, parameter :: kridge=4.
#endif
      real :: offset,flxdiv,factor, icedf
      real, parameter :: fice_min=0.05 ! Should be gotten from mod_common_ice
      include 'common_blocks.h'
C --- define step function;  step=1  if  xx1.ge.xx2 ; =0  if xx1.lt.xx2
      real :: xx1, xx2, step
      step(xx1,xx2)=(1.+sign(1.,xx1-xx2))*.5


#if defined (ICE)  

c --- Keep old margin
      oldmargin=margin
c
c --- Ice variables are now initialized on the full margin - we tile later
      margin=nbdy
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jdm+nbdy
         do i=1-nbdy,idm+nbdy
C
            ! NB, use a significant layer thickness , since advem uses pressure coords
            dummy1(i,j)=tenm
            dummy2(i,j)=dummy1(i,j)
            iuflx(i,j)=0.
            ivflx(i,j)=0.
            icevol(i,j)=0.
            snwvol(i,j)=0.
#if defined (TEST_ICE_AGE)
            farea (i,j)=0.
            fyarea(i,j)=0.
            fyage (i,j)=0.
            fyageold(i,j)=0.
#endif
C
         end do
      end do
C$OMP END PARALLEL DO

      ! Set up ice variables as "layered" - so that integrating variable 
      ! over a layer retrieves origial fields. Also correct for flux
      ! divergence. Note that this may be better stated in terms of ice
      ! volume rather than a fake layer thickness.
      margin=nbdy-1 ! due to i+1,i-1
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
C
         ! Set up advection variables
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
C
            ! Now these are "layered", eg farea*dummy1 gives ice conc.
            icevol(i,j)=ficem(i,j)*hicem(i,j)*tenm/dummy1(i,j)
            snwvol(i,j)=ficem(i,j)*hsnwm(i,j)*tenm/dummy1(i,j)
            farea (i,j)=ficem(i,j)*tenm/dummy1(i,j)
#if defined (TEST_ICE_AGE)
            fyarea   (i,j)=fy_frac (i,j)*tenm/dummy1(i,j)
            fyage    (i,j)=fy_age  (i,j)*fy_frac(i,j)*tenm/dummy1(i,j)
            fyageold (i,j)=fy_age  (i,j)
            ridgearea(i,j)=rdg_frac(i,j)*tenm/dummy1(i,j)
#endif /*TEST_ICE_AGE*/
         end do
         end do
C
         ! Set up u-fluxes
         do l=1,isu(j)
         do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            iuflx(i,j)=iceu(i,j)*scuy(i,j)*dummy1(i-1,j)
         end do
         end do
C
         ! Set up v-fluxes
         do l=1,isv(j)
         do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ivflx(i,j)=icev(i,j)*scvx(i,j)*dummy1(i,j-1)
         end do
         end do
C
      end do
C$OMP END PARALLEL DO
C
      ! Set up before and after fake layer thickness 
      margin=nbdy-1 ! due to i+1,i-1
C$OMP PARALLEL DO PRIVATE(j,l,i,flxdiv,offset) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         flxdiv=((iuflx(i+1,j)  -iuflx(i,j)  )
     &          +(ivflx(i,j+1)  -ivflx(i,j)  ))*baclin*scp2i(i,j)
         dummy2(i,j)=dummy1(i,j)-flxdiv
         offset=min(0.,dummy1(i,j),dummy2(i,j))
         dummy2(i,j)=dummy2(i,j)-offset
         dummy1(i,j)=dummy1(i,j)-offset
      end do
      end do
      end do
C$OMP END PARALLEL DO


      ! Tile Variables
      margin=nbdy 
      call xctilr(dummy1(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(dummy2(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(iuflx (1-nbdy,1-nbdy),1,1, 6,6, halo_uv)
      call xctilr(ivflx (1-nbdy,1-nbdy),1,1, 6,6, halo_vv)
      call xctilr(icevol(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(snwvol(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(farea (1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
#if defined (TEST_ICE_AGE)
      call xctilr(ridgearea(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(fyarea   (1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(fyage    (1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
#endif /*TEST_ICE_AGE*/



Cdiag tmp_old_icevol=ficem*hicem*scp2
Cdiag tmp_old_icearea=ficem*scp2
Cdiag call xcsum(old_vol_sum,tmp_old_icevol,ip)
Cdiag call xcsum(old_area_sum,tmp_old_icearea,ip)
       
                  

C      ! For now use the mpdata routine.
C      call advem(2,farea,iuflx,ivflx,
C     &           scp2,scp2i,baclin,dummy1,dummy2,iutil3,iutil4)
C      call advem(2,icevol,iuflx,ivflx,
C     &           scp2,scp2i,baclin,dummy1,dummy2,iutil3,iutil4)
C      call advem(2,snwvol,iuflx,ivflx,
C     &              scp2,scp2i,baclin,dummy1,dummy2,iutil3,iutil4)
C#if defined (TEST_ICE_AGE)
C      call advem(2,fyarea,iuflx,ivflx,
C     &              scp2,scp2i,baclin,dummy1,dummy2,iutil3,iutil4)
C      call advem(2,fyage ,iuflx,ivflx,
C     &              scp2,scp2i,baclin,dummy1,dummy2,iutil3,iutil4)
C      call advem(2,ridgearea,iuflx,ivflx,
C     &              scp2,scp2i,baclin,dummy1,dummy2,iutil3,iutil4)
C#endif /*TEST_ICE_AGE*/

      ! For now use the mpdata routine - which does not 
      ! use middle time step
      call advem(1,farea,farea,iuflx,ivflx,
     &           dummy1,dummy2,0.,scp2,scp2i,baclin )
      call advem(1,icevol,icevol,iuflx,ivflx,
     &           dummy1,dummy2,0.,scp2,scp2i,baclin )
      call advem(1,snwvol,snwvol,iuflx,ivflx,
     &           dummy1,dummy2,0.,scp2,scp2i,baclin )
#if defined (TEST_ICE_AGE)
      call advem(1,fyarea,iuflx,ivflx,
     &              scp2,scp2i,baclin,dummy1,dummy2,iutil3,iutil4)
      call advem(1,fyage ,iuflx,ivflx,
     &              scp2,scp2i,baclin,dummy1,dummy2,iutil3,iutil4)
      call advem(1,ridgearea,iuflx,ivflx,
     &              scp2,scp2i,baclin,dummy1,dummy2,iutil3,iutil4)
#endif /*TEST_ICE_AGE*/


      ! Integrate layered variable to recover original
      ficem=farea*dummy2/tenm
      icevol=icevol*dummy2/tenm
      snwvol=snwvol*dummy2/tenm
#if defined (TEST_ICE_AGE)
      ridgearea=ridgearea*dummy2/tenm
      fyarea   =fyarea   *dummy2/tenm
      fyage    =fyage    *dummy2/tenm
#endif /*TEST_ICE_AGE*/

      ! Artificial ridging
Cdiag tmp_new_icearea=ficem*scp2
      where (ficem > fice_max )
#if defined (TEST_ICE_AGE)

         ! This is the fraction to be ridged
         crdg=ficem-fice_max

         ! This is the level ice to be removed. 
         ! k*h_level is the thickness of ice ridged in this way.
         ! h_level is unknown. We only operate with mean ice thickness
         ! here.
         c_l=crdg*kridge/(kridge-1)

         ! Remove fy area. c_l is removed. max is needed in case
         ! my level ice (or ridged ice) is removed as well. 
         ! Note that both my and fy level
         ! ice is assumed to have the same thickness in this
         ! approximation.
         fyarea    =max(0.,fyarea -   c_l)

         ! increase in ridged area. min is needed in case ridged ice is
         ! re-ridged. (In this case ridged ice is not conserved)
         ridgearea =min(1.,ridgearea + c_l/kridge)
#endif
         ficem =fice_max
      else where (ficem<fice_min .and. ficem>0.)
#if defined (TEST_ICE_AGE)
         fyarea   =fyarea  *   (1. + fice_min-ficem)
         ridgearea=ridgearea * (1. + fice_min-ficem)
#endif /*TEST_ICE_AGE*/
         ficem=fice_min
      end where

      ! NB - Should check the logic here
      where (ficem >= fice_min )
         hicem=icevol/ficem
         hsnwm=snwvol/ficem
      end where

#if defined (TEST_ICE_AGE)
      fy_frac=fyarea
      rdg_frac=ridgearea
      where (fy_frac>1e-3)
         fy_age = fyage/fy_frac
      elsewhere
         fy_age = fyageold
      end where
#endif /*TEST_ICE_AGE*/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diffusion of variables starts here ..
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Use 1/5 "temdf2"  as diffusion velocity
      icedf=temdf2*0.2
      icevol=hicem*ficem
      snwvol=hsnwm*ficem
      ficem_old=ficem
      margin=nbdy-1 ! due to i+1, i-1
      do k=1,3
         if (k==1) then
            ifld=ficem
         elseif (k==2) then
            ifld=icevol
         elseif (k==3) then
            ifld=snwvol
         end if
!$OMP PARALLEL DO PRIVATE(j,l,i,factor)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
           do l=1,isu(j)
             do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
               factor=icedf*aspux(i,j)*scuy(i,j)*
     &                min(step(ficem_old(i  ,j),fice_min),
     &                    step(ficem_old(i-1,j),fice_min))
               iuflx(i,j)=factor*(ifld(i-1,j)-ifld(i,  j))
             enddo
           enddo
           do l=1,isv(j)
             do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
               factor=icedf*aspvy(i,j)*scvx(i,j)*
     &          min(step(ficem_old(i,j  ),fice_min),
     &              step(ficem_old(i,j-1),fice_min))
               ivflx(i,j)=factor*(ifld(i,j-1)-ifld(i,j))
             enddo
           enddo
         enddo
C$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(j,l,i,factor)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
           do l=1,isp(j)
             do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
               factor=-delt1/scp2(i,j)
               util1(i,j)=(iuflx(i+1,j)-iuflx(i,j)
     &                 +ivflx(i,j+1)-ivflx(i,j))*factor
               ifld(i,j)=ifld(i,j)+util1(i,j)
             enddo
           enddo
         enddo
C$OMP END PARALLEL DO

         if (k==1) then
            ficem=ifld
         elseif (k==2) then
            icevol=ifld
         elseif (k==3) then
            snwvol=ifld
         end if
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final adjustment of fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Adjustment for small numbers after diffusion
      !where (ficem<fice_min .and. ficem>fice_tol)
      !   ficem=fice_min
      !elsewhere (ficem<=fice_tol)
      !   ficem=0.
      !   hsnwm=0.
      !   hicem=0.
      !end where
      where (ficem<fice_min .and. ficem>0.)
         ficem=fice_min
      end where
      where (ficem >= fice_min )
         hicem=icevol/ficem
         hsnwm=snwvol/ficem
      end where





      ! Tile final variables
      call xctilr(ficem( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(hicem( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(hsnwm( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
#if defined (TEST_ICE_AGE)
      call xctilr(fy_age  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(fy_frac ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(rdg_frac( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
#endif /*TEST_ICE_AGE*/

      

      ! Restore old margin
      margin=oldmargin
#elif defined (ICESTATE) 
      if (mnproc==1) print *,'Advection of icestate vars not yet...'
      call xcstop('(icemodels_step)')
      stop '(icemodels_step)'
#endif /* ICE or ICESTATE */


      end subroutine icemodels_advect
      end module m_icemodels_advect
!!***
