module m_icestate_mechred


contains

   subroutine icestate_mechred
      use mod_icestate
      use mod_icestate_redist
      use mod_icestate_tools
      use mod_icestate_fluxes
      use mod_icestate_thermo
      use m_icestate_calcah
      implicit none

      type(t_ice), dimension(nthick) :: icem, psim
      REAL :: gh,harvest,sigf,ficetot,arealoss,toomuch, ah(0:nthick)
      real :: psi_f1m(nthick), divg, tml
      logical :: ridgm
      integer :: i,j,l,hk
      integer, dimension(nthick) :: kri_eff

#if defined (MPI)
      integer :: ierr
      real*4  :: r4harvest
      include 'mpif.h'
#endif

                                                  



      ! Calculation of a random ridging coefficient (between 2 and 6). It is
      ! assumed that ridging coefficients are distributed according to a  
      ! Gaussian law, centered in avkri and of standard deviation devkri.
      ! Statistics are performed in order to check the repartition of kri
      ! among different possible classes : stat.
      !harvest=.5 ! For testing if OMP is ok
      CALL random_number(harvest)
#if defined (MPI)
      r4harvest=harvest
#if defined (NOMPIR8)
      call mpi_bcast(r4harvest,1,MPI_REAL , 0, MPI_COMM_WORLD,ierr)
#else
      call mpi_bcast(r4harvest,1,MPI_REAL4, 0, MPI_COMM_WORLD,ierr)
#endif
      harvest=r4harvest
#endif

#ifdef MP_TEST_DUMP
      harvest=0.5
#endif
      !print *,harvest
      DO hk = 1, stepfunct-1
         IF (repfunct(hk) .LE. harvest .AND. harvest .LE. repfunct(hk+1)) EXIT 
      END DO
      stat(hk) = stat(hk) + 1
      kri     = 2. + 4. * (repfunct(hk) + repfunct(hk+1)) / 2.
      !print *,'kri is ',kri,h0




      !$OMP  PARALLEL DO PRIVATE(j,l,i,icem,divg,ficetot,sigf,hk,   & 
      !$OMP   ah,ridgm,gh,toomuch,psi_f1m,psim, arealoss,kri_eff,tml) &
      !$OMP SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))

         ! Set things from state etc....
         icem = icestate(i,j)%ice
         tml  = icestate(i,j)%tml
         divg    = divu(i,j) * dtt
         ficetot = SUM(icem%fice)


         ! This routine is called directly after advection. Make sure that 
         ! the thickness distribution is consistent
         psim = clear_ice(psim,nlay)
         call combine_ice(psim,icem,tml)
         icem=psim



         ! Calculation of probability that ice of thickness h should
         ! participate to rafting or ridging.  (ah)
         if ( (divg<epsil1.and.fice_max-ficetot<gstar) .or. (ficetot > fice_max) ) then
            ridgm=.true.
            ah=icestate_calcah(icem)
         else
            ah=0.
            ridgm=.false.
         end if


        
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !---------------START of redistribution calculations---------------!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! First attempt at redistribution. Use the amount specified in the 
         ! ah - distribution
         psi_f1m = 0.
         psim    = clear_ice(psim,nlay)
         if (ridgm) then
            where (icem%hice <= h0) 
               kri_eff = 2
            elsewhere
               kri_eff = kri
            endwhere

            ! Area removed from this iceclass ...
            psi_f1m     = - kri_eff / (kri_eff - 1.) * ah(1:nthick) * abs(divg)

            ! ... results in this iceclass
            psim        = icem
            psim%hice   = kri_eff * icem%hice
            psim%hsnw   = kri_eff * icem%hsnw
            psim%fice   = 1. / (kri_eff - 1.) * ah(1:nthick) * abs(divg)

            ! Check that we do not remove too much ice ...
            where (psi_f1m  < - icem%fice )  
               psi_f1m   = - icem%fice 
               psim%fice = - psi_f1m / kri_eff
            endwhere

         end if




         ! Now check wether the net area loss is sufficient to 
         ! keep ficetot below fice_max:
         arealoss      =  - sum( psi_f1m + psim%fice )  
         toomuch       = max(0.,ficetot+epsil1-fice_max)
         if (toomuch>epsil1.and.arealoss<toomuch) then

            ! Revert to traditional scheme, icestate properties are rebuilt later
            ridgm = .true.
            psi_f1m   = - icem%fice   ! All ice will first be removed in redice
            psim      = icem          ! Rebuilt ice based on icem
            psim%fice = psim%fice * fice_max / ficetot
            psim%hice = psim%hice * ficetot / fice_max
            psim%hsnw = psim%hsnw * ficetot / fice_max
         endif





        ! Finally, for totfice < 0.01 we adjust the thickness so that totfice = 0.01
        if ( ficetot > epsil1  .and. ficetot < fice_min .and. sum(icem%fice*icem%hice) > epsil1 ) then
           where (icem%fice > epsil1)
              psim      = icem          ! New "unridged" ice
              psi_f1m   = -icem%fice    ! All ice will be rebuilt 

              ! "Unridging" process
              psim%hice = psim%hice * psim%fice * ficetot / fice_min
              psim%hsnw = psim%hsnw * psim%fice * ficetot / fice_min
              psim%fice =             psim%fice * fice_min/ ficetot
           end where
           ridgm     = .true.        ! yes for ridging 

           ! Clear ice classes if necessary
           where (psim%fice<epsil1.or.psim%hice<epsil1) psim=clear_ice(psim,nlay)
        else if ( ficetot < fice_min .and. sum(icem%fice*icem%hice) < epsil1 ) then
           icem=clear_ice(icem,nlay)
        end if





        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !-------------------START of actual redistribution-----------------!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (ridgm) then
           
            ! Apply psi_f1, which indicates the quantity of each ice type that 
            ! should disappear :
            icem(:)%fice = icem(:)%fice + psi_f1m(:)

            ! Combine ice in icem and psim distribution to a final 
            ! distribution in icem
            call combine_ice(icem,psim,tml)


            ! Remove empty classes
            where (icem%fice<epsil1.or.icem%hice<epsil1)
               icem=clear_ice(icem,nlay)
            endwhere


            icestate(i,j)%ice = icem
        endif


      enddo
      enddo
      enddo
      !OMP END PARALLEL DO


   end subroutine icestate_mechred
end module m_icestate_mechred
