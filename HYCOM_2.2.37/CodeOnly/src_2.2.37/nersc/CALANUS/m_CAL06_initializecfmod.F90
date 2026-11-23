      module m_CAL06_initializecfmod
      use mod_CAL06_calnusparticles
      use mod_xc    ! HYCOM communication interface
      use mod_year_info

      contains


      subroutine initializecfmod(rt,nflt,rungen,fltflg)

      implicit none 

      type(year_info), intent(in) :: rt
      character(len=3), intent(in) :: rungen
      integer, intent(in) :: nflt,fltflg
      integer :: i,t,b,j,newpop,x,rpop,ye,cd(6,1992:2000)
      real	   :: nevn,ze,norm,expo,rand,nmax,cum,sd
      integer :: iostat, nfl
      real, parameter :: ER = 6378135.0

      real :: tmp
  
      iostat=0 

     if (fltflg.eq.1) then
        print*, 'New output files for calanus are generated'
        open(801,file='calanus_out',status='replace',form='formatted')
        close(801)
        open(801,file='mortality_out',status='replace',form='formatted')
        close(801)
        open(801,file='Biol_mortality.dat',status='replace',form='unformatted')
        close(801)
     endif

     call td(0,cmort)

   do nfl=1,nflt
      if (nfl.gt.mpop) then
         print*, 'Dimension for floats array is too small'
         print*, 'mpop:',mpop, 'nfl:',nfl
         call xcstop('(initializecfmod)')
      end if
  !fix    cfs(nfl)%physprop=flts(nfl) !done in mode_floats
      cfs(nfl)%mirr=0
      cfs(nfl)%fi=0          
   end do

   close(41) 
   newpop = mpop*0.4 !pop
!predation density
     call preddens  

   if (readpop==0) then ! initilize new population
!Create chromosomes for copepods

! Set to reasonable values
      do nfl=1,nflt
         call random_number(rand)
      !AS   cfs(nfl)%sv(1) = calasc + 2.0*(rand - 0.5)*calstd  ! Ascent day
         tmp=(2.0*(exp(-(rand-0.5)**2.0/(2.0*(0.42**2.0)))-0.5))
         call random_number(rand)
         cfs(nfl)%sv(1)=calasc + 2.0*(rand-0.5)*calstd*tmp
         cfs(nfl)%sv(2) = 0.4  ! Fat-soma-rato of descent
         cfs(nfl)%sv(3) = 144.0 ! Overwinter moulting
         cfs(nfl)%sv(4) = cfs(nfl)%pos(3) ! Overwintering depth AS0707, changed from 780
         cfs(nfl)%sv(5) =    13.0 ! Ny: 13.0  Old:  30.0 !Overwintering depth
         cfs(nfl)%sv(6) = 20500.0 ! Ny: 20500 Old: 158.0 !Overwintering depth
      enddo 

!Create attributes of super-individuals
      do nfl = 1,nflt
         cfs(nfl)%age=150.0

         cfs(nfl)%av(0) = 1.0
         cfs(nfl)%av(1) = 12.0
         cfs(nfl)%av(2) = maxegg*0.5
         cfs(nfl)%av(3) = 120.0  !in microgram C !120
         cfs(nfl)%av(4) = 150.0 !150.
         cfs(nfl)%av(5) = 20.
         cfs(nfl)%av(6) = 0.
         cfs(nfl)%av(7) = 0.
         cfs(nfl)%av(8) = cfs(nfl)%pos(3) !cfs(nfl)%sv(4)
         cfs(nfl)%av(9) = cfs(nfl)%pos(1)
         cfs(nfl)%av(10) = cfs(nfl)%pos(2)
         cfs(nfl)%av(11) = 0.
         cfs(nfl)%av(12) = 0.
         cfs(nfl)%av(13) = 0.
         cfs(nfl)%av(14) = cfs(nfl)%num
         cfs(nfl)%av(18) = cfs(nfl)%av(2)
         cfs(nfl)%av(19) = 0.0
         cfs(nfl)%nsip   = 0
      cfs(nfl)%generation = 0
! make sure gender is the same in both inner and outer model
         call random_number(rand)
         if (rand.gt.0.5) then
           cfs(nfl)%gender = 1
         else
           cfs(nfl)%gender = 0
         end if
         totbio = totbio + cfs(nfl)%av(2)*(cfs(nfl)%av(5) + cfs(nfl)%av(3)*4.)  !Convert to wet weight
      enddo    

      call print_calanus_info(cfs(1))
      totbio = totbio*1.E-18 !Convert to millions of tonnes
      print*, 'Totbio: ',totbio

!      call biomass
#if defined (NOMODEL)
         if (initfl==0) then 
!fix           call write_boundary_calanus(rt,rungen)
         end if
#endif ! /* NOMODEL */

      elseif(readpop==1) then 

!Read adapted population from file
         open(17, file = 'library/'//'sv_all.txt')	          
         open(188, file = 'library/'//'av_all.txt')	          
         open(199, file = 'library/'//'bsv_all.txt')	          
         read(188,*)rpop
         pop = min(rpop,mpop)
         print*,'Number of released super individuals:',pop
         newpop = pop
         do i = 1,pop
            read(17,*) (cfs(i)%sv(i), j = 1,maxs)     !sv_all.txt
            read(188,*)(cfs(i)%av(j), j = 0,maxa)	   !av_all.txt
!            read(199,*)(cfs(i)%bsv(j),j = 1,3) !bsv_all.txt
         enddo
         close(17)
         close(188)
         close(199)
      endif
      pop = newpop

!Cumulative normal distribution
      sd = 0.2*maxt
      nevn = sqrt(6.28318*sd)
      norm =	(1./nevn)
      nmax = 0.
      cum = 0.
      do t = 1,maxt
         ze = ((t-int(0.18*maxt))**2.)/(2*sd**2.)
         expo =  exp(-ze)
         npar(t) = norm*expo
         cum = cum + npar(t)
         cumn(t) = cum 
         if(npar(t)>nmax) nmax = npar(t)
      enddo

!Print
!open(27, file = 'out/'//'hormon.txt')
!   do t = 1,maxt
!      npar(t) = npar(t)/nmax
!      cumn(t) = cumn(t)/cum
!      !write(27,*)t,npar(t),cumn(t)
!   enddo
!close(27)

      do x = 0,13
         avpred(x) = 0.
         numpred(x) = 0.
         agemolt(x) = 0.
         molters(x) = 0.
      enddo
      pop=nflt
 
!Calculate mackerel predation probability
!   open(27, file = 'in/'//'makrelldatoer.txt')
!!   do ye = 1992,2000
!      read(27,*)cd(1:6,ye)
!   enddo
!   close(27)
!   do ye = 1992,2000
!      sd = 0.2*cd(4,ye)
!      nevn = sqrt(6.28318*sd)
!      norm =	(1./nevn)
!      nmax = 0.
!      do t = cd(2,ye),cd(3,ye)
!!         ze = ((t-cd(5,ye))**2.)/(2*sd**2.)
!         expo =  exp(-ze)
!         npar(t) = norm*expo
!         if(npar(t)>nmax) nmax = npar(t)
!      enddo
!      do t = 1,365
!         if(t<cd(2,ye)) then 
!!            mackprob(t,ye) = 0.
!         elseif(t<cd(3,ye)) then 
!            mackprob(t,ye) = npar(t)*cd(6,ye)/11.44 !Scale to max 1.
!         else 
!             mackprob(t,ye) = 0.
!         endif
!      enddo
!   enddo

!Read herring predation probability
!   open(27, file = 'in/'//'herprob.txt')
!   do t = 1,365 
!      read(27,*)herprob(t,1992),herprob(t,1993),herprob(t,1994),herprob(t,1995),&
!           herprob(t,1996),herprob(t,1997),herprob(t,1998),herprob(t,1999) !,herprob(t,2000)
!   enddo
!   close(27)
!   do ye = 1992,2000
!      do t = 1,365 
!         herprob(t,ye) = herprob(t,ye)/14.33  !Scale to max 1.
!       !if(ye==1993) print*,ye,t,herprob(t,ye)
!      enddo
!   enddo

!pen(27, file = 'out/'//'hormon.txt')
!do ye = 1992,1999
! do t = 1,365
!  write(27,*)ye,t,herprob(t,ye),mackprob(t,ye)
! enddo
!enddo
!close(27)
!stop
!fix   if (appflt==0) then
!      call write_calanus_to_file(rt)
!      call write_calanus_to_file_new(rt,rungen)
!   end if

      print*,'Finished initialising data!' 
   
      open(28,file=rungen//'calanus_param_and_initcond.txt',form='formatted')
      write(28,*) 'Parameters set in comcal.h'
      write(28,*)
      write(28,*) '1.  readpop  = ', readpop, ' If == 0 initialize new population, see also initfl'
      write(28,*) '                  If == 1 read adpted population - does not work'
      write(28,*)   

      write(28,*) '2.  mpop     = ', mpop, ' Maximum number of superindividuals'
      write(28,*)

      write(28,*) '3.  maxt     = ', maxt, ' Number of days in a year'
      write(28,*) '3.  offsp    = ', offsp, ' Maximum number of new superindividuals in on timestep'
      write(28,*)

      write(28,*) '4.  mthresh  = ', mthresh, ' If internal number of individuals are above this'
      write(28,*) '                             number mortality is done by reduction in numbers'
      write(28,*)
   
      write(28,*) '5.  maxs  = ', maxs, ' Maximun number of strategies'
      write(28,*) '5.  maxa  = ', maxa, ' Maximum number of attributes'
      write(28,*)

      write(28,*) '6.  mute    = ', mute, ' Chance of a mutation happening'
      write(28,*) '6.  effect  = ', effect, ' Magnitude of mutation'
      write(28,*) '6.  cmat    = ', cmat, ' Minimum carbon content for maturation'
      write(28,*) '6.  cs_cal  = ', cs_cal, ' Clutch size'
      write(28,*)

      write(28,*) '7.  eggw       = ', eggw, ' Egg weight (ugC)'
      write(28,*) '7.  diff_cal   = ', diff_cal, ' Diffuse light attenuation rate (1/m)'
      write(28,*)

      write(28,*) '8.  ksur       = ', ksur, ' Loss through surface (dim.less)'
      write(28,*) '8.  contr      = ', contr, ' Prey light absorption-reflectance coeff.'
      write(28,*) '8.  beam       = ', beam, ' Parameter in calculating visual range'
      write(28,*)

      write(28,*) '9.  sthresh    = ', sthresh, ' Sensitivity threshold for eye  (uE m-2 s-1)'
      write(28,*) '9.  red        = ', red, ' Redfield ration'
      write(28,*) '9.  card       = ', card, ' Density dependent feeding parameter (g Cm-2)'
      write(28,*)

      write(28,*) '10. pi_cal     = ', pi_cal, ' pi'
      write(28,*) '10. cchl       = ', cchl, ' Carbon to chloropgyll ratio'
      write(28,*) '10. mld        = ', mld, ' Mixed layer depth -- fixed only used in calrep'
      write(28,*)

      write(28,*) '11. tpdens     = ', tpdens, ' Max value for tactile predators'
      write(28,*) '11. tr         = ', tr, ' Tactile predator radius'
      write(28,*) '11. tpvel      = ', tpvel, ' Tactile predator swimming speed'
      write(28,*) '11. tacw       = ', tacw, ' Tactile predator weight? something'
      write(28,*)

      write(28,*) '12. mpvel      = ', mpvel, ' Mesopelagic predator swimming speed'
      write(28,*) '12. mesl       = ', mesl, ' Mesopelagic predator lenth'
      write(28,*) '12. mesw       = ', mesw, ' Mesopelagic predator weight? something'
      write(28,*) '12. mpdens     = ', mpdens, ' Max value for mesopelagic predators'
      write(28,*)
      write(28,*) ' ******************************************************** '
      write(28,*) 'Parameters set in mod_CAL06_calnusparticles.F90'
      write(28,*)
   
      write(28,*) 'Stage longivity parameter of Lynch et al. 1998   '
      do i=0,12
         write(28,*) 'Stage: ', i, ' a: ', a(i)
      end do
      write(28,*)

      write(28,*) 'Stage treshhold in ug C from Carlotti & Radach, 1996'
      do i=3,13
         write(28,*) 'Stage: ', i, ' a: ', st_wgt(i)
      end do
      write(28,*)

      write(28,*) 'Mortality data per stage (Aksnes & Blindheim 1996) -- not used'
      do i=0,12
         write(28,*) 'Stage: ', i, ' a: ', mort(i)
      end do
      write(28,*)
   
      write(28,*) 'Maxegg     = ', maxegg, ' Twice of the initial number of individuals in the superindividuals'
      write(28,*) 'biodt_cal = ', biodt_cal, ' Timestep (days)'
      write(28,*)
      write(28,*) 'pm_rid    = ', pm_rid, ' Productivity model rungen'
      write(28,*)
      write(28,*) 'idealized = ', idealized, ' Idealized or realistic food availability'
      write(28,*)
      write(28,*)

      write(28,*) 'Initial values set in m_CAL06_initializecfmod.F90'
      write(28,*)
      write(28,*) 'Initial strategy vector'
      write(28,*) 'calasc = ',calasc, ' Central ascent day'
      write(28,*) 'calstd = ',calstd, ' Spread in ascent day'
      write(28,*) 'sv(2)  = ',cfs(1)%sv(2), ' Fat-soma-rato of descent'
      write(28,*) 'sv(3)  = ',cfs(1)%sv(3), ' Overwinter moulting day'
      write(28,*) 'sv(4)  = ',cfs(1)%sv(4), ' Overwintering depth'
      write(28,*) 'sv(5)  = ',cfs(1)%sv(5), ' Parameter 1 for day-time depth'
      write(28,*) 'sv(6)  = ',cfs(1)%sv(6), ' Parameter 2 for day-time depth'

      write(28,*)
      write(28,*) 'Initial attrubute vector'
      write(28,*) 'Age    = ',cfs(1)%age, ' days'
      write(28,*) 'av(0)  = ',cfs(1)%av(0), ' 1=alive, 0=dead'
      write(28,*) 'av(1)  = ',cfs(1)%av(1), ' Stage'
      write(28,*) 'av(2)  = ',cfs(1)%av(2), ' Number of individuals'
      write(28,*) 'av(3)  = ',cfs(1)%av(3), ' Carbon content'
      write(28,*) 'av(4)  = ',cfs(1)%av(4), ' Age '
      write(28,*) 'av(5)  = ',cfs(1)%av(5), ' Fat content'
      write(28,*) 'av(6)  = ',cfs(1)%av(6), ' Development in stage '
      write(28,*) 'av(7)  = ',cfs(1)%av(7), ' Fitness - how many eggs can they produce '
      write(28,*) 'av(8)  = ',cfs(1)%av(8), ' Depth'
      write(28,*) 'av(9)  = ',cfs(1)%av(9), ' Longitude'
      write(28,*) 'av(10) = ',cfs(1)%av(10), ' Latitude'
      write(28,*) 'av(11) = ',cfs(1)%av(11), ' Flag for vertcal migration status'
      write(28,*) 'av(12) = ',cfs(1)%av(12), ' Food - feeding '
      write(28,*) 'av(13) = ',cfs(1)%av(13), ' Egestion'
      write(28,*) 'av(14) = ',cfs(1)%av(14), ' Uniq individual number'
      write(28,*) 'av(15) = ',cfs(1)%av(15), ' Ingestion'
      write(28,*) 'av(16) = ',cfs(1)%av(16), ' Growth rate ( growth/individual weight)'
      write(28,*) 'av(17) = ',cfs(1)%av(17), ' Not sure - somthing with reproduction'
      write(28,*) 'av(18) = ',cfs(1)%av(18), ' Not sure - somthing with reproduction'
      write(28,*) 'av(19) = ',cfs(1)%av(19), ' Not sure - somthing with reproduction'
      close(28)

      end subroutine initializecfmod

      end module m_CAL06_initializecfmod
