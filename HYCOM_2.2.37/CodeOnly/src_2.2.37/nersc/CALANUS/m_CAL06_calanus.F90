module m_CAL06_calanus
use mod_CAL06_calnusparticles
use mod_year_info
use mod_CAL06_read_norwecomcarbon

contains

subroutine calanus(deltfl,rt,nflt,rungen)
    implicit none

    type (year_info),intent (in) :: rt
    character(len=*),  intent(in) :: rungen
    integer, intent(inout) :: nflt
    real             :: deltfl
    integer          :: nfl


!!$!This program calculates movement,stagedynamics and reproduction of a Calanoid copepod
!!$. 
!!$!The model is based on super individual representation.
!!$! 
!!$!By Geir Huse, Institute of Marine Research, Bergen.
!!$!
!!$!-----------------------------------------------------------
!!$
!!$program calmod
!!$
!!$!-----------------------------------------------------------
!!$implicit none                                                                      
!!$integer  tid,month,maane(12),gen,pop,year,sim,rep,rf
!!$integer  tell,sumdelta,day,hour
!!$real     pops,cmort,maxegg,biom,herdep!!$!character (len=7) dir
!!$
!!$include  'comcal.txt' 
!!$
!!$print*,'Initialising calanus model!'
!!$
!!$do sim = 3,3 !1-11   
!!$do rep = 1,1
!!$pop = mpop*0.5
!!$call srand(65949.226)
!!$call initialise(maane,pop,maxegg)        !Initialise values for first Timestep
!!$call openfiles(sim)                                         
!!$call readstm  
!!$year = 1991
!!$do gen = start,horizon
!!$ year = year + 1   
!!$ if(year<1992) year = 1992
!!$ if(year>1999) year = 1992 
!!$ call td(gen,cmort)
!!$ !call simspec(sim)
!!$ tid = 0
!!$ write(*,1)sim,rep,gen,svm(1),svm(2),svm(3),svm(4),avm(2),pop,sumdelta  
!!$ write(4,1)sim,rep,gen,svm(1),svm(2),svm(3),svm(4),avm(2),pop,sumdelta  
!!$ do month = 1,12          !Simulation period from January to December
!!$  print*,year,month,pop,sumdelta,biom
!!$  sumdelta = 0
!!$  do day = 1,maane(month)    !Days of the month     
!!$   tid = tid + 1 
!!$   rf = rf + 1
!!$   tell = tell + 1
!!$   call preddens
!!$   do hour = 1,24                  !Timesteps within the day
!!$    call vertmov(tid,gen,pop,cmort,year,hour,herdep)  
!!$        call calsta(tid,pop,hour,year)  
!!$    call calmort(pop,gen,hour)  
!!$   enddo  !Hour loop
!!$   call popzip(pops,pop)
!!$   call printres(gen,tid,pop,day,biom)
!!$   call calrep(pop,pops,sumdelta,maxegg)   
!!$  enddo         !Day loop
!!$ enddo      !Month loop
!!$2 continue
!!$ enddo       !Year loop
!!$enddo
!!$
!!$call closefiles  
!!$
!!$enddo
!!$
!!$1 format(i2,x,i4,x,i4,4(f8.2,x),f15.1,x,i6,x,i4)
!!$
!!$end     program calmod! year = 1991

! update age
   cfs(1:nflt)%age=cfs(1:nflt)%age+deltfl
! DO ONLY ONCE A YEAR 
   if (rt%idd==0.and.rt%ihh==0) sumdelta=0

!Moved to mod_floats.F
!!READ IN DAILY CARBON FRM NORWECOM ONCE A DAY
!    if (.not.idealized) then
!    if (rt%ihh==0.and.rt%iss==0) then
!       call read_daily_carbon(rt)
!    end if
!    end if

!  Calculate biomass once a day
 !   if (rt%ihh==12.and.rt%iss==0) then
 !      call biomass
 !   end if

! DO EVERY HOUR

! do vertical migration
 
        if (.false.) then
        do nfl=1,nflt
        if (cfs(nfl)%num == 2) then
           call print_calanus_info(cfs(nfl))
        end if
        end do
        end if


    call vertmov(rt%idd,1,rt%iyy,rt%ihh) 

    call calsta(rt%idd,rt%ihh)

    do nfl=1,nflt
      cfs(nfl)%pos(1)=cfs(nfl)%av(9) 
      cfs(nfl)%pos(2)=cfs(nfl)%av(10)
      cfs(nfl)%pos(3)=cfs(nfl)%av(8) 
    end do

    call calmort(1,rt%idd)

      print*,'Mortality',cmort,rt%ihh, rt%idd, zm_cal(12,4),zn_cal(12,4),zb_cal(12,4)
      do nfl=1,nflt
        if (cfs(nfl)%av(0).eq.0.) then
           cfs(nfl)%pos(1:2)=-999.0 
        end if 
      end do

    if (rt%ihh==0)  call calrep(nflt)

    if (rt%ihh==0.and.rt%iss==0) then
      call write_calanus_to_file_new(rt,nflt,rungen)
    end if

#if defined (NOMODEL)
!If in NOMODEL write to file the individuals in the NWSNO model area
    if (rt%ihh==0.and.rt%iss==0) then
    print*, 'calling routine for writing boundary conditions'
    call write_boundary_calanus(rt,rungen)
    end if
#endif

#if defined (NWSNO)
!If in NWSNO read in individuals from the NOMODEL code.
    if (rt%ihh==0.and.rt%iss==0) then
    call read_boundary_calanus(rt,rungen)
    end if
#endif

end subroutine calanus
end module m_CAL06_calanus
   
