module mod_CAL06_calnusparticles

   use mod_xc, ONLY : idm, jdm, kdm, nbdy
  
   implicit none
   include 'comcal.h'
   type calfin
      integer             :: num       ! float identity number
      integer             :: type      ! type of float
      integer             :: layer     ! number of model layer
      real, dimension(3)  :: pos       ! float 3D position
      real                :: starttime ! initiation time
      real                :: endtime   ! termination time
      real                :: wdepth    ! water-depth
      real                :: temp,salt ! temperature and salinity 
      real                :: dens      ! density
      real                :: age
      real                :: mirr                ! irradiation at pos
      real                :: mixl                ! mixed layer depth
      real                :: phy_conc            ! phytoplanktopn concentration
      real                :: depth_max_phy_conc  ! depth of maximum phytoplankton concetration
      integer :: fi

!Strategy codes
!sv1 = wud
!sv2 = fsr
!sv3 = afd
!sv4 = owd
!sv5:maxs = ANN weights
!bsv1:nhid+1 = ANN weights
      real :: sv(6)!sv(1:maxs)
!      real :: bsv(1+2)

!Attribute codes
!av0 = 0=dead, 1=alive
!av1 = Stage, 0 is egg, 1-6 nauplia, 7-11 copepodits, 12 adult
!av2 = Internal number in individuals
!av3 = Structural weight  in micrograms
!av4 = Stage longevity
!av5 = Fat energy level in KJ
!av6 = Moult cycle fraction (Egg,N2,N3) | cum. egg number (adult)
!av7 = Total number of eggs in a super individual
!av8 = habd 
!av9 = habe 
!av10 = habn
!av11 = diapause,0 diapause,1 active,2 move down,3 move up
!av12 = Summing of daily growth  
!av13 = Summing of hourly egestion
!av14 = Unique identifier
!av15 = neutral overwintering depth
!av16 = growth rate
!av17 = antall avkom 
      real :: av(0:maxa)    !attribute vector
      
      real :: prisk(4)
      integer :: nsip                    ! number of superindividuals produced
      integer :: generation              ! generation 
      integer :: gender   !0=male, 1=female
   end type calfin

   integer, parameter :: maane(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/) ! Day in month
!stage longivity parameter of Lynch et al. 1998   
   real, parameter :: a(0:12)=(/595., 387., 582.,1387., 759.,&
                                716., 841., 966.,1137.,1428.,&
                               2166.,4083.,4083./) 
!stage longivity parameter of Campbell et al. 2001  
!   real, parameter :: a(0:12)=(/595.,  983., 1564.,2951.,3710.,&
!                               4426., 5267., 6233.,7370.,8798.,&
!                               10964.,15047.,15047./) 
!stage treshhold in ug C from Carlotti & Radach, 1996
   real, parameter :: st_wgt(3:13)=(/ 0.05, 0.20, 0.30,  0.45,  0.75, &
                                      1.10, 2.50, 7.00, 15.00,  90.00, &
                                      40.00/)
!Mortality data per stage (Aksnes & Blindheim 1996)
!Adult mortality from Miller et al. 1998 anf Ohman 2004 
!AS  real, parameter :: mort(0:13)=(/0.16, 0.08,  0.08,   0.08, 0.08, &
!AS                                   0.08, 0.08,  0.16,   0.18, 0.15, & 
!AS                                   0.04, 0.005, 0.0253, 0.0253/)
   real, parameter :: mort(0:13)=(/0.30,  0.15,  0.10,  0.02,  0.02, &
                                   0.02, 0.02,   0.025,   0.02, 0.02, & 
                                   0.025, 0.03,   0.03,   0.03/)
   real, parameter :: diapause_mort=0.001

!   real, parameter :: mort(0:13)=(/0.3,  0.3,    0.3,    0.3,  0.05, &
!                                   0.05, 0.05,   0.05,   0.05, 0.02, & 
!                                   0.02, 0.02,   0.02,   0.02/)

! from comcal.h
   integer, parameter :: mpop=49000  ! maximum number of sythetic floats
   integer,parameter:: offsp=mpop*0.25
   real :: prisk(4,mpop)

!CAS mort(13:14) satt inn 19.06.06 
!   real, parameter :: maxegg = totcop*2./(mpop*1.)
   real, parameter :: maxegg = 2.0e+12
!   real, parameter :: maxegg = 1.3333333e+12  ! for run with 30000 ind in NOMO

   real :: totbio
   real :: cmort

   type (calfin) :: cfs(mpop+1)
   integer       :: pop
   integer       :: sumdelta,pops
   real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: radmax
   real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: mod_carb
   real, parameter :: biodt_cal=1./24. ! 1 time

   character(len=3), parameter :: pm_rid='NW1'

!idealized or realistic chlorophyll concetration
   logical, parameter :: idealized=.false.  

   contains
!*********************************************************************
!!$   subroutine vertical_migration(cf, deltfl, rt)
!!$      use mod_year_info
!!$      implicit none
!!$      type (year_info) :: rt
!!$      type (calfin)    :: cf
!!$      real             :: deltfl
!!$      
!!$      if ((real(rt%idd).ge.ascent_time).and. &
!!$          (real(rt%idd).le.decent_time)) then
!!$         if (cf%pos(3).gt.cf%upper_depth) then
!!$           cf%pos(3)=cf%pos(3)-ascent_speed*deltfl
!!$        end if
!!$      end if
!!$      
!!$      if (real(rt%idd).gt.decent_time) then
!!$         if (cf%pos(3).lt.cf%lower_depth) then
!!$           cf%pos(3)=cf%pos(3)+decent_speed*deltfl
!!$        end if
!!$      end if
!!$      
!!$   end subroutine vertical_migration
!*********************************************************************
subroutine write_calanus_to_file_new(rt,nflt,rid)
    use mod_year_info

    implicit none
    type(year_info) :: rt
    character*3     :: rid
    character*70    :: calfile
    integer, intent(in) :: nflt
    integer         :: nfl

    calfile=rid//'CAL_y'//rt%cyy//'_d'//rt%cdd//'_h'//rt%chh//'.uf'

    open(802,file=trim(calfile),form='unformatted')
    write(802) nflt
    do nfl=1,nflt
       write(802) cfs(nfl), rt%iyy, rt%idd 
    end do
    close(802)

!!! for now (mortality file)    
    open(803,file='Biol_mortality.dat',form='unformatted',position='append')
    write(803) rt%iyy,rt%idd
    write(803) zm_cal
    write(803) zn_cal
    write(803) zb_cal
    close(803)

end subroutine write_calanus_to_file_new
!*********************************************************************
  subroutine print_calanus_info(calanus)
     type (calfin) :: calanus
 
201 format(a27,i6)
202 format(a27,f8.4,a8,f8.4, a8,f12.4)
203 format(a27,f12.4)
204 format(a22,f8.4,a22,f8.4, a22,f8.4)
205 format(a27,f8.4,a8,f8.4, a8,e12.4)
206 format(a27,e12.4)

     write(*,201) 'Calanus identity number      ', calanus%num
     write(*,201) 'Calanus type                 ', calanus%type
     write(*,201) 'Calanus layer                ', calanus%layer
     write(*,202) 'Calanus position: Lon        ', calanus%pos(1), &
                                  'Lat',   calanus%pos(2), &
                                  'Depth ',calanus%pos(3)
     write(*,203) 'Calanus start-time           ', calanus%starttime
     write(*,203) 'Calanus end-time             ', calanus%endtime
     write(*,203) 'Calanus position water depth:', calanus%wdepth
     write(*,204) 'Calanus temperature:', calanus%temp, &
             'Calanus salinity:', calanus%salt, &
             'Calanus density:', calanus%dens

     write(*,203) 'Calanus age                  ', calanus%age
     write(*,203) 'Max irradiance at position   ', calanus%mirr 
     write(*,203) 'Mix Layer Depth              ', calanus%mixl 
     write(*,203) 'Phytoplankton concentration  ', calanus%phy_conc 
     write(*,203) 'Depth of max phytopl. conc.  ', calanus%depth_max_phy_conc 
     write(*,201) 'Calanus fitness              ', calanus%fi
!sv1 = wud
!sv2 = fsr
!sv3 = afd
!sv4 = owd
!sv5:maxs = ANN weights
!bsv1:nhid+1 = ANN weights
     write(*,203) 'Wake up day                  ',calanus%sv(1)
     write(*,203) 'Fat soma ratio               ',calanus%sv(2)
     write(*,203) 'SV(3)                        ',calanus%sv(3)
     write(*,203) 'Overwintering depth          ',calanus%sv(4)
     write(*,203) 'Diurnal migration param. 1   ',calanus%sv(5)
     write(*,203) 'Diurnal migration param. 2   ',calanus%sv(6)
!      real :: bsv(1:nhid+2)

!Attribute codes
!av0 = 0=dead, 1=alive
!av1 = Stage, 0 is egg, 1-6 nauplia, 7-11 copepodits, 12 adult
!av2 = Internal number in individuals
!av3 = Structural weight  in micrograms
!av4 = Stage longevity
!av5 = Fat energy level in KJ
!av6 = Moult cycle fraction (Egg,N2,N3) | cum. egg number (adult)
!av7 = Total number of eggs in a super individual
!av8 = habd 
!av9 = habe 
!av10 = habn
!av11 = diapause,0 diapause,1 active,2 move down,3 move up
!av12 = Summing of daily growth  
!av13 = Summing of hourly egestion
!av14 = Unique identifier
!av15 = neutral overwintering depth
!av16 = growth rate
!av17 = antall avkom 
     write(*,203) 'dead or alive?               ',calanus%av(0)
     write(*,203) 'Stage                        ',calanus%av(1)
     write(*,206) 'Internal # in individuals    ',calanus%av(2)
     write(*,203) 'Structural weight  in ug     ',calanus%av(3)
     write(*,203) 'Stage longevity              ',calanus%av(4)
     write(*,203) 'Fat energy level in KJ       ',calanus%av(5)
     write(*,203) 'Moult cycle fraction         ',calanus%av(6)
     write(*,203) 'Total # of eggs in a SI      ',calanus%av(7)
     write(*,203) 'Depth                        ',calanus%av(8)
     write(*,203) 'Longitude                    ',calanus%av(9)
     write(*,203) 'Latitude                     ',calanus%av(10)
     write(*,203) '0 diap.,1 active,2 down,3 up ',calanus%av(11)
     write(*,203) 'Summing of daily growth      ',calanus%av(12)
     write(*,203) 'Summing of hourly egestion   ',calanus%av(13)
     write(*,203) 'Unique identifier            ',calanus%av(14)
     write(*,203) 'Neutral owerwintering depth  ',calanus%av(15)
     write(*,203) 'growth rate                  ',calanus%av(16)
     write(*,203) 'Number of eggs               ',calanus%av(17)
     write(*,203) 'AV(18)                       ',calanus%av(18)
     write(*,203) 'AV(19)                       ',calanus%av(19)


     write(*,203) 'prisk(1)                     ',calanus%prisk(1)
     write(*,203) 'prisk(2)                     ',calanus%prisk(1)
     write(*,203) 'prisk(3)                     ',calanus%prisk(1)
     write(*,203) 'prisk(4)                     ',calanus%prisk(1)
     write(*,201) 'Number of super-ind. prod.   ',calanus%nsip
     write(*,201) 'Generation                   ',calanus%generation
     write(*,201) 'Gender                       ',calanus%gender

  end subroutine print_calanus_info
!*********************************************************************

end module mod_CAL06_calnusparticles
!******************************************************************************

