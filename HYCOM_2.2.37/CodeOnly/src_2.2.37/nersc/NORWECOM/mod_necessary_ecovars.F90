module mod_necessary_ecovars

  use mod_xc
  implicit none 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  * * * I M P O R T A N T * * *
!
!  The variables in this module must be defined for ANY biochemical model coupled to HYCOM, 
!  i.e., same names, and with nbio set to the correct number of compartments in your model
!  and eco_tag set to a 5 character tag for your model.
!
!  Then, in mxkpp, tsadvc, regrid_tracer, read_write_restart_bio, and hycom 
!  (+ your specific ecosystem routines), 
!  make a use statement to the correct version of this module...
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef EVA85
  integer, PARAMETER :: nbio=3         ! number of biochemical compartments
  character(len=5) :: eco_tag='EVA85'      ! tag for the current ecomodel
#endif  
#ifdef FDM02
  integer, PARAMETER :: nbio=11         ! number of biochemical compartments
  character(len=5) :: eco_tag='FDM02'      ! tag for the current ecomodel
#endif  
#ifdef SCH02
  integer, PARAMETER :: nbio=15         ! number of biochemical compartments
  character(len=5) :: eco_tag='SCH02'      ! tag for the current ecomodel
#endif  
#if defined (NOR05)
#if defined (ZOOPL)
  integer, PARAMETER :: nbio=11         ! number of biochemical compartments
#elif defined (DETPHO)
  integer, PARAMETER :: nbio=9         ! number of biochemical compartments
#else
  integer, PARAMETER :: nbio=8         ! number of biochemical compartments
#endif /* DETPHO-ZOOPL */
  character(len=5) :: eco_tag='NOR05'      ! tag for the current ecomodel
  integer, PARAMETER :: init= 1, ipho = 2, isil = 3,  idet  = 4,  isis= 5, &
                        ifla= 6, idia = 7, ioxy = 8,  ised  = 12, iyel=13, &
                        idetp=9, imeso=10, imicro=11, icha=14;
!AS15062010: 
!The mortality cc(3) maqybe a bit high when ZOOPL is not defined
!   CC(1)    : FRACTION OF PHOSPHATE AND NITRATE IN A CELL
!   CC(2)    : FRACTION OF SILICATE AND NITRATE IN A CELL
!   CC(3)    : DEATH RATE
!   CC(4)    : RATE OF DECOMPOSITION OF DETRITUS
#if defined (ZOOPL)
  real, PARAMETER :: CC(1:4)=(/0.138,1.15,4.0E-7,1.52E-7/)
#else
  real, PARAMETER :: cc(1:4)=(/0.138,1.15,1.6E-6,1.52E-7/)
#endif
#endif

!previously defined in blkdat_eco
  integer, PARAMETER :: sinkfl=1
  integer, PARAMETER :: rivflg=1
  integer, PARAMETER :: sedflg=1


!!! just for debugging
  integer, PARAMETER :: susp_i=104
  integer, PARAMETER :: susp_j=6
  integer, PARAMETER :: susp_k=14
  integer, PARAMETER :: susp_var=idia

! time step-used for the biochem. sources and sinks: = baclin*nphys_sou:
  real :: biodt_sou                           


  real pp_m2                              ! volume of box-area for primp comp.
  

  real :: bot_layer(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2,nbio  )     ! Sediment layer
  real :: bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2,nbio)         ! biochemical compartments
  real :: rad(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)         ! biochemical compartments
  real :: rivnit(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                 ! river nitrogen
  real :: rivpho(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                 ! river phosphate
  real :: rivsil(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                 ! river silicate
  real :: rivdon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                 ! river DON
  real :: grosspp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)            ! gross primary production      
  real :: netpp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)              ! net primary production      
  real :: wvel(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)               ! veritcal velocity      
!  logical pp_area(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)               ! Area for primp   
           
end module mod_necessary_ecovars
