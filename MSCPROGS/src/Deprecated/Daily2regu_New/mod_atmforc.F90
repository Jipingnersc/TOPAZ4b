module mod_atmforc
   use mod_dimensions
   type forcing_fields
      character(len=5) tforce
!
!      real :: slp(nx,ny)         !  Sea level pressure
!      real :: taux(nx,ny)        !  wind stress in x direction
!      real :: tauy(nx,ny)        !  wind stress in y direction
!      real :: wndspd(nx,ny)      !  wind speed (tke source)
!      real :: airtmp(nx,ny)      !  pseudo air temperature
!      real :: relhum(nx,ny)      !  relative humidity
!      real :: clouds(nx,ny)      !  cloud cover
!      real :: precip(nx,ny)      !  precipitation
!      real :: sss(nx,ny)         !  SSS for relax 
!      real :: sst(nx,ny)         !  SST for relax 
!      real :: uwind(nx,ny)       !  u-component of wind
!      real :: vwind(nx,ny)       !  v-component of wind
!      real :: tauxice(nx,ny)     !  ice stress on water in x direction
!      real :: tauyice(nx,ny)     !  ice stress on water in y direction
!
      real,pointer :: slp(:,:)         !  Sea level pressure
      real,pointer :: taux(:,:)        !  wind stress in x direction
      real,pointer :: tauy(:,:)        !  wind stress in y direction
      real,pointer :: wndspd(:,:)      !  wind speed (tke source)
      real,pointer :: airtmp(:,:)      !  pseudo air temperature
      real,pointer :: relhum(:,:)      !  relative humidity
      real,pointer :: clouds(:,:)      !  cloud cover
      real,pointer :: precip(:,:)      !  precipitation
      real,pointer :: sss(:,:)         !  SSS for relax 
      real,pointer :: sst(:,:)         !  SST for relax 
      real,pointer :: uwind(:,:)       !  u-component of wind
      real,pointer :: vwind(:,:)       !  v-component of wind
      real,pointer :: tauxice(:,:)     !  ice stress on water in x direction
      real,pointer :: tauyice(:,:)     !  ice stress on water in y direction

   end type forcing_fields

   type forcing_variances
      real slp
      real taux          
      real tauy
      real wndspd
      real airtmp
      real relhum
      real clouds
      real precip
      real sss
      real sst
   end type forcing_variances

!   type (forcing_fields) force
!   type (forcing_fields) ran1,ran
!   type (forcing_variances) vars

   interface operator(+)
      module procedure add_force
   end interface

   interface operator(*)
      module procedure real_force_mult,&
                       variances_forcing_fields_mult
   end interface

   interface assignment(=)
      module procedure assign_force
      module procedure assign_vars
   end interface

   interface sqrt
      module procedure var_sqrt
   end interface
   

contains
   function add_force(A,B)
      type(forcing_fields) add_force
      type(forcing_fields), intent(in) :: A
      type(forcing_fields), intent(in) :: B
      add_force%taux = A%taux + B%taux
      add_force%tauy = A%tauy + B%tauy
      add_force%wndspd = A%wndspd + B%wndspd
      add_force%relhum = A%relhum + B%relhum
      add_force%airtmp = A%airtmp + B%airtmp
      add_force%clouds = A%clouds + B%clouds
      add_force%precip = A%precip + B%precip
      add_force%sss = A%sss + B%sss
      add_force%sst = A%sst + B%sst
   end function add_force

   function real_force_mult(A,B)
      type(forcing_fields) real_force_mult
      real, intent(in) :: A
      type(forcing_fields), intent(in) :: B
      real_force_mult%taux = A * B%taux
      real_force_mult%tauy = A * B%tauy
      real_force_mult%wndspd = A * B%wndspd
      real_force_mult%airtmp = A * B%airtmp
      real_force_mult%relhum = A * B%relhum
      real_force_mult%clouds = A * B%clouds
      real_force_mult%precip = A * B%precip
      real_force_mult%sss = A * B%sss
      real_force_mult%sst = A * B%sst
   end function real_force_mult

   function variances_forcing_fields_mult(A,B)
      type(forcing_fields) variances_forcing_fields_mult
      type(forcing_variances), intent(in) :: A
      type(forcing_fields), intent(in) :: B
      variances_forcing_fields_mult%taux = A%taux * B%taux
      variances_forcing_fields_mult%tauy = A%tauy * B%tauy
      variances_forcing_fields_mult%wndspd = A%wndspd * B%wndspd
      variances_forcing_fields_mult%airtmp = A%airtmp * B%airtmp
      variances_forcing_fields_mult%relhum = A%relhum * B%relhum
      variances_forcing_fields_mult%clouds = A%clouds * B%clouds
      variances_forcing_fields_mult%precip = A%precip * B%precip
      variances_forcing_fields_mult%sss = A%sss * B%sss
      variances_forcing_fields_mult%sst = A%sst * B%sst
   end function variances_forcing_fields_mult

   function var_sqrt(A)
      type(forcing_variances) var_sqrt
      type(forcing_variances), intent(in) :: A
      var_sqrt%taux = sqrt(A%taux)
      var_sqrt%tauy = sqrt(A%tauy)
      var_sqrt%wndspd = sqrt(A%wndspd)
      var_sqrt%airtmp = sqrt(A%airtmp)
      var_sqrt%relhum = sqrt(A%relhum)
      var_sqrt%clouds = sqrt(A%clouds)
      var_sqrt%precip = sqrt(A%precip)
      var_sqrt%sss = sqrt(A%sss)
      var_sqrt%sst = sqrt(A%sst)
   end function var_sqrt

   subroutine assign_force(A,r)
      type(forcing_fields), intent(out) :: A
      real, intent(in) :: r
      A%slp  = r
      A%taux = r
      A%tauy = r
      A%wndspd = r
      A%airtmp = r
      A%relhum = r
      A%clouds = r
      A%precip = r
      A%sss = r
      A%sst = r
      A%uwind = r
      A%vwind = r
      A%tauxice = r
      A%tauyice = r
   end subroutine assign_force


   subroutine assign_vars(A,r)
      type(forcing_variances), intent(out) :: A
      real, intent(in) :: r
      A%taux = r
      A%tauy = r
      A%wndspd = r
      A%airtmp = r
      A%relhum = r
      A%clouds = r
      A%precip = r
      A%sss = r
      A%sst = r
   end subroutine assign_vars

end module  mod_atmforc
