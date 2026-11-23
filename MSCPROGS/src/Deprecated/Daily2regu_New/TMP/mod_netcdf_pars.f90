

module mod_netcdf_pars

   ! Module containing variable info for netcdf dump. To add entry, put it 
   ! use the procedure below (See temperature for instance), remember to
   ! increase the "irec" counter....

   ! Only the name8, vname, and max/min val are important to specify. If you are
   ! unsure of what to put in other variable declarations, just put blanks. 
   
   ! COARDS are some definitions on how a netcdf file shuld look like. If
   ! you want to be a nice guy, look up the COARDS definitions for how to
   ! specify units. Also important for COARDS is the "standard_name" of your 
   ! variable. The standard_name makes it possible for software to know what
   ! variable it is looking at.

   ! KAL - 23.05.2005


   type netcdf_pars
      character(len=8 ) :: name8            ! name used in daily ave file
      character(len=20) :: vname            ! netcdf variable name
      character(len=20) :: standard_name    ! netcdf standard_name (COARDS conventions)
      character(len=20) :: long_name        ! A descriptive name
      character(len=20) :: units            ! units (COARDS)
      logical           :: process          ! Should we process this variable?
      real              :: max_val          ! max allowed value
      real              :: min_val          ! min allowed value
      real              :: scale_factor     ! calculated from max/min
      real              :: add_offset       ! calculated from max/min
      real              :: fillvalue        ! Value for no data
      real              :: num_dim          ! Number of dims of variable (2/3)
   end type netcdf_pars

   integer, parameter :: maxvar_info=50
   integer :: nvar_info

   type (netcdf_pars) :: var_info(maxvar_info)

contains

subroutine init_var_info()
   implicit none
   integer :: irec,i,ios,nbits
   logical :: ex
   character*8 :: char8

   irec=1

   ! Define fields - temperature
   var_info(irec)%name8        ='temp    '
   var_info(irec)%vname        ='temperature'
   var_info(irec)%standard_name='sea_water_potential_temperature'
   var_info(irec)%long_name    ='potential temperature'
   var_info(irec)%units        ='degrees_celsius'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =45.
   var_info(irec)%min_val      =-2.0
   var_info(irec)%num_dim      =3
   irec=irec+1

   ! Define fields - salinity
   var_info(irec)%name8        ='saln    '
   var_info(irec)%vname        ='salinity'
   var_info(irec)%standard_name='sea_water_salinity'
   var_info(irec)%long_name    ='salinity'
   var_info(irec)%units        ='psu'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =45.
   var_info(irec)%min_val      =10.
   var_info(irec)%num_dim      =3
   irec=irec+1


   ! Define fields - velocity
   var_info(irec)%name8        ='utot    '
   var_info(irec)%vname        ='u'
   var_info(irec)%standard_name='eastward_sea_water_velocity'
   var_info(irec)%long_name    ='Eastward Velocity'
   var_info(irec)%units        ='m s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =5.
   var_info(irec)%min_val      =-5.
   var_info(irec)%num_dim      =3
   irec=irec+1

   ! Define fields - velocity
   var_info(irec)%name8        ='vtot    '
   var_info(irec)%vname        ='v'
   var_info(irec)%standard_name='northward_sea_water_velocity'
   var_info(irec)%long_name    ='Northward Velocity'
   var_info(irec)%units        ='m s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =5.
   var_info(irec)%min_val      =-5.
   var_info(irec)%num_dim      =3
   irec=irec+1


   ! Define fields - velocity
   var_info(irec)%name8        ='ssh     '
   var_info(irec)%vname        ='ssh'
   var_info(irec)%standard_name='sea_surface_elevation'
   var_info(irec)%long_name    ='Sea Surface height'
   var_info(irec)%units        ='m'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =3.
   var_info(irec)%min_val      =-3.
   var_info(irec)%num_dim      =32
   irec=irec+1

   ! Define fields - mld1
   var_info(irec)%name8        ='mld1    ' !Placeholder
   var_info(irec)%vname        ='mld'
   var_info(irec)%standard_name='ocean_mixed_layer_thickness'
   var_info(irec)%long_name    ='Temperature Defined Mixed Layer'
   var_info(irec)%units        ='m'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =6000.
   var_info(irec)%min_val      =0.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - mld2
   var_info(irec)%name8        ='mld2    ' !Placeholder
   var_info(irec)%vname        ='mlp'
   var_info(irec)%standard_name='ocean_mixed_layer_thickness'
   var_info(irec)%long_name    ='Density Defined Mixed Layer'
   var_info(irec)%units        ='m'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =6000.
   var_info(irec)%min_val      =0.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - taux
   var_info(irec)%name8        ='taux    ' !Placeholder
   var_info(irec)%vname        ='taux'
   var_info(irec)%standard_name='surface_downward_eastward_stress'
   var_info(irec)%long_name    ='Atmospheric stress eastwards'
   var_info(irec)%units        ='pascal'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =2.
   var_info(irec)%min_val      =-2.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - tauy
   var_info(irec)%name8        ='tauy    ' 
   var_info(irec)%vname        ='tauy'
   var_info(irec)%standard_name='surface_downward_northward_stress'
   var_info(irec)%long_name    ='Atmospheric stress northwards'
   var_info(irec)%units        ='pascal'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =2.
   var_info(irec)%min_val      =-2.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - tauy
   var_info(irec)%name8        ='surflx  ' 
   var_info(irec)%vname        ='qtot'
   var_info(irec)%standard_name='surface_downward_heat_flux_in_air'
   var_info(irec)%long_name    ='Net surface heat flux into ocean'
   var_info(irec)%units        ='w m-2'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =3000.
   var_info(irec)%min_val      =-3000.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - swflx
   var_info(irec)%name8        ='swflx   ' 
   var_info(irec)%vname        ='qsw'
   var_info(irec)%standard_name='net_downward_shortwave_flux_in_air'
   var_info(irec)%long_name    ='shortwave heat flux into ocean'
   var_info(irec)%units        ='w m-2'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =3000.
   var_info(irec)%min_val      =0.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - emnp
   var_info(irec)%name8        ='emnp    ' 
   var_info(irec)%vname        ='emp'
   var_info(irec)%standard_name='water_flux_into_ocean'
   var_info(irec)%long_name    ='Freshwater flux'
   var_info(irec)%units        ='kg m-2 s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      = 1e-3
   var_info(irec)%min_val      =-1e-3
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - ub
   var_info(irec)%name8        ='ubavg   ' 
   var_info(irec)%vname        ='ub'
   var_info(irec)%standard_name='eastward_barotropic_sea_water_velocity'
   var_info(irec)%long_name    ='Eastward barotropic velocity'
   var_info(irec)%units        ='m s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =  5.
   var_info(irec)%min_val      = -5.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - vb
   var_info(irec)%name8        ='vbavg   ' 
   var_info(irec)%vname        ='vb'
   var_info(irec)%standard_name='northward_barotropic_sea_water_velocity'
   var_info(irec)%long_name    ='Northward barotropic velocity'
   var_info(irec)%units        ='m s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =  5.
   var_info(irec)%min_val      = -5.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - stream function
   var_info(irec)%name8        ='bsfd   '  ! Placeholder
   var_info(irec)%vname        ='bsfd'
   var_info(irec)%standard_name='ocean_barotropic_streamfunction'        
   var_info(irec)%long_name    ='Barotropic Streamfunction'     
   var_info(irec)%units        ='m3 s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      = 1e10
   var_info(irec)%min_val      =-1e10
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - ice concentration
   var_info(irec)%name8        ='fice   '
   var_info(irec)%vname        ='fice'
   var_info(irec)%standard_name='sea_ice_concentration'        
   var_info(irec)%long_name    ='Sea Ice Concentration'     
   var_info(irec)%units        =''
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =1.00
   var_info(irec)%min_val      =0.00
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - ice thickness
   var_info(irec)%name8        ='hice   '
   var_info(irec)%vname        ='hice'
   var_info(irec)%standard_name='sea_ice_thickness'        
   var_info(irec)%long_name    ='Sea Ice Thickness'     
   var_info(irec)%units        ='m'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =15.0
   var_info(irec)%min_val      =0.0
   var_info(irec)%num_dim      =2
   irec=irec+1


   ! Define fields - ice velocity
   var_info(irec)%name8        ='uice   '
   var_info(irec)%vname        ='uice'
   var_info(irec)%standard_name='eastward_sea_ice_velocity'        
   var_info(irec)%long_name    ='Eastward ice Velocity'     
   var_info(irec)%units        ='m s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =2.
   var_info(irec)%min_val      =-2.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - icevelocity
   var_info(irec)%name8        ='vice   '
   var_info(irec)%vname        ='vice'
   var_info(irec)%standard_name='northwardward_sea_ice_velocity'        
   var_info(irec)%long_name    ='Northward ice Velocity'     
   var_info(irec)%units        ='m s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =2.
   var_info(irec)%min_val      =-2.
   var_info(irec)%num_dim      =2
   irec=irec+1
   
   

   nvar_info=irec-1

   nbits=2**16

   do i=1,nvar_info
      ! Convert to 16 bits
      var_info(i)%add_offset   = (var_info(i)%max_val+var_info(i)%min_val)/2.
      var_info(i)%scale_factor = (var_info(i)%max_val-var_info(i)%min_val)/(nbits-10)
      var_info(i)%fillvalue=-(nbits-2)/2
   end do

   ! Check for file "fields.daily2regunc". If it exists, use the fields
   ! specified there to determine what should be processed
   inquire(exist=ex,file='fields.daily2regunc')
   if (ex) then
      do i=1,nvar_info
         var_info(i)%process=.false.
      end do
      open(10,file='fields.daily2regunc', status='old')
      ios=0
      do while(ios==0)
         read(10,'(a8)',iostat=ios) char8
         if (ios==0) then
            do i=1,nvar_info
               var_info(i)%process=var_info(i)%process.or. &
                  trim(char8)==trim(var_info(i)%name8)
            end do
         end if
      end do
      close(10)

      print *,'"fields.daily2regunc" found - will dump these vars (if present)'
      do i=1,nvar_info
         if (var_info(i)%process) then
            print *,var_info(i)%name8
         end if
      end do
   else
      print *,'No file "fields.daily2regunc" found - will dump all vars'
   end if

         
   end subroutine

   integer function var_info_getindex(inname)
   implicit none
   character(len=*), intent(in) :: inname
   integer :: i

   var_info_getindex=-1
   do i=1,nvar_info
      if (trim(var_info(i)%name8)==trim(inname)) var_info_getindex=i
   end do

   if (var_info_getindex==-1) then
      print *,'Warning: Could not find field '//trim(inname)
   end if
   end function



end module
