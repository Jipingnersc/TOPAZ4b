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
      character(len=80) :: standard_name    ! netcdf standard_name (COARDS conventions)
      character(len=80) :: standard_name_rot! netcdf standard_name for rotated vel(COARDS conventions)
      character(len=80) :: long_name        ! A descriptive name
      character(len=80) :: long_name_rot    ! A descriptive name
      character(len=20) :: units            ! units (COARDS)
      character(len=40) :: cell_methods     ! Info on averaging etc 
      logical           :: process          ! Should we process this variable?
      logical           :: gridrotate       ! Rotate velocities to output grid?
      logical           :: isvectorcomp     ! Is this a vector component?
      real              :: max_val          ! max allowed value
      real              :: min_val          ! min allowed value
      real              :: scale_factor     ! calculated from max/min
      real              :: add_offset       ! calculated from max/min
      integer*2         :: fillvalue        ! Value for no data
      real              :: num_dim          ! Number of dims of variable (2/3)
   end type netcdf_pars

   integer, parameter :: maxvar_info=50
   integer :: nvar_info

   type (netcdf_pars) :: var_info(maxvar_info)

contains

subroutine init_var_info(hycfile,gridrotate)
   implicit none
   integer, intent(in) :: hycfile
   logical, intent(in) :: gridrotate
   integer :: irec,i,ios,nbits
   logical :: ex
   character*8 :: char8


   var_info(:)%cell_methods=''
   if (hycfile==2 .or. hycfile==3) then
      var_info(:)%cell_methods='lat: lon: time: mean'
   elseif (hycfile==1) then
      var_info(:)%cell_methods='lat: lon: mean  time: point'
   end if

   ! Default
   var_info(:)%standard_name_rot=''
   var_info(:)%long_name_rot    =''
   var_info(:)%gridrotate=.false.
   var_info(:)%isvectorcomp=.false.

   irec=1

   ! Define fields - temperature
   var_info(irec)%name8        ='temp    '                        
   var_info(irec)%vname        ='temperature'
   var_info(irec)%standard_name='sea_water_potential_temperature'
   var_info(irec)%long_name    ='Potential temperature of sea water'
   var_info(irec)%units        ='degrees_celsius'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =45.
   var_info(irec)%min_val      =-3.0
   var_info(irec)%num_dim      =3
   irec=irec+1

   ! Define fields - salinity
   var_info(irec)%name8        ='saln    '
   var_info(irec)%vname        ='salinity'
   var_info(irec)%standard_name='sea_water_salinity'
   var_info(irec)%long_name    ='Salinity of sea water'
   var_info(irec)%units        ='psu'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =50.
   var_info(irec)%min_val      =0.
   var_info(irec)%num_dim      =3
   irec=irec+1


   ! Define fields - velocity
   var_info(irec)%name8        ='utot    '
   var_info(irec)%vname        ='u'
   var_info(irec)%standard_name='eastward_sea_water_velocity'
   var_info(irec)%standard_name_rot='x_sea_water_velocity'
   var_info(irec)%long_name    ='Eastward velocity of sea water'
   var_info(irec)%long_name_rot='X velocity component of sea water'
   var_info(irec)%units        ='m s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =10.
   var_info(irec)%min_val      =-10.
   var_info(irec)%num_dim      =3
   var_info(irec)%isvectorcomp =.true.
   var_info(irec)%gridrotate   =gridrotate
   irec=irec+1

   ! Define fields - velocity
   var_info(irec)%name8        ='vtot    '
   var_info(irec)%vname        ='v'
   var_info(irec)%standard_name='northward_sea_water_velocity'
   var_info(irec)%standard_name_rot ='x_sea_water_velocity'
   var_info(irec)%long_name    ='Northward velocity of sea water'
   var_info(irec)%long_name_rot='Y velocity component of sea water'
   var_info(irec)%units        ='m s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =10.
   var_info(irec)%min_val      =-10.
   var_info(irec)%num_dim      =3
   var_info(irec)%isvectorcomp =.true.
   var_info(irec)%gridrotate   =gridrotate
   irec=irec+1


   ! Define fields - ssh
   var_info(irec)%name8        ='ssh     '
   var_info(irec)%vname        ='ssh'
   var_info(irec)%standard_name='sea_surface_elevation'
   var_info(irec)%long_name    ='Sea Surface height'
   var_info(irec)%units        ='m'
   if (hycfile==1) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   end if
   var_info(irec)%max_val      =10.
   var_info(irec)%min_val      =-10.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - mld1
   var_info(irec)%name8        ='mld1    ' !Placeholder
   var_info(irec)%vname        ='mld'
   var_info(irec)%standard_name='ocean_mixed_layer_thickness'
   var_info(irec)%long_name    ='Temperature Defined Mixed Layer Depth'
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
   var_info(irec)%long_name    ='Density Defined Mixed Layer Depth'
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
   var_info(irec)%standard_name_rot='surface_downward_x_stress'
   var_info(irec)%long_name    ='Area averaged eastward atmospheric stress'
   var_info(irec)%long_name_rot='Area averaged x-component of atmospheric stress'
   var_info(irec)%units        ='pascal'
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   end if
   var_info(irec)%max_val      =2.
   var_info(irec)%min_val      =-2.
   var_info(irec)%num_dim      =2
   var_info(irec)%isvectorcomp =.true.
   var_info(irec)%gridrotate   =gridrotate
   irec=irec+1

   ! Define fields - tauxi
   var_info(irec)%name8        ='tauxi   ' 
   var_info(irec)%vname        ='tauxice'
   var_info(irec)%standard_name='downward_eastward_stress_at_sea_ice_base'
   var_info(irec)%standard_name_rot='downward_x_stress_at_sea_ice_base'
   var_info(irec)%long_name    ='Area averaged eastward ice stress on ocean'
   var_info(irec)%long_name_rot='Area averaged x-component of ice stress on ocean'
   var_info(irec)%units        ='pascal'
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   endif
   var_info(irec)%max_val      =2.
   var_info(irec)%min_val      =-2.
   var_info(irec)%num_dim      =2
   var_info(irec)%isvectorcomp =.true.
   var_info(irec)%gridrotate   =gridrotate
   irec=irec+1

   ! Define fields - tauy
   var_info(irec)%name8        ='tauy    ' 
   var_info(irec)%vname        ='tauy'
   var_info(irec)%standard_name='surface_downward_northward_stress'
   var_info(irec)%standard_name_rot='surface_downward_y_stress'
   var_info(irec)%long_name    ='Area averaged northward atmospheric stress'
   var_info(irec)%long_name_rot='Area averaged y-component of atmospheric stress'
   var_info(irec)%units        ='pascal'
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   endif
   var_info(irec)%max_val      =2.
   var_info(irec)%min_val      =-2.
   var_info(irec)%num_dim      =2
   var_info(irec)%isvectorcomp =.true.
   var_info(irec)%gridrotate   =gridrotate
   irec=irec+1


   ! Define fields - tauyi
   var_info(irec)%name8        ='tauyi   ' 
   var_info(irec)%vname        ='tauyice'
   var_info(irec)%standard_name='downward_northward_stress_at_sea_ice_base '
   var_info(irec)%standard_name_rot='downward_y_stress_at_sea_ice_base '
   var_info(irec)%long_name    ='Area averaged northwardward ice stress on ocean'
   var_info(irec)%long_name_rot='Area averaged y-component of ice stress on ocean'
   var_info(irec)%units        ='pascal'
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   endif
   var_info(irec)%max_val      =2.
   var_info(irec)%min_val      =-2.
   var_info(irec)%num_dim      =2
   var_info(irec)%isvectorcomp =.true.
   var_info(irec)%gridrotate   =gridrotate
   irec=irec+1


   ! Define fields - OCEAN heat flux
   var_info(irec)%name8        ='surflx  ' 
   var_info(irec)%vname        ='qtot'
   var_info(irec)%standard_name='surface_downward_heat_flux_in_sea_water'
   var_info(irec)%long_name    ='Net heat flux into ocean'
   var_info(irec)%units        ='w m-2'
   if (hycfile==1) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   endif
   var_info(irec)%max_val      =3000.
   var_info(irec)%min_val      =-3000.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - SURFACE heat flux (modified by ice)
   var_info(irec)%name8        ='surflxs ' 
   var_info(irec)%vname        ='qtotair'
   var_info(irec)%standard_name='surface_downward_heat_flux_in_air'
   var_info(irec)%long_name    ='Net surface atmospheric heat flux'
   var_info(irec)%units        ='w m-2'
   if (hycfile==1) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   endif
   var_info(irec)%max_val      =3000.
   var_info(irec)%min_val      =-3000.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - swflxs
   var_info(irec)%name8        ='swflx   ' 
   var_info(irec)%vname        ='qsw'
   var_info(irec)%standard_name='net_downward_shortwave_flux_in_air'
   var_info(irec)%long_name    ='Net surface shortwave heat flux'
   var_info(irec)%units        ='w m-2'
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   endif
   var_info(irec)%max_val      =2000.
   var_info(irec)%min_val      =0.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - surface albedo
   var_info(irec)%name8        ='albedo  ' 
   var_info(irec)%vname        ='albedo'
   var_info(irec)%standard_name='surface_albedo'
   var_info(irec)%long_name    ='Albedo of surface'
   var_info(irec)%units        =''
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   endif
   var_info(irec)%max_val      =1.01
   var_info(irec)%min_val      =0.
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - emnp
   var_info(irec)%name8        ='emnp    ' 
   var_info(irec)%vname        ='fwflux'
   var_info(irec)%standard_name='water_flux_into_ocean'
   var_info(irec)%long_name    ='Freshwater flux into ocean'
   var_info(irec)%units        ='kg m-2 s-1'
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   endif
   var_info(irec)%max_val      = 1e-3
   var_info(irec)%min_val      =-1e-3
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - ub
   var_info(irec)%name8        ='ubavg   ' 
   var_info(irec)%vname        ='ub'
   var_info(irec)%standard_name='eastward_barotropic_sea_water_velocity'
   var_info(irec)%standard_name_rot='x_barotropic_sea_water_velocity'
   var_info(irec)%long_name    ='Eastward barotropic ocean velocity'
   var_info(irec)%long_name_rot='X-component of barotropic ocean velocity'
   var_info(irec)%units        ='m s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =  5.
   var_info(irec)%min_val      = -5.
   var_info(irec)%num_dim      =2
   var_info(irec)%isvectorcomp =.true.
   var_info(irec)%gridrotate   =gridrotate
   irec=irec+1

   ! Define fields - vb
   var_info(irec)%name8        ='vbavg   ' 
   var_info(irec)%vname        ='vb'
   var_info(irec)%standard_name='northward_barotropic_sea_water_velocity'
   var_info(irec)%standard_name_rot='y_barotropic_sea_water_velocity'
   var_info(irec)%long_name    ='Northward barotropic ocean velocity'
   var_info(irec)%long_name_rot='Y-component of barotropic ocean velocity'
   var_info(irec)%units        ='m s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      =  5.
   var_info(irec)%min_val      = -5.
   var_info(irec)%num_dim      =2
   var_info(irec)%gridrotate   =gridrotate
   irec=irec+1

   ! Define fields - stream function
   var_info(irec)%name8        ='bsfd   '  ! Placeholder
   var_info(irec)%vname        ='bsfd'
   var_info(irec)%standard_name='ocean_barotropic_streamfunction'        
   var_info(irec)%long_name    ='Barotropic Streamfunction of the ocean'     
   var_info(irec)%units        ='m3 s-1'
   var_info(irec)%process      =.true.
   var_info(irec)%max_val      = 1e10
   var_info(irec)%min_val      =-1e10
   var_info(irec)%num_dim      =2
   var_info(irec)%isvectorcomp =.true.
   irec=irec+1

   ! Define fields - ice concentration
   var_info(irec)%name8        ='fice   '
   var_info(irec)%vname        ='fice'
   var_info(irec)%standard_name='sea_ice_area_fraction'        
   var_info(irec)%long_name    ='Sea Ice Concentration'     
   var_info(irec)%units        =''
   var_info(irec)%max_val      =1.01
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   end if
   var_info(irec)%min_val      =-0.01
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - ice thickness
   var_info(irec)%name8        ='hice   '
   var_info(irec)%vname        ='hice'
   var_info(irec)%standard_name='sea_ice_thickness'        
   var_info(irec)%long_name    ='Sea Ice Thickness (area average)'     
   var_info(irec)%units        ='m'
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   end if
   var_info(irec)%max_val      =15.0
   var_info(irec)%min_val      =-0.01
   var_info(irec)%num_dim      =2
   irec=irec+1

   ! Define fields - snow thickness
   var_info(irec)%name8        ='hsnw   '
   var_info(irec)%vname        ='hsnow'
   var_info(irec)%standard_name='surface_snow_thickness'        
   var_info(irec)%long_name    ='Thickness of snow on top of sea ice (area average)'
   var_info(irec)%units        ='m'
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   end if
   var_info(irec)%max_val      =15.0
   var_info(irec)%min_val      =-0.01
   var_info(irec)%num_dim      =2
   irec=irec+1


   ! Define fields - ice velocity
   var_info(irec)%name8        ='uice   '
   var_info(irec)%vname        ='uice'
   var_info(irec)%standard_name='eastward_sea_ice_velocity'        
   var_info(irec)%standard_name_rot='x_sea_ice_velocity'        
   var_info(irec)%long_name    ='Eastward velocity of sea ice (area average)'     
   var_info(irec)%long_name_rot='X velocity component of sea ice (area average)'     
   var_info(irec)%units        ='m s-1'
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   end if
   var_info(irec)%max_val      =2.
   var_info(irec)%min_val      =-2.
   var_info(irec)%num_dim      =2
   var_info(irec)%isvectorcomp =.true.
   var_info(irec)%gridrotate   =gridrotate
   irec=irec+1

   ! Define fields - icevelocity
   var_info(irec)%name8        ='vice   '
   var_info(irec)%vname        ='vice'
   var_info(irec)%standard_name='northward_sea_ice_velocity'        
   var_info(irec)%standard_name_rot='y_sea_ice_velocity'        
   var_info(irec)%long_name    ='Northward velocity of sea ice (area average)'     
   var_info(irec)%long_name_rot='y velocity component of sea ice (area average)'     
   var_info(irec)%units        ='m s-1'
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   end if
   var_info(irec)%max_val      =2.
   var_info(irec)%min_val      =-2.
   var_info(irec)%num_dim      =2
   var_info(irec)%isvectorcomp =.true.
   var_info(irec)%gridrotate   =gridrotate
   irec=irec+1

   ! Define fields - ice thickness tendency
   var_info(irec)%name8        ='htndncy'
   var_info(irec)%vname        ='htendency'
   var_info(irec)%standard_name='tendency_of_sea_ice_thickness_due_to_thermodynamics'
   var_info(irec)%long_name    ='Ice thickness tendency due to thermodynamics'
   var_info(irec)%units        ='m s-1'
   if (hycfile==1 .or. hycfile==3) then
      var_info(irec)%process      =.false.
   else
      var_info(irec)%process      =.true.
   end if
   var_info(irec)%max_val      =.00001
   var_info(irec)%min_val      =-.00001
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
   inquire(exist=ex,file='fields.hyc2proj')
   if (ex) then
      do i=1,nvar_info
         var_info(i)%process=.false.
      end do
      open(10,file='fields.hyc2proj', status='old')
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

      print *,'"fields.hyc2proj" found - will dump these vars (if present)'
      do i=1,nvar_info
         if (var_info(i)%process) then
            print *,var_info(i)%name8
         end if
      end do
   else
      print *,'No file "fields.hyc2proj" found - will dump all vars'
   end if

         
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   logical function var_gridrotate(varinfo)
   implicit none
   type (netcdf_pars), intent(in)  :: varinfo
   var_gridrotate=varinfo%gridrotate
   end function
   

   function var_longname(varinfo)
   implicit none
   type (netcdf_pars), intent(in)  :: varinfo
   character(len=80) :: var_longname
   if (varinfo%gridrotate) then
      var_longname=varinfo%long_name_rot
   else
      var_longname=varinfo%long_name
   end if
   end function

   function var_standard_name(varinfo)
   implicit none
   type (netcdf_pars), intent(in)  :: varinfo
   character(len=80) :: var_standard_name
   if (varinfo%gridrotate) then
      var_standard_name=varinfo%standard_name_rot
   else
      var_standard_name=varinfo%standard_name
   end if
   end function
   
   


end module
