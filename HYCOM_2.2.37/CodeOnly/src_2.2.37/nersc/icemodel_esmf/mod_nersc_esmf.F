      module mod_esmf_nersc
#if ! defined(USE_ESMF_NERSC)
      logical, parameter :: empty
#else
      use ESMF_Mod    ! ESMF  Framework
      use mod_xc      ! HYCOM communication interface
      use mod_za      ! HYCOM I/O interface
c
c --- -----------------------------------------
c --- MICOM-based HYbrid Coordinate Ocean Model
c ---               H Y C O M
c ---           v e r s i o n  2.2
c ---          ESMF  interface only 
c ---           to NERSC submodels
c --- -----------------------------------------
c
      implicit none
c
c
      public HYCOM_SetServices
c
      logical, save, public  :: put_export   !set in main program
      logical, save, public  :: get_import   !set in main program
      logical, save, public  :: end_of_run   !set in HYCOM_Run
      integer, save, public  :: nts_day      !set in HYCOM_init, timesteps/day
      integer, save, public  :: nts_ice      !set in HYCOM_init, timesteps/ice
c
c
c
c --- Data types for Import/Export array pointers
      type ArrayPtrReal2D
        real(ESMF_KIND_R4), dimension(:,:), pointer :: p
      end type ArrayPtrReal2D
c
c --- Attribute names for fields
      character(ESMF_MAXSTR), save ::
     &    attNameLongName = "long_name",
     &    attNameStdName  = "standard_name",
     &    attNameUnits    = "units",
     &    attNameSclFac   = "scale_factor",
     &    attNameAddOff   = "add_offset"
c
c --- Import Fields
      integer, parameter :: numImpFields=11
      character(ESMF_MAXSTR), save :: impFieldName(    numImpFields),
     &                                impFieldLongName(numImpFields),
     &                                impFieldStdName( numImpFields),
     &                                impFieldUnits(   numImpFields)
      real(ESMF_KIND_R4),     save :: impFieldSclFac(  numImpFields),
     &                                impFieldAddOff(  numImpFields)
c
c --- Export Fields
      integer, parameter :: numExpFields=7
      character(ESMF_MAXSTR), save :: expFieldName(    numExpFields),
     &                                expFieldLongName(numExpFields),
     &                                expFieldStdName( numExpFields),
     &                                expFieldUnits(   numExpFields)
      real(ESMF_KIND_R4),     save :: expFieldSclFac(  numExpFields),
     &                                expFieldAddOff(  numExpFields)
c
c --- ESMF related variables
      type(ESMF_Bundle),    save :: expBundle,
     &                              impBundle
      type(ESMF_Field),     save :: expField(numExpFields),
     &                              impField(numImpFields)
      type(ArrayPtrReal2D), save :: expData( numExpFields),
     &                              impData( numImpFields)
c
      type(ESMF_Clock),     save :: intClock
      type(ESMF_VM),        save :: vm
      type(ESMF_DELayout),  save :: deLayout
      integer,              save :: petCount, localPet, mpiCommunicator
      type(ESMF_Grid),      save :: grid2D
      type(ESMF_ArraySpec), save :: arraySpec2Dr

      real, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &       sic_import !Sea Ice Concentration
     &,     sitx_import !Sea Ice X-Stress
     &,     sity_import !Sea Ice Y-Stress
     &,     siqs_import !Solar Heat Flux thru Ice to Ocean
     &,     sifh_import !Ice Freezing/Melting Heat Flux
     &,     sifs_import !Ice Freezing/Melting Salt Flux
     &,     sifw_import !Ice Net Water Flux
     &,      sit_import !Sea Ice Temperature
     &,      sih_import !Sea Ice Thickness
     &,      siu_import !Sea Ice X-Velocity
     &,      siv_import !Sea Ice Y-Velocity
     &,      ocn_mask   !Ocean Currents Mask

      logical, save ::
     &       ocn_mask_init

#endif
      contains

      subroutine HYCOM_SetServices(gridComp, rc)
c
      type(ESMF_GridComp) :: gridComp
      integer             :: rc
c
      call ESMF_GridCompSetEntryPoint(
     &     gridComp,
     &     ESMF_SETINIT,
     &     HYCOM_Init,
     &     ESMF_SINGLEPHASE,
     &     rc)
      call ESMF_GridCompSetEntryPoint(
     &     gridComp,
     &     ESMF_SETRUN,
     &     HYCOM_Run,
     &     ESMF_SINGLEPHASE,
     &     rc)
      call ESMF_GridCompSetEntryPoint(
     &     gridComp,
     &     ESMF_SETFINAL,
     &     HYCOM_Final,
     &     ESMF_SINGLEPHASE,
     &     rc)
c
      end subroutine HYCOM_SetServices

      subroutine Setup_ESMF(gridComp, impState, expState, extClock, rc)
c
c --- Calling parameters
      type(ESMF_GridComp) :: gridComp
      type(ESMF_State)    :: impState
      type(ESMF_State)    :: expState
      type(ESMF_Clock)    :: extClock
      integer             :: rc
c
c --- set up ESMF data structures for HYCOM.
c
      integer                 :: i,j
      real(ESMF_KIND_R8)      :: coord1(itdm),coord2(jtdm) !cell centers
      real                    :: plonrl(itdm),platrl(jtdm)
      integer                 :: cnts(2)
      real(ESMF_KIND_R8)      :: mgcpd(2), dpd(2)
      character(10)           :: dimNames(2),dimUnits(2)
      type(ESMF_Logical)      :: periodic(2)
      integer(ESMF_KIND_I4)   :: year,month,day,hour,minute
      integer(ESMF_KIND_I4)   :: sec,msec,usec,nsec
      real(8)                 :: dsec,dmsec,dusec,dnsec
      type(ESMF_TimeInterval) :: timeStep, runDuration
      type(ESMF_Time)         :: startTime
      character(ESMF_MAXSTR)  :: msg, gridName
c
c --- Report
      call ESMF_LogWrite("HYCOM Setup routine called",
     &     ESMF_LOG_INFO)
!-----call ESMF_LogFlush
c
c  Attributes for import fields, identical to CICE export fields
      impFieldAddOff(:) = 0.0 !default is no offset
      impFieldSclFac(:) = 1.0 !default is no scale factor

      impFieldName(     1) = "sic"
      impFieldLongName( 1) = "Sea Ice Concentration"
      impFieldStdName(  1) = "sea_ice_area_fraction"
      impFieldUnits(    1) = "1"
      impFieldName(     2) = "sitx"
      impFieldLongName( 2) = "Sea Ice X-Stress"
      impFieldStdName(  2) = "downward_x_stress_at_sea_ice_base"
      impFieldSclFac(   2) = -1.0  !field is upward
      impFieldUnits(    2) = "Pa"
      impFieldName(     3) = "sity"
      impFieldLongName( 3) = "Sea Ice Y-Stress"
      impFieldStdName(  3) = "downward_y_stress_at_sea_ice_base"
      impFieldSclFac(   3) = -1.0  !field is upward
      impFieldUnits(    3) = "Pa"
      impFieldName(     4) = "siqs"
      impFieldLongName( 4) = "Solar Heat Flux thru Ice to Ocean"
      impFieldStdName(  4) = "downward_sea_ice_basal_solar_heat_flux"
      impFieldUnits(    4) = "W m-2"
      impFieldName(     5) = "sifh"
      impFieldLongName( 5) = "Ice Freezing/Melting Heat Flux"
      impFieldStdName(  5) = "upward_sea_ice_basal_heat_flux"
      impFieldSclFac(   5) = -1.0  !field is downward
      impFieldUnits(    5) = "W m-2"
      impFieldName(     6) = "sifs"
      impFieldLongName( 6) = "Ice Freezing/Melting Salt Flux"
      impFieldStdName(  6) = "downward_sea_ice_basal_salt_flux"
      impFieldUnits(    6) = "kg m-2 s-1"
      impFieldName(     7) = "sifw"
      impFieldLongName( 7) = "Ice Net Water Flux"
      impFieldStdName(  7) = "downward_sea_ice_basal_water_flux"
      impFieldUnits(    7) = "kg m-2 s-1"
      impFieldName(     8) = "sit" !diagnostic
      impFieldLongName( 8) = "Sea Ice Temperature"
      impFieldStdName(  8) = "sea_ice_temperature"
      impFieldAddOff(   8) = +273.15 !field is in degC
      impFieldUnits(    8) = "K"
      impFieldName(     9) = "sih" !diagnostic
      impFieldLongName( 9) = "Sea Ice Thickness"
      impFieldStdName(  9) = "sea_ice_thickness"
      impFieldUnits(    9) = "m"
      impFieldName(    10) = "siu" !diagnostic
      impFieldLongName(10) = "Sea Ice X-Velocity"
      impFieldStdName( 10) = "sea_ice_x_velocity"
      impFieldUnits(   10) = "m s-1"
      impFieldName(    11) = "siv" !diagnostic
      impFieldLongName(11) = "Sea Ice Y-Velocity"
      impFieldStdName( 11) = "sea_ice_y_velocity"
      impFieldUnits(   11) = "m s-1"
*
*     impFieldName(    12) = "patm"
*     impFieldLongName(12) = "Surface Air Pressure"
*     impFieldStdName( 12) = "surface_air_pressure"
*     impFieldUnits(   12) = "Pa"
*     impFieldName(    13) = "xwnd"
*     impFieldLongName(13) = "X-Wind"
*     impFieldStdName( 13) = "x_wind"
*     impFieldUnits(   13) = "m s-1"
*     impFieldName(    14) = "ywnd"
*     impFieldLongName(14) = "Y-Wind"
*     impFieldStdName( 14) = "y_wind"
*     impFieldUnits(   14) = "m s-1"
c
c  Attributes for export fields, identical to CICE import fields
      expFieldAddOff(:) = 0.0 !default is no offset
      expFieldSclFac(:) = 1.0 !default is no scale factor

      expFieldName(     1) = "sst"
      expFieldLongName( 1) = "Sea Surface Temperature"
      expFieldStdName(  1) = "sea_surface_temperature"
      expFieldAddOff(   1) = +273.15 !field is in degC
      expFieldUnits(    1) = "K"
      expFieldName(     2) = "sss"
      expFieldLongName( 2) = "Sea Surface Salinity"
      expFieldStdName(  2) = "sea_surface_salinity"
      expFieldUnits(    2) = "1e-3"
      expFieldName(     3) = "ssu"
      expFieldLongName( 3) = "Sea Surface X-Current"
      expFieldStdName(  3) = "sea_water_x_velocity"
      expFieldUnits(    3) = "m s-1"
      expFieldName(     4) = "ssv"
      expFieldLongName( 4) = "Sea Surface Y-Current"
      expFieldStdName(  4) = "sea_water_y_velocity"
      expFieldUnits(    4) = "m s-1"
      expFieldName(     5) = "ssh"
      expFieldLongName( 5) = "Sea Surface Height"
      expFieldStdName(  5) = "sea_surface_height_above_sea_level"
      expFieldUnits(    5) = "m"
      expFieldName(     6) = "ssfi"
      expFieldLongName( 6) = "Oceanic Heat Flux Available to Sea Ice"
      expFieldStdName(  6) = "upward_sea_ice_basal_available_heat_flux"
      expFieldSclFac(   6) = -1.0  !field is downward
      expFieldUnits(    6) = "W m-2"
      expFieldName(     7) = "mlt"  !diagnostic
      expFieldLongName( 7) = "Ocean Mixed Layer Thickness"
      expFieldStdName(  7) = "ocean_mixed_layer_thickness"
      expFieldUnits(    7) = "m"
c
c  Create a DE layout to match HYCOM layout
      deLayout = ESMF_DELayoutCreate(vm,
     &           deCountList=(/ipr, jpr/),
     &           rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &   "Setup_ESMF: DELayoutCreate failed", rc))
     &   call ESMF_Finalize(rc=rc)
c
c  Create array specifications
      call ESMF_ArraySpecSet(arraySpec2Dr,
     &     rank=2,
     &     type=ESMF_DATA_REAL,
     &     kind=ESMF_R4,
     &     rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &   "Setup_ESMF: ArraySpecSet failed", rc))
     &   call ESMF_Finalize(rc=rc)
c
c  Create an ESMF grid that matches the HYCOM 2D grid
#if defined(ESMF_CURVILINEAR)
c  Use indices as cartesian coordinates
c
      cnts(    1)=itdm;       cnts(    2)=jtdm;
      mgcpd(   1)=1.d0;       mgcpd(   2)=1.d0;
      dpd(     1)=1.d0;       dpd(     2)=1.d0;
      dimNames(1)="index_i";  dimNames(2)="index_j";
      dimUnits(1)="1";        dimUnits(2)="1";
      periodic(1)=ESMF_FALSE; periodic(2)=ESMF_FALSE;
      gridName="HYCOM X-Y Cell-Center Grid"
      grid2D = ESMF_GridCreateHorzXYUni(
     &               counts     =cnts,
     &      minGlobalCoordPerDim=mgcpd,
     &               deltaPerDim=dpd,
     &               horzStagger=ESMF_GRID_HORZ_STAGGER_A, !CELL_CENTER only
     &               dimNames   =dimNames,
     &               dimUnits   =dimUnits,
     &               coordOrder =ESMF_COORD_ORDER_XYZ,
     &               periodic   =periodic,
     &               name       =trim(gridName),
     &               rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &   "Setup_ESMF: GridCreateHorzXYUni", rc))
     &   call ESMF_Finalize(rc=rc)
#else
c  Rectilinear, lat-lon, coordinates
      call xclget(plonrl,itdm, plon,1,1,1,0, 0)
      call xclget(platrl,jtdm, plat,1,1,0,1, 0)
      coord1(1:itdm) = plonrl(1:itdm)
      coord2(1:jtdm) = platrl(1:jtdm)
      dimNames(1)="longitude";    dimNames(2)="latitude";
      dimUnits(1)="degrees_east"; dimUnits(2)="degrees_north";
      if     (plonrl(itdm)-plonrl(1) .gt. 350.0) then !probably global
        periodic(1)=ESMF_TRUE
        periodic(2)=ESMF_FALSE
      else
        periodic(1)=ESMF_FALSE
        periodic(2)=ESMF_FALSE
      endif
      gridName="HYCOM Lat-Lon Cell-Center Grid"
      grid2D = ESMF_GridCreateHorzLatLon(
     &               coord1     =coord1(1:itdm),
     &               coord2     =coord2(1:jtdm),
     &               horzStagger=ESMF_GRID_HORZ_STAGGER_A, !CELL_CENTER only
     &               dimNames   =dimNames,
     &               dimUnits   =dimUnits,
     &               coordOrder =ESMF_COORD_ORDER_XYZ,
     &               periodic   =periodic,
     &               name       =trim(gridName),
     &               rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &   "Setup_ESMF: GridCreateHorzLatLon", rc))
     &   call ESMF_Finalize(rc=rc)
#endif
c
c  Distribute the grid on the DE layout using HYCOM distribution
c  Only Separable 2-D block distributions with no land skipping are allowed
      call ESMF_GridDistribute(grid2D,
     &     deLayout       =deLayout,
     &     countsPerDEDim1=countde1(1:ipr),
     &     countsPerDEDim2=countde2(1:jpr),
     &     rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &   "Setup_ESMF: GridDistribute", rc))
     &   call ESMF_Finalize(rc=rc)
c
c  Associate grid with ESMF gridded component
      call ESMF_GridCompSet(gridComp, grid=grid2D, rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &   "Setup_ESMF: GridCompSet", rc))
     &   call ESMF_Finalize(rc=rc)
c  Setup export fields, bundles & state
      do i = 1,numExpFields
        expField(i) = ESMF_FieldCreate(grid2D, arraySpec2Dr,
     &                     allocFlag =ESMF_ALLOC,
     &                     horzRelloc=ESMF_CELL_CENTER,
     &                     haloWidth =0,
     &                     name      =trim(expFieldName(i)),
     &                     rc=rc)
        call ESMF_FieldSetAttribute(expField(i),
     &       trim(attNameLongName), trim(expFieldLongName(i)), rc=rc)
        call ESMF_FieldSetAttribute(expField(i),
     &       trim(attNameStdName),  trim(expFieldStdName(i)),  rc=rc)
        call ESMF_FieldSetAttribute(expField(i),
     &       trim(attNameUnits),    trim(expFieldUnits(i)),    rc=rc)
        call ESMF_FieldSetAttribute(expField(i),
     &       trim(attNameAddOff),        expFieldAddOff(i),    rc=rc)
        call ESMF_FieldSetAttribute(expField(i),
     &       trim(attNameSclFac),        expFieldSclFac(i),    rc=rc)
        call ESMF_FieldGetDataPointer(expField(i),
     &       expData(i)%p, copyFlag=ESMF_DATA_REF, rc=rc)
        expData(i)%p(:,:) = 0.0
      enddo
c
c  Create bundle from list of fields
      expBundle = ESMF_BundleCreate(numExpFields,
     &                  expField(:), name="HYCOM",
     &                  rc=rc)
c
c  Add bundle to the export state
      call ESMF_StateAddBundle(expState, expBundle, rc=rc)
c
c  Setup import fields, bundles & state
      do i = 1,numImpFields
        impField(i) = ESMF_FieldCreate(grid2D, arraySpec2Dr,
     &                     allocFlag =ESMF_ALLOC,
     &                     horzRelloc=ESMF_CELL_CENTER,
     &                     haloWidth =0,
     &                     name      =trim(impFieldName(i)),
     &                     rc=rc)
        call ESMF_FieldSetAttribute(impField(i),
     &       trim(attNameLongName), trim(impFieldLongName(i)), rc=rc)
        call ESMF_FieldSetAttribute(impField(i),
     &       trim(attNameStdName),  trim(impFieldStdName(i)),  rc=rc)
        call ESMF_FieldSetAttribute(impField(i),
     &       trim(attNameUnits),    trim(impFieldUnits(i)),    rc=rc)
        call ESMF_FieldSetAttribute(impField(i),
     &       trim(attNameAddOff),        impFieldAddOff(i),    rc=rc)
        call ESMF_FieldSetAttribute(impField(i),
     &       trim(attNameSclFac),        impFieldSclFac(i),    rc=rc)
        call ESMF_FieldGetDataPointer(impField(i),
     &       impData(i)%p, copyFlag=ESMF_DATA_REF, rc=rc)
        impData(i)%p(:,:) = 0.0
      enddo
       sic_import(:,:) = 0.0 !Sea Ice Concentration
      sitx_import(:,:) = 0.0 !Sea Ice X-Stress
      sity_import(:,:) = 0.0 !Sea Ice Y-Stress
      siqs_import(:,:) = 0.0 !Solar Heat Flux thru Ice to Ocean
      sifh_import(:,:) = 0.0 !Ice Freezing/Melting Heat Flux
      sifs_import(:,:) = 0.0 !Ice Freezing/Melting Salt Flux
      sifw_import(:,:) = 0.0 !Ice Net Water Flux
       sit_import(:,:) = 0.0 !Sea Ice Temperature
       sih_import(:,:) = 0.0 !Sea Ice Thickness
       siu_import(:,:) = 0.0 !Sea Ice X-Velocity
       siv_import(:,:) = 0.0 !Sea Ice Y-Velocity
c
c  Create bundle from list of fields
      impBundle = ESMF_BundleCreate(numImpFields,
     &                  impField(:), name="HYCOM",
     &                  rc=rc)
c
c  Add bundle to the import state
      call ESMF_StateAddBundle(impState, impBundle, rc=rc)
c
      ocn_mask_init = .true.  !still need to initialize ocn_mask
c
      end subroutine Setup_ESMF
      
      subroutine Export_ESMF
c
c --- Fill export state.
c --- Calculate ssfi "in place"
c
      integer i,j,k
      real    ssh2m
      real    tmxl,smxl,umxl,vmxl,hfrz,tfrz,t2f,ssfi
      real    dp1,usur1,vsur1,psur1,dp2,usur2,vsur2,psur2
c
c --- Report
      call ESMF_LogWrite("HYCOM Export routine called",
     &     ESMF_LOG_INFO)
!-----call ESMF_LogFlush
c
      margin = 0
c
      if     (ocn_mask_init) then  !very 1st call to this routine
        ocn_mask_init = .false.
c
        if     (iceflg.eq.4) then
          ocn_mask(:,:) = 0.0  !export ocean currents nowhere
        elseif (nestfq.ne.0.0) then
c         export ocean currents away from open boundaries
          do j= 1,jj
            do i= 1,ii
              if     (rmunv(i,j).ne.0.0) then
                ocn_mask(i,j) = 0.0
              else
                ocn_mask(i,j) = 1.0
              endif
            enddo !i
          enddo !j
          do i= 1,10
            call psmooth(ocn_mask,0)  !not efficient, but only done once
          enddo !i
        else
          ocn_mask(:,:) = 1.0  !export ocean currents everywhere
        endif
      endif !ocn_mask_init
c
c --- Assume Export State is as defined in Setup_ESMF
c --- Average two time levels since (the coupling frequency) icefrq >> 2
c
      ssh2m = 1.0/g
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
c ---       quantities for available freeze/melt heat flux
c ---       relax to tfrz with e-folding time of icefrq time steps
c ---       assuming the effective surface layer thickness is hfrz
c ---       multiply by dpbl(i,j)/hfrz to get the actual e-folding time
            hfrz = min( thkfrz*onem, dpbl(i,j) )
            t2f  = (spcifh*hfrz)/(baclin*icefrq*g)
            smxl = 0.5*(saln(i,j,1,n)+saln(i,j,1,m))
            tmxl = 0.5*(temp(i,j,1,n)+temp(i,j,1,m))
            tfrz = tfrz_0 + smxl*tfrz_s  !salinity dependent freezing point
            ssfi = (tfrz-tmxl)*t2f       !W/m^2 into ocean
c ---       average currents over top 10m
            usur1 = 0.0
            vsur1 = 0.0
            psur1 = 0.0
            usur2 = 0.0
            vsur2 = 0.0
            psur2 = 0.0
            do k= 1,kk
              dp1   = min( dp(i,j,k,1), max( 0.0, tenm - psur1 ) )
              usur1 = usur1 + dp1*(u(i,j,k,1)+u(i+1,j,k,1))
              vsur1 = vsur1 + dp1*(v(i,j,k,1)+v(i,j+1,k,1))
              psur1 = psur1 + dp1
              dp2   = min( dp(i,j,k,2), max( 0.0, tenm - psur2 ) )
              usur2 = usur2 + dp2*(u(i,j,k,2)+u(i+1,j,k,2))
              vsur2 = vsur2 + dp2*(v(i,j,k,2)+v(i,j+1,k,2))
              psur2 = psur2 + dp2
              if     (max(dp1,dp2).eq.0.0) then
                exit
              endif
            enddo
            umxl  = 0.25*( usur1/psur1 + ubavg(i,  j,1) +
     +                                   ubavg(i+1,j,1) +
     &                     usur2/psur2 + ubavg(i,  j,2) +
     &                                   ubavg(i+1,j,2)  )
            vmxl  = 0.25*( vsur1/psur1 + vbavg(i,j,  1) +
     &                                   vbavg(i,j+1,1) +
     &                     vsur2/psur2 + vbavg(i,j,  2) +
     &                                   vbavg(i,j+1,2)  )
            expData(1)%p(i,j) = tmxl
            expData(2)%p(i,j) = smxl
            expData(3)%p(i,j) = umxl*ocn_mask(i,j)
            expData(4)%p(i,j) = vmxl*ocn_mask(i,j)
            expData(5)%p(i,j) = ssh2m*srfhgt(i,j)  !ssh in m
            expData(6)%p(i,j) = max(-1000.0,min(1000.0,ssfi))  !as in CICE
            expData(7)%p(i,j) = dpbl(i,j)*qonem
          endif !ip
        enddo !i
      enddo !j
c
      end subroutine Export_ESMF
      
      subroutine Import_ESMF
c
c --- Extract import state.
c
      integer i,j
c
c --- Report
      call ESMF_LogWrite("HYCOM Import routine called",
     &     ESMF_LOG_INFO)
!-----call ESMF_LogFlush
c
c --- Assume Import State is as defined in Setup_ESMF
c
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
             sic_import(i,j) = impData( 1)%p(i,j) !Sea Ice Concentration
            sitx_import(i,j) = impData( 2)%p(i,j) !Sea Ice X-Stress
            sity_import(i,j) = impData( 3)%p(i,j) !Sea Ice Y-Stress
            siqs_import(i,j) = impData( 4)%p(i,j) !Solar Heat Flux thru Ice to Ocean
            sifh_import(i,j) = impData( 5)%p(i,j) !Ice Freezing/Melting Heat Flux
            sifs_import(i,j) = impData( 6)%p(i,j) !Ice Freezing/Melting Salt Flux
            sifw_import(i,j) = impData( 7)%p(i,j) !Ice Net Water Flux
             sit_import(i,j) = impData( 8)%p(i,j) !Sea Ice Temperature
             sih_import(i,j) = impData( 9)%p(i,j) !Sea Ice Thickness
             siu_import(i,j) = impData(10)%p(i,j) !Sea Ice X-Velocity
             siv_import(i,j) = impData(11)%p(i,j) !Sea Ice Y-Velocity
            if     (iceflg.ge.2 .and. icmflg.ne.3) then
              covice(i,j) = impData(1)%p(i,j) !Sea Ice Concentration
                si_c(i,j) = impData(1)%p(i,j) !Sea Ice Concentration
              if     (covice(i,j).gt.0.0) then
                 si_tx(i,j) = -impData( 2)%p(i,j) !Sea Ice X-Stress into ocean
                 si_ty(i,j) = -impData( 3)%p(i,j) !Sea Ice Y-Stress into ocean
                fswice(i,j) =  impData( 4)%p(i,j) !Solar Heat Flux thru Ice to Ocean
                flxice(i,j) =  fswice(i,j) +
     &                         impData( 5)%p(i,j) !Ice Freezing/Melting Heat Flux
                sflice(i,j) =  impData( 6)%p(i,j)*1.e3 -
     &                         impData( 7)%p(i,j)*saln(i,j,1,n)
                                                  !Ice Virtual Salt Flux
                temice(i,j) =  impData( 8)%p(i,j) !Sea Ice Temperature
                thkice(i,j) =  impData( 9)%p(i,j) !Sea Ice Thickness
                  si_u(i,j) =  impData(10)%p(i,j) !Sea Ice X-Velocity
                  si_v(i,j) =  impData(11)%p(i,j) !Sea Ice Y-Velocity
              else
                 si_tx(i,j) = 0.0
                 si_ty(i,j) = 0.0
                fswice(i,j) = 0.0
                flxice(i,j) = 0.0
                sflice(i,j) = 0.0
                temice(i,j) = 0.0
                thkice(i,j) = 0.0
                  si_u(i,j) = 0.0
                  si_v(i,j) = 0.0
              endif !covice
            elseif (iceflg.ge.2 .and. icmflg.eq.3) then
                  si_c(i,j) =  impData( 1)%p(i,j) !Sea Ice Concentration
              if     (si_c(i,j).gt.0.0) then
                 si_tx(i,j) = -impData( 2)%p(i,j) !Sea Ice X-Stress into ocean
                 si_ty(i,j) = -impData( 3)%p(i,j) !Sea Ice Y-Stress into ocean
                  si_h(i,j) =  impData( 9)%p(i,j) !Sea Ice Thickness
                  si_t(i,j) =  impData( 8)%p(i,j) !Sea Ice Temperature
                  si_u(i,j) =  impData(10)%p(i,j) !Sea Ice X-Velocity
                  si_v(i,j) =  impData(11)%p(i,j) !Sea Ice Y-Velocity
              else
                 si_tx(i,j) = 0.0
                 si_ty(i,j) = 0.0
                  si_h(i,j) = 0.0
                  si_t(i,j) = 0.0
                  si_u(i,j) = 0.0
                  si_v(i,j) = 0.0
              endif !covice
            endif !iceflg>=2 (icmflg)
          endif !ip
        enddo !i
      enddo !j
c
      end subroutine Import_ESMF

C      subroutine Archive_ESMF(iyear,jday,ihour)
C      integer   iyear,jday,ihour
Cc
Cc --- Create a HYCOM "archive-like" file from Import/Export state..
Cc --- Import state may not be at the same time as Export.
Cc
C      character*8  cname
C      character*80 cfile
C      integer      i,j,k,nop,nopa
C      real         coord,xmin,xmax
Cc
C      write(cfile,'(a,i4.4,a1,i3.3,a1,i2.2)')
C     &   'arche.',iyear,'_',jday,'_',ihour
C      nopa=13
C      nop =13+uoff
Cc
C      call zaiopf(trim(cfile)//'.a', 'new', nopa)
C      if     (mnproc.eq.1) then
C      open (unit=nop,file=trim(cfile)//'.b',status='new') !uoff+13
C      write(nop,116) ctitle,iversn,iexpt,yrflag,itdm,jtdm
C      call flush(nop)
C      endif !1st tile
C 116  format (a80/a80/a80/a80/
C     & i5,4x,'''iversn'' = hycom version number x10'/
C     & i5,4x,'''iexpt '' = experiment number x10'/
C     & i5,4x,'''yrflag'' = days in year flag'/
C     & i5,4x,'''idm   '' = longitudinal array size'/
C     & i5,4x,'''jdm   '' = latitudinal  array size'/
C     & 'field       time step  model day',
C     & '  k  dens        min              max')
Cc
Cc --- surface fields
Cc
C      coord=0.0
C      do k= 1,numExpFields
C        do j= 1,jj
C          do i= 1,ii
C            if     (ip(i,j).eq.1) then
C              util1(i,j) = expData(k)%p(i,j)
C            endif !ip
C          enddo !i
C        enddo !j
C        cname = expFieldName(k)(1:8)
C        call zaiowr(util1,ip,.true.,
C     &              xmin,xmax, nopa, .false.)
C        if     (mnproc.eq.1) then
C        write (nop,117) cname,nstep,time,0,coord,xmin,xmax
C        call flush(nop)
C        endif !1st tile
C      enddo !k
C      do j= 1,jj
C        do i= 1,ii
C          if     (ip(i,j).eq.1) then
C            util2(i,j) = impData(1)%p(i,j)  !ice concentration
C          else
C            util2(i,j) = 0.0
C          endif !ip
C        enddo !i
C      enddo !j
C      cname = impFieldName(1)(1:8)
C      call zaiowr(util2,ip,.true.,
C     &            xmin,xmax, nopa, .false.)
C      if     (mnproc.eq.1) then
C      write (nop,117) cname,nstep,time,0,coord,xmin,xmax
C      call flush(nop)
C      endif !1st tile
C      do k= 2,3 !si_tx,si_ty
C        do j= 1,jj
C          do i= 1,ii
C            if     (util2(i,j).ne.0.0) then
C              util1(i,j) = -impData(k)%p(i,j)  !into ocean
C            else
C              util1(i,j) = 0.0
C            endif !ice:no-ice
C          enddo !i
C        enddo !j
C        cname = impFieldName(k)(1:8)
C        call zaiowr(util1,ip,.true.,
C     &              xmin,xmax, nopa, .false.)
C        if     (mnproc.eq.1) then
C        write (nop,117) cname(1:4)//'down',
C     &                  nstep,time,0,coord,xmin,xmax
C        call flush(nop)
C        endif !1st tile
C      enddo !k
C      do k= 4,7 !fluxes
C        do j= 1,jj
C          do i= 1,ii
C            if     (util2(i,j).ne.0.0) then
C              util1(i,j) = util2(i,j)*impData(k)%p(i,j)
C            else
C              util1(i,j) = huge !mask where there is no ice
C            endif !ice:no-ice
C          enddo !i
C        enddo !j
C        cname = impFieldName(k)(1:8)
C        call zaiowr(util1,ip,.false.,  !mask on ice
C     &              xmin,xmax, nopa, .false.)
C        if     (mnproc.eq.1) then
C        write (nop,117) cname,nstep,time,0,coord,xmin,xmax
C        call flush(nop)
C        endif !1st tile
C      enddo !k
C      do k= 8,numImpFields
C        do j= 1,jj
C          do i= 1,ii
C            if     (util2(i,j).ne.0.0) then
C              util1(i,j) = impData(k)%p(i,j)
C            else
C              util1(i,j) = huge !mask where there is no ice
C            endif !ice:no-ice
C          enddo !i
C        enddo !j
C        cname = impFieldName(k)(1:8)
C        call zaiowr(util1,ip,.false.,  !mask on ice
C     &              xmin,xmax, nopa, .false.)
C        if     (mnproc.eq.1) then
C        write (nop,117) cname,nstep,time,0,coord,xmin,xmax
C        call flush(nop)
C        endif !1st tile
C      enddo !k
C      do j= 1,jj
C        do i= 1,ii
C          if     (util2(i,j).ne.0.0) then
C            util1(i,j) = impData( 6)%p(i,j)*1.e3 -
C     &                   impData( 7)%p(i,j)*saln(i,j,1,n)  !virtual salt flux
C          else
C            util1(i,j) = huge !mask where there is no ice
C          endif !ice:no-ice
C        enddo !i
C      enddo !j
C      cname = 'surtx   '
C      call zaiowr(surtx,ip,.true.,
C     &            xmin,xmax, nopa, .false.)
C      if     (mnproc.eq.1) then
C      write (nop,117) cname,nstep,time,0,coord,xmin,xmax
C      call flush(nop)
C      endif !1st tile
C      cname = 'surty   '
C      call zaiowr(surty,ip,.true.,
C     &            xmin,xmax, nopa, .false.)
C      if     (mnproc.eq.1) then
C      write (nop,117) cname,nstep,time,0,coord,xmin,xmax
C      call flush(nop)
C      endif !1st tile
C      cname = 'sflice  '
C      call zaiowr(util1,ip,.false.,  !mask on ice
C     &            xmin,xmax, nopa, .false.)
C      if     (mnproc.eq.1) then
C      write (nop,117) cname,nstep,time,0,coord,xmin,xmax
C      call flush(nop)
C      endif !1st tile
C 117  format (a8,' =',i11,f11.3,i3,f7.3,1p2e16.7)
Cc
C      close (unit=nop)
C      call zaiocl(nopa)
c
C      end subroutine Archive_ESMF
#endif  /* ! defined(USE_ESMF_NERSC) */
      end module mod_hycom
