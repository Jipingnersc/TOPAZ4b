      module mod_nerscice_esmf
      use ESMF_Mod    ! ESMF  Framework
      use mod_xc
c
c --- -----------------------------------------
c --- Coupling framework for the NERSC ice model
c ---        For now only set up with 
c ---               H Y C O M
c ---           v e r s i o n  2.2
c --- -----------------------------------------
c
      implicit none
c
      public NERSCICE_SetServices
      public NERSCICE_Init, NERSCICE_Run, NERSCICE_Final
c
      logical, save, public  :: put_export   !set in main program
      logical, save, public  :: get_import   !set in main program
      logical, save, public  :: end_of_run   !set in NERSCICE_Run
      integer, save, public  :: nts_day      !set in NERSCICE_init, timesteps/day
      integer, save, public  :: nts_ice      !set in NERSCICE_init, timesteps/ice
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
      integer, parameter :: numImpFields=12
      character(ESMF_MAXSTR), save :: impFieldName(    numImpFields),
     &                                impFieldLongName(numImpFields),
     &                                impFieldStdName( numImpFields),
     &                                impFieldUnits(   numImpFields)
      real(ESMF_KIND_R4),     save :: impFieldSclFac(  numImpFields),
     &                                impFieldAddOff(  numImpFields)
c
c --- Export Fields
      integer, parameter :: numExpFields=14
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
     &      ocn_sst   !Sea Ice Concentration
     &,     ocn_sss   !Sea Ice X-Stress
     &,     ocn_srfu  !Sea Ice Y-Stress
     &,     ocn_srfv  !Solar Heat Flux thru Ice to Ocean
     &,     ocn_ssh   !Ice Freezing/Melting Heat Flux
     &,     ocn_Qfi   !Ice Freezing/Melting Salt Flux
     &,     ocn_mld   !Ice Net Water Flux
     &,     air_cloud !Sea Ice Temperature
     &,     air_tair  !Sea Ice Thickness
     &,     air_slp   !Sea Ice X-Velocity
     &,     air_rlhum !Sea Ice Y-Velocity
     &,     air_uwind !Ocean Currents Mask
     &,     air_vwind !Ocean Currents Mask
     &,     air_shwflx!Ocean Currents Mask

      logical, save ::
     &       ocn_mask_init

      contains

      subroutine NERSCICE_SetServices(gridComp, rc)
c
      type(ESMF_GridComp) :: gridComp
      integer             :: rc
c
      call ESMF_GridCompSetEntryPoint(
     &     gridComp,
     &     ESMF_SETINIT,
     &     NERSCICE_Init,
     &     ESMF_SINGLEPHASE,
     &     rc)
      call ESMF_GridCompSetEntryPoint(
     &     gridComp,
     &     ESMF_SETRUN,
     &     NERSCICE_Run,
     &     ESMF_SINGLEPHASE,
     &     rc)
      call ESMF_GridCompSetEntryPoint(
     &     gridComp,
     &     ESMF_SETFINAL,
     &     NERSCICE_Final,
     &     ESMF_SINGLEPHASE,
     &     rc)
c
      end subroutine HYCOM_SetServices
C
      subroutine Setup_NERSCICE(gridComp, impState, 
                                expState, extClock, rc)
c
c --- Calling parameters
      type(ESMF_GridComp) :: gridComp
      type(ESMF_State)    :: impState
      type(ESMF_State)    :: expState
      type(ESMF_Clock)    :: extClock
      integer             :: rc
c
c --- set up ESMF data structures for NERSCICE
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
c  Attributes for export fields, identical to hycom import fields
      expFieldAddOff(:) = 0.0 !default is no offset
      expFieldSclFac(:) = 1.0 !default is no scale factor
c
c --- Later on we will remove some of these
c
      expFieldName(     1) = "sic"
      expFieldLongName( 1) = "Sea Ice Concentration"
      expFieldStdName(  1) = "sea_ice_area_fraction"
      expFieldUnits(    1) = "1"
      expFieldName(     2) = "sitx"
      expFieldLongName( 2) = "Sea Ice X-Stress"
      expFieldStdName(  2) = "downward_x_stress_at_sea_ice_base"
      expFieldSclFac(   2) = -1.0  !field is upward
      expFieldUnits(    2) = "Pa"
      expFieldName(     3) = "sity"
      expFieldLongName( 3) = "Sea Ice Y-Stress"
      expFieldStdName(  3) = "downward_y_stress_at_sea_ice_base"
      expFieldSclFac(   3) = -1.0  !field is upward
      expFieldUnits(    3) = "Pa"
      expFieldName(     4) = "siqs"
      expFieldLongName( 4) = "Solar Heat Flux thru Ice to Ocean"
      expFieldStdName(  4) = "downward_sea_ice_basal_solar_heat_flux"
      expFieldUnits(    4) = "W m-2"
      expFieldName(     5) = "sifh"
      expFieldLongName( 5) = "Ice Freezing/Melting Heat Flux"
      expFieldStdName(  5) = "upward_sea_ice_basal_heat_flux"
      expFieldSclFac(   5) = -1.0  !field is downward
      expFieldUnits(    5) = "W m-2"
      expFieldName(     6) = "sifs"
      expFieldLongName( 6) = "Ice Freezing/Melting Salt Flux"
      expFieldStdName(  6) = "downward_sea_ice_basal_salt_flux"
      expFieldUnits(    6) = "kg m-2 s-1"
      expFieldName(     7) = "sifw"
      expFieldLongName( 7) = "Ice Net Water Flux"
      expFieldStdName(  7) = "downward_sea_ice_basal_water_flux"
      expFieldUnits(    7) = "kg m-2 s-1"
      expFieldName(     8) = "sit" !diagnostic
      expFieldLongName( 8) = "Sea Ice Temperature"
      expFieldStdName(  8) = "sea_ice_temperature"
      expFieldAddOff(   8) = +273.15 !field is in degC
      expFieldUnits(    8) = "K"
      expFieldName(     9) = "sih" !diagnostic
      expFieldLongName( 9) = "Sea Ice Thickness"
      expFieldStdName(  9) = "sea_ice_thickness"
      expFieldUnits(    9) = "m"
      expFieldName(    10) = "siu" !diagnostic
      expFieldLongName(10) = "Sea Ice X-Velocity"
      expFieldStdName( 10) = "sea_ice_x_velocity"
      expFieldUnits(   10) = "m s-1"
      expFieldName(    11) = "siv" !diagnostic
      expFieldLongName(11) = "Sea Ice Y-Velocity"
      expFieldStdName( 11) = "sea_ice_y_velocity"
      expFieldUnits(   11) = "m s-1"
c
      expFieldName(    12) = "surflx"
      expFieldLongName(12) = "X"
      expFieldStdName( 12) = "X"
      expFieldUnits(   12) = "Wm-2"
      expFieldName(    13) = "salflx"
      expFieldLongName(13) = "X"
      expFieldStdName( 13) = "X"
      expFieldUnits(   13) = "g m-2 s-1"
      expFieldName(    14) = "shwflx"
      expFieldLongName(14) = "X"
      expFieldStdName( 14) = "X"
      expFieldUnits(   14) = "W m-2"
c
c  Attributes for export fields, identical to CICE import fields
      impFieldAddOff(:) = 0.0 !default is no offset
      impFieldSclFac(:) = 1.0 !default is no scale factor

      impFieldName(     1) = "sst"
      impFieldLongName( 1) = "Sea Surface Temperature"
      impFieldStdName(  1) = "sea_surface_temperature"
      impFieldAddOff(   1) = +273.15 !field is in degC
      impFieldUnits(    1) = "K"
      impFieldName(     2) = "sss"
      impFieldLongName( 2) = "Sea Surface Salinity"
      impFieldStdName(  2) = "sea_surface_salinity"
      impFieldUnits(    2) = "1e-3"
      impFieldName(     3) = "ssu"
      impFieldLongName( 3) = "Sea Surface X-Current"
      impFieldStdName(  3) = "sea_water_x_velocity"
      impFieldUnits(    3) = "m s-1"
      impFieldName(     4) = "ssv"
      impFieldLongName( 4) = "Sea Surface Y-Current"
      impFieldStdName(  4) = "sea_water_y_velocity"
      impFieldUnits(    4) = "m s-1"
      impFieldName(     5) = "ssh"
      impFieldLongName( 5) = "Sea Surface Height"
      impFieldStdName(  5) = "sea_surface_height_above_sea_level"
      impFieldUnits(    5) = "m"
      impFieldName(     6) = "ssfi"
      impFieldLongName( 6) = "Oceanic Heat Flux Available to Sea Ice"
      impFieldStdName(  6) = "upward_sea_ice_basal_available_heat_flux"
      impFieldSclFac(   6) = -1.0  !field is downward
      impFieldUnits(    6) = "W m-2"
      impFieldName(     7) = "mlt"  !diagnostic
      impFieldLongName( 7) = "Ocean Mixed Layer Thickness"
      impFieldStdName(  7) = "ocean_mixed_layer_thickness"
      impFieldUnits(    7) = "m"
c ---
c --- KAL - HYCOM reads the atmospheric stuff, so lets pass that too
c ---
      impFieldName(     8) = "clouds"  !diagnostic
      impFieldLongName( 8) = "X"
      impFieldStdName(  8) = "X"
      impFieldUnits(    8) = ""
      impFieldName(     9) = "tair"  !diagnostic
      impFieldLongName( 9) = "X"
      impFieldStdName(  9) = "X"
      impFieldUnits(    9) = ""
      impFieldName(    10) = "slp"  !diagnostic
      impFieldLongName(10) = "X"
      impFieldStdName( 10) = "X"
      impFieldUnits(   10) = ""
      impFieldName(    11) = "relhum"  !diagnostic
      impFieldLongName(11) = "X"
      impFieldStdName( 11) = "X"
      impFieldUnits(   11) = ""
      impFieldName(    12) = "uwind"  !diagnostic
      impFieldLongName(12) = ""
      impFieldStdName( 12) = ""
      impFieldUnits(   12) = ""
      impFieldName(    13) = "vwind"  !diagnostic
      impFieldLongName(12) = ""
      impFieldStdName( 12) = ""
      impFieldUnits(   12) = ""
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
     &                  expField(:), name="NERSCICE",
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
      ocn_sst   =0.0
      ocn_sss   =0.0
      ocn_srfu  =0.0
      ocn_srfv  =0.0
      ocn_ssh   =0.0
      ocn_Qfi   =0.0
      ocn_mld   =0.0
      air_cloud =0.0
      air_tair  =0.0
      air_slp   =0.0
      air_rlhum =0.0
      air_uwind =0.0
      air_vwind =0.0
      air_shwflx=0.0
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
C

c --- ----------------------------------------------------------
c --- 
c --- Subroutine exports variables from NERSCICE to HYCOM
c --- 
c --- ----------------------------------------------------------
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
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
c ---       quantities for available freeze/melt heat flux
c ---       relax to tfrz with e-folding time of icefrq time steps
c ---       assuming the effective surface layer thickness is hfrz
c ---       multiply by dpbl(i,j)/hfrz to get the actual e-folding time
            expData(1 )%p(i,j) = ficem(i,j)
            expData(2 )%p(i,j) = tauxice(i,j)
            expData(3 )%p(i,j) = tauyice(i,j)
            expData(4 )%p(i,j) = 0.0
            expData(5 )%p(i,j) = 0.0
            expData(6 )%p(i,j) = 0.0
            expData(7 )%p(i,j) = 0.0
            expData(8 )%p(i,j) = 0.0
            expData(9 )%p(i,j) = hicem(i,j)
            expData(10)%p(i,j) = uice(i,j)
            expData(11)%p(i,j) = vice(i,j)
            expData(12)%p(i,j) = surflx(i,j)
            expData(13)%p(i,j) = 0.0
            expData(14)%p(i,j) = 0.0
          endif !ip
        enddo !i
      enddo !j
c
      end subroutine Export_ESMF

c --- ----------------------------------------------------------
c --- 
c --- Subroutine exports variables from HYCOM to NERSCICE 
c --- 
c --- ----------------------------------------------------------
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
            ocn_sst   = impData( 1)%p(i,j) 
            ocn_sss   = impData( 2)%p(i,j) 
            ocn_srfu  = impData( 3)%p(i,j)
            ocn_srfv  = impData( 4)%p(i,j) 
            ocn_ssh   = impData( 5)%p(i,j) 
            ocn_Qfi   = impData( 6)%p(i,j) 
            ocn_mld   = impData( 7)%p(i,j) 
            air_cloud = impData( 8)%p(i,j) 
            air_tair  = impData( 9)%p(i,j) 
            air_slp   = impData(10)%p(i,j) 
            air_rlhum = impData(11)%p(i,j) 
            air_uwind = impData(12)%p(i,j) 
            air_vwind = impData(13)%p(i,j) 
            air_shwflx= impData(14)%p(i,j) 
          endif !ip
        enddo !i
      enddo !j
c
      end subroutine Import_ESMF
C
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
Cc
C      end subroutine Archive_ESMF
 

c --- ----------------------------------------------------------
c --- 
c --- Subroutine exports variables from HYCOM to NERSCICE 
c --- 
c --- ----------------------------------------------------------
      subroutine NERSCICE_Init 
         (gridComp, impState, expState, extClock, rc)
      implicit none
c
c --- Calling parameters
      type(ESMF_GridComp) :: gridComp
      type(ESMF_State)    :: impState
      type(ESMF_State)    :: expState
      type(ESMF_Clock)    :: extClock
      integer             :: rc
c
c --- Initialize (before the 1st time step).
c
      character(ESMF_MAXSTR) :: msg
c
c --- Report
      call ESMF_LogWrite("NERSCICE initialize routine called",
     &     ESMF_LOG_INFO)
!-----call ESMF_LogFlush
c
c --- Get VM from gridComp
      call ESMF_GridCompGet(gridComp, vm=vm, rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &   "NERSCICE_Init: GridCompGet failed", rc))
     &   call ESMF_Finalize(rc=rc)
c
c --- Get VM info
      call ESMF_VMGet(vm,
     &     petCount=petCount, localPET=localPet,
     &     mpiCommunicator=mpiCommunicator, rc=rc)
      if (ESMF_LogMsgFoundError(rc,
     &   "NERSCICE_Init: VMGet failed", rc))
     &   call ESMF_Finalize(rc=rc)
      write(msg,'(a,i4)') "NERSCICE_Init: petCount = ",petCount
      call ESMF_LogWrite(msg, ESMF_LOG_INFO)
!-----call ESMF_LogFlush
c
c
c --- set up ESMF data structures
c
      call Setup_ESMF(gridComp, impState, expState, extClock, rc)
      if (ESMF_LogMsgFoundError(rc,
     &   "NERSCICE_Init: Setup_ESMF failed", rc))
     &   call ESMF_Finalize(rc=rc)
c
c
c     nod=14 
c     nstep=nstep1
c     if (mnproc.eq.1) then
c     write (lp,'(/2(a,f8.1),2(a,i9),a/)') 'model starts at day',
c    &   time0,', goes to day',time0+day2-day1,'   (steps',nstep1,
c    &   ' --',nstep2,')'
c     open (unit=nod,file='summary_out',status='unknown')
c     write(nod,'(/2(a,f8.1),2(a,i9),a/)') 'model starts at day',
c    &   time0,', goes to day',time0+day2-day1,'   (steps',nstep1,
c    &   ' --',nstep2,')'
c     endif !1st tile
c
      call Export_ESMF
c
      call icemodels_init()
c
      end subroutine NERSCICE_Init
C
C

c --- ----------------------------------------------------------
c --- 
c --- Subroutine runs one step of ice model
c --- 
c --- ----------------------------------------------------------
      subroutine NERSCICE_Run(gridComp, impState, expState, extClock,rc)
      implicit none
c
c --- Calling parameters
      type(ESMF_GridComp) :: gridComp
      type(ESMF_State)    :: impState
      type(ESMF_State)    :: expState
      type(ESMF_Clock)    :: extClock
      integer             :: rc

      if     (get_import) then     ! new ESMF Import fields
        call Import_ESMF
        if     (nstep.eq.nstep0+1) then
          time=dtime0
          call forday(dtime0,yrflag, iyear,jday,ihour)
          call Archive_ESMF(iyear,jday,ihour)
          time=dtime
        endif !initial ESMF Archive
      endif !get_import
      call icemodels_step
      if     (put_export) then
c
c ---   fill the ESMF Export State.
c
        call Export_ESMF
        call forday(dtime,yrflag, iyear,jday,ihour)
        call Archive_ESMF(iyear,jday,ihour)
      endif
      end subroutine NERSCICE_Run




c --- ----------------------------------------------------------
c --- 
c --- Subroutine finalizes NERSCICE ESMF run
c --- 
c --- ----------------------------------------------------------
      subroutine NERSCICE_Final
     &           (gridComp, impState, expState, extClock, rc)
c
c --- Calling parameters
      type(ESMF_GridComp) :: gridComp
      type(ESMF_State)    :: impState
      type(ESMF_State)    :: expState
      type(ESMF_Clock)    :: extClock
      integer             :: rc
c
c --- Report
      call ESMF_LogWrite("NERSCICE finalize routine called",
     &                   ESMF_LOG_INFO)
!-----call ESMF_LogFlush
c
c --- Destroy internal ocean clock
*     call ESMF_ClockDestroy(intClock, rc=rc)
c
c
      if (mnproc .eq. 1) then
        write(nod,'(a)') 'normal stop'
        call flush(nod)
      endif
c
      end subroutine NERSCICE_Final
      end module mod_nerscice_esmf
