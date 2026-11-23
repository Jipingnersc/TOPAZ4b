module mod_measurement
   type measurement
      real d                       ! Measurement value
      real var                     ! Error variance of measurement
      character(len=5) id          ! Type of measurement ('SST', 'SLA', 'COL', 'HICE', 'FICE', 'SAL', 'TEM'
      real lon                     ! Longitude position
      real lat                     ! Latitude position
      real depths                  ! depths of position 
      integer ipiv                 ! i-pivot point in grid
      integer jpiv                 ! j-pivot point in grid
      integer ns                   ! representativity in mod cells (meas. support)
                                   ! ns=0 means: point measurements
                                   ! used in m_Generate_element_Sij.F90
      real a1                      ! bilinear coeffisients (if ni=0)
      real a2                      ! bilinear coeffisients
      real a3                      ! bilinear coeffisients
      real a4                      ! bilinear coeffisients
      logical status               ! active or not
   end type measurement
end module mod_measurement

