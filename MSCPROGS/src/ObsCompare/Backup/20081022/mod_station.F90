module mod_station


   !KAL  Format of the infile: Example 1  with 3 sections - 
   !KAL  station_group name interpreted as section name

   !#Station_group_name_1
   !lon lat
   !lon lat
   !lon lat
   !lon lat
   !#Station_group_name_2
   !lon lat
   !lon lat
   !lon lat
   !lon lat
   !#Station_group_name_3
   !lon lat
   !lon lat
   !lon lat
   !lon lat


   !KAL  Format of the infile: Example 1  with 3 moorings - 
   !KAL  station_group name interpreted as mooring name

   !#Station_group_name_1
   !lon lat
   !#Station_group_name_2
   !lon lat
   !#Station_group_name_3
   !lon lat


   character(len=*), parameter  :: infile_stations   = 'stations.in'

   type station
      real lon
      real lat
      integer ipiv
      integer jpiv
      real a(4)
      real depth ! water depth
   end type station

   type(station), save,  allocatable, dimension(:,:) :: stations
   character(len=40) , save, dimension(:), allocatable :: groupname
   integer , save :: ngroup, maxstat
   integer , save, dimension(:), allocatable :: stations_per_group


contains

subroutine get_stations()
   use mod_dimensions
   use m_initconfmap
   use m_pivotp
   use m_bilincoeff
   use m_oldtonew
   implicit none
   logical :: ex
   integer :: ios,counter,prevsec
   integer :: nstat,ipiv,jpiv
   character(len=80) :: c80
   real :: lon,lat,lon_n,lat_n,a(4)

   ! Get stations to process
   inquire(exist=ex,file=infile_stations)
   if (.not. ex) then
      print '(a)','You must specify stations to use in'
      print '(a)','the file called  '//infile_stations
      stop '(get_stations)'
   endif

   ! Start reading to get a) number of groups b) max number of stations per group
   open(10,file=infile_stations,form='formatted')
   ios=0
   ngroup=0
   nstat=0
   maxstat=0
   read(10,'(a)',iostat=ios) c80
   do while (ios==0)
      print *,'|'//trim(c80)//'|'
      if (c80(1:1)=='#') then
         if (ngroup/=0) then
            maxstat=max(maxstat,nstat)
         end if
         nstat=0
         ngroup=ngroup+1
      else
         nstat=nstat+1
      end if
      read(10,'(a80)',iostat=ios) c80
   end do
   print *,ngroup, nstat
   maxstat=max(maxstat,nstat)
   close(10)

   if (maxstat==0 .or. ngroup==0)  then
      print *,'No  stations or no groups'
      print *,'Groups      :',ngroup
      print *,'Max stations:',maxstat
      stop '(get_stations)'
   end if

   call initconfmap

   ! Allocate 
   allocate(stations_per_group(ngroup))
   allocate(groupname         (ngroup))
   allocate(stations          (ngroup,maxstat))

   ! Start reading (again)
   open(10,file=infile_stations,form='formatted')
   !rewind(10)
   ios=0
   ngroup=0
   nstat=0
   read(10,'(a80)',iostat=ios) c80
   print *,ios
   do while (ios==0)
      if (c80(1:1)=='#') then
         nstat=0
         ngroup=ngroup+1
         groupname(ngroup)=c80(2:len_trim(c80))
      else 
         read(c80,*) lon,lat
         nstat=nstat+1
         print *,ngroup,nstat,lon,lat
         stations_per_group(ngroup)=nstat
         stations(ngroup,nstat)%lon=lon
         stations(ngroup,nstat)%lat=lat
         print *,'ok'

         ! corresponding pivot points
         call oldtonew(lat,lon,lat_n,lon_n)
         call pivotp(lon_n,lat_n,ipiv,jpiv)

         stations(ngroup,nstat)%ipiv=ipiv
         stations(ngroup,nstat)%jpiv=jpiv
         print *,'ok2'

         ! Bilinear interpolation factors
         call bilincoeff(modlon,modlat,nx,ny,lon,lat,ipiv,jpiv,&
                          a(1),a(2),a(3),a(4))

         stations(ngroup,nstat)%a=a
         print *,'ok3'


      end if
      print *,trim(c80)
      read(10,'(a80)',iostat=ios) c80
   end do
   close(10)
end subroutine get_stations

real function fld2d_to_station(fld2d,station2)
use mod_dimensions
implicit none
real, intent(in) , dimension(nx,ny) :: fld2d
type(station), intent(in) :: station2

integer :: i,j

! Simple closest point for now
i=station2%ipiv
j=station2%jpiv
fld2d_to_station = fld2d(i,j)
end function


end module mod_station

