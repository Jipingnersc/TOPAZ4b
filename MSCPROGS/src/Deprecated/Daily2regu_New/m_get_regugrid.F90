module m_get_regugrid
contains
subroutine get_regugrid(lon_sw,lon_ne,dellon, &
                        lat_sw,lat_ne,dellat)
   implicit none
   real, intent(out) :: lon_sw,lon_ne,dellon
   real, intent(out) :: lat_sw,lat_ne,dellat

   logical :: ex
   integer :: ios
   character(len=*), parameter :: infile_regugrid='regugrid.in'
   character(len=40)  :: tmpchar


   ! Get Regular grid specification
   inquire(exist=ex,file=infile_regugrid)
   if (.not. ex) then
      print '(a)','You must specify grid in'
      print '(a)','the file '//infile_regugrid
      stop '(get_regugrid)'
   endif


   ! Start reading
   open(10,file=infile_regugrid,form='formatted')
   read(10,*,iostat=ios) lon_sw, tmpchar
   print '(a,f10.3)','SW Lon   :',lon_sw
   read(10,*,iostat=ios) lon_ne, tmpchar
   print '(a,f10.3)','NE lon   :',lon_ne
   read(10,*,iostat=ios) dellon, tmpchar
   print '(a,f10.3)','Delta lon:', dellon
   read(10,*,iostat=ios) lat_sw, tmpchar
   print '(a,f10.3)','SW lat   :',lat_sw
   read(10,*,iostat=ios) lat_ne, tmpchar
   print '(a,f10.3)','NE lat   :',lat_ne
   read(10,*,iostat=ios) dellat, tmpchar
   print '(a,f10.3)','Delta lat:', dellat
   close(10)

   if (ios/=0) then
      print *,'ERROR: occured when reading '//infile_regugrid
      stop '(get_regugrid)'
   endif

   end subroutine get_regugrid
end module m_get_regugrid
