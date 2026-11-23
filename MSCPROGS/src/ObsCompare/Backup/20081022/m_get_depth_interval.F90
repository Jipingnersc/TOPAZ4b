module m_get_depth_interval

contains

subroutine  get_depth_interval(zinterval,nzint)
   implicit none

   real, pointer, dimension (:,:) :: zinterval
   integer, intent(out) :: nzint



   character(len=*), parameter  :: infile_depths   = 'depthintervals.in'
   character(len=20) :: tmpchar

   logical :: ex
   integer :: ios,counter
   real :: v1,v2


   ! Get depth levels to use
   inquire(exist=ex,file=infile_depths)
   if (.not. ex) then
      print '(a)','You must specify depth levels to use in'
      print '(a)','the file called  '//infile_depths
      stop '(get_depth_interval)'
   endif


   ! Start reading
   print *
   ios=0
   counter=1
   open(10,file=infile_depths,form='formatted')
   read(10,*,iostat=ios) nzint
   print '(i3.3,a)',nzint,' Depth Intervals:'
   allocate(zinterval(nzint,2))
   do while (ios==0 .and. counter<=nzint)
      read(10,*,iostat=ios) v1,v2
      if (ios==0) then
         zinterval(counter,1)=v1
         zinterval(counter,2)=v2
         print '(a,i3,2f10.2)','Depth interval :',counter,zinterval(counter,:)

         if (zinterval(counter,1) > zinterval(counter,2) .or. any(zinterval(counter,:)<0) ) then
            print *,' Zintervals: depths must be  positive in '//infile_depths
            print *,' Zintervals: upper (1) level must be higher than higher(2) level '
            stop '(get_depth_interval)'
         end if
      end if
      counter=counter+1
   end do
   close(10)
   if (counter/=nzint+1) then
      print *,'ERROR: Not all depth levels set'
      stop '(get_depth_interval)'
   endif
   if (ios/=0) then
      print *,'ERROR: occured when reading '//infile_depths
      stop '(get_depth_interval)'
   endif

end subroutine get_depth_interval

end module m_get_depth_interval
