module m_get_depth_levels

contains

subroutine  get_depth_levels(zlevels,nzlevel)
   implicit none

   real, pointer, dimension (:) :: zlevels
   integer, intent(out) :: nzlevel



   character(len=*), parameter  :: infile_depths   = 'depthlevels.in'
   character(len=20) :: tmpchar

   logical :: ex
   integer :: ios,counter
   real :: v1,v2




   ! Get depth levels to use
   inquire(exist=ex,file=infile_depths)
   if (.not. ex) then
      print '(a)','You must specify depth levels to use in'
      print '(a)','the file called  '//infile_depths
      stop '(get_depth_levels)'
   endif

   ! Start reading
   ios=0
   counter=1
   open(10,file=infile_depths,form='formatted')
   read(10,*,iostat=ios) nzlevel,tmpchar
   print '(i3,a)',nzlevel,' Depth levels:'
   allocate(zlevels(nzlevel))
   do while (ios==0 .and. counter<=nzlevel)
      read(10,*,iostat=ios) zlevels(counter)
      print '(a,i3,f10.2)','Depth level :',counter,zlevels(counter)
      counter=counter+1
   end do
   close(10)
   if (counter/=nzlevel+1) then
      print *,'ERROR: Not all depth levels set'
      stop '(get_depth_levels)'
   endif
   if (ios/=0) then
      print *,'ERROR: occured when reading '//infile_depths
      stop '(get_depth_levels)'
   endif



end subroutine get_depth_levels

end module m_get_depth_levels
