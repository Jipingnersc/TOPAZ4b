module m_get_daily_infiles
contains
subroutine get_daily_infiles(daily_files,numfiles,maxfiles)
   implicit none 

   integer , intent(out) :: maxfiles, numfiles
   character(len=*), dimension(:),pointer :: daily_files


   character(len=*), parameter  :: infile_daily    = 'infiles.daily'
   integer         , parameter  :: l_maxfiles=1000

   logical :: ex
   integer :: ios, counter

   print *
   print '(a)','Reading infiles info:'

   ! Open infile for program -- files to process
   inquire(exist=ex,file=infile_daily)
   if (.not. ex) then
      print '(a)','You must specify files to process in'
      print '(a)','the file  '//infile_daily
      stop '(get_daily_infiles)'
   endif


   ! Cycle through infile to get files to process
   ios=0
   counter=1
   maxfiles=l_maxfiles
   allocate(daily_files(maxfiles))
   open(10,file=infile_daily,form='formatted')
   do while (ios==0 .and. counter<maxfiles)
      read(10,*,iostat=ios) daily_files(counter)
      if (ios==0) then
         print '(a)','INPUT FILE: '//trim(daily_files(counter))
         counter=counter+1
      end if
   end do
   numfiles=counter-1
   if (counter==maxfiles) then
      print *,'WARNING: max num of files reached !'
      stop '(get_daily_files)'
   end if
   print '(a,i4,a)','Read ',numfiles,' files'
   close(10)
end subroutine get_daily_infiles
end module m_get_daily_infiles
