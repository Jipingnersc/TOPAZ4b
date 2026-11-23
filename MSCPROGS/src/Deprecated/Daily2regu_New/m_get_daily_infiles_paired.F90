module m_get_daily_infiles_paired
contains
subroutine get_daily_infiles_paired(daily_files_paired,numfiles,maxfiles)
   implicit none 

   integer , intent(out) :: maxfiles, numfiles
   character(len=20), dimension(:,:),pointer :: daily_files_paired


   character(len=*), parameter  :: infile_daily    = 'infiles_paired.daily'
   integer         , parameter  :: l_maxfiles=500

   logical :: ex
   integer :: ios, counter

   print *
   print '(a)','Reading infiles info:'

   ! Open infile for program -- files to process
   inquire(exist=ex,file=infile_daily)
   if (.not. ex) then
      print '(a)','You must specify files to process in'
      print '(a)','the file  '//infile_daily
      stop '(get_daily_infiles_paired)'
   endif


   ! Cycle through infile to get files to process
   ios=0
   counter=1
   maxfiles=l_maxfiles
   allocate(daily_files_paired(maxfiles,2))
   open(10,file=infile_daily,form='formatted')
   do while (ios==0 .and. counter<maxfiles)
      read(10,*,iostat=ios) daily_files_paired(counter,1) , daily_files_paired(counter,2)
      if (ios==0) then
         print '(a)','INPUT FILES: '//trim(daily_files_paired(counter,1))//'  '//trim(daily_files_paired(counter,1))
         counter=counter+1
      end if
   end do
   numfiles=counter-1
   if (counter==maxfiles) then
      print *,'WARNING: max num of files reached !'
      stop '(get_daily_files_paired)'
   end if
   print '(a,i4,a)','Read ',numfiles,' files'
   close(10)
end subroutine get_daily_infiles_paired
end module m_get_daily_infiles_paired
