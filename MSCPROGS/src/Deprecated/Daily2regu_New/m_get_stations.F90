module m_get_stations
contains

subroutine get_stations(stat_lon,stat_lat,stat_sec_name,stat_sec_num, &
                        max_nstations, nstations, numsec, cformat)
   implicit none


   integer, intent(in) :: max_nstations
   integer, intent(out):: nstations, numsec
   real, dimension(max_nstations),intent(out) :: stat_lon
   real, dimension(max_nstations),intent(out) :: stat_lat
   character(len=*), dimension(max_nstations),intent(out) :: stat_sec_name
   integer, dimension(max_nstations),intent(out) :: stat_sec_num
   character(len=*), intent(in), optional :: cformat

   character(len=20) :: cformat2
   character(len=*), parameter  :: infile_secdef   = 'sectiondefs.in'
   character(len=*), parameter  :: infile_secdef2  = 'sectiondefsXBT.in'

   logical :: ex
   integer :: ios,counter,prevsec

   if (present(cformat)) cformat2=cformat

   ! Get stations to process
   inquire(exist=ex,file=infile_secdef)
   if (.not. ex) then
      print '(a)','You must specify stations to use in'
      print '(a)','the file called  '//infile_secdef
      stop '(get_stations)'
   endif

   ! Start reading
   ios=0
   counter=1
   prevsec=0
   !print *,'format :'
   !print 100,-10.0,10.0,'TEST',1
   if (.not.present(cformat)) then
      open(10,file=infile_secdef,form='formatted')
   elseif (trim(cformat2)=='NORMAL') then
      open(10,file=infile_secdef,form='formatted')
   elseif (trim(cformat2)=='XBT') then
      open(10,file=infile_secdef2,form='formatted')
   else
      print *,'Unknowm format specification : |'//trim(cformat2)//'|'
      stop '(get_stations (1))'
   end if
   do while (ios==0 .and. counter<=max_nstations)
      if (.not.present(cformat)) then
         read(10,100,iostat=ios)                  &
            stat_lon(counter), stat_lat(counter), &
            stat_sec_name(counter), stat_sec_num(counter)
      elseif (trim(cformat2)=='NORMAL') then
         read(10,100,iostat=ios)                  &
            stat_lon(counter), stat_lat(counter), &
            stat_sec_name(counter), stat_sec_num(counter)
      elseif (trim(cformat2)=='XBT') then
         read(10,200,iostat=ios)                  &
            stat_lon(counter), stat_lat(counter), &
            stat_sec_name(counter), stat_sec_num(counter)
            !print 200,                               &
            !   stat_lon(counter), stat_lat(counter), &
            !   stat_sec_name(counter), stat_sec_num(counter)
      else
         print *,'Unknowm format specification : |'//trim(cformat2)//'|'
         stop '(get_stations)'
      end if

      if (ios==0 ) then
         if ( prevsec/=stat_sec_num(counter)) then
            print 100,                               &
               stat_lon(counter), stat_lat(counter), &
               stat_sec_name(counter), stat_sec_num(counter)
         end if
         prevsec=stat_sec_num(counter)
         counter=counter+1
      end if
   end do
!100 format(f8.4,2x,f8.4,2x,a25,2x,i2,1x,i2)
100 format(f8.4,2x,f8.4,2x,a26,3x,i2)
200 format(f10.5,2x,f10.5,2x,a31,2x,i2)
   close(10)
   nstations=counter-1
   if (nstations>max_nstations) then
      print *,'ERROR: Too many stations'
      stop '(daily2sectionnc)'
   endif
   if (ios>0 ) then
      print *,'ERROR: occured when reading '//infile_secdef
      stop '(get_stations)'
   else if (nstations==0) then
      print *,'ERROR: No stations read from  '//infile_secdef
      stop '(get_stations)'
   endif
   
   ! NB: The following assumes each station from 1 to numsec
   ! has stations !
   numsec=maxval(stat_sec_num(1:nstations))
   print *,'A grand total of ',nstations,' over ',&
           numsec, ' sections'

end subroutine get_stations
end module m_get_stations
