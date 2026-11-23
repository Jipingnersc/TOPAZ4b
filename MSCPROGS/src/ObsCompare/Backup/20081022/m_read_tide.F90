module m_read_tide
contains

   !Reads SOFA data (SST) from a NetCDF file

   subroutine read_tide(datacat,time1,time2,fld,times,longitude,latitude,undef)
      use m_datetojulian
      implicit none

      integer, intent(in) :: time1,time2
      character(len=*),            intent(in)  :: datacat
      real, pointer, dimension(:) :: longitude,latitude
      integer, pointer, dimension(:) :: times
      real, pointer, dimension(:,:) :: fld
      real, intent(in) :: undef

      character(len=80) :: file1, cline
      character(len=14) :: clon, clat
      character(len= 8) :: stationname
      character(len= 1) :: hemilon,hemilat
      logical :: ex,lfound
      integer :: startyear, endyear, j2d, ios, ios2, tide0(12),  &
                 hourind, month, day, year, stationid
      integer :: ilatdeg, ilondeg
      real    :: flatmin, flonmin,flat,flon
      integer :: nrfiles,nrtimes, tmpnrtimes
      integer :: year1,year2,month1,month2,day1,day2
      integer :: itest, ifile,itestfile, itmp, itmp2
      integer :: istart, iend ,minind,maxind, istart2, iend2
      integer :: istartday, iendday, tmpnrdays
      integer :: day3,month3, year3,itime
      real, allocatable :: tide1(:),tide2(:),tide3(:)
      logical :: firstpass



      ! Startyear and endyear denote the time period used for intermediate
      ! calculations. The interval time1-time2 is within this time period
      startyear=(time1/10000) - 1
      endyear  =(time2/10000) + 1

      ! year, month and days for the time interval -- time1 and time2
      ! are integers of format yyyymmdd
      year1=(time1/10000) 
      year2=(time2/10000) 
      month1=(time1-year1*10000)/100
      month2=(time2-year2*10000)/100
      day1=(time1-year1*10000 - month1*100)
      day2=(time2-year2*10000 - month2*100)




      ! nrtimes is the array size for the time interval time1 to time2 (days)
      j2d=datetojulian(year2,month2,day2,year1,month1,day1)
      nrtimes=j2d+1

      ! tmpnrtimes is the array size for time interval 
      ! 1.1.startyear to 31.12.endyear. (hours)
      j2d=datetojulian(endyear,12,31,startyear,1,1)
      tmpnrtimes=(j2d+1)*24


      print *,year1,month1,day1
      print *,year2,month2,day2
      print *,'nrtimes=',nrtimes
      tmpnrdays=(j2d+1)
      print *,'tmpnrtimes=',tmpnrtimes
      print *,'tmpnrdays =',tmpnrdays

      

      ! List Data catalogue contents -- dump them in a file
      call system("ls "//trim(datacat)//"*.dat > "//trim(datacat)//".tidefiles")

      ! Count number of entrys in that file -- gives array size
      open(10,file=trim(datacat)//".tidefiles")
      ios=0
      nrfiles=0
      do while( ios==0)
         read(10,'(a80)',iostat=ios) cline
         print *,cline
         if (ios==0) nrfiles=nrfiles+1
      end do
      if (nrfiles==0) then
         print *,'No data files found .. stopping'
         stop '(read_tide)'
      end if
      rewind(10)
      print *,'nrfiles=',nrfiles

      ! Allocate arrays -- final data holders
      allocate(longitude(nrfiles))
      allocate(latitude (nrfiles))
      allocate(fld      (nrfiles,nrtimes))
      allocate(times    (nrtimes))

      ! Allocate arrays -- temporary data holders
      allocate(tide1    (tmpnrtimes))
      allocate(tide2    (tmpnrtimes))
      allocate(tide3    (tmpnrdays))

      ! The times are in yyyymmdd format, we can process them at this stage:
      j2d=datetojulian(year2,month2,day2,year1,month1,day1)
      nrtimes=j2d+1
      print *,time1,time2
      do itime=1,nrtimes
         call juliantodate(itime-1,year3,month3,day3,year1,month1,day1)
         times(itime)=year3*10000+month3*100+day3
      end do

         




      ! Loop through the contents of the file containing datafiles
      ! -- Initialize fields
      itestfile=2
      ios=0
      ifile=0
      fld=undef
      do while(ios==0)
         ifile=ifile+1

         ! Read name of this datafile -- initialize temporary tide holders
         read(10,'(a80)',iostat=ios) file1
         tide1=undef
         tide2=undef
         tide3=undef

         if (ios==0) then 
            print *,file1
            open(11,file=trim(file1))

            ! In file1, read everything from "startyear" onwards -- skip to that time
            year=-1
            ios2=0
            do while (year <= startyear)
               read(11,1010,iostat=ios2) stationid, stationname,year
            end do

            ! Now we should be in our "startyear"-"endyear" time interval
            flon=undef
            flat=undef
            if (ios2==0) then 

               ! Step one back
               backspace(11)

               ! We should be ready to read tide station data for our time
               ! period -- read until we reach endyear
               ios2=0
               maxind=-1
               minind=tmpnrtimes+1
               do while (year <= endyear .and. ios2==0)

                  ! Read data line 
                  read(11,'(a80)',iostat=ios2) cline

                  if (ios2==0) then

                     ! Try to read line as data line
                     read(cline,1030,iostat=ios2) stationid, stationname,year,month,day,hourind,tide0

                     ! If ok -- put data into "tide0" -- raw unsmoothed data -- hourly values
                     if (ios2==0) then
                        
                        j2d=datetojulian(year,month,day,startyear,1,1)+1
                        !print *,'year, month,day,hour index ',year,month,day,hourind
                        if (hourind==1) then
                           tide1((j2d-1)*24+1:(j2d-1)*24+12) = tide0
                           maxind=max(maxind,(j2d-1)*24+12)
                           minind=min(minind,(j2d-1)*24+12)
                        elseif (hourind==2) then
                           tide1((j2d-1)*24+13:(j2d-1)*24+24) = tide0
                           maxind=max(maxind,(j2d-1)*24+24)
                           minind=min(minind,(j2d-1)*24+24)
                        end if



                     ! On error, try to read data line as year header
                     elseif (ios2/=0) then
                        read(cline,1020,iostat=ios2) stationid, stationname,year,clat,clon
                        if (ios2==0) then

                           firstpass=abs(flon-undef)<1e-4 .and.  abs(flat-undef)<1e-4

                           ! Convert to decimal lon/lat
                           read(clat,'("LAT=",i3,f4.1,a1)') ilatdeg,flatmin,hemilat
                           if (hemilat=="S") then
                              ilatdeg=-ilatdeg
                              flatmin=-flatmin
                           end if
                           flat=ilatdeg+flatmin/60.
                           read(clon,'(" LONG=",i3,f4.1,a1)') ilondeg,flonmin,hemilon
                           if (hemilon=="W") then
                              ilondeg=-ilondeg
                              flonmin=-flonmin
                           end if
                           flon=ilondeg+flonmin/60.
                           if (firstpass) then
                              print *,'Station id, name, year, lat, lon',&
                                 stationid,stationname,year,clat,clon
                              print *,'Longitude -- deg,min,hemisphere,float',ilondeg,flonmin,hemilon,flon
                              print *,'Latitude  -- deg,min,hemisphere,float',ilatdeg,flatmin,hemilat,flat
                           end if
                        end if

                     end if
                  end if
               end do

               ! We should be sitting with a full time series now... but it
               ! might have gaps... (value=undef)
               itmp=minind
               do while (itmp<=maxind)

                  istart=maxind+1
                  iend=minind-1

                  ! 1) Find first and last valid hourly entry
                  lfound=.false.
                  do while (itmp<=maxind .and..not.lfound)
                     if (abs(tide1(itmp)-undef)<1e-4) then
                        itmp=itmp+1
                     else
                        lfound=.true.
                        istart=itmp
                     end if
                  end do

                  lfound=.false.
                  do while (itmp<=maxind .and..not.lfound)
                     if (abs(tide1(itmp)-undef)>1e-4) then
                        iend=itmp
                        itmp=itmp+1
                     else
                        lfound=.true.
                     end if
                  end do

                  if (ifile==itestfile) then
                     print *,istart,iend,minind,maxind
                     !print *,tide1(istart:iend)
                  end if


                  ! 2) Smooth the sucker to remove fluctuations with period < 3
                  ! days -- tide2 is now the array containing smoothed data
                  if (iend-istart+1 > 3*24+1) then
                     call rosame_filtre_demerliac(iend-istart+1,iend-istart+1, &
                        tide1(istart:iend),undef,tide2(istart:iend))

                     !  The tide filter renders the first and last few indices invalid
                     ! Again, find first and last valid entry
                     lfound=.false.
                     itmp2=istart
                     do while (itmp2<=iend .and..not.lfound)
                        if (abs(tide2(itmp2)-undef)<1e-4) then
                           itmp2=itmp2+1
                        else
                           lfound=.true.
                           istart2=itmp2
                        end if
                     end do

                     lfound=.false.
                     do while (itmp2<=iend .and..not.lfound)
                        if (abs(tide2(itmp2)-undef)>1e-4) then
                           iend2=itmp2
                           itmp2=itmp2+1
                        else
                           lfound=.true.
                        end if
                     end do

                     ! Ready to put data into daily average array. This is the
                     ! first and last  "day indexes" which are valid (incomplete
                     ! days are discarded)
                     istartday=ceiling(float(istart2)/24.)+1
                     iendday=  floor  (float(iend2  )/24.)
                     print *,istartday,iendday,istart2,iend2,(istartday-1)*24+1, &
                       (iendday-1)*24+24

                     ! Daily averages -- mean value of "tide2" over a day
                     do itmp2=istartday,iendday
                        !print *,itmp2,tide2((itmp2-1)*24+1), tide2((itmp2-1)*24+24-1)
                        tide3(itmp2)=sum(tide2((itmp2-1)*24+1:(itmp2-1)*24+24))/24
                        tide3(itmp2)= tide3(itmp2)*1e-3 ! - mm -> m
                     end do
                        


                  end if
               end do

            end if
            close(11)


         end if

         ! -- We are ready to put the data into the data field
         ! 
         do itmp=1,nrtimes

            ! itmp is index into "fld". This is the corresponding date  (3)
            call juliantodate(itmp,year3,month3,day3,year1,month1,day1)

            ! day index into tide3
            j2d=datetojulian(year3,month3,day3,startyear,1,1)+1

            ! Set data 
            fld(ifile,itmp) = tide3(j2d)
         end do

         ! -- Put lon, lat into the data field
         if (abs(flon-undef)>1e-4 .and. abs(flat-undef)>1e-4) then
            longitude(ifile)=flon
            latitude (ifile)=flat
         else
            print *,'Error:Could not get position for station'
            stop '(m_read_tide)'
         end if



         ! Test
         if (ifile==itestfile) then
            open(13,file='testtide.dat',status='replace')
            do itest=minind,maxind
               write(13,'(i5,4f10.2)') itest,tide1(itest),tide2(itest),tide3((itest/24)+1)
            end do
            close(13)
         end if
      end do

      !print *
      !print *,'Times:'
      !print *,times
      !print *
      !print *,'Longitude:'
      !print *,longitude
      !print *
      !print *,'Latitude:'
      !print *,latitude
      !print *
      !print *,'Tide data'
      !print *,fld

      ! Dump final data for testing
      open(13,file='testtide_final.dat',status='replace')
      do itime=1,nrtimes
         year3=(times(itime)/10000) 
         month3=(times(itime)-year3*10000)/100
         day3=(times(itime)-year3*10000 - month3*100)
         j2d=datetojulian(year3,month3,day3,startyear,1,1)
         write(13,'(i5,i6,3i5,i12,f10.2)') itime,j2d,year3,month3,day3,times(itime),fld(itestfile,itime)
      end do
      close(13)

            
! Common format -- all lines
1010 format(i3,a8,i4)
! Format for a year header
1020 format(i3,a8,i4,2x,2a14)
! Format for a data line
1030 format(i3,a8,i4,2i2,i1,12i5)



   end subroutine read_tide



end module m_read_tide
