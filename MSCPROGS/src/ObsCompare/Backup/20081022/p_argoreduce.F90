program p_argoreduce
   use mod_sphere_tools
   implicit none

   real, parameter :: undef=-999.

   integer, parameter :: maxfiles=20000
   integer :: ncorners
   integer :: ios
   character(len=80) :: c80
   character(len=30) :: cpos

   real   , dimension(maxfiles) :: lons,lats
   integer, dimension(maxfiles) :: dates
   character(len=80), dimension(maxfiles) :: infiles

   real, dimension(:), allocatable :: clons, clats,ptembias, psalbias
   integer :: psalnum, ptemnum
   integer :: depthrange(2)
   integer :: daterange(2)
   integer :: numfiles,numfiles2
   integer :: strindE,strindN,i,strindpos
   logical :: dateflag,regionflag
   integer :: numsalnpoints,numtemppoints
   real :: asaln, msaln, meanasaln,meanmsaln,depth, &
           biassaln, biastemp,rmssaln,rmstemp, atemp, &
           mtemp, meanatemp,meanmtemp
   integer :: firstyear, lastyear, firstmonth, lastmonth,imonth,iyear
   character(len=4) :: cyear
   character(len=2) :: cmonth,cdirct
   character(len=80) :: fulldir


   ! Input file lists 4D region of comparison
   ! 1) Corner lons
   ! 2) Corner lats
   ! 3) depth range
   ! 3) date range
   open(10,file='argoreduce.in',action='read')
   read(10,*) ncorners !; print *,ncorners
   allocate(clons(ncorners))
   allocate(clats(ncorners))
   read(10,*)clons !; print *,clons
   read(10,*)clats !; print *,clats
   read(10,*)depthrange !; print *,depthrange
   read(10,*)daterange !; print *,daterange
   close(10)

   ! List files in "ArgoCmp" directory/ies
   firstyear=daterange(1)/10000
   firstmonth=(daterange(1)-firstyear*10000)/100
   lastyear=daterange(2)/10000
   lastmonth=(daterange(2)-lastyear*10000)/100
   print *,'Firstyear, firstmonth:',firstyear,firstmonth
   print *,'Lastyear, lastmonth  :',lastyear,lastmonth
   !print *


   iyear=firstyear
   imonth=firstmonth
   cdirct=" >"
   do while (iyear<lastyear .or. (iyear==lastyear.and.imonth<=lastmonth)) 

      
     write(cyear,'(i4.4)') iyear
     write(cmonth,'(i2.2)') imonth
     fulldir="ArgoCmp/"//cyear//"/"//cmonth
     !print *,iyear,imonth, fulldir

     ! Augment month/year
     iyear=iyear+imonth/12
     imonth=mod(imonth,12)+1

     call system("find "//trim(fulldir)//" -type f"//cdirct//" ArgoCmpFiles")
      
      cdirct=">>"
   end do
      




   ! New -- list in subdirs w year/month
   !ymlow=daterange(1)/100
   !ymhi =daterange(2)/100
   !ym=ymlow
   !while (ym <= ymhi) do
   !   yy=ym/100
   !   mm=ym-yy*100000
   !   write(cyy,'(i4.4)') yy
   !   write(cmm,'(i2.2)') mm
   !   fullpath='ArgoCmp/'//cyy//cmm
   !   call system("dir ArgoCmp | grep -v profile_legend >> ArgoCmpFiles")
   !   yy = (mm+1)/12
   !   mm = mod(mm,12) + 1
   !   ym=yy*100000+mm
   !end do

   

   ! Go through list - deduce region and date from each
   open(10,file='ArgoCmpFiles',action='read')

   ios=0
   numfiles=0
   do while(ios==0)
      read(10,'(a80)',iostat=ios) c80
      !print *,c80

      strindpos=index(c80,'profile_date')
      !print *,strindpos
      !print *,c80(strindpos+12:strindpos+19)
      if (ios==0 .and. strindpos>0 ) then

         numfiles=numfiles+1
         ! Search for keyword "_profile_date" in file name
         read(c80(strindpos+12:strindpos+19),'(i8)') dates(numfiles)
         infiles(numfiles)=c80

         ! Search for keyword "_pos" in file name
         strindpos=index(c80,'_pos')
         cpos=c80(strindpos+4:len_trim(c80))


         ! Deduce lon/lat from remainder of string
         strindE=index(cpos,'E')
         strindN=index(cpos,'N')
         read(cpos(        1:strindE-1),*) lons(numfiles)
         read(cpos(strindE+2:strindN-1),*) lats(numfiles)

         !if (lons(numfiles)>0) then
         !   print *,trim(c80)//' '//cpos(1:strindE),' ',cpos(strindE+2:strindN)
         !   print *,cpos,dates(numfiles),lons(numfiles),lats(numfiles)
         !end if

         !print *,dates(numfiles),lons(numfiles),lats(numfiles)
      end if
   end do
   close(10)
   print '(a,i4)','Files present in ArgoCmp ',numfiles

   ! We have the files and positions - > Find candidat input files for reduction 4D
   ! box
   numfiles2=0
   do i=1,numfiles
      
      dateflag=dates(i)>=daterange(1) .and. dates(i)<=daterange(2)
      regionflag=inbox(clons,clats,ncorners,lons(i),lats(i))
      !print *,dateflag,regionflag

      if (dateflag.and.regionflag) then
         numfiles2=numfiles2+1
         dates(numfiles2)=dates(i)
         lons(numfiles2)=lons(i)
         lats(numfiles2)=lats(i)
         infiles(numfiles2)=infiles(i)

         !if (lons(numfiles2)>0) then
         !   print '(a,i10,2f10.2)',trim(infiles(numfiles2)),dates(numfiles2),lons(numfiles2),lats(numfiles2)
         !end if
      end if
   end do
   numfiles=numfiles2
   print '(a,i4)','Files present in ArgoCmp and matching criteria ',numfiles

   ! KAL -- for pointwise stats - only bias for now
   allocate(psalbias(numfiles)) 
   allocate(ptembias(numfiles))

   ! Now go through remaining files
   numsalnpoints=0
   numtemppoints=0
   meanasaln=0
   meanmsaln=0
   rmssaln=0
   rmstemp=0
   biassaln=0
   biastemp=0
   psalbias=0.
   ptembias=0.
   do i=1,numfiles

      psalnum=0
      ptemnum=0

      open(10,file=trim(infiles(i)),action='read')
      ios=0
      do while (ios==0)

         read(10,*,iostat=ios) ,depth,asaln,msaln,atemp,mtemp

         if (ios==0) then

            if(depth/=undef.and.asaln/=undef.and.msaln/=undef) then
               if (depth>=depthrange(1) .and. depth  <=depthrange(2)) then
                  meanasaln=meanasaln+asaln
                  meanmsaln=meanmsaln+msaln
                  rmssaln=rmssaln+(msaln-asaln)**2
                  biassaln=biassaln+(msaln-asaln)
                  numsalnpoints=numsalnpoints+1

                  ! Pointwise
                  psalbias(i)=psalbias(i)+(msaln-asaln) ; 
                  psalnum=psalnum+1
               end if
            end if

            if(depth/=undef.and.atemp/=undef.and.mtemp/=undef) then
               if (depth>=depthrange(1) .and. depth  <=depthrange(2)) then
                  meanatemp=meanatemp+atemp
                  meanmtemp=meanmtemp+mtemp
                  rmstemp=rmstemp+(mtemp-atemp)**2
                  biastemp=biastemp+(mtemp-atemp)
                  numtemppoints=numtemppoints+1

                  ! Pointwise
                  ptembias(i)=ptembias(i)+(mtemp-atemp) ; 
                  ptemnum=ptemnum+1
               end if
            end if

         end if



      end do
      close(10)
      if (ptemnum>0) then
         ptembias(i)=ptembias(i)/ptemnum
      else
         ptembias(i)=-999;
      end if
      if (psalnum>0) then
         psalbias(i)=psalbias(i)/psalnum
      else
         psalbias(i)=-999;
      end if
   end do
   print '(a,i6)','Number of salinity    data points (in depth range, for all locations) ',numsalnpoints
   print '(a,i6)','Number of temperature data points (in depth range, for all locations) ',numtemppoints

   print *

   if (numsalnpoints>0)  then
      print *,'Mean model salinity ',meanmsaln/numsalnpoints
      print *,'Model salinity bias ',biassaln/numsalnpoints
      print *,'Model salinity RMS  ',sqrt(rmssaln/numsalnpoints)
      rmssaln=sqrt(rmssaln/numsalnpoints)
      biassaln=biassaln/numsalnpoints
      meanmsaln=meanmsaln/numsalnpoints
   else
      meanmsaln=undef
      biassaln=undef
      numsalnpoints=0
      rmssaln=undef
   end if

   print *

   if (numtemppoints>0)  then
      print *,'Mean model temperature ',meanmtemp/numtemppoints
      print *,'Model temperature bias ',biastemp/numtemppoints
      print *,'Model temperature RMS  ',sqrt(rmstemp/numtemppoints)
      rmstemp=sqrt(rmstemp/numtemppoints)
      biastemp=biastemp/numtemppoints
      meanmtemp=meanmtemp/numtemppoints
   else
      meanmtemp=undef
      biastemp=undef
      numtemppoints=0
      rmstemp=undef
   end if

   ! Output regions involved to a data file
   open(10,file='argoreduce.locations',action='write',status='replace')
   do i=1,numfiles
      write(10,*) lons(i),lats(i), ptembias(i),psalbias(i)
   end do
   close(10)

   ! Output RMS values to "argoreduce.out" - post processing 
   ! should handle the rest
   open(10,file='argoreduce.out',action='write',status='replace')
   write(10,'(2i10,3f12.4,i6,3f12.4,i6)') &
              daterange(1),daterange(2), &
              meanmtemp, biastemp, rmstemp, &
              numtemppoints,&
              meanmsaln, biassaln, rmssaln, &
              numsalnpoints
   close(10)

 





end program
