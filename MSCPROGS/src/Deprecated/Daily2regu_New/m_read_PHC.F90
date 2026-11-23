module m_read_PHC
contains

   !Reads PHC data from ascii files 

   subroutine read_PHC(climid,fldname,fld,depths,reflon,dlon,reflat,dlat,maskval)
      implicit none

      integer, parameter :: cii=360,cjj=180,ckk_y=33,ckk_m=24,ctt=12
      real,    parameter :: data_reflon=.5    , data_dlon=1.
      real,    parameter :: data_reflat=-89.5 , data_dlat=1.
      real,    parameter :: mask_value=-99.

      real, dimension(ckk_y), parameter :: &
      dz=(/0,10,20,30,50,75,100,125,150,200,250,300,400,500,600,700, &
           800,900,1000,1100,1200,1300,1400,1500,1750,2000,2500,3000,&
           3500,4000,4500,5000,5500/)
 
      !real,parameter :: phcver=.98  
      real,parameter :: phcver=2.1   ! Current version from PHC


      ! Dummy variables
      character(len=*),                intent(in)  :: fldname,climid
      real, pointer,dimension(:)       :: depths
      real, pointer,dimension(:,:,:,:) :: fld
      real, intent(out) :: reflon,reflat,dlon,dlat,maskval

      ! Local variables
      character(len=80) :: fldfile
      character(len=40) :: fname,fnamegz,filebase,fileend,fileendgz
      character(len=2 ) :: cll
      logical :: ex,exgz
      real, allocatable, dimension(:) :: xvar,yvar,zvar,tvar
      real :: updiff,lwdiff,diff,data1(cii,cjj)
      integer :: l,ios,i,j,k


      if (trim(climid)=="PHC") then
         if (fldname=="temp") then
            filebase= './Climatology/PHC/Temp'
         else if (fldname=="sal") then
            filebase = './Climatology/PHC/Salt'
         else 
            print 'Unknown PHC field '//trim(fldname)
            stop '(read_phc)'
         end if
      else
         print *,'Only PHC Climatology supported'
         stop '(read_phc)'
      end if

      ! We have different PHC versions...
      if (phcver==0.98) then
         fileend='.98.p.obj'
         fileendgz='.98.p.obj.gz'
      else if (phcver==2.1) then
         fileend='_p2.1.obj'
         fileendgz='_p2.1.obj.gz'
      else
         print '(a,f6.2)','Unknown PHC version: ',phcver
         stop '(read_PHC)'
      end if






      allocate(fld(cii,cjj,ckk_y,ctt))
      allocate(depths(ckk_y))
      reflon=data_reflon
      reflat=data_reflat
      dlat  = data_dlat
      dlon  = data_dlon
      maskval = mask_value
      depths=dz

      !  First read ANNUAL data
      write(cll,'(i2.2)') 0
      fname=trim(filebase)//cll//trim(fileend)
      inquire(file=trim(fname),exist=ex) 

      fnamegz=trim(filebase)//cll//trim(fileendgz)
      inquire(file=trim(fnamegz),exist=exgz) 

      if (ex) then
         print *,'Reading Data from PHC file ',trim(fname)
      else if (exgz) then
         print *,'Unpacking ',trim(fnamegz)
         call system('gunzip '//trim(fnamegz))
      else
         print *,'Can not find ',trim(fname),' or ',trim(fnamegz)
         print *,'Make sure you have the right PHC version.'
         print *,'Current version is ',phcver
         stop '(read_PHC)'
      end if

      open(10,file=fname,form='formatted',status='old')
      do k=1,ckk_y
         read(10,1000,iostat=ios,end=100) ((data1(i,j),i=1,cii),j=1,cjj)
         fld(:,:,k,1)=data1
      end do
100   close(10)

      do l=2,ctt
         fld(:,:,:,l) = fld(:,:,:,1)
      end do



      !Read monthly products
      do l=1,ctt

         write(cll,'(i2.2)') l

         fname=trim(filebase)//cll//trim(fileend)
         inquire(file=trim(fname),exist=ex) 

         fnamegz=trim(filebase)//cll//trim(fileendgz)
         inquire(file=trim(fnamegz),exist=exgz) 


         if (ex) then
            print *,'Reading Data from PHC file ',trim(fname)
         else if (exgz) then
            print *,'Unpacking ',trim(fnamegz)
            call system('gunzip '//trim(fnamegz))
            print *,'Reading Data from PHC file ',trim(fname)
         else
            print *,'Can not find ',trim(fname),' or ',trim(fnamegz)
            print *,'Make sure you have the right PHC version.'
            print *,'Current version is ',phcver
            stop '(read_PHC)'
         end if

         open(10,file=fname,form='formatted',status='old')
         rewind(10)
         do k=1,ckk_m
            read(10,1000,iostat=ios,end=100) ((data1(i,j),i=1,cii),j=1,cjj)
            !print *,k,ios,cii,cjj
            fld(:,:,k,l)=data1
         end do
200      close(10)

         

         if (ios/=0) then 
            print *,'An eror occured when reading from '//trim(fname)
            print *,'IOSTAT was ',ios
            if (l==1) print 1000,fld(1:10,1,1,l)
            !deallocate(fld,depths)
            !stop '(read_PHC)'
         end if

         !call system('gzip '//trim(fname))

      end do


1000 format(10f8.4)

   end subroutine read_PHC
end module m_read_PHC
