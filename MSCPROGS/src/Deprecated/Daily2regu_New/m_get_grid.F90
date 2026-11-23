module m_get_grid 
contains

subroutine get_grid()
   use mod_common
   use m_bigrid
   implicit none

   character(len=80)  gridid
   character(len=50)  nfile
   character(len=7)   tag7

   logical ex,ex2,exuf
   integer :: i,j
   real*8, allocatable, dimension(:,:) ::  &
      qlon,qlat, &
      plon,plat, &
      ulon,ulat, &
      vlon,vlat,r8depths


   ! First, try getting from latlon.uf
   inquire(file='./Data/latlon.uf',exist=exuf)
   if (exuf) then
     write(*,*)'Load grid positions from file: latlon.uf'
     open(10,file='./Data/latlon.uf',form='unformatted')
     read(10)idm,jdm
     rewind(10)
     allocate(qlat(0:idm+1,0:jdm+1))
     allocate(qlon(0:idm+1,0:jdm+1))
     allocate(ulat(0:idm+1,0:jdm+1))
     allocate(ulon(0:idm+1,0:jdm+1))
     allocate(vlat(0:idm+1,0:jdm+1))
     allocate(vlon(0:idm+1,0:jdm+1))
     allocate(plat(0:idm+1,0:jdm+1))
     allocate(plon(0:idm+1,0:jdm+1))
     read(10)idm,jdm,qlat,qlon,plat,plon,ulat,ulon,vlat,vlon
     close(10)

   else ! Try getting from latlon.dat

      ! read lat-lon positions from grid
      inquire(file='./latlon.dat',exist=ex2)
      inquire(file='./Data/latlon.dat',exist=ex)
      if (ex) then
         nfile='./Data/latlon.dat'
        write(*,*)'Load grid positions from file: ./Data/latlon.dat'
      else if (ex2) then
         nfile='./latlon.dat'
        write(*,*)'Load grid positions from file: ./latlon.dat'
      else
         print *,'./Data/latlon.dat or ./latlon.dat does not exist'
         stop
      endif
      open(10,file=trim(nfile),form='formatted')
      READ(10,'(2i5)',ERR=250) idm,jdm
      250 continue
      allocate(qlat(0:idm+1,0:jdm+1))
      allocate(qlon(0:idm+1,0:jdm+1))
      allocate(ulat(0:idm+1,0:jdm+1))
      allocate(ulon(0:idm+1,0:jdm+1))
      allocate(vlat(0:idm+1,0:jdm+1))
      allocate(vlon(0:idm+1,0:jdm+1))
      allocate(plat(0:idm+1,0:jdm+1))
      allocate(plon(0:idm+1,0:jdm+1))
      read(10,'(a80)')gridid
      if (gridid(1:14) == 'CONFORMAL GRID') then
         print *,'Reading CONFORMAL GRID'
      else
         print *,'NBNBNB: This is not a conformal grid'
         print *,'Do you want to continue?'
         read (*,*)
      end if
      read(10,'(15e14.7)')((qlat(i,j),i=0,idm+1),j=0,jdm+1)
      read(10,'(15e14.7)')((qlon(i,j),i=0,idm+1),j=0,jdm+1)
      read(10,'(15e14.7)')((plat(i,j),i=0,idm+1),j=0,jdm+1)
      read(10,'(15e14.7)')((plon(i,j),i=0,idm+1),j=0,jdm+1)
      read(10,'(15e14.7)')((ulat(i,j),i=0,idm+1),j=0,jdm+1)
      read(10,'(15e14.7)')((ulon(i,j),i=0,idm+1),j=0,jdm+1)
      read(10,'(15e14.7)')((vlat(i,j),i=0,idm+1),j=0,jdm+1)
      read(10,'(15e14.7)')((vlon(i,j),i=0,idm+1),j=0,jdm+1)
      close(10)
   end if

   nx=idm ; ny=jdm
   ii=idm; jj=jdm ;
   ii1=ii-1; jj1=jj-1

   ! Some other fields -- get them from the ones above
   print *,'Horizontal size:',idm,jdm
   allocate(depths(idm,jdm))
   allocate(r8depths(idm,jdm))
   allocate(modlon(idm,jdm)) 
   allocate(modlat(idm,jdm))
   allocate(mqlon (0:idm+1,0:jdm+1))
   allocate(mqlat (0:idm+1,0:jdm+1))
   mqlon=qlon
   mqlat=qlat
   modlon=plon(1:idm,1:jdm)
   modlat=plat(1:idm,1:jdm)


   ! Get depth mask -- from depth file 
   write(tag7,'(i3.3,a,i3.3)')idm,'x',jdm
   inquire(file='depths'//tag7//'.uf',exist=ex)
   if (.not.ex) stop 'depths file does not exist'
   open (unit=10,file='depths'//tag7//'.uf',status='old',form='unformatted')
   read(10)r8depths
   close(10)
   depths=r8depths

   ! Allocate bigrid stuff
   allocate(ip(idm,jdm),iu(idm,jdm),iv(idm,jdm),iq(idm,jdm))
   allocate(isp(jdm),isu(jdm),isv(jdm),isq(jdm))
   allocate(jsp(idm),jsu(idm),jsv(idm),jsq(idm))
   allocate(ifp(jdm,ms),ilp(jdm,ms),ifu(jdm,ms),ilu(jdm,ms))
   allocate(ifv(jdm,ms),ilv(jdm,ms),ifq(jdm,ms),ilq(jdm,ms))
   allocate(jfp(idm,ms),jlp(idm,ms),jfu(idm,ms),jlu(idm,ms))
   allocate(jfv(idm,ms),jlv(idm,ms),jfq(idm,ms),jlq(idm,ms))

   call bigrid(depths)

   deallocate(qlat)
   deallocate(qlon)
   deallocate(ulat)
   deallocate(ulon)
   deallocate(vlat)
   deallocate(vlon)
   deallocate(plat)
   deallocate(plon)
   deallocate(r8depths)

end subroutine get_grid
end module m_get_grid
