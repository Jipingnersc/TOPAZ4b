      module m_NOR05_area_pp

      contains

      subroutine NOR05_area_pp
  
      use mod_xc
      use mod_necessary_ecovars

c --- CH: Subroutine to compute the area where you want to add up the
c ---     primary production. Called once in newhycom, as it slows 
c ---     down the code to half the  speed if called every timestep.
c ---     Based on the inbox-code from KAL, be careful when you define 
c ---     the corners, so you don't get the area outside the square...

      implicit none

      include 'common_blocks.h'

      real,dimension(4):: crnlon,crnlat
      real :: lon_i,lat_i
      real :: pp_m2
      integer :: npt
      integer :: i,j

      margin=0.
c --- CH: define area you want to add the primary production over      
      crnlon(1)= 5.3731
      crnlon(2)= -2.0583
      crnlon(3)= 1.13
      crnlon(4)= 9.843

      crnlat(1)= 62.279 
      crnlat(2)= 64.432 
      crnlat(3)= 66.639 
      crnlat(4)= 63.965 
      
      npt=4
      pp_m2=0.

      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         lat_i=plat(i,j)
         lon_i=plon(i,j)
         if(inbox(crnlon,crnlat,npt,lon_i,lat_i)) then
            if(ip(i,j)==1)then
               pp_area(i,j)=.true.
               pp_m2=pp_m2+scp2(i,j)
            else
               pp_area(i,j)=.false.
            endif
         else
            pp_area(i,j)=.false.
         endif
      enddo
      enddo
      print*,'pp_m2 is:',pp_m2
      end subroutine NOR05_area_pp
c --- Got this from Knut:)      
! Routine to calculate wether the point (plon,plat)  is in the box
! defined by crnlon,crnlat. Cnrlon/crnlat must be traversed so that 
! they form a convex polygon in 3D coordinates when following indices.
! It should work for for all regions defined by crnlon/crnlat...
      function inbox(crnlon,crnlat,npt,plon,plat)
      implicit none
      real, parameter :: rad=1.7453292519943295E-02,deg=57.29577951308232

      integer, intent(in) :: npt
      real, dimension(npt), intent(in) :: crnlon,crnlat
      real, intent(in) :: plon,plat
      logical :: inbox

      real, dimension(npt,3) :: cvec
      real, dimension(3)     :: pvec
      real, dimension(3,3)   :: rvec
      real, dimension(3)     :: nvec, nvec_prev, cprod
      integer :: i,im1,ip1
      logical :: lsign
      real    :: rotsign,old_rotsign


      ! point vector from origo
      pvec = geo2cart(plon,plat)


      ! vector to rectangle corner
      do i=1,npt
         cvec(i,:) = geo2cart(crnlon(i),crnlat(i))
      end do

      ! Traverse box boundaries -- Check that traversion is
      ! consistent and that point is in box
      lsign=.true.
      i=1
      old_rotsign=0.
      do while (i<npt+1 .and. lsign)
         
         ip1= mod(i      ,npt)+1
         im1= mod(i-2+npt,npt)+1


         ! Vectors used to span planes
         rvec(3,:) = cvec(ip1,:)
         rvec(2,:) = cvec(i  ,:)
         rvec(1,:) = cvec(im1,:)

         ! Normal vector to two spanning planes
         nvec      = cross_product(rvec(2,:),rvec(3,:))
         nvec_prev = cross_product(rvec(1,:),rvec(2,:))

         ! As we move to new planes, the cross product rotates in 
         ! a certain direction
         cprod = cross_product(nvec_prev,nvec)

         ! For anticlockwise rotation, this should be positive
         rotsign=sign(1.,dot_product(cprod,rvec(2,:)))

         ! Check that box is consistently traversed
         if (i>1 .and. rotsign * old_rotsign < 0) then
            print *,'Grid cell not consistently traversed'
            print *,'or polygon is not convex'
            stop '(inbox2)'
         end if
         old_rotsign=rotsign

         ! If this is true for all four planes, we are in grid box
         lsign = lsign .and. (dot_product(nvec,pvec)*rotsign)>0
         i=i+1

      end do


      inbox = lsign

      end function inbox
! Routine to get cartesian coordinates from geographical coordinates

      function geo2cart(lon,lat)
         real, parameter :: rad=1.7453292519943295E-02,deg=57.29577951308232
         real, intent(in) :: lon,lat
         real, dimension(3) :: geo2cart

         real :: lambda

         lambda=lat*rad
         theta=lon*rad
         geo2cart(1)=cos(lambda)*cos(theta)
         geo2cart(2)=cos(lambda)*sin(theta)
         geo2cart(3)=sin(lambda)
      end function geo2cart
! Routine to calculate cross product of two 3D vectors.
      function cross_product(v1,v2)
         implicit none
         real, intent(in), dimension(3) :: v1,v2
         real, dimension(3) :: cross_product

         cross_product(1) = v1(2)*v2(3) - v1(3)*v2(2)
         cross_product(2) = v1(3)*v2(1) - v1(1)*v2(3)
         cross_product(3) = v1(1)*v2(2) - v1(2)*v2(1)
      end function cross_product

      end module m_NOR05_area_pp
