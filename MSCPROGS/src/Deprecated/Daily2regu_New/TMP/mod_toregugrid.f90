

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! mod_toregugrid:
!
! A collection of routines and variables for converting
! the conformal mapping to a regular grid levels
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_toregugrid


   real, save :: firstlon = -999.
   real, save :: lastlon  = -999.
   real, save :: dlon     = -999.

   real, save :: firstlat = -999.
   real, save :: lastlat  = -999.
   real, save :: dlat     = -999.

   logical, save :: ini_called = .false.

   integer :: nlons
   integer :: nlats

   real, save, allocatable :: lons(:)
   real, save, allocatable :: lats(:)

   logical, save, dimension(:,:), allocatable :: ongrid
   real   , save, dimension(:,:), allocatable :: a1s, a2s, a3s, a4s
   integer, save, dimension(:,:), allocatable :: ipivs, jpivs

   logical,save :: testflag=.false.

   real, save :: undef_regu


   interface to_regugrid
      module procedure to_regugrid2d
      module procedure to_regugrid3d
   end interface



contains


subroutine regugrid_ini(flon,llon,deltalon, &
                        flat,llat,deltalat,undef)

   use mod_common
   use m_initconfmap
   use m_bilincoeff
   use m_oldtonew
   use m_pivotp
   use m_nearestpoint
   implicit none
   real, intent(in) :: flon,llon,deltalon
   real, intent(in) :: flat,llat,deltalat
   real, intent(in) :: undef

   integer :: ipiv, jpiv, ipib, jpib,ipsum,ilon,jlat
   real    :: t1,t2,a1,a2,a3,a4,lon_n, lat_n, asum
   logical :: gmsk(nx,ny),ass

   firstlon=flon;
   firstlat=flat;

   lastlon =llon;
   lastlat =llat;

   dlon =deltalon;
   dlat =deltalat;

   ! Error check
   if (dlon<0. .or. dlat<0.) then
      print *,'negative dlon or dlat!'
      stop '(regugrid_ini)'
   end if

   if (firstlon>lastlon .or. firstlat>lastlat) then
      print *,'last lon/lat less than first lon/lat!'
      stop '(regugrid_ini)'
   end if

   nlons=floor((lastlon-firstlon)/dlon) + 1
   nlats=floor((lastlat-firstlat)/dlat) + 1
   ini_called=.true.

   allocate(ongrid(nlons,nlats))
   allocate(ipivs (nlons,nlats))
   allocate(jpivs (nlons,nlats))
   allocate(a1s   (nlons,nlats))
   allocate(a2s   (nlons,nlats))
   allocate(a3s   (nlons,nlats))
   allocate(a4s   (nlons,nlats))
   allocate(lons(nlons))
   allocate(lats(nlats))

   do ilon=1,nlons
      lons(ilon)=firstlon+(ilon-1)*dlon
   end do
   do jlat=1,nlats
      lats(jlat)=firstlat+(jlat-1)*dlat
   end do

   if (.not. ini_called) then
      print *,'Initialize regugrid first!'
      stop '(to_regugrid)'
   end if

   ! Initialize conformal mapping
   call initconfmap

   do jlat=1,nlats
   do ilon=1,nlons

      ! Positions to interpolate to
      t1=lons(ilon)
      t2=lats(jlat)

      ! corresponding pivot points
      call oldtonew(t2,t1,lat_n,lon_n)
      call pivotp(lon_n,lat_n,ipiv,jpiv)

      ! Interpolation to regugrid
      if ((ipiv < 1).or.(ipiv > nx).or.(jpiv < 1).or.(jpiv > ny)) then
         ongrid(ilon,jlat)=.false.
      else
         ipib=min(ipiv+1,nx)
         JPib=min(jpiv+1,ny)
         ipsum=sum(ip(ipiv:ipiv+1,jpiv:jpiv+1))
         !if (ipsum >= 2) then
         if (ipsum >= 3) then
            call bilincoeff(modlon,modlat,nx,ny,t1,t2,ipiv,jpiv,&
                             a1,a2,a3,a4)
            
            a1=a1*ip(ipiv,jpiv)
            a2=a2*ip(ipib,jpiv)
            a3=a3*ip(ipib,jpib)
            a4=a4*ip(ipiv,jpib)

            asum = a1+a2+a3+a4
            a1=a1/asum
            a2=a2/asum
            a3=a3/asum
            a4=a4/asum

            !Call nearestpoint(modlon,modlat,nx,ny,t1,t2,ipiv,jpiv,&
            !                   A1,a2,a3,a4,gmsk,ass)

            a1s(ilon,jlat)=a1
            a2s(ilon,jlat)=a2
            a3s(ilon,jlat)=a3
            a4s(ilon,jlat)=a4
            ipivs(ilon,jlat) = ipiv
            jpivs(ilon,jlat) = jpiv
            ongrid(ilon,jlat)=.true.

         else
            ongrid(ilon,jlat)=.false.
         endif
      endif
   end do
   end do
   !$OMP END PARALLEL DO 

   undef_regu=undef
end subroutine regugrid_ini



subroutine to_regugrid2d(varin,varout)
   use mod_dimensions
   implicit none
   real, intent(in) :: varin (idm,jdm)
   real, intent(out) :: varout(nlons,nlats)

   integer :: ilon,jlat,ipiv,jpiv,jpib,ipib
   real a1,a2,a3,a4
   do jlat=1,nlats
   do ilon=1,nlons
   if (ongrid(ilon,jlat)) then
      a1=a1s(ilon,jlat)
      a2=a2s(ilon,jlat)
      a3=a3s(ilon,jlat)
      a4=a4s(ilon,jlat)
      ipiv=ipivs(ilon,jlat)
      jpiv=jpivs(ilon,jlat)
      ipib=min(ipiv+1,nx)
      JPib=min(jpiv+1,ny)


      varout(ilon,jlat) = a1*varin(ipiv,jpiv) &
                        + a2*varin(ipib,jpiv) &
                        + a3*varin(ipib,jpib) &
                        + a4*varin(ipiv,jpib)

      !if (testflag) print *,varout(ilon,jlat),varin(ipiv,jpiv)
   else
      varout(ilon,jlat)=undef_regu
   end if
   end do
   end do
end subroutine to_regugrid2d


subroutine to_regugrid3d(varin,varout)
   use mod_dimensions
   implicit none
   real, intent(in) :: varin (idm,jdm,kdm)
   real, intent(out) :: varout(nlons,nlats,kdm)


   integer :: k

   do k=1,kdm
      call to_regugrid2d(varin(:,:,k),varout(:,:,k))
   end do
end subroutine to_regugrid3d


end module mod_toregugrid
