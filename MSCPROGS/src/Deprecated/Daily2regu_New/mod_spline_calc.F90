module mod_spline_calc
   !integer, parameter :: ndeep=12
   !real, parameter :: Deeps(1:ndeep)=(/ 5.0, 30.0, 50.0,  100.0, 200.0, 400.0, 700.0, 1000.0, &
   !                    1500.0, 2000.0, 2500.0, 3000.0 /)
   ! 
   !
   integer, save :: ndeep
   real, save, allocatable :: deeps(:)
   real, save :: undef_spline=-1e9


contains

subroutine spline_calc_ini(z,nzlev,undef)
   !use m_initconfmap
   implicit none
   integer, intent(in) :: nzlev
   real, intent(in) :: z(nzlev)
   real, intent(in) :: undef

   allocate(deeps(nzlev))
   ndeep=nzlev
   deeps(:)=z(:)
   undef_spline=undef
   !print *,deeps

   !call initconfmap
end subroutine


subroutine spline_calc(fldin,depthin,nlons,nlats,mask,fldout,nzlev)
use mod_dimensions
use m_initconfmap
implicit none
integer, intent(in) :: nlons,nlats,nzlev
real, dimension(nlons,nlats,kdm)  ,intent(in) :: fldin,depthin
real, dimension(nlons,nlats,nzlev),intent(out) :: fldout !,depthout
logical, intent(in) :: mask(nlons,nlats)

integer :: ilon,jlat,sub_last,k
real, dimension(kdm) :: tmpfld,tmpd
real, dimension(kdm) :: subfld,subd
real, dimension(ndeep) :: sfld
real :: maxdin


!print *,'enter spline_calc'

if (nzlev/=ndeep) then
   print *,'Dimension mismatch in'
   print *,nzlev,ndeep
   stop '(spline_calc)'
end if

!print *,'loop begins'
do jlat=1,nlats
!if (mod(jlat,nlats/10)==0) print *,jlat,nlats
if (mod(jlat,nlats/10)==0) write(6,'(a1)',advance='no') '.'
do ilon=1,nlons
if(mask(ilon,jlat) .and. depthin(ilon,jlat,2)>0 ) then

      tmpfld=fldin  (ilon,jlat,:)
      tmpd  =depthin(ilon,jlat,:)
      maxdin=maxval(depthin(ilon,jlat,:))
      !print *,'tmpd  ',tmpd
      !print *,'tmpfld',tmpfld

      call ini_spline(tmpd,tmpfld,nz,subd,subfld,sub_last)
      !print *,'subd  ',subd
      !print *,'subfld',subfld

      call compute_spline(subd,subfld,nz,sub_last,sfld,deeps,ndeep,undef_spline)
      !print *,'deep  ',deeps
      !print *,'sfld  ',sfld

      fldout(ilon,jlat,:)=sfld
      !do k=1,ndeep
      !   depthout(ilon,jlat,k)=deeps(k)
      !   depthout(ilon,jlat,k)=min(depthout(ilon,jlat,k),maxdin)
      !end do

   else
      fldout(ilon,jlat,:)=undef_spline
   endif
enddo
enddO


end subroutine spline_calc



subroutine ini_spline(tmpd,tmpfld,nz,subd,subfld,sub_last)
! setting up vectors to interpolate from.
! neglecting thin layers
! For spline version

      implicit none
      integer, intent(in)  :: nz
      real, intent(in) , dimension(nz) :: tmpd, tmpfld
      real, intent(out), dimension(nz) :: subd, subfld
      integer, intent(out) :: sub_last

      integer ip(nz)
      integer l,k,nip, tmpip

      ! nip is lowest mass-containing layer
      ip(1)=1
      l=1
      do k=2,nz
         if ((tmpd(k)-tmpd(k-1)) > 1.0) then
            l=l+1
            ip(l)=k
         endif
      enddo
      nip=l

      ! Put into array
      l=0
      do k=1,nip
         l=l+1
         tmpip=ip(k) ! Vectors into arrays are bad...
         subfld(l)=tmpfld(tmpip)
         subd  (l)=tmpd  (tmpip)
      enddo
      sub_last=l

end subroutine ini_spline

subroutine compute_spline(subd,subfld,nz,sub_last,sfld,deeps2,ndepths,undef_spline)
   use m_calcu
   use m_calcconst
   implicit none
   integer, intent(in) :: sub_last,ndepths,nz
   real, intent(in), dimension(nz) :: subd, subfld
   real, intent(out),  dimension(ndepths) :: sfld
   real, intent(in),  dimension(ndepths) :: deeps2
   real, intent(in) :: undef_spline

   real, dimension(sub_last):: a, b
   real, dimension(0:sub_last):: c
   integer ic
   integer k

   call calcconst(sub_last,subfld,subd,a,b,c,'A','a')
   sfld=undef_spline
   do k=1,ndepths
      sfld(k)=calcu(sub_last,a,b,c,subd,deeps2(k),ic,undef_spline)
   enddo

end subroutine compute_spline


end module mod_spline_calc
