module mod_spline_calc
   !integer, parameter :: ndeep=12
   !real, parameter :: Deeps(1:ndeep)=(/ 5.0, 30.0, 50.0,  100.0, 200.0, 400.0, 700.0, 1000.0, &
   !                    1500.0, 2000.0, 2500.0, 3000.0 /)
   ! 
   !
   integer, save :: ndeep
   real, save, allocatable :: deeps(:)
   real, save :: undef_spline=-1e9

   integer, dimension(:,:,:), allocatable, save :: kindex
   logical,save :: splineflag

contains

subroutine spline_calc_ini(z,nzlev,undef,lsplinein)
   !use m_initconfmap
   implicit none
   integer, intent(in) :: nzlev
   real, intent(in) :: z(nzlev)
   real, intent(in) :: undef
   logical, intent(in) :: lsplinein

   allocate(deeps(nzlev))
   ndeep=nzlev
   deeps(:)=z(:)
   undef_spline=undef

   splineflag=lsplinein

end subroutine

subroutine spline_calc_release()
   deallocate(deeps)
end subroutine


subroutine spline_calc(fldin,depthin,nlons,nlats,mask,fldout,nzlev)
use mod_common
use m_initconfmap
implicit none
integer, intent(in) :: nlons,nlats,nzlev
real, dimension(nlons,nlats,kdm)  ,intent(in) :: fldin,depthin
real, dimension(nlons,nlats,nzlev),intent(out) :: fldout !,depthout
logical, intent(in) :: mask(nlons,nlats)

integer :: ilon,jlat,sub_last,k,klevel,kzlevel,klevel2
real, dimension(kdm) :: tmpfld,tmpd
real, dimension(kdm) :: subfld,subd
real, dimension(ndeep) :: sfld
real :: maxdin,stepk



!print *,'enter spline_calc'

if (nzlev/=ndeep) then
   print *,'Dimension mismatch in'
   print *,nzlev,ndeep
   stop '(spline_calc)'
end if

if  (splineflag) then
!print *,'loop begins'
do jlat=1,nlats
!if (mod(jlat,nlats/10)==0) print *,jlat,nlats
if (mod(jlat,nlats/10)==0) write(6,'(a1)',advance='no') '.'
do ilon=1,nlons
!print *,'beg',jlat,ilon,nlats,nlons,mask(ilon,jlat)
if(mask(ilon,jlat) .and. depthin(ilon,jlat,2)>0 ) then

      tmpfld=fldin  (ilon,jlat,:)
      tmpd  =depthin(ilon,jlat,:)
      maxdin=maxval(depthin(ilon,jlat,:))
      !print *,'tmpd  ',tmpd
      !print *,'tmpfld',tmpfld

      call ini_spline(tmpd,tmpfld,kdm,subd,subfld,sub_last)
      !print *,'subd  ',subd
      !print *,'subfld',subfld
      !print *,'sub_last',sub_last

      call compute_spline(subd,subfld,kdm,sub_last,sfld,deeps,ndeep,undef_spline)
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
   !print *,'end',jlat,ilon,nlats,nlons
enddo
enddo

else

   !print *,'splineflag',splineflag
   fldout=undef_spline

   ! Fast option - use layer value, can be improved a lot further
   if (.not.allocated(kindex)) then
      allocate(kindex(nlons,nlats,ndeep))
      kindex=-1
      !print *,'splineflag',splineflag

      do kzlevel=1,ndeep
      do klevel=1,kdm
      do jlat=1,nlats
      do ilon=1,nlons
         if (klevel==1) then
            if (depthin(ilon,jlat,klevel)>deeps(kzlevel)) then
               fldout(ilon,jlat,kzlevel)=fldin(ilon,jlat,klevel)
               kindex(ilon,jlat,kzlevel)=klevel
            end if
         else 
            if (depthin(ilon,jlat,klevel  )> deeps(kzlevel) .and. &
                depthin(ilon,jlat,klevel-1)<=deeps(kzlevel)) then
               fldout(ilon,jlat,kzlevel)=fldin(ilon,jlat,klevel)
               kindex(ilon,jlat,kzlevel)=klevel
            end if
         end if
      end do
      end do
      end do
      end do
      !print *,'End First kindex loop'

   else

      !print *,'2nd kindex loop'

      do kzlevel=1,ndeep
      do jlat=1,nlats
      do ilon=1,nlons
         klevel=kindex(ilon,jlat,kzlevel)
         klevel2=max(klevel,1)
         stepk=(1.+sign(1.,float(klevel)))*.5


         fldout(ilon,jlat,kzlevel)= &
            stepk*fldin(ilon,jlat,klevel2) + (1.-stepk)*undef_spline
         !if (klevel/=-1) then
         !   fldout(ilon,jlat,kzlevel)=fldin(ilon,jlat,klevel)
         !else
         !   fldout(ilon,jlat,kzlevel)=undef_spline
         !end if
      end do
      end do
      end do
   end if
end if



   


end subroutine spline_calc

subroutine spline_calc_1d(fldin,depthin,fldout,nzlev)
use mod_common
use m_initconfmap
implicit none
integer, intent(in) :: nzlev
real, dimension(kdm)  ,intent(in) :: fldin,depthin
real, dimension(nzlev),intent(out) :: fldout !,depthout

integer :: ilon,jlat,sub_last,k
real, dimension(kdm) :: tmpfld,tmpd
real, dimension(kdm) :: subfld,subd
real, dimension(ndeep) :: sfld
real :: maxdin


if (nzlev/=ndeep) then
   print *,'Dimension mismatch in'
   print *,nzlev,ndeep
   stop '(spline_calc)'
end if


tmpfld=fldin  
tmpd  =depthin
maxdin=maxval(depthin)

call ini_spline(tmpd,tmpfld,kdm,subd,subfld,sub_last)
!print *,'subd  ',subd
!print *,'subfld',subfld
!print *,'sub_last',sub_last

call compute_spline(subd,subfld,kdm,sub_last,sfld,deeps,ndeep,undef_spline)
!print *,'deep  ',deeps
!print *,'sfld  ',sfld

fldout=sfld
!do k=1,ndeep
!   depthout(ilon,jlat,k)=deeps(k)
!   depthout(ilon,jlat,k)=min(depthout(ilon,jlat,k),maxdin)
!end do
end subroutine spline_calc_1d


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
         !print *,'diff',(tmpd(k)-tmpd(k-1))
         if ((tmpd(k)-tmpd(k-1)) > 1.0) then
            l=l+1
            ip(l)=k
         endif
      enddo
      nip=l
      !print *,'nip,nz',nip,nz

      ! Put into array
      l=0
      do k=1,nip
         l=l+1
         tmpip=ip(k) 
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
