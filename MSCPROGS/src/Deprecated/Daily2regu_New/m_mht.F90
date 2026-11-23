module m_mht
contains

subroutine mht(u,v,sal,tem,dp,qlon,qlat,idm,jdm,kdm,merht,lats,nlats)
use m_spherdist
use m_tecfld
implicit none

   integer, intent(in) :: idm,jdm,kdm,nlats
   real,    intent(in) :: qlat(0:idm+1,0:jdm+1),qlon(0:idm+1,0:jdm+1)
   real,    intent(in), dimension(idm,jdm,kdm) :: u,v,dp,tem,sal
   real,    intent(out):: merht(nlats)
   real,    intent(in) :: lats(nlats)
      
   real depthu(idm,jdm),depthv(idm,jdm) ,scvx(idm,jdm),scuy(idm,jdm), &
        depth(idm,jdm),fumask(idm,jdm), fvmask(idm,jdm)
   real,dimension(idm,jdm,kdm) :: dpu, dpv, temu,temv,sigu,sigv
   real, dimension(idm,jdm) :: uint,vint,dptmp
   real, dimension(idm,jdm,kdm+1) :: sumdpu,sumdpv
   logical :: mask(idm,jdm)
   real :: salu,salv

   integer :: ideep, ilat, i, j, im1, ip1, jm1, jp1,k 
   logical, parameter :: masktest=.false.
   real, parameter :: cpsw=3987.
   real :: bartr,barvel,crssec

   real :: pi, radian, thref ! needed in stmt_funcs

   include 'stmt_funcs.h'

   pi = acos(0.)*2
   radian=180./pi
   thref=1e-3

!
! --- ------------------------------------------------------------------
! --- calculate scale factors and depths in u- and v-points
! --- ------------------------------------------------------------------
!
   !$OMP PARALLEL DO PRIVATE(i,j)
   do j=1,jdm
   do i=1,idm
      depth(i,j) = sum(dp(i,j,:))
   end do
   end do
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO PRIVATE(i,j,im1,ip1,jm1,jp1)
   do j=1,jdm
     jp1=j+1
     jm1=mod(j+jdm-2,jdm)+1
     do i=1,idm
       im1=max(1,i-1)
       ip1=i+1
       scvx(i,j)=spherdist(qlon(ip1,j),qlat(ip1,j), qlon(i  ,j),qlat(i  ,j))
       scuy(i,j)=spherdist(qlon(i,jp1),qlat(i,jp1), qlon(i,j  ),qlat(i,j  ))
       depthu(i,j)=min(depth(i,j),depth(im1,j))
       depthv(i,j)=min(depth(i,j),depth(i,jm1))
       if (depth(i,j).lt.1.) depth(i,j)=0.
       if (depthu(i,j).lt.1.) depthu(i,j)=0.
       if (depthv(i,j).lt.1.) depthv(i,j)=0.
     enddo
   enddo
   !$OMP END PARALLEL DO


!
! --- ------------------------------------------------------------------
! --- calculate dp in u and v points
! --- ------------------------------------------------------------------
!
   sumdpu=0.
   sumdpv=0.
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,im1,jm1)
   do j=1,jdm
     jm1=max(1,j-1)
     do i=1,idm
       im1=max(1,i-1)

       dpu(i,j,k)   = 0.5*(dp(i,j,k) + dp(im1,j,k))
       dpv(i,j,k)   = 0.5*(dp(i,j,k) + dp(i,jm1,k))

       sumdpu(i,j,k+1) = min(depthu(i,j),sumdpu(i,j,k)  + dpu(i,j,k))
       sumdpv(i,j,k+1) = min(depthv(i,j),sumdpv(i,j,k)  + dpv(i,j,k))

       dpu(i,j,k) = sumdpu(i,j,k+1)  - sumdpu(i,j,k)
       dpv(i,j,k) = sumdpv(i,j,k+1)  - sumdpv(i,j,k)
       
     end do
   end do
   !$OMP END PARALLEL DO
   end do
      
   do k=1,kdm
      !$OMP PARALLEL DO PRIVATE(i,j,im1,jm1,salu,salv)
      do j=1,jdm
        jm1=max(1,j-1)
        do i=1,idm
          im1=max(1,i-1)
          !temu(i,j,k) = (dp(i,j,k)*tem(i,j,k) + dp(im1,j,k)*tem(im1,j,k))/ &
          !              (2*max(1e-4,dpu(i,j,k)))
          !temv(i,j,k) = (dp(i,j,k)*tem(i,j,k) + dp(i,jm1,k)*tem(i,jm1,k))/ &
          !              (2*max(1e-4,dpv(i,j,k)))
          !salu        = (dp(i,j,k)*sal(i,j,k) + dp(im1,j,k)*sal(im1,j,k))/ &
          !              (2*max(1e-4,dpu(i,j,k)))
          !salv        = (dp(i,j,k)*sal(i,j,k) + dp(i,jm1,k)*sal(i,jm1,k))/ &
          !              (2*max(1e-4,dpv(i,j,k)))
          temu(i,j,k) = .5*(tem(i,j,k) + tem(im1,j,k))
          temv(i,j,k) = .5*(tem(i,j,k) + tem(i,jm1,k))
          salu        = .5*(sal(i,j,k) + sal(im1,j,k))
          salv        = .5*(sal(i,j,k) + sal(i,jm1,k))

          sigu(i,j,k) = sig(temu(i,j,k),salu)+1000.
          sigv(i,j,k) = sig(temv(i,j,k),salv)+1000.
        end do
      end do
      !$OMP END PARALLEL DO
   end do

!
! --- ------------------------------------------------------------------
! --- calculate  meridional overturning stream function
! --- ------------------------------------------------------------------
!



   ! For every zonal band:
   merht=0.
   do ilat=1,nlats

      ! Find points inside this region (North of lats(i)
      !mask=.false.
      !where (qlat(1:idm,1:jdm)>lats(ilat))
      !   mask=.true.
      !endwhere
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         mask(i,j)=qlat(i,j)>lats(ilat)
      end do
      end do
      !$OMP END PARALLEL DO


      ! This will only retain transport at the boundary
      ! of the region
      fumask=0.
      fvmask=0.
      !$OMP PARALLEL DO PRIVATE(i,j,ip1,jp1)
      do j=1,jdm
      do i=1,idm
      if (mask(i,j)) then

         jp1=min(jdm,j+1)
         ip1=min(idm,i+1)

         fumask(i,j)=fumask(i,j)+1.
         fvmask(i,j)=fvmask(i,j)+1.

         fumask(ip1,j)=fumask(ip1,j)-1.
         fvmask(i,jp1)=fumask(i,jp1)-1.
      end if
      end do
      end do
      !$OMP END PARALLEL DO

      ! Drop points over land
      !where (depth<10.0)
      !   fumask=0.
      !   fvmask=0.
      !endwhere

      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         if (depth(i,j)<1.) then
            fumask(i,j)=0.
            fvmask(i,j)=0.
         end if
         mask(i,j) = abs(fumask(i,j))>1e-4 .or. abs(fvmask(i,j))>1e-4
      end do
      end do
      !$OMP END PARALLEL DO

      !if (ilat==nlats/2.and.masktest) then
      !   call tecfld('fumask',fumask,idm,jdm,depth)
      !   call tecfld('fvmask',fvmask,idm,jdm,depth)
      !endif

      ! We should now have transport masks
      ! for the transport into the region
      ! north of lats(i). Let the streamfunction
      ! be zero at the surface

      ! Integrate over latitude band
      crssec = 0.
      bartr=0.
      barvel=0.

      ! Cross section area and barotropic velocity for transport
      do k=1, kdm
         !$OMP PARALLEL DO PRIVATE(i,j)  &
         !$OMP REDUCTION(+:crssec,bartr)
         do j=1,jdm
         do i=1,idm
         if (mask(i,j)) then

            ! cross section
            crssec= crssec + &
            abs(fumask(i,j))*dpu(i,j,k)*scuy(i,j) + &
            abs(fvmask(i,j))*dpv(i,j,k)*scvx(i,j)

            ! Section transport
            bartr= bartr + &
            fumask(i,j)*u(i,j,k)*dpu(i,j,k)*scuy(i,j) + &
            fvmask(i,j)*v(i,j,k)*dpv(i,j,k)*scvx(i,j)

         end if
         end do
         end do
         !$OMP END PARALLEL DO
      end do

      ! Section velocity component
      barvel=bartr/crssec
      !print *,ilat,barvel,bartr,crssec


      do k=1, kdm
         !$OMP PARALLEL DO PRIVATE(i,j)
         do j=1,jdm
         do i=1,idm
         if (mask(i,j)) then

            ! heat transport
            merht(ilat)= merht(ilat) + &
            fumask(i,j)*(u(i,j,k)-barvel)*temu(i,j,k)*dpu(i,j,k)*sigu(i,j,k)*scuy(i,j) + &
            fvmask(i,j)*(v(i,j,k)-barvel)*temv(i,j,k)*dpv(i,j,k)*sigv(i,j,k)*scvx(i,j)

         end if
         end do
         end do
         !$OMP END PARALLEL DO
      end do

      merht(ilat)=merht(ilat)*cpsw
   end do

end subroutine mht
end module m_mht


