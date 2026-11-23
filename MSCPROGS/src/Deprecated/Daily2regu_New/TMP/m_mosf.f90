

module m_mosf
contains

subroutine mosf(u,v,dp,qlon,qlat,idm,jdm,kdm,mostrf,lats,ds,nlats,nds,undef)
use m_spherdist
use m_tecfld
implicit none

   integer, intent(in) :: idm,jdm,nlats,nds,kdm
   real,    intent(in) :: qlat(0:idm+1,0:jdm+1),qlon(0:idm+1,0:jdm+1)
   real,    intent(in), dimension(idm,jdm,kdm) :: u,v,dp
   real,    intent(out) :: mostrf(nlats,nds)
   real,    intent(in) :: lats(nlats)
   real,    intent(in) :: ds(nds)
   real,    intent(in) :: undef
      
   real depthu(idm,jdm),depthv(idm,jdm) ,scvx(idm,jdm),scuy(idm,jdm), &
        depth(idm,jdm),fumask(idm,jdm), fvmask(idm,jdm)
   real,dimension(idm,jdm,kdm) :: dpu, dpv ,transu,transv
   real, dimension(idm,jdm) :: uint,vint,dptmp, inttransu,inttransv,ftstmask
   real, dimension(idm,jdm,kdm+1) :: sumdpu,sumdpv
   real, dimension(idm,jdm,nds) :: inttransu2,inttransv2,crssecu,crssecv
   logical :: mask(idm,jdm)
   real :: dptmpu, dptmpv
   real :: crssec(nlats,nds)

   integer :: ideep, ilat, i, j, im1, ip1, jm1, jp1,k 
   logical, parameter :: masktest=.false.
   real :: maxdepth, velresidual, trresidual



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

   !$OMP PARALLEL DO PRIVATE(i,j,ip1,im1,jp1,jm1)
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
     jm1=mod(j+jdm-2,jdm)+1
     do i=1,idm
       im1=max(1,i-1)
       dpu(i,j,k)   = 0.5*(dp(i,j,k) + dp(im1,j,k))
       dpv(i,j,k)   = 0.5*(dp(i,j,k) + dp(i,jm1,k))

       sumdpu(i,j,k+1)=min(sumdpu(i,j,k)+dpu(i,j,k),depthu(i,j))
       sumdpv(i,j,k+1)=min(sumdpv(i,j,k)+dpv(i,j,k),depthv(i,j))

       transu(i,j,k)=dpu(i,j,k)*u(i,j,k)*scuy(i,j)
       transv(i,j,k)=dpv(i,j,k)*v(i,j,k)*scvx(i,j)
     end do
   end do
   !$OMP END PARALLEL DO
   end do
   !print *,maxval(transu),minval(transu)
   !print *,maxval(transv),minval(transv)


   ! Top-down
   inttransu2=0.
   inttransv2=0.
   crssecu=0.
   crssecv=0.
   do ideep=1,nds
   !print *,ideep,nds,ds(ideep)
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,dptmpu,dptmpv)
   do j=1,jdm
   do i=1,idm
      dptmpu= min(sumdpu(i,j,k+1),ds(ideep)) - min(sumdpu(i,j,k),ds(ideep)) 
      dptmpv= min(sumdpv(i,j,k+1),ds(ideep)) - min(sumdpv(i,j,k),ds(ideep)) 
      inttransu2(i,j,ideep) = inttransu2(i,j,ideep)+dptmpu*transu(i,j,k)/(dpu(i,j,k)+1e-4)
      inttransv2(i,j,ideep) = inttransv2(i,j,ideep)+dptmpv*transv(i,j,k)/(dpv(i,j,k)+1e-4)


      ! For computing "cross section area" later
      crssecu(i,j,ideep)= crssecu(i,j,ideep)+dptmpu*scuy(i,j)
      crssecv(i,j,ideep)= crssecv(i,j,ideep)+dptmpv*scvx(i,j)


   end do
   end do
   !$OMP END PARALLEL DO
   end do
   end do
   !print *,maxval(inttransu2),minval(inttransu2)
   !print *,maxval(inttransv2),minval(inttransv2)


!
! --- ------------------------------------------------------------------
! --- calculate  meridional overturning stream function
! --- ------------------------------------------------------------------
!



   ! For every zonal band:
   mostrf=0.
   crssec=0.
   !oldmask=.false.
   !print *,maxval(dpv(:,:,1)),maxval(dpu(:,:,1)),maxval(dp(:,:,1))
   mask=.false.
   do ilat=1,nlats
      !print *,ilat,nlats,lats(ilat)

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
      ftstmask=0.
      !$OMP PARALLEL DO PRIVATE(i,j,ip1,jp1)
      do j=1,jdm
      do i=1,idm
      if (mask(i,j)) then

         jp1=min(jdm,j+1)
         ip1=min(idm,i+1)

         fumask(i,j)=fumask(i,j)+1.
         fvmask(i,j)=fvmask(i,j)+1.

         fumask(ip1,j)=fumask(ip1,j)-1.
         fvmask(i,jp1)=fvmask(i,jp1)-1.

         ftstmask(i,j)=1.
      end if
      end do
      end do
      !$OMP END PARALLEL DO

      ! Decrease mask to active points
      !mask = (abs(fumask)>1e-4 .or. abs(fvmask)>1e-4)
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
      !   call tecfld('fmask',ftstmask,idm,jdm,depth)
      !endif



       uint=0.
       vint=0.
       !if (ideep==nds.and.ilat==nlats/2) then
       !   call tecfld('uint',inttransu,idm,jdm,depth)
       !   call tecfld('vint',inttransv,idm,jdm,depth)
       !endif



       !$OMP PARALLEL DO PRIVATE(i,j,ideep)
       do ideep=1,nds

          !mostrf(ilat,ideep)= sum(fumask*inttransu2(:,:,ideep),mask=mask) + &
          !                    sum(fvmask*inttransv2(:,:,ideep),mask=mask) 


          ! Integrate over latitude band
          do j=1,jdm
          do i=1,idm
          if(mask(i,j)) then
             mostrf(ilat,ideep) = mostrf(ilat,ideep) + &
                   fumask(i,j)*inttransu2(i,j,ideep) + &
                   fvmask(i,j)*inttransv2(i,j,ideep)

            ! Cross section area
            crssec(ilat,ideep) = crssec(ilat,ideep) +  &
               crssecu(i,j,ideep)*abs(fumask(i,j))  +  &
               crssecv(i,j,ideep)*abs(fvmask(i,j))
          end if
          end do
          end do


          mostrf(ilat,ideep)= mostrf(ilat,ideep)*1e-6
          !if (ideep==1.and.ilat>1) then
          !   mostrf(ilat,ideep)= mostrf(ilat,ideep)+  mostrf(ilat-1,ideep)
          !end if

          !if (ds(ideep)>maxdepth) then
          !    mostrf(ilat,ideep)=undef
          !end if
      end do
      !$OMP END PARALLEL DO

      ! Bottom transport residual
      trresidual  = mostrf(ilat,nds)
      !print *,'residual ',mostrf(ilat,nds)

      ! Which translates into the following barotropic velocity (1e-6 m/s)
      velresidual = trresidual / crssec(ilat,nds)

      ! Make transport "nondivergent" -- Subtract transport residual
      do ideep=1,nds
         mostrf(ilat,ideep) = mostrf(ilat,ideep) - &
            crssec(ilat,ideep)*velresidual
      end do

      ! This should be zero now:
      !print *,'residual fixed ',mostrf(ilat,nds)

      ! Set depths below a maximum depth to undefined
      maxdepth=max( maxval(sumdpu(:,:,kdm+1),mask=abs(fumask)>1e-4), &
                    maxval(sumdpv(:,:,kdm+1),mask=abs(fvmask)>1e-4) )
      do ideep=1,nds
          if (ds(ideep)>maxdepth) then
              mostrf(ilat,ideep)=undef
          end if
      end do

   end do

   !call tecfld('mostrf',mostrf,nlats,nds,mostrf)
   !open(10,file='mostrf.tec',status="unknown")
   !write(10,*)'TITLE = "mostrf"'
   !write(10,*)'VARIABLES = "lat" "depth" "mostrf"'
   !write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',nlats,', J=',nds,', K=1'
   !write(10,'(10e15.5)')   ((lats(ilat)        ,ilat=1,nlats),ideep=1,nds)
   !write(10,'(10e15.5)')   ((ds(ideep)         ,ilat=1,nlats),ideep=1,nds)
   !write(10,'(10e15.5)')((mostrf(ilat,ideep),ilat=1,nlats),ideep=1,nds)
   !close(10)


end subroutine mosf
end module m_mosf


