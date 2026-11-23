module m_mosf_theta0
contains

subroutine mosf_theta0(u,v,tem,dp,qlon,qlat,idm,jdm,kdm,mostrf,lats,theta0array,nlats,ntheta0,undef)
use m_spherdist
use m_tecfld
implicit none

   integer, intent(in) :: idm,jdm,kdm,nlats,ntheta0
   real,    intent(in) :: qlat(0:idm+1,0:jdm+1),qlon(0:idm+1,0:jdm+1)
   real,    intent(in), dimension(idm,jdm,kdm) :: u,v,dp,tem
   real,    intent(out) :: mostrf(nlats,ntheta0)
   real,    intent(in) :: lats(nlats)
   real,    intent(in) :: theta0array(ntheta0)
   real, intent(in) :: undef
      
   real,    dimension(idm,jdm) ::  depthu,depthv,scvx,scuy, &
      depth,fumask, fvmask, ftstmask
   real,    dimension(idm,jdm,kdm)   :: dpu, dpv ,transu,transv
   real,    dimension(idm,jdm,kdm+1) :: sumdpu,sumdpv
   real,    dimension(idm,jdm,ntheta0) :: inttransu2,inttransv2, crssecu,crssecv
   integer, dimension(idm,jdm,ntheta0) :: ncnt_uthetaindex, ncnt_vthetaindex
   integer, dimension(idm,jdm,kdm) :: utheta0_index, vtheta0_index
   logical :: mask(idm,jdm)
   real :: thetadel
   real :: temu,temv
   real :: crssec(nlats,ntheta0)
   real :: trres, velres

   integer :: ideep, ilat, i, j, im1, ip1, jm1, jp1,k, itheta,ithetau,ithetav,nvalid
   logical, parameter :: masktest=.false.

   !include 'stmt_funcs.h'
   !pi = acos(0.)*2
   !radian=180/pi
   !thref=1e-3

   print *,'mosf_theta entry'
   thetadel = (maxval(theta0array)-minval(theta0array))/(ntheta0-1)
   !print *,'thetadel=',thetadel

!
! --- ------------------------------------------------------------------
! --- calculate scale factors and depths in u- and v-points
! --- ------------------------------------------------------------------
!
   print *,'mosf_theta depth, scuy,scvx'
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
   ncnt_uthetaindex=0
   ncnt_vthetaindex=0
   print *,'dpu,dpv'
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,im1,jm1temu,temv)
   do j=1,jdm
     jm1=max(j-1,1)
     do i=1,idm
       im1=max(1,i-1)

       dpu(i,j,k)   = 0.5*(dp(i,j,k) + dp(im1,j,k))
       dpv(i,j,k)   = 0.5*(dp(i,j,k) + dp(i,jm1,k))

       sumdpu(i,j,k+1) = min(depthu(i,j),sumdpu(i,j,k)  + dpu(i,j,k))
       sumdpv(i,j,k+1) = min(depthv(i,j),sumdpv(i,j,k)  + dpv(i,j,k))

       dpu(i,j,k) = sumdpu(i,j,k+1)  - sumdpu(i,j,k)
       dpv(i,j,k) = sumdpv(i,j,k+1)  - sumdpv(i,j,k)

       transu(i,j,k)=dpu(i,j,k)*u(i,j,k)*scuy(i,j)
       transv(i,j,k)=dpv(i,j,k)*v(i,j,k)*scvx(i,j)

       !temu = (dp(i,j,k)*tem(i,j,k)+tem(im1,j,k)*dp(im1,j,k)) / &
       !       (2*max(dpu(i,j,k),1e-4))
       !temv = (dp(i,j,k)*tem(i,j,k)+tem(i,jm1,k)*dp(i,jm1,k)) / &
       !       (2*max(dpv(i,j,k),1e-4))
       temu = .5*(tem(i,j,k)+tem(im1,j,k))
       temv = .5*(tem(i,j,k)+tem(i,jm1,k))

       ! Put into index array
       utheta0_index(i,j,k)= min(max(1,floor((temu-minval(theta0array))/thetadel)+1),ntheta0)
       vtheta0_index(i,j,k)= min(max(1,floor((temv-minval(theta0array))/thetadel)+1),ntheta0)
       !print *,i,j,k,temu,utheta0_index(i,j,k),sal(i,j,k)

       if ( dpu(i,j,k) <1e-4 ) then
          itheta=utheta0_index(i,j,k)
          ncnt_uthetaindex(i,j,itheta) =  ncnt_uthetaindex(i,j,itheta) +1
       end if
       if ( dpv(i,j,k) <1e-4 ) then
          itheta=vtheta0_index(i,j,k)
          ncnt_vthetaindex(i,j,itheta) =  ncnt_vthetaindex(i,j,itheta) +1
       end if

     end do
   end do
   !$OMP END PARALLEL DO
   end do


   ! Top-down
   print *,' layer transport'
   inttransu2=0.
   crssecu=0.
   crssecv=0.
   inttransv2=0.
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,ithetau,ithetav)
   do j=1,jdm
   do i=1,idm

      ithetau=utheta0_index(i,j,k)
      ithetav=vtheta0_index(i,j,k)
      
      inttransu2(i,j,ithetau) = inttransu2(i,j,ithetau)+ transu(i,j,k)
      inttransv2(i,j,ithetav) = inttransv2(i,j,ithetav)+ transv(i,j,k)

      ! Cross section area for theta layers
      crssecu(i,j,itheta) = crssecu(i,j,itheta) + scuy(i,j)*dpu(i,j,k)
      crssecv(i,j,itheta) = crssecv(i,j,itheta) + scvx(i,j)*dpv(i,j,k)

   end do
   end do
   !$OMP END PARALLEL DO
   end do

   !! The above gives the transport in each layer. Integrate.
   print *,' integrated transport'
   do itheta=ntheta0-1,1,-1
   !$OMP PARALLEL DO PRIVATE(i,j)
   do j=1,jdm
   do i=1,idm
      inttransu2(i,j,itheta) = inttransu2(i,j,itheta)+ inttransu2(i,j,itheta+1)
      inttransv2(i,j,itheta) = inttransv2(i,j,itheta)+ inttransv2(i,j,itheta+1)

      crssecu(i,j,itheta) = crssecu(i,j,itheta) + crssecu(i,j,itheta+1) 
      crssecv(i,j,itheta) = crssecv(i,j,itheta) + crssecv(i,j,itheta+1) 
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
   mostrf=0.
   crssec=0.
   !oldmask=.false.
   !print *,maxval(dpv(:,:,1)),maxval(dpu(:,:,1)),maxval(dp(:,:,1))
   mask=.false.
   do ilat=1,nlats
      !print *,ilat,nlats,lats(ilat)

      ! Find points inside this region (North of lats(i)
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


       !$OMP PARALLEL DO PRIVATE (i,j,itheta)
       do itheta=1,ntheta0

          do j=1,jdm
          do i=1,idm
          if (mask(i,j)) then
             mostrf(ilat,itheta) =  mostrf(ilat,itheta) + &
                fumask(i,j)*inttransu2(i,j,itheta)      + &
                fvmask(i,j)*inttransv2(i,j,itheta)     

             crssec(ilat,itheta) = crssec(ilat,itheta)  + &
                abs(fumask(i,j))*crssecu(i,j,itheta)    + &
                abs(fvmask(i,j))*crssecv(i,j,itheta)
          end if
          end do
          end do

          mostrf(ilat,itheta)= mostrf(ilat,itheta)*1e-6

          !nvalid= sum(abs(nint(fumask))*ncnt_uthetaindex(:,:,itheta)) + &
          !        sum(abs(nint(fvmask))*ncnt_vthetaindex(:,:,itheta))
          !if (nvalid==0) then
          !   mostrf(ilat,itheta)= undef
          !end if
          
      end do
      !$OMP END PARALLEL DO

      print *,ilat,mostrf(ilat,1),crssec(ilat,1)
      ! transport residual at "bottom"
      trres=mostrf(ilat,1)

      ! Velocity residual (1e-6 m/s)
      velres=trres/crssec(ilat,1)

      print *,ilat,trres,velres,crssec(ilat,1)

      ! Make streamfunction "nondivergent"
      do itheta=1,ntheta0
        mostrf(ilat,itheta) = mostrf(ilat,itheta) - velres*crssec(ilat,itheta)
     end do


   end do

   !print *,minval(mostrf),maxval(mostrf)

   !call tecfld('mostrf',mostrf,nlats,nds,mostrf)
   !open(10,file='mostrftheta0.tec',status="unknown")
   !write(10,*)'TITLE = "mostrf"'
   !write(10,*)'VARIABLES = "lat" "depth" "mostrf"'
   !write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',nlats,', J=',ntheta0,', K=1'
   !write(10,'(10e15.5)')   ((lats(ilat)        ,ilat=1,nlats),itheta=1,ntheta0)
   !write(10,'(10e15.5)')   ((theta0array(itheta)         ,ilat=1,nlats),itheta=1,ntheta0)
   !write(10,'(10e15.5)')((mostrf(ilat,itheta),ilat=1,nlats),itheta=1,ntheta0)
   !close(10)


end subroutine mosf_theta0
end module m_mosf_theta0


