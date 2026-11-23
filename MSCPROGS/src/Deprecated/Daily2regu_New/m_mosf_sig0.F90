module m_mosf_sig0
contains

subroutine mosf_sig0(u,v,sal,tem,dp,qlon,qlat,idm,jdm,kdm,mostrf,lats,sigarray,nlats,nsig,undef)
use m_spherdist
use m_tecfld
implicit none

   integer, intent(in) :: idm,jdm,nlats,nsig,kdm
   real,    intent(in) :: qlat(0:idm+1,0:jdm+1),qlon(0:idm+1,0:jdm+1)
   real,    intent(in), dimension(idm,jdm,kdm) :: u,v,dp,sal,tem
   real,    intent(out) :: mostrf(nlats,nsig)
   real,    intent(in) :: lats(nlats)
   real,    intent(in) :: sigarray(nsig)
   real, intent(in) :: undef
      
   real depthu(idm,jdm),depthv(idm,jdm) ,scvx(idm,jdm),scuy(idm,jdm), &
        depth(idm,jdm),fumask(idm,jdm), fvmask(idm,jdm)
   real,dimension(idm,jdm,kdm) :: dpu, dpv ,transu,transv
   real, dimension(idm,jdm) :: uint,vint,dptmp, inttransu,inttransv,ftstmask
   real, dimension(idm,jdm,kdm+1) :: sumdpu,sumdpv
   real, dimension(idm,jdm,nsig) :: inttransu2,inttransv2, crssecu,crssecv
   logical :: mask(idm,jdm)
   real :: dptmpu, dptmpv, sigdel
   real :: salu,salv,temu,temv,sigu,sigv
   real :: radian, pi, thref
   integer, dimension(idm,jdm,kdm)  :: usig0_index, vsig0_index
   integer, dimension(idm,jdm,nsig) :: ncnt_usigindex, ncnt_vsigindex
   real :: crssec(nlats,nsig)
   real :: velres, trres

   integer :: ideep, ilat, i, j, im1, ip1, jm1, jp1,k, isig,isigu,isigv,nvalid
   logical, parameter :: masktest=.false.

   include 'stmt_funcs.h'

   pi = acos(0.)*2
   radian=180/pi
   thref=1e-3

   sigdel = (maxval(sigarray)-minval(sigarray))/(nsig-1)
   !print *,'sigdel=',sigdel

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
   ncnt_usigindex=0
   ncnt_vsigindex=0
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,im1,jm1,salu,salv,temu,temv,sigu,sigv,isig)
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
       !salu = (dp(i,j,k)*sal(i,j,k)+sal(im1,j,k)*dp(im1,j,k)) / &
       !       (2*max(dpu(i,j,k),1e-4))
       !salv = (dp(i,j,k)*sal(i,j,k)+sal(i,jm1,k)*dp(i,jm1,k)) / &
       !       (2*max(dpv(i,j,k),1e-4))

       temu = .5*(tem(i,j,k)+tem(im1,j,k))
       temv = .5*(tem(i,j,k)+tem(i,jm1,k))
       salu = .5*(sal(i,j,k)+sal(im1,j,k))
       salv = .5*(sal(i,j,k)+sal(i,jm1,k))
   
       ! Densities
       sigu = sig(temu,salu)
       sigv = sig(temv,salv)

       ! Put into index array
       usig0_index(i,j,k)= max(1,min(floor((sigu-minval(sigarray))/sigdel)+1,nsig))
       vsig0_index(i,j,k)= max(1,min(floor((sigv-minval(sigarray))/sigdel)+1,nsig))
       !print *,i,j,k,temu,salu,sigu,usig0_index(i,j,k),sal(i,j,k)

       if ( dpu(i,j,k) <1e-4 ) then
          isig=usig0_index(i,j,k)
          ncnt_usigindex(i,j,isig) =  ncnt_usigindex(i,j,isig) +1
       end if
       if ( dpv(i,j,k) <1e-4 ) then
          isig=vsig0_index(i,j,k)
          ncnt_vsigindex(i,j,isig) =  ncnt_vsigindex(i,j,isig) +1
       end if

     end do
   end do
   !$OMP END PARALLEL DO
   end do
   !print *,maxval(transu),minval(transu)
   !print *,maxval(transv),minval(transv)
   !print *,maxval(usig0_index),minval(usig0_index)
   !print *,maxval(vsig0_index),minval(vsig0_index)
   !call tecfld('transu',transu(:,:,5),idm,jdm,depth)
   !call tecfld('usig0_index',float(usig0_index(:,:,5)),idm,jdm,depth)


   ! Top-down
   inttransu2=0.
   inttransv2=0.
   crssecu=0.
   crssecv=0.
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,isigu,isigv)
   do j=1,jdm
   do i=1,idm

      isigu=usig0_index(i,j,k)
      isigv=vsig0_index(i,j,k)
      
      inttransu2(i,j,isigu) = inttransu2(i,j,isigu)+ transu(i,j,k)
      inttransv2(i,j,isigv) = inttransv2(i,j,isigv)+ transv(i,j,k)

      ! Cross section area for layers
      crssecu(i,j,isigu) = crssecu(i,j,isigu) + dpu(i,j,k)*scuy(i,j)
      crssecv(i,j,isigu) = crssecv(i,j,isigv) + dpv(i,j,k)*scvx(i,j)

   end do
   end do
   !$OMP END PARALLEL DO
   end do
   !print *,'inttrans'
   !call tecfld('inttransu_1',inttransu2(:,:,5),idm,jdm,depth)
   !print *,maxval(inttransu2),minval(inttransu2)
   !print *,maxval(inttransv2),minval(inttransv2)

   ! The above gives the transport in each sigma layer. Integrate.
   do isig=2,nsig
   !$OMP PARALLEL DO PRIVATE(i,j)
   do j=1,jdm
   do i=1,idm
      inttransu2(i,j,isig) = inttransu2(i,j,isig)+ inttransu2(i,j,isig-1)
      inttransv2(i,j,isig) = inttransv2(i,j,isig)+ inttransv2(i,j,isig-1)

      crssecu(i,j,isig) = crssecu(i,j,isig) + crssecu(i,j,isig-1)
      crssecv(i,j,isig) = crssecv(i,j,isig) + crssecv(i,j,isig-1)
   end do
   end do
   !$OMP END PARALLEL DO
   end do
   !print *,'inttrans'
   !call tecfld('inttransu_2',inttransu2(:,:,5),idm,jdm,depth)
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
   do ilat=1,nlats
      !print *,ilat,nlats,lats(ilat)

      ! Find points inside this region (North of lats(i)
      mask=.false.
      !where (plat(1:idm,1:jdm)>lats(ilat).and.depth>0.)
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
      !mask = (abs(fumask)>1e-4 .or. abs(fvmask)>1e-4) !.and. depth(i,j)>1.
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



       !uint=0.
       !vint=0.
       !if (ideep==10) then
       !if (ilat==nlats/2) then
       !   call tecfld('uint',inttransu2(:,:,5),idm,jdm,depth)
       !   call tecfld('vint',inttransv2(:,:,5),idm,jdm,depth)
       !endif

       !$OMP PARALLEL DO PRIVATE(i,j,isig)
       do isig=1,nsig
          ! Integrate over latitude band
          !mostrf(ilat,isig)= sum(fumask*inttransu2(:,:,isig),mask=mask) + &
          !                   sum(fvmask*inttransv2(:,:,isig),mask=mask) 
          do j=1,jdm
          do i=1,idm
          if(mask(i,j)) then
             mostrf(ilat,isig)= mostrf(ilat,isig)+ &
             fumask(i,j)*inttransu2(i,j,isig)    + &
             fvmask(i,j)*inttransv2(i,j,isig)   

             crssec(ilat,isig) = crssec(ilat,isig) + &
                crssecu(i,j,isig)*abs(fumask(i,j)) + &
                crssecv(i,j,isig)*abs(fvmask(i,j))
          end if
          end do
          end do
             
             
          mostrf(ilat,isig)= mostrf(ilat,isig)*1e-6

          !nvalid= sum(abs(nint(fumask))*ncnt_usigindex(:,:,isig)) + &
          !        sum(abs(nint(fvmask))*ncnt_vsigindex(:,:,isig))
          !if (nvalid==0) then
          !   mostrf(ilat,isig)= undef
          !end if
             
      end do
      !$OMP END PARALLEL DO 

      ! Find transport residual
      trres=mostrf(ilat,nsig)

      ! velocity residual
      velres=trres/crssec(ilat,nsig)

      ! Make stream function "nondivergent"
      do isig=1,nsig
         mostrf(ilat,isig)= mostrf(ilat,isig) - velres*crssec(ilat,isig)
      end do

   end do

   !print *,minval(mostrf,mostrf/=undef),maxval(mostrf,mostrf/=undef)

   !call tecfld('mostrf',mostrf,nlats,nds,mostrf)
   !open(10,file='mostrfsig0.tec',status="unknown")
   !write(10,*)'TITLE = "mostrf"'
   !write(10,*)'VARIABLES = "lat" "depth" "mostrf"'
   !write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',nlats,', J=',nsig,', K=1'
   !write(10,'(10e15.5)')   ((lats(ilat)        ,ilat=1,nlats),isig=1,nsig)
   !write(10,'(10e15.5)')   ((sigarray(isig)         ,ilat=1,nlats),isig=1,nsig)
   !write(10,'(10e15.5)')((mostrf(ilat,isig),ilat=1,nlats),isig=1,nsig)
   !close(10)


end subroutine mosf_sig0
end module m_mosf_sig0


