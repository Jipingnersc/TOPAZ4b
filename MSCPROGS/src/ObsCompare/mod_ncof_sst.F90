module mod_ncof_sst
   integer :: linter=0
   integer,save :: Nlon,Nlat
   integer, allocatable,save :: Ipv(:,:),Jpv(:,:)
   
contains

   ! reading the observation grid information in the SST nc file
   subroutine ncof_sst_info(Fnc,mlon,mlat,idm,jdm)
   use netcdf
   use mod_confmap
   use mod_grid
   implicit none
   integer,parameter       :: Ti=100,Tj=400
   character(*),intent(in) :: Fnc 
   integer,intent(in) :: idm,jdm
   real,   intent(in),dimension(idm,jdm) :: mlon,mlat
   !------------------------------------------
   integer :: ncid,lonid,latid,varid
   real,allocatable :: Olon(:),Olat(:)
   !------------------------------------------
   integer :: i,j,ip0,jp0
   real    :: lon_n,lat_n
   !------------------------------------------

   print *,trim(Fnc) 
   call ncerr(nf90_open(trim(Fnc),nf90_nowrite,ncid)) 
   ! Get dimensions in netcdf file
   call ncerr(nf90_inq_dimid(ncid,"lon",lonid))
   call ncerr(nf90_inquire_dimension(ncid,lonid,len=Nlon))
   call ncerr(nf90_inq_dimid(ncid,"lat",latid))
   call ncerr(nf90_inquire_dimension(ncid,latid,len=Nlat))
   if (Nlon<1.or.Nlat<1) then
      print *, "There are some errors: ",trim(Fnc)
      return
   endif
   print *,Nlon,Nlat 
   ! Get variables in netcdf file
   if (linter.eq.0) then
     allocate(Olon(Nlon),Olat(Nlat),Ipv(Nlon,Nlat),Jpv(Nlon,Nlat))     
     call ncerr(nf90_inq_varid(ncid,"lon",varid))
     call ncerr(nf90_get_var(ncid,varid,Olon))
     call ncerr(nf90_inq_varid(ncid,"lat",varid))
     call ncerr(nf90_get_var(ncid,varid,Olat))
     print *, Olon(Ti),Olat(Tj)
     do i=1,Nlon
       do j=1,Nlat
         call oldtonew(Olat(j),Olon(i),lat_n,lon_n)
         call pivotp(lon_n,lat_n,ip0,jp0)
         if (ip0>1.and.jp0>1.and.ip0<idm.and.jp0<jdm) then
           Ipv(i,j)=ip0;  Jpv(i,j)=jp0; 
         end if
       end do
     end do
     deallocate(Olon,Olat)
     linter=linter+1
   endif
   end subroutine


   ! reading the SST nc file,and then interpoating on the model grid
   subroutine ncof_sst_read(Fnc,Nsst,Fice,Ersst,idm,jdm,depth)
   use netcdf
   use mod_confmap
   use mod_grid
   implicit none
   integer,parameter       :: Ti=100,Tj=400
   character(*),intent(in) :: Fnc 
   integer,intent(in) :: idm,jdm
   real,intent(in) :: depth(idm,jdm)
   real,  intent(out),dimension(idm,jdm) :: Nsst,Fice,Ersst 
   !------------------------------------------
   integer :: ncid,varid,varid1
   integer,dimension(nf90_max_var_dims) :: dimIDs
   real,  allocatable :: Var(:,:),Var2(:,:)
   integer*1,allocatable :: Var1(:,:)
   !------------------------------------------
   integer :: i,j,ip0,jp0
   !------------------------------------------

   print *,'linter=', linter
   call ncerr(nf90_open(trim(Fnc),nf90_nowrite,ncid))
   ! Get dimensions in netcdf file
   print *,Nlon,Nlat 
   ! Get variables in netcdf file
   allocate(Var(Nlon,Nlat),Var1(Nlon,Nlat),Var2(Nlon,Nlat));
   call ncerr(nf90_inq_varid(ncid,"analysed_sst",varid))
   call ncerr(nf90_inquire_variable(ncid,varid,dimids=dimIDs))
   print *, dimIDs(1),dimIDs(2),varid
   call ncerr(nf90_get_var(ncid,varid,Var))

   call ncerr(nf90_inq_varid(ncid,"analysis_error",varid))
   call ncerr(nf90_inquire_variable(ncid,varid,dimids=dimIDs))
   print *, dimIDs(1),dimIDs(2),varid
   call ncerr(nf90_get_var(ncid,varid,Var2))

   call ncerr(nf90_inq_varid(ncid,"sea_ice_fraction",varid1))
   print *, dimIDs(1),dimIDs(2),varid1
   call ncerr(nf90_get_var(ncid,varid1,Var1))
   print *, Var(Ti,Tj), Var(Tj,Ti) 
   ! Get new variables

   do j=1,Nlat
     do i=1,Nlon
       ip0=Ipv(i,j);  jp0=Jpv(i,j);
       if (ip0>0.and.jp0>0.and.ip0<idm.and.jp0<jdm) then
         if (abs(depth(ip0,jp0))>0.and.abs(depth(ip0,jp0))<20000) then
           Nsst(ip0,jp0)=Var(i,j)*0.01
           Ersst(ip0,jp0)=Var2(i,j)*0.01
           Fice(ip0,jp0)=Var1(i,j)*0.01
           if (abs(Nsst(ip0,jp0))>100) then
             Nsst(ip0,jp0)=-32768
             Ersst(ip0,jp0)=-32768
             Fice(ip0,jp0)=-32768
           endif
         else
           Nsst(ip0,jp0)=-32768
           Ersst(ip0,jp0)=-32768
           Fice(ip0,jp0)=-32768
         endif
       endif
     end do
   end do
   call ncerr(nf90_close(ncid))
   deallocate(Var,Var1)
   end subroutine

   subroutine ncof_sst_dump(Newnc,idm,jdm,mlon,mlat,Nsst0,Fice0,Ersst0)
   use netcdf
   implicit none
   character(*),intent(in) :: Newnc
   integer,     intent(in) :: idm,jdm
   real,   intent(in),dimension(idm,jdm) :: mlon,mlat,Nsst0,Fice0,Ersst0
   !---------------------
   integer :: rdimid, jdimid, idimid
   integer :: ncid,varid1,varid2,varid3,varid4,varid5
   integer :: dimms2d(2)
   !---------------------
   real*4, allocatable,dimension(:,:) :: Nsst,Fice,Ersst
   allocate(Nsst(idm,jdm),Fice(idm,jdm),Ersst(idm,jdm))
   Nsst=Nsst0
   Fice=Fice0
   Ersst=Ersst0
   
   print *,trim(Newnc)
   call ncerr(nf90_create(trim(Newnc),nf90_noclobber,ncid))
   call ncerr(nf90_def_dim(ncid,"idim",idm,idimid))
   call ncerr(nf90_def_dim(ncid,"jdim",jdm,jdimid))
   print *,idm
   dimms2d=(/idimid,jdimid/) 
   !call ncerr(nf90_redef(ncid))
   call ncerr(nf90_def_var(ncid,'lon',nf90_float,dimms2d,varid1))
   call ncerr(nf90_def_var(ncid,'lat',nf90_float,dimms2d,varid2))
   call ncerr(nf90_def_var(ncid,'ncof_sst',nf90_float,dimms2d,varid3))
   call ncerr(nf90_put_att(ncid,varid3,'_fillvalue',-32768))
   call ncerr(nf90_put_att(ncid,varid3,'missing_value',-32768))
   call ncerr(nf90_def_var(ncid,'ncof_fice',nf90_float,dimms2d,varid4))
   call ncerr(nf90_put_att(ncid,varid4,'_fillvalue',-32768))
   call ncerr(nf90_put_att(ncid,varid4,'missing_value',-32768))
   call ncerr(nf90_def_var(ncid,'sst_error',nf90_float,dimms2d,varid5))
   call ncerr(nf90_put_att(ncid,varid5,'_fillvalue',-32768))
   call ncerr(nf90_put_att(ncid,varid5,'missing_value',-32768))
   call ncerr(nf90_enddef(ncid))
   call ncerr(nf90_put_var(ncid,varid1,mlon))
   call ncerr(nf90_put_var(ncid,varid2,mlat))
   call ncerr(nf90_put_var(ncid,varid3,Nsst))
   call ncerr(nf90_put_var(ncid,varid4,Fice))
   call ncerr(nf90_put_var(ncid,varid5,Ersst))
   call ncerr(nf90_close(ncid))
   deallocate(Fice,Ersst,Nsst) 
   end subroutine

   ! Give netcdf error message
   subroutine ncerr(error)
      use netcdf
      implicit none
      integer, intent(in) :: error
      if (error/=nf90_noerr) then
         print *,'ncfile error '//nf90_strerror(error)
         stop
      end if
   end subroutine ncerr


end module 
