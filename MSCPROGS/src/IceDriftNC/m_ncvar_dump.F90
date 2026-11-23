module m_ncvar_dump


contains

subroutine ncvar_dump2D(nx,ny,mdx,mdy,imem)
   use netcdf
   use m_handle_err
   implicit none
   integer,                intent(in) :: nx,ny,imem
   real, dimension(:,:), intent(in) :: mdx,mdy
   character(len=80) :: filenc

   real*4, parameter :: undefr4=-1.e+10
   integer :: dimms3d(3),dimms2d(2)
   integer :: rdimid, jdimid, idimid,varid,ncid
   character*3 tmpX


   logical :: ex1,status0
   character(len=8)   :: dateinfo

   integer, dimension(nx,ny) :: mask0
   real,    dimension(nx,ny) :: tmpdx,tmpdy
   
   mask0=0
   where(mdx==0.and.mdy==0)  mask0=1
   tmpdx=mdx
   tmpdy=mdy
   where(mask0==1) tmpdx=undefr4   
   where(mask0==1) tmpdy=undefr4   

   write(tmpX,'(i3.3)') imem
   filenc='tmp'//trim(tmpX)//'.nc' 

   inquire(file=trim(filenc),exist=ex1)
   if (.not.ex1) then
      ! This also creates the netcdf file
      if (NF90_CREATE(trim(filenc),NF90_CLOBBER,ncid) /= NF90_NOERR) then
         print *,'An error occured when opening a netcdf file'
         stop '(ncvar_dump2D)'
      end if
     print *,'idm=',nx,'jdm=',ny,'rdm=',nf90_unlimited
     call handle_err(NF90_DEF_DIM(ncid,'idim',nx,idimid))
     call handle_err(NF90_DEF_DIM(ncid,'jdim',ny,jdimid))
     call handle_err(NF90_DEF_DIM(ncid,'rdim',nf90_unlimited,rdimid))
 
     call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'title','from daily-nc generated ice drifting fields marched with the OSISAF drift nc-file!'))
     call date_and_time(date=dateinfo)
     call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Dataset_generation_date', &
         dateinfo(1:8)))
     call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'history', &
         'Created '//dateinfo(1:8)//' by program icedrift'))
     call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Conventions','-'))


   else 
     if (NF90_OPEN(trim(filenc),NF90_WRITE,ncid) /= NF90_NOERR) then
       print *,"Error to open "//trim(filenc)
       stop '(ncvar_dump2D)'
     end if
     call handle_err (nf90_inq_dimid(ncid, 'idim', idimid))
     call handle_err (nf90_inq_dimid(ncid, 'jdim', jdimid))
     call handle_err (nf90_inq_dimid(ncid, 'rdim', rdimid))
     call handle_err(NF90_REDEF(ncid))
   endif

   print *,"open "//trim(filenc)
   dimms2d=(/idimid,jdimid/)
   dimms3d=(/idimid,jdimid,rdimid/)

   call handle_err(NF90_DEF_VAR(ncid,'dX-ice',NF90_Float,dimms2d,varid))
   call handle_err(NF90_PUT_ATT(ncid,varid,'units',"km"))
   call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',undefr4))
   call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
   call handle_err(NF90_PUT_ATT(ncid,varid,'long_name',"component of the displacement along the x axis of the grid"))
   call handle_err(NF90_ENDDEF(ncid))
   call handle_err(NF90_PUT_VAR(ncid,varid,tmpdx,start=(/1,1/)))

   call handle_err(NF90_REDEF(ncid))
   call handle_err(NF90_DEF_VAR(ncid,'dY-ice',NF90_Float,dimms2d,varid))
   call handle_err(NF90_PUT_ATT(ncid,varid,'units',"km"))
   call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',undefr4))
   call handle_err(NF90_PUT_ATT(ncid,varid,'missing_value',undefr4))
   call handle_err(NF90_PUT_ATT(ncid,varid,'long_name',"component of the displacement along the y axis of the grid"))
   call handle_err(NF90_ENDDEF(ncid))
   call handle_err(NF90_PUT_VAR(ncid,varid,tmpdy,start=(/1,1/)))


   call handle_err(NF90_CLOSE(ncid))


end subroutine ncvar_dump2D



end module m_ncvar_dump
