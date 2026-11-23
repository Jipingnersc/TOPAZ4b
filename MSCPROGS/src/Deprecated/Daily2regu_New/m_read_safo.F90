module m_read_safo
private :: errcond, getdimms,getvarstat
contains

   !Reads SOFA data (SST) from a NetCDF file

   subroutine read_safo(filename,fld,longitude,latitude,undef)
      use netcdf
      implicit none

      ! Dummy variables
      character(len=*),                intent(in)  :: filename
      real, pointer,   dimension(:)   :: longitude,latitude
      real, pointer,   dimension(:,:) :: fld
      real, intent(in) :: undef

      ! Local variables
      integer, dimension(nf90_max_dims) :: lonvdims,latvdims,sstvdims
      integer :: londimid,londimlen,lonvarid,lonvndim
      integer :: latdimid,latdimlen,latvarid,latvndim
      integer :: sstvarid,sstvndim
      integer :: ncid
      integer :: tvarid,fldndim,fldvarid,i,stat
      logical :: ex
      real, allocatable, dimension(:) :: xvar,yvar,zvar,tvar
      integer, allocatable, dimension(:,:) :: tmpfld
      real :: updiff,lwdiff,diff
      real :: add_offset,scale_factor, fillvalue


      ! Check for existence of netcdf file
      inquire(file=trim(filename),exist=ex)
      if (.not.ex) then
         print *,'Can not find file '//trim(filename)
         stop '(read_sofa)'
      end if

      stat = nf90_open(trim(filename),nf90_nowrite,ncid) 
      if (stat /= nf90_noerr) then
         call ncerr( stat)
      end if

      print *,'reading SOFA fields from '//trim(filename)

      ! Dimension IDs and lengths (assuming the names are right)
      call getdimms(ncid,'longitude',londimid,londimlen)
      call getdimms(ncid,'latitude', latdimid,latdimlen)

      ! Variable IDs and dimensions
      call getvarstat(ncid,'longitude'    ,lonvarid,lonvndim,lonvdims)
      call getvarstat(ncid,'latitude'     ,latvarid,latvndim,latvdims)
      call getvarstat(ncid,'analyzed_sst' ,sstvarid,sstvndim,sstvdims)


      ! Read grid and time variables.
      allocate(longitude(londimlen))
      allocate(latitude (latdimlen))
      allocate(fld   (londimlen,latdimlen))
      allocate(tmpfld(londimlen,latdimlen))
      call ncerr( nf90_get_var(ncid, lonvarid, longitude) )
      call ncerr( nf90_get_var(ncid, latvarid, latitude ) )
      call ncerr( nf90_get_var(ncid, sstvarid, tmpfld) )

      ! Get sst add_offset and scale_factor
      call ncerr( nf90_get_att(ncid, sstvarid, 'add_offset'  , add_offset  ) )
      call ncerr( nf90_get_att(ncid, sstvarid, 'scale_factor', scale_factor) )
      call ncerr( nf90_get_att(ncid, sstvarid, '_FillValue'  , fillvalue   ) )
      fld=tmpfld*scale_factor+add_offset

      where (tmpfld==fillvalue) fld=undef


   end subroutine read_safo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The following are auxillary routines 



   ! Give netcdf error message
   subroutine ncerr(error)
      use netcdf
      implicit none
      integer, intent(in) :: error
      if (error/=nf90_noerr) then
         print *,'Read_Safo: '//nf90_strerror(error)
         stop
      end if
   end subroutine ncerr




   subroutine getdimms(ncid,dname,dmid,dmlen)
      use netcdf
      implicit none
      integer,          intent(in)  :: ncid
      character(len=*), intent(in)  :: dname
      integer,          intent(out) :: dmid,dmlen
      call ncerr( nf90_inq_dimid(ncid,dname,dmid) )
      call ncerr( nf90_inquire_dimension(ncid,dmid,len=dmlen) )
   end subroutine




   subroutine getvarstat(ncid,varname,varid,var_ndim,var_dimids)
      use netcdf
      implicit none
      integer,                          intent(in)  :: ncid
      character(len=*),                 intent(in)  :: varname
      integer,                          intent(out) :: varid,var_ndim
      integer, dimension(nf90_max_dims),intent(out) :: var_dimids
      call ncerr( nf90_inq_varid(ncid,varname,varid) )
      call ncerr( nf90_inquire_variable(ncid,varid,ndims=var_ndim,dimids=var_dimids) )
   end subroutine getvarstat


end module m_read_safo
