module m_dump_coriolis
private :: errcond, getdimms,getvarstat
contains

   !Reads Coriolis data from a NetCDF file

   subroutine dump_coriolis(filename,varname,var,ndepths,nprof,undef)
      use netcdf
      implicit none

      ! Dummy variables
      integer, intent(in) :: nprof, ndepths
      character(len=*),             intent(in)  :: filename,varname
      real, dimension(ndepths,nprof) :: var
      real, intent(in) :: undef

      ! Local variables
      character(len=1), allocatable, dimension(:,:) :: &
         qcvar, qcpres,qcdeph
      character(len=1), allocatable, dimension(:) :: &
         qcpos
      integer, dimension(nf90_max_dims) ::                 &
         lonvdims,latvdims,vardims,qcvardims, qcpresvdims, &
         presvdims,juldvdims,qcposvdims,dephvdims,qcdephvdims, &
         pltfvdims
      integer ::    var_id,&
         lonvndim, latvndim, varndim, qcvarndim,       &
         qcpresvndim, presvndim, juldvndim,qcposvndim, &
         dephvndim, qcdephvndim,pltfvndim
      real, allocatable :: deph(:,:)

      integer :: profdimid,profdimlen,zlevdimid,zlevdimlen,s8dimid,s8dimlen
      integer :: &
         lonvarid,latvarid, presvarid, qcvarid, pltfvarid, &
         qcpresvarid, varid,juldvarid,qcposvarid,dephvarid,qcdephvarid
      integer :: sstvarid,sstvndim
      integer :: ncid
      integer :: i,stat,j
      logical :: ex
      real :: updiff,lwdiff,diff
      real :: add_offset,scale_factor, fillvalue
      logical :: lpres, ldeph

      ! Check if variable type is correct
      !if (varname/='PSAL' .and. varname/='TEMP') then
      !   print *,'Variable name '//varname//'is unknown'
      !   print *,'valid names are PSAL or TEMP'
      !   stop '(read_coriolis)'
      !end if


      ! Check for existence of netcdf file
      inquire(file=trim(filename),exist=ex)
      if (.not.ex) then
         print *,'Can not find file '//trim(filename)
      end if

      print *,'reading Coriolis fields from '//trim(filename)
      stat = nf90_open(trim(filename),nf90_write,ncid) 
      if (stat /= nf90_noerr) then
         call ncerr( stat)
      end if


      ! Dimension IDs and lengths (assuming the names are right)
      call getdimms(ncid,'N_PROF',profdimid,profdimlen)
      call getdimms(ncid,'N_LEVELS',zlevdimid,zlevdimlen)

      ! Naive security cehck
      if (nprof/=profdimlen .or. zlevdimlen/=ndepths) then
         print *,'Dimension mismatch in dump_coriolis'
         stop
      end if

      ! Define variable     
      call ncerr(nf90_redef(ncid))
      call ncerr(NF90_DEF_VAR(ncid,varname,NF90_FLOAT,(/zlevdimid,profdimid/),var_id))
      call ncerr(NF90_PUT_ATT(ncid,var_id,'_FillValue'   ,undef))
      call ncerr(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
      call ncerr(nf90_enddef(ncid))
      call ncerr(NF90_PUT_VAR(ncid,var_id,var))
      call ncerr(NF90_CLOSE(ncid))



   end subroutine dump_coriolis


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The following are auxillary routines 



   ! Give netcdf error message
   subroutine ncerr(error)
      use netcdf
      implicit none
      integer, intent(in) :: error
      if (error/=nf90_noerr) then
         print *,'read_coriolis: '//nf90_strerror(error)
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


end module m_dump_coriolis
