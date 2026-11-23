module m_read_coriolis
private :: errcond, getdimms,getvarstat
contains

   !Reads Coriolis data from a NetCDF file

   subroutine read_coriolis(filename,var,pres,longitude,latitude,juld,floatid,varname,undef)
      use netcdf
      implicit none

      ! Dummy variables
      character(len=*),                intent(in)  :: filename,varname
      real, pointer,   dimension(:)   :: longitude,latitude,juld
      real, pointer,   dimension(:,:) :: var,pres
      character(len=8),dimension(:)   , pointer    :: floatid
      real, intent(in) :: undef

      ! Local variables
      character(len=1), allocatable, dimension(:,:) :: &
         qcvar, qcpres,qcdeph,dstate,dmode
      character(len=1), allocatable, dimension(:) :: &
         qcpos
      integer, dimension(nf90_max_dims) ::                 &
         lonvdims,latvdims,vardims,qcvardims, qcpresvdims, &
         presvdims,juldvdims,qcposvdims,dephvdims,qcdephvdims, &
         pltfvdims,dstatvdims, dmodevdims
      integer ::    &
         lonvndim, latvndim, varndim, qcvarndim,       &
         qcpresvndim, presvndim, juldvndim,qcposvndim, &
         dephvndim, qcdephvndim,pltfvndim,dstatvarid,dstatvndim, &
         dmodevarid,dmodevndim
      real, allocatable :: deph(:,:)

      integer :: profdimid,profdimlen,zlevdimid,zlevdimlen,s8dimid,s8dimlen, &
         s4dimlen, s4dimid
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
      if (varname/='PSAL' .and. varname/='TEMP') then
         print *,'Variable name '//varname//'is unknown'
         print *,'valid names are PSAL or TEMP'
         stop '(read_coriolis)'
      end if


      ! Check for existence of netcdf file
      inquire(file=trim(filename),exist=ex)
      if (.not.ex) then
         print *,'Can not find file '//trim(filename)
         stop '(read_coriolis)'
      end if

      print *,'reading Coriolis fields from '//trim(filename)
      stat = nf90_open(trim(filename),nf90_nowrite,ncid) 
      if (stat /= nf90_noerr) then
         call ncerr( stat)
      end if

      !print *,'reading Coriolis fields from '//trim(filename)


      ! Dimension IDs and lengths (assuming the names are right)
      call getdimms(ncid,'N_PROF'  ,profdimid  ,profdimlen)
      call getdimms(ncid,'N_LEVELS',zlevdimid  ,zlevdimlen)
      call getdimms(ncid,'N_PARAM' ,nparamdimid,nparamdimlen)
      call getdimms(ncid,'STRING16',s16dimid   ,s16dimlen)
      call getdimms(ncid,'STRING8' ,s8dimid    ,s8dimlen)
      call getdimms(ncid,'STRING4' ,s4dimid    ,s4dimlen)
      print *,'# Profiles , # levels ',profdimlen, zlevdimlen
      !print *,profdimid, zlevdimid

      ! Better safe than sorry ...
      if (profdimlen<1 .or. zlevdimlen<1) then
         print *,'Dimensions are zero for this file ... '
         print *,'Exiting..'
         stop '(read_coriolis)'
      end if


      ! Get data state indicator, mode and station parameters
      allocate(dstate   (s4dimlen,profdimlen))
      allocate(dmode    (1       ,profdimlen))
      do j=1,profdimlen
         
         ! Data state indicator
         call ncerr( nf90_get_var(ncid, dstatvarid, dstate(:,j), &
             start=(/1,j/),count=(/s4dimlen,1/)))
         print *,dstate(:,j)

         ! Data mode (R)ealtime, (D)elayed or (A)djusted
         call ncerr( nf90_get_var(ncid, dmodevarid, dmode (:,j), &
             start=(/1,j/),count=(/1,1/)))
         print *,dmode(1,j)

         ! Station parameters
      end do



      ! Variable IDs and dimensions
      lpres=probe_var(ncid,'PRES')
      ldeph=probe_var(ncid,'DEPH')
      call getvarstat(ncid,'DATA_STATE_INDICATOR'  ,dstatvarid ,dstatvndim  ,dstatvdims)
      call getvarstat(ncid,'DATA_MODE'     ,dmodevarid ,dmodevndim  ,dmodevdims)
      call getvarstat(ncid,'LONGITUDE'    ,lonvarid ,lonvndim  ,lonvdims)
      call getvarstat(ncid,'LATITUDE'     ,latvarid ,latvndim  ,latvdims)
      call getvarstat(ncid,'POSITION_QC'  ,qcposvarid ,qcposvndim  ,qcposvdims)
      call getvarstat(ncid,'JULD'         ,juldvarid,juldvndim  ,juldvdims)
      if (lpres) then
         call getvarstat(ncid,'PRES'         ,presvarid,presvndim ,presvdims )
         call getvarstat(ncid,'PRES_QC'      ,qcpresvarid   ,qcpresvndim ,qcpresvdims )
      end if
      if (ldeph) then
         call getvarstat(ncid,'DEPH'         ,dephvarid,dephvndim ,dephvdims )
         call getvarstat(ncid,'DEPH_QC'      ,qcdephvarid   ,qcdephvndim ,qcdephvdims )
      end if
      call getvarstat(ncid,varname        ,varid    ,varndim   ,vardims )
      call getvarstat(ncid,varname//'_QC' ,qcvarid       ,qcvarndim   ,qcvardims )
      call getvarstat(ncid,'PLATFORM_NUMBER'  ,pltfvarid ,pltfvndim  ,pltfvdims)
      !print *,pltfvarid,pltfvndim


      ! Read variables
      allocate(floatid  (profdimlen))
      allocate(longitude(profdimlen))
      allocate(latitude (profdimlen))
      allocate(juld     (profdimlen))
      allocate(qcpos    (profdimlen))
      allocate(var      (zlevdimlen,profdimlen))
      allocate(pres     (zlevdimlen,profdimlen))
      allocate(deph     (zlevdimlen,profdimlen))
      allocate(qcvar    (zlevdimlen,profdimlen))
      allocate(qcpres   (zlevdimlen,profdimlen))
      allocate(qcdeph   (zlevdimlen,profdimlen))
      call ncerr( nf90_get_var(ncid, lonvarid, longitude) )
      call ncerr( nf90_get_var(ncid, latvarid, latitude ) )
      call ncerr( nf90_get_var(ncid, juldvarid, juld ) )
      if (lpres) call ncerr( nf90_get_var(ncid, presvarid, pres ) )
      if (ldeph) call ncerr( nf90_get_var(ncid, dephvarid, deph ) )
      call ncerr( nf90_get_var(ncid, varid, var) )





      ! New fast way
      do j=1,profdimlen
         if (lpres) then
            call ncerr( nf90_get_var(ncid, qcpresvarid, qcpres(:,j), &
                start=(/1,j/),count=(/zlevdimlen,1/)))
         endif
         if (ldeph) then
            call ncerr( nf90_get_var(ncid, qcdephvarid, qcdeph(:,j), &
                start=(/1,j/),count=(/zlevdimlen,1/)))
            !print *,qcdeph(:,j)
         endif
         call ncerr( nf90_get_var(ncid, qcvarid    , qcvar(:,j), &
                start=(/1,j/),count=(/zlevdimlen,1/)))
         call ncerr( nf90_get_var(ncid, pltfvarid    ,floatid(j)(:),   &
                start=(/1,j/),count=(/s8dimlen,1/)) )
         !print *,floatid(j)
      end do
      call ncerr( nf90_get_var(ncid, qcposvarid    , qcpos,   &
         start=(/1/),count=(/profdimlen/)) )



      ! Go through data and discard poor observations ("undef" them)
      do j=1,profdimlen
      do i=1,zlevdimlen
         if(qcvar (i,j) /= '1' .and. qcvar (i,j) /= '2' .and. qcvar(i,j) /='0') &
            var(i,j)=undef
         if (lpres) then
            if(qcpres(i,j) /= '1' .and. qcpres(i,j) /= '2' .and. qcpres(i,j) /='0') &
               pres(i,j)=undef
         end if
         if (ldeph) then
            if(qcdeph(i,j) /= '1' .and. qcdeph(i,j) /= '2' .and. qcdeph(i,j) /='0') &
               deph(i,j)=undef
         end if
      end do
      end do

      ! Select what to return of deph and pres .. 
      if (ldeph .and. .not. lpres) pres=deph
      if (ldeph.and.lpres) then
         do j=1,profdimlen
            if (all(abs(pres(:,j)-undef)<1e-4)) pres(:,j)=deph(:,j)
         end do
      end if
         

      do j=1,profdimlen
         if(qcpos (j) /= '1' .and. qcpos (j) /= '2' .and. qcpos(j) /='0') then
            longitude(j)=undef
            latitude (j)=undef
         end if
      end do


      ! Deallocate data
      deallocate(qcvar,qcpres,qcpos,deph,qcdeph)


   end subroutine read_coriolis


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

   logical function probe_var(ncid,varname)
      use netcdf
      implicit none
      character(len=*), intent(in) :: varname
      integer         , intent(in) :: ncid
      integer :: test, varid


      test=nf90_inq_varid(ncid,varname,varid) 
      probe_var = test==nf90_noerr
   end function probe_var



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


end module m_read_coriolis
