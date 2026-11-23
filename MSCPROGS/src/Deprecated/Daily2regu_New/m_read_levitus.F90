module m_read_levitus
private :: errcond, getdimms,getvarstat
contains

   !Reads Levitus data from a NetCDF file 

   subroutine read_levitus(climid,fldname,fld,depths,reflon,dlon,reflat,dlat)
      use netcdf
      implicit none

      ! Dummy variables
      character(len=*),                intent(in)  :: fldname,climid
      real, pointer,dimension(:)       :: depths
      real, pointer,dimension(:,:,:,:) :: fld
      real, intent(out) :: reflon,reflat,dlon,dlat

      ! Local variables
      character(len=80) :: fldfile
      integer, dimension(nf90_max_dims) :: xvdims,yvdims,zvdims,tvdims,flddims
      integer :: xdimid,ydimid,zdimid,tdimid,xdimlen,ydimlen,zdimlen,tdimlen
      integer :: xvndim,yvndim,zvndim,tvndim,ncid,xvarid,yvarid,zvarid
      integer :: tvarid,fldndim,fldvarid,i,stat
      logical :: ex
      real, allocatable, dimension(:) :: xvar,yvar,zvar,tvar
      real :: updiff,lwdiff,diff

      if (trim(climid)=="Levitus") then
         if (fldname=="sal") then
            fldfile = './Climatology/Levitus/season_salinity.cdf'
         else if (fldname=="temp") then
            fldfile = './Climatology/Levitus/season_temperature.cdf'
         else 
            print *,'Unknown Levitus field '//fldname
            stop '(read_Levitus)'
         end if
      else
         print *,'Only Levitus Climatology supported'
         stop '(read_Levitus)'
      end if


      ! Check for existence of netcdf file
      inquire(file=trim(fldfile),exist=ex)
      if (ex) then
         stat = nf90_open(trim(fldfile),nf90_nowrite,ncid) 
         if (stat == nf90_noerr) then

            print *,'reading Levitus fields from '//trim(fldfile)

            ! Dimension IDs and lengths (assuming the names are right)
            call getdimms(ncid,'X',xdimid,xdimlen)
            call getdimms(ncid,'Y',ydimid,ydimlen)
            call getdimms(ncid,'Z',zdimid,zdimlen)
            call getdimms(ncid,'T',tdimid,tdimlen)

            !print *,xdimlen,ydimlen,zdimlen,tdimlen

            ! Variable IDs and dimensions
            call getvarstat(ncid,'X',xvarid,xvndim,xvdims)
            call getvarstat(ncid,'Y',yvarid,yvndim,yvdims)
            call getvarstat(ncid,'Z',zvarid,zvndim,zvdims)
            call getvarstat(ncid,'T',tvarid,tvndim,tvdims)


            ! Error check -- array assumptions
            if (       xvndim > 1 .or. xvdims(1).ne.xdimid  .or. &
                       yvndim > 1 .or. yvdims(1).ne.ydimid  .or. &
                       zvndim > 1 .or. zvdims(1).ne.zdimid  .or. &
                       tvndim > 1 .or. tvdims(1).ne.tdimid ) then
                print *,'Assumed arrays does not match with data'
                stop
            end if

            ! Read grid and time variables.
            allocate(xvar(xdimlen))
            allocate(yvar(ydimlen))
            allocate(zvar(zdimlen))
            allocate(tvar(tdimlen))
            call ncerr( nf90_get_var(ncid, xvarid, xvar) )
            call ncerr( nf90_get_var(ncid, yvarid, yvar) )
            call ncerr( nf90_get_var(ncid, zvarid, zvar) )
            call ncerr( nf90_get_var(ncid, tvarid, tvar) )

            ! Error check -- dlon,reflon
            diff = xvar(2)-xvar(1) 
            updiff=diff ; lwdiff=diff
            do i=2,xdimlen
               diff   = xvar(i)-xvar(i-1)
               updiff = (max(updiff,diff))
               lwdiff = (min(lwdiff,diff))
            end do
            if (abs(updiff-lwdiff)>1e-8) then 
               stop 'Not constant longitude spacing'
            else
               dlon = diff
               reflon=xvar(1)
            end if

            ! Error check -- dlat,reflat
            diff = yvar(2)-yvar(1) 
            updiff=diff ; lwdiff=diff
            do i=2,ydimlen
               diff   = yvar(i)-yvar(i-1)
               updiff = (max(updiff,diff))
               lwdiff = (min(lwdiff,diff))
            end do
            if (abs(updiff-lwdiff)>1e-8) then 
               stop 'Not constant latitude spacing'
            else
               dlat = diff
               reflat=yvar(1)
            end if


            ! Read field stats
            call getvarstat(ncid,fldname,fldvarid,fldndim,flddims)

            ! Error check -- field dimensions
            if (fldndim==4 .and. flddims(1)==xdimid .and. &
                                 flddims(2)==ydimid .and. &
                                 flddims(3)==zdimid .and. &
                                 flddims(4)==tdimid ) then
               allocate(fld(xdimlen,ydimlen,zdimlen,tdimlen))
            else
               stop 'Error in field dimensions'
            end if

            ! Get field
            call ncerr( nf90_get_var(ncid,fldvarid,fld) )

            ! Depths
             allocate(depths(zdimlen))
             depths=zvar

         else
            
            call ncerr( stat)

         end if

      else

         print *,'No such file '//fldfile
         stop

      end if

   end subroutine read_levitus


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The following are auxillary routines 



   ! Give netcdf error message
   subroutine ncerr(error)
      use netcdf
      implicit none
      integer, intent(in) :: error
      if (error/=nf90_noerr) then
         print *,'Read_Levitus: '//nf90_strerror(error)
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


end module m_read_levitus
