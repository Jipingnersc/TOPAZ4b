module mod_read_netcdf
use netcdf


contains

   subroutine openNCfile(Fname,ncid)
     implicit none
     character(*), intent(in) :: Fname
     integer,     intent(out) :: ncid
     call handle_err(nf90_open(trim(Fname),nf90_nowrite,ncid),.true.)
   end subroutine

   subroutine readNCDim(ncid,Dimnam,Nlen)
     implicit none
     integer,             intent(in) :: ncid
     character(*),        intent(in) :: Dimnam
     integer,            intent(out) :: Nlen
     !-------------------------
     integer :: varid
     call handle_err(nf90_inq_dimid(ncid,trim(Dimnam),varid))
     call handle_err(nf90_inquire_dimension(ncid, varid, len = Nlen))
   end subroutine

   subroutine readNC1vec(ncid,Varnam,Vtype,Var,Nx)
     implicit none
     integer,             intent(in) :: ncid,Vtype
     character(*),        intent(in) :: Varnam
     real,dimension(:),  intent(out) :: Var
     integer,            intent(out) :: Nx
     !-------------------------
     integer :: varid
     integer,dimension(nf90_max_var_dims) :: vardims
     real,dimension(:),allocatable :: var0
     integer,dimension(:),allocatable :: var1
     call handle_err(nf90_inq_varid(ncid,trim(Varnam),varid),.true.)
     call handle_err(nf90_inquire_variable(ncid,varid,dimids=vardims),.true.)
     call handle_err(nf90_inquire_dimension(ncid,vardims(1),len=Nx))
     allocate(var0(Nx),var1(Nx))
     print *,Nx
     if(Nx>0.and.Vtype==2) then
       call handle_err(nf90_get_var(ncid,varid,var0),.true.)
       Var(1:Nx)=var0;
     elseif(Nx>0.and.Vtype==1) then
       call handle_err(nf90_get_var(ncid,varid,var1),.true.)
       Var(1:Nx)=real(var1);
     else
       print *, 'There are some error about '//trim(Varnam)//' in the nc file'
       Nx=0
     endif
     deallocate(var0,var1)
   end subroutine

   subroutine readNC3d(ncid,Varnam,Nx,Ny,Nz,Var)
     implicit none
     integer,     intent(in) :: ncid
     integer,     intent(in) :: Nx,Ny,Nz
     character(*),intent(in) :: Varnam
     real,       intent(out) :: Var(Nx,Ny,Nz)
     !-------------------------
     integer :: varid
     integer,dimension(nf90_max_var_dims) :: vardims
     integer :: N1,N2,N3
     call handle_err(nf90_inq_varid(ncid,trim(Varnam),varid),.true.)
     call handle_err(nf90_inquire_variable(ncid,varid,dimids=vardims),.true.)
!     print *,vardims(1:3)
     call handle_err(nf90_inquire_dimension(ncid,vardims(1),len=N1))
     call handle_err(nf90_inquire_dimension(ncid,vardims(2),len=N2))
     call handle_err(nf90_inquire_dimension(ncid,vardims(3),len=N3))
     print *,N1,N2,N3
     if(N1==Nx.and.N2==Ny.and.N3==Nz) then
       call handle_err(nf90_get_var(ncid,varid,Var),.true.)
     else
       print *, 'There are some error about '//trim(Varnam)//' in the nc file'
       return    
     endif
   end subroutine






   subroutine closeNCfile(ncid)
   implicit none
   integer,intent(in) :: ncid
   if (ncid==-1) then
      print *,'Cant close a closed file..'
      stop '(CloseNCFile)'
   end if
   if (NF90_CLOSE(ncid)/=NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(closeNCFile)'
   end if
   end subroutine

   subroutine handle_err(errcode,critical)
      use netcdf
      implicit none
      integer, intent(in) :: errcode
      logical, optional, intent(in) :: critical
      logical :: crit

      crit=.true.
      if (present(critical)) crit=critical
      if (errcode/=NF90_NOERR) then
         write(6,'(a)') NF90_STRERROR(errcode)
         if (crit) stop '(handle_err)'
      end if
   end subroutine
end module
