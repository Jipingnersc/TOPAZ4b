module m_def_netcdf_var

contains

subroutine def_netcdf_var(ncid,var_id,var_dims,ndim,varinfo1,dontpack,undefin)
use mod_netcdf_pars
use m_handle_err
use netcdf
implicit none
   integer, intent(in)  :: ncid, ndim
   integer, intent(out) ::var_id
   integer, intent(in), dimension(ndim) :: var_dims
   type(netcdf_pars), intent(in) :: varinfo1
   logical, optional, intent(in) :: dontpack
   real   , optional, intent(in) :: undefin

   logical :: dontpack2
   real    :: undef


   dontpack2=.false.
   if (present(dontpack)) dontpack2=dontpack
   call handle_err(nf90_redef(ncid))

   if (present(undefin)) then
      undef=undefin
   else
      undef=-1e14
   end if

   if (dontpack2) then
      call handle_err(NF90_DEF_VAR(ncid,trim(varinfo1%vname),NF90_FLOAT,var_dims,var_id))
      !call handle_err(NF90_PUT_ATT(ncid,var_id,'_FillValue'   ,varinfo1%fillvalue))
      !call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',varinfo1%fillvalue))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'_FillValue'   ,undef))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
   else
      call handle_err(NF90_DEF_VAR(ncid,trim(varinfo1%vname),NF90_SHORT,var_dims,var_id))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',0_2))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'_FillValue',0_2))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'add_offset',0.))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'scale_factor',1.))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'add_offset'   ,varinfo1%add_offset))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'scale_factor' ,varinfo1%scale_factor))
      !call handle_err(NF90_PUT_ATT(ncid,var_id,'_FillValue'   ,varinfo1%fillvalue))
      call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',varinfo1%fillvalue))
   end if
   call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name',     trim(var_longname(varinfo1))))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'units'    ,     trim(varinfo1%units))) 
   call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name', trim(var_standard_name(varinfo1))))

   call handle_err(nf90_enddef(ncid))
end subroutine
end module
