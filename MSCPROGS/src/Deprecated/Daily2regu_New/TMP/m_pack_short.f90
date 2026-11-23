

module m_pack_short

private :: get_predefined_range

contains



   subroutine pack_short(field3D,idm,jdm,kdm,undef,ncid,varid,varinfo1)
      use mod_netcdf_pars
      use netcdf
      use m_handle_err
      implicit none
      integer, intent(in)  :: idm,jdm,kdm
      real   , intent(in)  :: field3D(idm,jdm,kdm)
      real   , intent(in)  :: undef
      integer, intent(in)  :: ncid,varid
      !character(len=*), intent(in) :: varname
      type(netcdf_pars), intent(inout) :: varinfo1


      integer :: nbits,i,j,k
      integer :: field3Dint(idm,jdm,kdm)
      integer*2 :: fillvalue
      integer :: counthi,countlo

      ! Check for over/undershoot
      counthi=count(varinfo1%max_val<field3D)
      countlo=count(varinfo1%min_val>field3D)
      if (counthi>0 .and. countlo>0. ) then
         print *,'Over/Undershoot for variable '//varinfo1%vname
         print *,maxval(field3D)
         print *,minval(field3D)
         print *,countlo,counthi
      end if
      

      where (field3D/=undef)
         field3Dint = floor((max(min(field3D,varinfo1%max_val),varinfo1%min_val) - varinfo1%add_offset) / varinfo1%scale_factor)
      elsewhere
         field3Dint=varinfo1%fillvalue
      endwhere

      call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue',varinfo1%fillvalue))

      call handle_err(NF90_PUT_VAR(ncid,varid,field3dint))! ,start=(/1,1,1/)))
   end subroutine












   subroutine pack_short2D(field2D,idm,jdm,undef,ncid,varid,varinfo1)
      use mod_netcdf_pars
      use netcdf
      use m_handle_err
      implicit none
      integer, intent(in)  :: idm,jdm
      real   , intent(in)  :: field2D(idm,jdm)
      real   , intent(in)  :: undef
      integer, intent(in)  :: ncid,varid
      type(netcdf_pars), intent(inout) :: varinfo1


      integer :: nbits,i,j,k
      integer :: field2Dint(idm,jdm)
      integer*2 :: fillvalue
      integer :: counthi,countlo

      ! Check for over/undershoot
      counthi=count(varinfo1%max_val<field2D)
      countlo=count(varinfo1%min_val>field2D)
      if (counthi>0 .and. countlo>0. ) then
         print *,'Over/Undershoot for variable '//varinfo1%vname
         print *,countlo,counthi
      end if

      where (field2D/=undef)
         field2Dint = floor((max(min(field2D,varinfo1%max_val),varinfo1%min_val) - varinfo1%add_offset) / varinfo1%scale_factor)
      elsewhere
         field2Dint=varinfo1%fillvalue
      endwhere

      call handle_err(NF90_PUT_ATT(ncid,varid,'_FillValue'   ,varinfo1%fillvalue))
      call handle_err(NF90_PUT_VAR(ncid,varid,field2dint))! ,start=(/1,1,1/)))
   end subroutine









   subroutine get_predefined_range(varname,minfld,maxfld)
      implicit none
      character(len=*), intent(in) :: varname
      real, intent(out) :: minfld,maxfld

      if (varname=="temperature") then
         minfld=-2.0
         maxfld=45.0
      else if (varname=="salinity") then
         minfld=10.0
         maxfld=43.0
      else if (varname=="u"  .or. varname=="v" .or. &
               varname=="ub" .or. varname=="vb") then
         minfld=-5;
         maxfld=5.;
      else if (varname=="ssh") then
         minfld=-3;
         maxfld=3.;
      else if (varname=="bsfd") then
         minfld=-1e10
         maxfld=1e10
      else if (varname=="mosf" .or. varname == "mosfsig0" .or. &
               varname=="mosftheta0") then
         minfld=-200
         maxfld=200
      else if (varname=="mld") then
         minfld=0.
         maxfld=6000.
      else if (varname=="mld" .or. varname=="mlp") then
         minfld=0.
         maxfld=6000.
      else if (varname=="tauy" .or. varname=="taux") then
         minfld=-2.0
         maxfld= 2.0
      else if (varname=="qtot" .or. varname=="qsw") then
         minfld=-3000.
         maxfld= 3000.
      else if (varname=="emp") then
         minfld=-1e-3
         maxfld= 1e-3
      else if (varname=="interface" ) then
         minfld=0.
         maxfld=14000.
      else if (varname=="fice" ) then
         minfld=0.
         maxfld=1.001
      else if (varname=="hice" ) then
         minfld=0.
         maxfld=12.
      else if (varname=="uice" ) then
         minfld=-2.
         maxfld=2.
      else if (varname=="vice" ) then
         minfld=-2.
         maxfld=2.
      else if (varname=="wice" ) then
         minfld=0.
         maxfld=4.
      else 
         print *,'No packing method defined for variable '//varname
         stop '(get_predefined_range)'
      end if
   end subroutine

end module






