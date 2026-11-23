program p_ncofsstcmp
   use mod_year_info
   use mod_confmap
   use mod_hycomfile_io
   use netcdf
   use mod_grid
   use mod_xc
   use mod_za 
   use mod_ncof_sst
   implicit none

#if defined (IARGC)
   integer*4, external :: iargc
#endif
   character*80 :: Odir,Edir
   character*80 :: Ofile,Efile
   character(len=5)  :: Cdate

   ! Version info of this program

   ! Input file with days
   character(len=80) :: hyc_file
   integer            :: year,day,hour,find,ios,ifile, &
                         counter, refyear
   logical            :: ex
   logical            :: l_ncdump
   real,allocatable   :: Osst(:,:),fice(:,:),Esst(:,:)
   integer            :: k,k0
   logical            :: ex1
   
   ! Process input arguments
   if (iargc()/=3 ) then
      print *,'Usage:'
      print *,'  ncofsstcmp ncofsst-file hycom-file [ncdump]'
      stop
   else
      call getarg(1,Odir)
      call getarg(2,Edir)
      call getarg(3,hyc_file)

  !    if (iargc()==3) then
  !       call getarg(3,tmpchar)
  !       if (trim(tmpchar)=='-ncdump') then
  !          l_ncdump=.true.
  !       else
  !          print *,'Usage:'
  !          print *,'  ncofsstcmp ncofsst-file hycom-file [-ncdump]'
  !          stop
  !       end if
  !    end if
   end if
!   Odir="/work/xiejp/work_2015/OBS_TOPAZ/data0/sst/"
   Odir="/work/xiejp/work_2015/OBS_TOPAZ/downdata/"
   Edir="/work/xiejp/work_2015/OBS_TOPAZ/data0/sst/"

   ! Get bathymetry - lon/lat fields
   call xcspmd()
   call zaiost()
   call initconfmap(idm,jdm)
   call get_grid 
   allocate(Osst(idm,jdm),fice(idm,jdm),Esst(idm,jdm))

   k0=0;
!   do k=21915,23375
!   do k=17897,21914
   do k=14610,17897
     write(Cdate,'(i5)') k
     Ofile=trim(Odir)//"ncof_sst_"//trim(Cdate)//'.nc'
     Efile=trim(Edir)//"ncof_TP4_"//trim(Cdate)//'.nc'
     if (k0.eq.0) then
       ! Initialize the hycomfile module - and retrieve vertical levels
       call ncof_sst_info(Ofile,plon,plat,idm,jdm)
       k0=k0+1
     endif
     print *,trim(Efile)
     inquire(file=trim(Ofile),exist=ex1)
     if (ex1) then
       call ncof_sst_read(trim(Ofile),Osst,fice,Esst,idm,jdm,depths)
       call ncof_sst_dump(trim(Efile),idm,jdm,plon,plat,Osst,fice,Esst)
     endif

   end do

   print *,'end'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   ! Interpolation done - the rest is diagnostics
   deallocate(Osst,fice,Esst)


end program p_ncofsstcmp
