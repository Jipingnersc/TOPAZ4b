module m_getfileinfo
contains

      subroutine getfileinfo(fbase,hycfile,rti,rtd,nrmem)
      use mod_year_info
      use mod_read_weekab
      use mod_read_dailyab
      use mod_read_rstab
      implicit none
      type(year_info) , intent(out) :: rti,rtd
      character(len=*), intent(in ) :: fbase
      integer         , intent(in)  :: hycfile
      integer         , intent(out) :: nrmem
      integer :: nrmem, idm,jdm,kdm

      ! Read different fields - returns undef if not defined
      if (hycfile==2) then
         call daily_average_read_header(trim(fbase),rti,rtd,nrmem,idm,jdm,kdm)
      else if (hycfile==1) then
         call rst_read_header(trim(fbase),rti,nrmem,idm,jdm,kdm)
         rtd=rti
         nrmem=1
      else if (hycfile==3) then
         call week_read_header(trim(fbase),rti,idm,jdm,kdm)
         rtd=rti
         nrmem=1
      end if


      end subroutine getfileinfo

end module m_getfileinfo
