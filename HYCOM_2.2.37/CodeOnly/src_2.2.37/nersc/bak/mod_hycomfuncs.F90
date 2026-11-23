module mod_hycomfuncs
use mod_xc
use mod_year_info
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This retrieves the field corresponding to variable name, layer index and
! time level. Calls restart_getindex first
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine restart_getfield(filebase,varname,ilayer,itime,field,ierr)
   use mod_za
   implicit none
   character(len=*), intent(in)  :: filebase, varname
   integer,          intent(in)  :: ilayer, itime
   real,             intent(out) :: field(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
   integer,          intent(out) :: ierr

   integer :: indx,k
   real    :: bmin, bmax, amin, amax
   integer, parameter ::  nop=711

   ierr = 0
   call restart_getindex(filebase,varname,ilayer,itime,indx,bmin,bmax)
   if (indx<=0)  then
      if (mnproc==1) print '(a,2i5)','restart_getfield: Can not find field '// &
         trim(varname),ilayer, itime
      ierr=-1
      return
   end if
   if (mnproc==1) write(lp,'(a,i4,a,2i4,1x,a)'), &
      'restart_getfield: index is ',indx, ' for ',ilayer,itime,varname

   ! Safe to open
   call zaiopf(filebase(1:len_trim(filebase))//'.a','old', nop)

   ! Skip to indx "indx"
   do k=1,indx-1
      call zaiosk(nop)
   end do
   call zaiord(field, ip,.false., amin,amax, nop)

   if     (abs(amin-bmin).gt.abs(bmin)*1.e-4 .or. abs(amax-bmax).gt.abs(bmax)*1.e-4 ) then
      if     (mnproc.eq.1) then
         write(lp,'(/ a / a,2i3 / a,2i3 /  a,1p3e14.6 / a,1p3e14.6 /)') &
               'error - .a and .b files not consistent:',               &
               'iunit,indx = ',711,indx,                                &
               'vname,itime,ilayer = '//trim(varname),itime,ilayer,     &
               '.a,.b min = ',amin,bmin,amin-bmin,                      &
               '.a,.b max = ',amax,bmax,amax-bmax 
      endif !1st tile
      ierr=-1
      ! Deemed a critical error - something is way wrong
      call xcstop('(restart_getfield)')
      stop '(restart_getfield)'
   endif
   call zaiocl(nop)
   end subroutine




   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This appends the field with varname, ilayer, itime to the file with base
! name "filebase".  Retrieves last index of .a-file first
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine restart_appendfield(filebase,varname,ilayer,itime,field,ierr)
   use mod_za
   implicit none
   character(len=*), intent(in)    :: filebase
   character(len=8), intent(in)    :: varname
   integer,          intent(in)    :: ilayer, itime
   real,             intent(inout) :: field(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
   integer,          intent(out)   :: ierr

   real    :: xmin, xmax
   integer :: indx,k
   integer, parameter :: nop=711

   ierr=0 ! Do something more meaningful later on
   indx=restart_lastindex(trim(filebase))
   call xcsync(flush_lp)
   if (mnproc==1) print *,indx
   call zaiopf(filebase(1:len_trim(filebase))//'.a','unknown', nop)
   if (mnproc==1) then
      open (unit=nop,file=filebase(1:len_trim(filebase))//'.b', &
            status='unknown',action='readwrite',form='formatted')
   end if
   if (mnproc==1) write(lp,'(a,i4,a,2i4,1x,a)'), &
      'restart_appendfield: index is ',indx, ' for ',ilayer,itime,varname

   !skip first two lines + indx lines
   do k=1,2+indx
      if (mnproc==1) read(nop,*)
   end do
   do k=1,indx
      call zaiosk(nop)
   end do

   call zaiowr(field(1-nbdy,1-nbdy),ip,.false., xmin,xmax, nop, .false.)
   if (mnproc==1) write(nop,4100) varname,ilayer,itime,xmin,xmax
   if (mnproc==1) close(nop)
   call zaiocl(nop)
 4100 format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)
   end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This retrieves the .ab-index corresponding to variable name, layer index and
! time level. For now every thread opens the .b file which is not optimal.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine restart_getindex(filebase,varname,ilayer,itime,fldindx,bmin,bmax)
   use mod_za
   implicit none
   character(len=*), intent(in)  :: filebase, varname
   integer,          intent(in)  :: ilayer, itime
   integer,          intent(out) :: fldindx
   real,             intent(out) :: bmin, bmax

   integer :: i1,i2,k,ios, indx
   character(len=8) :: char8
   real :: xmin, xmax
   logical :: match
   integer, parameter :: nop=711

   open (unit=nop,file=filebase(1:len_trim(filebase))//'.b', &
                  status='unknown',action='read',form='formatted')

   ! Skip two first lines (header)
   fldindx=-1
   read(nop,*,iostat=ios) 
   read(nop,*,iostat=ios) 
   if (ios/=0) then
      close(nop)
      return
   end if

   ! Read until match
   ios=0
   indx=0
   match=.false.
   do while   (ios==0)
      read(nop,4100,iostat=ios) char8,i1,i2, xmin, xmax
      indx=indx+1
      if (trim(char8)==trim(varname) .and. i1==ilayer .and. i2==itime .and.  ios==0) then
         fldindx=indx
         bmin=xmin
         bmax=xmax
         match=.true.
      end if
   end do
   close(nop)
!  if (mnproc==1) write(lp,'(a,i4,a,2i4,1x,a)'), &
!     'restart_getindex is ',fldindx, ' for ',ilayer,itime,varname
4100 format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)
   end subroutine restart_getindex
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This retrieves the number of the last index in .ab files. Input is file
! base (without .a or .b ending)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function restart_lastindex(filebase)
   use mod_za
   implicit none
   character(len=*), intent(in) :: filebase
   integer :: i1,i2,k,ios, indx
   character(len=8) :: char8
   real :: xmin, xmax, rli
   integer, parameter :: nop=711
   integer :: iunitb, iunita

   restart_lastindex=-1
   rli=-1.
   if (mnproc==1)  then
      !inquire(file=filebase(1:len_trim(filebase))//'.b',number=iunitb)
      !inquire(file=filebase(1:len_trim(filebase))//'.a',number=iunita)
      open (unit=nop,file=filebase(1:len_trim(filebase))//'.b', &
                     status='unknown',action='read',form='formatted')

      ! Skip two first lines
      read(nop,*,iostat=ios) 
      read(nop,*,iostat=ios) 
      if (ios/=0) then
         close(nop)
         restart_lastindex=-1
         return
      end if
      ! Read until EOF
      ios=0
      indx=0
      do while   (ios==0)
         read(nop,4100,iostat=ios) char8,i1,i2, xmin, xmax
         !if (mnproc==1) print *,char8,i1,i2
         indx=indx+1
      end do
      close(nop)
      restart_lastindex=indx-1
      rli=restart_lastindex
   end if
   call xcmaxr(rli) ; restart_lastindex=rli
4100 format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)
   end function restart_lastindex


   ! Base name of restart file
   function restartfile(rungen,rt,imem)
      use mod_year_info
      implicit none
      character(len=80) restartfile
      character(len=3) , intent(in) :: rungen
      integer,           intent(in) :: imem
      type(year_info) rt
      character(len=3) :: cmem
      integer indx
      restartfile=''
      restartfile(1:3)=rungen(1:3)
      restartfile(4:10)='restart'
      restartfile(11:14)=rt%cyy(1:4)
      restartfile(15:15)='_'
      restartfile(16:18)=rt%cdd(1:3)
      restartfile(19:19)='_'
      restartfile(20:21)=rt%chh(1:2)
      ! If imem > 1 - append member index ...
      if (imem>1) then
         write(cmem,'(i3.3)') imem
         restartfile=restartfile(1:len_trim(restartfile))//'_mem'//cmem
      end if
   end function restartfile


!   function filemem(rungen,rt)
!      character(len=24) filemem
!      character(len=3) rungen
!      type(year_info) rt
!      filemem(1:3)=rungen(1:3)
!      filemem(4:10)='restart'
!      filemem(11:14)=rt%cyy(1:4)
!      filemem(15:15)='_'
!      filemem(16:18)=rt%cdd(1:3)
!      filemem(19:19)='_'
!      filemem(20:21)=rt%chh(1:2)
!      filemem(22:24)='.uf'
!   end function filemem

   function fileice(rungen,rt)
      character(len=27) fileice
      character(len=3) rungen
      type(year_info) rt
      fileice(1:3)=rungen(1:3)
      fileice(4:10)='restart'
      fileice(11:14)=rt%cyy(1:4)
      fileice(15:15)='_'
      fileice(16:18)=rt%cdd(1:3)
      fileice(19:19)='_'
      fileice(20:21)=rt%chh(1:2)
      fileice(22:27)='ICE.uf'
   end function fileice

!   function fileeco(rungen,rt)
!      character(len=27) fileeco
!      character(len=3) rungen
!      type(year_info) rt
!      fileeco(1:3)=rungen(1:3)
!      fileeco(4:10)='restart'
!      fileeco(11:14)=rt%cyy(1:4)
!      fileeco(15:15)='_'
!      fileeco(16:18)=rt%cdd(1:3)
!      fileeco(19:19)='_'
!      fileeco(20:21)=rt%chh(1:2)
!      fileeco(22:27)='ECO.uf'
!   end function fileeco
!
!   function fileave(rungen,rt)
!      character(len=19) fileave
!      type(year_info),  intent(in) :: rt
!      character(len=3), intent(in) :: rungen
!      fileave(1:3)=rungen(1:3)
!      fileave(4:9)='MULTI_'
!      fileave(10:13)=rt%cyy(1:4)
!      fileave(14:14)='_'
!      fileave(15:16)=rt%cmm(1:2)
!      fileave(17:19)='.uf'
!   end function fileave
!
!   function fileaveweek(rungen,rt)
!      implicit none
!      character(len=18) fileaveweek
!      type(year_info),  intent(in) :: rt
!      character(len=3), intent(in) :: rungen
!      character(len=2) tag2
!      write(tag2,'(i2.2)')int( (float(rt%idd)-0.1)/7. ) + 1
!      fileaveweek(1:3)=rungen(1:3)
!      fileaveweek(4:8)='WEEK_'
!      fileaveweek(9:12)=rt%cyy(1:4)
!      fileaveweek(13:13)='_'
!      fileaveweek(14:15)=tag2(1:2)
!      fileaveweek(16:18)='.uf'
!   end function fileaveweek

   function fileaveres(rungen,rt)
      character(len=27) fileaveres
      character(len=3) rungen
      type(year_info) rt
      fileaveres(1:3)=rungen(1:3)
      fileaveres(4:10)='restart'
      fileaveres(11:14)=rt%cyy(1:4)
      fileaveres(15:15)='_'
      fileaveres(16:18)=rt%cdd(1:3)
      fileaveres(19:19)='_'
      fileaveres(20:21)=rt%chh(1:2)
      fileaveres(22:27)='AVE.uf'
   end function fileaveres

!   function check_tsteps(baclin,batrop,lstep,dlt)
!      implicit none
!      logical check_tsteps
!      real baclin,batrop
!      integer, intent(out) :: lstep
!      real, intent(out) :: dlt
!!    'lstep' = number of barotropic time steps per baroclinic time step.
!!     lstep   m u s t   be even.
!      check_tsteps=.true.
!      if (mod(baclin,batrop) /= 0) then
!         print *,'baclin is not an integer times batrop'
!         print *,'baclin=',baclin
!         print *,'batrop=',batrop
!         check_tsteps=.false.
!      endif
!      lstep=nint(baclin/batrop)
!      if (mod(lstep,2) /= 0) then
!         print *,'lstep is not even',lstep
!         check_tsteps=.false.
!      endif
!      write (*,'(i4,'' barotropic steps per baroclinic time step'')')lstep
!      dlt=baclin/lstep
!
!   end function check_tsteps

!   subroutine set_random_seed(seed)
!      implicit none
!      integer, intent(in) :: seed        ! seed read from limits.dat
!      integer seedsize                   ! dimension of seedarray
!      integer, allocatable :: seedarray(:)
!      real ttt
!      call random_seed
!      call random_seed(size=seedsize)
!      allocate (seedarray(seedsize))
!      if (seedsize >= 2) seedarray(2)=seed
!      call random_seed(put = seedarray(1:seedsize))
!      call random_number(ttt)
!      print *,'ttt',ttt
!      deallocate(seedarray)
!   end subroutine set_random_seed

end module mod_hycomfuncs
