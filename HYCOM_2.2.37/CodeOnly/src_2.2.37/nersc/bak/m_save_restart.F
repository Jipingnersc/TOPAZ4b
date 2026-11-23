      module m_save_restart
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function save_restart_mem(rungen,rt,imem,nstep0,dtime0)
      use mod_xc
      use mod_za
      use mod_hycomfuncs
      implicit none

      logical :: save_restart_mem
      type(year_info),  intent(in) :: rt
      character(len=3), intent(in) :: rungen
      integer,          intent(in) :: imem
      integer,          intent(in) :: nstep0
      real*8,           intent(in) :: dtime0

      character(len=80) filename 
      include 'common_blocks.h'


      !! use hycom restart routine. Modify restart filename
      filename=trim(restartfile(rungen,rt,imem))
      flnmrso=trim(filename) 
      if (mnproc==1)  write(lp,*) 'Dumping to file '//
     &   trim(flnmrso)//'.[ab]'
      call restart_out(nstep0,dtime0)

      ! Hmm .. not closed automatically ...
      call zaiocl(12)
      if (mnproc==1)  then
         close(12)
      end if

      save_restart_mem=.true.
      return
      end function save_restart_mem






      function save_restart_ice(rungen,rt,imem)
      use mod_xc
      use mod_hycomfuncs
      use mod_year_info
      use mod_common_ice
      implicit none

      type(year_info),       intent(in)  :: rt
      character(len=3),      intent(in) :: rungen
      integer,               intent(in) :: imem


      ! Local arrays for file-saving
      real, dimension(itdm,jtdm)  :: 
     &   tficem,thicem,thsnwm,tticem,ttsrfm,ticeU,ticeV
      real*8, dimension(itdm,jtdm)  ::  ! Hard-coded real*8 for file IO
     &   io_ficem,io_hicem,io_hsnwm,io_ticem,io_tsrfm,io_iceU,io_iceV

      logical :: save_restart_ice

      character(len=27) icename
      integer i,j
#if ! defined (SHMEM) 
      integer iostat
#endif
      integer ioproc, ierr
      
#if defined(MPI)
      include 'mpif.h'
#elif defined (SHMEM)
      include 'mpp/shmem.fh'
      INTEGER, save :: PSYNC(SHMEM_REDUCE_SYNC_SIZE)
      INTEGER, save :: PWRK(MAX(1/2+1,SHMEM_REDUCE_MIN_WRKDATA_SIZE))
      integer, save :: isaveiceall, isaveice, alliostat, iostat
#endif


#if defined (SINGLE_RESTART) /* We always dump when this flag is set */
      save_restart_ice=save_restart_ice_ab(rungen,rt,imem)
#if defined (SINGLE_RESTART_ONLY) /* Short-circuit when this flag is set */
      return
#endif /* SINGLE_RESTART_ONLY */
#endif /* SINGLE_RESTART */

      ioproc=mnproc 
#if defined(MPI)
      ! Broadcast "io" mnproc (ioproc) from node with rank 0
      call mpi_bcast(ioproc,1,MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
#endif

      ! NB get stop - signals valid for MPI runs as well !

      ! Put ice data into non-tiled arrays
      call xcaget(tficem,ficem,0) ; io_ficem=tficem
      call xcaget(thicem,hicem,0) ; io_hicem=thicem
      call xcaget(thsnwm,hsnwm,0) ; io_hsnwm=thsnwm
      call xcaget(tticem,ticem,0) ; io_ticem=tticem
      call xcaget(ttsrfm,tsrfm,0) ; io_tsrfm=ttsrfm
      call xcaget(ticeU ,iceU ,0) ; io_iceu =ticeu 
      call xcaget(ticeV ,iceV ,0) ; io_icev =ticev

      ! Only one node is alloved to access this file (only important for MPI)
      save_restart_ice=.true.
      if (mnproc==ioproc) then
         icename=fileice(rungen,rt)
         write(lp,'(a,i4)',advance='yes')'Saving : '//
     &      trim(icename)//'                     record=',imem
         !inquire(iolength=j)tficem,thicem,thsnwm,tticem,ttsrfm  !,ticeU,ticeV
         inquire(iolength=j)
     &      io_ficem,io_hicem,io_hsnwm,io_ticem,io_tsrfm  !,ticeU,ticeV
         !print *,'iolength=',j ! Not important
         open(10,file=icename,status='unknown',form='unformatted',
     &      access='direct',recl=j)
         write(10,rec=imem,err=200,iostat=iostat)io_ficem,io_hicem,
     &      io_hsnwm,io_ticem,io_tsrfm   !,io_iceU,io_iceV
200      continue
         close(10)
         if (iostat==0) then
            !write(lp,*) '  ICE done'
            save_restart_ice=.true.
         else
            write(lp,'(a,i4)')
     &         'An Error encountered while writing record ',imem
            save_restart_ice=.false.
            write(lp,'(a,i4)') 'IOSTAT is :',iostat
         end if
      end if

      ! KAL : Should communicate result to all processes in MPI...
#if defined(MPI)
      call mpi_bcast(save_restart_ice,1,MPI_LOGICAL, 0,
     &               MPI_COMM_WORLD,ierr)
      call mpi_bcast(iostat,1,MPI_INTEGER, 0,  MPI_COMM_WORLD,ierr)
#elif defined (SHMEM)
      isaveice=0
      if (save_restart_ice) isaveice=1
      call SHMEM_INT8_AND_TO_ALL(isaveiceall,isaveice,1,0,0,N$PES, 
     &                           PWRK, PSYNC)
      save_restart_ice=isaveiceall==1
      call shmem_barrier_all()
      call shmem_int8_sum_to_all(alliostat,iostat,0,0,N$PES,PWRK,PSYNC)
      iostat=alliostat
#endif

      ! KAL - new: bail out on save errors.. Prevent disaster in the
      ! making
      if (iostat/=0) then
         if (mnproc==1) then
            write(lp,*) 'Exiting on iostat save  errors...'
            write(lp,*) 'IOSTAT is ',iostat
         end if
         call xcstop('(save_restart_ice)')
         stop '(save_restart_ice)'
      end if

      end function save_restart_ice





#if defined (WEEKLY_AVERAGE)
      function save_restart_ave(rungen,rt,imem)
      use mod_xc
      use mod_hycomfuncs
      use mod_year_info
      use mod_average
      implicit none

      logical :: save_restart_ave
      type(year_info),  intent(in) :: rt
      character(len=3), intent(in) :: rungen
      integer,          intent(in) :: imem

      !type(average_state_global) :: gave_week

      character(len=27) avename
      character(len=60) tmpchar
      integer i,j,find
      integer ierr
      
#if defined(MPI)
      include 'mpif.h'
#endif

      ! Dump with new version
      tmpchar=fileaveres(rungen,rt)
      if (mnproc==1) 
     &   write(lp,'(a)',advance='no') 'save_restart_ave: saving:'//
     &      trim(tmpchar)//'.[ab]'

      find=index(tmpchar,'.',back=.true.)
      call save_ave2(tmpchar(1:find-1),ave_week)
      if (mnproc==1)
     &    write(lp,'(a,i5,a)') '   counter= ',ave_week%counter, 
     &      '   AVE done'
         
      ! Always true .. for now ..
      save_restart_ave=.true.
      end function save_restart_ave
#endif


! Routine to save random forcing fields
! 22.11.2004 - KAL: Changed to use HYCOM type io

      logical function save_restart_rand(rungen,rt,ran,imem,
     &   nstep0,dtime0)
      use mod_forcing_nersc
      use mod_year_info
      use mod_za
      implicit none
      character(len=3)    , intent(in) :: rungen
      type(year_info)     , intent(in) :: rt
      type(forcing_fields), intent(inout) :: ran
      integer,              intent(in) :: imem
      integer,              intent(in) :: nstep0
      real*8,               intent(in) :: dtime0
      character(len=50) :: filename
      integer :: j,iostat,ios
      integer ioproc, ierr
      
#if defined(MPI)
      include 'mpif.h'
#endif

      integer :: nop,nfields,offset,k,l_nstep,l_imem,offsetb
      real :: xmin,xmax, l_time
      character(len=50) :: filebase
      character(len= 1) :: char1,carret
      character(len= 8) :: char8

      include 'common_blocks.h'


#if defined (SINGLE_RESTART) /* We always dump when this flag is set */
      save_restart_rand=save_restart_rand_ab(rungen,rt,ran,imem)
#if defined (SINGLE_RESTART_ONLY) /* Short-circuit when this flag is set */
      return
#endif /* SINGLE_RESTART_ONLY */
#endif /* SINGLE_RESTART */

      ! New output format -- uses hycom-style IO (zaio)
      nop=710
      filebase=rungen//'restRAND'//rt%cyy//'_'//rt%cdd//'_'//rt%chh
      if (mnproc==1) then 
         write(lp,'(a,i4)',advance='yes') 'save_restart_rand: saving '//
     &      trim(filebase)//'.[ab] record=',imem
      end if
      call zaiopf(trim(filebase)//'.a','unknown',nop)
      if (mnproc==1) then
         open(nop,file=trim(filebase)//'.b',status='unknown',
     &        form='formatted',access='direct',recl=80)
      end if


      ! Number of fields for one member 
      nfields=10

      ! Offset for IO
      offset=(imem-1)*nfields
      offsetb=offset+1 ! First record is model info 
      if (mnproc==1) then
         write(nop,118,rec=1) itdm,jtdm,kdm
         call flush(nop)
      end if

      ! Skip offset number of records
      do k=1,offset
         call zaiosk(nop)
      end do


      ! The rest ...
      call zaiowr(ran%slp(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+1)
     &      'slp     ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%taux(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+2)
     &      'taux    ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%tauy(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+3)
     &      'tauy    ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%wndspd(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+4)
     &      'wndspd  ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%airtmp(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+5)
     &      'airtmp  ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%relhum(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+6)
     &      'relhum  ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%clouds(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+7)
     &      'clouds  ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%precip(1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+8)
     &      'precip  ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%sss   (1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+9)
     &      'sss     ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiowr(ran%sst   (1-nbdy,1-nbdy),ip,.false.,
     &            xmin,xmax, nop, .false.)
      if (mnproc==1) then
         write(nop,117,rec=offsetb+10)
     &      'sst     ',nstep0,dtime0,imem,xmin,xmax
         call flush(nop)
      end if

      call zaiocl(nop)
      if (mnproc==1)  close(nop)
      save_restart_rand=.true.

      return

 117  format (a8,' =',i11,f11.2,i3,1p2e16.7)
 118  format ('itdm= ',i5, ' jtdm= ',i5,' kdm= ',i5)


      end function save_restart_rand



#if defined (SINGLE_RESTART)
      ! KAL - with the advent of ice types, fy fraction and ice age, 
      ! KAL - it is finally time to drop the old restart file setup for 
      ! KAL - ice variables. The .ab files are much more flexible for
      ! KAL - future changes. This version appends everything 
      ! KAL - to the HYCOM restart file
      logical function save_restart_ice_ab(rungen,rt,imem)
      use mod_xc
      use mod_za
      use mod_common_ice
      use mod_year_info
      use mod_hycomfuncs
      implicit none
      type(year_info) , intent(in) :: rt
      character(len=3), intent(in) :: rungen
      integer, intent(in) :: imem

      character filebase*80
      integer :: ierr

      save_restart_ice_ab=.true.
      filebase=trim(restartfile(rungen,rt,imem))
      if (mnproc==1) then
         write(lp,'(a)',advance='yes')'save_restart_ice_ab:  Saving  '//
     &      trim(filebase)//'.[ab]'
      end if
      call restart_appendfield(filebase,'ficem   ',0,1,ficem,ierr)
      call restart_appendfield(filebase,'hicem   ',0,1,hicem,ierr)
      call restart_appendfield(filebase,'hsnwm   ',0,1,hsnwm,ierr)
      call restart_appendfield(filebase,'tsrfm   ',0,1,tsrfm,ierr)
      call restart_appendfield(filebase,'ticem   ',0,1,ticem,ierr)
#if defined (TEST_ICE_AGE)
      call restart_appendfield(filebase,'fy_frac ',0,1,fy_frac,ierr)
      call restart_appendfield(filebase,'fy_age  ',0,1,fy_age ,ierr)
      call restart_appendfield(filebase,'rdg_frac',0,1,rdg_frac,ierr)
#endif
      if (ierr/=0) save_restart_ice_ab=.false.
      end function save_restart_ice_ab



      logical function save_restart_rand_ab(rungen,rt,ran,imem)
      use mod_forcing_nersc
      use mod_year_info
      use mod_za
      use mod_hycomfuncs
      implicit none
      character(len=3)    , intent(in) :: rungen
      type(year_info)     , intent(in) :: rt
      type(forcing_fields), intent(inout) :: ran
      integer,              intent(in) :: imem

      character(len=80)   :: filebase
      integer             :: ierr

      save_restart_rand_ab=.true.
      filebase=trim(restartfile(rungen,rt,imem))
      if (mnproc==1) then 
         write(lp,'(a)',advance='yes') 
     &      'save_restart_rand_ab: saving '//
     &      trim(filebase)//'.[ab]'
      end if
      call restart_appendfield(filebase,'rf_slp  ',0,1,ran%slp   ,ierr)
      call restart_appendfield(filebase,'rf_taux ',0,1,ran%taux  ,ierr)
      call restart_appendfield(filebase,'rf_tauy ',0,1,ran%tauy  ,ierr)
      call restart_appendfield(filebase,'rf_wspd ',0,1,ran%wndspd,ierr)
      call restart_appendfield(filebase,'rf_tair ',0,1,ran%airtmp,ierr)
      call restart_appendfield(filebase,'rf_rlhum',0,1,ran%relhum,ierr)
      call restart_appendfield(filebase,'rf_cloud',0,1,ran%clouds,ierr)
      call restart_appendfield(filebase,'rf_prcp ',0,1,ran%precip,ierr)
      call restart_appendfield(filebase,'rf_sss  ',0,1,ran%sss   ,ierr)
      call restart_appendfield(filebase,'rf_sst  ',0,1,ran%sst   ,ierr)
      if (ierr/=0) save_restart_rand_ab=.false.
      end function save_restart_rand_ab
#endif /* SINGLE_RESTART */





      end module m_save_restart
