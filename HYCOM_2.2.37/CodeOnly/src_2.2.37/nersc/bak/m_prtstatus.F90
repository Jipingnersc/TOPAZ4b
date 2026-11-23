module m_prtstatus
contains
subroutine prtstatus(rt,rungen,iday,iens,nstep)
   use mod_year_info
   use mod_xc
   implicit none
   real*8, external :: wtime
#if ! defined (MPI) && ! defined (SHMEM)
   integer, external :: omp_get_max_threads
#endif
   integer, intent(in) :: iday
   integer, intent(in) :: iens
   integer, intent(in) :: nstep
   type(year_info), intent(in) :: rt
   character(len=3), intent(in) :: rungen
   real :: par_time=0.0

   real*8, save :: time_last=0.0
   real*8, save :: time_now=0.0
   real*8, save :: tot_time=0.0
   real*8, save :: time_per_tstep=0.0
   logical, save :: first=.true.
   integer, save :: steps=0
   integer, save :: icpu = 1


   if (first) then
      first=.false.
#if ! defined (MPI) && ! defined (SHMEM)
      icpu = max(1,omp_get_max_threads())
#else
      icpu=1
#endif
      time_now=wtime()
      time_last=time_now
   else
      time_last=time_now
      time_now=wtime()
   endif

   steps=steps+1
   tot_time = tot_time + (time_now-time_last)
   time_per_tstep=tot_time/steps

   !print *,'prtstatus routine -- mnproc : ',mnproc,time_now,time_last
   if (mnproc==1) then
      open(10,file=rungen//'status',status='unknown')
      write(10,'(a,i4,a,i8,a,a4,a2,a3,a2,a2,a2,a4,3(a,f6.2))')&
            'iens=',iens,', nstep=',nstep,&
           &' --- y',rt%cyy,'_d',rt%cdd,'_h',rt%chh,'_s',rt%css,&
           &' --- CPU=',time_now-time_last,&
           &', Avg CPU=',time_per_tstep,&
           &', Tot CPU=',tot_time/3600.0
      close(10)

      write(lp,'(a,i4,a,i8,a,a4,a2,a3,a2,a2,a2,a4,3(a,f6.2))')&
         'iens=',iens,', nstep=',nstep,&
        &' --- y',rt%cyy,'_d',rt%cdd,'_h',rt%chh,'_s',rt%css,&
        &' --- CPU=',time_now-time_last,&
        &', Avg CPU=',time_per_tstep,&
        &', Tot CPU=',tot_time/3600.0
     call flush(lp)
   end if


end subroutine prtstatus
end module m_prtstatus
