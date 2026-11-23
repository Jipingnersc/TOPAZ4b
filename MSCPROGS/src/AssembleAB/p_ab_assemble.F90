program ab_assemble
! KAL -- to its own file named "analysis_procXXX.[ab]
! KAL --
! KAL -- This routine will gather the
! KAL -- analysis from the separate analyzed files into one complete restart file.
! KAL -- The final files produced by this routine are named "analysisXXX.[ab]".
! KAL -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
! KAL -- Input arguments:
! KAL --     number of files used in analysis
! KAL -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

   use mod_raw_io
   use m_parse_blkdat
   use m_put_mod_fld
   implicit none

   integer*4, external :: iargc

   character(len=80) :: template,cetemplate ! restart template

   character(len=80) :: afile,cfile
   integer          :: fnd, iafile
   logical          :: ex, allok, nomatch
   character(len=8) :: cfld, ctmp
   character(len=3) :: cproc
   integer          :: tlevel, vlevel, nproc
   real             :: bmin, bmax, rdummy
   integer          :: idm,jdm,kdm
   integer          :: rstind,rstindi,rstindj 
   real, allocatable:: fld(:,:)
   real*4, allocatable:: fldr4(:,:)
   real*4 :: spval,amin,amax


   integer,parameter :: numfields=2
   integer :: ios,ios2, reclICE,ifld
   character(len=8) :: fieldnames(numfields)
   integer :: fieldlevels(numfields)


   if (iargc()==1) then
      call getarg(1,ctmp)
      read(ctmp,*) nproc
   else
      print *,'usage: assemble_process nproc'
      call exit(1)
   endif

   ! Get dimensions from blkdat
   call parse_blkdat('idm   ','integer',rdummy,idm)
   call parse_blkdat('jdm   ','integer',rdummy,jdm)
   call parse_blkdat('kdm   ','integer',rdummy,kdm)

   if (idm>0 .and. idm < 1e4 .and. jdm>0 .and. jdm<1e4) then
      allocate(fld  (idm,jdm))
      allocate(fldr4(idm,jdm))
   else
      print *,'fld allocate error'
      stop '(assemble_postprocess)'
   end if


   ! Inquire for existence
!   inquire(exist=ex,file=template(1:fnd-1)//'.b')
!   if (.not.ex) then
!      write(*,*) 'Can not find '//template(1:fnd-1)//'.b'
!      stop '(assemble_postprocess)'
!   end if
   ! Loop over restart file

   ! KAL -- need list of analysis files here, for now we hardcode
    cfile='analysis'
    nomatch=.false.
    rstind=1
    do iafile=1,nproc ! List of procs used in analysis
       write(cproc,'(i3.3)') iafile-1
       afile='analysis_proc'//cproc
       inquire(exist=ex,file=trim(afile)//'.b')
       if (.not.ex) then
         write(*,*) 'Can not find '//trim(afile)//'.b'
         stop '(assemble_postprocess)'
       end if
       print *, trim(afile)
       allok=.true.
       rstindi=1
       do while (allok)
         ! Get header info from each .b file 
         call rst_header_from_index(trim(afile)//'.b', cfld,vlevel,tlevel,rstindi,bmin,bmax,.true.)
         if (nomatch) then
            print '(a8," .. time, layer:",3i5," - for test")',cfld,tlevel,vlevel,rstindi
         endif
         if (tlevel/=-1) then
            if (iafile>1) then
              call rst_index_from_header(trim(cfile)//'.b', cfld,vlevel,tlevel,rstindj,bmin,bmax,.false.)
              if (rstindj==-1) then
                ! Read field from analysed file 
                call READRAW(fldr4,amin,amax,idm,jdm,.false.,spval,trim(afile)//'.a',rstindi)
                fld=fldr4

                print '(a8," .. time, layer:",3i5," - involved from file ",i4)',cfld,tlevel,vlevel,rstindi,iafile
                call put_mod_fld(trim(cfile),fld,1,cfld,vlevel,tlevel,rstind,idm,jdm)
                rstind=rstind+1 
              endif       
            else
          !    if (rstind==1) then
          !      call system("head -n2 "//trim(afile)//".b > "//trim(cfile)//".b")
          !    endif
              ! Read field from analysed file 
              call READRAW(fldr4,amin,amax,idm,jdm,.false.,spval,trim(afile)//'.a',rstindi)
              fld=fldr4
              print '(a8," -- time, layer:",2i4," - convered from file ", i4)',cfld,tlevel,vlevel,iafile
              call put_mod_fld(trim(cfile),fld,1,cfld,vlevel,tlevel,rstind,idm,jdm)
              rstind=rstind+1 
            endif
            rstindi=rstindi+1
         else
            allok=.false.
         endif

      end do

   end do
   print *,'Normal exit of assemble_postprocess'
#if defined(AIX)
   call exit_(0)
#else
   call exit(0)
#endif



end program
