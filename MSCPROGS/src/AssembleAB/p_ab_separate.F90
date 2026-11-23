program ab_separate
! KAL -- to its own file named "TMPrst_2D.[ab] amd TMPrst_3D.[ab]
! KAL --
! KAL -- This routine will separate the source file into two files.
! KAL -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
! KAL -- Input arguments:
! KAL --     source file 
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
      call getarg(1,template)
   else
      print *,'usage: separate_process sourcefile'
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
      stop '(separate_process)'
   end if

   ! Remove postfix of template file
   fnd=max(index(template,'.a'),index(template,'.b'))

   ! Inquire for existence
   inquire(exist=ex,file=template(1:fnd-1)//'.b')
   if (.not.ex) then
      write(*,*) 'Can not find '//template(1:fnd-1)//'.b'
      stop '(separate_process)'
   end if

   ! KAL -- need list of analysis files here, for now we hardcode
   afile='TMPrst_2D'
   cfile='TMPrst_3D'

   nomatch=.false.
   allok=.true.
   rstind=1
   rstindi=1
   rstindj=1
   do while (allok)
     ! Get header info from each .b file 
     call rst_header_from_index(template(1:fnd-1)//'.b', cfld,vlevel,tlevel,rstind,bmin,bmax,.true.)
     if (nomatch) then
        print '(a8," .. time, layer:",3i5," - for source")',cfld,tlevel,vlevel,rstind
     endif
     if (tlevel/=-1) then
       ! Read field from source file 
       call READRAW(fldr4,amin,amax,idm,jdm,.false.,spval,template(1:fnd-1)//'.a',rstind)
       fld=fldr4
       if (vlevel==0) then
         print '(a8," .. time, layer:",3i5," - write into 2D afile")',cfld,tlevel,vlevel,rstindi
         call put_mod_fld(trim(afile),fld,1,cfld,vlevel,tlevel,rstindi,idm,jdm)
         rstindi=rstindi+1 
       elseif (vlevel>0.and.vlevel<=kdm) then
         print '(a8," .. time, layer:",3i5," - write into 3D afile")',cfld,tlevel,vlevel,rstindj
         call put_mod_fld(trim(cfile),fld,1,cfld,vlevel,tlevel,rstindj,idm,jdm)
         rstindj=rstindj+1 
       endif
       rstind=rstind+1
     else
       allok=.false.
     endif
   end do

   print *,'Normal exit of separate_postprocess'
#if defined(AIX)
   call exit_(0)
#else
   call exit(0)
#endif



end program
