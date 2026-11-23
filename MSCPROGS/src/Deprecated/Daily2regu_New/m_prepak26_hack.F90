module m_prepak26_hack
contains
subroutine prepak26_hack(rungen,fname,n2d,n3d,varlist,maxfld,idat,rt,idm,jdm,kdm,sigma)
!  Routine to do preprocessing of data packing
   !use mod_dimensions
   !use mod_common
   use mod_year_info
   implicit none
   integer, parameter :: maxhdr=1000


   integer :: maxfld

   character(len=*), intent(in) :: rungen
   character(len=*), intent(in) :: fname
   character(len=*), intent(in) :: varlist(maxfld)
   type(year_info), intent(in)  :: rt
   integer, intent(in) :: idm,jdm,kdm
   real, intent(in) :: sigma(kdm)

   integer :: n2d, n3d
   integer :: iw,it,nn,in,iy,n,icc,inc,loop

   integer, parameter :: xpivn=-999
   real, parameter :: gridn =-999.


   integer lines

   character*80 whole(maxhdr)

   integer idat(3)
   INTEGER count

!  Initialise Variables for file header  ???
   whole(:)=" "

!  Create file Header
      !nsday=(86400./baclin)+0.0001

      WRITE(whole(2),'(a,a3)')'Run: ',rungen 
!      WRITE(whole(3),'(a,a4,a,a2,a,a3,a,a2,a,a4)')'y=',rt%cyy,' m=',rt%cmm,' d=',rt%cdd,' h=',rt%chh,' s=',rt%css
      WRITE(whole(3),'(a)')fname
      WRITE(whole(4),103) n2d
      DO count=5,n2d+4
         WRITE (whole(count),'(a5)') varlist(count-4)
      enddo

      WRITE(whole(count),105) n3d
      DO loop=count+1,count+n3d
         WRITE(whole(loop),'(a5)') varlist(loop-count+n2d)
      enddo
      WRITE(whole(loop),106) IDM,JDM,KDM
      WRITE(whole(loop+1),107)

      inc=loop+2
      icc=1
      DO n=1,kdm,5
         in=n+4
         IF (in.GT.KDM) in=KDM

         DO nn=n,in
            WRITE(whole(inc)(icc:icc+10),108) nn,sigma(nn)
            icc=icc+11
         enddo

         inc=inc+1
         icc=1
      enddo

      WRITE(whole(inc),109) gridn,xpivn
      WRITE(whole(inc+1),113)
      WRITE(whole(inc+2),114)
      WRITE(whole(inc+3),115)
      WRITE(whole(inc+4),116)
      DO iw=maxhdr,2,-1
         IF (whole(iw).EQ." ") lines=iw
      enddo

      WRITE(whole(1),'(a3,i4)')'1.1',lines

 103  FORMAT("2D fields: ",I2)
 105  FORMAT("3D fields: ",I2)
 106  FORMAT("Domain: I=",I4," J=",I4," K=",I4)
 107  FORMAT("Layer Sigma-Theta:")
 108  FORMAT(I2,"=",F6.3)
 109  FORMAT("Resolution (gridn)= ",F6.2," Equator is at xpivn=",F6.2) 
 113  FORMAT("Matrix written K(J(I)),all 2D and Layer 1 fields first")
 114  FORMAT("Field separator is '++'.")
 115  FORMAT("Field header format 'Field Name,Layer,x,y,a,b,c,ASCII"," string length'")
 116  FORMAT("in fortran format('++ ',A5,3I2,3I4,I5) and data ","format (40(A2)).")



   open(10,file=trim(fname)//'.hdr')
      do it=1,lines
         iy=index(whole(it),'        ')
         write(10,'(a)') whole(it)(1:iy)
      enddo
   close(10)

end subroutine prepak26_hack
end module m_prepak26_hack
