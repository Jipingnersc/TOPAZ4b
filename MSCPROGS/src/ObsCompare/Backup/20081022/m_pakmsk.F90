module m_pakmsk
contains
subroutine pakmsk(mask,array,work,ii,jj,util,length)
! send 'array' to subr. pakk after deleting array elements outside 'mask'
   real array(ii,jj),work(ii,jj)
   integer mask(ii,jj)
   character*2 util(ii*jj+14)

   do j=1,jj
   do i=1,ii
      work(i,j)=0.
     if (mask(i,j).ne.0) work(i,j)=array(i,j)
   enddo
   enddo
   call pakk(work,ii,jj,util,length)
end subroutine
end module m_pakmsk
