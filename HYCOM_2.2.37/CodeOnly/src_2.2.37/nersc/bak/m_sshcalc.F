      module m_sshcalc
      ! computes ssh from micom variables
      contains
      subroutine sshcalc(m,ssh,slpin)
      use mod_xc
      use mod_forcing_nersc
      implicit none
      real, intent(out) :: ssh  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real, intent( in) :: slpin(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) 
      integer, intent(in) :: m
      integer i,j,l

      include 'common_blocks.h'

!$OMP PARALLEL DO PRIVATE(j,l,i  ) 
!$OMP&SCHEDULE(STATIC,jblk)
      do 81 j=1-margin,jj+margin
      do 81 l=1,isp(j)
      do 81 i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)
         ssh(i,j)=(montg1(i,j)/thref+pbavg(i,j,m))/onem + 
     &            (slp0-slpin(i,j))/100.0
81    continue
      end subroutine sshcalc
      end module m_sshcalc
