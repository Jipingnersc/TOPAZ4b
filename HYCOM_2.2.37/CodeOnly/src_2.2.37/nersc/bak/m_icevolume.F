      module m_icevolume
      contains
      real function icevolume()
      use mod_common_ice
      implicit none
      include 'common_blocks.h'

      integer i,j,l
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vol
      real*8 :: icevol8

      vol=0.0
!$OMP PARALLEL DO PRIVATE(j,l,i) 
!$OMP&SCHEDULE(STATIC,jblk)
      do j=1,jj
      do l=1,isp(j)
      do i=ifp(j,l),ilp(j,l)
         vol(i,j)=ficem(i,j)*hicem(i,j)*scpx(i,j)*scpy(i,j)
      end do
      end do
      end do
!$OMP END PARALLEL DO

      !Gather tot ice volume from tiles
      call xcsum(icevol8,vol,ip)
      icevolume=icevol8

      end function icevolume
      end module m_icevolume
