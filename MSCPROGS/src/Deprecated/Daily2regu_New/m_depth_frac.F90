module m_depth_frac

contains

   ! Routine to calculate fraction of a layer which is positioned between two
   ! depth levels
   

   subroutine depth_frac(upper,lower,interfaces,dfrac,idm,jdm,kdm)
      implicit none
      integer, intent(in) :: idm,jdm,kdm
      real,    intent(in) :: upper  ! Upper interface for integration
      real,    intent(in) :: lower  ! lower interface for integration  (< upper)
      real, dimension(idm,jdm,kdm+1), intent( in) :: interfaces ! Lower layer interfaces
      real, dimension(idm,jdm,kdm), intent(out) :: dfrac ! Lower layer interfaces

      integer :: i,j,k
      real    :: up_int,lw_int,dphere
      integer :: itest,jtest
      real, parameter :: epsiloon=1e-3

      ! Error check
      if (lower < upper .or. upper<0 ) then
         print *,'Error -- inconsistent upper or lower boundary'
         print *,'upper:', upper
         print *,'lower:', lower
         stop '(depth_integrate)'
      end if



      itest=20
      jtest=20
      do k=1,kdm
      do j=1,jdm
      do i=1,idm

         up_int= max( min(interfaces(i,j,k  ),lower) , upper)
         lw_int= min( max(interfaces(i,j,k+1),upper) , lower) 

         dphere= interfaces(i,j,k+1)-interfaces(i,j,k)

         dfrac(i,j,k)= (lw_int-up_int) / ( dphere + epsiloon)

         !if (i==itest.and.j==jtest) then
         !   print *,maxval(interfaces(:,:,k))
         !   print '(7f10.2)',up_int,upper,lw_int,lower,interfaces(i,j,k), &
         !      interfaces(i,j,k+1),lw_int-up_int,int_dp(i,j)
         !end if
      end do
      end do
      end do
         
   end subroutine depth_frac
end module m_depth_frac

