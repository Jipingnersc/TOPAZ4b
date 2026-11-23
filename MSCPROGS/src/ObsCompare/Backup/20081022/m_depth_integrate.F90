module m_depth_integrate

contains

   ! Routine to integrate with depth, upper interface given by "upper, and lower
   ! interface by "lower". surface is z=0, bottom z=depth.
   ! Ie positive direction downwards.

   ! Dimensions are arbitrary, as long as upper, lower and interfaces 
   ! have the same dimensions

   subroutine depth_integrate(upper,lower,int_val,int_dp,values,interfaces,idm,jdm,kdm,undef)
      implicit none
      integer, intent(in) :: idm,jdm,kdm
      real,    intent(in) :: undef
      real,    intent(in) :: upper        ! Upper interface for integration  (>0)
      real,    intent(in) :: lower        ! lower interface for integration  (>0,upper>)
      real, dimension(idm,jdm,kdm  ), intent( in) :: values     ! Variable to integrate
      real, dimension(idm,jdm,kdm+1), intent( in) :: interfaces ! Lower layer interfaces
      real, dimension(idm,jdm      ), intent(out) :: int_val    ! Integrated variable
      real, dimension(idm,jdm      ), intent(out) :: int_dp     ! Integrated thickness

      integer :: i,j,k
      real    :: up_int,lw_int
      integer :: itest,jtest
      real, parameter :: epsiloon=1e-3

      ! Error check
      if (lower < upper .or. upper<0 ) then
         print *,'Error -- inconsistent upper or lower boundary'
         print *,'upper:', upper
         print *,'lower:', lower
         stop '(depth_integrate)'
      end if



      !print *,'max p(2):',maxval(interfaces(:,:,2))
      !print *,'max p(3):',maxval(interfaces(:,:,3))
      int_val=0.
      int_dp =0.
      itest=20
      jtest=20
      do k=1,kdm
      do j=1,jdm
      do i=1,idm

         up_int= max( min(interfaces(i,j,k  ),lower) , upper)
         lw_int= min( max(interfaces(i,j,k+1),upper) , lower) 

         int_val(i,j)=int_val(i,j) + (lw_int-up_int)*values(i,j,k)
         int_dp (i,j)=int_dp (i,j) + (lw_int-up_int)

         !if (i==itest.and.j==jtest) then
         !   print *,maxval(interfaces(:,:,k))
         !   print '(7f10.2)',up_int,upper,lw_int,lower,interfaces(i,j,k), &
         !      interfaces(i,j,k+1),lw_int-up_int,int_dp(i,j)
         !end if
      end do
      end do
      end do
         
      do j=1,jdm
      do i=1,idm
         if ( int_dp(i,j)  > (lower-upper)*0.6 .and.  &
              lower-upper>epsiloon ) then
            int_val(i,j) = int_val(i,j) / int_dp(i,j)
         else
            int_dp (i,j)=undef
            int_val(i,j)=undef
         end if
         !if (i==itest.and.j==jtest)print *,int_dp(i,j),int_val(i,j),lower-upper,undef
         !print *,int_dp(i,j),int_val(i,j),lower-upper,undef
      end do
      end do

   end subroutine depth_integrate
end module m_depth_integrate

