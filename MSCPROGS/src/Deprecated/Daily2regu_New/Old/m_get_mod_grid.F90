module m_get_mod_grid
contains 
subroutine get_mod_grid(modlon,modlat,depths)
   use mod_dimensions
   implicit none
   real, dimension(nx,ny), intent(out) :: modlon,modlat,depths

   character(len=7) tag7    
   logical ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read position from model files
   inquire(file='newpos.uf',exist=ex)
   if (.not.ex) stop 'newpos.uf file does not exist'
   open(10,file='newpos.uf',form='unformatted',status='old')
     read(10)modlat,modlon
   close(10)
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read depths from model files
   write (tag7,'(i3.3,a,i3.3)') nx,'x',ny
   inquire(file='depths'//tag7//'.uf',exist=ex)
   if (.not.ex) stop 'depths.uf file does not exist'
   open(10,file='depths'//tag7//'.uf',form='unformatted',status='old')
      read(10)depths
   close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine  get_mod_grid
end module  m_get_mod_grid
