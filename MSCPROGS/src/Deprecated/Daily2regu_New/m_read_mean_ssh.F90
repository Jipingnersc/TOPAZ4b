 module m_read_mean_ssh
 contains
 subroutine read_mean_ssh(mean_ssh)
! use mod_meanssh
use mod_dimensions
   implicit none
real, intent(out):: mean_ssh(nx,ny)
logical lmeanssh
!Read in model SSH mean
  inquire(file='meanssh.uf',exist=lmeanssh)
  if (lmeanssh) then
  open (10,file='meanssh.uf', status='unknown',form='unformatted')
  read (10) mean_ssh
  close (10)
  else 
  print*,'m_read_meanssh: File not found'
  stop
  endif
print*,'m_read_mean_ssh: mean ssh data is read'       
end subroutine read_mean_ssh
end module m_read_mean_ssh
