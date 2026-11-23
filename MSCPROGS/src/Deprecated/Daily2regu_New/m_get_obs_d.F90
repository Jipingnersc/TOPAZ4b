module m_get_obs_d
contains
subroutine get_obs_d(obs,nrobs,string)
   use mod_measurement
   implicit none

   integer,           intent(in)  :: nrobs           ! Number of measurements
   type(measurement), intent(out) :: obs(nrobs)      ! measurements
   character(len=*),  intent(in)  :: string

   integer reclO,j
   logical :: ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read measurements and store in obs
   inquire(exist=ex,file='observations'//trim(string)//'.uf')
   if (.not.ex) then
      print *,'File "observations'//trim(string)//'.uf" does not exist'
      stop
   end if
   inquire(iolength=reclO)obs(1)
   open(10,file='observations'//trim(string)//'.uf',form='unformatted',access='direct',recl=reclO)
   do j=1,nrobs
      read(10,rec=j)obs(j)
   enddo
   close(10)

end subroutine get_obs_d
end module m_get_obs_d
