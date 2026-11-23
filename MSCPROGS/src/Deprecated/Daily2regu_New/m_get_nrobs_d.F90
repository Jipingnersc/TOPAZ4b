module m_get_nrobs_d
contains
integer function get_nrobs_d(string)
! Reads the observations to be used for assimilation from the file
! observation.uf. Each element is of type measurement and will be 
! stored in the vector d defined below.

   use mod_measurement
   implicit none

   character(len=*), intent(in) :: string
   integer reclO             ! Record length for an observation type
   type(measurement) Obs     ! measurements

   integer j
   logical ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find nb of  measurements stored
   Print *, 'Examining file observations'//trim(string)//'.uf'
   inquire(file='observations'//trim(string)//'.uf',exist=ex)
   if (.not.ex) then
      print *,'File "observations'//trim(string)//'.uf" does not exist'
      get_nrobs_d=0
      return
   end if
   inquire(iolength=reclO)Obs
   open(10,file='observations'//trim(string)//'.uf',form='unformatted',access='direct',recl=reclO)
   do j=1,1000000
      read(10,rec=j,err=200)Obs
      !print *,obs%id
   enddo
   200 get_nrobs_d=j-1
   close(10)
end function get_nrobs_d
end module m_get_nrobs_d
