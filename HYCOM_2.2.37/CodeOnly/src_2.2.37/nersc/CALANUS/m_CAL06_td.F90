!---------------------------------------------------------------------------------

subroutine td(gen,cmort)

!---------------------------------------------------------------------------------
  implicit   none
  integer	   gen,s,m
  real	   cmort
  include    'comcal.h'

!Perform training decomposition of 
  if(gen<=24) then
     cmort = 1. 
  elseif(gen<=48) then
     cmort = 1.2
  elseif(gen<72) then
     cmort = 1.4
  else
     cmort = 1.6
  endif

!Set mortality to zero
  do s = 0,13
     do m = 1,10
        zm_cal(s,m) = 0.
        zn_cal(s,m) = 0.
        zb_cal(s,m) = 0.
     enddo
  enddo

end subroutine td
