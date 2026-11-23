!-----------------------------------------------------
! surface light
subroutine surlig(hour,slig,tid,maxirr,dag)

!-----------------------------------------------------

   use mod_CAL06_calnusparticles!!!!include    'comcal.txt' 
   implicit   none
   integer, intent(in) ::  hour,tid
   real   , intent(in) ::  maxirr
   integer, intent(out)::  dag
   real   , intent(out)::  slig
   real    ::   peak,step,timestep
   real    ::	   delta,height,sunangle
   real    ::	   sst,twlight,bg,maxl

   bg = 66.
   twlight = 5.76
   timestep = 1.
   maxl = 0.
   peak = 0.
   step = hour						   
   slig = 0
   delta=0.3979*sin((pi_cal/180.)*(0.9856*(tid-80)+1.9171*(sin((pi_cal/180.)*(0.9856*tid))-0.98112)))		 
   sst=(15.*pi_cal/180.)*(step*24./24)				
   height=delta*sin(bg*pi_cal/180.)-sqrt(1.-delta**2.)*cos(bg*pi_cal/180.)*cos(sst)
   sunangle = asin(height)*57.295 

!Calculate surface light
   if(sunangle>0) then
      slig = maxirr*height + twlight
   elseif(sunangle>=-6.) then
      slig = ((twlight -.048)/6.)*(6.+sunangle)+.048 
   elseif(sunangle>=-12.) then
      slig = ((.048 - 1.15E-4)/6.)*(12.+sunangle)+1.15E-4
   elseif(sunangle>=-18) then
      slig = (((1.15E-4)-1.15E-5)/6.)*(18.+sunangle)+1.15E-5
   else 
      slig = 1.15E-5
   endif

!Set day
   if(slig<twlight) then
      dag = 0
   else
      dag = 1
   endif

!if(tid==160) write(27,*)tid,hour,dag

end subroutine surlig
