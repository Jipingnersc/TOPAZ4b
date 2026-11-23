!Calanus mortality dynamics model
!-------------------------------------------------------------------------
subroutine calmort(gen,day)
!-------------------------------------------------------------------------
   use mod_CAL06_calnusparticles ! include      'comcal.txt'
   implicit none
!
   integer, intent(in) :: gen, day
   integer ::	 i,intern,j,x
   real    ::     rand,pred! ab(0:13),zsum(0:13)
   integer ::     stagesum(0:13) 
!
   open(unit=811,file='mortality_out',status='old', &
             form='formatted',position='append')
 119 format(a12,i12,a8,i8,a45) 
  
   stagesum=0;
!Calculate mortality of calanusA								
   do i = 1,pop
      x = int(cfs(i)%av(1))
      stagesum(x)=stagesum(x)+1
      if(x>2) then 
         if(cfs(i)%av(3)<(st_wgt(x)*0.9)) then 
!AS 14082008            cfs(i)%av(0) = 0.	 !Kill individuals with too low structural weight
            cfs(i)%av(2) = cfs(i)%av(2)*exp(-0.1/24.0)	 !Impose a 10% mortality on individuals with low structural weight
!
write(811,119) 'Individual #', cfs(i)%num,'on day', day,&!
' killed - too low strucural weight'!
            zm_cal(x,5) = zm_cal(x,5) + cfs(i)%av(2)*(1.0-exp(-0.1/24.0))
!!            zn_cal(x,5) = zn_cal(x,5) + 1.
            zb_cal(x,5) = zb_cal(x,5) + cfs(i)%av(2)*(cfs(i)%av(3)+cfs(i)%av(5))*(1.0-exp(-0.1/24.0))
            goto 1
         endif 
      endif
      if(cfs(i)%av(2)<1.) then
         cfs(i)%av(0) = 0.
!
write(811,119) 'Individual #', cfs(i)%num, 'on day', day, &!
' killed - number of individuals less than 1'!
!
         zm_cal(x,7) = zm_cal(x,7) + cfs(i)%av(2)
         zn_cal(x,7) = zn_cal(x,7) + 1.
         zb_cal(x,7) = zb_cal(x,7) + cfs(i)%av(2)*(cfs(i)%av(3)+cfs(i)%av(5))
      endif				 
!CAS06      if(cfs(i)%av(4)>400.0) then
      if(cfs(i)%age>400.0) then
         cfs(i)%av(0) = 0. 
!
write(811,119) 'Individual #', cfs(i)%num,'on day', day, &!
' killed - age more than 400 days'!
!
         zm_cal(x,9) = zm_cal(x,9) + cfs(i)%av(2)
         zn_cal(x,9) = zn_cal(x,9) + 1.
         zb_cal(x,9) = zb_cal(x,9) + cfs(i)%av(2)*(cfs(i)%av(3)+cfs(i)%av(5))
      endif
      if(cfs(i)%av(6)>800.0) then
         cfs(i)%av(0) = 0.	!Kill worn out adults..
!
write(811,119) 'Individual #', cfs(i)%num,'on day', day, &!
' killed - worn out adults'!
!
         zm_cal(x,8) = zm_cal(x,8) + cfs(i)%av(2)
         zn_cal(x,8) = zn_cal(x,8) + 1.
         zb_cal(x,8) = zb_cal(x,8) + cfs(i)%av(2)*(cfs(i)%av(3)+cfs(i)%av(5))
      endif
      if (cfs(i)%av(0)<=0.) goto 1
      pred = 0.

! AS25012008 set constant predation mortality also for the older stages
! From Aksnes & Blindheim, Ohman & al 2004
! Mortalirty rate from From Ohman 1994
! Function modified with depth : AS 30012008, but minimum 0.001 (diapause_mort)  at depth      
      if (cfs(i)%av(11).ne.1.0) then
        pred = diapause_mort/24.0
      else
        pred = mort(x)/24.0
      end if

      zm_cal(x,10) = zm_cal(x,10) + cfs(i)%av(2)*(1.-exp(-pred)) 
 !!     zn_cal(x,10) = zn_cal(x,10) + 1.*(1.-exp(-pred))
      zb_cal(x,10) = zb_cal(x,10) + cfs(i)%av(2)*(cfs(i)%av(3)+cfs(i)%av(5))*(1.-exp(-pred)) 
      pred = pred + cfs(i)%prisk(4)

      if(cfs(i)%av(2)>mthresh) then            !Kill individuals by reduction of their internal number 
         zm_cal(x,1) = zm_cal(x,1) + cfs(i)%av(2)*(1.-exp(-cfs(i)%prisk(1)))
 !!        zn_cal(x,1) = zn_cal(x,1) + 1.*(1.-exp(-cfs(i)%prisk(1)))
         zb_cal(x,1) = zb_cal(x,1) + cfs(i)%av(2)*(1.-exp(-cfs(i)%prisk(1)))*&
                             (cfs(i)%av(3)+cfs(i)%av(5))

         zm_cal(x,2) = zm_cal(x,2) + cfs(i)%av(2)*(1.-exp(-cfs(i)%prisk(2)))
 !!        zn_cal(x,2) = zn_cal(x,2) + 1.*(1.-exp(-cfs(i)%prisk(2)))
         zb_cal(x,2) = zb_cal(x,2) + cfs(i)%av(2)*(1.-exp(-cfs(i)%prisk(2)))*&
                             (cfs(i)%av(3)+cfs(i)%av(5))

         zm_cal(x,3) = zm_cal(x,3) + cfs(i)%av(2)*(1.-exp(-cfs(i)%prisk(3)))
 !!        zn_cal(x,3) = zn_cal(x,3) + 1.*(1.-exp(-cfs(i)%prisk(3)))
         zb_cal(x,3) = zb_cal(x,3) + cfs(i)%av(2)*(1.-exp(-cfs(i)%prisk(3)))*&
                             (cfs(i)%av(3)+cfs(i)%av(5))

         zm_cal(x,4) = zm_cal(x,4) + cfs(i)%av(2)*(1.-exp(-cfs(i)%prisk(4)))
 !!        zn_cal(x,4) = zn_cal(x,4) + 1.*(1.-exp(-cfs(i)%prisk(4)))
         zb_cal(x,4) = zb_cal(x,4) + cfs(i)%av(2)*(1.-exp(-cfs(i)%prisk(4)))*&
                             (cfs(i)%av(3)+cfs(i)%av(5))

         cfs(i)%av(2) = cfs(i)%av(2)*exp(-pred)       !Predation mortality 

      else         !Kill individuals by Monte-Carlo simulations          
!
         intern = int(cfs(i)%av(2))
         do j = 1,intern	                 !Remove individuals
            call random_number(rand)
            if(rand<(1.-1.*exp(-pred))) then
               cfs(i)%av(2) =cfs(i)%av(2) - 1.0  
            endif	 
         enddo     
      endif
      cfs(i)%prisk(1) = 0.
      cfs(i)%prisk(2) = 0.
      cfs(i)%prisk(3) = 0.
      cfs(i)%prisk(4) = 0.

!
1   enddo
!

!do x=0,13
!print*,'BGM', zm_cal(x,10), zn_cal(x,10), zb_cal(x,10), x, stagesum(x)
!end do

   close(811)
end subroutine calmort
!
!
!Code for mortality cause z
!1 = herring predation
!2 = mackerel predation
!3 = mesopelagic predation
!4 = tactile predation
!5 = starvation
!7 = internal number less than 1
!8 = reproduced more than 800 eggs
!9 = age > 400
!10 = unspecified mortality - standard rates

!updates:
!AS25012008 set constant predation mortality also for the older stages
