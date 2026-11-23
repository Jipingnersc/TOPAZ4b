!Calanus mortality dynamics model
!-------------------------------------------------------------------------
subroutine vertmov(tid,gen,year,hour)
!subroutine vertmov(tid,gen,pop,cmort,year,hour,vertmov)
!
!-------------------------------------------------------------------------
!
   use mod_CAL06_calnusparticles
   implicit none
   integer, intent(in) :: tid, gen,year,hour
   integer 	   ::   d_cal,i,dag,y
   real            ::   nyd,bg,slig,rand,size,stay,visran
   real            ::   maxdepth,dyp,up,phyt,syn,dd,str
   real            ::   day_depth, night_depth
   real, parameter ::   dvm_speed=30.0 !AS: added (m/hour)
   real, parameter ::   cm1=0.7,cm2=0.278 
!   real, parameter ::   svm_upspeed=2.0, svm_downspeed=1.0 !seasonal vertical migration swimming speed
   real, parameter ::   svm_upspeed=0.833, svm_downspeed=0.833 !seasonal vertical migration swimming speed
   real :: annchl
!!!   include 'comcal.h'

!Set year to 2000 outside period 1992-1999 for the herring data
y = year
if((year>1999).or.(year<1992)) y = 2000
if((tid==1).and.(hour==1)) bh(1:12) = 0.

   open(unit=811,file='mortality_out',status='old', &
             form='formatted',position='append')
 119 format(a12,i12,a8,i8,a45) 

!Calculate vertical movement and environmental variables
   do i = 1,pop
      if(cfs(i)%av(0)==0.0) goto 1
      nyd = 0.
      stay = 0.
      up = 0.
      d_cal = int(cfs(i)%av(8))
      maxdepth = cfs(i)%wdepth
      bg = cfs(i)%av(10)
!Body length C&W 1998, convert from cm to m, divided by to to get radius of prey
      size = 0.001*cm1*(cfs(i)%av(3)+cfs(i)%av(5))**cm2   
!Find visual range at depth
      call surlig(hour,slig,tid,cfs(i)%mirr,dag)

       dd = 1.0 !exp(-cfs(i)%caldens_bio/card*1.)  !Density dependence in wet weight m-2

! chlorophyll based on artificial annual cycles
       if (idealized) then
         if (cfs(i)%wdepth.lt.400.0) then
           phyt=annchl(real(tid),1)*cchl*dd
         else
           phyt=annchl(real(tid),2)*cchl*dd
        end if

        if (cfs(i)%av(8)<200.0) then
           phyt=0.01*phyt
        elseif (cfs(i)%av(8)<150.0) then
           phyt=0.02*phyt
        elseif (cfs(i)%av(8)<100.0) then
           phyt=0.2*phyt
        elseif (cfs(i)%av(8)<50.0) then
           phyt=1.0*phyt
        endif
      else 
         phyt=cfs(i)%phy_conc*dd
      end if

!CAS
!Calculate movement for different diapause states
      if (cfs(i)%av(11)==0.0) then
         nyd = cfs(i)%av(8)
      elseif(cfs(i)%av(11)==1.0) then 
!Active movement of different stages 
         if(cfs(i)%av(1)<3.0) then         !<N3
            if(cfs(i)%av(8)<cfs(i)%mixl) then
               call random_number(rand)
               nyd = cfs(i)%av(8) + rand*2 - 1.
            else
               nyd = cfs(i)%av(8) + 1.      !Egg - is sinking
            endif
         elseif(cfs(i)%av(1)<9.0) then  !>N3 !CAS endret 10.6.2008: <C3
            if (idealized) then       
               nyd = min(20.0,cfs(i)%wdepth-5.0) 
            else   
               nyd = min(cfs(i)%depth_max_phy_conc,cfs(i)%wdepth-5.0,150.0) 
            end if
         else                   !Copepodites 
if (i==3) print*, 'Dag', dag, cfs(i)%av(14), ' Old depth: ', cfs(i)%av(8)
if (i==3) print*, 'Day:   ', cfs(i)%sv(5)+cfs(i)%sv(6)*size,cfs(i)%sv(5:6), size
if (i==3) print*, 'Night: ', cfs(i)%depth_max_phy_conc,cfs(i)%wdepth 
            day_depth=min(cfs(i)%sv(5) + cfs(i)%sv(6)*size,cfs(i)%wdepth-5.0)
            night_depth=min(day_depth,cfs(i)%depth_max_phy_conc)
            if(dag==0) then !Night
               if (idealized) then
                  nyd = min(30.0,cfs(i)%wdepth-5.0)
               else
!AS swim speed                 nyd = min(cfs(i)%depth_max_phy_conc,cfs(i)%wdepth-5.0)
                  if (cfs(i)%av(8) > night_depth) then
                    if (cfs(i)%av(8)-dvm_speed > night_depth) then
                      nyd = cfs(i)%av(8)-dvm_speed 
                    elseif (cfs(i)%av(8)-dvm_speed < night_depth) then
                      nyd = night_depth
                    endif
                  else
                    if (cfs(i)%av(8)+dvm_speed > night_depth) then
                      nyd = night_depth
                    elseif (cfs(i)%av(8)+dvm_speed < night_depth) then
                      nyd = cfs(i)%av(8)+dvm_speed 
                    endif
                  endif
                  nyd=min(nyd,cfs(i)%wdepth-5.0)
               end if 
            else           !Day
!AS swim speed               nyd = min(cfs(i)%sv(5) + cfs(i)%sv(6)*size,cfs(i)%wdepth-5.0)
               if (cfs(i)%av(8) > day_depth) then
                 if (cfs(i)%av(8)-dvm_speed > day_depth) then
                   nyd = cfs(i)%av(8)-dvm_speed 
                 elseif (cfs(i)%av(8)-dvm_speed < day_depth) then
                   nyd = day_depth
                 endif
               else
                 if (cfs(i)%av(8)+dvm_speed > day_depth) then
                   nyd = day_depth
                 elseif (cfs(i)%av(8)+dvm_speed < day_depth) then
                   nyd = cfs(i)%av(8)+dvm_speed 
                 endif
               endif
               nyd=min(nyd,cfs(i)%wdepth-5.0)

            endif    
if (i==3) Print*, cfs(i)%av(14), ' New depth: ', nyd, 'pos:', cfs(i)%av(8:10)
         endif !Stage boolean
 
      elseif (CFS(I)%av(11)==2.0) then       !Migrating down to diapause !Active movement of different stages 
!Swimming at 1 m hour-1 Carlotti & Wolf 1998:
          nyd = min(cfs(i)%av(8) + svm_downspeed,cfs(i)%wdepth-5.0)	          
          dyp = cfs(i)%sv(4)
          !dyp = cfs(i)%av(15)
          if (dyp>maxdepth) dyp = maxdepth
          if (nyd>=dyp) stay = 1           !Induce diapause, go from 2 to 0
      elseif (cfs(i)%av(11)==3.0) then	      !Migrating up from diapause
          if (cfs(i)%av(8)>cfs(i)%mixl*0.5) then 
             nyd = cfs(i)%av(8) - svm_upspeed
          else 
             nyd = cfs(i)%av(8) - svm_upspeed
             up = 1
          endif
       endif
      !Set depth
       if(nyd>maxdepth) then
          cfs(i)%av(8) = maxdepth
       elseif(nyd<1) then 
          cfs(i)%av(8) = 1.
       else
          cfs(i)%av(8) = nyd
       endif
 !Calculate phytoplankton biomass
       cfs(i)%av(12) = 0.
!       dd = 1.
! chlorophyll based on artificial annual cycles
       if (idealized) then
          if (cfs(i)%wdepth.lt.400.0) then
             cfs(i)%av(12)=annchl(real(tid),1)*cchl*dd
          else
             cfs(i)%av(12)=annchl(real(tid),2)*cchl*dd
          end if
          if (cfs(i)%av(8)>200.0) then
             cfs(i)%av(12)=0.01*cfs(i)%av(12)
          elseif (cfs(i)%av(8)>150.0) then
             cfs(i)%av(12)=0.02*cfs(i)%av(12)
          elseif (cfs(i)%av(8)>100.0) then
             cfs(i)%av(12)=0.2*cfs(i)%av(12)
          elseif (cfs(i)%av(8)>50.0) then
             cfs(i)%av(12)=1.0*cfs(i)%av(12)
         endif
      else
         cfs(i)%av(12)=cfs(i)%phy_conc*dd
      end if

if (i==1) then
!check phytoplankton concentration
  if (idealized) then
     print*, 'Phytoplankton concetration'!
     if (cfs(i)%wdepth.lt.400.0) then
        print*, cfs(i)%av(8),cfs(i)%av(12),annchl(tid,1),cchl*dd,cchl
     else
        print*, cfs(i)%phy_conc, cfs(i)%av(8),cfs(i)%av(12),annchl(tid,2),cchl*dd,cchl
     end if
  else
     print*, 'Phytoplankton concetration, depth:',cfs(i)%phy_conc,cfs(i)%av(8)
  end if
end if

!CAS

       d_cal = int(cfs(i)%av(8))
       if(cfs(i)%av(1)>8.0) then


          if(d_cal<=1000) then
             str = 0.00000141/((0.8+0.03*mesl)**2.)
             call getr(slig,d_cal,visran,size,diff_cal,str)
             visran = visran*0.01	 !Convert to m
!if (i==100) print*, 'Visran', visran,slig,d_cal,size,diff_cal,str
             syn =  max(pi_cal*(visran*0.5)**2.,0.00000000001)
!Calculate encounter rate with visual predators
             cfs(i)%prisk(3) = syn*mpvel*mesopred(d_cal,dag)*cmort/mesw*1.	       
          endif
       endif
!Calculate encounter rate with tactile predators  
       if (d_cal<=1000) cfs(i)%prisk(4) = pi_cal*tr*tr*tpvel*tacpred(d_cal,dag)*cmort/tacw*1.
! 
!Update diapause status
       if(stay==1.) cfs(i)%av(11) = 0.!
       if(up==1.) then
          cfs(i)%av(11) = 1.
          call random_number(rand)
!Kill the males of the population, assuming a sex ratio of 0.5
          if(cfs(i)%gender==0) then
             cfs(i)%av(2) = 0. 
             write(811,119) 'Individual #', cfs(i)%num, &
               'on day', tid,' killed - male .... 2'!
          endif 
       endif
     !  if(cfs(i)%av(8)>cfs(i)%sv(4)) cfs(i)%av(8) = cfs(i)%sv(4)
1  enddo
   close(811)
end subroutine vertmov

