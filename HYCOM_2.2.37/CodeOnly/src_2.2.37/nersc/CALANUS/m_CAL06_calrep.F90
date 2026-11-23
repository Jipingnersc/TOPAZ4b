!-----------------------------------------

!subroutine calrep(pops,sumdelta) 
subroutine calrep(nflt) 

!-----------------------------------------  
   use mod_CAL06_calnusparticles  !include     'comcal.txt'

implicit none

   integer     ind,i,newpop,l,r,comp,avk(mpop),x,nn_cal, ee_cal,e_cal,n_cal
   integer     j,delta,new,newsuper,spawners
   integer, intent(inout) :: nflt    
   real        rand,kum,kumf(mpop),lotto,cuminv,suminv
   real        invN(mpop),sumpop,sumremove(0:13),stageind(0:13),fitness,fit

!Calculate fitness  

   spawners = 0
   fitness = 0.
   do i = 1,pop
!Calculate reproduction
      cfs(i)%av(7) = 0.0
      cfs(i)%av(17) = 0.0
      if(cfs(i)%av(0)==0.0) goto 11
      if(cfs(i)%av(2)<1.0) goto 11
      if((cfs(i)%av(1)==13.0).and.(cfs(i)%av(3)>=cmat)) then
         if(cfs(i)%av(8)<mld) then
            !Need to be in mixed layer to spawn:  
            if((cfs(i)%av(5)/eggw*1.)>cs_cal) then                
               !Limit number of eggs laid  
               !cfs(i)%av(7) = min(cfs(i)%av(2)*cfs(i)%av(5)/eggw*1.,maxegg)
               !Limit number of eggs laid   
               cfs(i)%av(7) = cfs(i)%av(2)*cfs(i)%av(5)/eggw*1.!,maxegg)  
               cfs(i)%av(17) = 1.0 !superi
               spawners = spawners + 1
               fitness = fitness + cfs(i)%av(7)
            endif
         endif
      endif
11 continue 
   enddo  

goto 111
!Legg nye individ til eksisterende egg i samme rute
do i = 1,pop
 if(cfs(i)%av(7)>0.) then 
  e_cal = int(cfs(i)%av(9))
  n_cal = int(cfs(i)%av(10))
  do j = 1,pop
   if(i.ne.j) then
    ee_cal = int(cfs(j)%av(9))
    nn_cal = int(cfs(j)%av(10))
    r = cfs(j)%av(1)
    if((e_cal==ee_cal).and.(n_cal==nn_cal).and.(r==0)) then
     cfs(j)%av(2) = cfs(j)%av(2) + cfs(i)%av(7)
     cfs(i)%av(5) = max(0.,cfs(i)%av(5) - (cfs(i)%av(7)*eggw)/cfs(i)%av(2)*1.) !0. 
     cfs(i)%av(6) = cfs(i)%av(6) + cfs(i)%av(7)/(cfs(i)%av(2)+0.0001)
     cfs(i)%av(19) = cfs(i)%av(19) + cfs(i)%av(7)/cfs(i)%av(18)*1.
     cfs(i)%av(7) = 0.
     cfs(i)%av(17) = 0.
     goto 55
    endif
   endif
  enddo
 endif
55 continue
enddo
111 continue

!Count number of new super individuals
newsuper = 0
do i = 1,pop
 if(cfs(i)%av(7)>0.) newsuper = newsuper + 1
enddo

!Calculate inverse abundance of individuals
   suminv = 0.
   do i = 1,pop
      if(cfs(i)%av(2)>0.0) invN(i) = 1./(cfs(i)%av(2)+0.01)
      suminv = suminv + invN(i)
   enddo
   sumpop = 0.
   do i = 1,pop
      if(cfs(i)%av(2)>0.0) invN(i) = invN(i)/suminv*1.
      if(cfs(i)%av(0)==1.0) sumpop = sumpop + cfs(i)%av(2)
   enddo

!Perform routine for making room for new individuals
   new = min(newsuper,offsp)
   sumremove(0:13) = 0.
   delta = 0
   if((pop+new)>mpop) then 
      delta = pop + new - mpop
      do j = 1,delta+1
66 continue 
         call random_number(rand)
         lotto = rand								    
         cuminv = 0.
         do i = 1,pop
            cuminv = cuminv + invN(i) 
            if(lotto<=cuminv) then 
               if(cfs(i)%av(17)==0.0) then
                  sumremove(int(cfs(i)%av(1))) = sumremove(int(cfs(i)%av(1))) + cfs(i)%av(2)
       	          cfs(i)%av(0) = 0.0
                  goto 6
               else
	          goto 66
	       endif
            endif
         enddo
6     enddo
!!      call popzip(pops,pop)
   endif
   sumdelta = sumdelta + delta

!Calculate number of super inds in each stage
   do x = 0,13
      stageind(x)=0.001
   enddo
   do i = 1,pop
      x=int(cfs(i)%av(1))
      stageind(x) = stageind(x)+1
   enddo

!Divide the removed individuals among the survivors to maintain mass balance
   do i = 1,pop
      x=int(cfs(i)%av(1))
      cfs(i)%av(2) = cfs(i)%av(2) + sumremove(x)/stageind(x)*1. !if(av(0,i)==1.0) 
   enddo

!Sort parents according to their egg production
   do l = 1,pop
      r = 0
      do comp = 1,pop
         if(cfs(l)%av(7)<cfs(comp)%av(7)) r = r+1
         if(cfs(l)%av(7)==cfs(comp)%av(7)) then
            if(l>=comp) r = r+1
         endif
      enddo
      cfs(r)%fi = l
   enddo									  

!Enter spawning routine if there are spawners
   if(newsuper<=0) goto 9
   if(newsuper>offsp) then 
      kum = 0.
      do i = 1,spawners	  !Calculate cumulative fitness of parents
         avk(cfs(i)%fi) = 0
         kumf(cfs(i)%fi) = 0.
         fit = cfs(cfs(i)%fi)%av(7)/fitness*1.  !????????
         kum = kum + fit   
         kumf(cfs(i)%fi) = kum 
      enddo
      do j = 1,offsp
         call random_number(rand)
         lotto = rand
         do i = 1,spawners	  !Decide which individuals will reproduce
            if(lotto<kumf(cfs(i)%fi)) then
               avk(cfs(i)%fi) = avk(cfs(i)%fi) + 1
               goto 5
            endif  
         enddo
5 continue 
      enddo
   else
      do i = 1,spawners	  !Decide which individuals will reproduce
         avk(cfs(i)%fi) = cfs(cfs(i)%fi)%av(17)
      enddo
   endif

!Perform reproduction
   newpop = pop + 1
   do ind = 1,spawners 
 !Assign reduction of fat of parent for realised reproduction
      cfs(cfs(ind)%fi)%av(6) = cfs(cfs(ind)%fi)%av(6) + &
                               cfs(cfs(ind)%fi)%av(7)/(cfs(cfs(ind)%fi)%av(2)+0.0001)
      cfs(cfs(ind)%fi)%av(5) = max(0.,cfs(cfs(ind)%fi)%av(5) &
                             - cfs(cfs(ind)%fi)%av(7)*eggw/cfs(cfs(ind)%fi)%av(2)*1.) !0.
     cfs(cfs(ind)%fi)%av(19) = cfs(cfs(ind)%fi)%av(19) &
                             + cfs(cfs(ind)%fi)%av(7)/cfs(cfs(ind)%fi)%av(18)*1.
      do j = 1,1  !avk(fi(ind))
         i = newpop  
         if(i>mpop) then 
            pop = newpop-1
            goto 9
         endif
         do r = 1,maxs ! cas - only 6 attributes maxs
            cfs(i)%sv(r) = cfs(cfs(ind)%fi)%sv(r)
         enddo

         cfs(i)%av(0) = 1.
         cfs(i)%av(1) = 0.
         !Calculate internal number   NBNB:
         cfs(i)%av(2) = cfs(cfs(ind)%fi)%av(7)  !min(av(7,fi(ind))/avk(fi(ind))*1.,maxegg) !*(1.+sumremove/sumpop*1.)  
         cfs(i)%av(3) = 0.115  !in microgram C, reflects N3 (Carlotti and Wolf 1998) 
         cfs(i)%av(4) = 0.
         cfs(i)%av(5) = 0.115
         cfs(i)%av(6) = 0.
         cfs(i)%av(7) = 0.
         cfs(i)%av(8) = cfs(cfs(ind)%fi)%av(8)
         cfs(i)%av(9) = cfs(cfs(ind)%fi)%av(9)
         cfs(i)%av(10) = cfs(cfs(ind)%fi)%av(10)
         cfs(i)%av(11) = 1.  
         cfs(i)%av(12) = 0.  
         cfs(i)%av(13) = 0.  
         cfs(i)%av(14) = cfs(cfs(ind)%fi)%av(14) + real(cfs(cfs(ind)%fi)%nsip+1)*1000000.0 
         cfs(i)%av(18) = cfs(i)%av(2)
         cfs(i)%av(19) = 0.
         cfs(i)%nsip   = 100.  ! ti distinguish between generations
         cfs(i)%generation = cfs(cfs(ind)%fi)%generation+1
         call random_number(rand)
         if (rand.gt.0.5) then
           cfs(i)%gender = 1
         else
           cfs(i)%gender = 0
         end if
!Perform mutations

         do r = 1,maxs ! cas - only 6 attributes maxs
            call random_number(rand)	 
            if(rand<mute) then
               call random_number(rand)
               if(rand<=0.5) then 
                  call random_number(rand)
                  cfs(i)%sv(r) = cfs(i)%sv(r) +  cfs(i)%sv(r)*rand*effect
               else 
                  call random_number(rand)
                  cfs(i)%sv(r) = cfs(i)%sv(r) - cfs(i)%sv(r)*rand*effect
               endif 
            endif
         enddo

         if(cfs(i)%sv(1)>cfs(i)%sv(3)) cfs(i)%sv(1) = cfs(i)%sv(3)
         if(cfs(i)%sv(3)>365) cfs(i)%sv(3) = 365
         if(cfs(i)%sv(3)<1) cfs(i)%sv(3) = 1
         if(cfs(i)%sv(1)<1.) cfs(i)%sv(1) = 1

         if(cfs(i)%sv(2)>1) cfs(i)%sv(2) = 1
         if(cfs(i)%sv(2)<0.01) cfs(i)%sv(2) = 0.01
         if(cfs(i)%sv(4)<1.) cfs(i)%sv(4) = 1.
         if(cfs(i)%sv(4)<1.) cfs(i)%sv(4) = 1.
         if(cfs(i)%sv(4)>1600.) cfs(i)%sv(4) = 1600.
         if(cfs(i)%sv(5)<1.) cfs(i)%sv(5) = 1.
!CAS
         cfs(i)%type=cfs(cfs(ind)%fi)%type
         cfs(i)%layer=cfs(cfs(ind)%fi)%layer
         cfs(i)%pos=cfs(cfs(ind)%fi)%pos
         cfs(i)%starttime=cfs(cfs(ind)%fi)%starttime
         cfs(i)%endtime=cfs(cfs(ind)%fi)%endtime
         cfs(i)%wdepth=cfs(cfs(ind)%fi)%wdepth
         cfs(i)%temp=cfs(cfs(ind)%fi)%temp
         cfs(i)%salt=cfs(cfs(ind)%fi)%salt
         cfs(i)%dens=cfs(cfs(ind)%fi)%dens

         cfs(i)%num=int(cfs(i)%av(14))
         cfs(i)%age=0.0
         cfs(i)%mirr=cfs(cfs(ind)%fi)%mirr
         cfs(i)%mixl=cfs(cfs(ind)%fi)%mixl

         cfs(cfs(ind)%fi)%nsip=cfs(cfs(ind)%fi)%nsip+1
!CAS
         newpop = newpop + 1
      enddo
   enddo
8 continue
!???
  pop = newpop - 1
  nflt= newpop - 1

print*,cfs(pop-5:pop)%num
print*,cfs(pop-5:pop)%av(11)
print*, 'newpop', newpop, 'pop', pop

9 continue

end subroutine calrep
