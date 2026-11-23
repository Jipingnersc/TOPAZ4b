
!Calanus stage dynamics model based on Carlotti and Wolf 1998, as opposed to
!the previous version of Calsta which was based on Fiksen 2000. This routine 
!calculates growth from a regular bioenergetics model with ingestion, egestion
!and respiration processes for copepod stages N3 to Adult
!-------------------------------------------------------------------------

subroutine calsta(tid,hour)

!-------------------------------------------------------------------------
   use mod_CAL06_calnusparticles  !include     'comcal.txt'
   implicit none
   integer, intent(in) :: tid,hour

   integer 	 d,e,n,i,x,up,down,f
   real         dev,growth,t,met,phyto,rand,rloss,owmet
   real         Q10,QR10,AW,fmax,lamda,resp1,resp2
   real         i1,i2,i3,ingestion,r1,r2,r3,resp,egestion,alpha
!parameters used in the new grazing formulation:
   real, parameter :: imax_n=0.53/24.0 , &
                      imax_c=0.43/24.0 , &
                      Iv_n    =0.0324, &
                      Iv_c    =0.0219, &
                      phyto0_n=0.0, &
                      phyto0_c=0.0
 
   open(unit=811,file='mortality_out',status='old', &
             form='formatted',position='append')
 119 format(a12,i12,a8,i8,a45) 

!Bioenergetics parameters
   Q10 =  2.1
   QR10 = 3.4
   AW =   0.8
   fmax = 100.
   lamda = 0.3
   resp1 = 0.01/24.  !NB Changed from per day in C&W to hourly values
   resp2 = 0.2
!   imaxim = 1.50/24.  ! skal være 0.86 AS use value from C&W 1998  
!   imaxim = 0.7/24.    !Reduced maximum ingestion to accomodate development rates from Campbell 2001
   alpha = 0.8			 
   rloss = 0.1  !0.1
   owmet = 0.001/24.


!Calculate stage dynamics of calanoid copepods
   do i = 1,pop
      if(cfs(i)%av(0)==0.0) goto 1
      d = int(cfs(i)%av(8))
      e = int(cfs(i)%av(9))
      n = int(cfs(i)%av(10))
      up = 0
      down = 0
      ingestion=0

      t = cfs(i)%temp
      phyto = cfs(i)%av(12)
      
      if (i==1)   print*, 'Food conc', phyto, cfs(i)%av(12)
      x = int(cfs(i)%av(1))

      growth = 0.
!      avtemp(x) = avtemp(x) + t
!      numtemp(x) = numtemp(x) + 1.
      if(x<3) then         !Egg and nauplia 1 and 2

         if(hour==23) then    !Update egg and nauplia <3 only on a daily basis
            dev = a(x)*(t+9.11)**(-2.05)	  !Stage specific development
            cfs(i)%av(6) = cfs(i)%av(6) + 1./dev 
            if(cfs(i)%av(6)>=1.) then 
               cfs(i)%av(1) = cfs(i)%av(1) + 1.0     !Molting
               cfs(i)%av(6) = 0.
               agemolt(x) = agemolt(x) + cfs(i)%av(4)
               molters(x) = molters(x) + 1.0
               cfs(i)%av(4) = 0.0   !Set stage longevity to zero
            endif
         endif

      elseif(x<7) then    !Bioenergetics for N3-N6 

         i1 = imax_n*Q10**(t*0.1)
         i2 = cfs(i)%av(3)**AW
!AS         i3 = 1.1*phyto/(15.+phyto)                !Fitted to Campbell et al 2001, with 90% of max feeding at 70 mg C M-3    
         i3 = (1-exp(-Iv_n*(phyto-phyto0_n))) !Fitted to Campbell et al 2001
         ingestion = max(0.,i1*i2*i3) 		!Correct for average feeding 
         egestion = lamda*ingestion
         cfs(i)%av(15) = ingestion
         cfs(i)%av(13) = egestion
         r1 = resp1*cfs(i)%av(3)**AW 
         r2 = QR10**(t*0.1)
         r3 = resp2*ingestion	  !Active metabolism 20% of ingestion
         resp = r1*r2 + r3
         growth = ingestion - egestion - resp
         cfs(i)%av(16) = growth/cfs(i)%av(3)*1.
         cfs(i)%av(3) = cfs(i)%av(3) + growth 
         if((hour==23).and.(cfs(i)%av(3)>=st_wgt(x+1))) then
            cfs(i)%av(1) = cfs(i)%av(1) + 1.0 !Molting
            agemolt(x) = agemolt(x) + cfs(i)%av(4)
            molters(x) = molters(x) + 1.0
            cfs(i)%av(4) = 0.0   !Set stage longevity to zero
         endif 

      elseif(x<11) then    !Bioenergetics for C1-C4 

         i1 = imax_c*Q10**(t*0.1)
         i2 = cfs(i)%av(3)**AW
!AS         i3 = 1.12*phyto/(30.+phyto)                !Fitted to Campbell et al 2001, with 90% of max feeding at 71 mg C M-3    
         i3 = (1-exp(-Iv_c*(phyto-phyto0_c)))  !Fitted to Campbell et al 2001
         ingestion = max(0.,i1*i2*i3) 		!Correct for average feeding 
         egestion = lamda*ingestion
         cfs(i)%av(15) = ingestion
         cfs(i)%av(13) = egestion
         r1 = resp1*cfs(i)%av(3)**AW 
         r2 = QR10**(t*0.1)
         r3 = resp2*ingestion	  !Active metabolism 20% of ingestion
         resp = r1*r2 + r3
         growth = ingestion - egestion - resp
         cfs(i)%av(16) = growth/cfs(i)%av(3)*1.
         cfs(i)%av(3) = cfs(i)%av(3) + growth
         if((hour==23).and.(cfs(i)%av(3)>=st_wgt(x+1))) then
            cfs(i)%av(1) = cfs(i)%av(1) + 1.0 !Molting
            agemolt(x) = agemolt(x) + cfs(i)%av(4)
            molters(x) = molters(x) + 1.0
            cfs(i)%av(4) = 0.0   !Set stage longevity to zero
         endif 
  
      elseif(x==11) then     !Copepodid stage V

         i1 = imax_c*Q10**(t*0.1)
         i2 = cfs(i)%av(3)**AW
!AS         i3 = 1.12*phyto/(30.+phyto)                !Fitted to Campbell et al 2001, with 90% of max feeding at 71 mg C M-3    
         i3 = (1-exp(-Iv_c*(phyto-phyto0_c)))   !Fitted to Campbell et al 2001
         ingestion = max(0.,i1*i2*i3) 		!Correct for average feeding 
         egestion = lamda*ingestion
         cfs(i)%av(15) = ingestion
         cfs(i)%av(13) = egestion
         r1 = resp1*cfs(i)%av(3)**AW  
         r2 = QR10**(t*0.1)
         r3 = resp2*ingestion
         resp = r1*r2 + r3
         growth = ingestion - egestion - resp
         cfs(i)%av(16) = growth/cfs(i)%av(3)*1.
         cfs(i)%av(3) = cfs(i)%av(3) + growth
         if((hour==23).and.(cfs(i)%av(3)>=st_wgt(x+1))) then
            if(tid>cfs(i)%sv(3)) then 
               cfs(i)%av(1) = 12.0  !Moulting to overwinter
               agemolt(x) = agemolt(x) + cfs(i)%av(4)
               molters(x) = molters(x) + 1.0
               cfs(i)%av(4) = 0.0   !Set stage longevity to zero
            else 
       	       cfs(i)%av(1) = 13.0  !Moulting to mature
               agemolt(x) = agemolt(x) + cfs(i)%av(4)
               molters(x) = molters(x) + 1.0
               call random_number(rand)
               if(cfs(i)%gender==0) then
                !Kill the males of the population, assuming a sex ratio of 0.5
                 cfs(i)%av(2) = 0.  
                 write(811,119) 'Individual #', cfs(i)%num, 'on day', tid, &!
                                ' killed - male .... 1'!
               endif
            endif
         endif

      elseif(x==12) then     !Overwintering copepodid stage V
         if(cfs(i)%av(11)==0.0) then  !In diapause
            met = cfs(i)%av(3)*owmet  !Diapause metabolism, 0.001 in Fiksen 2001
            if(met<=cfs(i)%av(5)) then
               cfs(i)%av(5) = cfs(i)%av(5)-met
            else
               cfs(i)%av(3) = cfs(i)%av(3)-met
            endif
            if((tid==int(cfs(i)%sv(1))).and.(hour==23)) then !Wake up day from strategy vector
               up = 1                    !Wake up from diapause in spring, go from 0 to 3 
               cfs(i)%av(1) = 13.0              !Mature	
               agemolt(x) = agemolt(x) + cfs(i)%av(4)
               molters(x) = molters(x) + 1.0
               cfs(i)%av(4) = 0.0   !Set stage longevity to zero
            endif
         elseif(cfs(i)%av(11)==1.) then !Active
            i1 = imax_c*Q10**(t*0.1)
            i2 = cfs(i)%av(3)**AW
!AS            i3 = 1.12*phyto/(30.+phyto)                !Fitted to Campbell et al 2001, with 90% of max feeding at 71 mg C M-3    
            i3 = (1-exp(-Iv_c*(phyto-phyto0_c)))  !Fitted to Campbell et al 2001
            ingestion = max(0.,i1*i2*i3) 		!Correct for average feeding 
            egestion = lamda*ingestion
            cfs(i)%av(15) = ingestion
            cfs(i)%av(13) = egestion
            r1 = resp1*cfs(i)%av(3)**AW 
            r2 = QR10**(t*0.1)
            r3 = resp2*ingestion
            resp = r1*r2 + r3
            growth = ingestion - egestion - resp
            cfs(i)%av(16) = growth/cfs(i)%av(3)*1.
            if(growth>=0.) then
               cfs(i)%av(5) = cfs(i)%av(5) + growth
            elseif(growth<0.) then 
               if(cfs(i)%av(5)>abs(growth)) then  
                  cfs(i)%av(5) = cfs(i)%av(5) + growth
               else 
                  cfs(i)%av(3) = cfs(i)%av(3) + growth
               endif
            endif
            if((cfs(i)%av(5)/cfs(i)%av(3))>=cfs(i)%sv(2)) then
               down = 1  !Migrating down to start diapause at fsr
            endif
         elseif(cfs(i)%av(11)==2.) then    !Ascending
            met = cfs(i)%av(3)*owmet        !Diapause metabolism, 0.001 in Fiksen 2001
            cfs(i)%av(5) = cfs(i)%av(5)-met
         endif
 
      elseif(x==13) then           !C5->Adult females
         if(cfs(i)%av(11)==3.) then 
            met = cfs(i)%av(3)*owmet        !Diapause metabolism, 0.001 in Fiksen 2001
            if(met<=cfs(i)%av(5)) then
               cfs(i)%av(5) = cfs(i)%av(5)-met
            else
               cfs(i)%av(3) = cfs(i)%av(3)-met
            endif
         elseif(cfs(i)%av(11)==1.0) then
            i1 = imax_c*Q10**(t*0.1)
            i2 = cfs(i)%av(3)**AW
!AS            i3 = 1.12*phyto/(30.+phyto)                !Fitted to Campbell et al 2001, with 90% of max feeding at 71 mg C M-3    
            i3 = (1-exp(-Iv_c*(phyto-phyto0_c)))  !Fitted to Campbell et al 2001
            ingestion = max(0.,i1*i2*i3) 		      !Correct for average feeding 
            egestion = lamda*ingestion

            cfs(i)%av(15) = ingestion
            cfs(i)%av(13) = egestion
            r1 = resp1*cfs(i)%av(3)**AW 
            r2 = QR10**(t*0.1)
            r3 = resp2*ingestion
            resp = r1*r2 + r3
            growth = ingestion - egestion - resp
            cfs(i)%av(16) = growth/cfs(i)%av(3)*1.
            if(cfs(i)%av(3)<=cmat) then 
               if(growth>=0.) then
                  cfs(i)%av(3) = cfs(i)%av(3) + growth
               elseif(growth<0.) then 
                  if(cfs(i)%av(5)>abs(growth)) then  
                     cfs(i)%av(5) = cfs(i)%av(5) + growth
       	          else 
 	             cfs(i)%av(3) = cfs(i)%av(3) + growth
                  endif
               endif
            else						  !Allocation to fat
               if(growth>=0.) then
                  cfs(i)%av(5) = cfs(i)%av(5) + growth
               elseif(growth<0.) then 
                  if(cfs(i)%av(5)>abs(growth)) then  
	             cfs(i)%av(5) = cfs(i)%av(5) + growth
	          else 
	             cfs(i)%av(3) = cfs(i)%av(3) + growth
	          endif
               endif
            endif
         endif  !av11 booolean
      endif  !Stage boolean	
      if(cfs(i)%av(5)<0.) cfs(i)%av(5) = 0.  !Do not allow negative fat reserves
 !Update diapause variables and induce molting cost
      if(up==1.) cfs(i)%av(11) = 3.
      if(down==1.) cfs(i)%av(11) = 2.
      if(real(x)<cfs(i)%av(1)) cfs(i)%av(3) = cfs(i)%av(3)*(1. - rloss)
1  enddo
   close(811)

end subroutine calsta
