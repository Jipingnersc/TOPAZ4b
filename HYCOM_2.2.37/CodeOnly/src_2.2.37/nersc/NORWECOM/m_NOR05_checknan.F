      module m_NOR05_checknan
      contains 
      subroutine checknan(bioval,str)
         real,intent(in) :: bioval
         character(len=*), intent(in) :: str

         if(bioval .ne. bioval) then
            print*,'bio is NaNQ:',bioval
            print*,str
            stop 
         endif   
      end subroutine checknan                            
      subroutine checknan_pos(bioval,i,j,k,kbio,str,n)
         use mod_necessary_ecovars
         use mod_xc
         real,intent(in) :: bioval
         character(len=*), intent(in) :: str
         integer,intent(in) :: i,j,k,kbio,n
         include 'common_blocks.h'

         if(bioval .ne. bioval) then
            print*,'bio is NaNQ:',bioval
            print*,str
            print*,'and posision is:',i,j,k,kbio
            print*,'depth:',depths(i,j)
            print*,'nit,pho,sil,det:',bio(i,j,k,n,init),
     &             bio(i,j,k,n,ipho),bio(i,j,k,n,isil),bio(i,j,k,n,idet)   
            print*,'sis,fla,dia,oxy:',bio(i,j,k,n,isis),
     &             bio(i,j,k,n,ifla),bio(i,j,k,n,idia),bio(i,j,k,n,ioxy)   
#if defined (YELLOWSED)
            print*,'sed,yel:',bio(i,j,k,n,ised),bio(i,j,k,n,iyel)
#endif /*YELLOWSED*/
            stop 
         endif   
      end subroutine checknan_pos                            

      subroutine checkneg(bioval,str)
         real,intent(in) :: bioval
         character(len=*), intent(in) :: str
         if(bioval .lt. 0.) then
            print*,'bio is negative:',bioval
            print*,str
            stop 
         endif  
      end subroutine checkneg
      subroutine checkneg_pos(bioval,i,j,k,kbio,k_list,str)
         real,intent(in) :: bioval
         character(len=*), intent(in) :: str
         integer,intent(in) :: i,j,k,kbio,k_list
         if(bioval .lt. 0.) then
            print*,'bio is negative:',bioval
            print*,str
            print*,'and position is:',i,j,k,kbio
            print*,'number of  layers:',k_list
            stop 
         endif  
      end subroutine checkneg_pos
      !To be able to compile EVA85, have to remove this..
      !subroutine print_bio(i,j,k,n,str)
      !   use mod_necessary_ecovars
      !   use mod_xc
      !   integer,intent(in) :: i,j,k,n 
      !   character(len=*), intent(in) :: str
      !   print *,str
      !   print *,'DIN, DIC, ALK',bio(i,j,k,n,1:3)
      !   print *,'PHYN,PHYC, CHL',bio(i,j,k,n,4:6)
      !   print *,'DETN,DETC',bio(i,j,k,n,7:8)
      !   print *,'HETN,HETC',bio(i,j,k,n,9:10)
      !   print *,'EON,EOC',bio(i,j,k,n,11:12)
      !end subroutine print_bio
      subroutine checkzero(bioval,str)
         real,intent(in) :: bioval
         character(len=*), intent(in) :: str

         if(bioval .eq. 0.) then
            print*,'bioval is zero:',bioval
            print*,str
         endif
      end subroutine checkzero

      subroutine checkzero_pos(bioval,str,i,j,k)
         real,intent(in) :: bioval
         character(len=*), intent(in) :: str

         if(bioval .eq. 0.) then
            print*,'bioval is zero,pos:',bioval,i,j,k
            print*,str
         endif
      end subroutine checkzero_pos

      subroutine checkinf(bioval,i,j,k,kbio, str)
         real,intent(in) :: bioval
         character(len=*), intent(in) :: str
       
         if (abs(bioval) .gt. 1.0e6) then
            print *, 'bioval is inf:', bioval
            print *, 'pos. and kbio is:',i,j,k,kbio
            print *,str
            stop
         endif   
      end subroutine checkinf

      subroutine checkall(bioval,str,ldepth_box,i,j,k,n,m,k_max)
         use mod_necessary_ecovars
         use mod_xc
         implicit none
         
         real, intent(in)            :: bioval
         integer, intent(in)         :: i,j,k,n,m,k_max
         real,intent(in),dimension(i-1:i+1,j-1:j+1,1:kdm+1) :: 
     &      ldepth_box
                         
         character(len=*),intent(in) :: str
         
         integer kbio,ib,jb,kb 

         if (abs(bioval) .gt. 1.0e6 .or. bioval .lt. 0.) then
            open(39,file='negative_bio.dat',status='unknown',
     &           form='formatted') 
            do kbio=1,nbio
            do ib=i-1,i+1
            do jb=j-1,j+1
            do kb=1,k_max
             write(39,*) bio(ib,jb,kb,n,kbio),bio(ib,jb,kb,m,kbio),
     &       ldepth_box(ib,jb,kb+1),ib,jb,kb,kbio 
            enddo 
            enddo 
            enddo 
            enddo 
            close(39)
            print *, str
            stop 
         endif   
      end subroutine checkall
      end module m_NOR05_checknan



