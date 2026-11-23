      module mod_checknan
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

c****************************************************************

#if defined (NOR05)
      subroutine checknan_pos(bioval,i,j,k,kbio,k_list,str,n)
      use mod_necessary_ecovars
      use mod_xc
      real,intent(in) :: bioval
      character(len=*), intent(in) :: str
      integer,intent(in) :: i,j,k,kbio,k_list,n
      include 'common_blocks.h'

      if(bioval .ne. bioval) then
        print*,'bio is NaNQ:',bioval
        print*,str
        print*,'and posision is:',i,j,k,kbio
        print*,'klist is:',k_list
        print*,'depth:',depths(i,j)
        print*,'nit,pho,sil,det:',bio(i,j,k,n,init),bio(i,j,k,n,ipho),
     &                         bio(i,j,k,n,isil),bio(i,j,k,n,idet)   
        print*,'sis,fla,dia,oxy:',bio(i,j,k,n,isis),bio(i,j,k,n,ifla),
     &                         bio(i,j,k,n,idia),bio(i,j,k,n,ioxy)
        print*,'dp is:',dp(i,j,k,n)
        stop 
      endif   
      end subroutine checknan_pos                            
#endif

c***********************************************************
                         
      subroutine checkneg(bioval,i,j,k,n,kbio,str)
      use mod_xc , only : kdm
      implicit none
      real,intent(in) :: bioval
      integer, intent(in) :: i,j,k,n,kbio
      character(len=*), intent(in) :: str
      integer :: icount

      if(bioval .lt. 0.) then
        print*,'bio is negative:',bioval
        print*,'position and component:',i,j,k,kbio
c      do icount=1,kdm
c         call print_bio(i,j,icount,n,'')
c      end do
        print*,str
        stop 
      endif  
      end subroutine checkneg

c********************************************************

      subroutine checkneg_pos(bioval,i,j,k,kbio,str)
      real,intent(in) :: bioval
      character(len=*), intent(in) :: str
      integer,intent(in) :: i,j,k,kbio
      if(bioval .lt. 0.) then
        print*,'bio is negative:',bioval
        print*,str
        print*,'and position is:',i,j,k,kbio
        stop 
      endif  
      end subroutine checkneg_pos

c********************************************************

#if defined (SCH02)
      subroutine print_bio(i,j,k,n,str)
      use mod_necessary_ecovars
      use mod_xc
      implicit none
      integer,intent(in) :: i,j,k,n
      integer :: icount
      character(len=*), intent(in) :: str
      character*25 :: biostr(15)
   
      biostr(1)='DIN'
      biostr(2)='DIC'
      biostr(3)='ALK'
      biostr(4)='PHYN'
      biostr(5)='PHYC'
      biostr(6)='CHLA'
      biostr(7)='DETN'
      biostr(8)='DETC'
      biostr(9)='HETN'
      biostr(10)='HETC'
      biostr(11)='EON'
      biostr(12)='EOC'
      biostr(13)='PCHO'
      biostr(14)='TEPC'
      biostr(15)='PRIMARY PRODUCTIVITY'
   
      print *,str
      do icount= 1,nbio
        print *,biostr(icount),bio(i,j,k,n,icount)
      end do
      end subroutine print_bio
#endif

c************************************************

#if defined (EVA85)
      subroutine print_bio(i,j,k,n,str)
      use mod_necessary_ecovars
      use mod_xc
      include 'common_blocks.h'
      implicit none
      integer,intent(in) :: i,j,k,n 
      character(len=*), intent(in) :: str
      print *,str,'NIT, PHY, ZOO',bio(i,j,k,n,1:3), & 
           'level', k,p(i,j,k)/onem
      end subroutine print_bio
#endif

c************************************************

      subroutine checkzero(bioval,str)
      real,intent(in) :: bioval
      character(len=*), intent(in) :: str

      if(bioval .eq. 0.) then
        print*,'bioval is zero:',bioval
        print*,str
      endif
      end subroutine checkzero

c***************************************************

      subroutine checkzero_pos(bioval,str,i,j,k)
      real,intent(in) :: bioval
      character(len=*), intent(in) :: str

      if(bioval .eq. 0.) then
        print*,'bioval is zero,pos:',bioval,i,j,k
        print*,str
      endif
      end subroutine checkzero_pos

c***************************************************

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

c****************************************************

      subroutine checkall(bioval,str,ldepth_box,i,j,k,n,m,k_max)
      use mod_necessary_ecovars
      use mod_xc
      implicit none
   
      real, intent(in)            :: bioval
      integer, intent(in)         :: i,j,k,n,m,k_max
      real,intent(in),dimension(i-1:i+1,j-1:j+1,1:kdm+1) :: ldepth_box
                   
      character(len=*),intent(in) :: str
   
      integer kbio,ib,jb,kb 

      if (abs(bioval) .gt. 1.0e6 .or. bioval .lt. 0.) then
        open(39,file='negative_bio.dat',
     &       status='unknown',form='formatted') 
        do kbio=1,nbio
        do ib=i-1,i+1
        do jb=j-1,j+1
        do kb=1,k_max
        write(39,*) bio(ib,jb,kb,n,kbio),bio(ib,jb,kb,m,kbio),
     &             ldepth_box(ib,jb,kb+1),ib,jb,kb,kbio 
        enddo 
        enddo 
        enddo 
        enddo 
        close(39)
        print *, str
        stop 
      endif   
      end subroutine checkall

c*******************************************************************

      end module mod_checknan



