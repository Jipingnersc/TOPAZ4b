      module m_NOR05_initialize_sed
       
      contains

      subroutine NOR05_initialize_sed()
      use mod_xc
      use mod_necessary_ecovars

      implicit none
      include 'common_blocks.h'
      
      
      integer :: i,j,kbio,kmax
     
      if (mnproc==1) print*,'initializing bot_layer'
      
      do j=1-nbdy,jdm+nbdy
      do i=1-nbdy,idm+nbdy
            bot_layer(i,j,1,init)=40.0
            bot_layer(i,j,2,init)=40.0
            bot_layer(i,j,1,ipho)=4.0  ! AS20111026: 20.0
            bot_layer(i,j,2,ipho)=4.0  ! AS20111026: 20.0
            bot_layer(i,j,1,isil)=27.0 ! AS20111026: 200.0
            bot_layer(i,j,2,isil)=27.0 ! AS20111026: 200.0
            bot_layer(i,j,1,idet)=30.0
            bot_layer(i,j,2,idet)=30.0
            bot_layer(i,j,1,isis)=20.0
            bot_layer(i,j,2,isis)=20.0
            bot_layer(i,j,1,ifla)=0.0
            bot_layer(i,j,2,ifla)=0.0
            bot_layer(i,j,1,idia)=0.0
            bot_layer(i,j,2,idia)=0.0
            bot_layer(i,j,1,ioxy)=3000.0
            bot_layer(i,j,2,ioxy)=3000.0
#ifdef DETPHO
            bot_layer(i,j,1,idetp)=3.0
            bot_layer(i,j,2,idetp)=3.0
#endif
      enddo
      enddo
      end subroutine NOR05_initialize_sed
      end module m_NOR05_initialize_sed
