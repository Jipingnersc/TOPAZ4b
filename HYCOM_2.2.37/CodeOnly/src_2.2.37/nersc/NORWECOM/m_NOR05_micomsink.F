      module m_NOR05_micomsink

      contains
!-----------------------------------------------------
!  subroutine that, if given vertical distribution of:
!    -variable
!    -vertical  velocity
!    -thickness of layer
!    calculates the sinking
!  from Karen A. 26.05.06
!-----------------------------------------------------

      subroutine NOR05_micomsink(n,timestep)
      use mod_xc
      use mod_necessary_ecovars

      implicit none

      include 'common_blocks.h'
      include 'sedcom.h'
      include 'biocom.h'

      integer,intent(in) :: n 
      real,intent(in) :: timestep 
      
      real    :: layer_thickn
      real    :: srdia
      integer :: i,j,k,kdonor
      
      
      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point 
      !margin=0
!--- CH: sinking velocities
!--- surface layer      
      k=1 
      do i=1-margin,ii+margin
      do j=1-margin,jj+margin
         layer_thickn=dp(i,j,k,n)/onem
         if(depths(i,j).gt.0. .and. layer_thickn .gt.0.) then
            bio(i,j,k,n,idet)=(bio(i,j,k,n,idet)*layer_thickn)/
     &                     (layer_thickn+srdet*dayinv*timestep)
            bio(i,j,k,n,isis)=(bio(i,j,k,n,isis)*layer_thickn)/
     &                     (layer_thickn+srsis*dayinv*timestep)
            if(bio(i,j,k,n,ifla).gt.0.1)then
               bio(i,j,k,n,ifla)=(bio(i,j,k,n,ifla)*layer_thickn)/
     &                        (layer_thickn+srfla*dayinv*timestep)
            endif
!use the variable sining rate
            if(bio(i,j,k,n,isil)/csil.lt.sib)then
              srdia = srdiamax
            else
!              srdia = srdiamin+(srdiamax-srdiamin)/(bio(i,j,k,n,isil)/csil)
              srdia = (srdiamin+a3l/(bio(i,j,k,n,isil)/csil))
            endif
            if(bio(i,j,k,n,idia).gt.0.1)then
               bio(i,j,k,n,idia)=(bio(i,j,k,n,idia)*layer_thickn)/
     &                        (layer_thickn+srdia*dayinv*timestep)
            endif
#if defined (YELLOWSED)
            bio(i,j,k,n,ised)=(bio(i,j,k,n,ised)*layer_thickn)/
     &                     (layer_thickn+srsed*dayinv*timestep)
            bio(i,j,k,n,iyel)=(bio(i,j,k,n,iyel)*layer_thickn)/
     &                     (layer_thickn+sryel*dayinv*timestep)
#endif /*YELLOWSED*/
#if defined DETPHO
            bio(i,j,k,n,idetp)=(bio(i,j,k,n,idetp)*layer_thickn)/
     &                     (layer_thickn+srdetp*dayinv*timestep)
#endif
         endif
      enddo 
      enddo
c      
!--- layer k=2..kdm      
c
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
        kdonor=1
        do k=2,kdm
          layer_thickn=dp(i,j,k,n)/onem
          if(depths(i,j).gt.1.) then
          if(layer_thickn .gt.0.) then
             bio(i,j,k,n,idet)=(bio(i,j,k,n,idet)*layer_thickn
     &                    +bio(i,j,kdonor,n,idet)*srdet*dayinv*timestep)
     &                    /(layer_thickn+srdet*dayinv*timestep)        
             bio(i,j,k,n,isis)=(bio(i,j,k,n,isis)*layer_thickn
     &                    +bio(i,j,kdonor,n,isis)*srsis*dayinv*timestep)
     &                    /(layer_thickn+srsis*dayinv*timestep)        
             if(bio(i,j,kdonor,n,ifla).gt.0.1)then
                bio(i,j,k,n,ifla)=(bio(i,j,k,n,ifla)*layer_thickn
     &                    +bio(i,j,kdonor,n,ifla)*srfla*dayinv*timestep)
     &                    /(layer_thickn+srfla*dayinv*timestep)        
             endif
!use the variable sining rate
            if(bio(i,j,k,n,isil)/csil.lt.sib)then
              srdia = srdiamax
            else
!              srdia = srdiamin+(srdiamax-srdiamin)/(bio(i,j,k,n,isil)/csil)
              srdia = (srdiamin+a3l/(bio(i,j,k,n,isil)/csil))
             endif
             if(bio(i,j,kdonor,n,idia).gt.0.1)then
                bio(i,j,k,n,idia)=(bio(i,j,k,n,idia)*layer_thickn
     &                         +bio(i,j,kdonor,n,idia)*srdia*dayinv
     &                         *timestep)
     &                         /(layer_thickn+srdia*dayinv*timestep)        
             endif
#if defined (YELLOWSED)
             bio(i,j,k,n,ised)=(bio(i,j,k,n,ised)*layer_thickn
     &                    +bio(i,j,kdonor,n,ised)*srsed*dayinv*timestep)
     &                    /(layer_thickn+srsed*dayinv*timestep)        
             bio(i,j,k,n,iyel)=(bio(i,j,k,n,iyel)*layer_thickn
     &                    +bio(i,j,kdonor,n,iyel)*sryel*dayinv*timestep)
     &                    /(layer_thickn+sryel*dayinv*timestep) 
#endif /*YELLOWSED*/
#if defined DETPHO       
             bio(i,j,k,n,idetp)=(bio(i,j,k,n,idetp)*layer_thickn
     &                 +bio(i,j,kdonor,n,idetp)*srdetp*dayinv*timestep)
     &                 /(layer_thickn+srdetp*dayinv*timestep)        
#endif     
             kdonor=k
          else
             bio(i,j,k,n,idet)=bio(i,j,kdonor,n,idet)
             bio(i,j,k,n,isis)=bio(i,j,kdonor,n,isis)
             bio(i,j,k,n,ifla)=bio(i,j,kdonor,n,ifla)
             bio(i,j,k,n,idia)=bio(i,j,kdonor,n,idia)
#if defined (YELLOWSED)
             bio(i,j,k,n,ised)=bio(i,j,kdonor,n,ised)
             bio(i,j,k,n,iyel)=bio(i,j,kdonor,n,iyel)
#endif /*YELLOWSED*/
#if defined DETPHO
             bio(i,j,k,n,idetp)=bio(i,j,kdonor,n,idetp)
#endif
          endif
          endif
        enddo
      enddo
      enddo
c --- in case of small concentrations, no sinking of diatoms and
c --- flagellates      
c      
      end subroutine NOR05_micomsink     
      end module m_NOR05_micomsink
