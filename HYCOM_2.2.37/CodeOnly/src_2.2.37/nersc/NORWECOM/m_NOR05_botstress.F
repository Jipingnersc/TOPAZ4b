      module m_NOR05_botstress
      contains
      subroutine NOR05_botstress(n,tbs)
      !-----------------------------------------------------------------!
      ! Routine to calculate total bottom  stress, to be used in the    !
      ! sedimentation part of the NORWECOM ecosystem model.             !
      ! Last modified: 22.06.06, Cecilie Hansen                         ! 
      ! KAL -- changed to fixed-source format....
      !-----------------------------------------------------------------! 
      use mod_xc
      implicit none

      integer,intent(in) :: n

      real,intent(out),dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: tbs

      real    :: p_h,p_l,p_b,s_l,d_sl
      real    :: factor
      real    :: u_bs,v_bs
      real    :: tbs_min,tbs_max
      integer :: i,j,k

      include 'common_blocks.h'
   
      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point.
      !KAL - botstress uses hor. derivation so eff. margin is reduced
      !margin=0
      tbs_min=1000.
      tbs_max=-1000.
      d_sl=10. !depth boundary layer

!      if (mnproc .eq. 1 .and. i.eq.4 .and. j.eq.57 .and. k.eq.5) &
!     & print*,'      RB: botstress 50: v(4,58,5,1) =',v(4,58,5,1)

      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      if(ip(i,j)==1)then
         u_bs=0.
         v_bs=0.
         do k=1,kdm
            p_h=p(i,j,k)
            p_l=p(i,j,k+1)
            p_b=p(i,j,kk+1)
            s_l=p_b-d_sl*onem
            if(p_h-s_l.gt.0.)then
               factor=1.
            elseif(p_h-s_l.lt.0 .and.p_l-s_l.gt.0)then
               factor=(p_l-s_l)/(p_l-p_h) 
            else
               factor=0.
            endif
            u_bs=u_bs+(u(i,j,k,n)*iu(i,j)+u(i+1,j,k,n)*iu(i+1,j))*factor
     &          /(max(iu(i,j)+iu(i+1,j),1)) 

!      if (mnproc .eq. 1 .and. i.eq.4 .and. j.eq.57 .and. k.eq.5) &
!     & print*,'      RB: botstress 55: v(4,58,5,1) =',v(4,58,5,1)

            v_bs=v_bs+(v(i,j,k,n)*iv(i,j)+v(i,j+1,k,n)*iv(i,j+1))*factor
     &          /(max(iv(i,j)+iv(i,j+1),1)) 
         enddo

         u_bs=u_bs/d_sl+(ubavg(i,j,n)*iu(i,j)+ubavg(i+1,j,n)*iu(i+1,j))
     &                 /(max(iu(i,j)+iu(i+1,j),1))

!      if (mnproc .eq. 1 .and. i.eq.4 .and. j.eq. 57) &
!     & print*,'      RB: botstress 56: v_bs =',v_bs,' d_sl =',d_sl

         v_bs=v_bs/d_sl+(vbavg(i,j,n)*iv(i,j)+ubavg(i,j+1,n)*iv(i,j+1))
     &                 /(max(iv(i,j)+iv(i,j+1),1))

!      if (mnproc .eq. 1 .and. i.eq.4 .and. j.eq. 57) &
!     & print*,'      RB: botstress 57: tbs(4,57) =',tbs(4,57)
!      if (mnproc .eq. 1 .and. i.eq.4 .and. j.eq. 57) &
!     & print*,'      RB: botstress 57: v_bs =',v_bs

         tbs(i,j)=cb*(sqrt(u_bs**2+v_bs**2)+cbar)/d_sl              

!      if (mnproc .eq. 1 .and. i.eq.4 .and. j.eq. 57) &
!     & print*,'      RB: botstress 60: tbs(4,57) =',tbs(4,57)
      endif
      enddo
      enddo

!      if (mnproc .eq. 1) &
!     & print*,'      RB: botstress 62: tbs(4,57) =',tbs(4,57)

      end subroutine NOR05_botstress
      end module m_NOR05_botstress
