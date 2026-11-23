module m_calcconst
contains
subroutine calcconst(kk,U,deep,a,b,c,cflag,bflag)
   implicit none 
   integer, intent(in) :: kk
   character(len=1), intent(in) :: cflag
   character(len=1), intent(in) :: bflag
   real, intent(in), dimension(kk)::U,deep
   real, intent(out), dimension(kk)::a,b
   real, intent(out), dimension(0:kk)::c

   integer i,j,k, dimen, eq, info
   real*8, dimension(3*kk+1)::x,rhs,Z
   real*8, dimension(4*(3*kk+1))::lapwork
   integer, dimension(3*kk+1)::IPVT
   real RCOND
   real*8, dimension(3*kk+1,3*kk+1)::Amatr
     
   dimen = 3*kk+1

! Generate Amatr and rhs

   rhs = 0.
   Amatr = 0.

! First layer
   Amatr(1,1) = -1.                ! c_0
   Amatr(1,2) = 1.                 ! c_1
   Amatr(1,3) = .5*deep(1)         ! b_1
   Amatr(1,4) = (.5*deep(1))**2    ! a_1

   Amatr(2,3) = 1.                 ! b_1
   Amatr(2,4) = deep(1)              ! a_1

   Amatr(3,1) = .5                   ! c_0
   Amatr(3,2) = .5                   ! c_1
   Amatr(3,3) = 3./8.*deep(1)      ! b_1
   Amatr(3,4) = 7./24.*deep(1)**2  ! a_1
   rhs(3) = U(1)

   eq = 3
! Interior layers
! 
   do i=2,kk 
! Cont i function 

      eq = eq + 1
      Amatr(eq,3*(i-1)-1) =-1.            ! C(i-1)
      Amatr(eq,3*(i-1)  ) =-deep(i-1)       ! b(i-1)
      Amatr(eq,3*(i-1)+1) =-deep(i-1)**2    ! a(i-1)

      Amatr(eq,3*(i  )-1) = 1.            ! C(i)
      Amatr(eq,3*(i  )  ) = deep(i-1)       ! b(i)
      Amatr(eq,3*(i  )+1) = deep(i-1)**2    ! a(i)

! Cont derivative

      eq = eq + 1
      if (cflag == 'A') then
         ! Continuity in derivative
         Amatr(eq,3*(i-1)  ) =-1.            ! b(i-1)
         Amatr(eq,3*(i-1)+1) =-deep(i-1)*2.  ! a(i-1)
         Amatr(eq,3*(i  )  ) = 1.            ! b(i)
         Amatr(eq,3*(i  )+1) = deep(i-1)*2.  ! a(i)
      elseif (cflag == 'B') then
         ! interpolates average interface velocity (break continuity in derivative)
         Amatr(eq,3*(i  )-1) = 1.            ! C(i)
         Amatr(eq,3*(i  )  ) = deep(i-1)       ! b(i)
         Amatr(eq,3*(i  )+1) = deep(i-1)**2    ! a(i)
         rhs(eq)=0.5*(U(i-1) + U(i))
      elseif (cflag == 'C') then
         ! interpolates weighted average interface velocity (break continuity in derivative)
         Amatr(eq,3*(i  )-1) = 1.            ! C(i)
         Amatr(eq,3*(i  )  ) = deep(i-1)       ! b(i)
         Amatr(eq,3*(i  )+1) = deep(i-1)**2    ! a(i)
         if (i == 2) then
            rhs(eq)=0.5*(deep(i-1)*U(i-1) +  (deep(i)-deep(i-1))*U(i) )/deep(i)
         else
            rhs(eq)=0.5*( (deep(i-1)-deep(i-2))*U(i-1) + (deep(i)-deep(i-1))* U(i))  / &
                                  (deep(i)-deep(i-2))
         endif
      else
         stop 'invalid cflag in calcconst' 
      endif
 

! Mean 
      eq = eq + 1
      Amatr(eq,3*(i  )-1) = 1.                                                 ! C(i)
      Amatr(eq,3*(i  )  ) = 0.5*(deep(i-1)+deep(i))                            ! b(i)
      Amatr(eq,3*(i  )+1) = (deep(i-1)**2 + deep(i)*deep(i-1)+deep(i)**2)/3.0  ! a(i)
      rhs(eq) = U(i)
   enddo

! Closure

   eq = eq + 1 
   i = kk  
   if (bflag == 'a') then
      Amatr(eq,3*(i  )-1) = 1.          ! C(kk)
      Amatr(eq,3*(i  )  ) = deep(i)       ! b(kk)
      Amatr(eq,3*(i  )+1) = deep(i)**2    ! a(kk)
      rhs(eq)=U(i)
   elseif (bflag == 'b') then
      Amatr(eq,3*(i  )  ) = 1.            ! b(i)
      Amatr(eq,3*(i  )+1) = deep(i)*2.    ! a(i)
   else
       stop 'invalid bflag in calcconst'
   endif
          

! Solve

!   write(*,'(13g10.2)')transpose(Amatr)

   ! Arch definition time ...







   !call DGECON('1',dimen,Amatr,dimen,'I',RCOND,lapwork,ipvt,info)
   call DGESV(dimen,1,Amatr,dimen,ipvt,rhs,dimen,info)







   c(0) = rhs(1)
   c(1:kk) = rhs(1+1:3*kk+1-2:3)
   b(1:kk) = rhs(1+2:3*kk+1-1:3)
   a(1:kk) = rhs(1+3:3*kk+1-0:3)
       
         
end subroutine  calcconst
end module  m_calcconst
