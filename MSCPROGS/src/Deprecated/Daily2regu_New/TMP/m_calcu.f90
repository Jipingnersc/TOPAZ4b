

module m_calcu
contains
   function calcu(k,a,b,c,d,z,ic,undef)
   implicit none
   integer, intent(in)::k
   real, intent(in), dimension(k)::a,b,d
   real, intent(in), dimension(0:k)::c
   real, intent(in)::z
   real, intent(in)::undef
   real calcu
   integer, intent(out):: ic
   integer i

   if (z .le. d(1)) then
      if (z.le.0.5*d(1)) then
         calcu = c(0)
         ic = 1
         return
      else
         calcu = a(1)*z**2+b(1)*z+c(1)
         ic = 1
         return
      endif
   endif

   do i = 2,k
      if (z .le. d(i)) then
         calcu = a(i)*z**2+b(i)*z+c(i)
         ic = i
         return
      endif
   enddo
   calcu=undef
!   stop 'missed in calcu'
end function calcu
end module m_calcu
