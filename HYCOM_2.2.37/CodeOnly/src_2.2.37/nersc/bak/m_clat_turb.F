      module m_clat_turb
      contains
      subroutine clat_turb
c
c --- -------------------------------------------------------------------
c --- read numerical values of the latent transfer coefficient based
c --- on the tabelled values in isemer et al, j clim., p. 1180, 1989
c ---  n o t e :
c --- i-index gives the wind at 10m height
c ---            from 0 to 30 m/s in intervals of 2 m/s
c --- j-index gives the virtual air-sea temperature difference at 10m height
c ---            from -7 to +7 deg c in intervals of .5 deg c
c --- for all but the equatorial and sub-equatorial waters, virtual temp is
c --- close to real temp (see gill, 1982, p. 41)
c --- -------------------------------------------------------------------
c


      use mod_xc
      use mod_common_ice
      implicit none

      REAL d_wind,d_temp,clat1(16,29)
      integer skip
      integer i,j
      logical :: ex
c
c      dimension clat1(16,29)
c
c --- read data file
      inquire(exist=ex,file='iwh_tabulated.dat')
      if (.not.ex) then
         if (mnproc==1) 
     &   write(lp,*) 'clat_turb: iwh_tabulated.dat does not exist'
         call xcstop ('(clat_turb)')
         stop '(clat_turb)'
      end if
      open(19,file='iwh_tabulated.dat',STATUS= 'OLD')
      if (mnproc==1) print *,'reading Data/iwh_tabulated.dat'

      do 100 j=1,29       
      do 100 i=1,16
         read(19,*) d_wind,d_temp,clat1(i,j)
         clat1(i,j)=clat1(i,j)*1.e-3
  100 continue
      close(19)
c
c --- let  i=1  represents wind speeds in the interval [0, 2) m/s
c ---      i=2  represents wind speeds in the interval [2, 4) m/s ...
c
      do 111 j=1,29
      do 110 i=1,15
  110 clat(i,j)=(clat1(i,j)+clat1(i+1,j))*.5
  111 clat(16,j)=clat(15,j)
c
c
c --- let  j=1  represents temp differences in the interval [-7, -6.5) deg c
c ---      j=2  represents temp differences in the interval [-6.5, -6) deg c...
c
      do 121 i=1,16
      do 120 j=1,28
  120 clat(i,j)=(clat1(i,j)+clat1(i,j+1))*.5
  121 clat(i,29)=clat(i,28)
c
      end subroutine clat_turb
      end module m_clat_turb
