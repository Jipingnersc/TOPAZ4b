module mod_common
   use mod_xc
   !use mod_dimensions



   integer, save :: nx,ny,nz

   integer, save :: kdm
   !integer, save :: idm,jdm,kdm
   integer, save :: ii,jj,ii1,jj1,kk
   integer, parameter :: ms=30
   !integer, save :: lp

   integer, save, allocatable, dimension(:  ) :: &
      isp,jsp,isu,jsu,isv,jsv,isq,jsq

   integer, save, allocatable, dimension(:,:) :: &
      ip,iu,iv,iq

   integer, save, allocatable, dimension(:,:) :: &
      ifq,ilq,ifp,ilp,ifu,ilu,ifv,ilv
   integer, save, allocatable, dimension(:,:) :: &
      jfq,jlq,jfp,jlp,jfu,jlu,jfv,jlv
   real, save, allocatable, dimension(:,:)   ::  &
      depths,modlon,modlat,mqlon,mqlat

   real, parameter :: onem=9806.

end module mod_common
