      module m_poflat
      contains
      function poflat(theta,xlat)
c
      real theta,xlat
c
c --- this routine returns pressure as function of density (sigma) and latitude
c
      integer        lp
      common/linepr/ lp
      save  /linepr/
c
      integer ix,kz
      real    p1,p2,pinthi,pintlo,thet,x,xla,z
c
      integer    kdpth,klat
      parameter (kdpth=15,klat=20)
      real pdat(kdpth,klat)
c
      real onem,thet1,thet2,dthet,xlat1,xlat2,dlat
c
      data onem/9806./  ! SI units
      data thet1,thet2,dthet/21.0,28.0,0.5/
      data xlat1,xlat2,dlat/-30.,65.,5./
c
c---  depth (m) of isopycnals of potential density 21.0, 21.5, ... , 28.0
c---  at latitudes  30s ... 65n  (source: levitus atlas)
c
      data pdat/
     . 0.,0.,0.,0., 0., 0., 0., 0., 0., 33., 87.,260.,587.,1219.,5000.,  ! 30s
     . 0.,0.,0.,0., 0., 0., 0., 0., 0., 53.,133.,273.,540.,1125.,5000.,  ! 25s
     . 0.,0.,0.,0., 0., 0., 0., 0., 0., 67.,147.,267.,480.,1063.,5000.,  ! 20s
     . 0.,0.,0.,0., 0., 0., 0., 0.,27., 73.,140.,213.,420.,1063.,5000.,  ! 15s
     . 0.,0.,0.,0., 0., 0., 0.,35.,60., 80.,127.,180.,360.,1031.,5000.,  ! 10s
     . 0.,0.,0.,0., 0., 0.,24.,47.,60., 80.,107.,160.,381.,1000.,5000.,  !  5s
     . 0.,0.,0.,0., 0.,20.,35.,52.,60., 80.,107.,160.,381.,1000.,5000.,  !  0
     . 0.,0.,0.,0., 0.,30.,40.,59.,60., 80.,107.,160.,381.,1000.,5000.,  !  5n
     . 0.,0.,0.,0., 0.,20.,35.,59.,73., 92.,108.,144.,344., 980.,5000.,  ! 10n
     . 0.,0.,0.,0., 0., 0.,41.,64.,83.,111.,141.,200.,387., 924.,5000.,  ! 15n
     . 0.,0.,0.,0., 0., 0., 0.,47.,87.,123.,163.,235.,472., 920.,5000.,  ! 20n
     . 0.,0.,0.,0., 0., 0., 0., 0.,47., 87.,140.,227.,520., 906.,5000.,  ! 25n
     . 0.,0.,0.,0., 0., 0., 0., 0., 0., 40., 87.,200.,600., 926.,5000.,  ! 30n
     . 0.,0.,0.,0., 0., 0., 0., 0., 0.,  0., 60.,147.,552., 912.,5000.,  ! 35n
     . 0.,0.,0.,0., 0., 0., 0., 0., 0.,  0., 33., 83.,333., 780.,5000.,  ! 40n
     . 0.,0.,0.,0., 0., 0., 0., 0., 0.,  0.,  0., 45.,147., 746.,5000.,  ! 45n
     . 0.,0.,0.,0., 0., 0., 0., 0., 0.,  0.,  0., 27.,104., 493.,5000.,  ! 50n
     . 0.,0.,0.,0., 0., 0., 0., 0., 0.,  0.,  0.,  0., 33., 287.,5000.,  ! 55n
     . 0.,0.,0.,0., 0., 0., 0., 0., 0.,  0.,  0.,  0.,  0., 163.,5000.,  ! 60n
     . 0.,0.,0.,0., 0., 0., 0., 0., 0.,  0.,  0.,  0.,  0., 107.,1106./  ! 65n
c
c---  quasi-hermite interpolation function (0 < xx < 1)
c
      real parabl,xx,a,b,c
ccc   real qsihmt,d
      parabl(xx,a,b,c)=b+.5*xx*(c-a+xx*(a+c-b-b))
ccc   qsihmt(xx,a,b,c,d)=(1.-xx)*parabl(xx,a,b,c)+xx*parabl(1.-xx,d,c,b)
c
      xla=(xlat*57.29578-xlat1)/dlat+1.
      ix=max(2,min(klat-2,int(xla)))
      x=max(0.,min(1.,xla-float(ix)))
      thet=(theta-thet1)/dthet+1.
      kz=max(1,min(kdpth-1,int(thet)))
      z=max(0.,min(1.,thet-float(kz)))
c
c --- horizontal interpolation: quasi-hermite. vertical interpol.: linear
c
      p1=parabl(   x,pdat(kz  ,ix-1),pdat(kz  ,ix  ),pdat(kz  ,ix+1))
      p2=parabl(1.-x,pdat(kz  ,ix+2),pdat(kz  ,ix+1),pdat(kz  ,ix  ))
      pintlo=p1*(1.-x)+p2*x
      p1=parabl(   x,pdat(kz+1,ix-1),pdat(kz+1,ix  ),pdat(kz+1,ix+1))
      p2=parabl(1.-x,pdat(kz+1,ix+2),pdat(kz+1,ix+1),pdat(kz+1,ix  ))
      pinthi=p1*(1.-x)+p2*x
      poflat=(pintlo*(1.-z)+pinthi*z)*onem
ccc   write (lp,'('' poflat'',2f7.2,2i5,2f7.2,f7.1)')
ccc  .  theta,xlat*57.29578,ix,kz,x,z,poflat/onem
      end function poflat
      end module m_poflat
c
c> Revision history
c>
c> May  2000 - conversion to SI units
