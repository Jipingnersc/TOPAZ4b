      module m_bigrid
      contains

      subroutine bigrid(depth)
      use mod_common
      use m_indxi 
      use m_indxj 
c
c --- set loop bounds for irregular basin in c-grid configuration
c --- q,u,v,p are vorticity, u-velocity, v-velocity, and mass points, resp.
c --- 'depth' = basin depth array, zero values indicate land
c
      implicit none
c
      real depth(idm,jdm)
c
      integer i,j,nfill,nzero,isec,ifrst,ilast
      character char3*3
c
      integer nchar
      character fmt*13
      data nchar/120/
      data fmt/'(i4,1x,120i1)'/
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 17 j=1,jj
      do 17 i=1,ii
      ip(i,j)=0
      iq(i,j)=0
      iu(i,j)=0
 17   iv(i,j)=0
c
c --- detect (and abort on) single-width inlets and 1-point seas.
 16   nfill=0
      do 15 j=1,jj1
      do 15 i=1,ii1
      nzero=0
      if (depth(i,j).gt.0.) then
        if (i.eq.1  .or.depth(i-1,j).le.0.) nzero=nzero+1
        if (i.eq.ii1.or.depth(i+1,j).le.0.) nzero=nzero+1
        if (j.eq.1  .or.depth(i,j-1).le.0.) nzero=nzero+1
        if (j.eq.jj1.or.depth(i,j+1).le.0.) nzero=nzero+1
        if (nzero.ge.3) then
          write (lp,'(a,i4,a,i4,a,i1,a)')
     +      'warning - dh(',i,',',j,') has ',
     +      nzero,' land neighbours'
          depth(i,j)=0.
          nfill=nfill+1
        end if
      end if
 15   continue
      if (nfill.gt.0) then
         goto 16
!        write(lp,'(/a/)') 
!     .   'Must correct bathymetry before running HYCOM'
!        call flush(lp)
!        stop '(bigrid)'
      endif
c
c --- mass points are defined where water depth is greater than zero
c --- u,v points are located halfway between any 2 adjoining mass points
c --- 'interior' q points require water on all 4 sides.
c --- 'promontory' q points require water on 3 (or at least 2
c --- diametrically opposed) sides
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj1
        do i=1,ii1
          if (depth(i,j).gt.0.) then
            ip(i,j)=1
          endif
        enddo
        do i=2,ii1
          if (ip(i-1,j).gt.0.and.ip(i,j).gt.0) then
            iu(i,j)=1
          endif
        enddo
      enddo
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=2,jj1
        do i=1,ii1
          if (ip(i,j-1).gt.0.and.ip(i,j).gt.0) then
            iv(i,j)=1
          endif
        enddo
        do i=2,ii1
          if (min(ip(i,j),ip(i-1,j),ip(i,j-1),ip(i-1,j-1)).gt.0) then
            iq(i,j)=1
          elseif ((ip(i  ,j).gt.0.and.ip(i-1,j-1).gt.0).or.
     .            (ip(i-1,j).gt.0.and.ip(i  ,j-1).gt.0)    ) then
            iq(i,j)=1
          endif
        enddo
      enddo
c
c --- determine loop bounds for vorticity points, including interior and
c --- promontory points
      call indxi(iq,ifq,ilq,isq)
      call indxj(iq,jfq,jlq,jsq)
c
c --- vorticity points on borders
c --- 'left' border (j increasing to right)
ccc      do 6 j=1,jj1
ccc      do 6 i=2,ii1
ccc      if (iq(i,j).ne.0) go to 6
ccc      if (ip(i,j).gt.0.and.ip(i-1,j).gt.0) iq(i,j)=2
ccc  6   continue
c --- 'right' border
ccc      do 7 j=2,jj
ccc      do 7 i=2,ii1
ccc      if (iq(i,j).ne.0) go to 7
ccc      if (ip(i,j-1).gt.0.and.ip(i-1,j-1).gt.0) iq(i,j)=3
ccc  7   continue
c --- 'top' border (i increasing downward)
ccc      do 8 j=2,jj1
ccc      do 8 i=1,ii1
ccc      if (iq(i,j).ne.0) go to 8
ccc      if (ip(i,j).gt.0.and.ip(i,j-1).gt.0) iq(i,j)=4
ccc  8   continue
c --- 'bottom' border
ccc      do 9 j=2,jj1
ccc      do 9 i=2,ii
ccc      if (iq(i,j).ne.0) go to 9
ccc      if (ip(i-1,j).gt.0.and.ip(i-1,j-1).gt.0) iq(i,j)=5
ccc  9   continue
c
c --- determine loop indices for mass and velocity points
      call indxi(ip,ifp,ilp,isp)
      call indxj(ip,jfp,jlp,jsp)
      call indxi(iu,ifu,ilu,isu)
      call indxj(iu,jfu,jlu,jsu)
      call indxi(iv,ifv,ilv,isv)
      call indxj(iv,jfv,jlv,jsv)
c
c --- velocity points on borders
c --- 'left' and 'right' border points for v (y pointing right)
ccc      do 13 i=1,ii1
ccc      do 13 j=1,jj1
ccc      if (ip(i,j).ne.1) go to 13
ccc      if (iv(i,j).eq.0) iv(i,j)=2
ccc      if (iv(i,j+1).eq.0) iv(i,j+1)=2
ccc 13   continue
c --- 'top' and 'bottom' border points for u (x pointing down)
ccc      do 14 i=1,ii1
ccc      do 14 j=1,jj1
ccc      if (ip(i,j).ne.1) go to 14
ccc      if (iu(i,j).eq.0) iu(i,j)=2
ccc      if (iu(i+1,j).eq.0) iu(i+1,j)=2
ccc 14   continue
c
c --- write out  -ip-  array
c --- data are written in strips nchar points wide
      return
      isec=(ii-1)/nchar
      do 9 ifrst=0,nchar*isec,nchar
      ilast=min(ii,ifrst+nchar)
      write (char3,'(i3)') ilast-ifrst
      fmt(8:10)=char3
      write (lp,'(''ip array, cols'',i5,'' --'',i5)') ifrst+1,ilast
 9    write (lp,fmt) (j,(10*ip(i,j),i=ifrst+1,ilast),j=jj,1,-1)
c
      return
      end subroutine bigrid
      end module m_bigrid
