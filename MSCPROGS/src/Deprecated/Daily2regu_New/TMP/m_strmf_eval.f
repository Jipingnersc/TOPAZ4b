

      module m_strmf_eval
      contains

      subroutine strmf_eval(idm,jdm,strmf,ub,vb,depth,qlat,qlon)
      use m_spherdist
      implicit none
c
      integer idm,jdm
      real strmf(idm,jdm)
     .    ,ub(idm,jdm),vb(idm,jdm)
     .    ,depth(idm,jdm)
     .    ,qlat(0:idm+1,0:jdm+1),qlon(0:idm+1,0:jdm+1)
c
      integer i,j,im1,jm1,ip1,jp1
      real depthu(idm,jdm),depthv(idm,jdm) 
     .    ,scvx(idm,jdm),scuy(idm,jdm)
c
c --- ------------------------------------------------------------------
c --- calculate scale factors and depths in u- and v-points
c --- ------------------------------------------------------------------
c
      do j=1,jdm
        jp1=j+1
        jm1=mod(j+jdm-2,jdm)+1
        do i=1,idm
          im1=max(1,i-1)
          ip1=i+1
          scvx(i,j)=spherdist(qlon(ip1,j),qlat(ip1,j),
     .                        qlon(i  ,j),qlat(i  ,j))
          scuy(i,j)=spherdist(qlon(i,jp1),qlat(i,jp1),
     .                        qlon(i,j  ),qlat(i,j  ))
          depthu(i,j)=min(depth(i,j),depth(im1,j))
          depthv(i,j)=min(depth(i,j),depth(i,jm1))
          if (depth(i,j).lt.1.) depth(i,j)=0.
          if (depthu(i,j).lt.1.) depthu(i,j)=0.
          if (depthv(i,j).lt.1.) depthv(i,j)=0.
        enddo
      enddo
c
c --- ------------------------------------------------------------------
c --- integrate the stream function with boundary condition strmf(1,1)=0
c --- ------------------------------------------------------------------
c
      strmf(1,1)=0.
c
      do j=2,jdm
        jm1=j-1
        strmf(1,j)=strmf(1,jm1)
     .            -ub(1,jm1)*depthu(1,jm1)*scuy(1,jm1)
      enddo
c
      do i=2,idm
        im1=i-1
        do j=1,jdm
          strmf(i,j)=strmf(im1,j)
     .              +vb(im1,j)*depthv(im1,j)*scvx(im1,j)
        enddo
      enddo
c
c --- ------------------------------------------------------------------
c --- interpolate the streamfunction to the p-point (also smooths)
c --- ------------------------------------------------------------------
c
      do j=1,jdm
        jp1=mod(j,jdm)+1
        do i=1,idm-1
          ip1=i+1
          strmf(i,j)=.25*(strmf(i,j  )+strmf(ip1,j  )
     .                   +strmf(i,jp1)+strmf(ip1,jp1))
        enddo
      enddo
c
      return
c
      end subroutine

      end module m_strmf_eval
