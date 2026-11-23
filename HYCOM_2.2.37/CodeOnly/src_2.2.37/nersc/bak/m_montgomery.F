      module m_montgomery
      ! computes montgomery potential 
      contains

      subroutine montgomery(m)
      use mod_xc
      use mod_tides
      implicit none
      integer i,j,l,k,m, oldmargin
      real :: oneta(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real :: sumdp, sumth
      include 'common_blocks.h'
      include 'stmt_fns.h'
c
!$OMP PARALLEL DO PRIVATE(j,l,k,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            sumdp = 0.0
            sumth = 0.0
            do k=1,kk
              if (kapref.ne.0) then  !thermobaric
c
c ---   sigma-star is virtual potential density, as defined in 
c ---   Sun et.al. (1999), 'Inclusion of thermobaricity in 
c ---   isopycnic-coordinate ocean models', JPO 29 pp 2719-2729.
c
c ---   use upper interface pressure in converting sigma to sigma-star.
c ---   to avoid density variations in layers intersected by bottom
c
                if     (kapref.gt.0) then
                  thstar(i,j,k,1)=th3d(i,j,k,m)
     &                             +kappaf(temp(i,j,k,m),
     &                                     saln(i,j,k,m),
     &                                        p(i,j,k),
     &                                     kapref)
                else
                  thstar(i,j,k,1)=th3d(i,j,k,m)
     &                             +kappaf(temp(i,j,k,m),
     &                                     saln(i,j,k,m),
     &                                        p(i,j,k),
     &                                     2)
                  thstar(i,j,k,2)=th3d(i,j,k,m)
     &                             +kappaf(temp(i,j,k,m),
     &                                     saln(i,j,k,m),
     &                                        p(i,j,k),
     &                                     kapi(i,j))
                endif !kapref
              else  !non-thermobaric
                thstar(i,j,k,1)=th3d(i,j,k,m)
              endif  !thermobaric:else
c
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,m)
c
              if     (sshflg.ne.0) then
                sumth = sumth + dp(i,j,k,m)*th3d(i,j,k,m)
                sumdp = sumdp + dp(i,j,k,m)
              endif !sshflg
            enddo !k
c
            if     (sshflg.ne.0) then
              sumth = sumth / max( sumdp, onemm )  !vertical mean of th3d
              sumdp = sumdp*qonem * g              !depth(m) * g
              steric(i,j) =  sshgmn(i,j) +
     &                      (sshgmn(i,j) + sumdp) *
     &                      (thmean(i,j) - sumth) /
     &                      (1000.0+thbase+sumth)
            endif !sshflg
c
c ---       store (1+eta) (= p_total/p_prime) in -oneta-
            oneta(i,j)=1.+pbavg(i,j,m)/p(i,j,kk+1)
c
c ---       m_prime in lowest layer:
            montg(i,j,kk,1)=psikk(i,j,1)+
     &        ( p(i,j,kk+1)*(thkk(i,j,1)-thstar(i,j,kk,1))
     &          -pbavg(i,j,m)*(thstar(i,j,kk,1)+thbase) )*thref**2
            if     (kapref.eq.-1) then
              montg(i,j,kk,2)=psikk(i,j,2)+
     &          ( p(i,j,kk+1)*(thkk(i,j,2)-thstar(i,j,kk,2))
     &            -pbavg(i,j,m)*(thstar(i,j,kk,2)+thbase) )*thref**2
            endif !kapref.eq.-1
c
c ---       m_prime in remaining layers:
            do k=kk-1,1,-1
              montg(i,j,k,1)=montg(i,j,k+1,1)+p(i,j,k+1)*oneta(i,j)
     &            *(thstar(i,j,k+1,1)-thstar(i,j,k,1))*thref**2
              if     (kapref.eq.-1) then
                montg(i,j,k,2)=montg(i,j,k+1,2)+p(i,j,k+1)*oneta(i,j)
     &              *(thstar(i,j,k+1,2)-thstar(i,j,k,2))*thref**2
              endif !kapref.eq.-1
            enddo !k
c
c ---       srfhgt (used diagnostically, in mxmyaij and for tidal SAL).
            if     (kapref.ne.-1) then
              montg1(i,j) =                 montg(i,j,1,1)
            else
              montg1(i,j) =      skap(i,j) *montg(i,j,1,1) +
     &                      (1.0-skap(i,j))*montg(i,j,1,2)
            endif !kapref
            srfhgt(i,j) = montg1(i,j) + thref*pbavg(i,j,m)
c
cdiag       if     (sshflg.ne.0) then
cdiag         if     (itest.gt.0 .and. jtest.gt.0) then
cdiag           write (lp,'(i9,2i5,3x,a,2f12.6,f12.2)')
cdiag&            nstep,itest+i0,jtest+j0,
cdiag&              'sssh =',
cdiag&              steric(i,j),sshgmn(i,j),sumdp
cdiag           write (lp,'(i9,2i5,3x,a,3f12.6)')
cdiag&            nstep,itest+i0,jtest+j0,
cdiag&              'thmn =',
cdiag&              sumth,thmean(i,j),1000.0+thbase+sumth
cdiag           write (lp,'(i9,2i5,3x,a,3f12.6)')
cdiag&            nstep,itest+i0,jtest+j0,
cdiag&              'ssh  =',
cdiag&              srfhgt(i,j),steric(i,j),srfhgt(i,j)-steric(i,j)
cdiag         endif !test
cdiag       endif !sshflg
c
c ---       tidal corrections, note that these are not part of montg
c ---       but are included here to simplify the pressure gradient
            if     (tidflg.gt.0 .and. sshflg.eq.0) then !tides
              do k=1,kk
                montg(i,j,k,1)=montg(i,j,k,1)
     &                            -g*etide(i,j)-tidsal*srfhgt(i,j)
                if     (kapref.eq.-1) then
                  montg(i,j,k,2)=montg(i,j,k,2)
     &                              -g*etide(i,j)-tidsal*srfhgt(i,j)
                endif !kapref.eq.-1
              enddo !k
            elseif (tidflg.gt.0 .and. sshflg.ne.0) then !tides
              do k=1,kk
                montg(i,j,k,1)=montg(i,j,k,1)
     &                            -g*etide(i,j)
     &                            -tidsal*(srfhgt(i,j)-steric(i,j))
                if     (kapref.eq.-1) then
                  montg(i,j,k,2)=montg(i,j,k,2)
     &                              -g*etide(i,j)
     &                              -tidsal*(srfhgt(i,j)-steric(i,j))
                endif !kapref.eq.-1
              enddo !k
            endif !tides (sshflg)
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c





      end subroutine montgomery
      end module m_montgomery
