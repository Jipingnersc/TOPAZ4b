!---------------------------------------------------------------------------------

subroutine preddens

!---------------------------------------------------------------------------------
  implicit         none
  integer	         d,skrivut
  real	         nevn,ze,norm,expo,sd,sumnp
  include          'comcal.h'!

!AS
  real, parameter :: backgmp=0.015, backgtp=0.015

!Vertical distribution of mesopelagic fish density
!Day distribution
  sd = 0.1*1000
  nevn = sqrt(6.28318*sd)
  norm =	(1./nevn)
  sumnp = 0.
  do d = 1,1000
     ze = ((d-75)**2.)/(2*sd**2.)
     expo =  exp(-ze)
     npar(d) = norm*expo
     sumnp = sumnp + npar(d)
  enddo
  do d = 1,1000
     mesopred(d,1) = backgmp + npar(d)*mpdens/sumnp*1. ! nmax*1.
  enddo

!Vertical distribution of mesopelagic fish density
!Night distribution
  sd = 0.08*1000
  nevn = sqrt(6.28318*sd)
  norm =	(1./nevn)
  sumnp = 0.
  do d = 1,1000
     ze = ((d-25)**2.)/(2*sd**2.)
     expo =  exp(-ze)
     npar(d) = norm*expo
     sumnp = sumnp + npar(d)
  enddo
  do d = 1,1000
     mesopred(d,0) = backgmp + npar(d)*mpdens/sumnp*1. !nmax*1. 
  enddo

!Vertical distribution of tactile predator density
!Day distribution
  sd = 0.2*1000
  nevn = sqrt(6.28318*sd)
  norm =	(1./nevn)
  sumnp = 0.
  do d = 1,1000
     ze = ((d-250)**2.)/(2*sd**2.)
     expo =  exp(-ze)
     npar(d) = norm*expo
     sumnp = sumnp + npar(d)
  enddo
  do d = 1,1000
     tacpred(d,1) = backgtp + npar(d)*tpdens/sumnp*1.
  enddo

!Vertical distribution of tactile predator density
!Night distribution
  sd = 0.2*1000
  nevn = sqrt(6.28318*sd)
  norm =	(1./nevn)
  sumnp = 0.
  do d = 1,1000
     ze = ((d-50)**2.)/(2*sd**2.)
     expo =  exp(-ze)
     npar(d) = norm*expo
     sumnp = sumnp + npar(d)
  enddo
  do d = 1,1000
     tacpred(d,0) = backgtp + npar(d)*tpdens/sumnp*1.
  enddo

  skrivut = 0
  if(skrivut==1) then 
     open(27, file = 'out/'//'hormon.txt')
     do d = 1,1000
        write(27,2)d,tacpred(d,0),mesopred(d,0),tacpred(d,1),mesopred(d,1)
     enddo
     close(27)
     stop
  2 format(i4,x,4(f10.3,x))
  endif


end subroutine preddens
