!
!-----------------------------------------------------
!
SUBROUTINE GETR(slig,D,VISRAN,size,diff,str)
!
!-----------------------------------------------------
!
!Obtain visual range by solving the non-linear equation 
!by means of Newton-Raphson iteration and derivation in 
!subroutine DERIV.
!
IMPLICIT  NONE
INTEGER   IEND,AS,IER,I 
INTEGER, INTENT(IN) :: D
REAL, INTENT(IN) :: slig,size,diff,str
REAL, INTENT(OUT) ::  VISRAN
REAL      EPS,RST,TOL,TOLF,F1,FDER,DX
include 'comcal.h'  
!
!Input parameters 
!RST    : start value of VISRAN (m) 
!KLIG   : diffuse light attenuation coefficient (m-1) 
!CLIG   : beam attenuation koefficient (m-1) 
!ksur   : loss through surface (dim.less) 
!SLIG   : above-surface light intensity (uE m-2 s-1) 
!D		: current Depth (m) 
!contr  : prey light absorption-reflectance coeff.
!SIZE   : prey size (m) 
!sthresh  : sensitivity threshold for eye  (uE m-2 s-1) 
!
!Output parameters 
!F1     : function value of equation in DERIV 
!FDER   : the derivative of the function 
!VISRAN : the predator's visual range (when F1=0) 
!
!Initial guess of visual range (RST)
CALL EASYR(D,RST,slig,size,diff,str)
!
!Upper boundary of allowed error of visual range 
 EPS =.0000001 
!Maximum number of iteration steps 
 IEND = 100  
 IER = 0 
!Prepare iteration 
 VISRAN = RST 
 TOL = VISRAN
!
 CALL DERIV(D,VISRAN,F1,FDER,slig,size,diff,str)
! 
 TOLF = 100.* EPS 
! 
!Start iteration loop 
 DO 6 I = 1, IEND 
  IF (F1) 1, 7, 1 
! 
!Equation is not satisfied by VISRAN 
1  IF (FDER) 2, 8, 2 
! 
!Iteration is possible 
2 DX = F1/FDER 
  VISRAN = VISRAN-DX
! 
!Test on allowed range 
 IF (VISRAN.LT.0.) GOTO 9 
! 
 TOL = VISRAN 
 CALL DERIV(D,VISRAN,F1,FDER,slig,size,diff,str)
! 
!Test on satisfactory accuracy 
   TOL = EPS 
   AS = ABS(VISRAN) 
   IF (AS-1.) 4, 4, 3 
3   TOL = TOL*AS 
4  IF (ABS(DX)-TOL) 5, 5, 6 
5  IF (ABS(F1)-TOLF) 7, 7, 6 
6  CONTINUE 
! 
!goto 15		   !NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNNBNNBNBNNBNBNBN
! 
!     No convergence after IEND steps.Error return.
      IER = 1
7     RETURN 
!     Return in case of zero divisor 
8     IER = 2
      RETURN 
!     VISRAN out of allowed range (negative) 
9     IER = 3
!       
!15 continue	   !NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNNBNNBNBNNBNBNBN
!
END SUBROUTINE GETR
!
!-----------------------------------------------------
!
SUBROUTINE EASYR(D,RST,slig,size,diff,str)
!
!-----------------------------------------------------
!
!Obtain a first estimate of visual range by using a linear expression 
! 
IMPLICIT NONE
INTEGER     D
REAL        VISRAN,R2
REAL        RST,slig,size,diff,str
include  'comcal.h'  
!
!KLIG   : diffuse light attenuation coefficient (m-1) 
!ksur   : loss through surface (dim.less) 
!SLIG   : above-surface light intensity (uE m-2 s-1) 
!D		: current Depth (m) 
!contr  : prey light absorption-reflectance coeff.
!SIZE   : prey size (m) 
!sthresh  : sensitivity threshold for eye  (uE m-2 s-1) 
!VISRAN : the predator's visual range (when F1=0) 
!R2     : the square of the visual range (m2) 
!
!R2 = EXP(-D*DIFF)*ksur*slig*ABS(contr)*PI_CAL*SIZE*SIZE/sthresh*1. 
R2 = EXP(-D*DIFF)*ksur*slig*ABS(contr)*PI_CAL*SIZE*SIZE/STR*1. 
VISRAN = SQRT(R2) 
RST = VISRAN
!
END SUBROUTINE EASYR 
!
!-----------------------------------------------------
!
SUBROUTINE DERIV(D,VISRAN,F1,FDER,slig,size,diff,str)
!
!-----------------------------------------------------
!
!Derivation of equation for visual range of a predator 
! 
IMPLICIT NONE
INTEGER   D 
REAL      VISRAN,FR1,FR2,EB,EMAX,KE 
REAL      F1,FDER,FR11,V1,slig,size,diff,str
include 'comcal.h'  
!  
!Input parameters 
!RST    : start value of VISRAN (m) 
!KLIG   : diffuse light attenuation coefficient (m-1) 
!CLIG   : beam attenuation koefficient (m-1) 
!ksur   : loss through surface (dim.less) 
!SLIG   : above-surface light intensity (uE m-2 s-1) 
!D  : current D (m) 
!contr  : prey light absorption-reflectance coeff.
!SIZE   : prey size (m) 
!sthresh : sensitivity threshold for eye  (uE m-2 s-1) 
! 
!Output parameters 
!F1     : function value of equation in DERIV 
!FDER   : the derivative of the function 
!VISRAN : the predator's visual range (when F1=0) 
! 
V1 = VISRAN 
Ke = 5. !uE m-2 s-1	 3.23E-6 
Emax = 50. !0 !org. = 50
!
EB = ksur*slig*EXP(-DIFF*D)
FR1 = 2*LOG(V1)+VISRAN*beam
!FR11 = (EB/(Ke+EB)*1.)*Emax*ABS(contr)*PI_CAL*SIZE*SIZE/sthresh*1. 
FR11 = (EB/(Ke+EB)*1.)*Emax*ABS(contr)*PI_CAL*SIZE*SIZE/STR*1. 
FR2 = LOG(FR11) 
F1 = FR2 - FR1 
FDER = -beam - 2./VISRAN 
!  
END SUBROUTINE DERIV
