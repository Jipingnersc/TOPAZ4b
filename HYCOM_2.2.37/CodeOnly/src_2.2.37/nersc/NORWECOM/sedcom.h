!
! Sediment common declarations
!
      REAL TTAU(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      LOGICAL WAVES,SEDIMENT
      COMMON/COMLOG/WAVES,SEDIMENT
      COMMON/COMSED/TTAU
!
      REAL,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm) :: SIS,OXY
     & ,YEL,SED
      COMMON/BIOSED/SIS,OXY,YEL,SED
!      REAL DETMAX,DETBUL,DENITFR,DETBU(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),DENIT(IM,JM)
!      COMMON/BOTPRO/DETMAX,DETBUL,DENITFR,DETBU,DENIT 
      REAL DETMAX,DETBUL,DENITFR
      REAL,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: DETBU,DENIT
      REAL SISBUL,SISMAX
      REAL,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)::SISBU
      COMMON/BOTPRO/DETBU,DENIT,SISBU
      REAL,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: TOTSED
      REAL FACTW,DETW,SISW,FLAW,DIAW,NIW,SIW
      COMMON/BOTW/TOTSED
      DATA  FACTW/1000./
      PARAMETER (DETW=3555.,FLAW=3555.,DIAW=4395.,NIW=224.,SIW=28.    
     &    ,SISW=60)
!#if defined (DETPHO)
      REAL DETPW, SRDETP, DETPMAX, DETPBUL
      DATA DETPW/355./,SRDETP/3./
      DATA DETPMAX/10000./,DETPBUL/63./
      REAL,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: DETPBU
      COMMON DETPBU 
!#endif
!
! Atomic weights
!   COXY = oxygen
!   COXY2 = oxygen gas (O2)
!   CCAR = carbone
!
      REAL COXY,COXY2,CCAR
      DATA COXY/16.00/,COXY2/32.00/,CCAR/12.01/
!
!   SCC(1)    : {O}/{N} FOR EACH CELL PRODUCED (weigth) 
!   SCC(2)    : {C}/{N} IN A CELL (weigth)
!   SCC(3)    : {O}/{N} FOR THE DENITRIFICATION
!   SCC(4)    : RATE OF DECOMPOSITION OF SILICATE SHELL (SIS) (180 DAYS) 
!   SCC(4)    : RATE OF DECOMPOSITION OF SILICATE SHELL (SIS) (40 DAYS) 
!   SCC(5)    :  MIN RATE OF LEAKAGE OF NUTRIENTS FROM SEDIMENT (10 DAYS)
!   SCC(6)    :  MAX RATE OF LEAKAGE OF NUTRIENTS FROM SEDIMENT (6 hours)
!   SCC(7)    :  DENITRIFICATION RATE (5 DAYS) (NOT USED)
!   SCC(8)    :  BURIAL RATE 120 DAYS
!
!
      REAL      SCC(8)
      DATA      SCC(1)/19.71/,SCC(2)/5.68/,SCC(3)/3.42/,                
     &          SCC(4)/1.45E-7/,                                        
     &          SCC(5)/1.157E-6/,SCC(6)/4.63E-5/,SCC(7)/2.3E-6/,        
     &          SCC(8)/9.65E-8/
      DATA      DENITFR/0.3/,DETMAX/100000./,DETBUL/630./
      DATA      SISMAX/100000./,SISBUL/2570./
!
!
!
!AS      INTEGER DAY, PARAMETER (DAY=86400)
      REAL SRDET,SRFLA,SRDIAMAX,SRDIAMIN,SRSIS,A3L,SIB
      REAL SRSED,SRYEL
      REAL SEDMIN,YELMIN
      DATA SEDMIN/0.0/,YELMIN/0.0/
      PARAMETER(A3L=2.7,SIB=1.0,SRDET=3.,SRFLA=0.25,          
     &          SRSED=1.0,SRDIAMAX=3.0,SRDIAMIN=0.3,SRSIS=3.0,          
     &          SRYEL=7.5)
!#ifdef CHATONELLA
!      REAL SRCHA
!      PARAMETER(SRCHA=12.0)
!#endif
