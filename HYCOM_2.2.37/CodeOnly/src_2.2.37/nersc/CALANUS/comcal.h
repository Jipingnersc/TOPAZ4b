      integer,parameter::  readpop=0
!      integer,parameter::  mpop=50000
      integer,parameter::  maxt=365!,offsp=mpop*0.25
      integer,parameter::  mthresh=10
      integer,parameter::  maxs=6,maxa=19
      real,parameter::     mute=0.05,effect=0.2,cmat=100.,cs_cal=40.
      real,parameter::     eggw=0.23,diff_cal=0.1
      real,parameter::	   ksur=0.5,contr=0.5,beam=0.40
      real,parameter::	   sthresh=7.5E-07,red=6.62,card=50. !g Cm-2
      real,parameter::	   pi_cal=3.1415,cchl=40., mld=50.
      real,parameter::     tpdens=363.,tr=0.1,tpvel=0.01,tacw=2.
      real,parameter::     mpvel=0.4,mesl=30.,mesw=50.,mpdens=185
!      real,parameter::     totcop=2.70E+17 
!C Total number adult 50 m tonn wet weight, 1 mg per adult, side 462 Norskehavsboka
      real	zm_cal(0:13,10),zn_cal(0:13,10),zb_cal(0:13,10)
      real      cumn(0:365),copd(100,2),naupd(100)
      real      npar(1500),tacpred(1500,0:1),popn(0:13)
      real      avgraz(365)
      real      numpred(0:13),avpred(0:13),agemolt(0:13),molters(0:13)
      real      numtemp(0:13),avtemp(0:13),bh(12),avegest(365)
      real      mesopred(1000,0:1),ns(12)!,prisk(4,mpop),ns(12)

      real, parameter :: calasc=60.0, calstd=15.0 ! ascent day - is set in blkdat_floats

      common // popn,tacpred,zm_cal
      common // numpred,copd,ns
      common // avpred
      common // cumn,npar!,nc,ec
      common // agemolt,molters,numtemp,avtemp,bh
      common // mesopred,naupd,zb_cal
      common // zn_cal,avegest,avgraz

!CAS Attribute codes
!CAS av0 = 0=dead, 1=alive
!CAS av1 = Stage, 0 is egg, 1-6 nauplia, 7-11 copepodits, 12 adult
!CAS av2 = Internal number in individuals
!CAS av3 = Structural weight  in micrograms
!CAS av4 = Stage longevity
!CAS av5 = Fat energy level in KJ
!CAS av6 = Moult cycle fraction (Egg,N2,N3) | cum. egg number (adult)
!CAS av7 = Total number of eggs in a super individual
!CAS av8 = habd 
!CAS av9 = habe 
!CAS av10 = habn
!CAS av11 = diapause,0 diapause,1 active,2 move down,3 move up
!CAS av12 = Summing of daily growth  
!CAS av13 = Summing of hourly egestion
!CAS av14 = Unique identifier
!CAS av15 = neutral overwintering depth
!CAS av16 = growth rate
!CAS av17 = antall avkom 


!CAS Strategy codes
!CAS sv1 = wud
!CAS sv2 = fsr
!CAS sv3 = afd
!CAS sv4 = owd
!CAS sv5:maxs = ANN weights
!CAS bsv1:nhid+1 = ANN weights
