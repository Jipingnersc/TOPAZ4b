#
set -x
C
C --- Create model interpolated winds for HYCOM. Use cubic cpline interpolation.
C --- Zero offset file.
C ---
C --- prebuild this script similar to a model script
# --- awk -f 284.awk y01=098 ab=a 284W.com > 284w098a.com
C
C --- Preamble, script keys on O/S name.
C
C
C --- pget, pput "copy" files between scratch and permanent storage.
C --- Can both be cp if the permanent filesystem is mounted locally.
C
export pget=cp
export pput=cp

C
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C --- N is data-set name, e.g. ec10m-reanal
C --- W is permanent native pcip directory
C
#setenv E 016
#setenv P hycom/ATLb2.00/expt_01.6/data
#setenv D ~/$P
#setenv N ec10m-oper
#
#KAL - these should be accessible from parent script - check this
if [ "$P" = "" ] ; then 
   echo "Experiment data not sourced"
   exit 1
fi
exit 
mkdir -p $S/wind
cd       $S/wind
C
C --- For whole year runs.
C ---   ymx number of years per model run.
C ---   Y01 initial model year of this run.
C ---   YXX is the last model year of this run, and the first of the next run.
C ---   A and B are identical, typically blank.
C --- For part year runs.
C ---   A is this part of the year, B is next part of the year.
C ---   Y01 initial model year of this run.
C ---   YXX year at end of this part year run.
C ---   ymx is 1.
C --- Note that these variables and the .awk generating script must
C ---  be consistant to get a good run script.
C
C --- For winds, only Y01 and A are used.
C
C --- One year spin-up run.
C
#@ ymx =  1
#C
#setenv A "a"
#setenv B "b"
#setenv Y01 "099"
#C
#switch ("${B}")
#case "${A}":
#    setenv YXX `echo $Y01 $ymx | awk '{printf("%03d", $1+$2)}'`
#    breaksw
##case "a":
#    setenv YXX `echo $Y01 | awk '{printf("%03d", $1+1)}'`
#    breaksw
#default:
#    setenv YXX $Y01
#endsw
#C
#echo "Y01 =" $Y01 "YXX = " $YXX  "A =" ${A} "B =" ${B}

C
C --- time limits.
C --- use "LIMITI" when starting a run after day zero.
C
#setenv TS `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $1}'`
#setenv TM `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $2}'`
#C
#echo "TS =" $TS "TM =" $TM


C
C --- input files from file server.
C
touch      fort.44 fort.44a fort.45 fort.45a
#/bin/rm -f fort.44 fort.44a fort.45 fort.45a
if [ ! -s  fort.44 ]; then
  ${pget} ${D}/../../force/offset/tauewd_zero.b fort.44  &
fi    
if [ ! -s  fort.44a ]; then
  ${pget} ${D}/../../force/offset/tauewd_zero.a fort.44a &
fi
if [ ! -s  fort.45 ]; then
  ${pget} ${D}/../../force/offset/taunwd_zero.b fort.45  &
fi
if [ ! -s  fort.45a ]; then
  ${pget} ${D}/../../force/offset/taunwd_zero.a fort.45a &
fi
C
touch regional.grid.a regional.grid.b
if [ ! -s regional.grid.b ]; then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
fi
if [ ! -s regional.grid.a ]; then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
fi

C
C KAL - retrieve forcing info from infile.in
C
FRC=$(cat infile.in | head -n 6 | tail -n1 | cut -c1-5)
CLM=$(cat infile.in | head -n 6 | tail -n1 | cut -c7-11)

if [ "$FRC" == "era40" -a "$CLM" = "era40" ] ; then
   C
   C KAL - generate just-in-time forcing from era40 (NERSC)
   C
   export FORFUN_WIND=yes
   ln -s ${NERSC_ERA40} .
   ln -s ${NERSC_ERA40_CLIM} .
   ~knutali/Progs/Setup/Synoptic_Forcing/Build/forfun_nersc $FRC $CLM
else
   C
   C KAL - generate just-in-time forcing from hycom ec10m
   C
   echo "For now only nersc HF ERA40 forcing is supported"
   C
   touch      fort.71 fort.72
   if [ ! -s fort.71 ]; then
     ${pget} ${W}/ds111.1_ec10m-oper_9801_06hr_ws.d fort.71 &
   fi
   if [ ! -s fort.72 ]; then
     ${pget} ${W}/ds111.1_ec10m-oper_0205_06hr_ws.d fort.72 &
   fi

   C
   C --- executable
   C
   /bin/cp /u/home/wallcraf/hycom/ALL/force/src/wi . &
   wait
   chmod ug+rx wi
   ls -laFq
   C
   C --- NAMELIST input.
   C
   cat <<E-o-D  >! fort.05i
    &WWTITL
     CTITLE = '1234567890123456789012345678901234567890',
     CTITLE = 'ECMWF, 10m DS111.1, 6-hr, MKS',
    /
    &WWTIME
     WSCALE =   1.0,  !already in MKS
     SPDMIN =   0.0,
     WSTART = ${TS},
     TSTART = ${TS},
     TMAX   = ${TM},
    /
    &WWFLAG
     IGRID  = 1,  !0=p; 1=u&v; 2=p;
     ISPEED = 2,  !0:none; 1:const; 2:kara; 3:coare
     INTERP = 1,  !0:bilinear; 1:cubic spline
     INTMSK = 0,  !0:no mask; 1:land/sea=0/1; 2:l/s=1/0;
     IFILL  = 3,  !0,1:tx&ty; 2,3:magnitude; 1,3:smooth;
     IOFILE = 0,
     IWFILE = 4,
    /
E-o-D
   C
   C --- run the wind interpolation.
   C
   touch fort.10 fort.10a
   /bin/rm -f fort.1[012] fort.1[012]a
   C
   export FOR010A=fort.10a
   export FOR011A=fort.11a
   export FOR012A=fort.12a
   C
   export FOR044A=fort.44a
   export FOR045A=fort.45a
   C
   /bin/rm -f core
   touch core
   ./wi < fort.05i
   C
   C --- Output.
   C
   /bin/mv fort.10  tauewd_${Y01}${A}.b
   /bin/mv fort.10a tauewd_${Y01}${A}.a
   /bin/mv fort.11  taunwd_${Y01}${A}.b
   /bin/mv fort.11a taunwd_${Y01}${A}.a
   /bin/mv fort.12  wndspd_${Y01}${A}.b
   /bin/mv fort.12a wndspd_${Y01}${A}.a
   C
   if [ -d ./SAVE ]; then
     ln tauewd_${Y01}${A}.a ./SAVE/tauewd_${Y01}${A}.a
     ln tauewd_${Y01}${A}.b ./SAVE/tauewd_${Y01}${A}.b
     ln taunwd_${Y01}${A}.a ./SAVE/taunwd_${Y01}${A}.a
     ln taunwd_${Y01}${A}.b ./SAVE/taunwd_${Y01}${A}.b
     ln wndspd_${Y01}${A}.a ./SAVE/wndspd_${Y01}${A}.a
     ln wndspd_${Y01}${A}.b ./SAVE/wndspd_${Y01}${A}.b
   fi
   fi
fi
C
C  --- END OF JUST IN TIME WIND GENERATION SCRIPT.
C
