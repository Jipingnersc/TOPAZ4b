#!/bin/bash
#
echo " Script needs setup"
exit
pget=cp
pput=cp
set -x
export X="NOEXPERIMENT" # dummy experiment

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd -

source ${BASEDIR}/EXPT.src $X || { echo "Could not source ${BASEDIR}/EXPT.src" ; exit 1 ; }

# Check that pointer to HYCOM_ALL is set (from EXPT.src)
if [ -z ${HYCOM_ALL} ] ; then
   echo "Environment not set "
   exit
else
   if [ ! -d ${HYCOM_ALL} ] ; then
      echo "HYCOM_ALL not properly set up"
      echo "HYCOM_ALL not a directory at ${HYCOM_ALL}"
      exit
   fi
fi

#
# --- Create a wind stress data file.
#
#
# --- S is scratch directory,
# --- D is permanent directory,
#
D=$BASEDIR/force/ecmwf-reanal_7993-mn6hr/
S=$D/SCRATCH
mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}



#
# --- Input.
#
touch      fort.44 fort.44A fort.45 fort.45A fort.71
/bin/rm -f fort.44 fort.44A fort.45 fort.45A fort.71
${pget} $BASEDIR/force/offset/tauewd_zero.b fort.44   || { echo "Can not get wind offset" ; exit 1; }
${pget} $BASEDIR/force/offset/tauewd_zero.a fort.44A  || { echo "Can not get wind offset" ; exit 1; }
${pget} $BASEDIR/force/offset/taunwd_zero.b fort.45   || { echo "Can not get wind offset" ; exit 1; }
${pget} $BASEDIR/force/offset/taunwd_zero.a fort.45A  || { echo "Can not get wind offset" ; exit 1; }
${pget} $ERA15/era15-sea_1979-1993_mon_strblk.D fort.71   || { echo "Can not get era15 stress data" ; exit 1; }

touch regional.grid.a regional.grid.b
rm regional.grid.a regional.grid.b
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b  ||  { echo "Can not get regional.grid.b" ; exit 1; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a  ||  { echo "Can not get regional.grid.a" ; exit 1; }

#cp ${D}/../../../ALL/force/src/wi_100_co . &
cp ${HYCOM_ALL}/force/src/wi .  ||  { echo "Can not get wi executable " ; exit 1; }
wait
#chmod a+rx wi_100_co
chmod a+rx wi
/bin/rm -f fort.1[0123]*

export FOR010A=fort.10A
export FOR011A=fort.11A
export FOR012A=fort.12A
export FOR013A=fort.13A
export FOR044A=fort.44A
export FOR045A=fort.45A

/bin/rm -f core
touch core
/bin/rm -f fort.10 fort.11 fort.12 fort.13 fort.46 fort.47 fort.48
ISPEED=2
IUSTAR=1
WSCALE=1.0
echo "KAL -- routine disabled until we have a clear understanding of wscale " && exit 1
#./wi_100_co <<E-o-D
./wi <<E-o-D
 &WWTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'ERA15 monthly, MKS',
 /
 &WWTIME
  SPDMIN = 0.0,
  WSCALE = $WSCALE,
 /
 &WWFLAG
  IGRID  = 1,
  ISPEED = $ISPEED,
  IUSTAR = $IUSTAR,
  INTERP = 1,
  IWFILE = 1,
 /
E-o-D
#
# --- Output.
#
${pput} fort.10  ${D}/tauewd.b
${pput} fort.10A ${D}/tauewd.a
${pput} fort.11  ${D}/taunwd.b
${pput} fort.11A ${D}/taunwd.a
${pput} fort.12  ${D}/wndspd.b
${pput} fort.12A ${D}/wndspd.a
${pput} fort.13  ${D}/ustar.b
${pput} fort.13A ${D}/ustar.a
#
# --- Delete all files.
#
#/bin/rm -f *
echo
echo "Wind climatology copied to $D"
