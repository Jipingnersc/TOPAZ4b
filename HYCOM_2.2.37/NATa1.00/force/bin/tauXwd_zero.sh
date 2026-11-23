#!/bin/bash
#
echo " Script needs further setup"
exit
set -x
pget=cp
pput=cp
#
# --- create HYCOM zero wind stress files.
#

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd -
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }

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
# --- S is scratch directory,
# --- D is permanent directory,
#
D=${BASEDIR}/force/offset/
S=$D/SCRATCH
mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}
#
# --- Input.
#
touch regional.grid.a regional.grid.b
rm    regional.grid.a regional.grid.b
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b || { echo "Could not get regional.grid.a" ; exit 1 ;}
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a || { echo "Could not get regional.grid.a" ; exit 1 ;}
cp ${HYCOM_ALL}/force/src/off_zero .  || { echo "Could not get executable off_zero" ; exit 1 ;}
wait
chmod a+rx off_zero

/bin/rm -f fort.10 fort.10A
/bin/rm -f fort.11 fort.11A
export FOR010A=fort.10A
export FOR011A=fort.11A
/bin/rm -f core
touch core
./off_zero
#
# --- Output.
#
${pput} fort.10  ${D}/tauewd_zero.b
${pput} fort.10A ${D}/tauewd_zero.a
${pput} fort.11  ${D}/taunwd_zero.b
${pput} fort.11A ${D}/taunwd_zero.a
#
# --- Delete all files.
#
#/bin/rm -f *
