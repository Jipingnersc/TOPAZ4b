#!/bin/bash
#
# Routine sets up hycom tidal data (CSRobc_elev.dat) in experiment directory $E
# Programs used is "csr2mod_GE"
#
# Final files are placed in $BASEDIR/tides_nersc/$E/
#
# REQUIRES: a working topography for experiment $E (regional.grid/regional.depth)
#           CSR_PATH to be set, points to location of CSR dataset and input.CSR
#           If input.CSR is present in dir experiment $E, that will be preferred 
#           so for local changes modify that file.


pget=cp
pput=cp

# Experiment also needs topography
if [ $# -ne 1 ] ; then
   echo "Script will set up boundary tidal data for the model using the"
   echo "CSR tidal database. Input is experiment number"
   echo
   echo "Usage:"
   echo " $(basename $0) experiment_version"
   echo "Example:"
   echo " $(basename $0) 01.0"
   exit
fi
export X=$1

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
HERE=$(pwd)
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd $HERE



source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }

# Check that CSR_PATH is set (usually done in REGION.src)
if [ -z ${CSR_PATH} ] ; then
   echo "environment variable CSR_PATH not set "
   exit
fi
if [ -z ${MSCPROGS} ] ; then
   echo "environment variable MSCPROGS not set "
   exit
fi



#
#
# --- Create CSR tidal forcing fields for nersc version 
#
#
# --- S is scratch directory,
# --- D is permanent directory,
#
D=$BASEDIR/tides_nersc/
S=$D/SCRATCH
TARGET=$D/$E
mkdir -p $S
mkdir -p $TARGET
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

#
# --- Input.
#
touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b input.CSR
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b input.CSR
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b    || { echo "Could not get regional.grid.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a    || { echo "Could not get regional.grid.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a || { echo "Could not get regional.depth.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b || { echo "Could not get regional.depth.b file " ; exit 1 ; }

# For CSR input 
${pget} ${CSR_PATH}/input.CSR  .  || { echo "Could not get input.CSR " ; exit 1 ; }
echo "Using input.CSR from $CSR_PATH"


#Run csr2mod
$MSCPROGS/src/Tides_CSR/csr2mod_GE

cp CSRobc_elev.dat $TARGET
cp input.CSR $TARGET
cp csrtide.nc $TARGET

echo
echo "CSRobc_elev.dat  in $TARGET"
echo "Active input.CSR in $TARGET"
echo "Diagnostics csrtide.nc in $TARGET"
