#!/bin/bash
#
# Routine sets up hycom tidal data (FESobc_elev.dat) in experiment directory $E
# Programs used is "fes2mod" and "fes2nc" (latter for diagnostics
#
# Final files are placed in $BASEDIR/tides_nersc/$E/
#
# REQUIRES: a working topography for experiment $E (regional.grid/regional.depth)
#           FES_PATH to be set, points to location of FES dataset. 
#           Number of constituents (8) is hardcoded in fes2mod/fes2nc


pget=cp
pput=cp

# Experiment also needs topography
if [ $# -ne 1 ] ; then
   echo "Script will set up boundary tidal data for the model using the"
   echo "FES tidal database. Input is experiment number"
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

# Check that FES_PATH is set (usually done in EXPT.src)
if [ -z ${FES_PATH} ] ; then
   echo "environment variable FES_PATH not set "
   exit
fi
if [ -z ${MSCPROGS} ] ; then
   echo "environment variable MSCPROGS not set "
   exit
fi



#
#
# --- Create tidal forcing fields for nersc version 
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
touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b 
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b 
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b    || { echo "Could not get regional.grid.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a    || { echo "Could not get regional.grid.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a || { echo "Could not get regional.depth.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b || { echo "Could not get regional.depth.b file " ; exit 1 ; }

#Run fes2mod 
#$MSCPROGS/src/Tides_FES/nersc_src/fes2mod
#Fanf:you need to use gnu and not pgi to compile fes2mod
#Fanf:module swap PrgEnv-pgi/2.0.62 PrgEnv-gnu 
$MSCPROGS/src/Tides_FES/nersc_src/fes2mod

#Run fes2nc for diagnostics
echo "Creating diagnostics in $S/fes2004.nc"
echo "(may take some time - note dimensions are transposed)"
$MSCPROGS/src/Tides_FES/nersc_src/fes2nc > /dev/null

cp FESobc_elev.dat fes2004.nc $TARGET

echo
echo "Tidal input data for hycom and diagnostics file placed in $TARGET"
