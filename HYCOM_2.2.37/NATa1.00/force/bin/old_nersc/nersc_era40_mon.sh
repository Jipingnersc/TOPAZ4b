#!/bin/bash
#
set -x
pget=cp
pput=cp

# Experiment number needed
if [ $# -ne 1 ] ; then
   echo " $(basename $0) needs experiment version as input (ex 01.0)"
   exit
fi
export X=$1


# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd -
source ${BASEDIR}/EXPT.src $X || { echo "Could not source ${BASEDIR}/EXPT.src" ; exit 1 ; }

# Check that pointer to MSCPROGS is set
if [ -z ${MSCPROGS} ] ; then
   echo "MSCPROGS Environment not set "
   exit
else
   if [ ! -d ${MSCPROGS} ] ; then
      echo "MSCPROGS not properly set up"
      echo "MSCPROGS not a directory at ${MSCPROGS}"
      exit
   fi
fi


#
#
# --- Create HYCOM atmospheric forcing (era40 climatology)
#
#
# --- S is scratch directory,
# --- D is permanent directory,
#
D=$BASEDIR/force/nersc_era40/
S=$D/SCRATCH
mkdir -p $S
mkdir -p $S/Data
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

#
# --- Input.
#
touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b #Data/sss_nodc.dat Data/sst_nodc.dat
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b #Data/sss_nodc.dat Data/sst_nodc.dat
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b     || { echo "Could not get regional.grid.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a     || { echo "Could not get regional.grid.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a  || { echo "Could not get regional.depth.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b  || { echo "Could not get regional.depth.b file " ; exit 1 ; }
#${pget} ${BASEDIR}/expt_${X}/blkdat.input blkdat.input      || { echo "Could not get file blkdat.input " ; exit 1 ; }
#${pget} ${BASEDIR}/expt_${X}/infile.in infile.in            || { echo "Could not get file infile.in " ; exit 1 ; }


$MSCPROGS/bin_setup/forfun_nersc-2.2 month era40 ||  { echo "Error running forfun_nersc-2.2 " ; exit 1 ; }

# The nersc era40 forcing is region-independent 
mkdir -p $D/$E  || { echo " Could not create data  destination dir $D/$E" ; exit 1;}
for i in forcing.*.[ab] ; do
   new=$(echo $i | sed "s/^forcing\.//")
   mv $i $D/$E/$new
done




