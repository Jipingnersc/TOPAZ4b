#!/bin/bash

# KAL - get X from input
if [ $# -ne 1 ] ; then
   echo "This script will set up some relaxation files needed when running"
   echo "with the old option for climatology. Files for SSS and SST relaxation are created"
   echo "Input is experiment number (ex 01.0)"
   exit
fi
export X=$1


# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd -
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }


# Check that pointer to MSCPROGS is set (from EXPT.src)
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


CDIR=${BASEDIR}/relax/$E/
SCRATCH=$CDIR/SCRATCH

mkdir -p $SCRATCH
mkdir -p $SCRATCH/Data
cd       $SCRATCH || { echo " Could not descend scratch dir $SCRATCH" ; exit 1;}

touch blkdat.input regional.grid.a regional.grid.b regional.depth.a regional.depth.b  Data/sst_nodc.dat Data/sss_nodc.dat
rm    blkdat.input regional.grid.a regional.grid.b regional.depth.a regional.depth.b  Data/sst_nodc.dat Data/sss_nodc.dat

cp ${BASEDIR}/expt_${X}/blkdat.input . || { echo "Could not copy blkdat.input " ; exit 1 ;}

cp ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b || { echo "Couldnt get depth_${R}_${T}.b " ; exit 1 ;}
cp ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a  || { echo "Couldnt get depth_${R}_${T}.a " ; exit 1 ;}

cp ${BASEDIR}/topo/regional.grid.b regional.grid.b || { echo "Couldnt get regional.grid.b " ; exit 1 ;}
cp ${BASEDIR}/topo/regional.grid.a regional.grid.a || { echo "Couldnt get regional.grid.b " ; exit 1 ;}

cp $MSCPROGS/data/Data/saln_nodc_srf_1x1deg.dat  $SCRATCH/Data || { echo "Couldnt get old levitus sss " ; exit 1 ;}
cp $MSCPROGS/data/Data/temp_nodc_srf_1x1deg.dat  $SCRATCH/Data || { echo "Couldnt get old levitus sst " ; exit 1 ;}

touch old_levini ;  rm old_levini
cp $MSCPROGS/bin_setup/old_levini  $SCRATCH || { echo "Couldnt get old_levitus " ; exit 1 ;}

./old_levini

# Move data to experiment directory
mv Data/sss_nodc.dat $CDIR
mv Data/sst_nodc.dat $CDIR
