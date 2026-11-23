#!/bin/bash

set -x

# Experiment also needs topography
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


CDIR=$BASEDIR/force/nersc_old/
SCRATCH=${CDIR}/SCRATCH
mkdir -p $SCRATCH
mkdir -p $SCRATCH/Data
cd       $SCRATCH || { echo " Could not descend scratch dir $SCRATCH" ; exit 1;}

touch blkdat.input Data/sss_nodc.dat Data/sst_nodc.dat infile.in regional.depth.a regional.depth.b regional.grid.a regional.grid.b
rm    blkdat.input Data/sss_nodc.dat Data/sst_nodc.dat infile.in regional.depth.a regional.depth.b regional.grid.a regional.grid.b

#cp ../../../expt_${X}/blkdat.input .            || { echo "Could not copy blkdat.input " ; exit 1 ;}
cp ${CDIR}/../../relax/${E}/sss_nodc.dat Data/  || { echo "Could not get file sss_nodc.dat (need to run old_lev_nersc.sh in relax first)" ; exit 1 ; }
cp ${CDIR}/../../relax/${E}/sst_nodc.dat Data/  || { echo "Could not get file sst_nodc.dat (need to rub old_lev_nersc.sh in relax first)" ; exit 1 ; }
#cp ${CDIR}/../../expt_${X}/infile.in infile.in  || { echo "Could not get file infile.in " ; exit 1 ; }
cp ${CDIR}/../../topo/regional.grid.b regional.grid.b     || { echo "Could not get regional.grid.b file " ; exit 1 ; }
cp ${CDIR}/../../topo/regional.grid.a regional.grid.a     || { echo "Could not get regional.grid.a file " ; exit 1 ; }
cp ${CDIR}/../../topo/depth_${R}_${T}.a regional.depth.a  || { echo "Could not get regional.depth.a file " ; exit 1 ; }
cp ${CDIR}/../../topo/depth_${R}_${T}.b regional.depth.b  || { echo "Could not get regional.depth.b file " ; exit 1 ; }



echo "Copying old Forcing to $SCRATCH"
cp -r $MSCPROGS/data/Data $SCRATCH
cp -r $MSCPROGS/data/Data/fields.inp* $SCRATCH
$MSCPROGS/bin_setup/old_newf;

# The previous file creates "old" hycom input file   convert them to .ab - type
$MSCPROGS/bin_setup/forfun_nersc-2.2 month old


# For now the "old" nersc forcing is region-dependent due to the region-dependent 
# levitus forcing
mkdir $CDIR/${E}
for i in forcing.*.[ab] ; do
   new=$(echo $i | sed "s/^forcing\.//")
   mv $i $CDIR/$E/$new
done

echo
echo " Forcing data in $CDIR/$E"
