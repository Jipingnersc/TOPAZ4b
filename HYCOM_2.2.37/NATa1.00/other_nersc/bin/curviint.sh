#!/bin/bash
#set -x
pget=cp
pput=cp

# Experiment needs topography
if [ $# -ne 4 ] ; then
   echo "This script will set up and run the curviint routine. This"
   echo "routine interpolates from a global model onto a local grid"
   echo "The vertical structure (layers) of the outer model are assumed"
   echo "to be the same as in the inner model"
   echo
   echo "Input to this routine is:"
   echo "1)Version number of local/nested model"
   echo "2)Version number of outer/gobal model"
   echo "3)Path to the top-level directory of global model"
   echo "4)full path to a hycom restart file of the global model"
   echo
   echo "Example:"
   echo "  $(basename $0) 01.1 01.0 /work/knutali/h22_tst/TP3a0.12 /work/knutali/h22_tst/TP3a0.12/expt_01.0/data/TP3restart0003_020_00.a"
   echo "This will run curviint on TP3restart0003_020_00.a - which is from experiment 01.0 of the global model. The top-level"
   echo "region directory of the global model is /work/knutali/h22_tst/TP3a0.12" 
   exit
fi

ODIR=$3 #"Outer model top-level directory
RPATH=$4 #Full path to restart file

# Source outer model first, and keep outer model setup in new variable names
echo "Sourcing Outer model"
OX=$2
source ${ODIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${ODIR}/expt_$OX/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$OX/EXPT.src" ; exit 1 ; }
OE=$E
OT=$T
OR=$R


# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd -
echo "Sourcing Inner model"
X=$1
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }

# Check that pointer to MSCPROGS is set
echo $MSCPROGS
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



# Data and scratch paths
D=$BASEDIR/curviint/$E
S=$D/SCRATCH
mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}


# Set up scratch dir with input data for curviint - this is for the inner model
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b     || { echo "Could not get regional.grid.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a     || { echo "Could not get regional.grid.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a  || { echo "Could not get regional.depth.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b  || { echo "Could not get regional.depth.b file " ; exit 1 ; }
${pget} ${BASEDIR}/expt_${X}/blkdat.input blkdat.input      || { echo "Could not get file blkdat.input " ; exit 1 ; }


# Set up scratch dir with input data from GLOBAL model 
cp $ODIR/topo/depth_${OR}_${OT}.a global.depth.a || { echo "Could not get global file " ; exit 1 ; } 
cp $ODIR/topo/depth_${OR}_${OT}.b global.depth.b || { echo "Could not get global file " ; exit 1 ; } 
cp $ODIR/topo/regional.grid.a  global.grid.a || { echo "Could not get global file " ; exit 1 ; } 
cp $ODIR/topo/regional.grid.b  global.grid.b || { echo "Could not get global file " ; exit 1 ; } 
cp $ODIR/topo/grid.info grid.info || { echo "Could not get global file " ; exit 1 ; } 

# Set up scratch dir with restart file from GLOBAL model 
[ -f "$RPATH" ] ||  { echo "HYCOM restart file $RPATH does not exist..";  exit 1 ; }
RDIR=$(dirname $RPATH)                #directory with restart file
RFILE=$(basename $RPATH)              #restart file
RBASE=$(echo $RFILE | sed "s/\.a$//") #restart file without appendix
cp $RDIR/${RBASE}*  ./


$MSCPROGS/bin_setup/curviint-2.2  ${RBASE}.a ${RBASE}ICE.uf
#$MSCPROGS/bin_setup/curviint-2.2  <<EOF
#
#${RBASE}.a
#${RBASE}ICE.uf
#EOF

# Copy to experment dir
TARGET=$BASEDIR/curviint/$E/
mkdir -p $TARGET
for i in curviint.*; do
   nname=$(echo $i | sed "s/curviint\.//")
   cp $i $TARGET/$nname
done


