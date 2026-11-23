#!/bin/bash
pget=cp
pput=cp
# Experiment also needs topography (through experiment number)
if [ $# -ne 5 ] ; then
   echo "This routine will set up nesting bathymetry for the inner model. "
   echo "You need to specify a topography version of the inner "
   echo "model (this model). This topography will be used as a STARTING"
   echo "POINT(!) for the new, modified topography. You will also need "
   echo "to specify the NEW topography version of the 'nestbatted' topography"
   echo "In addition you need to specify the path to the outer model top-"
   echo "level directory along with the experiment number of the outer model"
   echo "you wish to use for nesting"
   echo
   echo "Example:"
   echo "  $(basename $0) 01.0 00 01 /work/knutali/h22_tst/TP3a0.12/ 01.1"
   echo
   echo "where 01.0 is the experiment number of the inner model and you wish "
   echo "to use as a starting point, 00 is the correponding bathy"
   echo ", and 01 is the 'nestbatted' topo version "
   echo "01.1 is the experiment number of the outer model. The path"
   echo "to the top-level directory of outer model is in this example "
   echo "/work/knutali/h22_tst/TP3a0.12/"
   exit
fi

X=$1
TOPO_ORIG=$2
TOPO_FINAL=$3
PATHOUTER=$4
OX=$5
 
# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd -
echo "Sourcing Inner model"
echo $E $X

# Source inner model setup
# Run with noexperiment because we just want to modify the topography at this point
source ${BASEDIR}/expt_${X}/EXPT.src $X || { echo "Could not source ${BASEDIR}/EXPT.src" ; exit 1 ; }
echo X=$X E=$E T=$T R=$R
IR=$R

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

# Now source inner model setup
echo
X=$OX
echo "Sourcing Outer model"
source ${PATHOUTER}/REGION.src $OX || { echo "Could not source ${PATHOUTER}/REGION.src" ; exit 1 ; }
source ${PATHOUTER}/expt_${OX}/EXPT.src $OX || { echo "Could not source ${PATHOUTER}/EXPT.src" ; exit 1 ; }
echo X=$X E=$E T=$T R=$R

# Rename these variables, cause we need to get the same from the outer model
OX=$X
OE=$E
OT=$T
OR=$R

#Unset these variables to avoid unwanted effects
unset  -v X E T R


#
#
# --- Create nesting fields for outer model (should be compatible with standard hycom)
#
#
# --- S is scratch directory,
# --- D is permanent directory,
#
D=$BASEDIR/topo/
S=$D/SCRATCH
mkdir -p $S
mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

# Get regional files and grid.info from inner model
touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b grid.info
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b grid.info
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b    || \
   { echo "Could not get inner regional.grid.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a    || \
   { echo "Could not get inner regional.grid.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${IR}_${TOPO_ORIG}.a regional.depth.a || \
   { echo "Could not get inner regional.depth.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${IR}_${TOPO_ORIG}.b regional.depth.b || \
   { echo "Could not get inner regional.depth.b file " ; exit 1 ; }

# Get regional files from outer model - note rename !
${pget} ${PATHOUTER}/topo/regional.grid.b global.grid.b    || \
   { echo "Could not get outer regional.grid.b file " ; exit 1 ; }
${pget} ${PATHOUTER}/topo/regional.grid.a global.grid.a    || \
   { echo "Could not get outer regional.grid.a file " ; exit 1 ; }
${pget} ${PATHOUTER}/topo/depth_${OR}_${OT}.a global.depth.a || \
   { echo "Could not get outer regional.depth.a file " ; exit 1 ; }
${pget} ${PATHOUTER}/topo/depth_${OR}_${OT}.b global.depth.b || \
   { echo "Could not get outer regional.depth.b file " ; exit 1 ; }
${pget} ${PATHOUTER}/topo/grid.info grid.info || \
   { echo "Could not get outer grid.info file " ; exit 1 ; }

# Run nestbat 
${MSCPROGS}/src/Nestbat/nestbat-2.2 

# Check if old topo version exists
if [ -f ${BASEDIR}/topo/depth_${IR}_${TOPO_FINAL}.a  -o  -f ${BASEDIR}/topo/depth_${IR}_${TOPO_FINAL}.b ] ; then
   echo
   echo
   echo "ERROR: topography file(s) depth_${IR}_${TOPO_FINAL}.[ab] already exist"
   echo "       Cowardly refusing to overwrite existing files....."
else
   echo
   echo
   cp -i nestbat.depth.a ${BASEDIR}/topo/depth_${IR}_${TOPO_FINAL}.a  
   cp -i nestbat.depth.b ${BASEDIR}/topo/depth_${IR}_${TOPO_FINAL}.b  
   echo "All looks fine, copying depth_${IR}_${TOPO_FINAL}.[ab] to $D"
   echo "Diagnostic files can be found in $S"
   echo "You are now ready to set up a nesting experiment.."
fi
