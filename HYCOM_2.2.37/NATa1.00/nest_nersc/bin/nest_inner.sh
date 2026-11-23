#!/bin/bash
#


#set -x
pget=cp
pput=cp
# Experiment also needs topography (through experiment number)
if [ $# -ne 1 ] ; then
   echo "This routine will set up nesting for the inner model. "
   echo "You need to specify the experiment number of the inner "
   echo "model (this model)"
   echo 
   echo "Routine sets up the nesting relaxation masks (rmu.[ab]),"
   echo "and a ports file needed to do nesting. Ports are placed"
   echo "whereever ocean is present in the 2nd grid cell from eastern,"
   echo "western, northern and southern edges"
   echo 
   echo "NB: nestbat.sh should be rune before this routine"
   echo
   echo "Usage: "
   echo "$(basename $0) experiment_number"
   exit
fi
export X=$1

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd -
echo "Sourcing Inner model"
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }
echo X=$X E=$E T=$T R=$R

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
# --- Create nesting fields for outer model (should be compatible with standard hycom)
#
#
# --- S is scratch directory,
# --- D is permanent directory,
#
D=$BASEDIR/nest_nersc/$E/inner/
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
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a || \
   { echo "Could not get inner regional.depth.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b || \
   { echo "Could not get inner regional.depth.b file " ; exit 1 ; }
${pget} ${BASEDIR}/expt_${X}/blkdat.input blkdat.input || \
   { echo "Could not get blkdat.input file " ; exit 1 ; }

#Create nesting relaxation mask
echo "Setting up relaxaion mask:"
echo "--------------------------"
${MSCPROGS}/src/Nesting-2.2/nestini
echo
echo
echo

#Create nesting ports (ports along boundary where applicable)
echo "Setting up boundary ports:"
echo "--------------------------"
${MSCPROGS}/src/Nesting-2.2/nestport
echo
echo
echo


# Move to experiment directory - prepend with id and topo version of inner model
for i in rmu.[ab]  ports.input.nest nestrel.nc ; do
   mv $i $D/
done

echo
echo
echo
echo "All looks ok - ports file and relaxation file (rmu.[ab]) in $D"
echo "Diagnostic file nestrel.nc also in $D"
