#!/bin/bash
if [ $# -ne 1 ] ; then
   echo "Input to this routine must be TOPO version T (e.g 01)"
   exit 1
fi
T=$1



# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd -
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }

echo
echo "Topo names for topography version $T are:"
echo depth_${R}_${T}.a  depth_${R}_${T}.b

