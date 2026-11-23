#!/bin/bash
# Set up 2d tiling for region

# Check that environment is set
if [ -z ${BASEDIR} ] ; then
   echo "Environment not set "
   exit
fi


# Input args
[ $# -ne 2 ] && [ $# -ne 3 ] \
             && echo "Usage:" \
             && echo "   $(basename $0) [1st tile dimension] [2nd tile dimension]" -noplot \
             && echo "  tip  - two negative tile dimensions gives a uniform grid " \
             && exit
iqr=$1
jqr=$2
echo "iqr is  $iqr  , jqr is $jqr"
export nm=`echo $iqr $jqr | awk '{printf("%4.4d\n",$1*$2)}'`

source $BASEDIR/EXPT.src
PROGDIR="$BASEDIR/../ALL/"
export PARTIT=$PROGDIR/topo/src/partit
export TOPO_PPM=$PROGDIR/topo/src/topo_ppm

export FOR051=$BASEDIR/topo/depth_${R}_${T}.b
export FOR051A=$BASEDIR/topo/depth_${R}_${T}.a
export FOR021=fort.21

touch regional.grid.a regional.grid.b
rm    regional.grid.a regional.grid.b
cp  $BASEDIR/topo/regional.grid.a .
cp  $BASEDIR/topo/regional.grid.b .



# Calculate "load-balancing" patches
echo "$iqr $jqr 0.75" | $PARTIT
cp fort.21 depth_${R}_${T}.$nm

# Create graphics
$TOPO_PPM
mv fort.31 depth_${R}_${T}.$nm.ppm
rm fort.*

# Show it ...
if [ "$3" != "-noplot" ] ; then
   display  depth_${R}_${T}.$nm.ppm
fi

