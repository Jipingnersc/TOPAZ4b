#!/bin/bash
#
set -x
source ../../EXPT.src
N=$T
#
# --- force/plot directory softlinks for topography and source.
#
#

touch   regional.grid.a regional.depth.a
/bin/rm *.[ab]
#
/bin/ln -s ../../topo/regional.grid.a   regional.grid.a
/bin/ln -s ../../topo/regional.grid.b   regional.grid.b
#
/bin/ln -s ../../topo/depth_${R}_${N}.a regional.depth.a
/bin/ln -s ../../topo/depth_${R}_${N}.b regional.depth.b
#
if [ -s  ../../topo/landsea_${R}.a ] ; then
  /bin/ln -s ../../topo/landsea_${R}.a regional.mask.a
fi
#
#if (! -e src) then
#  /bin/ln -s ../../../ALL/plot/src_2.1.00 src
#endif

