#!/bin/bash
#######################################################################
# Start post-processing

#[ -r ./data ] && mv data data_KEEP
#[ -r ./data_results ] && ln -sf data_results data

#Check environment
cd $(dirname 0)/../
export BASEDIR=/cluster/work/users/xiejp/TP4b0.12
echo BASEDIR=$BASEDIR
if [ -z ${BASEDIR} ] ; then
   tellerror "BASEDIR Environment not set "
   exit
fi

# Initialize environment (sets Scratch dir ($S), Data dir $D ++ ) - must be in
# experiment dir for this to work
cd ${BASEDIR}/expt_01.1

source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src "; exit 1; }


# Print error if we cannot descend scratch dir $S
cd       $S || { cd $P ; echo "BADRUN" > log/hycom.stop ;  exit 1 ;}

touch   PIPE_DEBUG
/bin/rm PIPE_DEBUG


#
# Move restart and daily files to data location $D
# first test if its a link, if so, then remove it
for i in ???restart*
do
  if [ -L "${i}" ]
    then
       rm $i
  fi
done

#mv ???restart* $D
mv TP4restart* TP4DAILY* TP4AVE*  $D
#mv TP4restart* TP4DAILY* TP4AVE* ???archv.* ???GP[0-9][0-9][0-9][0-9]_[0-9][0-9]/  $D
#Tim: change permissions of these files so nersc group can read them:
cd $D
chmod 640 ???restart* ???DAILY* ???AVE* ???archv.* ???GP[0-9][0-9][0-9][0-9]_[0-9][0-9]
cd $S
chmod 640 wave*

#
# Copy some files useful for analysis
#cp regional.grid.* regional.depth.* blkdat.input archv.* $D

# Move the nesting output files to the data directory.
#mv nest_out_* $D



#
# --- HYCOM error stop is implied by the absence of a normal stop.
#
touch summary_out
tail -1 summary_out
tail -1 summary_out | grep -c "^normal stop"
if  [ `tail -1 summary_out | grep -c "^normal stop"` == 0 ] ; then
  cd $P
  echo "BADRUN"  > log/hycom.stop
else
  cd $P
  echo "GOODRUN" > log/hycom.stop
fi
