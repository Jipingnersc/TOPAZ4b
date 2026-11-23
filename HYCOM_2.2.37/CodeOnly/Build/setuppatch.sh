#!/bin/bash
# Script to set up dimensions for compiling hycom
#
# Retrieves
# itdm - Total grid size from patch.input
# jtdm - Total grid size from patch.input
# mpe  - number of tiles in i-dir from patch.input
# npe  - number of tiles in j-dir from patch.input
# idm  - Tile  grid size from patch.input
# jdm  - Tile  grid size from patch.input
# kdm  - vertical layers from blkdat.input

# Must input number of CPUs
if [ $# -eq 1 ] ; then
   NMPI=$1
else
   echo "Usage: $(basename $0) NMPI"
   echo "where NMPI is the number of MPI patches to use"
   echo "(Set this to 0 for OpenMP-only)"
   exit 1
fi
echo "Number of CPUs set to $NMPI" 

# MUST be run in build directory
PWD=$(pwd)
X=$(echo $PWD | sed "s/.*\///" | sed "s/.*_X//")
source ../REGION.src
source ../expt_$X/EXPT.src


# Look for blkdat.input file in ../expt_$X/blkdat.input
# Where X is experiment number 
blkdatfile="../expt_${X}/blkdat.input"
if [ -f $blkdatfile ] ; then
   echo "found blkdat-file $blkdatfile"
else
   echo "Could not find blkdat-file $blkdatfile"
   exit
fi

# Some requirements
# Info from blkdat.input
echo
iversn=$(cat  $blkdatfile | grep "'iversn'" | sed "s/^[ ]*//" | tr  \t " " |  cut -f1 -d" ")
iexpt=$(cat  $blkdatfile | grep "'iexpt '" | sed "s/^[ ]*//" | tr  \t " " |  cut -f1 -d" ")
itdm2=$(cat  $blkdatfile | grep "'idm   '" | sed "s/^[ ]*//" | tr  \t " " |  cut -f1 -d" ")
jtdm2=$(cat  $blkdatfile | grep "'jdm   '" | sed "s/^[ ]*//" | tr  \t " " |  cut -f1 -d" ")
kdm=$(cat  $blkdatfile | grep "'kdm   '" | sed "s/^[ ]*//" | tr  \t " " |  cut -f1 -d" ")
echo "blkdat.input says:"
echo "itdm   = $itdm2"
echo "jtdm   = $jtdm2"
echo "kdm    = $kdm"
echo "iexpt  = $iexpt"
echo "iversn = $iversn"

[ $iversn -ne 22 ] && echo "warning: $blkdatfile has iversn = $iversn. Change to version 22 for hycom 2.2"
if [ $iexpt  -ne $E ] ; then
   echo "ERROR: $blkdatfile has iexpt = $iexpt - change this to the current experiment ($E)"
   echo "ERROR: no dimensions set up"
   exit 1
fi

if [ $NMPI -gt 1 ] ; then
   # Look for partition file in ../topo/partit/depth_${R}_${T}.${NMPI}
   # Where R is top level dir, T is topography version (depends on X) and 
   # NMPI is number of CPUs
   partitfile=../topo/partit/depth_${R}_${T}.`echo 000$NMPI | tail -5c`
   #echo $partitfile
   if [ -f $partitfile ] ; then
      echo "found partit-file $partitfile"
   else
      echo "ERROR:Could not find partit-file $partitfile - "
      echo "ERROR:you need to set it up if you want to use NMPI=$NMPI"
      exit
   fi

   # Retrieve dim sizes from patch.input
   input=$(head -n2 $partitfile | tail -n1 | sed "s/^[ ]*//" | sed "s/[ ]+/ /" )
   #echo $input
   npes=$(echo $input | cut -f1 -d" ")
   mpe=$(echo $input  | cut -f2 -d" ")
   npe=$(echo $input  | cut -f3 -d" ")
   itdm=$(echo $input | cut -f4 -d" ")
   jtdm=$(echo $input | cut -f5 -d" ")
   idm=$(echo $input  | cut -f6 -d" ")
   jdm=$(echo $input  | cut -f7 -d" ")

   # Info from patch.input
   echo
   echo "patch.input says:"
   echo "NCPUS = $npes"
   echo "mpe   = $mpe"
   echo "npe   = $npe"
   echo "itdm  = $itdm"
   echo "jtdm  = $jtdm"
   echo "idm   = $idm"
   echo "jdm   = $jdm"
else
   echo "NMPI set to $NMPI - patch.input not needed"
   npes=1
   mpe=1
   npe=1
   itdm=$itdm2
   jtdm=$jtdm2
fi


# Consistency check
[ $itdm -ne $itdm2 ] && { echo "blkdat.input and patch.input disagree on itdm" ; exit 3; }
[ $jtdm -ne $jtdm2 ] && { echo "blkdat.input and patch.input disagree on jtdm" ; exit 4; }

filnavn="dimensions_nersc.h"
echo "CKAL -- script-generated dimensions from blkdat.input and patch.input
      integer, parameter :: itdm=$itdm
      integer, parameter :: jtdm=$jtdm
      integer, parameter :: kdm=$kdm
      integer, parameter :: idm=$idm
      integer, parameter :: jdm=$jdm
      integer, parameter :: iqr=$mpe
      integer, parameter :: jqr=$npe " > $filnavn


echo
echo "$filnavn is created, set info in flags and you can compile"



