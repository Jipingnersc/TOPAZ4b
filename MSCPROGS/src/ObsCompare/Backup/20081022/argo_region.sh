#!/bin/bash
# This routine looks in the catalogue ArgoCmp previously created profiles
# by the routine argocmp. It will look in that catalogue for files matching
# region, depth and time restrictions specified. Region specification templates
# are present in the $REGIONDIR catalogue. depth and time speifications 
# are given as arguments to this routine
#
#NB - time format is yyyymmdd

usage="Usage   $(basename $0) [opional args]    region depth1 depth2 time1 time2 time3 .... \n
Optional arguments: \n
   -p  : plots points in data region "

while getopts "p" options; do
  case $options in
     p ) plot=1 ;;
     * ) echo -e $usage
     exit 1;;
  esac
done
shift $(($OPTIND - 1))

# Check number of args - must be greater equal to 5
if [ $# -lt 5  ] ; then 
   echo -e $usage
   exit 1
fi

# This location contain info on regions. If
# arguments are supplied, the routine will look in this 
# catalogue for argoreduce.in.${REGION} files. If no arguments a
# are given, it will look for argoreduce.in in the present directory
REGIONDIR=/home/fimm/nersc/knutali/Data/Argo/Regions
REGION=$1
D1=$2
D2=$3
shift 3

[ -f argoreduce.timeseries ] && rm argoreduce.timeseries 

while [ $# -ge 2  ] ; do

   T1=$1
   T2=$2
   [ -f ${REGIONDIR}/argoreduce.in.${REGION} ] &&  cp ${REGIONDIR}/argoreduce.in.${REGION} .
   [ ! -f ${REGIONDIR}/argoreduce.in.${REGION} ] && \
      { echo "No region file argoreduce.in.${REGION}  found ! " ; exit ;}

   # Produce argoreduce.in
   cat argoreduce.in.${REGION} | sed "s/D1/$D1/" | sed "s/D2/$D2/"  \
                               | sed "s/T1/$T1/" | sed "s/T2/$T2/" > argoreduce.in

   #Name to append to output files
   POST=$REGION.$T1.$T2

   # Reduction routine which calculate statistics over regions and
   # depth intervals. Reads argoreduce.in
   [ -f argoreduce.out ] && rm argoreduce.out
   [ -f argoreduce.locations ] && rm argoreduce.locations
   ./argoreduce

   # Produce a plot of the locations used in the calculation
   if [ $plot==1 -a -s argoreduce.locations ] ; then
matlab -nojvm -nodisplay << EOF
plot_argoreduce_locations();
exit;
EOF
      cp argoreduce.locations.eps argoreduce.locations.${POST}.eps
      cp argoreduce.locations.png argoreduce.locations.${POST}.png
   elif [  ! -s argoreduce.locations ] ; then
      echo "ERROR: No argo locations match criteria ... "
   fi

   # Move the argoreduce.out to a more descriptive file
   #???
   cat argoreduce.out >> argoreduce.timeseries 

   shift
done
