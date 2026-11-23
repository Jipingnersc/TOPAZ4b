#!/bin/bash
if [ $# -ne 3 -a $# -ne 4 ] ; then  
   echo "Usage : $(basename $0) rungen byear bday "; 
   echo "Usage : $(basename $0) rungen byear bday -skipargoprof"; 
   exit; 
fi

# Test argument..
if [ "$4" == "-skipargoprof" ] ; then
   doargoprof=0
else
   doargoprof=1
fi

# Run id
rungen=$1

# Retrieve bulletin year and day from input args
byear=$2
bjuld=$3
bjuld2=$(echo 00$bjuld | tail -c4)


# Clear ArgoCmp directory
if [ $doargoprof -eq 1 ] ; then
   [ -d ArgoCmp ] && rm ./Argocmp/*
fi

# file list
files=$(ls ${rungen}DAILY_${byear}_${bjuld2}_*.a)
#echo $files

# Go through file list to get max and min dates
mindate=30000000
maxdate=0
for k in $files ; do 
   iday=$(echo $k | sed "s/${rungen}DAILY_${byear}_${bjuld2}_[0-9]\+_//" | sed "s/\..*//")
   iyear=$(echo $k | sed "s/${rungen}DAILY_${byear}_${bjuld2}_//" | sed "s/\_[0-9]\+.*//")
   echo $iyear $iday

   # convert to dates
   datestr=$(jultodate $iday $iyear 1 1)

   [ $datestr -gt $maxdate ] && let maxdate=$datestr
   [ $datestr -lt $mindate ] && let mindate=$datestr
done
echo "Max and min date for bulletin:  $maxdate $mindate"

# Floating point year
newyear=$(datetojul $byear 12 31 $byear 1 1)
fyear=$(echo $byear $bjuld $newyear | perl -lane 'print $F[0] + $F[1]/$F[2]')
#echo $fyear
#exit


# Extract profiles - finds argo profiles which match file dates, 
# interpolates model to obs points. Result is dumped in directory 
# ArgoCmp
if [ $doargoprof -eq 1 ] ; then
   ./argo_extract_profiles.sh $files
fi

# Depth intervals to process
declare -a dlevels
dlevels=(0 100 300 500 1000 3000)
numlevels=${#dlevels[@]}

#Regions to process
declare -a regions
regions=(NAT NORDIC LABSEA)
numregions=${#regions[@]}

#Start processing sub - regions
regcnt=1
while [ $regcnt -le $numregions ] ; do
   levcnt=1
   while [ $levcnt -lt $numlevels ] ; do

      D1=${dlevels[$levcnt-1]}
      D2=${dlevels[$levcnt]}
      REG=${regions[$regcnt-1]}

      echo $D1 $D2 $REG

      # The remaining routines use files in "ArgoCmp"
      # Process "NAT" region - 0 - 100 m - mindate to maxdate
      if [ $levcnt -eq 1 ] ; then
         ./argo_region.sh $REG $D1 $D2 $mindate $maxdate -plotdata
      else
         ./argo_region.sh $REG $D1 $D2 $mindate $maxdate 
      fi
      POST=$REG.$D1.$D2
      echo $fyear $(cat argoreduce.out) >> argoreduce.${POST}.timeseries

      #exit

      let levcnt=$levcnt+1
   done
   let regcnt=$regcnt+1
done
