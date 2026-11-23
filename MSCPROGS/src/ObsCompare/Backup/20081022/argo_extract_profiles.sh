#!/bin/bash
# Routine extracts profiles for input files. it Finds correct dates by checking
# dates of the files you specify. This involves checking locally (see DATADIR below), 
# as well as looking in the FTP web site.
#
# Finally the vertical/horizontal interpolation is done, using the argocmp routine.
# 

# Check number of args - must be greater than or equal to 1
if [ $# -lt 1 ] ; then 
   echo "Usage:"
   echo "    $(basename $0) files " 
   exit
fi

ncdump=""
if [ "$1" == "-ncdump" ] ;then
   ncdump="-ncdump"
fi

# This assumes routines are in the same location as this routine
basedir=$(dirname $0)


files=$@



# This location contins data cache - when not obtained from ifremer
DATADIR=/home/fimm/nersc/knutali/Data/Argo 

# Where argo data can be found - FTPDIR contains month/year as well (see below)
FTPARGO=ftp.ifremer.fr
FTPDIRARGO=ifremer/argo/geo/atlantic_ocean/


# For each file check that it is present in my Data dir
for i in $files; do

   # Check that this is a daily file. To be changed later on
   if echo $i | grep "DAILY" ; then

   rungen=$(echo $i | cut -c1-3)

   # Get date info from file name FORDAILY_2006_108_2006_108.b
   #thisrefyear=$(echo $i | sed "s/${rungen}DAILY_[0-9]+_[0-9]+_//" | \
   #              sed  "s/_.*//")
   thisrefyear=$(echo $i | sed "s/${rungen}DAILY_[0-9]\+_[0-9]\+_//" | sed "s/_.*//")
   thisjday=$(echo $i |    sed "s/${rungen}DAILY_[0-9]\+_[0-9]\+_[0-9]\+_//" |  sed  "s/\..*//")

   #echo $rungen  $thisrefyear $thisjday
   #exit

   # Get actual date 
   act_date=$(jultodate $thisjday $thisrefyear 1 1)
   act_year=$(echo $act_date | cut -c1-4)
   act_month=$(echo $act_date | cut -c5-6)

   # This is the argo file name
   argofile=${act_date}_prof.nc
   argofilegz=${act_date}_prof.nc.gz

   echo $thisrefyear $thisjday $act_date $act_year $act_month
   #exit

   # Is the corresponding  file present in my data dir ?
   tmpgz=${DATADIR}/${act_year}/${act_month}/${argofilegz}
   tmpdir=${DATADIR}/${act_year}/${act_month}/
   if [ -f $argofile  ] ; then
      # Copy to data to use later
      [ ! -d $tmpdir ] && mkdir -p $tmpdir
      [ ! -f $tmpgz  ] && gzip < $argofile > $tmpgz
      echo "$argofile already present"
   elif [ -f $tmpgz ] ; then
      gunzip < ${tmpgz} > $argofile
      echo "$argofile Retrieved from local backup"
   elif [ ! -f $tmpgz ] ; then

      # Try to retrieve from ifremer
ftp << EOF
open $FTPARGO
cd $FTPDIRARGO/${act_year}/${act_month}/
bin
hash
prompt
get $argofile
bye
EOF

      echo "$argofile:"
      if [ ! -f $argofile  ] ; then
         echo "Retrieve failed !"
      else
         echo "Retrieve successful"
      fi
   fi


   # The file SHOULD be here now
   if [ -f $argofile  ] ; then
      $basedir/argocmp $argofile $i $ncdump
   fi

   fi # Check for daily
done



# At this stage, we have a catalog Argocmp, with profiles from model and data.
