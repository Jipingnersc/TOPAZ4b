#!/bin/bash
if [ $# -ne 2 ] ; then
   echo
   echo "$(basename $0) - a routine to set up and initialize a region, "
   echo "based on existing regional.* files (++).  Input to this routine must"
   echo "be topography version T (e.g 01) and location of existing"
   echo "topo/grid files."
   echo " Example:"
   echo "   $(basename $0) 01 /work/$USER/TOPAZ3/"
   exit 1
fi
TIN=$1
DIR=$2
[ ! -d $DIR ] && { echo "input dir $DIR does not exist" ; exit 1 ; }



# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd -
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
T=$TIN

fbase=depth_${R}_${T}
echo
mkdir -p $BASEDIR/topo/SCRATCH/
rm $BASEDIR/topo/SCRATCH/*


# First check if topo dir contains regional.grid.[ab]. Copy the files in there if not present
if [ ! -f $BASEDIR/topo/regional.grid.a -o ! -f $BASEDIR/topo/regional.grid.b  ] ; then
   echo "Copying regional.grid files into $BASEDIR/topo"
   cp $DIR/regional.grid* $BASEDIR/topo/
fi
if [ ! -f $BASEDIR/topo/grid.info ] ; then
   echo "Copying grid.info file into $BASEDIR/topo"
   cp $DIR/grid.info $BASEDIR/topo/
fi
if [ -f $DIR/latlon.dat ] ; then
   if [ ! -f $BASEDIR/topo/latlon.dat ] ; then
      echo "Copying latlon.dat file into $BASEDIR/topo"
      cp $DIR/latlon.dat $BASEDIR/topo/
   else
      echo "File $BASEDIR/topo/latlon.dat already present"
      echo "Cowardly refusing to overwrite it!"
   fi
fi
if [ -f $DIR/newpos.uf ] ; then
   if [ ! -f $BASEDIR/topo/newpos.uf ] ; then
      echo "Copying newpos.uf file into $BASEDIR/topo"
      cp $DIR/newpos.uf $BASEDIR/topo/
   else
      echo "File $BASEDIR/topo/newpos.uf already present"
      echo "Cowardly refusing to overwrite it!"
   fi
fi


# copy files to scratch
cp $DIR/regional.grid.a $BASEDIR/topo/SCRATCH/ || \
   { echo "could not copy regional.grid.a files from region topo dir - " ; \
      echo "Make sure they are present!"; exit 1; }
cp $DIR/regional.grid.b $BASEDIR/topo/SCRATCH/ || \
   { echo "could not copy regional.grid.b files from region topo dir - " ; \
      echo "Make sure they are present!"; exit 1; }
if [ -f $DIR/regional.depth.a ] ; then
   cp $DIR/regional.depth.* $BASEDIR/topo/SCRATCH/
fi
if [ -f $DIR/depths*.uf ] ; then
   cp $DIR/depths*.uf       $BASEDIR/topo/SCRATCH/
fi
if [ -f $DIR/ndepths*.uf ] ; then
   cp $DIR/ndepths*.uf      $BASEDIR/topo/SCRATCH/ 
fi
cd $BASEDIR/topo/SCRATCH/ || { echo "Cant descend $BASEDIR/topo/SCRATCH/"; exit 1 ;}


# Get idm, jdm from regional.grid.[ab]
rgidm=$(head -n1 regional.grid.b | sed "s/'idm.*//" | tr -d [:space:] )
rgjdm=$(head -n2  regional.grid.b | tail -n1 | sed "s/'jdm.*//" |tr -d [:space:] ) 
rgidm=$(head -n1 regional.grid.b | sed "s/^[ ]*//" | sed "s/[ \t].*//")
rgjdm=$(head -n2 regional.grid.b | tail -n1 | sed "s/^[ ]*//" | sed "s/[ \t].*//")
echo
echo
echo
echo "Grid size from regional.grid.b $rgidm $rgjdm "

# first attempt is to use existing regional.depth files
#if [ 0 -eq 1 ] ; then
if [ -f regional.depth.a -a -f regional.depth.b ] ; then
   # copy addhuge utility to scratch - run program
   echo "New topography files will be based on  regional.depth.[ab]"
   cp $BASEDIR/topo/bin/addhuge.py $BASEDIR/topo/SCRATCH/
   ./addhuge.py regional.depth $fbase || { echo "addhuge failed .." ; exit 1; }

# second attempt is to use existing ndepths file 
elif ls ndepths*.uf &> /dev/null ; then
   fdepths=$(ls ndepths*.uf)

   idm=$(echo $fdepths | sed "s/ndepths//" | sed "s/x.*//")
   jdm=$(echo $fdepths | sed "s/ndepths.*x//" | sed "s/\.uf.*//")

   if [ $idm -ne $rgidm  -o $jdm -ne $rgjdm ] ; then
      echo "Mismatch betweem depths*.uf dimension and regional.grid.b dimension"
      exit 1
   fi

   echo "New topography files will be based on  $fdepths  with dimensions $idm $jdm"
   cp $BASEDIR/topo/bin/depthtoregional.py $BASEDIR/topo/SCRATCH/
   ./depthtoregional.py $fdepths $fbase $idm $jdm || \
      { echo "depthtoregional failed .." ; exit 1; }

   echo
   echo
   echo "NB! File is based on ndepths*.uf file - "
   echo "NB! make sure it shouldnt be based on depths*.uf!"

# last attempt is to use existing depths file
elif ls depths*.uf &> /dev/null ; then
   fdepths=$(ls depths*.uf)
   idm=$(echo $fdepths | sed "s/depths//" | sed "s/x.*//")
   jdm=$(echo $fdepths | sed "s/depths.*x//" | sed "s/\.uf.*//")

   if [ $idm -ne $rgidm  -o $jdm -ne $rgjdm ] ; then
      echo "Mismatch betweem depths*.uf dimension and regional.grid.b dimension"
      exit 1
   fi

   echo "New topography files will be based on  $fdepths  with dimensions $idm $jdm"
   cp $BASEDIR/topo/bin/depthtoregional.py $BASEDIR/topo/SCRATCH/
   ./depthtoregional.py $fdepths $fbase $idm $jdm || \
      { echo "depthtoregional failed .." ; exit 1; }

   echo
   echo
   echo "NB! File is based on depths*.uf file - "
   echo "NB! make sure it shouldnt be based on ndepths*.uf!"

else
   echo "Could not find any topography files (tried regional.depth.[ab],"
   echo "depths*.uf, ndepths*.uf"
   exit 1
fi


echo
[ -f $BASEDIR/topo/$fbase.a ] && \
   { echo "$fbase.a already exists in topo dir - I will not copy.." ; \
     echo "(copy by hand if you feel its right  ....)" ; exit 1; }
[ -f $BASEDIR/topo/$fbase.b ] && \
   { echo "$fbase.b already exists in topo dir - I will not copy.." ; \
     echo "(copy by hand if you feel its right  ....)" ; exit 1; }

# Safe to copy now
echo "Copying new files to  $BASEDIR/topo"
cp $fbase.a $fbase.b $BASEDIR/topo/




