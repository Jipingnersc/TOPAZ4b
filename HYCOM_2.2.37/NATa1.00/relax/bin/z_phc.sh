#!/bin/bash

# KAL - get X from input
if [ $# -ne 1 ] ; then
   echo "This script will set up files with the PHC climatology. The "
   echo "climatology will be given on z-levels, and interpolated in "
   echo "the horizontal - also onto land points. The resulting files"
   echo "are needed for setting up the relaxation files to be used in"
   echo "HYCOM."
   echo
   echo "The only input is ksigma (0 or 2)"
   exit
fi
export KSIGMA=$1


# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd -
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }

# Check that pointer to HYCOM_ALL is set
if [ -z ${HYCOM_ALL} ] ; then
   echo "Environment not set "
   exit
else
   if [ ! -d ${HYCOM_ALL} ] ; then
      echo "HYCOM_ALL not properly set up"
      echo "HYCOM_ALL not a directory at ${HYCOM_ALL}"
      exit
   fi
fi

#set -x
pget=cp
pput=cp

#
# --- Interpolate Levitus climtology to a HYCOM region.
# --- KSIGMA=0 for _sig0_ fields
#
#
# --- P is primary path,
# --- S is scratch   directory,
# --- D is permanent directory,
# --- PHC_PATH is PHC directory
#
D=$BASEDIR/relax/phc/
S=$D/SCRATCH


[  "${PHC_PATH}" == ""  ] && { echo "PHC_PATH is not set " ; exit 1; }
[ ! -d ${PHC_PATH} ] && { echo "Cant find PHC_PATH directory at ${PHC_PATH}" ; exit 1; }

mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}


touch regional.grid.a regional.grid.b
rm    regional.grid.a regional.grid.b
touch   z_levitus
/bin/rm z_levitus
#
# --- 12 months
#
for MM  in  01 02 03 04 05 06 07 08 09 10 11 12 ; do
   #
   # --- Input.
   #
   # --- Download Levitus from:
   # --- ftp://obelix.rsmas.miami.edu/awall/hycom/levitus_for_hycom.tar.gz
   #
   touch      fort.71 fort.73
   /bin/rm -f fort.71 fort.73
   ${pget} ${PHC_PATH}/r_m${MM}.d fort.71 || { echo "Cant copy ${PHC_PATH}/r_m${MM}.d "; exit 1 ; }
   ${pget} ${PHC_PATH}/s_m${MM}.d fort.73 || { echo "Cant copy ${PHC_PATH}/r_m${MM}.d "; exit 1 ; }
   
   if [ ! -s regional.grid.b ] ; then
     ${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b || { echo "Cant copy regional.grid.b "; exit 1 ; }
   fi
   if [ ! -s regional.grid.a ] ; then
     ${pget} ${BASEDIR}//topo/regional.grid.a regional.grid.a || { echo "Cant copy regional.grid.a "; exit 1 ; }
   fi
   touch z_levitus
   if [ ! -s  z_levitus ] ; then
     ${pget} ${HYCOM_ALL}/relax/src/z_levitus .  || { echo "Cant copy z_levitus "; exit 1 ; }
   fi
   wait

   chmod a+rx z_levitus
   /bin/rm -f core
   touch core
   export FOR010A=fort.10A
   export FOR011A=fort.11A
   export FOR012A=fort.12A
   /bin/rm -f fort.10 fort.10A fort.11 fort.11A fort.12 fort.12A
   ./z_levitus <<E-o-D
    &AFTITL
     CTITLE = '1234567890123456789012345678901234567890',
     CTITLE = 'PHC     monthly',
    /
    &AFFLAG
     ICTYPE =   2,
     KSIGMA =   ${KSIGMA},
     INTERP =   1,
     ITEST  =  40,
     JTEST  =  17,
     MONTH  = $MM,
    /
E-o-D
   #
   # --- Required Output, potential density and temperature.
   #
   ${pput} fort.10  ${D}/temp_sig${KSIGMA}_m${MM}.b
   ${pput} fort.10A ${D}/temp_sig${KSIGMA}_m${MM}.a
   ${pput} fort.12  ${D}/dens_sig${KSIGMA}_m${MM}.b
   ${pput} fort.12A ${D}/dens_sig${KSIGMA}_m${MM}.a
   #
   # --- Optional Output.
   #
   ${pput} fort.11  ${D}/saln_sig${KSIGMA}_m${MM}.b
   ${pput} fort.11A ${D}/saln_sig${KSIGMA}_m${MM}.a
   #
   # --- end of month foreach loop
   #
   /bin/rm fort.7[123]
done
#
# --- Delete all files.
#
#/bin/rm -f *
