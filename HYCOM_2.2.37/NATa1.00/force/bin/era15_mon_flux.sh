#!/bin/bash
#

echo " Script needs setup"
exit
set -x
pget=cp
pput=cp
export X="NOEXPERIMENT" # dummy experiment


# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
cd -

source ${BASEDIR}/EXPT.src $X || { echo "Could not source ${BASEDIR}/EXPT.src" ; exit 1 ; }

# Check that pointer to HYCOM_ALL is set (from EXPT.src)
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

#
#
# --- Create HYCOM heat flux files.
# --- Zero W/m**2 radiation correction, +ve means gained by ocean.
#
#
# --- P is primary path,
# --- S is scratch directory,
# --- D is permanent directory,
#
D=$BASEDIR/force/ecmwf-reanal_7993-mn6hr/ # TODO: Must be fixed
S=$D/SCRATCH
mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

echo "KAL -- routine disabled for testing(e.g use at own risk) " && exit 1

#
# --- Input.
#
touch      fort.71
/bin/rm -f fort.71
${pget} $ERA15/era15-sea_1979-1993_mon_TaqaQrQp.D fort.71 || { echo "Could not get era15 file " ; exit 1 ; }

#
touch regional.grid.a regional.grid.b
rm regional.grid.a regional.grid.b
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b || { echo "Could not get regional.grid.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a  || { echo "Could not get regional.grid.a file " ; exit 1 ; }
#cp ${D}/../../../ALL/force/src/ap_100_co . || { echo "Could not get routine ap_100_co " ; exit 1 ; }
cp ${HYCOM_ALL}/force/src/ap . || { echo "Could not get routine ap " ; exit 1 ; }
wait

#chmod a+rx ap_100_co
chmod a+rx ap
/bin/rm -f fort.1[01234]*

export FOR010A=fort.10A
export FOR011A=fort.11A
export FOR012A=fort.12A
export FOR013A=fort.13A
export FOR014A=fort.14A

/bin/rm -f core
touch core

#./ap_100_co <<E-o-D
./ap <<E-o-D
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'ERA15 monthly, MKS',
 /
 &AFTIME
  BIASPC =   0.0,
  BIASRD =   0.0,
 /
 &AFFLAG
  IFTYPE =   4,
  INTERP =   1,
 /
E-o-D
#
# --- Output.
#
${pput} fort.10  ${D}/airtmp.b
${pput} fort.10A ${D}/airtmp.a
${pput} fort.11  ${D}/vapmix.b
${pput} fort.11A ${D}/vapmix.a
${pput} fort.12  ${D}/radflx.b
${pput} fort.12A ${D}/radflx.a
${pput} fort.13  ${D}/shwflx.b
${pput} fort.13A ${D}/shwflx.a
#${pput} fort.14  ${D}/precip.b
#${pput} fort.14A ${D}/precip.a

#
# --- Input.
#
/bin/rm -f fort.1[01234]*
touch      fort.71
/bin/rm -f fort.71
${pget} $ERA15/era15-sea_1979-1993_mon_ttlpcp.D fort.71 || { echo "Could not get era15 file " ; exit 1 ; }

#./ap_100_co <<E-o-D
./ap <<E-o-D
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'ERA15 monthly, MKS',
 /
 &AFTIME
  BIASPC =   0.0,
  BIASRD =   0.0,
 /
 &AFFLAG
  IFTYPE =   1,
  INTERP =   1,
 /
E-o-D
#
# --- Output.
#
${pput} fort.14  ${D}/precip.b
${pput} fort.14A ${D}/precip.a



#
# --- Delete all files.
#
#/bin/rm -f *
echo
echo "Copied forcing files to $D"
