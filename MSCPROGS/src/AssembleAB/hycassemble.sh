#!/bin/bash
usage="####################################################################\n\
This routine collects the analysed fields, and assembles \n\
them into files of name analysis.[ab] used as restart files by HYCOM. 
The unassembled fields have names of type\n\
analysis_procXXX.[ab] \n\n\
Usage:   $(basename $0)  nproc \n\
Where: \n\
\t\"nproc\"            is the number of the combined files  \n\
\n\
Example:\n\
\t hycassemble.sh  3 
\n\n\
NB:\n\
####################################################################\n"


prog="$(dirname $0)/hycassemble"

[ $# -ne 1 ] && { echo -e $usage ; exit 1 ; }

# order of the restart fields
$prog $@
[ $? -ne 0 ]  && { echo hyc_rst_assemble failed ; exit 1 ; }

myrandom="$RANDOM_$RANDOM"
head -n2 analysis_proc000.b > tmp${myrandom}.b

cat analysis.b >> tmp${myrandom}.b
mv tmp${myrandom}.b analysis.b

echo $?
