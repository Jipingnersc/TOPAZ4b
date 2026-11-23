
# Need to define the source files for combining and their location
# 
# Check their permittions for these file
# 
Sdir=/work/cagyum/TOPAZ_DAY092/
hycrst_phy=TP4restart2017_092_00_mem002
hycrst_eco=TP0restart2017_092_00



# clean the possible temperary files
rm TMPrst_*.? analysis*.?


# check the file existing
Fini1=${Sdir}${hycrst_phy}
Fini2=${Sdir}${hycrst_eco}
if [ ! -s ${Fini1}.a -o ! -s ${Fini1}.b ]; then
  echo "Cannot find " ${Fini1}
  exit 0
elif [ ! -s ${Fini2}.a -o ! -s ${Fini2}.b ]; then
  echo "Cannot find " ${Fini2}
  exit 0
fi

myrandom="$RANDOM_$RANDOM"
# 0. prepare the header lines for the target .b files
head -n2 ${Fini1}.b > tmp${myrandom}.b

# 1. separate the restart file into 2/3D files
Ftmpsur=analysis_proc
Ftmplast=analysis
ii=0
hycseparate ${Fini1}.a
for ityp in 3D 2D; do
  if [ -s TMPrst_${ityp}.a -a -s TMPrst_${ityp}.b ]; then
    CN=$(echo 00${ii}|tail -4c)
    mv TMPrst_${ityp}.a ${Ftmpsur}${CN}.a 
    cat tmp${myrandom}.b > ${Ftmpsur}${CN}.b
    cat TMPrst_${ityp}.b >> ${Ftmpsur}${CN}.b
    let ii=ii+2
    rm TMPrst_${ityp}.?
  fi
done

ii=1
hycseparate ${Fini2}.a
for ityp in 3D 2D; do
  if [ -s TMPrst_${ityp}.a -a -s TMPrst_${ityp}.b ]; then
    CN=$(echo 00${ii}|tail -4c)
    mv TMPrst_${ityp}.a ${Ftmpsur}${CN}.a 
    cat tmp${myrandom}.b > ${Ftmpsur}${CN}.b
    cat TMPrst_${ityp}.b >> ${Ftmpsur}${CN}.b
    let ii=ii+2
    rm TMPrst_${ityp}.?
  fi
done


# 2. combine the separated files into one
hycassemble.sh 4

# clean the temperary files
if [ -s ${Ftmplast}.a -a -s ${Ftmplast}.b ]; then
  rm ${Ftmpsur}*.?
  rm tmp*.?
fi

echo $?

