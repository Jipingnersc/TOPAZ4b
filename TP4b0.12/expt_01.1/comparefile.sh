

Fdir2='/home/pr2n0108/TOPAZ/TP4b0.12/expt_02.1'

Fdir1='/scratch-shared/xiejp/TP4b0.12/expt_02.1'

cd ${Fdir1}

for ifile in $(ls ${Fdir1}/*.*) ; do
  #11    if [ ${ifile}]
   #echo ${ifile%\/*}
   tmpstr=${ifile%\/*}
   filestr=${ifile#${tmpstr}\/}
   if [ -s ${filestr} -a -s ${Fdir2}/${filestr} ]; then
      echo
      echo 'Comparing two files: ' ${filestr}
      diff ${filestr} ${Fdir2}/${filestr} 
   fi


done
