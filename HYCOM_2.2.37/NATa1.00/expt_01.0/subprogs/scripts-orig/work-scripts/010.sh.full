#!/bin/bash
#
# --- check that the C comment command is available.
#

suffix () {
   echo $i | sed "s/.*\.//g"
}

if  ! C >& /dev/null ;  then
  if [ -f ../../ALL/bin/C ] ; then
     cd ../../ALL/bin
     export PATH=$PATH:$(pwd)
     cd -
  else
    echo "Please put the command hycom/ALL/bin/C in your path"
    exit 1
  fi
fi
#
set -x
C
C --- Experiment GOMd0.08 - 17.X
C --- 20 layer HYCOM on 0.08 degree Gulf of Mexico region
C
C --- 06.1 - FCT2, KPP, tracer, 2.1.34,   nested inside ATLd0.08 2.6.
C --- 14.0 - FCT2, KPP, tracer, 2.2.06oo, nested inside ATLd0.08 2.6.
C --- 14.1 - Twin of 14.0 with            btr2nd=1
C --- 14.2 - Twin of 14.1 with  2.2.06q3, advflg=hybflg=0.0, temdfc=1.0.
C --- 14.7 - Twin of 14.2 with  2.2.08.
C --- 16.5 - Twin of 14.2 with  2.2.09,   hybiso=0.01 and hybmap=2.
C --- 17.2 - Twin of 14.2 with  2.2.12,   hybiso=0.01 and hybmap=2.
C
C --- Preamble, script keys on O/S name.
C
C --- Set parallel configuration, see ../README/README.expt_parallel.
C --- NOMP = number of OpenMP threads, 0 for no OpenMP, 1 for inactive OpenMP
C --- NMPI = number of MPI    tasks,   0 for no MPI
C
export NOMP=0
export NMPI=4
C
C --- modify NOMP and NMPI based on batch limits A
C --- KAL - skipped
C
C
C --- R is region name.
C --- V is source code version number.
C --- T is topography number.
C --- K is number of layers.
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C
source ../EXPT.src  # TODO - need better setup method, but use for now
export V=2.2.12
export K=`grep "'kdm   ' =" blk* | awk '{printf("%03d", $1)}'`
export P=$(pwd)/data
export D=$P
export S=/work/knutali/hycom/$R/expt_$X/SCRATCH
mkdir -p $D
mkdir -p $S
cd       $S

C
#C --- For whole year runs.
#C ---   Y01 initial model year of this run.
#C ---   YXX is the last model year of this run, and the first of the next run.
#C ---   A and B are identical, typically blank.
#C --- For part year runs.
#C ---   A is this part of the year, B is next part of the year.
#C ---   Y01 is the start model year of this run.
#C ---   YXX is the end   model year of this run, usually Y01.
#C --- For a few hour/day run
#C ---   A   is the start day and hour, of form "dDDDhHH".
#C ---   B   is the end   day and hour, of form "dXXXhYY".
#C ---   Y01 is the start model year of this run.
#C ---   YXX is the end   model year of this run, usually Y01.
#C --- Note that these variables are set by the .awk generating script.
#C
#setenv A "i"
#setenv B "j"
#setenv Y01 "099"
#setenv YXX "099"
#C
#echo "Y01 =" $Y01 "YXX = " $YXX  "A =" ${A} "B =" ${B}
#C
#C --- local input files.
#C
#if (-e ${D}/../${E}y${Y01}${A}.limits) then
#  /bin/cp ${D}/../${E}y${Y01}${A}.limits limits
#else
##  use "LIMITI"  when starting a run after day zero.
##  use "LIMITS9" (say) for a 9-day run
#  echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} >! limits
#endif
#cat limits
#C
# KAL - The above is for standard hycom - we use a slightly different approach
C
C
C --- Set up limits for this segment
C
C
export NUM_SEGMENT=1   # TODO set by python script or other
export THIS_SEGMENT=1  # TODO set by python script or other
if [ -f  ${D}/../limits_$(printf "%3.3d" ${THIS_SEGMENT}) ]; then
   cp  ${D}/../limits_$(printf "%3.3d" ${THIS_SEGMENT}) limits
else
   print " Can not find segment file  ${D}/../limits_$(printf "%3.3d" ${THIS_SEGMENT}) "
   exit 1
fi
if [ -s  ${D}/../ports.input_${THIS_SEGMENT} ] ; then
  /bin/cp ${D}/../ports.input_${THIS_SEGMENT} ports.input
else
  /bin/cp ${D}/../ports.input ports.input
fi
if [ -s  ${D}/../tracer.input_${THIS_SEGMENT} ] ; then
  /bin/cp ${D}/../tracer.input_${THIS_SEGMENT} tracer.input
else
  /bin/cp ${D}/../tracer.input tracer.input
fi
if [ -s  ${D}/../blkdat.input_${THIS_SEGMENT} ] ; then
  /bin/cp ${D}/../blkdat.input_${THIS_SEGMENT} tracer.input
else
  /bin/cp ${D}/../blkdat.input blkdat.input
fi

C
C
C --- Set up MPI partition for this segment
C
C
if [ $NMPI != 0 ] ; then
  export NPATCH=`echo $NMPI | awk '{printf("%04d", $1)}'`
  /bin/rm -f patch.input
  /bin/cp ${D}/../../topo/partit/depth_${R}_${T}.${NPATCH}  patch.input
# /bin/cp ${D}/../../topo/partit/depth_${R}_${T}.${NPATCH}u patch.input
fi

C
C --- check that iexpt from blkdat.input agrees with E from this script.
C
export EB=`grep "'iexpt ' =" blk* | awk '{printf("%03d", $1)}'`
if [ "$EB" != "$E" ] ; then
  cd $D/..
  /bin/mv LIST LIST_BADRUN
  echo "BADRUN" > LIST
  exit
fi

#C
#C --- turn on detailed debugging.
#C
#touch PIPE_DEBUG

C
C --- pget, pput "copy" files between scratch and permanent storage.
C --- Can both be cp if the permanent filesystem is mounted locally.
C --- KAL - we use pget=pput=cp
C
export pget=/bin/cp
export pput=/bin/cp


C
C --- input files from file server.
C
for i in regional.depth.a regional.depth.b ; do
   touch $i
   if [ ! -s  $i ] ; then
      ${pget} ${D}/../../topo/$i $i &
   fi
done
for i in regional.grid.a regional.grid.b ; do
   touch $i
   if [ ! -s  $i ] ; then
      ${pget} ${D}/../../topo/$i $i &
   fi
done

if [ ! -d  ./wind ]; then
   C
   C
   C --- Climatology atmospheric forcing 
   C
   C
   # FN Retrieved from EXPT.src
   #setenv FN ecmwf-reanal_7993-mn6hr
   for i in tauewd taunwd wndspd radflx shwflx vapmix \
      airtmp seatmp surtmp precip ; do
      # TODO - doesnt handle offsets
      touch forcing.$i.a forcing.$i.b
      if [  ! -s  forcing.$i.a  ] ; then
         ${pget} ${D}/../../force/${FN}/$i.a      forcing.$i.a &
      fi
      if [  ! -s  forcing.$i.b  ] ; then
         ${pget} ${D}/../../force/${FN}/$i.b      forcing.$i.b &
      fi
   done
fi

C
C --- time-invarent heat flux offset
C
#KAL - imported from source file
#setenv OF ""
#setenv OF "_413"
if [  "$OF" != "" ];  then
  touch  forcing.offlux.a  forcing.offlux.b
  if [ ! -s  forcing.offlux.a ] ; then
     ${pget} ${D}/../../force/${FN}/offlux${OF}.a forcing.offlux.a &
  fi
  if [ ! -s  forcing.offlux.b ] ; then
     ${pget} ${D}/../../force/${FN}/offlux${OF}.b forcing.offlux.b &
  fi 
fi 
C

C
C --- river forcing
C
touch  forcing.rivers.a  forcing.rivers.b
if [ ! -s  forcing.rivers.a ] ; then
   ${pget} ${D}/../../force/rivers/rivers_${T}.a forcing.rivers.a &
fi
if [ ! -s  forcing.rivers.b ] ; then
   ${pget} ${D}/../../force/rivers/rivers_${T}.b forcing.rivers.b &
fi

C
C --- kpar forcing
C
touch  forcing.kpar.a forcing.kpar.b
if [ ! -s forcing.kpar.a ] ; then
   ${pget} ${D}/../../force/seawifs/kpar.a forcing.kpar.a &
fi
if [ ! -s forcing.kpar.b ] ; then
   ${pget} ${D}/../../force/seawifs/kpar.b forcing.kpar.b &
fi


C
C --- relaxation
C
for i in rmu saln temp intf ; do
   touch relax.$i.a relax.$i.b
   j=$(echo $i | head -c3)
   if [ ! -s relax.$i.a ] ; then
      ${pget} ${D}/../../relax/${E}/relax_$j.a relax.$i.a  &
   fi
   if [ ! -s relax.$i.b ] ; then
      ${pget} ${D}/../../relax/${E}/relax_$j.b relax.$i.b  &
   fi
done

C
C - thermobaric reference state?
C
touch tbaric.a tbaric.b
if [ ! -s tbaric.a ]; then
   ${pget} ${D}/../../relax/${E}/tbaric.a tbaric.a  &
fi
if [ ! -s tbaric.b ]; then
   ${pget} ${D}/../../relax/${E}/tbaric.b tbaric.b  &
fi

C
C - Mystery field ....
C
touch iso.sigma.a iso.sigma.b
if [ ! -s iso.sigma.a ] ; then
   ${pget} ${D}/../../relax/${E}/iso_sigma.a iso.sigma.a  &
fi
if [ ! -s iso.sigma.b ] ; then
   ${pget} ${D}/../../relax/${E}/iso_sigma.b iso.sigma.b  &
fi

#C
#touch iso.top.a
#touch iso.top.b
#if (-z iso.top.a) then
#   ${pget} ${D}/../../relax/${E}/iso_top.a iso.top.a  &
#endif
#if (-z iso.top.b) then
#   ${pget} ${D}/../../relax/${E}/iso_top.b iso.top.b  &
#endif
#C
#touch thkdf4.a
#touch thkdf4.b
#if (-z thkdf4.a) then
#   ${pget} ${D}/../../relax/${E}/thkdf4.a thkdf4.a  &
#endif
#if (-z thkdf4.b) then
#   ${pget} ${D}/../../relax/${E}/thkdf4.b thkdf4.b  &
#endif



C
C --- restart input
C
touch   restart_in.a restart_in.b restart_out.a restart_out.b restart_out1.a restart_out1.b
if [ ! -s restart_in.b ]; then
  export RI="       0.00"
else
  export RI=`head -2 restart_in.b | tail -1 | awk  '{printf("%11.2f\n", $5)}'`
fi
if [ ! -s restart_out.b ] ;  then
  export RO="       0.00"
else
  export RO=`head -2 restart_out.b | tail -1 | awk  '{printf("%11.2f\n", $5)}'`
fi
if [ ! -s restart_out1.b ]; then
  export R1="       0.00"
else
  export R1=`head -2 restart_out1.b | tail -1 | awk  '{printf("%11.2f\n", $5)}'`
fi
export LI=`awk  '{printf("%11.2f\n", $1)}' limits`
C
if [ `echo $LI | awk '{if ($1 <= 0.0) print 1; else print 0}'` -eq 1 ];  then
   C --- no restart needed
   /bin/rm restart_in.a   restart_in.b
   /bin/rm restart_out.a  restart_out.b
   /bin/rm restart_out1.a restart_out1.b
elif [ `echo $LI $RI | awk '{if ($1-0.1 < $2 && $1+0.1 > $2) print 1; else print 0}'` -eq 1 ]; then
   C --- restart is already in restart_in
   /bin/rm restart_out.a  restart_out.b
   /bin/rm restart_out1.a restart_out1.b
elif [ `echo $LI $RO | awk '{if ($1-0.1 < $2 && $1+0.1 > $2) print 1; else print 0}'` -eq 1 ] ; then
   C --- restart is in restart_out
   /bin/mv restart_out.a  restart_in.a
   /bin/mv restart_out.b  restart_in.b
   /bin/rm restart_out1.a restart_out1.b
elif [ `echo $LI $R1 | awk '{if ($1-0.1 < $2 && $1+0.1 > $2) print 1; else print 0}'` -eq 1 ]; then
   C ---   restart is in restart_out1
   /bin/mv restart_out1.a restart_in.a
   /bin/mv restart_out1.b restart_in.b
   /bin/rm restart_out.a  restart_out.b
else
   C ---   get restart from permenant storage
   /bin/rm restart_in.a   restart_in.b
   /bin/rm restart_out.a  restart_out.b
   /bin/rm restart_out1.a restart_out1.b
   ${pget} ${D}/restart_in_${LAST_SEGMENT}.a restart_in.a &
   ${pget} ${D}/restart_in_${LAST_SEGMENT}.b restart_in.b &
fi

C
C --- model executable
C
if  [ $NMPI -eq 0 -a $NOMP -eq 0 ]; then
   export TYPE=one
elif [ $NMPI -eq 0 ];  then
   export TYPE=omp
elif [ $NOMP -eq 0 ]; then
   export TYPE=mpi
else
   export TYPE=ompi
fi
/bin/cp ${D}/../../Build_V${V}_X${X}/hycom . &

C
C --- summary printout
C
touch   summary_out
/bin/mv summary_out summary_old
C
C --- heat transport output
C
touch   flxdp_out.a flxdp_out.b
/bin/mv flxdp_out.a flxdp_old.a
/bin/mv flxdp_out.b flxdp_old.b
C
touch   ovrtn_out
/bin/mv ovrtn_out ovrtn_old
C
C --- clean up old archive files, typically from batch system rerun.
C
mkdir KEEP
touch archv.dummy.b
for f  in arch* ; do
  /bin/mv $f KEEP/$f
done
#TODO - check out nesting approach
#C
#C --- Nesting input archive files.
#C
#if (-e ./nest) then
#  setenv NA ${Y01}${A}
#  setenv NB ${YXX}${B}
#  cd ./nest
#  touch rmu.a rmu.b
#  if (-z rmu.a) then
#     ${pget} ${D}/../../relax/${E}/nest_rmu.a rmu.a &
#  endif
#  if (-z rmu.b) then
#     ${pget} ${D}/../../relax/${E}/nest_rmu.b rmu.b &
#  endif
#  touch   arch.dummy.b
#  /bin/rm arch*.[ab]
#  touch archv_${NA}.tar
#  if (-z archv_${NA}.tar) then
#    ${pget} ${D}/nest/archv_${NA}.tar archv_${NA}.tar
#  endif
#  tar xvf archv_${NA}.tar
#  cd ..
#endif

C
C --- let all file copies complete.
C
wait
C
C --- zero file length means no rivers.
C
if [ ! -s forcing.rivers.a ] ; then
   /bin/rm forcing.rivers.[ab]
fi

C
C --- Just in time atmospheric forcing.
C
if [ -d  ./wind ]; then
   if [ ! -d ./flux ] ; then
     echo './flux must exist if ./wind does'
     exit 1
   fi   
   C
   C --- Check to see if wind and flux files exist, if not make them and wait.
   C
   /bin/rm -f forcing.tauewd.a forcing.taunwd.a forcing.wndspd.a
   /bin/rm -f forcing.tauewd.b forcing.taunwd.b forcing.wndspd.b
   /bin/rm -f ./wind/${E}_${THIS_SEGMENT}
   if [ -s ./wind/tauewd_${THIS_SEGMENT}.a -a \
        -s ./wind/tauewd_${THIS_SEGMENT}.b -a \
        -s ./wind/taunwd_${THIS_SEGMENT}.a -a \
        -s ./wind/taunwd_${THIS_SEGMENT}.b -a \
        -s ./wind/wndspd_${THIS_SEGMENT}.a -a \
        -s ./wind/wndspd_${THIS_SEGMENT}.b    ] ; then
     /bin/ln ./wind/tauewd_${THIS_SEGMENT}.a forcing.tauewd.a
     /bin/ln ./wind/taunwd_${THIS_SEGMENT}.a forcing.taunwd.a
     /bin/ln ./wind/wndspd_${THIS_SEGMENT}.a forcing.wndspd.a
     /bin/ln ./wind/tauewd_${THIS_SEGMENT}.b forcing.tauewd.b
     /bin/ln ./wind/taunwd_${THIS_SEGMENT}.b forcing.taunwd.b
     /bin/ln ./wind/wndspd_${THIS_SEGMENT}.b forcing.wndspd.b
   else
     cd ./wind
     touch ${E}_${THIS_SEGMENT}
     #    touch ${E}w${Y01}${A}
     #    /bin/rm -f ${E}w${Y01}${A}.com ${E}w${Y01}${A}.log
     #    awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}W.com > ${E}w${Y01}${A}.com
     #    csh ${E}w${Y01}${A}.com >& ${E}w${Y01}${A}.log &
     # TODO - create wind forcing here ...
     cd ..
   fi    
   if [ -d ./ssta ] ; then
     /bin/rm -f forcing.surtmp.a forcing.surtmp.b ./ssta/${E}_segment${THIS_SEGMENT}
     if [ -s ./ssta/surtmp_${THIS_SEGMENT}.a -a \
          -s ./ssta/surtmp_${THIS_SEGMENT}.b    ] ;  then
       /bin/ln ./ssta/surtmp_${THIS_SEGMENT}.a forcing.surtmp.a
       /bin/ln ./ssta/surtmp_${THIS_SEGMENT}.b forcing.surtmp.b
     else
       cd ./ssta
       touch ${E}_${THIS_SEGMENT}
       #touch ${E}t${Y01}${A}
       #/bin/rm -f ${E}t${Y01}${A}.com ${E}t${Y01}${A}.log
       #awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}T.com > ${E}t${Y01}${A}.com
       #csh ${E}t${Y01}${A}.com >& ${E}t${Y01}${A}.log &
       # TODO - create surtmp forcing here ...
       cd ..
     fi
  fi

  if [ -d ./ssto ] ; then
     /bin/rm -f forcing.seatmp.a forcing.seatmp.b  ./ssto/${E}_${THIS_SEGMENT}
     if [ -s ./ssto/seatmp_${THIS_SEGMENT}.a -a \
          -s ./ssto/seatmp_${THIS_SEGMENT}.b    ]; then
       /bin/ln ./ssto/seatmp_${THIS_SEGMENT}.a forcing.seatmp.a
       /bin/ln ./ssto/seatmp_${THIS_SEGMENT}.b forcing.seatmp.b
     else
       cd ./ssto
       touch ${E}_${THIS_SEGMENT}
       #/bin/rm -f ${E}o${Y01}${A}.com ${E}o${Y01}${A}.log
       #awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}O.com > ${E}o${Y01}${A}.com
       #csh ${E}o${Y01}${A}.com >& ${E}o${Y01}${A}.log &
       # TODO - create seatmp forcing here ...
       cd ..
     fi
   fi

   if [  -d ./pcip ]; then
     /bin/rm -f forcing.precip.a forcing.precip.b  ./pcip/${E}p${Y01}${A}
     if [ -s ./pcip/precip_${THIS_SEGMENT}.a -a \
          -s ./pcip/precip_${THIS_SEGMENT}.b    ] ; then
       /bin/ln ./pcip/precip_${THIS_SEGMENT}.a forcing.precip.a
       /bin/ln ./pcip/precip_${THIS_SEGMENT}.b forcing.precip.b
     else
       cd ./pcip
       touch ${E}_${THIS_SEGMENT}
       #/bin/rm -f ${E}p${Y01}${A}.com ${E}p${Y01}${A}.log
       #awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}P.com > ${E}p${Y01}${A}.com
       #csh ${E}p${Y01}${A}.com >& ${E}p${Y01}${A}.log &
       # TODO - create precip forcing here ...
       cd ..
     fi
  fi

  if [ -d ./flxt ]; then
     /bin/rm -f forcing.radflx.a forcing.radflx.b ./flxt/${E}_${THIS_SEGMENT}
     if [ -s ./flxt/totflx_${Y01}${A}.a -a \
          -s ./flxt/totflx_${Y01}${A}.b    ] ; then
       /bin/ln ./flxt/totflx_${Y01}${A}.a forcing.radflx.a
       /bin/ln ./flxt/totflx_${Y01}${A}.b forcing.radflx.b
     else
       cd ./flxt
       touch ${E}_${THIS_SEGMENT}
       #touch ${E}q${Y01}${A}
       #/bin/rm -f ${E}p${Y01}${A}.com ${E}p${Y01}${A}.log
       #awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}Q.com > ${E}q${Y01}${A}.com
       #csh ${E}q${Y01}${A}.com >& ${E}q${Y01}${A}.log &
       # TODO - create flxt forcing here ...
       cd ..
     fi
   fi

   if [ -d ./flux ] ; then
     /bin/rm -f forcing.airtmp.a forcing.shwflx.a forcing.vapmix.a
     /bin/rm -f forcing.airtmp.b forcing.shwflx.b forcing.vapmix.b
     if [ ! -d ./pcip ]; then
       /bin/rm -f forcing.precip.a forcing.precip.b
       touch ./flux/precip_${THIS_SEGMENT}.a
       touch ./flux/precip_${THIS_SEGMENT}.b
     fi
     if [ ! -d ./flxt] ; then
       /bin/rm -f forcing.radflx.a forcing.radflx.b
       touch ./flux/radflx_${THIS_SEGMENT}.a
       touch ./flux/radflx_${THIS_SEGMENT}.b
     fi
     /bin/rm -f ./flux/${E}_${THIS_SEGMENT}
     if [-s ./flux/airtmp_${THIS_SEGMENT}.a -a \
         -s ./flux/airtmp_${THIS_SEGMENT}.b -a \
         -s ./flux/precip_${THIS_SEGMENT}.a -a \
         -s ./flux/precip_${THIS_SEGMENT}.b -a \
         -s ./flux/radflx_${THIS_SEGMENT}.a -a \
         -s ./flux/radflx_${THIS_SEGMENT}.b -a \
         -s ./flux/shwflx_${THIS_SEGMENT}.a -a \
         -s ./flux/shwflx_${THIS_SEGMENT}.b -a \
         -s ./flux/vapmix_${THIS_SEGMENT}.a -a \
         -s ./flux/vapmix_${THIS_SEGMENT}.b    ] ; then
       /bin/ln ./flux/airtmp_${THIS_SEGMENT}.a forcing.airtmp.a
       /bin/ln ./flux/airtmp_${THIS_SEGMENT}.b forcing.airtmp.b
       /bin/ln ./flux/shwflx_${THIS_SEGMENT}.a forcing.shwflx.a
       /bin/ln ./flux/shwflx_${THIS_SEGMENT}.b forcing.shwflx.b
       /bin/ln ./flux/vapmix_${THIS_SEGMENT}.a forcing.vapmix.a
       /bin/ln ./flux/vapmix_${THIS_SEGMENT}.b forcing.vapmix.b
       if [ ! -d ./pcip ]; then
         /bin/ln ./flux/precip_${THIS_SEGMENT}.b forcing.precip.b
         /bin/ln ./flux/precip_${THIS_SEGMENT}.a forcing.precip.a
       fi
       if [ ! -d ./flxt ]; then
         /bin/ln ./flux/radflx_${THIS_SEGMENT}.a forcing.radflx.a
         /bin/ln ./flux/radflx_${THIS_SEGMENT}.b forcing.radflx.b
       fi
     else
       cd ./flux
       touch ${E}_${THIS_SEGMENT}
       #/bin/rm -f ${E}f${Y01}${A}.com ${E}f${Y01}${A}.log
       #awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}F.com > ${E}f${Y01}${A}.com
       #csh ${E}f${Y01}${A}.com >& ${E}f${Y01}${A}.log &
       # TODO - create flux forcing here ...
       cd ..
     fi
   fi
   wait
C
   if [ -f ./wind/${E}_${THIS_SEGMENT} ]; then
     /bin/ln ./wind/tauewd_${THIS_SEGMENT}.a forcing.tauewd.a
     /bin/ln ./wind/taunwd_${THIS_SEGMENT}.a forcing.taunwd.a
     /bin/ln ./wind/wndspd_${THIS_SEGMENT}.a forcing.wndspd.a
     /bin/ln ./wind/tauewd_${THIS_SEGMENT}.b forcing.tauewd.b
     /bin/ln ./wind/taunwd_${THIS_SEGMENT}.b forcing.taunwd.b
     /bin/ln ./wind/wndspd_${THIS_SEGMENT}.b forcing.wndspd.b
   fi
C
   if [ -d ./ssta ]; then
     if [ -f ./ssta/${E}_${THIS_SEGMENT} ]; then
       /bin/ln ./ssta/surtmp_${THIS_SEGMENT}.a forcing.surtmp.a
       /bin/ln ./ssta/surtmp_${THIS_SEGMENT}.b forcing.surtmp.b
     fi
   fi
C
   if [ -d ./ssto ];  then
     if [ -f ./ssto/${E}_${THIS_SEGMENT} ]; then
       /bin/ln ./ssto/seatmp_${THIS_SEGMENT}.a forcing.seatmp.a
       /bin/ln ./ssto/seatmp_${THIS_SEGMENT}.b forcing.seatmp.b
     fi
   fi
C
   if [ -d ./pcip ] ; then
     if [ -f ./pcip/${E}_${THIS_SEGMENT} ]; then
       /bin/ln ./pcip/precip_${THIS_SEGMENT}.a forcing.precip.a
       /bin/ln ./pcip/precip_${THIS_SEGMENT}.b forcing.precip.b
     fi
   fi
C
   if [ -d ./flxt ]; then
     if [ -f ./flxt/${E}_${THIS_SEGMENT} ]; then
       /bin/ln ./flxt/totflx_${THIS_SEGMENT}.a forcing.radflx.a
       /bin/ln ./flxt/totflx_${THIS_SEGMENT}.b forcing.radflx.b
     fi
   fi
C
   if [ -d ./flux ]; then
     if [ -f ./flux/${E}_${THIS_SEGMENT} ]; then
       /bin/ln ./flux/airtmp_${THIS_SEGMENT}.a forcing.airtmp.a
       /bin/ln ./flux/airtmp_${THIS_SEGMENT}.b forcing.airtmp.b
       /bin/ln ./flux/shwflx_${THIS_SEGMENT}.a forcing.shwflx.a
       /bin/ln ./flux/shwflx_${THIS_SEGMENT}.b forcing.shwflx.b
       /bin/ln ./flux/vapmix_${THIS_SEGMENT}.a forcing.vapmix.a
       /bin/ln ./flux/vapmix_${THIS_SEGMENT}.b forcing.vapmix.b
       if [ ! -d ./pcip ]; then
         /bin/ln ./flux/precip_${THIS_SEGMENT}.b forcing.precip.b
         /bin/ln ./flux/precip_${THIS_SEGMENT}.a forcing.precip.a
       fi
       if [ ! -d ./flxt ]; then
         /bin/ln ./flux/radflx_${THIS_SEGMENT}.a forcing.radflx.a
         /bin/ln ./flux/radflx_${THIS_SEGMENT}.b forcing.radflx.b
       fi
     fi
   fi
C
C
C --- If the winds or fluxes for the next segment dont exist, 
C --- interpolate them to model grid when current segment is running.
C
   if [ -s ./wind/tauewd_${NEXT_SEGMENT}.a -a \
        -s ./wind/tauewd_${NEXT_SEGMENT}.b -a \
        -s ./wind/taunwd_${NEXT_SEGMENT}.a -a \
        -s ./wind/taunwd_${NEXT_SEGMENT}.b -a \
        -s ./wind/wndspd_${NEXT_SEGMENT}.a -a \
        -s ./wind/wndspd_${NEXT_SEGMENT}.b    ] ; then
C
C --- next segments winds already exist.
C
   else
     cd ./wind
     #/bin/rm -f ${E}w${YXX}${B}.com ${E}w${YXX}${B}.log
     #awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}W.com > ${E}w${YXX}${B}.com
     #csh ${E}w${YXX}${B}.com >& ${E}w${YXX}${B}.log &
     # TODO - generate forcing for next segment
     cd ..
   fi

   if [ -d ./ssta ];  then
     if [ -s ./ssta/surtmp_${NEXT_SEGMENT}.a -a \
          -s ./ssta/surtmp_${NEXT_SEGMENT}.b    ] ; then
C
C ---   next segments surtmp already exist.
C
     else
       cd ./ssta
       #/bin/rm -f ${E}t${YXX}${B}.com ${E}t${YXX}${B}.log
       #awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}T.com > ${E}t${YXX}${B}.com
       #csh ${E}t${YXX}${B}.com >& ${E}t${YXX}${B}.log &
       # TODO - generate ssta forcing for next segment
       cd ..
     fi
   fi    
   if [ -d ./ssto]; then
     if [-s ./ssto/seatmp_${NEXT_SEGMENT}.a -a \
         -s ./ssto/seatmp_${NEXT_SEGMENT}.b    ]; then
C
C ---   next segments seatmp already exist.
C
     else
       cd ./ssto
       #/bin/rm -f ${E}o${YXX}${B}.com ${E}o${YXX}${B}.log
       #awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}O.com > ${E}o${YXX}${B}.com
       #csh ${E}o${YXX}${B}.com >& ${E}o${YXX}${B}.log &
       # TODO - generate seatmp forcing for next segment
       cd ..
     fi
   fi
   if [ -d ./pcip ]; then
     if [ -s ./pcip/precip_${NEXT_SEGMENT}.a -a \
          -s ./pcip/precip_${NEXT_SEGMENT}.b    ]; then
C
C ---   next segments pcip already exist.
C
     else
       cd ./pcip
       #/bin/rm -f ${E}p${YXX}${B}.com ${E}p${YXX}${B}.log
       #awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}P.com > ${E}p${YXX}${B}.com
       #csh ${E}p${YXX}${B}.com >& ${E}p${YXX}${B}.log &
       # TODO - generate seatmp forcing for next segment
       cd ..
     fi
   fi
   if [ -d ./flxt ]; then
     if [ -s ./flxt/totflx_${NEXT_SEGMENT}.a -a \
          -s ./flxt/totflx_${NEXT_SEGMENT}.b    ];  then
C
C ---   next segments flxt already exist.
C
     else
       cd ./flxt
       #/bin/rm -f ${E}q${YXX}${B}.com ${E}q${YXX}${B}.log
       #awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}Q.com > ${E}q${YXX}${B}.com
       #csh ${E}q${YXX}${B}.com >& ${E}q${YXX}${B}.log &
       # TODO - generate flxt forcing for next segment
       cd ..
     fi    
   fi

   if [ -d ./flux ]; then
     if [ ! -d ./pcip ]; then
       touch ./flux/precip_${YXX}${B}.a
       touch ./flux/precip_${YXX}${B}.b
     fi
     if [ ! -d ./flxt ];  then
       touch ./flux/radflx_${YXX}${B}.a
       touch ./flux/radflx_${YXX}${B}.b
     fi
     if [ -s ./flux/airtmp_${YXX}${B}.a -a \
          -s ./flux/airtmp_${YXX}${B}.b -a \
          -s ./flux/precip_${YXX}${B}.a -a \
          -s ./flux/precip_${YXX}${B}.b -a \
          -s ./flux/radflx_${YXX}${B}.a -a \
          -s ./flux/radflx_${YXX}${B}.b -a \
          -s ./flux/shwflx_${YXX}${B}.a -a \
          -s ./flux/shwflx_${YXX}${B}.b -a \
          -s ./flux/vapmix_${YXX}${B}.a -a \
          -s ./flux/vapmix_${YXX}${B}.b    ]; then
C
C ---   next segments fluxes already exist.
C
     else
       cd ./flux
       #/bin/rm -f ${E}f${YXX}${B}.com ${E}f${YXX}${B}.log
       #awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}F.com > ${E}f${YXX}${B}.com
       #csh ${E}f${YXX}${B}.com >& ${E}f${YXX}${B}.log &
       # TODO - generate flux forcing for next segment
       cd ..
     fi
  fi
fi # -d ./wind

#C
#C --- Nesting input archive files for next segment.
#C
#if (-e ./nest) then
#  cd ./nest
#  touch archv_${NB}.tar
#  if (-z archv_${NB}.tar) then
#    ${pget} ${D}/nest/archv_${NB}.tar archv_${NB}.tar &
#  endif
#  cd ..
#endif
#C
#chmod ug+x hycom
#/bin/ls -laFq
#C
#if (-e ./nest) then
#  ls -laFq nest
#endif

C
C ---  Check to make sure restart file is there
C
if [ `echo $LI | awk '{print ($1 > 0.0)}'` == 1   -a ! -s restart_in.a ];  then
  cd $D/..
  /bin/mv LIST LIST_BADRUN
  echo "BADRUN" > LIST
  exit
fi
C

if [ $NMPI -eq  0 ] ;  then
C
C --- run the model, without MPI or SHMEM
C
   if [ $NOMP -eq 0 ]; then
     setenv NOMP 1
   fi
   C
   C - Removed OS switches for OMP
   C
else
   C
   C --- run the model, with MPI or SHMEM and perhaps also with OpenMP.
   C
   touch patch.input

   if [ ! -s  patch.input ] ; then
     C
     C --- patch.input is always required when using MPI or SHMEM.
     C
     cd $D/..
     /bin/mv LIST LIST_BADRUN
     echo "BADRUN" > LIST
     exit
   fi
fi





#######################################################################
# Start post-processing


touch   PIPE_DEBUG
/bin/rm PIPE_DEBUG

C
C --- archive output in a separate tar directory
C
touch archv.dummy.a archv.dummy.b archv.dummy.txt
C
if [ -d ./SAVE ] ; then
  for f  in archv.*.a  archv.*.b archv.*.txt ; do
    /bin/ln ${f} SAVE/${f}
  done
fi
C
mkdir ./tar_${THIS_SEGMENT}
C
for f in archv.*.a archv.*.b archv.*.txt ; do
  /bin/mv ${f} ./tar_${THIS_SEGMENT}/${E}_${f}
done
C
C --- build and run or submit the tar script
C
#awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}A.com >! tar_${Y01}${A}.com
#csh tar_${Y01}${A}.com >&! tar_${Y01}${A}.log &
##llsubmit ./tar_${Y01}${A}.com
C
C --- heat transport statistics output
C
if [ -s flxdp_out.a ] ; then
  ${pput} flxdp_out.a ${S}/flxdp_${THIS_SEGMENT}.a
fi
if [ -s flxdp_out.b ] ; then
  ${pput} flxdp_out.b ${S}/flxdp_${THIS_SEGMENT}.b
fi
if [ -s ovrtn_out ] ; then
  ${pput} ovrtn_out ${S}/ovrtn_out_${THIS_SEGMENT}.a
fi

C
C --- restart output
C
if [ -s  restart_out.a ] ; then
  ${pput} restart_out.a ${S}/restart_${NEXT_SEGMENT}.a
fi
if [ -s restart_out.b ] ; then
  ${pput} restart_out.b ${S}/restart_${NEXT_SEGMENT}.b
fi

C
if [ -d ./wind ] ; then
   C
   C --- Delete just in time wind and flux files.
   C
   touch summary_out
   tail -1 summary_out
   tail -1 summary_out | grep -c "^normal stop"
   if ( `tail -1 summary_out | grep -c "^normal stop"` == 1 ) then
     /bin/rm -f ./wind/*_${THIS_SEGMENT}*
     /bin/rm -f ./flux/*_${THIS_SEGMENT}*
     /bin/rm -f ./flxt/*_${THIS_SEGMENT}*
     /bin/rm -f ./pcip/*_${THIS_SEGMENT}*
     /bin/rm -f ./ssta/*_${THIS_SEGMENT}*
     /bin/rm -f ./ssto/*_${THIS_SEGMENT}*
   fi
   C
   C
   C --- Wait for wind and flux interpolation of next segment.
   C
   wait
fi

#C
#C --- wait for nesting .tar file.
#C
#if (-e ./nest) then
#  wait
#endif

C
C --- HYCOM error stop is implied by the absence of a normal stop.
C
touch summary_out
tail -1 summary_out
tail -1 summary_out | grep -c "^normal stop"
if  [ `tail -1 summary_out | grep -c "^normal stop"` == 0 ] ; then
  cd $D/..
  /bin/mv LIST LIST_BADRUN
  echo "BADRUN" > LIST
fi
C
C --- wait for tar bundles to complete
C
wait
C
C  --- END OF MODEL RUN SCRIPT
C
