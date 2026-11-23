#
#@$-lt 0:04:00
#@$-lT 0:04:00
#@$-lm 4mw
#@$-lM 4mw
#@$-s  /bin/csh
#@$
#
set echo
set time = 1
C
setenv OS `uname`
switch ($OS)
case 'SunOS':
    breaksw
case 'Linux':
    breaksw
case 'OSF1':
    breaksw
case 'IRIX64':
    breaksw
case 'AIX':
    breaksw
case 'unicosmk':
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
C
C --- Create HYCOM kpar files.
C
switch ($OS)
case 'SunOS':
    if (-e /scr) then
      setenv pget ~wallcraf/bin/pget
      setenv pput ~wallcraf/bin/pput
    else
      setenv pget cp
      setenv pput cp
    endif
    breaksw
case 'XXX':
    setenv pget ~wallcraf/bin/pget
    setenv pput ~wallcraf/bin/pput
    breaksw
default:
    setenv pget cp
    setenv pput cp
endsw
C
C --- P is primary path,
C --- S is scratch directory,
C --- D is permanent directory,
C
setenv P hycom/GLBt0.72/force/seawifs
C
switch ($OS)
case 'SunOS':
    if (-e /scr) then
#                  NAVO MSRC
      setenv S /scr/${user}/$P
      setenv D            ~/$P
    else if (-e /net/ajax/data) then
#                  NRLSSC
      setenv S /net/ajax/data/${user}/$P
      setenv D                      ~/$P
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv D ~/$P
    endif
    breaksw
case 'Linux':
    if (-e /external/fast) then
#                  NRLSSC
      setenv S /external/fast/${user}/$P
      setenv D                      ~/$P
    else if (-e /net/ajax/data) then
#                  NRLSSC
      setenv S /net/ajax/data/${user}/$P
      setenv D                      ~/$P
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv D ~/$P
    endif
    breaksw
default:
#                  Single Disk
    setenv S ~/$P/SCRATCH
    setenv D ~/$P
endsw
C
mkdir -p $S
cd       $S
C
C --- Input.
C
touch      fort.71
/bin/rm -f fort.71
${pget} ~wallcraf/flux_ieee/seawifs/seawifs_mon_kpar.d fort.71 &
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C
#${pget} ${D}/../../../ALL/force/src/kp . &
cp ${D}/../../../ALL/force/src/kp . &
wait
chmod a+rx kp
C
/bin/rm -f fort.1[01234]*
C
setenv FOR010A fort.10A
C
switch ($OS)
case 'SunOS':
    /bin/rm -f core kp.trace
    touch core
default:
    /bin/rm -f core
    touch core
endsw
C
./kp <<'E-o-D'
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'SeaWifs 1997:2001 monthly clim, 1/m',
 /
 &AFTIME
  PARMIN =    0.04, !minimum kpar (default 0.04)
  PARMAX =    0.20  !maximum kpar (default 0.20)
 /
 &AFFLAG
  INTERP =    0,  !0:piecewise-linear; 1:cubic-spline;
  INTMSK =    0,  !0:no mask; 1:land/sea=0/1; 2:land/sea=1/0;
 /
'E-o-D'
C
switch ($OS)
case 'SunOS':
    if (-e kp.trace) cat kp.trace
    breaksw
default:
endsw
C
C --- Output.
C
mv fort.10  kpar.b
mv fort.10A kpar.a
#${pput} kpar.b ${D}/kpar.b
#${pput} kpar.a ${D}/kpar.a
C
C --- Delete all files.
C
#/bin/rm -f *
