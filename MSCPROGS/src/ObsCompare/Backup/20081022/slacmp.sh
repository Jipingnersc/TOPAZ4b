#!/bin/bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/fimm/nersc/knutali/opt/libnc-dap/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/fimm/nersc/knutali/opt/libdap-3.7.10/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/fimm/nersc/knutali/opt/curl-7.18.0/lib/


${0%.sh}  $@
