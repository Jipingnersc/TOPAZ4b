#!/bin/bash

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../
export BASEDIR=$(pwd)/  
cd -


#find $BASEDIR -type d -name SCRATCH -exec rm -r {} \;
find $BASEDIR -type d -name SCRATCH -exec rm -r {} \;
