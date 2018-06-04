#! /bin/bash

HERE=`pwd`
#=================================================================================================
#PLEASE, provide MANUALLY the following PATHS:
#=================================================================================================

#PATHS FOR PYTHON/CTOOLS PIPELINES

DATA_PATH='/home/queenmab/DATA/LMC/LMC/' #Folder to store heavy data files (Observations and datacubes)
CTOOLS_PATH="/home/queenmab/anaconda3/envs/ctools/" #Folder where ctools is installed ("gamma folder must be contained here"

#PATHS FOR C++/ROOT PIPES
#This tools require the Math repository

MATH_PATH='/home/queenmab/GitHub/Math/' #Define the path of your Math repository

======================================================================

var=`grep -r -h '/Math/' $HERE/pipelines/pipes_in_C/ | head -1`
var=`echo $var | awk '{print $2}' | sed -e 's/^.//' -e 's/.$//'`
OLD_MATH_PATH=`echo "${var%/*}"`

find $HERE -type f | xargs sed -i "s%${OLD_MATH_PATH}%${MATH_PATH}%g"

var=`grep -r -h 'DATA_PATH = ' $HERE/pipelines/pipes_in_py/ | head -1`
OLD_DATA_PATH=`echo $var | awk '{print $3}' | sed -e 's/^.//' -e 's/.$//'`

find $HERE -type f | xargs sed -i "s%${OLD_DATA_PATH}%${DATA_PATH}%g"

var=`grep -r -h 'CTOOLS_PATH = ' $HERE/pipelines/pipes_in_py/ | head -1`
OLD_CTOOLS_PATH=`echo $var | awk '{print $3}' | sed -e 's/^.//' -e 's/.$//'`

find $HERE -type f | xargs sed -i "s%${OLD_CTOOLS_PATH}%${CTOOLS_PATH}%g"


#=================================================================================================
#Change the path for LMC folder (THIS IS AUTOMATC, NO NEED TO CHANGE MANUALLY) 

LMC_PATH=`pwd` #THIS FOLDER IS LMC FOLDER

var=`grep -r -h 'LMC_PATH = ' $LMC_PATH/tools/tools_in_py/ | head -1`
#tr -d $var 'LMC_PATH'
OLD_LMC_PATH=`echo $var | awk '{print $3}' | sed -e 's/^.//' -e 's/.$//'`
#OLD_LMC_PATH='/home/queenmab/GitHub/LMC'

find $HERE -type f | xargs sed -i "s%${OLD_LMC_PATH}%${LMC_PATH}%g"


