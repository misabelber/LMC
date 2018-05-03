#! /bin/bash

HERE=`pwd`

#PLEASE, provide manually the following PATHS:

DATA_PATH='/scratch/bernardos/LMC/' #Folder to store heavy data files (Observations and datacubes)
CTOOLS_PATH="/scratch/bernardos/ctools/" #Folder where ctools is installed ("gamma folder must be contained here"

var=`grep -r -h 'DATA_PATH = ' $HERE/pipelines/pipes_in_py/ | head -1`
OLD_DATA_PATH=`echo $var | awk '{print $3}' | sed -e 's/^.//' -e 's/.$//'`

find $HERE -type f | xargs sed -i "s%${OLD_DATA_PATH}%${DATA_PATH}%g"

var=`grep -r -h 'CTOOLS_PATH = ' $HERE/pipelines/pipes_in_py/ | head -1`
OLD_CTOOLS_PATH=`echo $var | awk '{print $3}' | sed -e 's/^.//' -e 's/.$//'`

find $HERE -type f | xargs sed -i "s%${OLD_CTOOLS_PATH}%${CTOOLS_PATH}%g"


#Change the path for LMC folder 

LMC_PATH=`pwd` #THIS FOLDER IS LMC FOLDER

var=`grep -r -h 'LMC_PATH = ' $LMC_PATH/tools/tools_in_py/ | head -1`
#tr -d $var 'LMC_PATH'
OLD_LMC_PATH=`echo $var | awk '{print $3}' | sed -e 's/^.//' -e 's/.$//'`
OLD_LMC_PATH='/afs/ciemat.es/user/b/bernardos/Documentos/LMC'

find $HERE -type f | xargs sed -i "s%${OLD_LMC_PATH}%${LMC_PATH}%g"


