#!/bin/sh

echo `pwd`
dir=/afs/cern.ch/user/b/byates/TopAnalysis
cd $dir
echo `pwd`
eval `scramv1 runtime -sh`

#cd ${_CONDOR_SCRATCH_DIR}
#echo `pwd`

#cp -R $dir/scripts .

#export CMSSW_SEARCH_PATH=${CMSSW_SEARCH_PATH}:`pwd`

echo "Executing: ${3}"

$3
#mv *.root $dir/condor
#rm *
