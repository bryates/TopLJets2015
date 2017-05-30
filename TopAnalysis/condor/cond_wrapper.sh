#!/bin/sh

#determine CMSSW config
SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`
ARCH=${SCRIPTPATH##/*/}
WORKDIR=${SCRIPTPATH}/../

#configure environment
cd $WORKDIR
export SCRAM_ARCH=$ARCH
eval `scram r -sh`

#echo `pwd`
#dir=/afs/cern.ch/user/b/byates/TopAnalysis
#cd $dir
#echo `pwd`
#eval `scramv1 runtime -sh`

#cd ${_CONDOR_SCRATCH_DIR}
#echo `pwd`

#cp -R $dir/scripts .

#export CMSSW_SEARCH_PATH=${CMSSW_SEARCH_PATH}:`pwd`

echo "Executing: ${*}"

$*
#mv *.root $dir/condor
#rm *
