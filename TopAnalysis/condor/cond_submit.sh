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

tag=$3
num=$2
period=$4

echo "Running: python scripts/runLocalAnalysis.py -i /store/user/byates/LJets2015/8db9ad6/${tag}/MergedMiniEvents_${num}.root -o LJets2015/2016/${tag}_${num}.root --tag ${tag} --method TOP::RunTop --era era2016 --runPeriod ${period}"

python scripts/runLocalAnalysis.py -i /store/user/byates/LJets2015/8db9ad6/${tag}/MergedMiniEvents_${num}.root -o LJets2015/2016/${tag}_${num}.root --tag ${tag} --method TOP::RunTop --era era2016 --runPeriod ${period}
#mv *.root $dir/condor
#rm *
