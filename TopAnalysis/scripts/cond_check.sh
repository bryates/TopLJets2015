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

in=$3
out=$4
pub=$5

#echo "Running: python scripts/runLocalAnalysis.py -i /store/user/byates/LJets2015/8db9ad6/${tag}/MergedMiniEvents_${num}.root -o LJets2015/2016/${tag}_${num}.root --tag ${tag} --method TOP::RunTop --era era2016 --runPeriod ${period}"

#python scripts/runLocalAnalysis.py -i /store/user/byates/LJets2015/8db9ad6/${tag}/MergedMiniEvents_${num}.root -o LJets2015/2016/${tag}_${num}.root --tag ${tag} --method TOP::RunTop --era era2016 --runPeriod ${period}

echo "python scripts/checkProductionIntegrity.py -i ${in} -o ${out} --nocheck --only ${pub}"

python scripts/checkProductionIntegrity.py -i ${in} -o ${out} --nocheck --only ${pub}
#mv *.root $dir/condor
#rm *
