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
tag=$5
num=$2
period=$6
method=$7
fit=$8
syst=$9
append=""
if [ $syst == 1 ] 
then
  append="_up"
elif [ $syst == -1 ]
then
  append="_down"
fi

#echo "Running: python scripts/runLocalAnalysis.py -i /store/user/byates/LJets2015/8db9ad6/${tag}/MergedMiniEvents_${num}.root -o LJets2015/2016/${tag}_${num}.root --tag ${tag} --method TOP::RunTop --era era2016 --runPeriod ${period}"

#python scripts/runLocalAnalysis.py -i /store/user/byates/LJets2015/8db9ad6/${tag}/MergedMiniEvents_${num}.root -o LJets2015/2016/${tag}_${num}.root --tag ${tag} --method TOP::RunTop --era era2016 --runPeriod ${period}

echo "Running: python scripts/runLocalAnalysis.py -i ${in}/${tag}/MergedMiniEvents_${num}_ext.root -o ${out}/${tag}${append}_ext_${num}.root --tag ${tag} --method ${method} --era era2016 --runPeriod ${period} --rbFit ${fit} --runSysts ${syst}"

python scripts/runLocalAnalysis.py -i ${in}/${tag}/MergedMiniEvents_${num}_ext.root -o ${out}/${tag}${append}_ext_${num}.root --tag ${tag} --method ${method} --era era2016 --runPeriod ${period} --rbFit ${fit} --runSysts ${syst}
#mv *.root $dir/condor
#rm *
