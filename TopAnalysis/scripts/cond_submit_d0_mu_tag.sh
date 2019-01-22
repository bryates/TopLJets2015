#!/bin/sh

echo `pwd`
dir=/afs/cern.ch/user/b/byates/TopAnalysis
cd $dir
echo `pwd`
eval `scramv1 runtime -sh`


cd LJets2015/2016/mtop/
#cd ${_CONDOR_SCRATCH_DIR}
#echo `pwd`

#cp -R $dir/scripts .

#export CMSSW_SEARCH_PATH=${CMSSW_SEARCH_PATH}:`pwd`

name=$3
isData=$4
epoch=$5
frag=$6

echo "root -q -b '/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/mtop/splot_d0_mu_tag.C(\"${name}\",${isData},\"${frag}\",${epoch})'"

root -q -b '/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/mtop/splot_d0_mu_tag.C("'${name}'",'${isData},'"'${frag}'",'${epoch}')'
