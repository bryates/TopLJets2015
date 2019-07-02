# TopLJets2015

## Analysis twiki
Please keep 
https://twiki.cern.ch/twiki/bin/view/Main/TopLJ2016Analysis
up to date with the on-going tasks and results

## Installation instructions
To execute in your lxplus work area.
```
cmsrel CMSSW_8_0_26
cd CMSSW_8_0_26/src
cmsenv
git clone -b 80x_rereco git@github.com:bryates/TopLJets2015.git
#For BFragmentationAnalyzer
mkdir TopQuarkAnalysis
cd TopQuarkAnalysis
git clone https://git@gitlab.cern.ch:8443/CMS-TOPPAG/BFragmentationAnalyzer.git (Kereros 5 on lxplus)
#git clone ssh://git@gitlab.cern.ch:7999/CMS-TOPPAG/BFragmentationAnalyzer.git (ssh)
cd -
scram b -j8
```

## Running ntuple creation
First time create a symbolic link to the jet energy corrections files
```
ln -s data/era2016/Spring16_25nsV3_DATA.db
ln -s data/era2016/Spring16_25nsV3_MC.db
ln -s data/era2016/RoccoR_13tev.txt 
```
To run locally the ntuplizer, for testing purposes
```
cmsRun test/runMiniAnalyzer_cfg.py runOnData=False/True outFilename=MiniEvents.root
```
To submit a list of samples, described in a json file to the grid you can use the following script.
```
python scripts/submitToGrid.py -j data/era2016/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --lfn /store/group/phys_top/byates -s
python scripts/submitToGrid.py -j data/era2016/syst_samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --lfn /store/group/phys_top/byates -s
```
Partial submission can be made adding "-o csv_list" as an option
Don't forget to init the environment for crab3 (e.g. https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3CheatSheet#Environment_setup)
```
source /cvmfs/cms.cern.ch/crab3/crab.sh
```
As soon as ntuple production starts to finish, to move from crab output directories to a simpler directory structure which can be easily parsed by the local analysis runThe merging can be run locally if needed by using the checkProductionIntegrity.py script
```
python scripts/submitCheckProductionIntegrity.py -i /store/group/phys_top/byates/92cadc6 -o /store/user/byates/LJets2016/8db9ad6
```

## Preparing the analysis 

Correction and uncertainty files are stored under data by era directories (e.g. data/era2015, data/era2016) in order no to mix different periods.
After ntuples are processed start by creating the json files with the list of runs/luminosity sections processed, e.g. as:
```
for i in $(find grid/ -maxdepth 1 | grep crab_Data);
do
    crab report ${i}; 
done
``` 
Then you can merge the json files for the same dataset to get the full list of run/lumi sections to analyse
```
mergeJSON.py grid/crab_Data13TeV_SingleMuon_2016B/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016C/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016D/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016E/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016F/results/processedLumis.json --output data/era2016/Data13TeV_SingleMuon_lumis_BCDEF.json
mergeJSON.py grid/crab_Data13TeV_SingleMuon_2016G/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016H_v2/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016H_v3/results/processedLumis.json --output data/era2016/Data13TeV_SingleMuon_lumis_GH.json

mergeJSON.py grid/crab_Data13TeV_SingleMuon_2016B/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016C/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016D/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016E/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016F/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016G/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016H_v2/results/processedLumis.json grid/crab_Data13TeV_SingleMuon_2016H_v3/results/processedLumis.json --output data/era2016/Data13TeV_SingleMuon_lumis.json
```
You can then run the brilcalc tool to get the integrated luminosity in total and per run (see https://twiki.cern.ch/twiki/bin/view/CMS/2015LumiNormtag for more details).
```
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH
brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/Normtags/normtag_DATACERT.json -u /pb -i data/era2016/Data13TeV_SingleMuon_lumis.json --hltpath="HLT_Iso*24_v*"
```
Use the table which is printed out to update the "lumiPerRun" method in ReadTree.cc.
That will be used to monitor the event yields per run in order to identify outlier runs.
* Pileup weighting. To update the pileup distributions run the script below. It will store the data pileup distributions for different min.bias cross section in data/pileupWgts.root
```
python scripts/runPileupEstimation.py --json data/era2016/Data13TeV_SingleMuon_lumis_BCDEF.json --out data/era2016/pileupWgtsBCDEF.root
python scripts/runPileupEstimation.py --json data/era2016/Data13TeV_SingleMuon_lumis_GH.json --out data/era2016/pileupWgtsGH.root
```
* B-tagging. To apply corrections to the simulation one needs the expected efficiencies stored somwewhere. The script below will project the jet pT spectrum from the TTbar sample before and after applying b-tagging, to compute the expecte efficiencies. The result will be stored in data/expTageff.root
```
python scripts/saveExpectedBtagEff.py -i /store/group/phys_top/byates/LJets2016//8db9ad6/MC13TeV_TTJets_powheg -o data/era2016/expTageff.root
```
* MC normalization. This will loop over all the samples available in EOS and produce a normalization cache (weights to normalize MC). The file will be available in data/genweights.root and data/genweights_syst.root
```
python scripts/produceNormalizationCache.py -i /store/group/phys_top/byates/LJets2016/8db9ad6 -o data/era2016/genweights.root
python scripts/produceNormalizationCache.py -i /store/group/phys_top/byates/syst_samples/ -o data/era2016/genweights_syst.root -j data/era2016/syst_samples.json
```
You're now ready to start locally the analysis.


## Running locally the analysis for testing

The analysis (histogram filling, final selection) is in src/ReadTree.cc.
Recompile (scram b) everytime you change it so that you can test the new features.
To test the code on a single file to produce plots.
```
python scripts/runLocalAnalysis.py -i MiniEvents.root
```
To run the code on a set of samples stored in EOS you can run it as shown below.
If "-q queue_name" is appended the jobs are submitted to the batch system instead of running locally. 
To check the status of your jobs run "bjobs" and then "bpeek job_number" if you want to inspect how the job is running in the cluster.
If "-n n_jobs" is passed the script runs locally using "n_jobs" parallel threads.
```
python scripts/runLocalAnalysis.py -i /store/group/phys_top/byates/LJets2016/8db9ad6/ -o LJets2015/2016/ --method TOP::RunTopKalman --era era2016 --runPeriod BCDEFGH
```
After the jobs have run you can merge the outputs with
```
python scripts/mergeOutputs.py LJets2016/8db9ad6/ True
```
The True flag merges the histograms only.
To plot the output of the local analysis you can run the following:
```
python scripts/plotter.py -i LJets2015/2016/ --puNormSF puwgtctr -j data/era2016/samples.json -l data/era2016/lumi.json --run BCDEFGH
```

## Submitting the full analysis to the batch system

A script wraps up the above procedure for all the signal and control regions used in the analyis.
To use it you can use the following script
```
sh scripts/steerTOPWidthAnalysis.sh <DISTS/MERGE/PLOT/BKG>
```

## Updating the code

Commit your changes regularly with
```
git commit -a -m'comment on the changes made'
```
Push to your forked repository
```
git push git@github.com:MYGITHUBLOGIN/TopLJets2015.git
```
From the github area of the repository cleak on the green button "Compare,review and create a pull request"
to create the PR to merge with your colleagues.
