# TopLJets2015

## Installation instructions
To execute in your lxplus work area.
```
cmsrel CMSSW_7_4_14
cd CMSSW_7_4_14/src
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git checkout 74x-root6
scram b -j 9
cd ..
git clone git@github.com:pfs/TopLJets2015.git
scram b -j 9
```

## Running ntuple creation
To run locally the ntuplizer, for testing purposes
```
cmsRun test/runMiniAnalyzer_cfg.py runOnData=False/True outFilename=MiniEvents.root
```
To submit a list of samples, described in a json file to the grid you can use the following script.
```
python submitToGrid.py -j data/samples_Run2015.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --lfn my_output_directory_in_eos -s
```
Partial submission can be made adding "-o csv_list" as an option
Don't forget to init the environment for crab3
(e.g. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial)

As soon as ntuple production starts to finish, to move from crab output directories to a simpler directory structure which can be easily parsed by the local analysis run 
```
python checkProductionIntegrity.py -i my_output_directory_in_eos
```
If "--cleanup" is passed, the original crab directories in EOS are removed.

## Running local analysis
