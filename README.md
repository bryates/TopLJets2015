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

## Running local analysis
