B1;2c# TopLJets2015

## Analysis twiki
Please keep 
https://twiki.cern.ch/twiki/bin/view/Main/TopLJ2015Analysis
up to date with the on-going tasks and results

## Installation instructions
To execute in your lxplus work area.
```
cmsrel CMSSW_7_4_14
cd CMSSW_7_4_14/src
cmsenv
git cms-merge-topic ikrav:egm_id_7.4.12_v1
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
python scripts/submitToGrid.py -j data/samples_Run2015.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --lfn my_output_directory_in_eos -s
```
Partial submission can be made adding "-o csv_list" as an option
Don't forget to init the environment for crab3
(e.g. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial)

As soon as ntuple production starts to finish, to move from crab output directories to a simpler directory structure which can be easily parsed by the local analysis run 
```
python scripts/checkProductionIntegrity.py -i /store/group/phys_top/psilva/1caa016 -o /store/cmst3/user/psilva/LJets2015/7bae03e
```
If "--cleanup" is passed, the original crab directories in EOS are removed.
For data, don't forget to create the json files with the list of runs/luminosity sections processed, e.g. as:
```
crab report grid/crab_Data13TeV_SingleElectron_2015D_v3
``` 
Then you can merge the json files for the same dataset to get the full list of run/lumi sections to analyse
```
mergeJSON.py grid/crab_Data13TeV_SingleElectron_2015D_v3/results/lumiSummary.json grid/crab_Data13TeV_SingleElectron_2015D_v4/results/lumiSummary.json --output SingleElectron_lumi.json
```
You can then run the brilcalc tool to get the integrated luminosity in total and per run.
(see https://twiki.cern.ch/twiki/bin/view/CMS/2015LumiNormtag for more details)
```
brilcalc lumi --normtag /afs/cern.ch/user/c/cmsbril/public/normtag_json/OfflineNormtagV1.json -i SingleElectron_lumi.json
```
Use the table which is printed out to update the "lumiPerRun" method in ReadTree.cc.
That will be used to monitor the event yields per run in order to identify outlier runs.
To update the pileup distributions run
```
python scripts/runPileupEstimation.py --json SingleElectron_lumi.json
```
It will store the data pileup distributions for different min.bias cross section values under data.
Before running it's better to save the expected b-tagging efficiency with
```
python scripts/saveExpectedBtagEff.py 
```
This will project the jet pT spectrum from the TTbar sample before and after applying b-tagging,
and save the expected efficiencies in data/expTageff.root.
You're now ready to start locally the analysis.

## Running locally the analysis

The analysis (histogram filling, final selection) is in src/ReadTree.cc.
Recompile (scram b) everytime you change it so that you can test the new features.
To test the code on a single file to produce plots.
```
python scripts/runLocalAnalysis.py -i MiniEvents.root
```
To run the code on a set of samples, listed in a json file you can run it as follows:
```
python scripts/runLocalAnalysis.py -i /store/cmst3/user/psilva/LJets2015/7bae03e -j data/samples_Run2015.json -n 8 -o analysis_muplus   --ch 13   --charge 1
python scripts/runLocalAnalysis.py -i /store/cmst3/user/psilva/LJets2015/7bae03e -j data/samples_Run2015.json -n 8 -o analysis_muminus  --ch 13   --charge -1
python scripts/runLocalAnalysis.py -i /store/cmst3/user/psilva/LJets2015/7bae03e -j data/samples_Run2015.json -n 8 -o analysis_munoniso --ch 1300
python scripts/runLocalAnalysis.py -i /store/cmst3/user/psilva/LJets2015/7bae03e -j data/syst_samples_Run2015.json -n 8 -o analysis_muplus   --ch 13   --charge 1
python scripts/runLocalAnalysis.py -i /store/cmst3/user/psilva/LJets2015/7bae03e -j data/syst_samples_Run2015.json -n 8 -o analysis_muminus  --ch 13   --charge -1
```
The first time it runs over the directory it will compute the normalization factor for MC
such that the distributions will correspond to 1/pb of data.
The normalization factor is given by (xsec / N generated events)
where xsec is stored in the json file, and N generated events is summed up
from the "counter" histogram stored in the the files to process for each process.
The first time it also computes the pileup weights on a sample-by-sample basis
by taking the ratio of the of the putrue distribution to the pileup distribution estimated in data.
Both the normalization factors and the pileup weights are stored under the "analysis" directory
in a cache file called ".xsecweights.pck".
To plot the output of the local analysis you can run the following:
```
python scripts/plotter.py -i analysis_muplus/ -j data/samples_Run2015.json -l 1262
python scripts/plotter.py -i analysis_muminus/ -j data/samples_Run2015.json -l 1262
python scripts/plotter.py -i analysis_munoniso/ -j data/samples_Run2015.json -l 1262
```
After the plotters are created one can run the QCD estimation normalization, by fitting the MET distribution.
The script will also produce the QCD templates using the data from the sideband region. It runs as
```
python scripts/runQCDEstimation.py --iso analysis_muplus/plots/plotter.root --noniso analysis_munoniso/plots/plotter.root --out analysis_muplus/ --norm metpt
python scripts/runQCDEstimation.py --iso analysis_muminus/plots/plotter.root --noniso analysis_munoniso/plots/plotter.root --out analysis_muminus/ --norm metpt
```
The output is a ROOT file called Data_QCDMultijets.root which can now be used in addition to the predictions of all the other backgrounds.
To include it in the final plots you can run the plotter script again (see instructions above).

## Cross section fitting
We use the Higgs combination tool to perform the fit of the production cross section.
(cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit for details of the release to use).
To create the datacard you can run the following script
```
python scripts/createDataCard.py -i analysis_muplus/plots/plotter.root -o analysis_muplus/datacard -d njetsnbtags
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
