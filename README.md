# TopLJets2015

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
First time create a symbolic link to the jet energy corrections files
```
ln -s data/Summer15_25nsV6_DATA.db
ln -s data/Summer15_25nsV6_MC.db
```
To run locally the ntuplizer, for testing purposes
```
cmsRun test/runMiniAnalyzer_cfg.py runOnData=False/True outFilename=MiniEvents.root
```
To submit a list of samples, described in a json file to the grid you can use the following script.
```
python scripts/submitToGrid.py -j data/samples_Run2015.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --lfn my_output_directory_in_eos -s
```
Partial submission can be made adding "-o csv_list" as an option
Don't forget to init the environment for crab3 (e.g. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial).

As soon as ntuple production starts to finish, to move from crab output directories to a simpler directory structure which can be easily parsed by the local analysis run 
```
python scripts/checkProductionIntegrity.py -i /store/group/phys_top/psilva/b18c191 -o /store/cmst3/user/psilva/LJets2015/b18c191
```
If "--cleanup" is passed, the original crab directories in EOS are removed.

## Preparing the analysis 

After ntuples are processed start by creating the json files with the list of runs/luminosity sections processed, e.g. as:
```
crab report grid/crab_Data13TeV_SingleElectron_2015D_v3
``` 
Then you can merge the json files for the same dataset to get the full list of run/lumi sections to analyse
```
mergeJSON.py grid/crab_Data13TeV_SingleElectron_2015D_v3/results/lumiSummary.json grid/crab_Data13TeV_SingleElectron_2015D_v4/results/lumiSummary.json --output data/SingleElectron_lumiSummary.json
```
You can then run the brilcalc tool to get the integrated luminosity in total and per run (see https://twiki.cern.ch/twiki/bin/view/CMS/2015LumiNormtag for more details).
```
brilcalc lumi --normtag /afs/cern.ch/user/c/cmsbril/public/normtag_json/OfflineNormtagV1.json -i data/SingleElectron_lumiSummary.json
```
Use the table which is printed out to update the "lumiPerRun" method in ReadTree.cc.
That will be used to monitor the event yields per run in order to identify outlier runs.
* Pileup weighting. To update the pileup distributions run the script below. It will store the data pileup distributions for different min.bias cross section in data/pileupWgts.root
```
python scripts/runPileupEstimation.py --json data/SingleElectron_lumiSummary.json
```
* B-tagging. To apply corrections to the simulation one needs the expected efficiencies stored somwewhere. The script below will project the jet pT spectrum from the TTbar sample before and after applying b-tagging, to compute the expecte efficiencies. The result will be stored in data/expTageff.root
```
python scripts/saveExpectedBtagEff.py 
```
* MC normalization. This will loop over all the samples available in EOS and produce a normalization cache (weights to normalize MC). The file will be available in data/genweights.pck
```
python scripts/produceNormalizationCache.py -i /store/cmst3/user/psilva/LJets2015/5736a2c
```
You're now ready to start locally the analysis.


## Running locally the analysis

The analysis (histogram filling, final selection) is in src/ReadTree.cc.
Recompile (scram b) everytime you change it so that you can test the new features.
To test the code on a single file to produce plots.
```
python scripts/runLocalAnalysis.py -i MiniEvents.root
```
To run the code on a set of samples stored in EOS you can run it as follows:
```
python scripts/runLocalAnalysis.py -i /store/cmst3/user/psilva/LJets2015/5736a2c -q 8nh --runSysts -o analysis_muplus   --ch 13   --charge 1
python scripts/runLocalAnalysis.py -i /store/cmst3/user/psilva/LJets2015/5736a2c -q 8nh --runSysts -o analysis_muminus  --ch 13   --charge -1
python scripts/runLocalAnalysis.py -i /store/cmst3/user/psilva/LJets2015/5736a2c -q 8nh            -o analysis_munoniso --ch 1300
```
If "-q queue_name" is appended the jobs are submitted to the batch system instead of running locally. 
To check the status of your jobs run "bjobs" and then "bpeek job_number" if you want to inspect how the job is running in the cluster.
If instead "-n n_jobs" is passed the script runs locally using "n_jobs" parallel threads.
After the jobs have run you can merge the outputs with
```
./scripts/mergeOutputs.py analysis_muplus
./scripts/mergeOutputs.py analysis_muminus
./scripts/mergeOutputs.py analysis_munoniso/
```
To plot the output of the local analysis you can run the following:
```
python scripts/plotter.py -i analysis_muplus/   -j data/samples_Run2015.json                           -l 2093.6
python scripts/plotter.py -i analysis_muminus/  -j data/samples_Run2015.json                           -l 2093.6
python scripts/plotter.py -i analysis_muplus/   -j data/syst_samples_Run2015.json -o syst_plotter.root -l 2093.6
python scripts/plotter.py -i analysis_muminus/  -j data/syst_samples_Run2015.json -o syst_plotter.root -l 2093.6
python scripts/plotter.py -i analysis_munoniso/ -j data/samples_Run2015.json                           -l 2093.6
```
After the plotters are created one can run the QCD estimation normalization, by fitting the MET distribution.
The script will also produce the QCD templates using the data from the sideband region. It runs as
```
python scripts/runQCDEstimation.py --iso analysis_muplus/plots/plotter.root --noniso analysis_munoniso/plots/plotter.root --out analysis_muplus/
python scripts/runQCDEstimation.py --iso analysis_muminus/plots/plotter.root --noniso analysis_munoniso/plots/plotter.root --out analysis_muminus/
```
The output is a ROOT file called Data_QCDMultijets.root which can now be used in addition to the predictions of all the other backgrounds.
To include it in the final plots you can run the plotter script again (see instructions above).

## Cross section fitting
We use the Higgs combination tool to perform the fit of the production cross section.
(cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit for details of the release to use).
It currently has to be run from a CMSSW_7_1_5 release.

### Cut in categories

To create the datacard you can run the following script
```
python scripts/createDataCard.py -i analysis_muplus/plots/plotter.root -o  analysis_muplus/datacard  -q analysis_muplus/.qcdscalefactors.pck -d nbtags
python scripts/createDataCard.py -i analysis_muminus/plots/plotter.root -o analysis_muminus/datacard -q analysis_muminus/.qcdscalefactors.pck -d nbtags
```
Combine the datacards per category above into the final one per channel
```
#charge combinations
a=(muplus muminus)
for i in ${a[@]}; do
    cd analysis_${i}/datacard;
    combineCards.py ${i}1j=datacard_1j.dat ${i}2j=datacard_2j.dat ${i}3j=datacard_3j.dat ${i}4j=datacard_4j.dat > datacard.dat
    cd -;
done

#final combination
mkdir -p analysis_mu/datacard
cd analysis_mu/datacard
combineCards.py muplus=../../analysis_muplus/datacard/datacard.dat muminus=../../analysis_muminus/datacard/datacard.dat > datacard.dat
cd -
```
Run the fits and show the results
```
python scripts/fitCrossSection.py "#mu^{+}"=analysis_muplus/datacard/datacard.dat -o analysis_muplus/datacard &
python scripts/fitCrossSection.py "#mu^{-}"=analysis_muminus/datacard/datacard.dat -o analysis_muminus/datacard &
python scripts/fitCrossSection.py "#mu^{#pm}"=analysis_mu/datacard/datacard.dat -o analysis_mu/datacard &
```
After all is run you can also compare the results with
```
python scripts/fitCrossSection.py "#mu^{+}"=analysis_muplus/datacard/datacard.dat  "#mu^{-}"=analysis_muminus/datacard/datacard.dat  "#mu^{#pm}"=analysis_mu/datacard/datacard.dat --noFit
```

### Full shape analysis

To create the datacard you can run the following script
```
a=(muplus muminus)
for i in ${a[@]}; do 
    python scripts/createDataCard.py -i analysis_${i}/plots/plotter.root --systInput analysis_${i}/plots/syst_plotter.root -o  analysis_${i}/datacard_shape  -q analysis_${i}/.qcdscalefactors.pck -d mt     -c 1j0t,2j0t,3j0t,4j0t;
    python scripts/createDataCard.py -i analysis_${i}/plots/plotter.root --systInput analysis_${i}/plots/syst_plotter.root -o  analysis_${i}/datacard_shape  -q analysis_${i}/.qcdscalefactors.pck -d minmlb -c 1j1t,2j1t,2j2t,3j1t,3j2t,4j1t,4j2t;
done
```
Combine the datacards per category above into the final one per channel
```
#charge combinations
a=(muplus muminus)
cats=(1j0t 1j1t 2j0t 2j1t 2j2t 3j0t 3j1t 3j2t 4j0t 4j1t 4j2t)
for i in ${a[@]}; do
    cd analysis_${i}/datacard_shape;
    tocombine=""
    for c in ${cats[@]}; do
    	tocombine="${i}${c}=datacard_${c}.dat ${tocombine}"
    done
    combineCards.py ${tocombine} > datacard.dat
    cd -;
done

#final combination
mkdir -p analysis_mu/datacard_shape
cd analysis_mu/datacard_shape
combineCards.py muplus=../../analysis_muplus/datacard_shape/datacard.dat muminus=../../analysis_muminus/datacard_shape/datacard.dat > datacard.dat
cd -
```
Run the fits and show the results
```
python scripts/fitCrossSection.py "#mu^{+}"=analysis_muplus/datacard_shape/datacard.dat  -o analysis_muplus/datacard_shape --POIs r,Mtop &
python scripts/fitCrossSection.py "#mu^{-}"=analysis_muminus/datacard_shape/datacard.dat -o analysis_muminus/datacard_shape --POIs r,Mtop &
python scripts/fitCrossSection.py "#mu^{#pm}"=analysis_mu/datacard_shape/datacard.dat    -o analysis_mu/datacard_shape  --POIs r,Mtop &
```
After all is run you can also compare the results with
```
python scripts/fitCrossSection.py "#mu^{+}"=analysis_muplus/datacard_shape/datacard.dat  "#mu^{-}"=analysis_muminus/datacard_shape/datacard.dat  "#mu^{#pm}"=analysis_mu/datacard_shape/datacard.dat --noFit
python scripts/fitCrossSection.py "#mu(c&c)"=analysis_muplus/datacard/datacard.dat  "#mu(shape)"=analysis_mu/datacard_shape/datacard.dat --noFit
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
