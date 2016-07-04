#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerTOP5TeVAnalysis.sh <SEL/MERGE/PLOT/WWW>";
    echo "        SEL          - selects data and MC";
    echo "        MERGE        - merge the output of the jobs";
    echo "        PLOT         - runs the plotter tool on the selection";
    echo "        WWW          - moves the plots to an afs-web-based area";    
    echo "        The previous step can be called with *TEST to run a simplified version of the analysis without charge selection/systematics"
    exit 1; 
fi

#suppress batch notifications by mail
export LSB_JOB_REPORT_MAIL=N

queue=8nh
sourcedir=/store/cmst3/group/hintt/LJets5TeV/
outdir=~/work/LJets-5TeV
wwwdir=~/www/LJets-5TeV
lumi=26.0
data=/store/cmst3/group/hintt/mverweij/PP5TeV/data/SingleMuHighPt/crab_FilteredSingleMuHighPt_v3/160425_163333/merge/HiForest_0.root

RED='\e[31m'
NC='\e[0m'
case $WHAT in
    PREPARE )
	python scripts/produceNormalizationCache.py -i ${sourcedir}                -o data/era5TeV/genweights.root  --HiForest;
	python scripts/saveExpectedBtagEff.py       -i ${sourcedir}/MCTTNominal_v2 -o data/era5TeV/expTageff.root   --HiForest; 
	;;
    SELTEST )
	queue=local
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_mutest   --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 13 --only Pythia8;
	#python scripts/runLocalAnalysis.py -i ${data}      -q ${queue} -o ${outdir}/analysis_mutest/FilteredSingleMuHighPt_v3.root --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 13;
	;;
    SEL )
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_mu   --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 13 --runSysts;
	python scripts/runLocalAnalysis.py -i ${data}      -q ${queue} -o ${outdir}/analysis_mu/FilteredSingleMuHighPt_v3.root --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 13;
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_munoniso --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 1300; 
	python scripts/runLocalAnalysis.py -i ${data}      -q ${queue} -o ${outdir}/analysis_munoniso/FilteredSingleMuHighPt_v3.root --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 1300;
	;;
    MERGETEST )
	./scripts/mergeOutputs.py ${outdir}/analysis_mutest;
	;;
    MERGE )
	a=(mu munoniso)
	for i in ${a[@]}; do
	    ./scripts/mergeOutputs.py ${outdir}/analysis_${i};
	done
	;;
    PLOTTEST )
	#python scripts/plotter.py -i ${outdir}/analysis_mutest -j data/era5TeV/samples.json      -l ${lumi} --saveLog;	
	python scripts/plotter.py -i ${outdir}/analysis_mutest -j data/era5TeV/Wsamples.json      -l ${lumi} --saveLog --noStack;	
	;;
    PLOT )
	a=(mu munoniso)
	for i in ${a[@]}; do
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/samples.json      -l ${lumi};	
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/syst_samples.json -l ${lumi} -o syst_ploter.root --silent;	
	done
	;;
    WWWTEST )
	mkdir -p ${wwwdir}/analysis_mutest
        cp ${outdir}/analysis_mutest/plots/*.{png,pdf} ${wwwdir}/analysis_mutest
        cp test/index.php ${wwwdir}/analysis_mutest
	;;
    WWW )
	a=(mu munoniso)
	for i in ${a[@]}; do
	    mkdir -p ${wwwdir}/analysis_${i}
	    cp ${outdir}/analysis_${i}/plots/*.{png,pdf} ${wwwdir}/analysis_${i}
	    cp test/index.php ${wwwdir}/analysis_${i}
	done
	;;
esac