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


queue=8nh
sourcedir=/store/cmst3/group/hintt/LJets5TeV/
outdir=~/work/TopLJets5TeV
wwwdir=~/www/TopLJets5TeV
lumi=26.0

RED='\e[31m'
NC='\e[0m'
case $WHAT in
    SELTEST )
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_mu   --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 13;
	;;
    SEL )
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_muplus   --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 13 --charge 1  --runSysts;
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_muminus  --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 13 --charge -1 --runSysts;
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_munoniso --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 1300; 
	;;
    MERGETEST )
	./scripts/mergeOutputs.py ${outdir}/analysis_mu;
	;;
    MERGE )
	a=(muplus muminus munoniso )
	for i in ${a[@]}; do
	    ./scripts/mergeOutputs.py ${outdir}/analysis_${i};
	done
	;;
    PLOTTEST )
	python scripts/plotter.py -i ${outdir}/analysis_mu -j data/era5TeV/samples.json      -l ${lumi};	
	;;
    PLOT )
	a=(muplus muminus munoniso )
	for i in ${a[@]}; do
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/samples.json      -l ${lumi};	
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/syst_samples.json -l ${lumi} -o syst_ploter.root --silent;	
	done
	;;
    WWWTEST )
	mkdir -p ${wwwdir}/analysis_mu
        cp ${outdir}/analysis_mu/plots/*.{png,pdf} ${wwwdir}/analysis_mu
        cp test/index.php ${wwwdir}/analysis_mu
	;;
    WWW )
	a=(muplus muminus munoniso )
	for i in ${a[@]}; do
	    mkdir -p ${wwwdir}/analysis_${i}
	    cp ${outdir}/analysis_${i}/plots/*.{png,pdf} ${wwwdir}/analysis_${i}
	    cp test/index.php ${wwwdir}/analysis_${i}
	done
	;;
esac