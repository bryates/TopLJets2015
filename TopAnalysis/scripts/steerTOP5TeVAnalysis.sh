#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerTOP5TeVAnalysis.sh <SEL/MERGE/BKG/PLOT/WWW>";
    echo "        SEL          - selects data and MC";
    echo "        MERGE        - merge the output of the jobs";
    echo "        BKG          - runs the background estimation from sidebands";
    echo "        PLOT         - runs the plotter tool on the selection";
    echo "        WWW          - moves the plots to an afs-web-based area";    
    echo "        FIT          - run the cross section fit"
    exit 1; 
fi

#suppress batch notifications by mail
export LSB_JOB_REPORT_MAIL=N

queue=8nh
sourcedir=/store/cmst3/group/hintt/LJets5TeV/
outdir=~/work/LJets-5TeV
wwwdir=~/www/LJets-5TeV
lumi=27.9
data=/store/cmst3/group/hintt/mverweij/PP5TeV/data/SingleMuHighPt/crab_FilteredSingleMuHighPt_v3/160425_163333/merge/HiForest_0.root

RED='\e[31m'
NC='\e[0m'
case $WHAT in
    PREPARE )
	python scripts/produceNormalizationCache.py -i ${sourcedir}                -o data/era5TeV/genweights.root  --HiForest;
	python scripts/saveExpectedBtagEff.py       -i ${sourcedir}/MCTTNominal_v2 -o data/era5TeV/expTageff.root   --HiForest; 
	;;
    SEL )
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_mu                                --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 13 --runSysts;
	python scripts/runLocalAnalysis.py -i ${data}      -q ${queue} -o ${outdir}/analysis_mu/FilteredSingleMuHighPt_v3.root --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 13;
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_munoniso                          --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 1300;
	python scripts/runLocalAnalysis.py -i ${data}      -q ${queue} -o ${outdir}/analysis_munoniso/FilteredSingleMuHighPt_v3.root --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 1300;
	;;
    MERGE )
	a=(mu munoniso)
	for i in ${a[@]}; do
	    ./scripts/mergeOutputs.py ${outdir}/analysis_${i};
	done
	;;
    BKG )
	a=(mu munoniso)
        for i in ${a[@]}; do
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/samples.json      -l ${lumi} --silent;
	done
	python scripts/runQCDEstimation.py \
	    --iso    ${outdir}/analysis_mu/plots/plotter.root \
	    --noniso ${outdir}/analysis_munoniso/plots/plotter.root \
	    --out    ${outdir}/analysis_mu/ \
	    --sels  ,0b,1b,2b;
	;;
    PLOT )
	a=(mu munoniso)
	for i in ${a[@]}; do
	    #python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/Wsamples.json     -l ${lumi} --saveLog --noStack;	
	    #mkdir ~/${outdir}/analysis_${i}/wplots;
	    #mv ~/${outdir}/analysis_${i}/plots/* ~/${outdir}/analysis_${i}/wplots/;
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/samples.json      -l ${lumi} --saveLog;	
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/syst_samples.json -l ${lumi} -o syst_plotter.root --silent;	
	done
	;;

    WWW )
	a=(mu munoniso)
	for i in ${a[@]}; do
	    mkdir -p ${wwwdir}/analysis_${i}
	    cp ${outdir}/analysis_${i}/plots/*.{png,pdf} ${wwwdir}/analysis_${i}
	    cp test/index.php ${wwwdir}/analysis_${i}
	done
	;;
    FIT )
	echo -e "[ ${RED} Creating datacards ${NC} ]"
	python scripts/createDataCard.py \
	    -i ${outdir}/analysis_mu/plots/plotter.root \
	    --systInput ${outdir}/analysis_mu/plots/syst_plotter.root \
            -q ${outdir}/analysis_mu/.qcdscalefactors.pck \
	    -o ${outdir}/analysis_mu/datacard \
	    --specs TOP-16-015 \
	    --signal tbart \
            -d mjj \
	    -c 0b,1b,2b;
	
	a=(0b 1b 2b)
	for i in ${a[@]}; do
	    python scripts/projectShapeUncs.py ${outdir}/analysis_mu/datacard/shapes_${i}.root btag,othertag,jes,jer;
	    python scripts/projectShapeUncs.py ${outdir}/analysis_mu/datacard/shapes_${i}.root ttPartonShower,Hadronizer,ttFactScale,ttRenScale,ttCombScale;
	    python scripts/projectShapeUncs.py ${outdir}/analysis_mu/datacard/shapes_${i}.root wFactScale,wRenScale,wCombScale W;
	done
	mkdir -p ${wwwdir}/shapes
	mv *.{png,pdf} ${wwwdir}/shapes;
	cp test/index.php ${wwwdir}/shapes;
	;;
esac