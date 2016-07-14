#!/bin/bash

WHAT=$1; 
UNBLIND=$2
if [ "$#" -lt 1 ]; then 
    echo "steerTOP5TeVAnalysis.sh <SEL/MERGE/BKG/PLOT/WWW/PREPAREFIT/FIT>";
    echo "        SEL          - selects data and MC";
    echo "        MERGE        - merge the output of the jobs";
    echo "        BKG          - runs the background estimation from sidebands";
    echo "        PLOT         - runs the plotter tool on the selection";
    echo "        WWW          - moves the plots to an afs-web-based area";    
    echo "        PREPAREFIT   - create datacards for the fit"
    echo "        FIT          - run the cross section fit (may need a special CMSSW release to use combine) if 1 is passed as well it will unblind"
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
	echo -e "[ ${RED} computing total number of eff. events available for analysis ${NC} ]"
	python scripts/produceNormalizationCache.py -i ${sourcedir}                -o data/era5TeV/genweights.root  --HiForest;
	echo -e "[ ${RED} projecting the b-tagging efficiency ${NC} ]"
	python scripts/saveExpectedBtagEff.py       -i ${sourcedir}/MCTTNominal_v2 -o data/era5TeV/expTageff.root   --HiForest; 
	;;
    SEL )
	echo -e "[ ${RED} Sending out jobs to batch ${NC} ]"
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_mu                                --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis       --ch 13 --runSysts --only MC;
	python scripts/runLocalAnalysis.py -i ${data}      -q ${queue} -o ${outdir}/analysis_mu/FilteredSingleMuHighPt_v3.root --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis       --ch 13;
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_munoniso                          --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis       --ch 1300 --only MC;
	python scripts/runLocalAnalysis.py -i ${data}      -q ${queue} -o ${outdir}/analysis_munoniso/FilteredSingleMuHighPt_v3.root --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 1300;
	;;
    MERGE )
	echo -e "[ ${RED} Merging job output ${NC} ]"
	a=(mu munoniso)
	for i in ${a[@]}; do
	    ./scripts/mergeOutputs.py ${outdir}/analysis_${i};
	done
	;;
    BKG )
	echo -e "[ ${RED} Running QCD estimation from non-isolated side-band ${NC} ]"
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
	echo -e "[ ${RED} Running plotter ${NC} ]"
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
	echo -e "[ ${RED} Moving plots to ${outdir} ${NC} ]"
	a=(mu munoniso)
	for i in ${a[@]}; do
	    mkdir -p ${wwwdir}/analysis_${i}
	    cp ${outdir}/analysis_${i}/plots/*.{png,pdf} ${wwwdir}/analysis_${i}
	    cp test/index.php ${wwwdir}/analysis_${i}
	done
	;;
    PREPAREFIT )
	echo -e "[ ${RED} Creating datacards ${NC} ]"
	python scripts/createDataCard.py \
	    -i ${outdir}/analysis_mu/plots/plotter.root \
	    --systInput ${outdir}/analysis_mu/plots/syst_plotter.root \
            -q ${outdir}/analysis_mu/.qcdscalefactors.pck \
	    -o ${outdir}/analysis_mu/datacard \
	    --specs TOP-16-015 \
	    --signal tbart \
            -d mjj \
	    -c 0b,1b,2b \
	    --addBinByBin 0.3 \
	    --rebin 2;
	
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
    FIT )
	echo -e "[ ${RED} $CMSSW_BASE will be used - make sure combine is compatible and is installed ${NC} ]"
	cd ${outdir}/analysis_mu/datacard;
	combineCards.py m0b=datacard_0b.dat m1b=datacard_1b.dat m2b=datacard_2b.dat > datacard.dat;
	text2workspace.py datacard.dat -m 0 -o workspace.root
        
        #expected
	commonOpts="-t -1 --expectSignal=1 --setPhysicsModelParameterRanges btag=-2,2:r=0,2 -m 0";
	combine workspace.root -M MaxLikelihoodFit ${commonOpts};
	mv mlfit.root mlfit_exp.root
	combine workspace.root -M MultiDimFit ${commonOpts} --algo=grid --points=100;
	mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_r.root
	combine workspace.root -M MultiDimFit ${commonOpts} --algo=grid --points=100 -S 0;
	mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_stat_r.root

	commonOpts="--redefineSignalPOIs btag -P btag --expectSignal=1 --algo=grid --points=100 --setPhysicsModelParameterRanges btag=-2,2:r=0,2 -m 0";	
	combine workspace.root -M MultiDimFit ${commonOpts};
	mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_btag.root
	combine workspace.root -M MultiDimFit ${commonOpts} -S 0;
	mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_stat_btag.root
	
        combine workspace.root -M MultiDimFit --algo=grid --points=2500 -m 0 -t -1 \
	    --redefineSignalPOIs r,btag -P r -P btag  --setPhysicsModelParameterRanges btag=-2,2:r=0,2 \
	    --expectSignal=1;
	mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_rvsbtag.root;

	#observed...
	if [ "${UNBLIND}" == "1" ]; then
	    echo -e "[ ${RED} will unblind the results now ${NC}]";

	    combine workspace.root -M MultiDimFit --redefineSignalPOIs btag -P btag --algo=grid --points=100 --setPhysicsModelParameterRanges btag=-2,2:r=0,2 -m 0;
	    mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_btag.root

	    combine workspace.root -M MultiDimFit -P r --algo=grid --points=100 --setPhysicsModelParameterRanges btag=-2,2:r=0,2 -m 0;
	    mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_r.root
	    
            combine workspace.root -M MultiDimFit --algo=grid --points=2500 -m 0 \
		--redefineSignalPOIs r,btag -P r -P btag  --setPhysicsModelParameterRanges btag=-2,2:r=0,2;
	    mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_rvsbtag.root;
	fi
	cd -
	;;
    SHOWFIT )
	cardsDir=${outdir}/analysis_mu/datacard;
	python scripts/fitSummaryPlots.py "#mu"=${cardsDir}/datacard.dat --POIs r,btag --label "27.9 pb^{-1} (5.02 TeV)" -o ${cardsDir};
	;;
esac