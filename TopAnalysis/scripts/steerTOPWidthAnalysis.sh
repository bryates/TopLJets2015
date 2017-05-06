#!/bin/bash

WHAT=$1; 
ERA=$2;
RUN=$3;
if [ "$#" -ne 3 ]; then 
    echo "steerTOPWidthAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/BKG/PLOT/WWW> <ERA> <RUN.";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGESEL     - merge the output of the jobs";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    echo "        ANA          - analyze the selected events";
    echo "        MERGE        - merge the output of the analysis jobs";
    echo "        BKG          - estimate DY scale factor from data";
    echo "        PLOT         - runs the plotter tool on the analysis outputs";
    echo "        WWW          - moves the analysis plots to an afs-web-based area";
    echo " "
    echo "        ERA          - era2015/era2016";
    echo "        RUN          - run range (e.g. BCDEF)";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N

queue=8nh
githash=8db9ad6
#lumi=12868.66
#lumi=35862.452 #DoubleMuon
lumi=35740.161 #SingleMuon
lumiUnc=0.062
#eosdir=/store/cmst3/user/psilva/LJets2016/${githash}
eosdir=/store/user/byates/LJets2015/${githash}
case $ERA in
    era2015)
	githash=8c1e7c9;
	lumi=2267.84;
	eosdir=/store/cmst3/user/psilva/LJets2015/${githash}
	;;
esac
lumi=`jq -r '.Data13TeV_SingleMuon["'$RUN'"]' data/era2016/lumi.json`

summaryeosdir=/store/cmst3/group/top/summer2016/TopWidth_${ERA}
#outdir=~/work/TopWidth_${ERA}
outdir=/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016
wwwdir=~/www/Top2016/2016


RED='\e[31m'
NC='\e[0m'
case $WHAT in
    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} -n 8 -q ${queue} -o ${outdir} --era ${ERA} --runPeriod ${RUN} -m TOP::RunTop --ch 0;
	;;
    MERGESEL )
	./scripts/mergeOutputs.py ${outdir} True;	
	;;
    PLOTSEL )
	#python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/${ERA}/samples.json -l ${lumi} --saveLog;# --mcUnc ${lumiUnc};	
	#python scripts/plotter.py -i ${outdir} -j data/${ERA}/samples.json -l ${lumi} --saveLog --run ${RUN};# --mcUnc ${lumiUnc};	
	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/${ERA}/samples.json -l ${lumi} --saveLog --run ${RUN};# --mcUnc ${lumiUnc};	
        #python scripts/plotter.py -i ${outdir} -j data/${ERA}/samples.json -l ${lumi};
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/test/eff_vtx/${RUN}
	mkdir -p ${wwwdir}/test/eff_vtx/${RUN}/log
	mkdir -p ${wwwdir}/test/eff_vtx/${RUN}/ee
	mkdir -p ${wwwdir}/test/eff_vtx/${RUN}/e
	mkdir -p ${wwwdir}/test/eff_vtx/${RUN}/em
	mkdir -p ${wwwdir}/test/eff_vtx/${RUN}/mumu
	mkdir -p ${wwwdir}/test/eff_vtx/${RUN}/mu
	cp -p  ${outdir}/plots/*_${RUN}.{png,pdf} ${wwwdir}/test/eff_vtx/${RUN}/
	cp -p  ${outdir}/plots/*_${RUN}_log.{png,pdf} ${wwwdir}/test/eff_vtx/${RUN}/
        rename _${RUN}. . ${wwwdir}/test/eff_vtx/${RUN}/*_${RUN}.*
        rename _${RUN}_log. _log. ${wwwdir}/test/eff_vtx/${RUN}/*_${RUN}_log.*
	rm ${wwwdir}/test/eff_vtx/${RUN}/*no_weight*_${RUN}.{png,pdf}
	mv ${wwwdir}/test/eff_vtx/${RUN}/*_log*.{png,pdf} ${wwwdir}/test/eff_vtx/${RUN}/log/
	mv ${wwwdir}/test/eff_vtx/${RUN}/*_ee_*.{png,pdf} ${wwwdir}/test/eff_vtx/${RUN}/ee/
	mv ${wwwdir}/test/eff_vtx/${RUN}/*_em_*.{png,pdf} ${wwwdir}/test/eff_vtx/${RUN}/em/
	mv ${wwwdir}/test/eff_vtx/${RUN}/*_e_*.{png,pdf} ${wwwdir}/test/eff_vtx/${RUN}/e/
	mv ${wwwdir}/test/eff_vtx/${RUN}/*_mm_*.{png,pdf} ${wwwdir}/test/eff_vtx/${RUN}/mumu/
	mv ${wwwdir}/test/eff_vtx/${RUN}/*_m_*.{png,pdf} ${wwwdir}/test/eff_vtx/${RUN}/mu/
	cp -p test/eff_vtx/index.php ${wwwdir}/test/eff_vtx/${RUN}/
	cp -p test/eff_vtx/index.php ${wwwdir}/test/eff_vtx/${RUN}/log/
	cp -p test/eff_vtx/index.php ${wwwdir}/test/eff_vtx/${RUN}/ee/
	cp -p test/eff_vtx/index.php ${wwwdir}/test/eff_vtx/${RUN}/e/
	cp -p test/eff_vtx/index.php ${wwwdir}/test/eff_vtx/${RUN}/em/
	cp -p test/eff_vtx/index.php ${wwwdir}/test/eff_vtx/${RUN}/mumu/
	cp -p test/eff_vtx/index.php ${wwwdir}/test/eff_vtx/${RUN}/mu/
	;;
    ANA )
	python scripts/runTopWidthAnalysis.py -i ${summaryeosdir} -o ${outdir}/analysis -q ${queue};
	;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir}/analysis;
	;;
    BKG )
	python scripts/plotter.py      -i ${outdir}/analysis  -j data/${ERA}/samples.json  -l ${lumi} --onlyData --only mll -o dy_plotter.root;        
	python scripts/runDYRinRout.py --in ${outdir}/analysis/plots/dy_plotter.root --categs 1b,2b --out ${outdir}/analysis/plots/;
	;;
    PLOT )
	python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json      -l ${lumi} --only count --saveTeX -o count_plotter.root --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck; 
        python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json      -l ${lumi} --onlyData --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck;        
	python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/syst_samples.json -l ${lumi} --silent -o syst_plotter.root;        
        ;;
    WWW )
        mkdir -p ${wwwdir}/ana
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana        
        cp test/index.php ${wwwdir}/ana
	;;
esac
