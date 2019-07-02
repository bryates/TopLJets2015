#!/bin/bash

WHAT=$1; 
ERA=$2;
RUN=$3;
if [ "$#" -lt 3 ]; then 
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
lumiFile=/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/data/era2016/lumi.json
#eosdir=/store/cmst3/user/psilva/LJets2016/${githash}
eosdir=/store/user/byates/LJets2015/${githash}
extdir=/store/group/phys_top/byates/ext
case $ERA in
    era2015)
	githash=8c1e7c9;
	lumi=2267.84;
	eosdir=/store/cmst3/user/psilva/LJets2015/${githash}
	;;
esac
lumi=`jq -r '.Data13TeV_SingleMuon["'$RUN'"]' ${lumiFile}`

summaryeosdir=/store/cmst3/group/top/summer2016/TopWidth_${ERA}
#outdir=~/work/TopWidth_${ERA}
outdir=/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/etaPiK/
subdir='/test'
if [ ! -z $4 ]; then
    subdir=$4
fi
wwwdir=~/www/Top2016/2016/$subdir
echo $wwwdir


RED='\e[31m'
NC='\e[0m'
case $WHAT in
    SEL )
        if [ "${RUN}" == "BCDEFGH" ]; then
	  python scripts/runLocalAnalysis.py -i ${eosdir} -n 8 -q ${queue} -o ${outdir} --era ${ERA} --runPeriod ${RUN} --dataOnly -m TOP::RunTop --ch 0;
	  python scripts/runLocalAnalysis.py -i ${eosdir} -n 8 -q ${queue} -o ${outdir} --era ${ERA} --runPeriod B --MCOnly -m TOP::RunTop --ch 0;
	  python scripts/runLocalAnalysis.py -i ${extdir} -n 8 -q ${queue} -o ${outdir} --era ${ERA} --runPeriod B --MCOnly -m TOP::RunTop --ch 0;
        else
  	  python scripts/runLocalAnalysis.py -i ${eosdir} -n 8 -q ${queue} -o ${outdir} --era ${ERA} --runPeriod ${RUN} -m TOP::RunTop --ch 0;
  	  python scripts/runLocalAnalysis.py -i ${extdir} -n 8 -q ${queue} -o ${outdir} --era ${ERA} --runPeriod ${RUN} -m TOP::RunTop --ch 0;
        fi
        #python scripts/runLocalAnalysis.py -i /store/group/phys_top/byates/LJets2016/8db9ad6/ -o LJets2015/2016/ --method TOP::RunTopKalman --era era2016 --runPeriod BCDEFGH -q 8nh
	;;
    MERGESEL )
	./scripts/mergeOutputs.py ${outdir} True;	
	;;
    PLOTSEL )
	#python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/${ERA}/samples.json -l ${lumi} --saveLog;# --mcUnc ${lumiUnc};	
	#python scripts/plotter.py -i ${outdir} -j data/${ERA}/samples.json -l ${lumi} --saveLog --run ${RUN};# --mcUnc ${lumiUnc};	
        echo "Processing with lumi=${lumi} pb^-1"
	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr_${RUN}  -j data/${ERA}/samples.json -l ${lumiFile} --saveLog --run ${RUN};# --mcUnc ${lumiUnc};	
	#python scripts/plotter.py -i LJets2015/2016 --puNormSF puwgtctr -j data/era2016/samples.json -l data/era2016/lumi.json --saveLog --run BCDEFGH
        #python scripts/plotter.py -i ${outdir} -j data/${ERA}/samples.json -l ${lumi};
	;;
    WWWSEL )
	mkdir -p ${outdir}/plots/${RUN}
	mkdir -p ${outdir}/plots/${RUN}/log
	mkdir -p ${outdir}/plots/${RUN}/ee
	mkdir -p ${outdir}/plots/${RUN}/e
	mkdir -p ${outdir}/plots/${RUN}/em
	mkdir -p ${outdir}/plots/${RUN}/mumu
	mkdir -p ${outdir}/plots/${RUN}/mu
	mkdir -p ${wwwdir}/${RUN}
	mkdir -p ${wwwdir}/${RUN}/log
	mkdir -p ${wwwdir}/${RUN}/ee
	mkdir -p ${wwwdir}/${RUN}/e
	mkdir -p ${wwwdir}/${RUN}/em
	mkdir -p ${wwwdir}/${RUN}/mumu
	mkdir -p ${wwwdir}/${RUN}/mu
	mv  ${outdir}/plots/*_${RUN}.{png,pdf} ${outdir}/plots/${RUN}/
	mv  ${outdir}/plots/*_${RUN}_log.{png,pdf} ${outdir}/plots/${RUN}/
        rename _${RUN}. . ${outdir}/plots/${RUN}/*_${RUN}.*
        rename _${RUN}_log. _log. ${outdir}/plots/${RUN}/*_${RUN}_log.*
        mv ${outdir}/plots/${RUN}/*_log*.{png,pdf} ${outdir}/plots/${RUN}/log/
        mv ${outdir}/plots/${RUN}/*_ee_*.{png,pdf} ${outdir}/plots/${RUN}/ee/
        mv ${outdir}/plots/${RUN}/*_em_*.{png,pdf} ${outdir}/plots/${RUN}/em/
        mv ${outdir}/plots/${RUN}/*_mm_*.{png,pdf} ${outdir}/plots/${RUN}/mumu/
        mv ${outdir}/plots/${RUN}/*_e_*.{png,pdf} ${outdir}/plots/${RUN}/e/
        mv ${outdir}/plots/${RUN}/*_m_*.{png,pdf} ${outdir}/plots/${RUN}/mu/
	cp -p ${outdir}/plots/${RUN}/*.{png,pdf} ${wwwdir}/${RUN}/
	cp -p ${outdir}/plots/${RUN}/log/*.{png,pdf} ${wwwdir}/${RUN}/log/
	cp -p ${outdir}/plots/${RUN}/ee/*.{png,pdf} ${wwwdir}/${RUN}/ee/
	cp -p ${outdir}/plots/${RUN}/em/*.{png,pdf} ${wwwdir}/${RUN}/em/
	cp -p ${outdir}/plots/${RUN}/e/*.{png,pdf} ${wwwdir}/${RUN}/e/
	cp -p ${outdir}/plots/${RUN}/mumu/*.{png,pdf} ${wwwdir}/${RUN}/mumu/
	cp -p ${outdir}/plots/${RUN}/mu/*.{png,pdf} ${wwwdir}/${RUN}/mu/
	cp -p test/index.php ${wwwdir}/
	cp -p test/index.php ${wwwdir}/${RUN}/
	cp -p test/index.php ${wwwdir}/${RUN}/log/
	cp -p test/index.php ${wwwdir}/${RUN}/ee/
	cp -p test/index.php ${wwwdir}/${RUN}/e/
	cp -p test/index.php ${wwwdir}/${RUN}/em/
	cp -p test/index.php ${wwwdir}/${RUN}/mumu/
	cp -p test/index.php ${wwwdir}/${RUN}/mu/
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
