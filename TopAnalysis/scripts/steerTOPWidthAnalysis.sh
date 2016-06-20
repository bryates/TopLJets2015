#!/bin/bash

WHAT=$1; 
ERA=$2
if [ "$#" -ne 2 ]; then 
    echo "steerTOPWidthAnalysis.sh <SELDATA/SELMC/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/PLOT/WWW> <ERA>";
    echo "        SEL{DATA,MC} - launches selection jobs to the batch on DATA or MC, output will contain summary trees and control plots"; 
    echo "        MERGESEL     - merge the output of the jobs";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    echo "        ANA          - analyze the selected events";
    echo "        MERGE        - merge the output of the analysis jobs";
    echo "        PLOT         - runs the plotter tool on the analysis outputs";
    echo "        WWW          - moves the analysis plots to an afs-web-based area";
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N

queue=2nd
githash=7e62835
lumi=589
case $ERA in
    era2015)
	githash=8c1e7c9;
	lumi=2267.84
	;;
esac

eosdir=/store/cmst3/user/psilva/LJets2015/${githash}
summaryeosdir=/store/cmst3/group/top/summer2016/TopWidth_${ERA}
outdir=~/work/TopWidth_${ERA}
wwwdir=~/www/TopWidth_${ERA}


RED='\e[31m'
NC='\e[0m'

case $WHAT in
    SELDATA )
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir} --era ${ERA} -m TOPWidth::RunTopWidth --ch 0 --only Data;
	;;
    SELMC )
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir} --era ${ERA} -m TOPWidth::RunTopWidth --ch 0 --only ^Data;
	;;
    MERGESEL )
	./scripts/mergeOutputs.py ${outdir} True;	
	;;
    PLOTSEL )
	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/${ERA}/samples.json -l ${lumi};	
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
    ANA )
	python scripts/runTopWidthAnalysis.py -i ${summaryeosdir} -o ${outdir}/analysis -q ${queue};
	;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir}/analysis;
	;;
    PLOT )
        python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json      -l ${lumi} --onlyData;        
	python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/syst_samples.json -l ${lumi} --silent;        
        ;;
    WWW )
        mkdir -p ${wwwdir}/ana
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana        
        cp test/index.php ${wwwdir}/ana
	;;
esac