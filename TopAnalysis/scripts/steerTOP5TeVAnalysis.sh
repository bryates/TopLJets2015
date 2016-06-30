#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerTOP5TeVAnalysis.sh <SEL/MERGE/PLOT/WWW>";
    echo "        SEL          - selects data and MC";
    echo "        MERGE        - merge the output of the jobs";
    echo "        PLOT         - runs the plotter tool on the selection";
    echo "        WWW          - moves the plots to an afs-web-based area";    
    exit 1; 
fi


queue=local
sourcedir=/store/cmst3/user/psilva/LJets5TeV
outdir=~/work/TopLJets5TeV
wwwdir=~/www/TopLJets5TeV

RED='\e[31m'
NC='\e[0m'
case $WHAT in
    SEL )
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir} --era era5TeV -m Run5TeVAnalysis::Run5TeVAnalysis --ch 13;
	;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir} True;	
	;;
    PLOT )
	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/${ERA}/samples.json -l ${lumi};	
	;;
    WWW )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
esac