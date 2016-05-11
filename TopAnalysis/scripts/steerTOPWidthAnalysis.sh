#!/bin/bash

WHAT=$1; 
if [[ "$1" == "" ]]; then 
    echo "steerTOPWidthAnalysis.sh <SEL/MERGE/PLOT>";
    echo "        SEL     - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGE   - merge the output of the jobs";
    echo "        PLOT    - runs the plotter tool"
    echo "        WWW     - moves the plots to an afs-web-based area"
    exit 1; 
fi

queue=8nh
eosdir=/store/cmst3/user/psilva/LJets2015/8c1e7c9
outdir=~/work/TopWidth
wwwdir=~/www/TopWidth
lumi=2267.84

RED='\e[31m'
NC='\e[0m'

case $WHAT in
    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir} -m TOPWidth::RunTopWidth --ch 0;
	;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir};	
	;;
    PLOT )
	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/samples_Run2015.json -l ${lumi};	
	;;
    WWW )
	mkdir -p ${wwwdir}
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}
	cp test/index.php ${wwwdir}
	;;
esac