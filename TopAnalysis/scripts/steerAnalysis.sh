#!/bin/bash

WHAT=$1; if [[ "$1" == "" ]]; then echo "steerAnalysis.sh <SEL/MERGE/PLOT/BKG>"; exit 1; fi

queue=8nh
eosdir=/store/cmst3/user/psilva/LJets2015/64217e8
outdir=~/work/LJets2015/
lumi=2134

RED='\e[31m'
NC='\e[0m'

case $WHAT in
    SEL )
	echo -e "[ ${RED} Submitting the selection for the signal regions ${NC} ]"
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_muplus   --ch 13   --charge 1
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_muminus  --ch 13   --charge -1
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_eplus   --ch 11   --charge 1
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_eminus  --ch 11   --charge -1
	
	echo -e "[ ${RED} Submitting the selection for the control regions ${NC} ]"
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue}            -o ${outdir}/analysis_munoniso --ch 1300
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue}            -o ${outdir}/analysis_enoniso --ch 1100
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_z --ch 21
	;;
    MERGE )
	a=(muplus muminus eplus eminus munoniso enoniso z)
	for i in ${a[@]}; do
	    echo -e "[ ${RED} Merging ${i} ${NC} ]"
	    ./scripts/mergeOutputs.py ${outdir}/analysis_${i};
	done
	;;
    PLOT )
	a=(muplus muminus eplus eminus munoniso enoniso z)
	for i in ${a[@]}; do
	    echo -e "[ ${RED} Creating plotter for ${i} ${NC} ]"
	    python scripts/plotter.py -i ${outdir}/analysis_${i}/ --puNormSF puwgtctr  -j data/samples_Run2015.json -l ${lumi} --saveLog
	done
	a=(muplus muminus eplus eminus)
	for i in ${a[@]}; do
	    echo -e "[ ${RED} Creating plotter for ${i} ${NC} ]"
	    python scripts/plotter.py -i ${outdir}/analysis_${i}/ --puNormSF puwgtctr  -j data/syst_samples_Run2015.json -l ${lumi} -o syst_plotter.root --silent
	done
	;;
    BKG )
	a=(mu e)
	b=(plus minus)
	for i in ${a[@]}; do
	    for j in ${b[@]}; do
		python scripts/runQCDEstimation.py --iso ${outdir}/analysis_${i}${j}/plots/plotter.root --noniso ${outdir}/analysis_${i}noniso/plots/plotter.root --out ${outdir}/analysis_${i}${j}/
	    done
	done
	;;
esac