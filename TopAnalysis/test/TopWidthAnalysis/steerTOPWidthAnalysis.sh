#!/bin/bash

WHAT=$1; 
ERA=$2;
if [ "$#" -ne 2 ]; then 
    echo "steerTOPWidthAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/PLOT/WWW>";
    echo "        WORKSPACE       - create a workspace with the TopHypoTest model";
    echo "        SCAN            - scans the likelihood";
    echo "        CLs             - calculates CLs using higgs_combine";
    echo "        TOYS            - plots toys and saves locally";
    echo "        MERGE_DATACARDS - merge all datacards for one LFS/WID into one WID datacard";
    echo "        QUANTILES       - plots the quantiles chart for all distributions, by lepton final state";
    echo "        <above>_M       - with <above> one of the commands above (past 'PLOT'), run over merged datacards"; 
    echo "        WWW             - moves the analysis plots to an afs-web-based area";
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

outdir=~/work/TopWidth_${ERA}/
cardsdir=~/work/TopWidth_${ERA}/datacards
wwwdir=~/www/TopWidth_${ERA}/
lumi=2267.84

lfs=(EE EM MM)
wid=(0p5w 1p0w 1p5w 2p0w 2p5w 3p0w 3p5w 4p0w 4p5w 5p0w)
cat=(1b 2b)
lbCat=(highpt lowpt)

RED='\e[31m'
NC='\e[0m'

CMSSW_7_4_7dir=~/CMSSW_7_4_7/src/
CMSSW_7_6_3dir=~/CMSSW_8_0_8_patch1/src/

case $WHAT in
    SHAPES )
        cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
        python test/TopWidthAnalysis/createShapesFromPlotter.py \
                -i ${outdir}/analysis/plots/plotter.root \
                -s tbart,tW -o ${outdir}/datacards/ \
                --systInput ${outdir}/analysis/plots/syst_plotter.root
    ;;
    MERGE_DATACARDS )
        cd ${CMSSW_7_4_7dir}
        for twid in ${wid[*]}
        do
            # for a given width, merge all
            allcmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
            for tlbCat in ${lbCat[*]}
            do
                # for a given width and lbcat, merge all
                lbccmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
                for tlfs in ${lfs[*]}
                do
                    # for a given width, lfs, and lbcat, merge all
                    lfscmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
                    for tCat in ${cat[*]}
                    do
                        allcmd="${allcmd} ${tlbCat}${tlfs}${tCat}=${cardsdir}/datacard__${twid}_${tlbCat}${tlfs}${tCat}.dat "
                        lbccmd="${lbccmd} ${tlbCat}${tlfs}${tCat}=${cardsdir}/datacard__${twid}_${tlbCat}${tlfs}${tCat}.dat "
                        lfscmd="${lfscmd} ${tlbCat}${tlfs}${tCat}=${cardsdir}/datacard__${twid}_${tlbCat}${tlfs}${tCat}.dat "
                    done
                    echo $lfscmd
                    echo " "
                    lfscmd="${lfscmd} > ${cardsdir}/datacard__${twid}_${tlbCat}${tlfs}.dat"
                done
                echo $lbccmd
                echo " "
                lbccmd="${lbccmd} > ${cardsdir}/datacard__${twid}_${tlbCat}.dat"
            done

            echo $allcmd
            echo " "
            allcmd="${allcmd} > ${cardsdir}/datacard__${twid}.dat"
            lfscmd="${lfscmd} > ${cardsdir}/datacard__${twid}.dat"
            
            eval $allcmd
            eval $lbccmd
            eval $lfscmd
        done
    ;;
    WORKSPACE )
        cd ${CMSSW_7_4_7dir}
        mkdir ${outdir}
        for twid in ${wid[*]}
        do
            for tlfs in ${lfs[*]}
            do
                for tlbCat in ${lbCat[*]}
                do
                    for tcat in ${cat[*]}
                    do
                        echo "Creating workspace for ${twid}${tlfs}" 
                        text2workspace.py ${cardsdir}/datacard__${twid}_${tlfs}.dat -P \
                            HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest \
                            -m 172.5 --PO verbose --PO altSignal=${twid} --PO muFloating \
                            -o ${outdir}/${twid}_${tlfs}.root
                    done
                done
            done
            
            # All datacards
            echo "Creating workspace for ${twid}" 
            text2workspace.py ${cardsdir}/datacard__${twid}.dat -P \
                HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest \
                -m 172.5 --PO verbose --PO altSignal=${twid} --PO muFloating \
                -o ${outdir}/${twid}.root
        done
    ;;
    SCAN )
        cd ${outdir}
        for twid in ${wid[*]}
        do
            for tlfs in ${lfs[*]}
            do
                for tlbCat in ${lbCat[*]}
                do
                    for tcat in ${cat[*]}
                    do
                        combine ${outdir}/${twid}_${tlfs}.root -M MultiDimFit \
                            -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200 \
                            -t -1 --expectSignal=1 --setPhysicsModelParameters x=0,r=1 \
                            -n x0_scan_Asimov_${twid}_${tlfs}
                    done
                done
            done

            # All datacards
            combine ${outdir}/${twid}.root -M MultiDimFit \
                -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200 \
                -t -1 --expectSignal=1 --setPhysicsModelParameters x=0,r=1 \
                -n x0_scan_Asimov_${twid}
        done
    ;;
    CLs )
        cd ${outdir}
        for twid in ${wid[*]}
        do
            for tlfs in ${lfs[*]}
            do
                for tlbCat in ${lbCat[*]}
                do
                    for tcat in ${cat[*]}
                    do
                        combine ${outdir}/${twid}_${tlfs}.root -M HybridNew --seed 8192 --saveHybridResult \
                            -m 172.5  --testStat=TEV --singlePoint 1 -T 500 -i 2 --fork 6 \
                            --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0 \
                            --expectedFromGrid 0.5 -n x_pre-fit_exp_${twid}_${tlfs} \
                            &> ${outdir}/x_pre-fit_exp_${twid}_${tlfs}.log
                    done
                done
            done

            # All datacards
            combine ${twid}.root -M HybridNew --seed 8192 --saveHybridResult \
                -m 172.5  --testStat=TEV --singlePoint 1 -T 500 -i 2 --fork 6 \
                --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0 \
                --expectedFromGrid 0.5 -n x_pre-fit_exp_${twid} \
                &> ${outdir}/x_pre-fit_exp__${twid}.log
        done
    ;;
    TOYS )
        cd ${CMSSW_7_4_7dir}
        cd ${outdir}
        for twid in ${wid[*]}
        do
            for tlfs in ${lfs[*]}
            do
                for tlbCat in ${lbCat[*]}
                do
                    for tcat in ${cat[*]}
                    do
                        cd ${outdir}
                        # each lfs datacard
                        root -l -q -b higgsCombinex_pre-fit_exp_${twid}_${tlfs}.HybridNew.mH172.5.8192.quant0.500.root \
                            "${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTreeTopWid.cxx(\"x_pre-fit_exp_${twid}_${tlfs}.qvals.root\",172.5,1,\"x\",1000,\"${tlfs}\",\"${twid}\")"
                    done
                done
            done
            
            # All datacards
            root -l -q -b higgsCombinex_pre-fit_exp_${twid}.HybridNew.mH172.5.8192.quant0.500.root \
                "${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTreeTopWid.cxx(\"x_pre-fit_exp__${twid}.qvals.root\",172.5,1,\"x\",1000,\"\",\"${twid}\")"
        done
    ;;
    QUANTILES )
        cd ${outdir}
        lfsStr=""
        for tlfs in ${lfs[*]}
        do
          if [[ "${lfsStr}" == "" ]];
          then
            lfsStr="${tlfs}"
          else
            lfsStr="${lfsStr},${tlfs}"
          fi
        done
        
        python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidAnalysis/getQuantilesPlot.py \
               -i ${outdir}/ -o ${outdir}/ --lfs ${lfsStr} 
        
        # All datacards
        python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidAnalysis/getQuantilesPlot.py \ 
               -i ${outdir}/ -o ${outdir}/
    ;;
    WWW )
        mkdir -p ${wwwdir}/ana
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana        
        cp ${outdir}/*.{png,pdf} ${wwwdir}/ana        
        cp test/index.php ${wwwdir}/ana
	;;
esac
