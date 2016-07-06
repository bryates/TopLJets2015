#!/bin/bash

ERA=era2015
queue=2nd
eosdir=/store/cmst3/user/psilva/LJets2015/8c1e7c9
eosdatadir=/store/cmst3/user/psilva/LJets2015/13c9e4c
summaryeosdir=/store/cmst3/group/top/summer2016/TopWidth_${ERA}
outdir=~/work/TopWidth_era2015/
cardsdir=~/work/TopWidth_era2015/datacards
wwwdir=~/www/
lumi=2267.84

lfs=(EE EM MM)
wid=(0p5w 1p0w 1p5w 2p0w 2p5w 3p0w 3p5w 4p0w 4p5w 5p0w)
cat=(1b 2b)
lbCat=(highpt lowpt)

RED='\e[31m'
NC='\e[0m'

CMSSW_7_4_7dir=~/CMSSW_7_4_7/src/
CMSSW_7_6_3dir=~/CMSSW_8_0_8_patch1/src/

# highptEM2b 
# --> +++++ highptEM1b
# --> +++++ lowptEM2b + lowptEM1b
# --> +++++ highptEE2b + highptEE1b + lowptEE2b + lowptEE1b  
# --> +++++ highptMM2b + highptMM1b + lowptMM2b + lowptMM1b  


for twid in ${wid[*]}
do
    echo "width: ${twid}"
    echo " - merging cards..."
    # high1b + high2b EM
    python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
        highptEM2b=${cardsdir}/datacard__${twid}_highptEM2b.dat \
        highptEM1b=${cardsdir}/datacard__${twid}_highptEM1b.dat \
        > ${cardsdir}/datacard__${twid}_step2.dat


    # full EM ch
    python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
        highptEM2b=${cardsdir}/datacard__${twid}_highptEM2b.dat \
        highptEM1b=${cardsdir}/datacard__${twid}_highptEM1b.dat \
        lowptEM2b=${cardsdir}/datacard__${twid}_lowptEM2b.dat \
        lowptEM1b=${cardsdir}/datacard__${twid}_lowptEM1b.dat \
        > ${cardsdir}/datacard__${twid}_step3.dat

    # 2 ch of 3
    python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
        highptEM2b=${cardsdir}/datacard__${twid}_highptEM2b.dat \
        highptEM1b=${cardsdir}/datacard__${twid}_highptEM1b.dat \
        lowptEM2b=${cardsdir}/datacard__${twid}_lowptEM2b.dat \
        lowptEM1b=${cardsdir}/datacard__${twid}_lowptEM1b.dat \
        highptEE2b=${cardsdir}/datacard__${twid}_highptEE2b.dat \
        highptEE1b=${cardsdir}/datacard__${twid}_highptEE1b.dat \
        lowptEE2b=${cardsdir}/datacard__${twid}_lowptEE2b.dat \
        lowptEE1b=${cardsdir}/datacard__${twid}_lowptEE1b.dat \
        > ${cardsdir}/datacard__${twid}_step4.dat

    # all ch
    python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
        highptEM2b=${cardsdir}/datacard__${twid}_highptEM2b.dat \
        highptEM1b=${cardsdir}/datacard__${twid}_highptEM1b.dat \
        lowptEM2b=${cardsdir}/datacard__${twid}_lowptEM2b.dat \
        lowptEM1b=${cardsdir}/datacard__${twid}_lowptEM1b.dat \
        highptEE2b=${cardsdir}/datacard__${twid}_highptEE2b.dat \
        highptEE1b=${cardsdir}/datacard__${twid}_highptEE1b.dat \
        lowptEE2b=${cardsdir}/datacard__${twid}_lowptEE2b.dat \
        lowptEE1b=${cardsdir}/datacard__${twid}_lowptEE1b.dat \
        highptMM2b=${cardsdir}/datacard__${twid}_highptMM2b.dat \
        highptMM1b=${cardsdir}/datacard__${twid}_highptMM1b.dat \
        lowptMM2b=${cardsdir}/datacard__${twid}_lowptMM2b.dat \
        lowptMM1b=${cardsdir}/datacard__${twid}_lowptMM1b.dat \
        > ${cardsdir}/datacard__${twid}_step5.dat

    echo " - making workspaces..."

    #//cd ${CMSSW_7_4_7dir}
    #//eval `scramv1 runtime -sh`

    text2workspace.py ${cardsdir}/datacard__${twid}_highptEM2b.dat -P \
        HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest \
        -m 172.5 --PO verbose --PO altSignal=${twid} --PO muFloating \
        -o ${outdir}/${twid}_step1.root

    for i in 2 3 4 5
    do
        text2workspace.py ${cardsdir}/datacard__${twid}_step${i}.dat -P \
            HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest \
            -m 172.5 --PO verbose --PO altSignal=${twid} --PO muFloating \
            -o ${outdir}/${twid}_step${i}.root
    done

    echo " - running combine"

    cd ${outdir}
    for i in 1 2 3 4 5
    do
        for j in 0 1
        do
            combine ${outdir}/${twid}_step${i}.root -M MultiDimFit \
                -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200 \
                -t -1 --expectSignal=1 --setPhysicsModelParameters x=${j},r=1 \
                -n x${j}_scan_Asimov_${twid}_step${i}
        done
    done

    echo " - producing CLs stats"

    cd ${outdir}
    for i in 1 2 3 4 5
    do
        combine ${twid}_step${i}.root -M HybridNew --seed 8192 --saveHybridResult \
            -m 172.5  --testStat=TEV --singlePoint 1 -T 500 -i 2 --fork 6 \
            --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0 \
            --expectedFromGrid 0.5 -n x_pre-fit_exp_${twid}_step${i} \
            &> ${outdir}/x_pre-fit_exp__${twid}_step${i}.log
    done

    echo " - throwing toys"

    cd ${outdir}
    for i in 1 2 3 4 5
    do
        root -l -q -b higgsCombinex_pre-fit_exp_${twid}_step${i}.HybridNew.mH172.5.8192.quant0.500.root \
            "${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/hypoTestResultTreeTopWid.cxx(\"x_pre-fit_exp__${twid}_step${i}.qvals.root\",172.5,1,\"x\",1000,\"step${i}\",\"${twid}\")"
    done

    echo " - getting quantiles plot"

    python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
        -i ${outdir}/ -o ${outdir}/${twid}_ --lfs step1,step2,step3,step4,step5 --wid ${twid}   \
        --axisOverwrite hipt_EM_2b,+hipt_EM_1b,+_lopt_EM_*b,+_all_EE,+_all_MM


done

# get likelihood scan plots
for wid in ${wid[*]}
do
    python getLikelihoodScans.py -i ${outdir}/ -o ${outdir}/
done

# copy to web directory
mkdir -p ${wwwdir}/ana07032016
cp ${outdir}/attempt1/plots/*.{png,pdf} ${wwwdir}/ana07032016/       
cp ${outdir}/*.{png,pdf} ${wwwdir}/ana07032016/        
cp ${wwwdir}/ana/index.php ${wwwdir}/ana07032016/
