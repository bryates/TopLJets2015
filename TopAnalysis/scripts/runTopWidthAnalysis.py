#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist

"""
Take the ratio of two Breit-Wigner functions at fixed mass as a reweighting factor
"""
def weightTopWidth(tmassList,bwigner,targetWidth,origWidth=1.324):
    bwigner.FixParameter(2,origWidth)
    origNorm=bwigner.Integral(0,300)

    bwigner.FixParameter(2,targetWidth)
    targetNorm=bwigner.Integral(0,300)

    wgt=1.0
    for m in tmassList:
        bwigner.FixParameter(2,origWidth)
        origVal=bwigner.Eval(m)
        bwigner.FixParameter(2,targetWidth)
        targetVal=bwigner.Eval(m)
        wgt *= (targetVal/targetNorm) / (origVal/origNorm)
    return wgt

"""
Analysis loop
"""
def runTopWidthAnalysis(fileName,
                        outFileName,
                        widthList=[0.5,1,2,4],
                        smMass=172.5,
                        smWidth=1.324,
                        systs=['','puup','pudn','btagup','btagdn','jerup','jerdn','jesup','jesdn','lesup','lesdn']):
        
    print '....analysing',fileName,'with output @',outFileName

    #check if this is data beforehand
    isData=False if 'MC13TeV' in fileName else True
    if isData:
        widthList=[1.0]
        systs=['']

    #define the relativistic Breit-Wigner function
    bwigner=ROOT.TF1('bwigner',
                     '[0]*([1]*[2]*sqrt([1]*[1]*([1]*[1]+[2]*[2]))/sqrt([1]*[1]+sqrt([1]*[1]*([1]*[1]+[2]*[2]))))/(TMath::Power(x*x-[1]*[1],2)+TMath::Power([1]*[2],2))',
                     0,300)
    bwigner.SetParName(0,"N")
    bwigner.SetParameter(0,1.0)
    bwigner.SetParName(1,"m_{t}")
    bwigner.FixParameter(1,smMass)
    bwigner.SetParName(2,"#Gamma_{t}")
    bwigner.FixParameter(2,smWidth)

    #book histograms
    observablesH={}    
    for s in systs:
        for i in ['lowpt','highpt']:
            for j in ['E','M','EE','MM','EM']:
                for b in ['1b','2b']:
                    for w in widthList:
                        var=s+i+j+b+'_mlb_%3.1fw'%w
                        observablesH[var]=ROOT.TH1F(var,';Mass(lepton,jet) [GeV];l+j pairs',50,0,300)
                        if w!=1.0 or len(s)>0 : continue
                        var=i+j+b+'_pairing'
                        observablesH[var]=ROOT.TH1F(var,';Pairing;l+j pairs',2,0,2)
                        observablesH[var].GetXaxis().SetBinLabel(1,'correct')
                        observablesH[var].GetXaxis().SetBinLabel(2,'wrong')

    for var in observablesH:
        observablesH[var].SetDirectory(0)
        observablesH[var].Sumw2()

    #open file
    puNormSF=1.0
    if isData:
        fIn=ROOT.TFile.Open(fileName)
        puCorrH=fIn.Get('puwgtctr')
        nonWgt=puCorrH.GetBinContent(1)
        wgt=puCorrH.GetBinContent(2)
        if wgt>0 : puNormSF=nonWgt/wgt
        fIn.Close()

    tree=ROOT.TChain('twev')
    tree.AddFile(fileName)

    #loop over events in the tree and fill histos
    totalEntries=tree.GetEntries()
    for i in xrange(0,totalEntries):

        tree.GetEntry(i)

        if i%100==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )

        #determine weighting factors for the width
        tmassList=[]
        for it in xrange(0,tree.nt): tmassList.append( tree.t_m[it] )
        widthWeight={}
        for w in widthList: widthWeight[w]=weightTopWidth(tmassList,bwigner,w*smWidth,smWidth)

        evcat='E'
        if abs(tree.cat)==13 : evcat='M'
        if abs(tree.cat)==11*11 : evcat='EE'
        if abs(tree.cat)==11*13 : evcat='EM'
        if abs(tree.cat)==13*13 : evcat='MM'

        #preselect the b-jets (central b-tag, b-tag up, b-tag dn, jer up, jer dn, jes up, jes dn)
        bjets=( [], [], [], [], [], [], [] )
        for ij in xrange(0,tree.nj):
                        
            jp4=ROOT.TLorentzVector()
            jp4.SetPtEtaPhiM(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij])
            
            #jres=tree.j_jer[ij]
            #jscale=tree.j_jscale[ij]

            for ibit in xrange(0,3):
                btagVal=((tree.j_btag[ij] >> ibit) & 0x1)
                if btagVal!=0:
                    bjets[ibit].append( (ij,jp4) )
                

        
        #pair with the leptons
        for il in xrange(0,tree.nl):

            stdlp4=ROOT.TLorentzVector()
            stdlp4.SetPtEtaPhiM(tree.l_pt[il],tree.l_eta[il],tree.l_phi[il],tree.l_m[il])
            lesScale=1.0
            if not isData: lesScale=ROOT.getLeptonEnergyScaleUncertainty(tree.l_id[il],lp4.Pt(),lpt.Eta())

            for s in systs:
                
                #event weight
                evWeight=puNormSF*tree.weight[0]

                lp4=ROOT.TLorentzVector(stdlp4)

                #experimental uncertainties
                ijhyp=0
                if s=='btagup' : ijhyp=1
                if s=='btagdn' : ijhyp=2
                if s=='jerup'  : ijhyp=3
                if s=='jerdn'  : ijhyp=4
                if s=='jesup'  : ijhyp=5
                if s=='jesdn'  : ijhyp=6
                if s=='lesup'  : lp4 *= (1.0+lesScale)
                if s=='lesdn'  : lp4 *= (1.0-lesScale)
                
                #btag hypothesis
                nbtags=len(bjets[ijhyp])
                if nbtags<1 : continue
                if nbtags>2 : nbtags=2
                btagcat='1b' if nbtags==1 else '2b'
                
                #pileup
                if s=='puup' : evWeight=puNormSF*tree.weight[1]
                if s=='pudn' : evWeight=puNormSF*tree.weight[2]

                for ib in xrange(0,nbtags):

                    ij,jp4 = bjets[ijhyp][ib]

                    #check if assignment is correct or not
                    assignmentType=1
                    if tree.gl_id[il]!=0 and tree.gj_flav[ij]!=0 and tree.nt>0:
                        if tree.gl_id[il]*tree.gj_flav[ij]<0 : assignmentType=0
                
                    #kinematics of the l,b system
                    mlb=(lp4+jp4).M()
                    ptlb=(lp4+jp4).Pt()
                    ptCat='lowpt' if ptlb<100 else 'highpt'
                        
                    #fill histos
                    for w in widthList:
                        var=s+ptCat+evcat+btagcat+'_mlb_%3.1fw'%w
                        observablesH[var].Fill(mlb,evWeight*widthWeight[w])

                        #only for standard width and syst variations
                        if w!=1.0 or len(s)>0 : continue
                        var=s+ptCat+evcat+btagcat+'_pairing'
                        observablesH[var].Fill(assignmentType,evWeight*widthWeight[w])

         
    #save results
    fOut=ROOT.TFile.Open(outFileName,'RECREATE')
    for var in observablesH: observablesH[var].Write()
    fOut.Close()

 
"""
Wrapper for when the analysis is run in parallel
"""
def runTopWidthAnalysisPacked(args):
    try:
        fileNames,outFileName=args
        runTopWidthAnalysis(fileNames,outFileName)
    except : # ReferenceError:
        print 50*'<'
        print "  Problem with", name, "continuing without"
        print 50*'<'
        return False
    
"""
Create analysis tasks
"""
def createAnalysisTasks(opt):

    onlyList=opt.only.split('v')

    ## Local directory
    file_list=[]
    if os.path.isdir(opt.input):
        for file_path in os.listdir(opt.input):
            if file_path.endswith('.root'):
                file_list.append(os.path.join(opt.input,file_path))
    elif opt.input.startswith('/store/'):
        file_list = getEOSlslist(opt.input)
    elif '.root' in opt.input:
        file_list.append(opt.input)

    #list of files to analyse
    tasklist=[]
    for filename in file_list:
        baseFileName=os.path.basename(filename)      
        tag,ext=os.path.splitext(baseFileName)
        if len(onlyList)>0:
            processThis=False
            for filtTag in onlyList:
                if filtTag in tag:
                    processThis=True
            if not processThis : continue
        tasklist.append((filename,'%s/%s'%(opt.output,baseFileName)))

    #loop over tasks
    if opt.queue=='local':
        if opt.jobs>1:
            print ' Submitting jobs in %d threads' % opt.jobs
            import multiprocessing as MP
            pool = MP.Pool(opt.jobs)
            pool.map(runTopWidthAnalysisPacked,tasklist)
        else:
            for fileName,outFileName in tasklist:
                runTopWidthAnalysis(fileName,outFileName)
    else:
        cmsswBase=os.environ['CMSSW_BASE']
        for fileName,_ in tasklist:
            localRun='python %s/src/TopLJets2015/TopAnalysis/scripts/runTopWidthAnalysis.py -i %s -o %s -q local'%(cmsswBase,fileName,opt.output)
            cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localRun)
            print cmd
            os.system(cmd)

"""
steer
"""
def main():
	usage = 'usage: %prog [options]'
	parser = optparse.OptionParser(usage)
	parser.add_option('-i', '--input',
                          dest='input',   
                          default='/afs/cern.ch/user/p/psilva/work/TopWidth',
                          help='input directory with the files [default: %default]')
	parser.add_option('--jobs',
                          dest='jobs', 
                          default=1,
                          type=int,
                          help='# of jobs to process in parallel the trees [default: %default]')
	parser.add_option('--only',
                          dest='only', 
                          default='',
                          type='string',
                          help='csv list of tags to process')
	parser.add_option('-o', '--output',
                          dest='output', 
                          default='analysis',
                          help='Output directory [default: %default]')
	parser.add_option('-q', '--queue',
                          dest='queue',
                          default='local',
                          help='Batch queue to use [default: %default]')
	(opt, args) = parser.parse_args()

        ROOT.FWLiteEnabler.enable() 
	os.system('mkdir -p %s' % opt.output)

        createAnalysisTasks(opt)
        

if __name__ == "__main__":
	sys.exit(main())
