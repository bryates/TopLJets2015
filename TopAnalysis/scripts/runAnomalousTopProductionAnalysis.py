#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import numpy
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist
from TopLJets2015.TopAnalysis.nuSolutions import *

"""
a dummy converter
"""
def convertToPtEtaPhiM(lVec,xyz,m=0.):    
    en=ROOT.TMath.Sqrt(xyz[0]**2+xyz[1]**2+xyz[2]**2)
    p4=ROOT.TLorentzVector(xyz[0],xyz[1],xyz[2],en)
    return lVec(p4.Pt(),p4.Eta(),p4.Phi(),p4.M())

"""
stop quark/chi0 filter the arguments are the masses to filter out
"""
def stopChiFilter(tree,filtArgs):
    weight=-1.0
    mstop,mchi0=0,0
    for i in xrange(0,tree.nt):
        if abs(tree.t_id[i])==100006 : mstop=tree.t_id[i]
        if abs(tree.t_id[i])==100022 : mchi0=tree.t_id[i]
    if mstop==float(filtArgs[0]) and mstop==float(filtArgs[1]): weight=1.0
    return weight


"""
Analysis loop
"""
def runAnomalousTopProductionAnalysis(fileName,outFileName,filterName):
        
    print '....analysing',fileName,'with output @',outFileName

    #book histograms
    observablesH={}
    for k in ['emu','ll']:
        observablesH['ht_'+k]=ROOT.TH1F('ht_'+k,';H_{T} [GeV];Events',20,0,1000)

        observablesH['dphibb_'+k]=ROOT.TH1F('dphibb_'+k,';#Delta#phi(b,#bar{b}) [rad];Events',20,0,3.15)
        observablesH['cosbstar_'+k]=ROOT.TH1F('cosbstar_'+k,';cos(#theta*_{b});Events',20,-1,1)
        observablesH['cosbstarprod_'+k]=ROOT.TH1F('cosbstarprod_'+k,';cos(#theta*_{b_{1}})cos(#theta*_{b_{2}});Events',20,-1,1)

        observablesH['dphill_'+k]=ROOT.TH1F('dphill_'+k,';#Delta#phi(l,l) [rad];Events',20,0,3.15)
        observablesH['coslstar_'+k]=ROOT.TH1F('coslstar_'+k,';cos(#theta*_{l});Events',20,-1,1)
        observablesH['coslstarprod_'+k]=ROOT.TH1F('coslstarprod_'+k,';cos(#theta*_{l_{1}})cos(#theta*_{l_{2}});Events',20,-1,1)

        observablesH['dphijj_'+k]=ROOT.TH1F('dphijj_'+k,';#Delta#phi(j,j) [rad];Events',20,0,3.15)
    for var in observablesH:
        observablesH[var].SetDirectory(0)
        observablesH[var].Sumw2()


    #open file
    puNormSF=1.0
    if 'MC13TeV' in fileName:
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
    lVec = ROOT.Math.LorentzVector(ROOT.Math.PtEtaPhiM4D('double')) 
    for i in xrange(0,totalEntries):

        tree.GetEntry(i)

        if i%100==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )

        filtWeight=1.0
        if filterName:
            filtFunc,filtArgs=filterName.split('=')
            filtWeight=globals()[filtFunc](tree,filtArgs.split(','))
            if filtWeight<0 : continue

        if abs(tree.cat)<100 : continue
        evcat = 'emu' if abs(tree.cat)==11*13 else 'll'

        evWeight=puNormSF*tree.weight[0]

        #leptons
        leptons=[]
        for il in xrange(0,tree.nl):
            leptons.append( lVec(tree.l_pt[il],tree.l_eta[il],tree.l_phi[il],tree.l_m[il]) )
        if len(leptons)<2 : continue

        #preselect the b-jets (save always the jet and the gen jet)
        bjets,otherjets=[],[]
        for ij in xrange(0,tree.nj):

            btagVal=(tree.j_btag[ij] & 0x1)
            if btagVal!=0 and len(bjets)<2: 
                bjets.append( lVec(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij]) )
            else:
                otherjets.append( lVec(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij]) )
        if len(bjets)!=2: continue

        #met
        metx,mety=tree.met_pt*ROOT.TMath.Cos(tree.met_phi),tree.met_pt*ROOT.TMath.Sin(tree.met_phi)

        #try to solve the kinematics (need to swap bl assignments)
        allSols=[]
        try:
            sols=doubleNeutrinoSolutions( (bjets[0],   bjets[1]), 
                                          (leptons[0], leptons[1]),
                                          (metx,mety) )
            for isol in xrange(0,len(sols.nunu_s)):               
                top=bjets[0]+leptons[0]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][0],0.)
                top_=bjets[1]+leptons[1]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][1],0.)
                allSols.append( (0,top,top_) )
        except numpy.linalg.linalg.LinAlgError:
            pass        
        try:
            sols=doubleNeutrinoSolutions( (bjets[0],   bjets[1]), 
                                          (leptons[1], leptons[0]),
                                          (metx,mety) )
            for isol in xrange(0,len(sols.nunu_s)):
                top=bjets[0]+leptons[1]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][0],0.)
                top_=bjets[1]+leptons[0]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][1],0.)
                allSols.append( (1,top,top_) )
        except numpy.linalg.linalg.LinAlgError :
            pass

        #sort solutions by increasing m(ttbar)
        if len(allSols)==0: continue
        allSols=sorted(allSols, key=lambda sol: (sol[1]+sol[2]).mass() )        
        l1idx=0 if allSols[0][0]==0 else 1
        l2idx=1 if allSols[0][0]==0 else 0
 
        #setup the Lorentz transformations to the top/anti-top rest frames
        topBoost, top_Boost = allSols[0][1].BoostToCM(), allSols[0][2].BoostToCM()

        #measure b-jet angles
        cosb1 = ROOT.TMath.Cos( ROOT.Math.VectorUtil.Angle( ROOT.Math.VectorUtil.boost(bjets[0],topBoost), allSols[0][1] ) )
        cosb2 = ROOT.TMath.Cos( ROOT.Math.VectorUtil.Angle( ROOT.Math.VectorUtil.boost(bjets[1],top_Boost), allSols[0][2] ) )
        observablesH['dphibb_'+evcat].Fill(ROOT.Math.VectorUtil.DeltaPhi(bjets[0],bjets[1]),evWeight)
        observablesH['cosbstar_'+evcat].Fill(cosb1,evWeight)
        observablesH['cosbstar_'+evcat].Fill(cosb2,evWeight)
        observablesH['cosbstarprod_'+evcat].Fill(cosb1*cosb2,evWeight)

        #measure leptonic angles
        cosl1 = ROOT.TMath.Cos( ROOT.Math.VectorUtil.Angle( ROOT.Math.VectorUtil.boost(leptons[l1idx],topBoost), allSols[0][1] ) )
        cosl2 = ROOT.TMath.Cos( ROOT.Math.VectorUtil.Angle( ROOT.Math.VectorUtil.boost(leptons[l2idx],top_Boost), allSols[0][2] ) )
        observablesH['dphill_'+evcat].Fill(ROOT.Math.VectorUtil.DeltaPhi(leptons[l1idx],leptons[l2idx]),evWeight)
        observablesH['coslstar_'+evcat].Fill(cosl1,evWeight)
        observablesH['coslstar_'+evcat].Fill(cosl2,evWeight)
        observablesH['coslstarprod_'+evcat].Fill(cosl1*cosl2,evWeight)

        if len(otherjets)>=2:
            observablesH['dphijj_'+evcat].Fill(ROOT.Math.VectorUtil.DeltaPhi(otherjets[0],otherjets[1]),evWeight)
         
        #other control variables
        observablesH['ht_'+evcat].Fill(bjets[0].pt()+bjets[1].pt()+leptons[0].pt()+leptons[1].pt()+tree.met_pt,evWeight)

    #save results
    fOut=ROOT.TFile.Open(outFileName,'RECREATE')
    for var in observablesH: 
        observablesH[var].Write()
    fOut.Close()

 
"""
Wrapper for when the analysis is run in parallel
"""
def runAnomalousTopProductionAnalysisPacked(args):
    try:
        fileNames,outFileName,filterName=args
        runAnomalousTopProductionAnalysis(fileNames,outFileName,filterName)
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
        tasklist.append((filename,'%s/%s'%(opt.output,baseFileName),opt.filter))

    #loop over tasks
    if opt.queue=='local':
        if opt.jobs>1:
            print ' Submitting jobs in %d threads' % opt.jobs
            import multiprocessing as MP
            pool = MP.Pool(opt.jobs)
            pool.map(runAnomalousTopProductionAnalysisPacked,tasklist)
        else:
            for fileName,outFileName,filterName in tasklist:
                runAnomalousTopProductionAnalysis(fileName,outFileName,filterName)
    else:
        cmsswBase=os.environ['CMSSW_BASE']
        for fileName,_,filterName in tasklist:
            localRun='python %s/src/TopLJets2015/TopAnalysis/scripts/runAnomalousTopProductionAnalysis.py -i %s -o %s -q local'%(cmsswBase,fileName,opt.output)
            if filterName : localRun+=' --filter %s'%filterName
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
                          default='/afs/cern.ch/user/p/psilva/work/Stop',
                          help='input directory with the files [default: %default]')
	parser.add_option('-f', '--filter',
                          dest='filter',   
                          default=None,
                          help='apply this filter function')
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
