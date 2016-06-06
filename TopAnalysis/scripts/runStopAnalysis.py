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
Analysis loop
"""
def runStopAnalysis(fileName,outFileName):
        
    print '....analysing',fileName,'with output @',outFileName

    #book histograms
    observablesH={}
    for k in ['emu','ll']:
        observablesH['mlb_'+k]=ROOT.TH1F('mlb_'+k,';Mass(lepton,b) [GeV];Events',20,0,250)
        observablesH['dphibb_'+k]=ROOT.TH1F('dphibb_'+k,';#Delta#phi(b,#bar{b}) [rad];Events',20,0,3.15)
        observablesH['dphill_'+k]=ROOT.TH1F('dphill_'+k,';#Delta#phi(l,l) [rad];Events',20,0,3.15)
        observablesH['dphijj_'+k]=ROOT.TH1F('dphijj_'+k,';#Delta#phi(l,l) [rad];Events',20,0,3.15)
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

        if abs(tree.cat)<100 : continue
        evcat = 'emu' if abs(tree.cat)==11*13 else 'll'

        evWeight=puNormSF*tree.weight[0]

        #preselect the b-jets
        bjets,otherjets=[],[]
        for ij in xrange(0,tree.nj):

            btagVal=(tree.j_btag[ij] & 0x1)
            if btagVal!=0 and len(bjets)<2: 
                bjets.append( (lVec(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij]),
                               lVec(tree.gj_pt[ij],tree.gj_eta[ij],tree.gj_phi[ij],tree.gj_m[ij]))
                              )
            else:
                otherjets.append( (lVec(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij]),
                                   lVec(tree.gj_pt[ij],tree.gj_eta[ij],tree.gj_phi[ij],tree.gj_m[ij]))
                                  )        
        if len(bjets)!=2: continue

        #leptons
        leptons=[]
        for il in xrange(0,tree.nl):
            leptons.append( (lVec(tree.l_pt[il],tree.l_eta[il],tree.l_phi[il],tree.l_m[il]),
                             lVec(tree.gl_pt[il],tree.gl_eta[il],tree.gl_phi[il],tree.gl_m[il])) )


        #met
        metx,mety=tree.met_pt*ROOT.TMath.Cos(tree.met_phi),tree.met_pt*ROOT.TMath.Sin(tree.met_phi)

        #try to solve the kinematics (need to swap bl assignments)
        allSols=[]
        try:
            sols=doubleNeutrinoSolutions( (bjets[0][0],   bjets[1][0]), 
                                          (leptons[0][0], leptons[1][0]),
                                          (metx,mety) )
            for isol in xrange(0,len(sols.nunu_s)):               
                top=bjets[0][0]+leptons[0][0]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][0],0.)
                top_=bjets[1][0]+leptons[1][0]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][1],0.)
                allSols.append( (0,top,top_) )
        except numpy.linalg.linalg.LinAlgError:
            pass        
        try:
            sols=doubleNeutrinoSolutions( (bjets[0][0],   bjets[1][0]), 
                                          (leptons[1][0], leptons[0][0]),
                                          (metx,mety) )
            for isol in xrange(0,len(sols.nunu_s)):
                top=bjets[0][0]+leptons[1][0]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][0],0.)
                top_=bjets[1][0]+leptons[0][0]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][0],0.)
                allSols.append( (0,top,top_) )
        except numpy.linalg.linalg.LinAlgError :
            pass

        #sort solutions by increasing m(ttbar)
        if len(allSols)==0: continue
        allSols=sorted(allSols, key=lambda sol: (sol[1]+sol[2]).mass() )
        l1idx=0 if allSols[0][0]==0 else 1
        l2idx=1 if allSols[0][0]==0 else 0
        observablesH['mlb_'+evcat].Fill((bjets[0][0]+leptons[l1idx][0]).mass(),evWeight)
        observablesH['mlb_'+evcat].Fill((bjets[0][0]+leptons[l2idx][0]).mass(),evWeight)
        observablesH['dphibb_'+evcat].Fill(ROOT.Math.VectorUtil.DeltaPhi(bjets[0][0],bjets[1][0]),evWeight)
        observablesH['dphill_'+evcat].Fill(ROOT.Math.VectorUtil.DeltaPhi(leptons[l1idx][0],leptons[l2idx][1]),evWeight)
        if len(otherjets)>=2:
            observablesH['dphijj_'+evcat].Fill(ROOT.Math.VectorUtil.DeltaPhi(otherjets[0][0],otherjets[1][0]),evWeight)
         
    #save results
    fOut=ROOT.TFile.Open(outFileName,'RECREATE')
    for var in observablesH: 
        observablesH[var].Write()
    fOut.Close()

 
"""
Wrapper for when the analysis is run in parallel
"""
def runStopAnalysisPacked(args):
    try:
        fileNames,outFileName=args
        runStopAnalysis(fileNames,outFileName)
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
            pool.map(runStopAnalysisPacked,tasklist)
        else:
            for fileName,outFileName in tasklist:
                runStopAnalysis(fileName,outFileName)
    else:
        cmsswBase=os.environ['CMSSW_BASE']
        for fileName,_ in tasklist:
            localRun='python %s/src/TopLJets2015/TopAnalysis/scripts/runStopAnalysis.py -i %s -o %s -q local'%(cmsswBase,fileName,opt.output)
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
