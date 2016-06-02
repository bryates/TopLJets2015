#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist

"""
Analysis loop
"""
def runTopRadiusAnalysis(fileName,outFileName):
        
    print '....analysing',fileName,'with output @',outFileName

    #book histograms
    observablesH={}
    for k in ['emu','ll']:
        observablesH['ht_'+k]=ROOT.TH1F('ht_'+k,';H_{T} [GeV];Events',20,0,1000)
        observablesH['dphill_'+k]=ROOT.TH1F('dphill_'+k,';#Delta#phi(l,l) [rad];Events',20,0,2*ROOT.TMath.Pi())
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
    for i in xrange(0,totalEntries):

        tree.GetEntry(i)

        if i%100==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )

        if abs(tree.cat)<100 : continue
        evcat = 'emu' if abs(tree.cat)==11*13 else 'll'

        evWeight=puNormSF*tree.weight[0]

        #preselect the b-jets
        bjets=[]
        for ij in xrange(0,tree.nj):

            btagVal=(tree.j_btag[ij] & 0x1)
            if btagVal==0: continue
            
            jp4=ROOT.TLorentzVector()
            jp4.SetPtEtaPhiM(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij])
            gjp4=ROOT.TLorentzVector()
            gjp4.SetPtEtaPhiM(tree.gj_pt[ij],tree.gj_eta[ij],tree.gj_phi[ij],tree.gj_m[ij])
            bjets.append( (jp4,gjp4) )
            if len(bjets)==2 : break
        if len(bjets)!=2: continue

        #leptons
        leptons=[]
        for il in xrange(0,tree.nl):

            lp4=ROOT.TLorentzVector()
            lp4.SetPtEtaPhiM(tree.l_pt[il],tree.l_eta[il],tree.l_phi[il],tree.l_m[il])
            glp4=ROOT.TLorentzVector()
            glp4.SetPtEtaPhiM(tree.gl_pt[il],tree.gl_eta[il],tree.gl_phi[il],tree.gl_m[il])
            leptons.append( (lp4,glp4) )

        #met
        met=tree.met_pt

        ht=bjets[0][0].Pt()+bjets[1][0].Pt()+leptons[0][0].Pt()+leptons[1][0].Pt()+met
        observablesH['ht_'+evcat].Fill(ht,evWeight)

        dphill=leptons[0][0].DeltaPhi(leptons[1][0])
        if ht>500 : dphill += ROOT.TMath.Pi()
        observablesH['dphill_'+evcat].Fill(dphill,evWeight)

         
    #save results
    fOut=ROOT.TFile.Open(outFileName,'RECREATE')
    for var in observablesH: 
        observablesH[var].Write()
    fOut.Close()

 
"""
Wrapper for when the analysis is run in parallel
"""
def runTopRadiusAnalysisPacked(args):
    try:
        fileNames,outFileName=args
        runTopRadiusAnalysis(fileNames,outFileName)
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
            pool.map(runTopRadiusAnalysisPacked,tasklist)
        else:
            for fileName,outFileName in tasklist:
                runTopRadiusAnalysis(fileName,outFileName)
    else:
        cmsswBase=os.environ['CMSSW_BASE']
        for fileName,_ in tasklist:
            localRun='python %s/src/TopLJets2015/TopAnalysis/scripts/runTopRadiusAnalysis.py -i %s -o %s -q local'%(cmsswBase,fileName,opt.output)
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
                          default='/afs/cern.ch/user/p/psilva/work/TopRadius',
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
