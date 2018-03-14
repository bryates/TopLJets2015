#!/usr/bin/env python

import optparse
import os,sys
import json
import ROOT
from subprocess import Popen, PIPE
from collections import OrderedDict

def sampleLoop(inDir,sample,HiForest,weightCounter,wgtCounter,labelH):
   for f in os.listdir('eos/cms/%s/%s' % (inDir,sample ) ):
       fIn=ROOT.TFile.Open('eos/cms/%s/%s/%s' % (inDir,sample,f ) )
       if not HiForest:
           if wgtCounter is None:
               try:
                   wgtCounter=fIn.Get('weightCounter/Event_weight').Clone('genwgts')
                   weightCounter=True
               except:
                   weightCounter=False
               try:
                   if not weightCounter:
                       wgtCounter=fIn.Get('analysis/fidcounter0').Clone('genwgts')
               except:
                   print 'Check eos/cms/%s/%s/%s probably corrupted?' % (inDir,sample,f )
                   continue
               wgtCounter.SetDirectory(0)
               wgtCounter.Reset('ICE')
           labelH=fIn.Get('analysis/generator_initrwgt')
           if labelH : labelH.SetDirectory(0)                
           #if 'MC13TeV' in sample:
               #tree=fIn.Get('analysis/data')
               #Nall=tree.Draw('ttbar_w[0]','','goff')
               #Nneg=tree.Draw('ttbar_w[0]','ttbar_w[0]<0','goff')
               #Nnet=Nall - 2*Nneg
               #print Nnet
           if weightCounter:
               wgtCounter.Add(fIn.Get('weightCounter/Event_weight'))
           else:
               wgtCounter.Add(fIn.Get('analysis/fidcounter0'))
       else:
           if wgtCounter is None:
               wgtCounter=ROOT.TH1F('genwgts','genwgts',500,0,500)
               wgtCounter.SetDirectory(0)
           hiTree=fIn.Get('hiEvtAnalyzer/HiTree')
           for i in xrange(0,hiTree.GetEntriesFast()):
               hiTree.GetEntry(i)
               try:
                   ttbar_w=getattr(hiTree,'ttbar_w')
                   for ibin in xrange(0,ttbar_w.size()):
                       wgtCounter.Fill(ibin,ttbar_w[ibin])
                   if ttbar_w.size()==0: raise ValueError('simple count required')
               except:
                   wgtCounter.Fill(0,1)
       fIn.Close()

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',       dest='inDir',       help='input directory with files',   default='/store/cmst3/user/psilva/LJets2015/5736a2c',        type='string')
    parser.add_option('-e', '--extDir',      dest='extDir',      help='ext directory with files',     default='/store/group/phys_top/byates/ext/',        type='string')
    parser.add_option(      '--HiForest',    dest='HiForest',    help='flag if these are HiForest',   default=False, action='store_true')
    parser.add_option('-o', '--output',      dest='cache',       help='output file',                  default='data/era2016/genweights.root',                      type='string')
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default='data/era2016/samples.json',              type='string')
    (opt, args) = parser.parse_args()

    #read list of samples
    jsonFile = open(opt.json,'r')
    samplesList=json.load(jsonFile,encoding='utf-8',object_pairs_hook=OrderedDict).items()
    #samplesList=list(reversed(samplesList))
    jsonFile.close()

    #mount locally EOS
    eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'
    Popen([eos_cmd, ' -b fuse mount', 'eos'],stdout=PIPE).communicate()

    #loop over samples available
    genweights={}
    for sample in os.listdir('eos/cms/%s' % opt.inDir):

        #sum weight generator level weights
        wgtCounter=None
        weightCounter=False
        labelH=None
        Nnet=None
        #sampleLoop(opt.inDir,sample,opt.HiForest,weightCounter,wgtCounter,labelH)
        for f in os.listdir('eos/cms/%s/%s' % (opt.inDir,sample ) ):
            fIn=ROOT.TFile.Open('eos/cms/%s/%s/%s' % (opt.inDir,sample,f ) )
            if not opt.HiForest:
                if wgtCounter is None:
                    try:
                        wgtCounter=fIn.Get('weightCounter/Event_weight').Clone('genwgts')
                        weightCounter=True
                    except:
                        weightCounter=False
                    try:
                        if not weightCounter:
                            wgtCounter=fIn.Get('analysis/fidcounter0').Clone('genwgts')
                    except:
                        print 'Check eos/cms/%s/%s/%s probably corrupted?' % (opt.inDir,sample,f )
                        continue
                    wgtCounter.SetDirectory(0)
                    wgtCounter.Reset('ICE')
                labelH=fIn.Get('analysis/generator_initrwgt')
                if labelH : labelH.SetDirectory(0)                
                #if 'MC13TeV' in sample:
                    #tree=fIn.Get('analysis/data')
                    #Nall=tree.Draw('ttbar_w[0]','','goff')
                    #Nneg=tree.Draw('ttbar_w[0]','ttbar_w[0]<0','goff')
                    #Nnet=Nall - 2*Nneg
                    #print Nnet
                if weightCounter:
                    wgtCounter.Add(fIn.Get('weightCounter/Event_weight'))
                else:
                    wgtCounter.Add(fIn.Get('analysis/fidcounter0'))
            else:
                if wgtCounter is None:
                    wgtCounter=ROOT.TH1F('genwgts','genwgts',500,0,500)
                    wgtCounter.SetDirectory(0)
                hiTree=fIn.Get('hiEvtAnalyzer/HiTree')
                for i in xrange(0,hiTree.GetEntriesFast()):
                    hiTree.GetEntry(i)
                    try:
                        ttbar_w=getattr(hiTree,'ttbar_w')
                        for ibin in xrange(0,ttbar_w.size()):
                            wgtCounter.Fill(ibin,ttbar_w[ibin])
                        if ttbar_w.size()==0: raise ValueError('simple count required')
                    except:
                        wgtCounter.Fill(0,1)
            fIn.Close()
        if os.path.isdir('eos/cms/%s/%s' % (opt.extDir,sample)):
            for f in os.listdir('eos/cms/%s/%s' % (opt.inDir,sample ) ):
                fIn=ROOT.TFile.Open('eos/cms/%s/%s/%s' % (opt.inDir,sample,f ) )
                #sampleLoop(opt.extDir,sample,opt.HiForest,weightCounter,wgtCounter,labelH)
                if wgtCounter is None:
                    try:
                        wgtCounter=fIn.Get('weightCounter/Event_weight').Clone('genwgts')
                        weightCounter=True
                    except:
                        weightCounter=False
                    try:
                        if not weightCounter:
                            wgtCounter=fIn.Get('analysis/fidcounter0').Clone('genwgts')
                    except:
                        print 'Check eos/cms/%s/%s/%s probably corrupted?' % (opt.extDir,sample,f )
                        continue
                    wgtCounter.SetDirectory(0)
                    wgtCounter.Reset('ICE')
                labelH=fIn.Get('analysis/generator_initrwgt')
                if labelH : labelH.SetDirectory(0)                
                #if 'MC13TeV' in sample:
                    #tree=fIn.Get('analysis/data')
                    #Nall=tree.Draw('ttbar_w[0]','','goff')
                    #Nneg=tree.Draw('ttbar_w[0]','ttbar_w[0]<0','goff')
                    #Nnet=Nall - 2*Nneg
                    #print Nnet
                if weightCounter:
                    wgtCounter.Add(fIn.Get('weightCounter/Event_weight'))
                else:
                    wgtCounter.Add(fIn.Get('analysis/fidcounter0'))
                fIn.Close()

        if wgtCounter is None: continue
        if labelH:
            for xbin in range(1,labelH.GetNbinsX()):
                label=labelH.GetXaxis().GetBinLabel(xbin)
                for tkn in ['<','>',' ','\"','/','weight','=']: label=label.replace(tkn,'')
                wgtCounter.GetXaxis().SetBinLabel(xbin+1,label)

        #invert to set normalization
        print sample,' initial sum of weights=',wgtCounter.GetBinContent(1)
        #print Nnet
        for xbin in xrange(1,wgtCounter.GetNbinsX()+1):
            val=wgtCounter.GetBinContent(xbin)
            if val==0: continue
            wgtCounter.SetBinContent(xbin,1./val)
            wgtCounter.SetBinError(xbin,0.)
        for tag,samples in samplesList:
            if tag[0] is None: continue
            if sample in tag: print tag,' xsec=',samples[0]
            if sample in tag: wgtCounter.SetBinContent(2, samples[0])
       
        genweights[sample]=wgtCounter

    #unmount locally EOS
    Popen([eos_cmd, ' -b fuse umount', 'eos'],stdout=PIPE).communicate()

    #dump to ROOT file    
    cachefile=ROOT.TFile.Open(opt.cache,'RECREATE')
    for sample in genweights:
        genweights[sample].SetDirectory(cachefile)
        genweights[sample].Write(sample)
    cachefile.Close()
    print 'Produced normalization cache @ %s'%opt.cache

    #all done here
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
