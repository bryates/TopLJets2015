import os
import sys
import optparse
import ROOT
import pickle
import json
from TopLJets2015.TopAnalysis.storeTools import *

"""
Get data and sum of MC from file
"""
def getTemplates(fIn,dist,tag,rebin=False):
    data,sumMC=None,None
    for key in fIn.Get(dist).GetListOfKeys():
        keyName=key.GetName()
        if 'Multijets' in keyName : continue
        if 'Graph' in keyName : continue
        h=fIn.Get('%s/%s'%(dist,keyName))
        if keyName==dist:
            data=h.Clone('data_'+tag)
            data.SetDirectory(0)
        else:
            if sumMC:
                sumMC.Add(h)
            else:
                sumMC=h.Clone('mcsum_'+tag)
                sumMC.SetDirectory(0)

    for xbin in xrange(1,sumMC.GetNbinsX()+1):
        val=sumMC.GetBinContent(xbin)
        if val==0:
            sumMC.SetBinContent(xbin,1e-3)

    #if rebin:
    #    data.Rebin()
    #    sumMC.Rebin()

    return data,sumMC


"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--iso',          dest='iso',          help='plotter file with the iso selection',        default=None,       type='string')
    parser.add_option('--noniso',       dest='noniso',       help='plotter file with the non iso selection',    default=None,       type='string')
    parser.add_option('--out',          dest='outdir',       help='output directory',                           default='./',       type='string')
    parser.add_option('--norm',         dest='norm',         help='distribution to be used for normalization',  default='metpt',    type='string')
    (opt, args) = parser.parse_args()

    #prepare output
    os.system('mkdir -p %s'%opt.outdir)

    #prepare fitter
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    ROOT.AutoLibraryLoader.enable()
    ROOT.gSystem.Load('libTopLJets2015TopAnalysis.so')
    ROOT.gROOT.LoadMacro('src/TemplatedFitTools.cc+')
    from ROOT import TemplatedFitTools
    fracFitter=ROOT.TemplatedFitTools()
        
    #open inputs
    fNonIso=ROOT.TFile.Open(opt.noniso)
    fIso=ROOT.TFile.Open(opt.iso)

    #perform a fit to a variable of interest
    nonIsoTemplateSF={}
    for sel in ['1j','2j','3j','4j']:

        #data in the sideband
        dataNonIso, sumMCNonIso = getTemplates(fIn=fNonIso, dist='%s_%s'%(opt.norm,sel), tag='noniso',rebin=True)
        dataNonIso.Add(sumMCNonIso,-1)
        dataNonIso.SetTitle('QCD multijets (data)')
        nsideband=dataNonIso.Integral()

        #data in the signal region
        dataIso,    sumMCIso    = getTemplates(fIn=fIso,    dist='%s_%s'%(opt.norm,sel), tag='iso',rebin=True)
        dataIso.SetTitle('data')
        sumMCIso.SetTitle('other processes')

        #initial estimate in signal region = data-sum MC
        nini=dataIso.Integral()-sumMCIso.Integral()
        dataNonIso.Scale(nini/dataNonIso.Integral())

        #perform the fit to adjust initial estimate
        templates=ROOT.TObjArray()
        templates.Add(dataNonIso)
        templates.Add(sumMCIso)
        res=fracFitter.fit(templates,dataIso,0,'%s/%s_%s'%(opt.outdir,opt.norm,sel))
        nonIsoTemplateSF[sel]=(res.sf*nini/nsideband,res.sfUnc*nini/nsideband)

        #free mem
        dataNonIso.Delete()
        sumMCNonIso.Delete()
        dataIso.Delete()
        sumMCIso.Delete()

    fIso.Close()

    #produce the QCD templates
    fOut=ROOT.TFile.Open('%s/Data_QCDMultijets.root'%(opt.outdir),'RECREATE')
    for distKey in fNonIso.GetListOfKeys():
        dist=distKey.GetName()
        if 'iso' in dist or 'ratevsrun' in dist: continue
        dataNonIso, sumMCNonIso = getTemplates(fIn=fNonIso, dist=dist, tag=dist)
        dataNonIso.Add(sumMCNonIso,-1)        
        dataNonIso.SetTitle('QCD multijets (data)')
        sel='inc'
        if '_1j' in dist : sel='1j'
        if '_2j' in dist : sel='2j'
        if '_3j' in dist : sel='3j'
        if '_4j' in dist : sel='4j'
        try:
            dataNonIso.Scale( nonIsoTemplateSF[sel][0] )
        except:            
            if 'catcountSecVtx' in dist:
                for xbin in xrange(1,dataNonIso.GetNbinsX()+1):
                    val,unc=dataNonIso.GetBinContent(xbin),dataNonIso.GetBinError(xbin)
                    njets=int((xbin-1)/21)+1
                    scale=nonIsoTemplateSF['%dj'%njets][0]
                    val=val*scale
                    unc=unc*scale
                    dataNonIso.SetBinContent(xbin,val)
                    dataNonIso.SetBinError(xbin,unc)
            elif 'catcount' in dist:
                for xbin in xrange(1,dataNonIso.GetNbinsX()+1):
                    val,unc=dataNonIso.GetBinContent(xbin),dataNonIso.GetBinError(xbin)
                    njets=int((xbin-1)/3)+1
                    scale=nonIsoTemplateSF['%dj'%njets][0]
                    val=val*scale
                    unc=unc*scale
                    dataNonIso.SetBinContent(xbin,val)
                    dataNonIso.SetBinError(xbin,unc)
            else:
                print 'unable to normalize',dist
        sumMCNonIso.Delete()
        fOut.cd()
        dataNonIso.SetDirectory(fOut)
        dataNonIso.Write(dist)

    #all done
    fOut.Close()
    fNonIso.Close()

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
