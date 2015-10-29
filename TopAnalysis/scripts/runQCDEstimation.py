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
    parser.add_option('--norm',         dest='norm',         help='SIP3d cut to normalized QCD',                default=4,          type=float)
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
    for sel in ['1j0t','1j1t','2j0t','2j1t','2j2t','3j0t','3j1t','3j2t','4j0t','4j1t','4j2t']:

        #data in the sideband
        dataNonIso, sumMCNonIso = getTemplates(fIn=fNonIso, dist='lsip3d_%s'%sel, tag='noniso',rebin=True)
        dataNonIso.Add(sumMCNonIso,-1)
        dataNonIso.SetTitle('QCD multijets (data)')
       

        #data in the signal region
        dataIso,    sumMCIso    = getTemplates(fIn=fIso,    dist='lsip3d_%s'%sel, tag='iso',rebin=True)
        dataIso.SetTitle('data')
        sumMCIso.SetTitle('other processes')

        #normalized QCD template above the cut in SIP3d
        xbin=dataNonIso.GetXaxis().FindBin(opt.norm)
        nxbins=dataNonIso.GetNbinsX()        
        niso=dataIso.Integral(xbin,nxbins)
        nmciso=sumMCIso.Integral(xbin,nxbins)
        nnoniso=dataNonIso.Integral(xbin,nxbins)
        nonIsoTemplateSF[sel]=(ROOT.TMath.Max(niso-nmciso,0.)/nnoniso,nnoniso)

        #free mem
        dataNonIso.Delete()
        sumMCNonIso.Delete()
        dataIso.Delete()
        sumMCIso.Delete()

    #combined categories
    for combCat in ['1j','2j','3j','4j']:
        totalIni,totalFinal=0,0
        for btags in xrange(0,3):
            key='%s%dt'%(combCat,btags)
            if key in nonIsoTemplateSF:
                totalIni   += nonIsoTemplateSF[key][1]
                totalFinal += nonIsoTemplateSF[key][0]*nonIsoTemplateSF[key][1]
        nonIsoTemplateSF[combCat]=(totalFinal/totalIni,totalIni)

    fIso.Close()

    #produce the QCD templates
    fOut=ROOT.TFile.Open('%s/Data_QCDMultijets.root'%(opt.outdir),'RECREATE')
    for distKey in fNonIso.GetListOfKeys():
        dist=distKey.GetName()
        if 'iso' in dist or 'ratevsrun' in dist: continue
        dataNonIso, sumMCNonIso = getTemplates(fIn=fNonIso, dist=dist, tag=dist)
        dataNonIso.Add(sumMCNonIso,-1)        
        dataNonIso.SetTitle('QCD multijets (data)')
        category=dist.split('_')[-1]
        try:
            dataNonIso.Scale( nonIsoTemplateSF[category][0] )
        except:            
            if 'njetsnbtags' in dist:
                for xbin in xrange(1,dataNonIso.GetNbinsX()+1):
                    val,unc=dataNonIso.GetBinContent(xbin),dataNonIso.GetBinError(xbin)
                    njets=int((xbin-1)/3)+1
                    nbtags=int(xbin%3)
                    key='%dj%dt'%(njets,nbtags)
                    if not key in nonIsoTemplateSF: continue
                    scale=nonIsoTemplateSF[key][0]
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
