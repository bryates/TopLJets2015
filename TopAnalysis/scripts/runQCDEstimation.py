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
    parser.add_option('--norm',         dest='norm',         help='SIP3d cut to normalized QCD',                default=3.0,        type=float)
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
    qcdNormUnc={}
    #for sel in ['1j0t','1j1t','2j0t','2j1t','2j2t','3j0t','3j1t','3j2t','4j0t','4j1t','4j2t']:  #too aggressive given QCD is small after 2t
    for sel in ['1j','2j','3j','4j']: 

        #data in the sideband
        dataNonIso, sumMCNonIso = getTemplates(fIn=fNonIso, dist='lsip3d_%s'%sel, tag='noniso',rebin=True)

        #data in the signal region
        dataIso,    sumMCIso    = getTemplates(fIn=fIso,    dist='lsip3d_%s'%sel, tag='iso',rebin=True)

        #normalized QCD template above the cut in SIP3d (include overflow)
        xbin=dataNonIso.GetXaxis().FindBin(opt.norm)
        nxbins=dataNonIso.GetNbinsX()+1        

        #signal region
        niso=dataIso.Integral(xbin,nxbins)
        nmciso=sumMCIso.Integral(xbin,nxbins)
        niso_max,niso_cen,niso_min=niso,ROOT.TMath.Max(niso-nmciso,0.),ROOT.TMath.Max(niso-2*nmciso,0.)

        #control region
        nmcnoniso=sumMCNonIso.Integral(xbin,nxbins)
        nnoniso=dataNonIso.Integral(xbin,nxbins)
        nnoniso_max,nnoniso_cen,nnoniso_min=nnoniso,ROOT.TMath.Max(nnoniso-nmcnoniso,0.),ROOT.TMath.Max(nnoniso-2*nmcnoniso,0.)

        #scale factors to apply and relative uncertainty (maximised)
        nonIsoTemplateSF[sel]=niso_cen/nnoniso_cen
        qcdNormUnc[sel]=( (niso_min/nnoniso_max)/nonIsoTemplateSF[sel],
                          (niso_max/nnoniso_min)/nonIsoTemplateSF[sel] )

        #free mem
        dataNonIso.Delete()
        sumMCNonIso.Delete()
        dataIso.Delete()
        sumMCIso.Delete()


    #combined categories
    #for combCat in ['1j','2j','3j','4j']:
    #    totalIni,totalFinal=0,0
    #    for btags in xrange(0,3):
    #        key='%s%dt'%(combCat,btags)
    #        if key in nonIsoTemplateSF:
    #            totalIni   += nonIsoTemplateSF[key][1]
    #            totalFinal += nonIsoTemplateSF[key][0]*nonIsoTemplateSF[key][1]
    #    nonIsoTemplateSF[combCat]=(totalFinal/totalIni,totalIni)

    fIso.Close()

    #produce the QCD templates
    fOut=ROOT.TFile.Open('%s/Data_QCDMultijets.root'%(opt.outdir),'RECREATE')
    for distKey in fNonIso.GetListOfKeys():
        dist=distKey.GetName()

        if 'iso' in dist or 'ratevsrun' in dist: continue
        category=dist.split('_')[-1][:2]    
        dataNonIso, sumMCNonIso = getTemplates(fIn=fNonIso, dist=dist, tag=dist)

        #do not subtract anything in the CR
        dataNonIsoUp=dataNonIso.Clone('%s_QCD%sUp'%(dist,category))
        dataNonIsoUp.SetTitle('QCD multijets (data) - up')

        #subtract 2xMC in the CR
        dataNonIsoDown=dataNonIso.Clone('%s_QCD%sDown'%(dist,category))
        dataNonIsoDown.Add(sumMCNonIso,-2)       
        dataNonIsoDown.SetTitle('QCD multijets (data) - down')

        #subtract 1xMC in the CR
        dataNonIso.Add(sumMCNonIso,-1)  
        dataNonIso.SetName(dist)
        dataNonIso.SetTitle('QCD multijets (data)')

        try:
            dataNonIso.Scale( nonIsoTemplateSF[category] )
            totalQCD=dataNonIso.Integral()
            if totalQCD>0:
                for xbin in xrange(0,dataNonIsoUp.GetNbinsX()):
                    val=dataNonIsoUp.GetBinContent(xbin+1)
                    if val<0 : dataNonIsoUp.SetBinContent(xbin+1,1e-5)
                    val=dataNonIsoDown.GetBinContent(xbin+1)
                    if val<0 : dataNonIsoDown.SetBinContent(xbin+1,1e-5)
                dataNonIsoUp.Scale(totalQCD/dataNonIsoUp.Integral())
                dataNonIsoDown.Scale(totalQCD/dataNonIsoDown.Integral())
        except:            
            print 'unable to normalize ',dist,' with category=',category
        sumMCNonIso.Delete()

        #dump to file
        fOut.cd()
        for h in [dataNonIso,dataNonIsoUp,dataNonIsoDown]:
            h.SetDirectory(fOut)
            h.Write()

        
    #dump to file    
    print 'QCD scale factors CR->SR'
    print nonIsoTemplateSF
    print 'QCD normalization uncertainties'
    print qcdNormUnc
    cachefile=open('%s/.qcdscalefactors.pck'%opt.outdir,'w')
    pickle.dump(nonIsoTemplateSF, cachefile, pickle.HIGHEST_PROTOCOL)
    pickle.dump(qcdNormUnc, cachefile, pickle.HIGHEST_PROTOCOL)
    cachefile.close()

    #all done
    print 'QCD templates and systematic variations stored in %s'%fOut.GetName()
    fOut.Close()
    fNonIso.Close()

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
