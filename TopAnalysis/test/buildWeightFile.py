#!/usr/bin/env python

import ROOT
from runWeightCreationLoop import TUNES,MAXEVENTS

BRs=[
    (511, 'B^{0}',        0.23845, 0.1043, 0.1033, 0.0028),
    (521, 'B^{+}',        0.25579, 0.1129, 0.1099, 0.0028),
    (531, 'B^{0}_{s}',    0.21920, 0.0930, 0.0960, 0.008),
    (5122, '#Lambda_{b}', 0.17870, 0.0770, 0.104, 0.022)
    ]

"""
Interpolate extremes and then derive the weights based on a 2nd order spline for the remaining nodes
"""
def smoothWeights(gr):

    #interpolate for low xb
    gr.Fit('pol9','QR+','',0,0.7)
    lowxb=gr.GetFunction('pol9')

    #flatten tail for xb>1
    gr.Fit('pol0','QR+','',1.03,2)
    highxb=gr.GetFunction('pol0')

    smoothExtremesGr=ROOT.TGraph()
    x,y=ROOT.Double(0),ROOT.Double(0)
    for i in xrange(0,gr.GetN()):
        gr.GetPoint(i,x,y)
        if x<0.55:
            smoothExtremesGr.SetPoint(i,x,lowxb.Eval(x))
        elif x>1.03:
            smoothExtremesGr.SetPoint(i,x,highxb.Eval(x))
        else:
            smoothExtremesGr.SetPoint(i,x,y)

    #smooth the weights
    tSpline=ROOT.TMVA.TSpline2("spline",smoothExtremesGr)
    smoothGr=ROOT.TGraph()
    for i in xrange(0,1000):
        x=2.0*i/1000
        smoothGr.SetPoint(i,x,tSpline.Eval(x))
    return smoothGr


def main():

    outf='${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/era2016/bfragweights.root'
    
    #derive the weights
    xb={}
    for tag,_,_ in TUNES:
        fName='LJets2015/2016/bfrag/xb_%s_numEvent%d.root'%(tag,MAXEVENTS) if MAXEVENTS>0 else 'LJets2015/2016/bfrag/xb_%s.root'%tag
        fIn=ROOT.TFile.Open(fName)
        xb[tag]=fIn.Get('bfragAnalysis/xb_semilepinc').Clone(tag)
        #xb[tag].Rebin(2);
        xb[tag].SetDirectory(0)
        fIn.Close()

    #save to file
    fOut=ROOT.TFile.Open(outf,'RECREATE')
    #for tag in ['up','uup','uuup','central','ccentral','cccentral','down','ddown','dddown']:
    #for tag in [ ('sup','BL',0.865), ('scentral','BL',0.854), ('sdown','BL',0.845) ]
    for tag,frag,param in TUNES:
    #for tag in ['up','central','down']:
    #for tag in ['up','central','down','Peterson']:
        if tag == 'cuetp8m2t4': continue
        print tag
        xb[tag].Scale(1./xb[tag].Integral())
        xb['cuetp8m2t4'].Scale(1./xb['cuetp8m2t4'].Integral())
        xb[tag].Divide(xb['cuetp8m2t4'])
        #xb[tag].Smooth()
        gr=ROOT.TGraphErrors(xb[tag])
        gr.SetMarkerStyle(20)

        sgr=smoothWeights(gr)
        sgr.SetName(tag+'Frag')
        sgr.SetLineColor(ROOT.kBlue)
        sgr.Write()

        #gr.Draw('ap')
        #sgr.Draw('l')
        #raw_input()

    #semi-leptonic BRs
    semilepbrUp=ROOT.TGraph()
    semilepbrUp.SetName("semilepbrUp")
    semilepbrDown=ROOT.TGraph()
    semilepbrDown.SetName("semilepbrDown")

    for entry in BRs:

        i=semilepbrUp.GetN()
        pid,_,py8inc,py8exc,pdg,pdgUnc=entry

        brUp   = py8inc*(1+ROOT.TMath.Max((pdg+pdgUnc)-py8exc,0.)/py8exc)
        semilepbrUp.SetPoint(i,     pid*(-1), (1-brUp)/(1-py8inc))
        semilepbrUp.SetPoint(i+1,   pid,      brUp/py8inc)

        brDown = py8inc*(1-ROOT.TMath.Max(py8exc-(pdg-pdgUnc),0.)/py8exc)
        semilepbrDown.SetPoint(i,   pid*(-1), (1-brDown)/(1-py8inc))
        semilepbrDown.SetPoint(i+1, pid,      brDown/py8inc)

    semilepbrUp.Write()
    semilepbrDown.Write()
        

    fOut.Close()
    print 'Fragmentation been saved to',outf

if __name__ == "__main__":
    main()
