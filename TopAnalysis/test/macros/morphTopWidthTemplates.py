#!/usr/bin/env python

import ROOT

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

ws=ROOT.RooWorkspace('sigws')
ws.factory('x[0,300]')

widths=[1.0,2.0,3.0,4.0]
masses=[(169.5,'~/work/TopWidth_era2015/analysis/plots/syst_plotter.root','t#bar{t} m=169.5'),
        (172.5,'~/work/TopWidth_era2015/analysis/plots/plotter.root',     't#bar{t}'),
        (175.5,'~/work/TopWidth_era2015/analysis/plots/syst_plotter.root','t#bar{t} m=175.5')]

widthDim=ROOT.RooBinning(len(widths)-1,0,len(widths)-1)
massDim=ROOT.RooBinning(len(masses)-1,0,len(masses)-1)
refGrid=ROOT.RooMomentMorphND.Grid(widthDim,massDim)

for ch in ['EM']:
    for cat in ['2b']:
        for lbCat in ['highpt']:
            for imass in xrange(0,len(masses)):
                mass,url,proc=masses[imass]

                fIn=ROOT.TFile.Open(url)
                for iwid in xrange(0,len(widths)):

                    #get histogram and convert to a PDF
                    dirname='%s%s%s_mlb_%3.1fw'%(lbCat,ch,cat,widths[iwid])
                    h=fIn.Get(dirname+"/"+dirname+'_'+proc)
                    name='%s_m%d'%(dirname,int(10*mass))
                    data=ROOT.RooDataHist(name,name,ROOT.RooArgList(ws.var("x")),h)
                    pdf=ROOT.RooHistPdf(name+"_pdf",name+"_pdf",ROOT.RooArgSet(ws.var("x")),data)
                    getattr(ws,'import')(pdf,ROOT.RooCmdArg())

                    #add pdf to the grid
                    print 'Adding',pdf.GetName(),'@ (',iwid,imass,')'
                    refGrid.addPdf(ws.pdf(pdf.GetName()),iwid,imass)
                fIn.Close()

ws.factory('alpha[0,3]')
ws.factory('beta[0,3]')
pdf=ROOT.RooMomentMorphND('widmorphpdf','widmorphpdf',
                          ROOT.RooArgList( ws.var('alpha'), ws.var('beta') ),
                          ROOT.RooArgList( ws.var('x') ),
                          refGrid,
                          ROOT.RooMomentMorphND.Linear)
pdf.useHorizontalMorphing(False)
getattr(ws,'import')(pdf,ROOT.RooCmdArg())

dpdfdalpha=pdf.derivative(ws.var('alpha'))
dpdfdbeta=pdf.derivative(ws.var('beta'))


ws.var('alpha').setVal(0.0)
ws.var('beta').setVal(1.0)

outfile=ROOT.TFile("test.root","RECREATE")
outfile.cd()
ws.Write()
outfile.Close()

lsSensitivities={'alpha':ROOT.TGraph(),
                 'beta':ROOT.TGraph()}
for x in xrange(0,330,30):
    ws.var('x').setVal(x)
    f=pdf.getVal()
    dfda=dpdfdalpha.getVal()
    dfdb=dpdfdbeta.getVal()
    if f==0 : continue
    np=lsSensitivities['alpha'].GetN()
    lsSensitivities['alpha'].SetPoint(np,x,(1/f)*(dfda**2))
    lsSensitivities['beta'].SetPoint(np,x,(1/f)*(dfdb**2))

c= ROOT.TCanvas("c","c",500,500)
frame = ws.var("x").frame()
ws.pdf("highptEM2b_mlb_1.0w_m1725_pdf").plotOn(frame, ROOT.RooFit.LineColor(ROOT.kGray),   ROOT.RooFit.LineStyle(ROOT.kSolid))
ws.pdf('widmorphpdf').plotOn(frame,ROOT.RooFit.LineStyle(ROOT.kDashed))
frame.Draw()
lsSensitivities['alpha'].SetLineColor(ROOT.kRed)
lsSensitivities['alpha'].SetLineWidth(2)
lsSensitivities['alpha'].Draw('c')
lsSensitivities['beta'].SetLineColor(ROOT.kAzure+3)
lsSensitivities['beta'].SetLineWidth(2)
lsSensitivities['beta'].Draw('c')
c.Modified()
c.Update()
raw_input()
