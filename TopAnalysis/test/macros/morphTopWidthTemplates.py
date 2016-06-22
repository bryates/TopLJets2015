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

c= ROOT.TCanvas("c","c",500,500)
frame = ws.var("x").frame()

ws.pdf("highptEM2b_mlb_1.0w_m1755_pdf").plotOn(frame, ROOT.RooFit.LineColor(ROOT.kGray),   ROOT.RooFit.LineStyle(ROOT.kSolid))
ws.var('alpha').setVal((1.0-widths[0])/(widths[-1]-widths[0]))
ws.var('beta').setVal((172.5-masses[0][0])/(masses[-1][0]-masses[0][0]))

ws.var('alpha').setVal(0)
ws.var('beta').setVal(2.0)

ws.pdf('widmorphpdf').plotOn(frame,ROOT.RooFit.LineStyle(ROOT.kDashed))

frame.Draw()
c.Modified()
c.Update()
raw_input()
