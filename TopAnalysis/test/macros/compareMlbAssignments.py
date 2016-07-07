#!/usr/bin/env python

import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

fIn=ROOT.TFile.Open('~/work/TopWidth_era2015/analysis/MC13TeV_TTJets.root')

c=ROOT.TCanvas('c','c',500,500)
for ptCat in ['highpt','lowpt']:
    for bCat in ['2b','1b']:
        h=None
        for ch in ['EE','MM','EM']:
            hnew=fIn.Get(ptCat+ch+bCat+'_pairing')
            try:
                h.Add(hnew)
            except:
                h=hnew.Clone(ptCat+bCat)
                h.SetDirectory(0)

        pie=ROOT.TPie(h)
        pie.SetRadius(.3)
        pie.SetLabelsOffset(.01)
        pie.SetLabelFormat("#splitline{%txt}{%perc}")
        for i in xrange(0,3):
            pie.SetEntryLineColor(i,ROOT.kGray+i);
            pie.SetEntryFillColor(i,ROOT.kGray+i);
        pie.Draw("nol <")
        txt=ROOT.TLatex()
        txt.SetNDC()
        txt.SetTextFont(42)
        txt.DrawLatex(0.12,0.95,'#bf{CMS} #it{simulation} 13 TeV')
        catTitle=bCat+', '
        catTitle += 'boosted' if ptCat=='highpt' else 'non-boosted'
        txt.DrawLatex(0.12,0.9, catTitle)
        for ext in ['.png','.pdf']:
            c.SaveAs(ptCat+bCat+ext)

fIn.Close()

