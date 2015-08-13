#!/usr/bin/env python

import math
import ROOT
from ROOT import THStack, TLatex
from ROOT import TCanvas, TPad, TLegend

class Plot:
    """
    A wrapper to store data and MC histograms for comparison
    """

    def __init__(self,name):
        self.name = name
        self.mc = []
        self.data = None
        self.garbageList = []
        self.loadedBaseTools=False

    def info(self):
        print self.name
        print len(self.mc),' mc processes', ' data=', self.data

    def add(self, h, title, color, isData):
        self.garbageList.append(h)
        h.SetTitle(title)

        if "GeV" in h.GetXaxis().GetTitle():
            h.GetYaxis().SetTitle("Events / %3.1f GeV" % h.GetBinWidth(1))
        else:
            h.GetYaxis().SetTitle("Events / %3.1f" % h.GetBinWidth(1))

        if isData:
            h.SetMarkerStyle(20)
            h.SetMarkerColor(color)
            h.SetLineColor(color)
            h.SetLineWidth(2)
            h.SetFillColor(0)
            h.SetFillStyle(0)
            self.data = h
        else:
            h.SetMarkerStyle(1)
            h.SetMarkerColor(color)
            h.SetLineColor(color)
            h.SetLineWidth(1)
            h.SetFillColor(color)
            h.SetFillStyle(1001)
            self.mc.append(h)

    def reset(self):
        for o in self.garbageList: o.Delete()

    def showTable(self, outDir, firstBin=1, lastBin=-1):
        if len(self.mc)==0:
            print '%s is empty' % self.name
            return

        if not self.loadedBaseTools:
            #ROOT.gSystem.Load("libUserCodellvv_fwk.so")
            self.loadedBaseTools=True

        if firstBin<1: firstBin = 1

        f = open(outDir+'/'+self.name+'.dat','w')
        f.write('------------------------------------------\n')
        f.write("Process".ljust(20),)
        f.write("Events after each cut\n")
        f.write('------------------------------------------\n')

        tot ={}
        err = {}
        f.write(' '.ljust(20),)
        try:
            for xbin in xrange(1,self.mc[0].GetXaxis().GetNbins()+1):
                pcut=self.mc[0].GetXaxis().GetBinLabel(xbin)
                f.write(pcut.ljust(40),)
                tot[xbin]=0
                err[xbin]=0
        except:
            pass
        f.write('\n')
        f.write('------------------------------------------\n')

        for h in self.mc:
            pname = h.GetTitle()
            f.write(pname.ljust(20),)

            for xbin in xrange(1,h.GetXaxis().GetNbins()+1):
                itot=h.GetBinContent(xbin)
                ierr=h.GetBinError(xbin)
                #pval=' & %s'%ROOT.toLatexRounded(itot,ierr,-1,True)
                pval=' & %f +/- %f'%(itot,ierr)
                f.write(pval.ljust(40),)
                tot[xbin] = tot[xbin]+itot
                err[xbin] = err[xbin]+ierr*ierr
            f.write('\n')

        f.write('------------------------------------------\n')
        f.write('Total'.ljust(20),)
        for xbin in tot:
            pval=' & %f +/- %f'%(itot,ierr)
            #pval=' & %s'%ROOT.toLatexRounded(tot[xbin],math.sqrt(err[xbin]),-1,True)
            f.write(pval.ljust(40),)
        f.write('\n')

        if self.data is None: return
        f.write('------------------------------------------\n')
        f.write('Data'.ljust(20),)
        for xbin in xrange(1,self.data.GetXaxis().GetNbins()+1):
            itot=self.data.GetBinContent(xbin)
            pval=' & %d'%itot
            f.write(pval.ljust(40))
        f.write('\n')
        f.write('------------------------------------------\n')
        f.close()


    def show(self, outDir, genPseudoData,logY=False):
        if len(self.mc)==0:
            print '%s is empty' % self.name
            return
        canvas = TCanvas('c_'+self.name,'C',600,600)
        canvas.cd()
        t1 = TPad("t1","t1", 0.0, 0.20, 1.0, 1.0)
        t1.SetBottomMargin(0)
        t1.Draw()
        t1.cd()
        t1.SetLogy(logY)
        self.garbageList.append(t1)

        frame = None
        # leg = TLegend(0.15,0.9,0.9,0.95)
        leg = TLegend(0.6,0.7,0.92,0.89)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        nlegCols = 0

        maxY = 1.0
        if self.data is not None:
            leg.AddEntry( self.data, self.data.GetTitle(),'p')
            frame = self.data.Clone('frame')
            self.garbageList.append(frame)
            maxY = self.data.GetMaximum()*1.1
            frame.Reset('ICE')
        elif genPseudoData:
            leg.AddEntry('pseudodata', 'Pseudo-data','p')

        totalMC = None
        stack = THStack('mc','mc')
        for h in self.mc:
            stack.Add(h,'hist')
            leg.AddEntry(h,h.GetTitle(),'f')
            nlegCols = nlegCols+1
            if totalMC is None:
                totalMC = h.Clone('totalmc')
                self.garbageList.append(totalMC)
                totalMC.SetDirectory(0)
            else:
                totalMC.Add(h)

        if totalMC is not None:
            maxY = max(totalMC.GetMaximum(),maxY)
            if frame is None:
                frame = totalMC.Clone('frame')
                frame.Reset('ICE')
                self.garbageList.append(frame)
            if genPseudoData and self.data is None and totalMC.ClassName().find('TH1')>=0:
                self.data=totalMC.Clone('pseudodata')
                self.data.Sumw2()
                self.data.SetLineColor(1)
                self.data.SetMarkerStyle(20)
                self.data.SetMarkerColor(1)
                self.data.SetFillStyle(0)
                self.data.SetFillColor(0)
                if self.name.find('flow')<0:
                    self.data.Reset('ICE')
                    for i in xrange(0,int(totalMC.Integral())): self.data.Fill( totalMC.GetRandom())

        if self.data is not None: nlegCols = nlegCols+1
        if nlegCols == 0:
            print '%s is empty'%self.name
            return 

        frame.GetYaxis().SetRangeUser(1e-2,1.2*maxY)
        frame.SetDirectory(0)
        frame.Draw()
        frame.GetYaxis().SetTitleOffset(1.6)
        stack.Draw('hist same')
        if self.data is not None: self.data.Draw('P same')
        # leg.SetNColumns(nlegCols)
        leg.SetNColumns(2)
        leg.Draw()

        ## Draw CMS Preliminary label
        tlat = TLatex()
        tlat.SetNDC()
        tlat.SetTextFont(62)
        tlat.SetTextSize(0.04)
        tlat.SetTextAlign(31)
#        prelim_text = 'CMS Preliminary, #sqrt{s} = 8 TeV'
        prelim_text = 'CMS Preliminary, #sqrt{s} = 13 TeV'
        tlat.DrawLatex(0.92, 0.95, prelim_text)


        if totalMC is None or self.data is None:
            t1.SetPad(0,0,1,1)
            t1.SetBottomMargin(0.12)
        else:
            canvas.cd()
            t2 = TPad("t2","t2", 0.0, 0.0, 1.0, 0.2)
            self.garbageList.append(t2)
            t2.SetTopMargin(0)
            t2.SetBottomMargin(0.4)
            t2.SetGridy()
            t2.Draw()
            t2.cd()
            ratio = self.data.Clone('ratio')
            self.garbageList.append(ratio)
            ratio.Divide(totalMC)
            ratio.SetDirectory(0)
            ratio.Draw('e1')
            ratio.GetYaxis().SetRangeUser(0.62,1.38)
            ratio.GetYaxis().SetTitle('Data/#SigmaBkg')
            ratio.GetYaxis().SetNdivisions(5)
            ratio.GetYaxis().SetLabelSize(0.15)
            ratio.GetYaxis().SetTitleSize(0.18)
            ratio.GetXaxis().SetLabelSize(0.18)
            ratio.GetXaxis().SetTitleSize(0.2)
            ratio.GetYaxis().SetTitleOffset(0.3)
            ratio.GetXaxis().SetTitleOffset(0.8)

        canvas.cd()
        canvas.Modified()
        canvas.Update()
        for ext in ['pdf','png'] : canvas.SaveAs(outDir+'/'+self.name+'.'+ext)


def customROOTstyle():
    """
    Loads TDR style
    """
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasDefH(928) #Height of canvas
    ROOT.gStyle.SetCanvasDefW(904) #Width of canvas
    ROOT.gStyle.SetCanvasDefX(1320) #POsition on screen
    ROOT.gStyle.SetCanvasDefY(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetPadColor(ROOT.kWhite)
    ROOT.gStyle.SetPadGridX(False)
    ROOT.gStyle.SetPadGridY(False)
    ROOT.gStyle.SetGridColor(0)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridWidth(1)
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetFrameBorderSize(1)
    ROOT.gStyle.SetFrameFillColor(0)
    ROOT.gStyle.SetFrameFillStyle(0)
    ROOT.gStyle.SetFrameLineColor(1)
    ROOT.gStyle.SetFrameLineStyle(1)
    ROOT.gStyle.SetFrameLineWidth(1)
    ROOT.gStyle.SetHistLineColor(1)
    ROOT.gStyle.SetHistLineStyle(0)
    ROOT.gStyle.SetHistLineWidth(1)
    # ROOT.gStyle.SetLegoInnerR(Float_t rad = 0.5)
    # ROOT.gStyle.SetNumberContours(Int_t number = 20)
    ROOT.gStyle.SetEndErrorSize(2)
    # ROOT.gStyle.SetErrorMarker(20)
    ROOT.gStyle.SetErrorX(0.)
    ROOT.gStyle.SetMarkerStyle(20)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetFitFormat("5.4g")
    ROOT.gStyle.SetFuncColor(2)
    ROOT.gStyle.SetFuncStyle(1)
    ROOT.gStyle.SetFuncWidth(1)
    ROOT.gStyle.SetOptDate(0)
    # ROOT.gStyle.SetDateX(Float_t x = 0.01)
    # ROOT.gStyle.SetDateY(Float_t y = 0.01)
    # For the statistics box:
    ROOT.gStyle.SetOptFile(0)
    ROOT.gStyle.SetOptStat(0) # To display the mean and RMS: SetOptStat("mr")
    ROOT.gStyle.SetStatColor(ROOT.kWhite)
    ROOT.gStyle.SetStatFont(42)
    ROOT.gStyle.SetStatFontSize(0.025)
    ROOT.gStyle.SetStatTextColor(1)
    ROOT.gStyle.SetStatFormat("6.4g")
    ROOT.gStyle.SetStatBorderSize(1)
    ROOT.gStyle.SetStatH(0.1)
    ROOT.gStyle.SetStatW(0.15)
    # ROOT.gStyle.SetStatStyle(Style_t style = 1001)
    # ROOT.gStyle.SetStatX(Float_t x = 0)
    # ROOT.gStyle.SetStatY(Float_t y = 0)
    # Margins:
    ROOT.gStyle.SetPadTopMargin(0.07)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.17)
    ROOT.gStyle.SetPadRightMargin(0.03)
    # For the Global title:
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleColor(1)
    ROOT.gStyle.SetTitleTextColor(1)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleFontSize(0.065)
    ROOT.gStyle.SetTitleH(0.07) # Set the height of the title box
    ROOT.gStyle.SetTitleW(0.80) # Set the width of the title box
    ROOT.gStyle.SetTitleX(0.15) # Set the position of the title box
    ROOT.gStyle.SetTitleY(1.00) # Set the position of the title box
    # ROOT.gStyle.SetTitleStyle(Style_t style = 1001)
    ROOT.gStyle.SetTitleBorderSize(1)
    # For the axis titles:
    ROOT.gStyle.SetTitleColor(1, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetTitleSize(0.06, "XYZ")
    # ROOT.gStyle.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
    # ROOT.gStyle.SetTitleYSize(Float_t size = 0.02)
    ROOT.gStyle.SetTitleXOffset(0.95)
    ROOT.gStyle.SetTitleYOffset(1.3)
    # ROOT.gStyle.SetTitleOffset(1.1, "Y") # Another way to set the Offset
    # For the axis labels:
    ROOT.gStyle.SetLabelColor(1, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
    ROOT.gStyle.SetLabelSize(0.044, "XYZ")
    # For the axis:
    ROOT.gStyle.SetAxisColor(1, "XYZ")
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetTickLength(0.03, "XYZ")
    ROOT.gStyle.SetNdivisions(510, "XYZ")
    ROOT.gStyle.SetPadTickX(1) # To get tick marks on the opposite side of the frame
    ROOT.gStyle.SetPadTickY(1)
    # Change for log plots:
    ROOT.gStyle.SetOptLogx(0)
    ROOT.gStyle.SetOptLogy(0)
    ROOT.gStyle.SetOptLogz(0)
    # Postscript options:
    ROOT.gStyle.SetPaperSize(20.,20.)


