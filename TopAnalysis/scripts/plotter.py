import optparse
import os,sys
import json
import ROOT
import math
import pickle
from collections import OrderedDict

from TopLJets2015.TopAnalysis.rounding import *

"""
increments the first and the last bin to show the under- and over-flows
"""
def fixExtremities(h,addOverflow=True,addUnderflow=True):
    if addUnderflow and h.Integral() > h.GetBinContent(0):
        fbin  = h.GetBinContent(0) + h.GetBinContent(1)
	fbine = ROOT.TMath.Sqrt(h.GetBinError(0)**2 + h.GetBinError(1)**2)
	h.SetBinContent(1,fbin)
	h.SetBinError(1,fbine)
	h.SetBinContent(0,0)
	h.SetBinError(0,0)
    if addOverflow and h.Integral() > h.GetBinContent(h.GetNbinsX()+1):
        nbins = h.GetNbinsX();
	fbin  = h.GetBinContent(nbins) + h.GetBinContent(nbins+1)
	fbine = ROOT.TMath.Sqrt(h.GetBinError(nbins)**2  + h.GetBinError(nbins+1)**2)
	h.SetBinContent(nbins,fbin)
	h.SetBinError(nbins,fbine)
	h.SetBinContent(nbins+1,0)
	h.SetBinError(nbins+1,0)


"""
A wrapper to store data and MC histograms for comparison
"""
class Plot(object):

    def __init__(self,name):
        self.name = name
        self.wideCanvas = True if 'ratevsrun' in self.name else False
        self.mc = OrderedDict()
        self.mcsyst = {}
        self.rbList = []
        self.drawUnc = True
        #self.mc = {}
        self.dataH = None
        self.data = None
        self.totalMCUnc = None
        self.totalMCUncShape = None
        self._garbageList = []
        self.plotformats = ['pdf','png']
        self.savelog = False
        #self.noPU = False
        self.ratiorange = (0.76,1.24)
        if "_jpsi" in name or "_meson" in name:
            self.ratiorange = (0.47,1.57)
        self.ratiorange = (0.76,1.24)

    def add(self, h, title, color, isData, isSyst, rbList):
        ## hack to fix impact parameter range (cm -> um) ##
        if "pf_d" in h.GetName() and "significance" not in h.GetName():
            nbins = h.GetNbinsX()
            h.GetXaxis().Set(nbins,0,1000)
            tmptitle = h.GetXaxis().GetTitle()
            tmptitle = tmptitle.replace("cm","#mum")
            h.GetXaxis().SetTitle(tmptitle)
        if "massJPsi" in h.GetName() and "massJPsi_l_" not in h.GetName():
            h.GetXaxis().SetRangeUser(2.5,3.4)
        if "D0_mu_tag_mu" in h.GetName():
            h.GetXaxis().SetTitle("P_{T}(D^{0}_{#mu}+#mu_{tag})/#Sigma p_{T}^{ch}")
        if "JPsioJet" in h.GetName():
            h.GetXaxis().SetTitle("P_{T}(J/#Psi)/#Sigma p_{T}^{ch}")

        h.SetTitle(title)
        if isData:
            try:
                self.dataH.Add(h)
            except:
                self.dataH=h
                self.dataH.SetDirectory(0)
                self.dataH.SetMarkerStyle(20)
                self.dataH.SetMarkerSize(0.9)
                #self.dataH.SetMarkerSize(1.4)
                self.dataH.SetMarkerColor(color)
                self.dataH.SetLineColor(ROOT.kBlack)
                self.dataH.SetLineWidth(2)
                self.dataH.SetFillColor(0)
                self.dataH.SetFillStyle(0)
                self._garbageList.append(h)
        elif isSyst:
            try:
                self.mcsyst[title].Add(h)
            except:
                self.mcsyst[title]=h
                self.mcsyst[title].SetName('%s_%s' % (h.GetName(), title ) )
                self.mcsyst[title].SetDirectory(0)
                self.mcsyst[title].SetMarkerStyle(1)
                self.mcsyst[title].SetMarkerColor(color)
                self.mcsyst[title].SetLineColor(ROOT.kBlack)
                self.mcsyst[title].SetLineWidth(1)
                self.mcsyst[title].SetFillColor(color)
                self.mcsyst[title].SetFillStyle(1001)
                self._garbageList.append(h)
        else:
            try:
                self.mc[title].Add(h)
            except:
                self.mc[title]=h
                self.mc[title].SetName('%s_%s' % (self.mc[title].GetName(), title ) )
                self.mc[title].SetDirectory(0)
                self.mc[title].SetMarkerStyle(1)
                self.mc[title].SetMarkerColor(color)
                self.mc[title].SetLineColor(ROOT.kBlack)
                self.mc[title].SetLineWidth(1)
                self.mc[title].SetFillColor(color)
                self.mc[title].SetFillStyle(1001)
                self._garbageList.append(h)
        self.rbList=rbList

    def finalize(self):
        self.data = convertToPoissonErrorGr(self.dataH)

    def appendTo(self,outUrl):
        outF = ROOT.TFile.Open(outUrl,'UPDATE')
        if not outF.cd(self.name):
            outDir = outF.mkdir(self.name)
            outDir.cd()
        for m in self.mcsyst :
            self.mcsyst[m].Write(self.mcsyst[m].GetName(), ROOT.TObject.kOverwrite)
        for m in self.mc :
            self.mc[m].Write(self.mc[m].GetName(), ROOT.TObject.kOverwrite)
        if self.dataH :
            self.dataH.Write(self.dataH.GetName(), ROOT.TObject.kOverwrite)
        if self.data :
            self.data.Write(self.data.GetName(), ROOT.TObject.kOverwrite)
        outF.Close()

    def reset(self):
        for o in self._garbageList:
            try:
                o.Delete()
            except:
                pass

    def show(self, outDir,lumi,noStack=False,saveTeX=False):

        if len(self.mc)<2 and self.dataH is None:
            print '%s has 0 or 1 MC!' % self.name
            return

        if len(self.mc)>0 and self.mc.values()[0].InheritsFrom('TH2') :
            print 'Skipping TH2'
            return

        cwid=1000 if self.wideCanvas else 500
        c = ROOT.TCanvas('c','c',cwid,500)
        c.SetBottomMargin(0.0)
        c.SetLeftMargin(0.0)
        c.SetTopMargin(0)
        c.SetRightMargin(0.00)
        c.SetFillStyle(4000)
        c.SetFillColor(400)
        c.SetBorderSize(0)
        c.SetLineWidth(0)

        #holds the main plot
        c.cd()
        p1 = None
        if self.dataH:
            p1=ROOT.TPad('p1','p1',0.0,0.2,1.0,1.0) if cwid!=1000 else ROOT.TPad('p1','p1',0.0,0.0,1.0,1.0)
            p1.SetRightMargin(0.05)
            p1.SetLeftMargin(0.12)
            p1.SetTopMargin(0.1)
            p1.SetBottomMargin(0.01)
        else:
            p1=ROOT.TPad('p1','p1',0.0,0.0,1.0,1.0)
            p1.SetRightMargin(0.05)
            p1.SetLeftMargin(0.12)
            p1.SetTopMargin(0.1)
            p1.SetBottomMargin(0.12)
        p1.Draw()

        p1.SetGridx(False)
        p1.SetGridy(False) #True)
        self._garbageList.append(p1)
        p1.cd()

        # legend
        iniy=0.9 if self.wideCanvas or noStack else 0.85
        dy=0.1 if noStack else 0.02
        ndy=len(self.mc) if noStack else max(len(self.mc)-2,0)
        leg = ROOT.TLegend(0.45, iniy-dy*ndy, 0.95, iniy+0.05)

        leg.SetBorderSize(0)
        leg.SetFillStyle(0)        
        leg.SetTextFont(43)
        leg.SetTextSize(12)
        nlegCols = 0

        if self.dataH is not None:
            if self.data is None: self.finalize()
            leg.AddEntry( self.data, self.data.GetTitle(),'p')
            nlegCols += 1
        for h in self.mc:
            
            #compare
            if noStack:
                refH=self.mc.values()[0]
                if refH!=self.mc[h]:
                    chi2=refH.Chi2Test( self.mc[h], 'WW CHI2')
                    pval=refH.Chi2Test( self.mc[h], 'WW')     
                    self.mc[h].SetTitle('#splitline{%s}{#chi^{2}=%3.1f (p-val: %3.3f)}'%(self.mc[h].GetTitle(),chi2,pval))
                else:
                    refH.SetLineWidth(2)

            leg.AddEntry(self.mc[h], self.mc[h].GetTitle(), 'f')
            nlegCols += 1
        if nlegCols ==0 :
            print '%s is empty'%self.name
            return

        if not noStack:
            leg.SetNColumns(ROOT.TMath.Min(nlegCols/2,3))

        # Build the stack to plot from all backgrounds
        totalMC = None
        nominalTTbar=None
        totalMCUnc=None
        stack = ROOT.THStack('mc','mc')
        for h in self.mc:

            if noStack:
                self.mc[h].SetFillStyle(0)
                self.mc[h].SetLineColor(self.mc[h].GetFillColor())
                
            stack.Add(self.mc[h],'hist')
            
            try:
                totalMC.Add(self.mc[h])
            except:
                totalMC = self.mc[h].Clone('totalmc')
                self._garbageList.append(totalMC)
                totalMC.SetDirectory(0)
        for h in self.mc:
            if 't#bar{t}' not in h: continue
            if h=='t#bar{t}':
                nominalTTbar = self.mc[h].Clone('nomttbar')
                self._garbageList.append(nominalTTbar)
                nominalTTbar.SetDirectory(0)

        #systematics
        nominalIntegral = 1.
        #if totalMC and nominalTTbar:# and 'FSR' in h:
        if totalMC and nominalTTbar and len(self.mcsyst)>0:
            nominalIntegral = nominalTTbar.Integral()
            #for hname in self.mcsyst:#.iteritems():
                #if self.mcsyst[hname].Integral()>0: print "%s: %s" % (hname, nominalIntegral/self.mcsyst[hname].Integral())
                #if(self.mcsyst[hname].Integral()>0): self.mcsyst[hname].Scale(nominalIntegral/self.mcsyst[hname].Integral())
            systUp=[0.]
            systDown=[0.]
            for xbin in xrange(1,nominalTTbar.GetNbinsX()+1):
                systUp.append(0.)
                systDown.append(0.)
                for hname in self.mcsyst:#.iteritems():
                    rbIdx=-1
                    mIdx=-1
                    if 'FSR-down' in hname: rbIdx=2
                    if 'FSR-up' in hname: rbIdx=4
                    if 'UEdown' in hname: rbIdx=6
                    if 'UEdown' in hname: rbIdx=8
                    if 'CR' in hname: rbIdx=5
                    #if opt.rbJson is not '': rbIdx=-1
                    if len(self.rbList)<1: rbIdx=-1
                    if 'meson' in self.name and 'mu_tag' in self.name: mIdx=0
                    elif 'meson' in self.name: mIdx=1
                    elif 'jpsi' in self.name: mIdx=2
                    #if totalMCUnc is None:
                        #totalMCUnc = self.mcsyst[hname].Clone()
                    #else: totalMCUnc.Add(self.mcsyst[hname].Clone())
                    #totalMCUnc = self.mcsyst[hname].Clone()
                    diff = self.mcsyst[hname].GetBinContent(xbin) - nominalTTbar.GetBinContent(xbin)
                    if rbIdx > 0 and diff != 0.0:
                        diff = abs(diff)
                        if(self.rbList[mIdx][1][rbIdx] < self.rbList[0][mIdx][0]): diff = -diff
                        print self.rbList[mIdx][1][rbIdx],self.rbList[0][mIdx][0]
                    print self.mcsyst[hname].GetBinContent(xbin),nominalTTbar.GetBinContent(xbin),diff
                    if (diff > 0):
                        systUp[xbin] = math.sqrt(systUp[xbin]**2 + diff**2)
                    else:
                        systDown[xbin] = math.sqrt(systDown[xbin]**2 + diff**2)
            totalMCUnc = totalMC.Clone('totalmcunc')
            self._garbageList.append(totalMCUnc)
            #if(totalMCUnc.Integral()>0): totalMCUnc.Scale(nominalIntegral/totalMCUnc.Integral())
            #totalMCUnc.Add(totalMC) #1) add total MC
            #totalMCUnc.Add(nominalTTbar,-1) #2) subtract nominal MC to add background only
            totalMCUnc.SetDirectory(0)
            totalMCUnc.SetFillColor(0)
            #totalMCUnc.SetLineColor(ROOT.kRed)
            #totalMCUnc.SetMarkerColor(ROOT.kRed)
            totalMCUnc.SetFillColor(ROOT.TColor.GetColor('#99d8c9'))
            totalMCUnc.SetFillColor(ROOT.kBlue-3)
            ROOT.gStyle.SetHatchesLineWidth(1)
            totalMCUnc.SetFillStyle(400)
            totalMCUnc.SetFillStyle(3245)
            for xbin in xrange(1,nominalTTbar.GetNbinsX()+1):
                #scale=nominalIntegral/totalMCUnc.Integral()
                #scale=len(systUp)
                scale=1.
                totalMCUnc.SetBinContent(xbin, totalMCUnc.GetBinContent(xbin) + (systUp[xbin]-systDown[xbin])/2.)
                totalMCUnc.SetBinError(xbin, math.sqrt(totalMCUnc.GetBinError(xbin)**2 + ((systUp[xbin]+systDown[xbin])/2.)**2) / math.sqrt(scale))
            #if(totalMCUnc.Integral()>0): totalMCUnc.Scale(nominalIntegral/totalMCUnc.Integral())
            for hname,h in self.mcsyst.iteritems():
                if(h.Integral>0): h.Scale(nominalIntegral/h.Integral())
            systUpShape=[0.]
            systDownShape=[0.]
            for xbin in xrange(1,nominalTTbar.GetNbinsX()+1):
                systUpShape.append(0.)
                systDownShape.append(0.)
                for hname in self.mcsyst:
                    rbIdx=-1
                    mIdx=-1
                    if 'FSR-down' in hname: rbIdx=2
                    if 'FSR-up' in hname: rbIdx=4
                    if 'UEdown' in hname: rbIdx=6
                    if 'UEdown' in hname: rbIdx=8
                    if 'CR' in hname: rbIdx=5
                    #if opt.rbJson is not '': rbIdx=-1
                    if len(self.rbList)<1: rbIdx=-1
                    if 'meson' in self.name and 'mu_tag' in self.name: mIdx=0
                    elif 'meson' in self.name: mIdx=1
                    elif 'jpsi' in self.name: mIdx=2
                    diff = self.mcsyst[hname].GetBinContent(xbin) - nominalTTbar.GetBinContent(xbin)
                    if rbIdx > 0 and diff != 0.0:
                        diff = abs(diff)
                        if(self.rbList[mIdx][1][rbIdx] < self.rbList[0][mIdx][0]): diff = -diff
                        print self.rbList[mIdx][1][rbIdx],self.rbList[0][mIdx][0]
                    if(diff > 0):
                        systUpShape[xbin] = math.sqrt(systUpShape[xbin]**2 + diff**2)
                    else:
                        systDownShape[xbin] = math.sqrt(systDownShape[xbin]**2 + diff**2)
            totalMCUncShape = totalMC.Clone("totaluncshape")
            self._garbageList.append(totalMCUncShape)
            totalMCUncShape.SetFillColor(ROOT.TColor.GetColor('#d73027'))
            totalMCUncShape.SetFillStyle(3254)
            for xbin in xrange(1,nominalTTbar.GetNbinsX()+1):
                totalMCUncShape.SetBinContent(xbin, totalMCUncShape.GetBinContent(xbin) + (systUpShape[xbin]-systDownShape[xbin])/2.)
                totalMCUncShape.SetBinError(xbin, math.sqrt(totalMCUncShape.GetBinError(xbin)**2 + ((systUpShape[xbin]+systDownShape[xbin])/2.)**2))
            self.totalMCUnc = totalMCUnc
            self.totalMCUncShape = totalMCUncShape

        #test for null plots
        if totalMC :
            if totalMC.Integral()==0:
                if self.dataH is None : return
                if self.dataH.Integral()==0: return
        elif self.dataH is None : return
        elif self.dataH.Integral()==0 : return 


        frame = totalMC.Clone('frame') if totalMC is not None else self.dataH.Clone('frame')
        frame.Reset('ICE')
        if noStack:
            maxY=stack.GetStack().At(0).GetMaximum()/1.25
        elif totalMC:
            maxY = totalMC.GetMaximum() 
            if self.dataH:
                if maxY<self.dataH.GetMaximum():
                    maxY=self.dataH.GetMaximum()
        else:
            maxY=self.dataH.GetMaximum()

        frame.GetYaxis().SetRangeUser(0.1,maxY*1.45)

        frame.SetDirectory(0)
        frame.Reset('ICE')
        self._garbageList.append(frame)
        frame.GetYaxis().SetTitleSize(0.045)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetNoExponent()
        frame.GetYaxis().SetTitleOffset(1.3)
        frame.GetXaxis().SetTitleSize(0.0)
        frame.GetXaxis().SetLabelSize(0.0)
        frame.Draw()
        if totalMC is not None   : 
            if noStack: stack.Draw('nostack same')
            else      : 
               stack.Draw('hist same')
               if len(self.mcsyst)>0:
                   #self.totalMCUnc.Draw("hist same")
                   #self.totalMCUnc.Draw("e2 same")
                   #leg.AddEntry(totalMCUnc, "UE-Down", 'f')
                   #leg.AddEntry(totalMCUnc, "Total unc.", 'f')
                   self.totalMCUncShape.Draw("e2 same")
                   leg.AddEntry(totalMCUncShape, "Total shape unc.", 'f')
        if self.data is not None : self.data.Draw('p')


        leg.Draw()
        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(43)
        txt.SetTextSize(16)
        txt.SetTextAlign(12)
        iniy=0.9 if self.wideCanvas else 0.95
        inix=0.12 if noStack else 0.12
        if lumi<100:
            txt.DrawLatex(inix,iniy,'#bf{CMS} #it{Preliminary} %3.1f pb^{-1} (13 TeV)' % (lumi) )
        else:
            txt.DrawLatex(inix,iniy,'#bf{CMS} #it{Preliminary} %3.1f fb^{-1} (13 TeV)' % (lumi/1000.) )

        #holds the ratio
        c.cd()
        if len(self.mc)>0 and self.dataH:
            p2 = ROOT.TPad('p2','p2',0.0,0.0,1.0,0.2)
            p2.Draw()
            p2.SetBottomMargin(0.4)
            p2.SetRightMargin(0.05)
            p2.SetLeftMargin(0.12)
            p2.SetTopMargin(0.01)
            p2.SetGridx(False)
            p2.SetGridy(True)
            self._garbageList.append(p2)
            p2.cd()
            ratioframe=frame.Clone('ratioframe')
            ratioframe.GetYaxis().SetTitle('Data/MC')
            #ratioframe.GetYaxis().SetTitle('Ratio')
            ratioframe.GetYaxis().SetRangeUser(self.ratiorange[0], self.ratiorange[1])
            self._garbageList.append(frame)
            ratioframe.GetYaxis().SetNdivisions(5)
            ratioframe.GetYaxis().SetLabelSize(0.18)        
            ratioframe.GetYaxis().SetTitleSize(0.2)
            ratioframe.GetYaxis().SetTitleOffset(0.25)
            ratioframe.GetXaxis().SetLabelSize(0.15)
            ratioframe.GetXaxis().SetTitleSize(0.2)
            ratioframe.GetXaxis().SetTitleOffset(0.8)
            ratioframeshape=ratioframe.Clone('ratioframeshape')
            self._garbageList.append(ratioframeshape)
            ratioframeshape.SetFillColor(ROOT.kRed)
            #ratioframeshape.SetMarkerColor(ROOT.kRed)
            #ratioframeshape.SetMarkerStyle(20)
            ratioframeshape.SetFillColor(ROOT.TColor.GetColor('#99d8c9'))
            ratioframeunc=ratioframeshape.Clone('ratioframeunc')
            ratioframeunc.SetFillColor(ROOT.kBlue-3)
            ratioframeunc.SetFillStyle(3245)
            ratioframeshape.SetFillColor(ROOT.TColor.GetColor('#d73027'))
            ratioframeshape.SetFillStyle(3254)
            ratio=self.dataH.Clone('ratio')
            if len(self.mcsyst)>0:
                for xbin in xrange(1,totalMC.GetNbinsX()+1):
                    val=totalMC.GetBinContent(xbin)
                    unc=self.totalMCUnc.GetBinError(xbin)
                    uncs=self.totalMCUncShape.GetBinError(xbin)
                    if val>0:
                        totalMCUnc=ROOT.TMath.Sqrt((unc/val)**2)# + unc**2)
                        totalMCUncShape=ROOT.TMath.Sqrt((uncs/val)**2)# + unc**2)
                        #if self.totalMCUnc.GetBinContent(xbin) > 0:
                            #ratioframeshape.SetBinContent(xbin,self.dataH.GetBinContent(xbin)/self.totalMCUnc.GetBinContent(xbin))
                        ratioframeshape.SetBinContent(xbin,self.totalMCUncShape.GetBinContent(xbin)/val)
                        ratioframeshape.SetBinError(xbin,totalMCUncShape)
                        ratioframeunc.SetBinContent(xbin,self.totalMCUnc.GetBinContent(xbin)/val)
                        ratioframeunc.SetBinError(xbin,totalMCUnc)
                        #ratioframeunc.SetBinError(xbin, ratio.GetBinError(xbin)/self.totalMCUnc.GetBinContent(xbin))
                        #if(ratio.GetBinContent(xbin)>0): ratioframeunc.SetBinError(xbin, math.sqrt(ratio.GetBinError(xbin)**2/ratio.GetBinContent(xbin)**2 + self.totalMCUnc.GetBinError(xbin)**2/self.totalMCUnc.GetBinContent(xbin)**2))
                        #ratioframeunc.SetBinContent(xbin, ratio.GetBinContent(xbin)/self.totalMCUnc.GetBinContent(xbin))
                        #ratioframeunc=self.totalMCUnc.Clone('ratioframunc')
                        #ratioframeunc.Divide(self.totalMCUnc,self.totalMCUnc)
                        #ratioframeunc.SetFillColor(ROOT.kBlue-3)
                        #ratioframeunc.SetFillStyle(3245)


            ratioframe.Draw()
            #if len(self.mcsyst)>0: ratioframeshape.Draw("p same")
            #if len(self.mcsyst)>0: ratioframeunc.Draw("e2 same")
            if len(self.mcsyst)>0: ratioframeshape.Draw("e2 same")

            try:
                #ratio=nominalTTbar.Clone('ratio') #use for purity plots (S/(S+B))
                ratio=self.dataH.Clone('ratio')
                ratio.SetDirectory(0)
                self._garbageList.append(ratio)
                ratio.Divide(totalMC)
                #for xbin in xrange(1,ratio.GetNbinsX()+1):
                    #if totalMC.GetBinContent(xbin) > 0.:
                        #ratio.SetBinError(xbin, ratio.GetBinError(xbin)/totalMC.GetBinContent(xbin))
                        #ratio.SetBinContent(xbin, ratio.GetBinContent(xbin)/totalMC.GetBinContent(xbin))
                    #else:
                        #ratio.SetBinError  (xbin, 0.)
                        #ratio.SetBinContent(xbin, 0.)
                gr=ROOT.TGraphAsymmErrors(ratio)
                gr.SetMarkerStyle(self.data.GetMarkerStyle())
                gr.SetMarkerSize(self.data.GetMarkerSize())
                gr.SetMarkerColor(self.data.GetMarkerColor())
                gr.SetLineColor(self.data.GetLineColor())
                gr.SetLineWidth(self.data.GetLineWidth())
                gr.Draw('p')
            except:
                pass

        #all done
        c.cd()
        c.Modified()
        c.Update()

        #save
        #if self.noPU: 
            #self.name += "_noPU"
        for ext in self.plotformats : c.SaveAs(os.path.join(outDir, self.name+'.'+ext))
        if self.savelog:
            p1.cd()
            frame.GetYaxis().SetRangeUser(1,maxY*50)
            p1.SetLogy()
            c.cd()
            c.Modified()
            c.Update()
            for ext in self.plotformats : c.SaveAs(os.path.join(outDir, self.name+'_log.'+ext))

        if saveTeX : self.convertToTeX(outDir=outDir)


    def convertToTeX(self, outDir):
        if len(self.mc)==0:
            print '%s is empty' % self.name
            return

        f = open(outDir+'/'+self.name+'.dat','w')
        f.write('------------------------------------------\n')
        f.write("Process".ljust(20),)
        f.write("Events after each cut\n")
        f.write('------------------------------------------\n')

        tot ={}
        err = {}
        f.write(' '.ljust(20),)
        try:
            for xbin in xrange(1,self.mc.values()[0].GetXaxis().GetNbins()+1):
                pcut=self.mc.values()[0].GetXaxis().GetBinLabel(xbin)
                f.write(pcut.ljust(40),)
                tot[xbin]=0
                err[xbin]=0
        except:
            pass
        f.write('\n')
        f.write('------------------------------------------\n')

        for pname in self.mc:
            h = self.mc[pname]
            f.write(pname.ljust(20),)

            for xbin in xrange(1,h.GetXaxis().GetNbins()+1):
                itot=h.GetBinContent(xbin)
                ierr=h.GetBinError(xbin)
                pval=' & %s'%toLatexRounded(itot,ierr)
                f.write(pval.ljust(40),)
                tot[xbin] = tot[xbin]+itot
                err[xbin] = err[xbin]+ierr*ierr
            f.write('\n')

        f.write('------------------------------------------\n')
        f.write('Total'.ljust(20),)
        for xbin in tot:
            pval=' & %s'%toLatexRounded(tot[xbin],math.sqrt(err[xbin]))
            f.write(pval.ljust(40),)
        f.write('\n')

        if self.dataH is None: return
        f.write('------------------------------------------\n')
        f.write('Data'.ljust(20),)
        for xbin in xrange(1,self.dataH.GetXaxis().GetNbins()+1):
            itot=self.dataH.GetBinContent(xbin)
            pval=' & %d'%itot
            f.write(pval.ljust(40))
        f.write('\n')
        f.write('------------------------------------------\n')
        f.close()



"""
converts a histogram to a graph with Poisson error bars
"""
def convertToPoissonErrorGr(h):

    htype=h.ClassName()
    if htype.find('TH1')<0 : return None

    #check https://twiki.cern.ch/twiki/bin/view/CMS/PoissonErrorBars
    alpha = 1 - 0.6827;
    grpois = ROOT.TGraphAsymmErrors(h);
    for i in xrange(0,grpois.GetN()+1) :
        N = grpois.GetY()[i]
        if N<200 :
            L = 0
            if N>0 : L = ROOT.Math.gamma_quantile(alpha/2,N,1.)
            U = ROOT.Math.gamma_quantile_c(alpha/2,N+1,1)
            grpois.SetPointEYlow(i, N-L)
            grpois.SetPointEYhigh(i, U-N)
        else:
            grpois.SetPointEYlow(i, math.sqrt(N))
            grpois.SetPointEYhigh(i,math.sqrt(N))
    return grpois


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default=None,              type='string')
    parser.add_option(      '--systJson',    dest='systJson',    help='json with list of systematics',  default=None,              type='string')
    parser.add_option(      '--rbJson',      dest='rbJson',      help='json with list of systematics',  default=None,              type='string')
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default=None,              type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--noStack',     dest='noStack',     help='don\'t stack distributions',     default=False,             action='store_true')
    parser.add_option(      '--saveLog',     dest='saveLog' ,    help='save log versions of the plots', default=False,             action='store_true')
    parser.add_option(      '--saveNorm',    dest='saveNorm' ,   help='save normalised versions of the plots', default=False,             action='store_true')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
    parser.add_option(      '--rebin',       dest='rebin',       help='rebin factor',                   default=1,                 type=int)
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi to print out',              default='data/era2016/lumi.json',        type='string')
    parser.add_option(      '--only',        dest='only',        help='plot only these (csv)',          default='',                type='string')
    parser.add_option(      '--run',         dest='run',         help='plot only in run',               default="BCDEFGH",         type='string')
    parser.add_option(      '--puNormSF',    dest='puNormSF',    help='Use this histogram to correct pu weight normalization', default=None, type='string')
    parser.add_option(      '--procSF',      dest='procSF',      help='Use this to scale a given process component e.g. "W":.wjetscalefactors.pck,"DY":dyscalefactors.pck', default=None, type='string')
    parser.add_option(      '--rbFit',       dest='rbFit' ,      help='Adjust normalization for r_B fit', default=False,            action='store_true')
    (opt, args) = parser.parse_args()

    #read list of samples
    jsonFile = open(opt.json,'r')
    samplesList=json.load(jsonFile,encoding='utf-8',object_pairs_hook=OrderedDict).items()
    #samplesList=list(reversed(samplesList))
    jsonFile.close()

    systs = ['/TRIGGER', '/TRK', '/LEP', '/PU', '/PI', '/JER' ]
    inDir = [opt.inDir]
    for syst in systs:
        for s in ['/up','/down']:
            inDir+=[opt.inDir+s+syst]

    #read list of syst samples
    systSamplesList = []
    if opt.systJson:
        systJsonList = opt.systJson.split(',')
        for jsonPath in systJsonList:
            systJsonFile = open(jsonPath,'r')
            systSamplesList += json.load(systJsonFile,encoding='utf-8').items()
            systJsonFile.close()

    #read list of r_B systs
    rbSamplesList = []
    if opt.rbJson:
        rbJsonList = opt.rbJson.split(',')
        for jsonPath in rbJsonList:
            rbJsonFile = open(jsonPath,'r')
            rbSamplesList += json.load(rbJsonFile,encoding='utf-8').items()
            rbJsonFile.close()

    #read list of lumis
    jsonFile = open(opt.lumi,'r')
    lumiList=json.load(jsonFile,encoding='utf-8')
    lumiList=lumiList["Data13TeV_SingleMuon"]
    jsonFile.close()

    #proc SF
    procSF={}
    if opt.procSF:
        procList=opt.procSF.split(',')
        for newProc in procList:
            proc,cacheUrl=newProc.split(':')
            if not os.path.isfile(cacheUrl) : continue
            cache=open(cacheUrl,'r')
            procSF[proc]=pickle.load(cache)
            cache.close()
            print 'Scale factors added for',proc

    onlyList=opt.only.split(',')

    #read plots 
    plots=OrderedDict()
    #plots={}
    report=''
    if opt.puNormSF.split('_')[-1] not in "BCDEFGH":
        opt.puNormSF += "_"
        opt.puNormSF += opt.run
    print "Using PU SF " + opt.puNormSF
    for slist,isSyst in [(reversed(samplesList),False),(systSamplesList,True)]:
        if slist is None: continue
        for tag,sample in slist:
            if isSyst and "t#bar{t}" not in sample[3]: continue
            if tag[0] is None: continue
            splitRun = lambda x: ["2016" + x[i] for i in range(0, len(x), 1)]
            split_run = splitRun( opt.run )
            if 'Data' in tag:
                if not any(run in tag for run in split_run): continue

            xsec=sample[0]
            isData=sample[1]
            doFlavourSplitting=sample[6]
            subProcs=[(tag,sample[3],sample[4])]
            if doFlavourSplitting:
                subProcs=[]
                for flav in [(1,sample[3]+'+l'),(4,sample[3]+'+c'),(5,sample[3]+'+b',sample[4])]:
                    subProcs.append(('%d_%s'%(flav[0],tag),flav[1],sample[4]+3*len(subProcs)))
            for sp in subProcs:
                for Dir in inDir:

                    tmp=''
                    if '/up/' in Dir:
                        tmp = Dir.split('/')[-1]
                        if tmp not in sp[0]: continue
                        tmp = 'up/' + tmp
                        if 'powheg' not in sp[0]: continue
                        #sp[0]+=tmp
                    if '/down/' in Dir:
                        tmp = Dir.split('/')[-1]
                        if tmp not in sp[0]: continue
                        if 'powheg' not in sp[0]: continue
                        tmp = 'down/' + tmp
                        #sp[0]+=tmp
                    if 'up' in Dir and '_up_' not in sp[0]: continue
                    if 'down' in Dir and '_down_' not in sp[0]: continue
                    if 'up' not in Dir and '_up_' in sp[0]: continue
                    if 'down' not in Dir and '_down_' in sp[0]: continue
                    fIn=ROOT.TFile.Open('%s/%s.root' % ( Dir, sp[0]) )
                    #fIn=ROOT.TFile.Open('%s/%s.root' % ( opt.inDir, sp[0]) )
                    if not fIn : continue
                    print 'Loading file: %s' % sp[0]

                    #fix top pT weighting normalization
                    topPtNorm=1
                    if not isData and "TTJets" in tag:
                        try:
                            topPtBCDEF=fIn.Get("topptwgt_BCDEF")
                            nonTopWgtBCDEF=topPtBCDEF.GetBinContent(1)
                            topWgtBCDEF=topPtBCDEF.GetBinContent(2)
                        except: pass
                        try:
                            topPtGH=fIn.Get("topptwgt_GH")
                            nonTopWgtGH=topPtGH.GetBinContent(1)
                            topWgtGH=topPtGH.GetBinContent(2)
                        except: pass
                        if topWgtBCDEF>0 :
                            topPtNormBCDEF=topWgtBCDEF/nonTopWgtBCDEF
                            report += '%s was scaled by %3.3f for top pT epoch %s reweighting\n' % (sp[0],topPtNormBCDEF,"BCDEF")
                        if topWgtGH>0 :
                            topPtNormGH=topWgtGH/nonTopWgtGH
                            report += '%s was scaled by %3.3f for top pT epoch %s reweighting\n' % (sp[0],topPtNormGH,"GH")

                    #fix pileup weighting normalization
                    puNormSF=1
                    if opt.puNormSF and not isData:
                        if(opt.run == "BCDEFGH"):
                          puBCDEF = opt.puNormSF.split("_")
                          if puBCDEF[-1] in opt.run:
                              puBCDEF[-1] = "BCDEF"
                          else:
                              puBCDEF.append("BCDEF")
                          puBCDEF = "_".join(puBCDEF)
                          try:
                              puCorrHBCDEF=fIn.Get(puBCDEF)
                              nonWgtBCDEF=puCorrHBCDEF.GetBinContent(1)
                              wgtBCDEF=puCorrHBCDEF.GetBinContent(2)
                          except: pass
                          puGH = opt.puNormSF.split("_")
                          if puGH[-1] in opt.run:
                              puGH[-1] = "GH"
                          else:
                              puGH.append("GH")
                          puGH = "_".join(puGH)
                          try:
                              puCorrHGH=fIn.Get(puGH)
                              nonWgtGH=puCorrHGH.GetBinContent(1)
                              wgtGH=puCorrHGH.GetBinContent(2)
                          except: pass
                          if wgtBCDEF>0 :
                              puNormSFBCDEF=nonWgtBCDEF/wgtBCDEF
                              if puNormSFBCDEF>1.3 or puNormSF<0.7 : 
                                  puNormSFBCDEF=1
                                  report += '%s wasn\'t be scaled as too large SF was found (probably low stats)\n' % sp[0]
                              else :
                                  report += '%s was scaled by %3.3f for pileup epoch %s normalization\n' % (sp[0],puNormSFBCDEF,"BCDEF")
                          if wgtGH>0 :
                              puNormSFGH=nonWgtGH/wgtGH
                              if puNormSF>1.3 or puNormSF<0.7 : 
                                  puNormSFGH=1
                                  report += '%s wasn\'t be scaled as too large SF was found (probably low stats)\n' % sp[0]
                              else :
                                  report += '%s was scaled by %3.3f for pileup epoch %s normalization\n' % (sp[0],puNormSFGH,"GH")

                        else:
                            puCorrH=fIn.Get(opt.puNormSF)
                            nonWgt=puCorrH.GetBinContent(1)
                            wgt=puCorrH.GetBinContent(2)
                            if wgt>0 :
                                puNormSF=nonWgt/wgt
                                if puNormSF>1.3 or puNormSF<0.7 : 
                                    puNormSF=1
                                    report += '%s wasn\'t be scaled as too large SF was found (probably low stats)\n' % sp[0]
                                else :
                                    report += '%s was scaled by %3.3f for pileup normalization\n' % (sp[0],puNormSF)

                    #fix r_B fit weighting normalization
                    rbFitSF=1
                    rbFitSFs=[[1,1],[1,1],[1,1]]
                    if opt.rbFit and not isData:
                        if(opt.run == "BCDEFGH"):
                          for num,rbFit in enumerate(["rbwgt_jpsi_BCDEFGH","rbwgt_d0_BCDEFGH","rbwgt_d0mu_BCDEFGH"]):
                              rbBCDEF = rbFit.split("_")
                              if rbBCDEF[-1] in opt.run:
                                  rbBCDEF[-1] = "BCDEF"
                              else:
                                  rbBCDEF.append("BCDEF")
                              rbBCDEF = "_".join(rbBCDEF)
                              try:
                                  rbCorrHBCDEF=fIn.Get(rbBCDEF)
                                  nonWgtBCDEF=rbCorrHBCDEF.GetBinContent(1)
                                  wgtBCDEF=rbCorrHBCDEF.GetBinContent(2)
                              except: pass
                              rbGH = rbFit.split("_")
                              if rbGH[-1] in opt.run:
                                  rbGH[-1] = "GH"
                              else:
                                  rbGH.append("GH")
                              rbGH = "_".join(rbGH)
                              try:
                                  rbCorrHGH=fIn.Get(rbGH)
                                  nonWgtGH=rbCorrHGH.GetBinContent(1)
                                  wgtGH=rbCorrHGH.GetBinContent(2)
                              except: pass
                              if wgtBCDEF>0 :
                                  rbFitSFBCDEF=nonWgtBCDEF/wgtBCDEF
                                  if rbFitSFBCDEF>1.3 or rbFitSFBCDEF<0.7 : 
                                      rbFitSFBCDEF=1
                                      report += '%s wasn\'t be scaled as too large SF was found (probably low stats)\n' % sp[0]
                                  else :
                                      report += '%s was scaled by %3.3f for rB fit %s epoch %s normalization\n' % (sp[0],rbFitSFBCDEF,rbFit.split("_")[1],"BCDEF")
                                  rbFitSFs[num][0]=wgtBCDEF
                              if wgtGH>0 :
                                  rbFitSFGH=nonWgtGH/wgtGH
                                  if rbFitSFGH>1.3 or rbFitSFGH<0.7 : 
                                      rbFitSFGH=1
                                      report += '%s wasn\'t be scaled as too large SF was found (probably low stats)\n' % sp[0]
                                  else :
                                      report += '%s was scaled by %3.3f for rB fit %s epoch %s normalization\n' % (sp[0],rbFitSFGH,rbFit.split("_")[1],"GH")
                                  rbFitSFs[num][1]=wgtGH
                    try:
                        for tkey in fIn.GetListOfKeys():

                            key=tkey.GetName()
                            keep=True
                            keep=False if len(onlyList)>0 else True
                            for pname in onlyList:
                                if pname in key: keep=True
                                #if "_ljet_" in pname and "_m_" in key: keep=True
                                #if "_ljet_" in pname and "_m_" in key: keep=True
                            #hack to ignore WJets in D meson mass plots FIXME
                            #if "massD" in key and "WJets" in tag:
                                #keep=False
                            #hack to ignoe WJets and DY in JPSi plots
                              #Single event with large weight
                            #if "_jpsi" in key and "WJets" in tag: continue
                            #if "_jpsi" in key and "DY" in tag: continue
                            #if "_tag" in key and "WJets" in tag: continue
                            #if "mass" in key and "_meson" in key and "WJets" in tag: continue
                            #if "_"+opt.run not in key and opt.run != "BCDEFGH" :
                            if opt.run not in key and opt.run != "BCDEFGH" :
                                 keep=False
                            #if key.split("_")[-1] != opt.run: continue
                            if "_BCDEFGH" in key: continue
                            if not keep: continue
                            obj=fIn.Get(key)
                            if not obj.InheritsFrom('TH1') : continue
                            if not isData and not '(data)' in sp[1]: 
                                sfVal=1.0
                                for procToScale in procSF:
                                    if sp[1]==procToScale:
                                    #if procToScale in sp[1]:
                                        for pcat in procSF[procToScale]:                                    
                                            if pcat not in key: continue
                                            sfVal=procSF[procToScale][pcat][0]
                                        #print 'Applying scale factor for ',sp[1],key,sfVal
                                lumi=obj.GetName().split("_")[-1]
                                if(lumi not in opt.run): lumi=opt.run
                                if(opt.run == "BCDEFGH" and lumi == "BCDEF"): puNormSF=puNormSFBCDEF
                                elif(opt.run == "BCDEFGH" and lumi == "GH"): puNormSF=puNormSFGH
                                if(opt.run == "BCDEFGH" and lumi == "BCDEF" and "TTJets" in tag): topPtNorm=topPtNormBCDEF
                                elif(opt.run == "BCDEFGH" and lumi == "GH" and "TTJets" in tag): topPtNorm=topPtNormGH
                                lumi=lumiList[lumi]
                                #if("TTJets" in tag): lumi=lumi*0.733417
                                #if "fsr" not in sp[0] and "isr" not in sp[0]:
                                    #xsec=1. #now stored in normH
                                
                                #if("_jpsi_" in obj.GetName()): 
                                    #if(opt.run == "BCDEFGH" and lumi == "BCDEF"): rbFitSF=rbFitSFs[0][0]
                                    #elif(opt.run == "BCDEFGH" and lumi == "GH"): rbFitSF=rbFitSFs[0][1]
                                #elif("_meson_" in obj.GetName() and "mu_tag" in obj.GetName()): 
                                    #if(opt.run == "BCDEFGH" and lumi == "BCDEF"): rbFitSF=rbFitSFs[2][0]
                                    #elif(opt.run == "BCDEFGH" and lumi == "GH"): rbFitSF=rbFitSFs[2][1]
                                #elif("_meson_" in obj.GetName()): 
                                    #if(opt.run == "BCDEFGH" and lumi == "BCDEF"): rbFitSF=rbFitSFs[1][0]
                                    #elif(opt.run == "BCDEFGH" and lumi == "GH"): rbFitSF=rbFitSFs[1][1]
                                if(opt.run == "BCDEFGH" and lumi == "BCDEF"): rbFitSF=rbFitSF
                                elif(opt.run == "BCDEFGH" and lumi == "GH"): rbFitSF=rbFitSF
                                obj.Scale(xsec*lumi*puNormSF*sfVal*topPtNorm*rbFitSF)
                                #obj.Scale(lumi*puNormSF*sfVal*topPtNorm)
                            if("JPsioJet" in obj.GetName()): obj.Rebin(2)
                            over=True
                            under=True
                            if "meson" in key: over=False
                            if "meson" in key: under=False
                            if "l3d" in key: over=False
                            if "l3d" in key: under=False
                            if "mass" in key: over=False
                            if "mass" in key: under=False
                            if "oJet" in key: over=False
                            fixExtremities(h=obj,addOverflow=over,addUnderflow=under)
                            if opt.rebin>1:  obj.Rebin(opt.rebin)
                            if opt.run != "BCDEFGH":
                                if not key in plots : plots[key]=Plot(key)
                                plots[key].add(h=obj,title=sp[1],color=sp[2],isData=sample[1],isSyst=isSyst)
                            else:
                                if key.split("_")[-1] in ("BCDEF","GH"):
                                    tmpkey = key.split("_")
                                    tmpkey[-1]="BCDEFGH"
                                    tmpkey="_".join(tmpkey)
                                    tmp = obj.GetName().split("_")
                                    tmp[-1]="BCDEFGH"
                                    tname="_".join(tmp)
                                    obj.SetName(tname)
                                    if not tmpkey in plots : plots[tmpkey]=Plot(tmpkey)
                                    plots[tmpkey].add(h=obj,title=sp[1],color=sp[2],isData=sample[1],isSyst=isSyst,rbList=rbSamplesList)
                    except:
                        pass

    #show plots
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    outDir=opt.inDir+'/plots'
    os.system('mkdir -p %s' % outDir)
    outName = opt.outName.replace(".root","_"+opt.run+".root")
    os.system('rm %s/%s'%(outDir,outName))
    for p in plots : 
        if opt.saveLog    : plots[p].savelog=True
        #if not opt.puNormSF    : plots[p].noPU=True
        lumiTotal=lumiList[opt.run]
        if not opt.silent : plots[p].show(outDir=outDir,lumi=lumiTotal,noStack=opt.noStack,saveTeX=opt.saveTeX)
        outName = opt.outName.replace(".root","_"+opt.run+".root")
        plots[p].appendTo('%s/%s'%(outDir,outName))
        plots[p].reset()

    print '-'*50
    print 'Plots and summary ROOT file can be found in %s' % outDir
    if len(report) : print report
    print '-'*50

        
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

