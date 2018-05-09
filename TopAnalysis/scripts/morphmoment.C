#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TPaveText.h"
#include <iostream>
#include <vector>
#include <RooFit.h>
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooAddition.h"
#include "RooArgSet.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooMomentMorph.h"
#include "RooNLLVar.h"
#include "RooBinning.h"

using namespace RooFit;

std::vector<TString> samples = { "_d0_mu_tag_mu", "_d0", "_jpsi" };
//std::vector<TString> samples = { "_d0_mu_tag_mu", "_jpsi" };
//std::vector<TString> samples = { "_d0_mu_tag_mu" };
//std::vector<TString> samples = { "_d0" };
//TUNES=[ ('up','BL',1.055),  ('uup','BL',1.000), ('uuup','BL',0.975), ('central','BL',0.955), ('ccentral','BL',0.900), ('cccentral','BL',0.875), ('cuetp8m2t4','BL',0.855), ('ddown','BL',0.800), ('dddown','BL',0.775), ('down','BL',0.755)]
/*
std::vector<TString> tunes = {"_down", "_ddown", "_dddown", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.755, 0.775, 0.800, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055};
*/
std::vector<TString> tunes = {"_down", "", "_central", "_up" };
std::vector<float> param = {0.755, 0.855, 0.955, 1.055};
std::vector<RooDataHist*> ptfrac_mc_hist, ptfrac_data_hist;
std::vector<RooHistPdf*> ptfrac_mc_pdf, ptfrac_data_pdf;
std::vector<RooHistPdf*> histMC, histData;
RooRealVar ptfrac;
RooPlot *frame;
int bins(22);

void morphmoment() {
/*
RooBinning bins(0,1.01);
bins.addBoundary(0.);
bins.addBoundary(0.2);
bins.addBoundary(0.4);
bins.addBoundary(0.5);
bins.addBoundary(0.6);
bins.addBoundary(0.7);
bins.addBoundary(0.8);
bins.addBoundary(0.9);
bins.addBoundary(1.0);
*/
  std::vector<TFile*> filesdata;
  std::vector<RooWorkspace*> wdata;
  for(auto it : samples) {
    TString fname = TString::Format("TopMass_Data_sPlot%s.root",it.Data());
    std::cout << "Opening " << fname << std::endl;
    filesdata.push_back(TFile::Open(fname));
  }
  //load Data into RooDataHist
  for(auto &file : filesdata) {
    wdata.push_back((RooWorkspace*)file->Get("w"));
    wdata.back()->var("ptfrac")->setBins(bins);
  }
  ptfrac=*(RooRealVar*)(wdata[0]->var("ptfrac")->Clone());
  frame = ptfrac.frame();
  for(auto &w : wdata) {
    RooDataSet *sigData = (RooDataSet*)w->data("sigData");
    int pos = &w - &wdata[0];
    std::cout << pos << std::endl;
    //if(samples[pos] == "_d0")
    sigData = (RooDataSet*)sigData;
    //sigData = (RooDataSet*)sigData->reduce(cut);
    TString title(TString::Format("ptfrac_hist_Data%s",samples[pos].Data()));
    ptfrac_data_hist.push_back(new RooDataHist(title, title, ptfrac, *sigData));
    //ptfrac_data_hist.push_back(new RooDataHist(title, title, ptfrac, *w->data("sigData")));
    title.ReplaceAll("hist","pdf");
    std::cout << title << " " << ptfrac_data_hist.back()->sumEntries() << std::endl;
    //ptfrac_data_pdf.push_back(new RooHistPdf(title, title, ptfrac, *ptfrac_data_hist.back()));//*mcData);
    //histData.push_back(new RooHistPdf(title, title, ptfrac, *ptfrac_data_hist.back()));
  }
  for(auto tune : tunes) {
    std::vector<TFile*> filesmc;//, filesdata;
    std::vector<RooWorkspace*> wmc;//, wdata;
    std::vector<double> tuneWgts;
    for(auto it : samples) {
      TString fname = TString::Format("TopMass_172v5%s_sPlot%s.root",tune.Data(),it.Data());
      std::cout << "Opening " << fname << std::endl;
      filesmc.push_back(TFile::Open(fname));
      //store N_evt weights for alternate tunes to scale N_evt (shape only)
      TH1F *tuneWgt = (TH1F*)filesmc.back()->Get("tuneWgt");
      //if (tuneWgt->GetBinContent(2))
      if (it != "")
        tuneWgts.push_back(tuneWgt->GetBinContent(1)/tuneWgt->GetBinContent(2));
      //nominal sample, no weight needed
      else tuneWgts.push_back(1.);
      //loading data is not optimal
      /*
      fname = TString::Format("TopMass_Data_sPlot%s.root",it.Data());
      std::cout << "Opening " << fname << std::endl;
      filesdata.push_back(TFile::Open(fname));
      */
    }
    
    //load RooWorkspace and set binning
    for(auto &file : filesmc) {
      wmc.push_back((RooWorkspace*)file->Get("w"));
      wmc.back()->var("ptfrac")->setBins(bins);
    }
    if(tune == tunes[0]) {
      ptfrac=*(RooRealVar*)(wmc[0]->var("ptfrac")->Clone());
      frame = ptfrac.frame();
    }

    //loading data is not optimal
    /*
    for(auto &file : filesdata) {
      wdata.push_back((RooWorkspace*)file->Get("w"));
      wdata.back()->var("ptfrac")->setBins(bins);
    }
    */
    
    //load MC into RooDataHist
    std::cout << "building DataHists" << std::endl;
    std::cout << "building HistPdfs" << std::endl;
    int num(0);
    TString cut = "j_pt_ch>50 && j_pt_ch<100";
    for(auto &w : wmc) {
      RooDataSet *sigData = (RooDataSet*)w->data("sigData");
      int pos = &w - &wmc[0];
      //if(samples[pos] == "_d0")
      sigData = (RooDataSet*)sigData;
      //sigData = (RooDataSet*)sigData->reduce(cut);
      TString title(TString::Format("ptfrac_hist_MC%s%s",tune.Data(),samples[pos].Data()));
      //tuneWgts[pos] to scale hist back to N_no_weight (shape only)
      ptfrac_mc_hist.push_back(new RooDataHist(title, title, ptfrac, *sigData));//, tuneWgts[pos]));
      std::cout << tuneWgts[pos] << std::endl;
      //ptfrac_mc_hist.push_back(new RooDataHist(title, title, ptfrac, *w->data("sigData")));//*mcData);
      title.ReplaceAll("hist","pdf");
      std::cout << title << " " << ptfrac_mc_hist.back()->sumEntries() << std::endl;
      ptfrac_mc_pdf.push_back(new RooHistPdf(title, title, ptfrac, *ptfrac_mc_hist.back()));//*mcData);
      histMC.push_back(new RooHistPdf(title, title, ptfrac, *ptfrac_mc_hist.back()));
    }

    //load Data into RooDataHist
    //loading data is not optimal
    /*
    for(auto &w : wdata) {
      RooDataSet *sigData = (RooDataSet*)w->data("sigData");
      int pos = &w - &wdata[0];
      //if(samples[pos] == "_d0")
      sigData = (RooDataSet*)sigData;
      //sigData = (RooDataSet*)sigData->reduce(cut);
      TString title(TString::Format("ptfrac_hist_Data%s",samples[pos].Data()));
      ptfrac_data_hist.push_back(new RooDataHist(title, title, ptfrac, *sigData));
      //ptfrac_data_hist.push_back(new RooDataHist(title, title, ptfrac, *w->data("sigData")));
      title.ReplaceAll("hist","pdf");
      std::cout << title << " " << ptfrac_data_hist.back()->sumEntries() << std::endl;
      ptfrac_data_pdf.push_back(new RooHistPdf(title, title, ptfrac, *ptfrac_data_hist.back()));//x*mcData);
      histData.push_back(new RooHistPdf(title, title, ptfrac, *ptfrac_data_hist.back()));
    }
    */
  }


  //RooRealVar rB = RooRealVar("rB", "r_{B}", 0.855, 0.845, 0.865);
  RooRealVar rB = RooRealVar("rB", "r_{B}", 0.855, 0.755, 1.055);
  std::vector<RooArgList> pdfs;
  for(size_t i = 0; i < samples.size(); i++)
    pdfs.push_back(RooArgList());
  RooArgList varlist;
  //pds = new RooArgList[samples.size()]; //one for each sample
  TVectorD paramVec = TVectorD(tunes.size());
  for(size_t i = 0; i < tunes.size(); i++) {
    for(size_t j = 0; j < samples.size(); j++) {
      int pos = j + i * samples.size();
      std::cout << pos << std::endl;
      std::cout << histMC.size() << std::endl;
      histMC[pos]->Print();
      histMC[pos]->plotOn(frame, RooFit::Binning(bins));
      pdfs[j].add(*histMC[pos]);
      if(j == 0)
        paramVec[i] = param[i];
    }
  }
  pdfs[0].Print();
  varlist.add(ptfrac);
  paramVec.Print();
  //frame->Draw();

  std::vector<RooMomentMorph> morph;
  //std::vector<std::pair<float,float>> rBval;
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  TString report;
  for(auto &it : ptfrac_data_hist)
    std::cout << it->GetName() << std::endl;
  for(size_t i = 0; i < samples.size(); i++) {
    morph.push_back(RooMomentMorph(Form("morph%s",samples[i].Data()),Form("morph%s",samples[i].Data()), rB, varlist, pdfs[i], paramVec,RooMomentMorph::Linear));
    TH1* hh = morph[i].createHistogram("hh",ptfrac,Binning(22),RooFit::YVar(rB,Binning(30)));
    //TH1* hh = morph[i].createHistogram("hh",ptfrac,Binning(22),RooFit::YVar(rB,Binning(30)));
    TString name(Form("%s vs p_{T} and r_{B}",samples[i].Data()));
    name.ReplaceAll("_d0_mu_tag_mu","D^{0}_{#mu}+#mu_{tag}");
    name.ReplaceAll("_d0","D^{0}");
    name.ReplaceAll("_jpsi","J/#Psi");
    hh->SetTitle(name);
    hh->GetYaxis()->SetNoExponent(kTRUE);
    hh->GetXaxis()->SetTitleOffset(2.0);
    hh->GetYaxis()->SetTitleOffset(2.0);
    hh->GetZaxis()->SetTitleOffset(1.5);
    hh->Draw("lego");
    c1->SaveAs(Form("TemplateMorph_172v5%s_rB_%d-%d.pdf",samples[i].Data(),(int)(param.front()*1000),(int)(param.back()*1000)));
    c1->SaveAs(Form("TemplateMorph_172v5%s_rB_%d-%d.png",samples[i].Data(),(int)(param.front()*1000),(int)(param.back()*1000)));

    //draw overlays
    //rB = RooRealVar("rB", "r_{B}", 0.855, 0.755, 1.055);
    /*
    frame = ptfrac.frame();
    std::vector<int> colors = {kBlue, kBlue+3, kGreen, kGreen+3, kOrange, kOrange+3, kCyan, kCyan+3, kMagenta, kMagenta+3};
    int loops(colors.size());
    for(int j = 0; j < loops; j++) {
      float step = (rB.getMax() - rB.getMin())/loops;
      rB.setVal(rB.getMin() + step);
      rB.Print();
      morph[i].plotOn(frame, RooFit::Binning(bins), RooFit::LineColor(colors[j]));
    }
    frame->Draw();
    c1->SaveAs(TString::Format("MorphOverlay_rB%s.pdf", samples[i].Data()));
    c1->SaveAs(TString::Format("MorphOverlay_rB%s.png", samples[i].Data()));
    */

    //fit morph to data
    //rB = RooRealVar("rB", "r_{B}", 0.855, 0.755, 1.055);
    RooFitResult *rBFit = morph[i].fitTo(*ptfrac_data_hist[i]);//,"rse");
    //std::cout << "Best fit rB = " << rBFit->getVal() << std::endl;
    frame = ptfrac.frame();
    ptfrac_data_hist[i]->plotOn(frame, RooFit::Binning(bins));
    report += TString::Format("%s r_B = %.3f +/- %.2e\n", samples[i].Data(), rB.getVal(), rB.getError());
    //rB->setVal(rBFit->getVal());
    morph[i].plotOn(frame, RooFit::Binning(bins), FillColor(38));
    float up = rB.getVal()+rB.getError();
    float down = rB.getVal()-rB.getError();
    /*
    rB.setVal(up);
    morph[i].plotOn(frame, FillColor(38), LineStyle(kDashed));
    rB.setVal(down);
    morph[i].plotOn(frame, FillColor(38), LineStyle(kDashed));
    */
    TPaveText *pt = new TPaveText(0.12,0.85,0.3,0.65,"NDC"); //NB blNDC
    pt->SetFillStyle(0);
    pt->SetTextAlign(11);
    pt->SetBorderSize(0);
    pt->SetTextFont(42);
    pt->SetTextSize(0.046);
    TString text = TString::Format("r_{B}= %.3f #pm %.2e",rB.getVal(), rB.getError());
    pt->AddText(text);
    frame->addObject(pt);
    frame->Draw();
    c1->SaveAs(TString::Format("BestFit_rB%s.pdf", samples[i].Data()));
    c1->SaveAs(TString::Format("BestFit_rB%s.png", samples[i].Data()));
    frame = rB.frame(Bins(100), Range(rB.getMin(), rB.getMax()-0.1));
    //frame = rB.frame(Bins(100), Range(0.755,0.955));
    RooNLLVar nll("nll", "nll", morph[i], *ptfrac_data_hist[i]);
    nll.plotOn(frame, ShiftToZero());
    frame->GetYaxis()->SetRangeUser(-1,20);
    if(samples[i] == "_jpsi")
    frame->GetYaxis()->SetRangeUser(-1,5);
    name = samples[i];
    name.ReplaceAll("_"," ");
    name.ReplaceAll(" mu tag mu","_{#mu} + #mu_{tag} ");
    name.ReplaceAll("d0","D^{0}");
    name.ReplaceAll("jpsi","J/#Psi");
    frame->SetTitle(Form("Log Likelihood Scan of r_{B}%s", name.Data()));
    //TLine *line = new TLine(rB.getMin(),0,rB.getMax(),0);
    TLine *line = new TLine(0.755,0,0.955,0);
    line->SetLineColor(kRed);
    frame->addObject(pt);
    frame->Draw();
    line->Draw();
    c1->SaveAs(TString::Format("NLL-Scan_rB%s.pdf", samples[i].Data()));
    c1->SaveAs(TString::Format("NLL-Scan_rB%s.png", samples[i].Data()));
  }
  std::cout << report << std::endl;

  //rB = RooRealVar("rB", "r_{B}", 0.855, 0.755, 1.055);
  RooCategory sample("sample","sample");
  /*
  for(int i = 0; i < samples.size(); i++) {
    sample.defineType(samples[i]);
  }
  */
  sample.defineType("d0_mu_tag");
  sample.defineType("d0");
  sample.defineType("jpsi");
  RooDataHist combined("combined", "combined", RooArgSet(ptfrac), Index(sample), Import("d0_mu_tag", *ptfrac_data_hist[0]), Import("d0", *ptfrac_data_hist[1]), Import("jpsi", *ptfrac_data_hist[2]));

  rB.setConstant(false);
  RooSimultaneous simPdf("simPdf", "simPdf", sample);
  /*
  for(int i = 0; i < samples.size(); i++) {
    simPdf.addPdf(morph[i], samples[i]);
  }
  */
  simPdf.addPdf(morph[0], "d0_mu_tag");
  simPdf.addPdf(morph[1], "d0");
  simPdf.addPdf(morph[2], "jpsi");
  simPdf.fitTo(combined, Minimizer("Minuit2","Migrad"), "e");
  simPdf.Print("v");

  frame = rB.frame(Bins(100), Range(0.755,0.955));
  /*
  combined.plotOn(frame, Cut("sample==sample::d0_mu_tag"));
  simPdf.plotOn(frame, Slice(sample, "d0_mu_tag"), ProjWData(sample, combined));
  */
  
  RooNLLVar nll("nll", "nll", simPdf, combined);
  nll.plotOn(frame, ShiftToZero());
  TLine *line = new TLine(0.755,0,0.955,0);
  line->SetLineColor(kRed);
  frame->GetYaxis()->SetRangeUser(-1,20);//frame->GetYaxis()->GetXmax());
  TPaveText *pt = new TPaveText(0.12,0.85,0.3,0.65,"NDC"); //NB blNDC
  pt->SetFillStyle(0);
  pt->SetTextAlign(11);
  pt->SetBorderSize(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.046);
  TString text = TString::Format("r_{B}= %.3f #pm %.2e",rB.getVal(), rB.getError());
  pt->AddText(text);
  frame->addObject(pt);
  frame->SetTitle("Log Likelihood Scan of r_{B} Simultaneous");
  frame->Draw();
  line->Draw();
  c1->SaveAs("NLL-Scan_rB_Simultaneous.pdf");
  c1->SaveAs("NLL-Scan_rB_Simultaneous.png");
  std::cout << rB.getVal() << " +/- " << rB.getError() << std::endl;
  //mg.Draw();
  //frame->Draw();
  /*
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  rB.setVal(0.855);
  morph.plotOn(frame, RooFit::Binning(22));
  rB.setVal(0.925);
  morph.plotOn(frame, LineColor(kBlue+2), RooFit::Binning(22));
  frame->Draw();
  */

}
