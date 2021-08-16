#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TPad.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooRealVar.h>
#include "convert.h"
#include "tdr.C"

void plotJPsi_fsr() {
std::vector<TString> tune = {"", "fsr-up", "fsr-down", "_fit"};
std::vector<float> param = {0.855, 0.855, 0.855, 0.80};
/*
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up", "_fit" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055, 0.80};
*/
std::vector<int> color = {kBlue, kRed, kOrange+3, kBlack}; //kOrange, kYellow-7, kViolet-4, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet, kBlack};

//RooHist *tmp;
std::vector<TH1F*> tmp;
TString cut("j_pt_ch<150");
int bin(22);
TCanvas *c1 = new TCanvas("c1","c1");
c1->cd();
for(size_t i=0; i<tune.size(); i++){
TString name(TString::Format("LJets2015/2016/mtop/TopMass_172v5%s_sPlot_jpsi.root",tune[i].Data()));
if(tune[i].Contains("fsr"))
  name=TString::Format("LJets2015/2016/mtop/TopMass_%s_sPlot_jpsi.root",tune[i].Data());
std::cout << name << std::endl;
TFile *fin = TFile::Open(name);
RooWorkspace *w = (RooWorkspace*)fin->Get("w");
w->var("ptfrac")->setBins(bin);
RooRealVar ptfrac=*w->var("ptfrac");
ptfrac.setBins(22);
RooDataSet *sigData = (RooDataSet*)w->data("sigData")->reduce(cut);
RooDataHist *ptfrac_hist = new RooDataHist("ptfrac_hist", "ptfrac_hist", *w->var("ptfrac"), *sigData);//*mcData);
RooPlot *ptfrac_mc = w->var("ptfrac")->frame();
ptfrac_hist->plotOn(ptfrac_mc);

  //tmp = ((RooPlot*)fin->Get("ptfrac_signal"))->getHist();
  tmp.push_back((TH1F*)convert(ptfrac_mc,false,0,1.1));
  //tmp.push_back(convert((RooPlot*)fin->Get("ptfrac_signal")));
  std::cout << "converted" << std::endl;
  tmp.back()->SetLineColor(color[i]);
  std::cout << "line color" << std::endl;
  //if(i==0) tmp->SetTitle("J/#Psi p_{T} / #Sigma p_{T}^{ch}");
  tmp.back()->SetTitle(TString::Format("%s",tune[i].Data()));
  if(tune[i] == "") tmp.back()->SetTitle("Nominal");
  if(tune[i] == "_fit") tmp.back()->SetTitle("Fit");
  if(tune[i] == TString("_fit")) tmp.back()->SetLineStyle(10);
  std::cout << "Title" << std::endl;
  tmp.back()->Scale(1./tmp.back()->Integral());
  if(i==0) {
    tmp.back()->Draw("hist");
    tmp.back()->GetXaxis()->SetRangeUser(0,1.1);
    //tmp.back()->GetYaxis()->SetRangeUser(0,.17);
  }
  else tmp.back()->Draw("same hist");
}
TFile *fin = TFile::Open("LJets2015/2016/mtop/TopMass_Data_sPlot_jpsi.root");
RooWorkspace *w = (RooWorkspace*)fin->Get("w");
w->var("ptfrac")->setBins(bin);
RooRealVar ptfrac=*w->var("ptfrac");
RooDataSet *sigData = (RooDataSet*)w->data("sigData")->reduce(cut);
RooDataHist *ptfrac_data_hist = new RooDataHist("ptfrac_hist", "ptfrac_hist", *w->var("ptfrac"), *sigData);//*dataData);
RooPlot *ptfrac_data = w->var("ptfrac")->frame();
ptfrac.setBins(22);
ptfrac_data_hist->plotOn(ptfrac_data);
tmp.push_back((TH1F*)convert(ptfrac_data,false,0,1.1));
tmp.back()->Scale(1./tmp.back()->Integral());
//tmp.push_back(convert((RooPlot*)fin->Get("ptfrac_signal")));
//tmp.back()->Rebin(2);
tmp.back()->SetLineColor(1);
tmp.back()->SetMarkerStyle(20);
tmp.back()->SetTitle("Data");
tmp.back()->Draw("same");
gStyle->SetOptStat(0);

TLegend *leg = c1->BuildLegend();
tmp.front()->GetXaxis()->SetTitle("J/#Psi p_{T} / #Sigma p_{T}^{ch}");
tmp.front()->GetYaxis()->SetTitle("Jets / 0.05");
tdr(tmp.front());
}
