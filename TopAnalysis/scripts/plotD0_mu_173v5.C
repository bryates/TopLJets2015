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
#include "LJets2015/2016/mtop/tdr.C"

void plotD0_mu_173v5() {
std::vector<TString> tune = {"_700", "_ddown", "", "_ccentral", "_fit"};
std::vector<float> param = {0.700, 0.800, 0.855,0.900, 0.80};
/*
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up", "_fit" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055, 0.80};
*/
std::vector<int> color = {kBlue, kCyan-3, kOrange+3, kRed+1, kBlack}; //kOrange, kYellow-7, kViolet-4, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet, kBlack};

//RooHist *tmp;
std::vector<TH1F*> tmp;
TString cut("j_pt_ch<150");
int bin(22);
TCanvas *c1 = setupCanvas();
TPad *p1 = setupPad();
p1->cd();
float iniy=0.83;
float dy=0.03;
float ndy=param.size();
TLegend *leg = new TLegend(0.2, iniy-dy*ndy, 0.29, iniy+0.05);
leg->SetBorderSize(1);
leg->SetFillStyle(0);      
leg->SetTextFont(43);
leg->SetTextSize(12);
for(size_t i=0; i<tune.size(); i++){
TString name(TString::Format("LJets2015/2016/mtop/sPlot/sPlot/TopMass_m173v5%s_sPlot_d0_mu_tag_mu.root",tune[i].Data()));
//TString name(TString::Format("LJets2015/2016/mtop/TopMass_fsr-up_%s_sPlot_d0_mu_tag_mu.root",tune[i].Data()));
std::cout << name << std::endl;
TFile *fin = TFile::Open(name);
/*
RooWorkspace *w = (RooWorkspace*)fin->Get("w");
w->var("ptfrac")->setBins(bin);
RooRealVar ptfrac=*w->var("ptfrac");
ptfrac.setBins(22);
RooDataSet *sigData = (RooDataSet*)w->data("sigData");//->reduce(cut);
RooDataHist *ptfrac_hist = new RooDataHist("ptfrac_hist", "ptfrac_hist", *w->var("ptfrac"), *sigData);//*mcData);
RooPlot *ptfrac_mc = w->var("ptfrac")->frame();
ptfrac_hist->plotOn(ptfrac_mc);
*/

  //tmp.push_back((TH1F*)convert(ptfrac_mc,false,0,1.1));
  tmp.push_back(convert((RooPlot*)fin->Get("ptfrac_mu_tag_signal")));
  std::cout << "converted" << std::endl;
  tmp.back()->SetLineColor(color[i]);
  std::cout << "line color" << std::endl;
  //if(i==0) tmp->SetTitle("J/#Psi p_{T} / #Sigma p_{T}^{ch}");
  tmp.back()->SetTitle(TString::Format("%0.3f",param[i]));
  if(tune[i] == TString("_fit")) tmp.back()->SetLineStyle(10);
  std::cout << "Title" << std::endl;
  tmp.back()->Scale(1./tmp.back()->Integral());
  if(i==0) {
    tmp.front()->GetXaxis()->SetTitle("(D^{0} p_{T} + #mu p_{T}) / #Sigma p_{T}^{ch}");
    tmp.front()->GetYaxis()->SetTitle("Jets / 0.05");
    tmp.back()->Draw("hist");
    tdr(tmp.back());
    tmp.back()->SetTitle(TString::Format("%0.3f",param[0]));
    tmp.back()->GetXaxis()->SetRangeUser(0,1.1);
    //tmp.back()->GetYaxis()->SetRangeUser(0,240);
  }
  else tmp.back()->Draw("same hist");
  leg->AddEntry( tmp.back(), tmp.back()->GetTitle(),"lp");
}
TFile *fin = TFile::Open("LJets2015/2016/mtop/sPlot/sPlot/TopMass_Data_sPlot_d0_mu_tag_mu.root");
/*
RooWorkspace *w = (RooWorkspace*)fin->Get("w");
w->var("ptfrac")->setBins(bin);
RooRealVar ptfrac=*w->var("ptfrac");
RooDataSet *sigData = (RooDataSet*)w->data("sigData")->reduce(cut);
RooDataHist *ptfrac_data_hist = new RooDataHist("ptfrac_hist", "ptfrac_hist", *w->var("ptfrac"), *sigData);//*dataData);
RooPlot *ptfrac_data = w->var("ptfrac")->frame();
ptfrac.setBins(22);
ptfrac_data_hist->plotOn(ptfrac_data);
tmp.push_back((TH1F*)convert(ptfrac_data,false,0,1.1));
*/
tmp.push_back(convert((RooPlot*)fin->Get("ptfrac_mu_tag_signal")));
tmp.back()->Scale(1./tmp.back()->Integral());
//tmp.push_back(convert((RooPlot*)fin->Get("ptfrac_signal")));
//tmp.back()->Rebin(2);
tmp.back()->SetLineColor(1);
tmp.back()->SetMarkerStyle(20);
tmp.back()->SetTitle("Data");
leg->AddEntry( tmp.back(), tmp.back()->GetTitle(),"l");
tmp.back()->Draw("same");
gStyle->SetOptStat(0);
/*

tmp.front()->SetTitle("p_{T}(D^{0}_{#mu} + #mu) / #Sigma p_{T}^{ch}");
tmp.front()->GetXaxis()->SetTitle("(D^{0}_{#mu} p_{T} + #mu p_{T}) / #Sigma p_{T}^{ch}");
tmp.front()->GetYaxis()->SetTitle("Jets / 0.05");
tdr(tmp.front());
*/
leg->Draw();
tmp.front()->SetTitle("");
/*
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_d0_mu_tag_m173v5.png");
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_d0_mu_tag_m173v5.pdf");
*/
}
