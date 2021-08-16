#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TPad.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooRealVar.h>
#include "LJets2015/2016/mtop/convert.h"
#include "LJets2015/2016/mtop/tdr.C"

void plotD0_mu(TString sample="172v5", int epoch=0, bool fin=false, float best=0.) {
std::vector<TString> tune = {"_700", "_ddown", "", "_ccentral", "_fit"};
std::vector<float> param = {0.700, 0.800, 0.855,0.900, 0.863};
std::vector<int> color = {kBlue, kCyan-3, kOrange+3, kRed+1, kBlack}; //kOrange, kYellow-7, kViolet-4, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet, kBlack};
/*
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up", "_fit" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055, 0};
std::vector<int> color = {kBlue, kCyan-3, kOrange+3, kRed+1, kOrange, kYellow-7, kViolet-4, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet, kBlack};
*/

//RooHist *tmp;
std::vector<TH1F*> tmp;
TString cut("j_pt_ch<150");
//TString sample("172v5");
std::vector<float> bin = {0, 0.2, 0.4, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0}; 
TCanvas *c1 = setupCanvas();
TPad *p1 = setupPad();
p1->cd();
//int epoch=1;
float iniy=0.83;
float dy=0.03;
float ndy=param.size();
TLegend *leg = new TLegend(0.15, iniy-dy*ndy, 0.24, iniy+0.05);
leg->SetBorderSize(1);
leg->SetFillStyle(0);      
leg->SetTextFont(43);
leg->SetTextSize(12);
for(size_t i=0; i<tune.size(); i++){
TString name("");
if(best>0. && param[i] != best) continue;
name = TString::Format("LJets2015/2016/mtop/sPlot/sPlot/TopMass_%s%s_sPlot_d0_mu_tag_mu.root",sample.Data(),tune[i].Data());
if(name.Contains("fit")) name.ReplaceAll(TString::Format("_%s%s_", sample.Data(), tune[i].Data()), "_172v5_");
name.ReplaceAll("d0_mu_tag_mu.root","d0_mu_tag_mu1.root");
TFile *f = TFile::Open(name);
std::cout << name << std::endl;
RooPlot *ptfrac_mc = (RooPlot*)f->Get("ptfrac_mu_tag_signal");
std::cout << name << std::endl;
TH1F *h1 = (TH1F*)convert(ptfrac_mc,false,bin);
h1->SetDirectory(0);
std::cout << name << std::endl;
f->Close();
name.ReplaceAll("d0_mu_tag_mu1.root","d0_mu_tag_mu2.root");
f = TFile::Open(name);
ptfrac_mc = (RooPlot*)f->Get("ptfrac_mu_tag_signal");
std::cout << name << std::endl;
TH1F *h2 = (TH1F*)convert(ptfrac_mc,false,bin);
h1->SetDirectory(0);
h1->Add(h2);
h1->Scale(1./h1->Integral());
f->Close();
//ptfrac_hist->plotOn(ptfrac_mc);

  //tmp = ((RooPlot*)f->Get("ptfrac_mu_tag_signal"))->getHist();
  //tmp.push_back((TH1F*)convert(ptfrac_mc,true,0,1.1));
  tmp.push_back(h1);
  //tmp.push_back(convert((RooPlot*)f->Get("ptfrac_mu_tag_signal")));
  tmp.back()->SetLineColor(color[i]);
  if(best>0.) tmp.back()->SetLineColor(kRed);
  //if(i==0) tmp->SetTitle("J/#Psi p_{T} / #Sigma p_{T}^{ch}");
  tmp.back()->SetTitle(TString::Format("%0.3f",param[i]));
  if(param[i]==0)
    tmp.back()->SetTitle("Fit");
  if(tune[i] == TString("_fit")) tmp.back()->SetLineStyle(10);
  tmp.back()->Scale(1./tmp.back()->Integral());
  if(i==0) {
    tmp.front()->GetXaxis()->SetTitle("(D^{0} p_{T} + #mu p_{T}) / #Sigma p_{T}^{ch}");
    tmp.front()->GetYaxis()->SetTitle("Jets / 0.05");
    tmp.back()->Draw("hist");
    tdr(tmp.back(),epoch,fin);
    tmp.back()->SetTitle(TString::Format("%0.3f",param[0]));
    tmp.back()->GetXaxis()->SetRangeUser(0,1.1);
    tmp.back()->GetYaxis()->SetRangeUser(0,0.45);
    /*
    if(jpT) tmp.back()->GetYaxis()->SetRangeUser(0,0.16);
    else tmp.back()->GetYaxis()->SetRangeUser(0,0.24);
    */
  }
  else tmp.back()->Draw("same hist");
  tmp.back()->Draw("same e");
  leg->AddEntry( tmp.back(), tmp.back()->GetTitle(),"lp");
}
TString name(TString::Format("LJets2015/2016/mtop/sPlot/sPlot/TopMass_Data_sPlot_d0_mu_tag_mu1.root"));
f = TFile::Open(name);
std::cout << name << std::endl;
/*
RooWorkspace *w = (RooWorkspace*)f->Get("w");
w->var("ptfrac")->setBins(bin);
RooRealVar ptfrac=*w->var("ptfrac");
ptfrac.setBins(22);
RooDataSet *sigData = (RooDataSet*)w->data("sigData")->reduce(cut);
RooDataHist *ptfrac_data_hist = new RooDataHist("ptfrac_hist", "ptfrac_hist", *w->var("ptfrac"), *sigData);//*dataData);
RooPlot *ptfrac_data = w->var("ptfrac")->frame();
*/
RooPlot *ptfrac_data = (RooPlot*)f->Get("ptfrac_mu_tag_signal");
TH1F *d1 = (TH1F*)convert(ptfrac_data,false,bin);
d1->SetDirectory(0);
f->Close();
name.ReplaceAll("d0_mu_tag_mu1.root","d0_mu_tag_mu2.root");
f = TFile::Open(name);
ptfrac_data = (RooPlot*)f->Get("ptfrac_mu_tag_signal");
std::cout << name << std::endl;
TH1F *d2 = (TH1F*)convert(ptfrac_data,false,bin);
d2->SetDirectory(0);
d1->Add(d2);
d1->Scale(1./d1->Integral());
f->Close();
//ptfrac_data_hist->plotOn(ptfrac_data);
//tmp.push_back((TH1F*)convert(ptfrac_data,true,0,1.1));
tmp.push_back(d1);
//tmp.push_back(convert((RooPlot*)f->Get("ptfrac_mu_tag_signal")));
tmp.back()->Scale(1./tmp.back()->Integral());
tmp.back()->SetMarkerStyle(20);
tmp.back()->SetLineColor(1);
tmp.back()->SetTitle("Data");
leg->AddEntry( tmp.back(), tmp.back()->GetTitle(),"l");
tmp.back()->Draw("same");
gStyle->SetOptStat(0);

for(auto &plot : tmp) {
  if(plot == tmp.back()) continue;
  std::cout << plot->GetTitle() << ": ";
  std::cout << tmp.back()->Chi2Test(plot, "CHI2 WW") << std::endl;
}
/*

tmp.front()->GetXaxis()->SetTitle("J/#Psi p_{T} / #Sigma p_{T}^{ch}");
tmp.front()->GetYaxis()->SetTitle("Jets / 0.05");
tdr(tmp.front());
*/
leg->Draw();
tmp.front()->SetTitle("");
TString b("");
if(best>0.) b = TString::Format("_%d", (int)(best*1000));
TString ext(TString::Format("_%s",sample.Data()));
if(fin) {
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_d0_mu_tag_mu" + ext + b + "_d0_mu_tag_mu_final.png");
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_d0_mu_tag_mu" + ext + b + "_d0_mu_tag_mu_final.pdf");
}
else {
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_d0mu" + ext + b + "_d0mu.png");
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_d0mu" + ext + b + "_d0mu.pdf");
}

}
