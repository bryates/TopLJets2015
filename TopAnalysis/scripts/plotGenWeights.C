#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TPad.h>
#include <RooPlot.h>
#include <RooHist.h>
#include "convert.h"
#include "/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/mtop/tdr_sim.C"

void plotGenWeights(TString proc="d0", bool fin=false) {
gStyle->SetOptStat(0);
std::vector<TString> tune = {TString::Format("sdowngen_%sFrag", proc.Data()), TString::Format("700gen_%sFrag", proc.Data()), TString::Format("725gen_%sFrag", proc.Data()), TString::Format("downgen_%sFrag", proc.Data()), TString::Format("dddowngen_%sFrag", proc.Data()), TString::Format("ddowngen_%sFrag", proc.Data()), TString::Format("scentralgen_%sFrag", proc.Data()), TString::Format("cccentralgen_%sFrag", proc.Data()), TString::Format("ccentralgen_%sFrag", proc.Data()), TString::Format("925gen_%sFrag", proc.Data()), TString::Format("centralgen_%sFrag", proc.Data()), TString::Format("uuupgen_%sFrag", proc.Data()), TString::Format("uupgen_%sFrag", proc.Data()), TString::Format("upgen_%sFrag", proc.Data()), TString::Format("Dzuchargedincgen_%sFrag", proc.Data()), TString::Format("Dzdchargedincgen_%sFrag", proc.Data()), TString::Format("Dzbuchargedincgen_%sFrag", proc.Data()), TString::Format("Dzbdchargedincgen_%sFrag", proc.Data())};
std::vector<TString> param = {"0.655", "0.700", "0.725", "0.755", "0.775", "0.800", "0.825", "0.875", "0.900", "0.925", "0.955", "0.975", "1.000", "1.055", "D^{0} up", "D^{0} down", "D^{0} up", "D^{0} down"} ;
tune = {TString::Format("sdowngen_%sFrag", proc.Data()), TString::Format("725gen_%sFrag", proc.Data()), TString::Format("downgen_%sFrag", proc.Data()), TString::Format("dddowngen_%sFrag", proc.Data()), TString::Format("scentralgen_%sFrag", proc.Data()), "",TString::Format("cccentralgen_%sFrag", proc.Data()), TString::Format("ccentralgen_%sFrag", proc.Data()), TString::Format("centralgen_%sFrag", proc.Data()), TString::Format("uuupgen_%sFrag", proc.Data()),TString::Format("B0Dzugen_%sFrag", proc.Data()), TString::Format("B0Dzdgen_%sFrag", proc.Data()), TString::Format("B0Dzbugen_%sFrag", proc.Data()), TString::Format("B0Dzbdgen_%sFrag", proc.Data()), TString::Format("BpmDzugen_%sFrag", proc.Data()), TString::Format("BpmDzdgen_%sFrag", proc.Data()), TString::Format("BpmDzbugen_%sFrag", proc.Data()), TString::Format("BpmDzbdgen_%sFrag", proc.Data()), TString::Format("BsDzugen_%sFrag", proc.Data()), TString::Format("BsDzdgen_%sFrag", proc.Data()), TString::Format("BsDzbugen_%sFrag", proc.Data()), TString::Format("BsDzbdgen_%sFrag", proc.Data()), TString::Format("BLbDzugen_%sFrag", proc.Data()), TString::Format("BLbDzdgen_%sFrag", proc.Data()), TString::Format("BLbDzbugen_%sFrag", proc.Data()), TString::Format("BLbDzbdgen_%sFrag", proc.Data())};
std::cout << tune.size() << std::endl;
param = {"0.655", "0.725", "0.755", "0.775", "0.825", "0.855", "0.875", "0.9", "0.955", "0.975", "B^{0} #rightarrow D^{0} up", "B^{0} #rightarrow D^{0} down", "B^{0} #rightarrow #bar{D^{0}} up", "B^{0} #rightarrow #bar{D^{0}} down", "B^{#pm} #rightarrow D^{0} up", "B^{#pm} #rightarrow D^{0} down", "B^{#pm} #rightarrow #bar{D^{0}} up", "B^{#pm} #rightarrow #bar{}D^{0}} down", "B* #rightarrow D^{0} up", "B* #rightarrow D^{0} down", "B* #rightarrow #bar{D^{0}} up", "B* #rightarrow #bar{D^{0}} down", "#Lambda_{b} #rightarrow D^{0} up", "#Lambda_{b} #rightarrow D^{0} down", "#Lambda_{b} #rightarrow #bar{D^{0}} up", "#Lambda_{b} #rightarrow #bar{D^{0}} down"};
//std::vector<TString> tune = {TString::Format("sdowngen_%sFrag", proc.Data()), TString::Format("700gen_%sFrag", proc.Data()), TString::Format("725gen_%sFrag", proc.Data()), TString::Format("scentralgen_%sFrag", proc.Data()), TString::Format("ccentralgen_%sFrag", proc.Data()), TString::Format("925gen_%sFrag", proc.Data()), TString::Format("centralgen_%sFrag", proc.Data())};//, TString::Format("uuupgen_%sFrag", proc.Data())};
//std::vector<TString> param = {"0.655", "0.700", "0.725", "0.825", "0.900", "0.925", "0.955"};//, "0.975"};
std::vector<int> color = {kBlue, kCyan-3, kOrange+3, kRed+1, kOrange, kYellow-7, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet+3, kRed-5, kBlack, kRed, kOrange, kBlue+5, kGreen+2, kRed+4, kViolet+2, kCyan+2, kMagenta-3};
//std::vector<int> color = {kBlue, kCyan-3, kOrange+3, kRed+1, kOrange, kYellow-7, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet+3, kRed-5, kBlack, kRed, kRed, kCyan, kCyan};
//std::vector<int> color = {kBlue, kCyan-3, kOrange+3, kRed+1, kOrange, kYellow-7, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet+3, kRed-5, kBlack, kBlue, kGreen, kRed, kRed, kCyan, kCyan};

//TFile *f = TFile::Open("/eos/cms/store/user/byates/top18012/bfragweights.root");
TFile *f = TFile::Open("data/era2016/bfragweights.root");
std::vector<TH1F*> tmp;
//TCanvas *c1 = new TCanvas("c1","c1",500,450);
//c1->cd();
TCanvas *c1 = setupCanvas();
c1->SetWindowSize(550,500);
TPad *p1 = new TPad("p1","p1",0.0,0.0,1.0,1.0);
p1->Range(-0.1536387,0.5421622,1.141653,1.438919);
//p1->SetRightMargin(0.05);
p1->SetLeftMargin(0.12);
p1->SetTopMargin(0.1);
p1->SetBottomMargin(0.12);
p1->Draw();
p1->cd();
//TPad *p1 = setupPad();
//p1->cd();
gPad->cd();
float iniy=0.4;
float dy=0.03;
float ndy=param.size()/3;
TLegend *leg = new TLegend(0.21, iniy-dy*ndy, 0.79, iniy+0.05);
//TLegend *leg = new TLegend(0.91, iniy-dy*ndy, 0.99, iniy+0.05);
leg->SetBorderSize(0);
leg->SetFillStyle(1001);
leg->SetTextFont(43);
leg->SetTextSize(12);
//leg->SetHeader("#it{r}_{b} value", "C");
leg->SetNColumns(3);
//leg->SetHeader("#it{r}_{b} value", "C");
for(size_t i=0; i<tune.size(); i++){

  //tmp = ((RooPlot*)f->Get("ptfrac_signal"))->getHist();
  if(tune[i] == TString("fitgenFrag")) continue; //Fit is now 0.855, no need to draw
  if((TH1F*)f->Get(tune[i]) == nullptr){
    std::cout << tune[i] << " not found!" << std::endl;
    continue;
  }
  std::cout << tune[i] << " " << param[i] << " " << i << "/" << tune.size() << std::endl;
  tmp.push_back((TH1F*)f->Get(tune[i]));
  std::cout << tmp.back()->GetName() << std::endl;
  /*
  for(int i = 0; i < tmp.back()->GetN(); i++) {
    double b = tmp.back()->GetFunction("pol9")->GetParameter(0);
    double x,y;
    tmp.back()->GetPoint(i, x, y);
    //tmp.back()->SetPoint(i, x, y/tmp.back()->Integral());
    tmp.back()->GetY()[i] *= 1./tmp.back()->Integral();
  }
  */
  tmp.back()->SetLineColor(color[i]);
  if(TString(param[i]).Contains("B") || TString(param[i]).Contains("Lambda")) {
    tmp.back()->SetLineColor(color[i-10]);
    tmp.back()->SetLineStyle(kDashDotted);
  }
  //tmp.back()->SetFillColor(color[i]);
  //tmp.back()->SetFillColor(0);
  //if(i==0) tmp->SetTitle("J/#Psi #it{p}_{T} / #Sigma #it{p}_{T}^{ch}");
  tmp.back()->SetTitle(TString::Format("%s",param[i].Data()));
  if(param[i]==0)
    tmp.back()->SetTitle("Fit");
  if(i==0 || tmp.size()==1) {
    tmp.back()->GetXaxis()->SetRangeUser(0.,1.0);
    tmp.back()->GetYaxis()->SetRangeUser(0.4,1.3);
    tmp.front()->GetXaxis()->SetTitle("B #it{p}_{T}^{#it{r}_{b}} / jet #it{p}_{T}");
    tmp.front()->GetXaxis()->SetTitle("#it{x}_{b}");
    //tmp.front()->GetXaxis()->SetTitle("(B #it{p}_{T}^{#it{r}_{b}} / jet #it{p}_{T}) / (B #it{p}_{T}^{#it{r}_{b}=0.855} / jet #it{p}_{T})");
    tmp.front()->GetYaxis()->SetTitle("Weight relative to #it{r}_{b} = 0.855");
    //tmp.front()->GetYaxis()->SetTitle("Jets / 0.05");
    tmp.back()->Draw("hist");
    //tmp.back()->GetYaxis()->SetRangeUser(0.55,1.35);
    tmp.back()->Draw("hist same");
    tdr(tmp.back(), 0, fin);
    tmp.back()->SetTitle(TString::Format("%s",param[0].Data()));
  }
  else tmp.back()->Draw("hist same");
  leg->AddEntry( tmp.back(), tmp.back()->GetTitle(),"f");
}

leg->Draw();
TPaveText *rb = new TPaveText(0.21,0.53,0.59,0.95,"NDC"); //NB blNDC
rb->SetFillStyle(0);
rb->SetTextAlign(11);
rb->SetBorderSize(0);
rb->SetTextFont(42);
rb->SetTextSize(0.036);
TString text = TString::Format("Weight relative to #it{r}_{b} = 0.855");
rb->AddText(text);
//rb->Draw();
rb = new TPaveText(0.41,iniy,0.73,iniy+0.15,"NDC"); //NB blNDC
rb->SetFillStyle(0);
rb->SetTextAlign(11);
rb->SetBorderSize(0);
rb->SetTextFont(42);
rb->SetTextSize(0.036);
text = TString::Format("#it{r}_{b} value");
rb->AddText(text);
rb->Draw();
tmp.front()->SetTitle("");
if(fin) {
c1->SaveAs(TString::Format("LJets2015/2016/mtop/www/meson/morph/genXbweights_%s_final.png", proc.Data()));
c1->SaveAs(TString::Format("LJets2015/2016/mtop/www/meson/morph/genXbweights_%s_final.pdf", proc.Data()));
}
else {
c1->SaveAs(TString::Format("LJets2015/2016/mtop/www/meson/morph/genXbweights_%s.png", proc.Data()));
c1->SaveAs(TString::Format("LJets2015/2016/mtop/www/meson/morph/genXbweights_%s.pdf", proc.Data()));
}
}
