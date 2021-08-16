#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TPad.h>
#include <RooPlot.h>
#include <RooHist.h>
#include "convert.h"
#include "LJets2015/2016/mtop/tdr_sim.C"

void plotDssWeights(bool fin=false) {
std::vector<TString> tune = {"DssuFrag", "DssdFrag", "DssDzuFrag", "DssDzdFrag", "DssDzTenUpFrag", "DssDzTenUpFrag" };
std::vector<float> param = {0.906, 0.810, 1.05, 0.95, 1.1, 0.9};
std::vector<int> color = {kBlue, kCyan-3, kOrange+3, kRed+1, kOrange, kYellow-7, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet+3, kRed-5, kBlack};

//TFile *f = TFile::Open("/eos/cms/store/user/byates/top18012/bfragweights.root");
TFile *f = TFile::Open("data/era2016/bfragweights.root");
std::vector<TGraph*> tmp;
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
float iniy=0.33;
float dy=0.03;
float ndy=param.size()/3;
TLegend *leg = new TLegend(0.21, iniy-dy*ndy, 0.79, iniy+0.05);
//TLegend *leg = new TLegend(0.91, iniy-dy*ndy, 0.99, iniy+0.05);
leg->SetBorderSize(0);
leg->SetFillStyle(1001);      
leg->SetTextFont(43);
leg->SetTextSize(12);
//leg->SetHeader("#it{r}_{b} value", "C");
leg->SetNColumns(7);
for(size_t i=0; i<tune.size(); i++){

  //tmp = ((RooPlot*)f->Get("ptfrac_signal"))->getHist();
  if(tune[i] == TString("fitFrag")) continue; //Fit is now 0.855, no need to draw
  tmp.push_back((TGraph*)f->Get(tune[i]));
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
  tmp.back()->SetFillColor(color[i]);
  //tmp.back()->SetFillColor(0);
  //if(i==0) tmp->SetTitle("J/#Psi #it{p}_{T} / #Sigma #it{p}_{T}^{ch}");
  tmp.back()->SetTitle(TString::Format("%0.3f",param[i]));
  if(param[i]==0)
    tmp.back()->SetTitle("Fit");
  if(tune[i] == TString("fitFrag")) tmp.back()->SetLineStyle(10);
  if(i==0) {
    tmp.back()->GetXaxis()->SetRangeUser(0.,1.0);
    tmp.back()->GetYaxis()->SetRangeUser(0.2,1.35);
    tmp.front()->GetXaxis()->SetTitle("(B #it{p}_{T}^{#it{r}_{b}} / jet #it{p}_{T}) / (B #it{p}_{T}^{#it{r}_{b}=0.855} / jet #it{p}_{T})");
    tmp.front()->GetYaxis()->SetTitle("Weight relative to #it{r}_{b} = 0.855");
    //tmp.front()->GetYaxis()->SetTitle("Jets / 0.05");
    tmp.back()->Draw("ACI same");
    tmp.back()->GetYaxis()->SetRangeUser(0.8,1.2);
    tmp.back()->Draw("ACI same");
    tdr(tmp.back(), 0, fin);
    tmp.back()->SetTitle(TString::Format("%0.3f",param[0]));
  }
  else tmp.back()->Draw("CI same");
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
rb = new TPaveText(0.26,iniy,0.58,iniy+0.15,"NDC"); //NB blNDC
rb->SetFillStyle(0);
rb->SetTextAlign(11);
rb->SetBorderSize(0);
rb->SetTextFont(42);
rb->SetTextSize(0.036);
text = TString::Format("#it{r}_{b} value -- BR_{var} / BR_{nom}");
rb->AddText(text);
rb->Draw();
tmp.front()->SetTitle("");
if(fin) {
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_final.png");
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_final.pdf");
}
else {
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights.png");
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights.pdf");
}
}
