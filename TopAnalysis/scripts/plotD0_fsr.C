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

void plotD0_fsr() {
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
TString name(TString::Format("LJets2015/2016/mtop/sPlot/sPlot/TopMass_172v5%s_sPlot_d0.root",tune[i].Data()));
if(tune[i].Contains("fsr"))
  name=TString::Format("LJets2015/2016/mtop/sPlot/sPlot/TopMass_%s_sPlot_d0.root",tune[i].Data());
std::cout << name << std::endl;
TFile *fin = TFile::Open(name);

  //tmp = ((RooPlot*)fin->Get("ptfrac_signal"))->getHist();
  tmp.push_back(convert((RooPlot*)fin->Get("ptfrac_signal"),false,0,1.1));
  std::cout << "converted" << std::endl;
  tmp.back()->SetLineColor(color[i]);
  std::cout << "line color" << std::endl;
  //if(i==0) tmp->SetTitle("J/#Psi p_{T} / #Sigma p_{T}^{ch}");
  tmp.back()->SetTitle(TString::Format("%0.3f",param[i]));
  std::cout << "Title" << std::endl;
  if(tune[i] == "") tmp.back()->SetTitle("Nominal");
  if(tune[i] == "_fit") tmp.back()->SetTitle("Fit");
  if(tune[i] == TString("_fit")) tmp.back()->SetLineStyle(10);
  tmp.back()->Scale(1./tmp.back()->Integral());
  if(i==0) {
    tmp.back()->GetXaxis()->SetRangeUser(0.,1.1);
    tmp.back()->GetYaxis()->SetRangeUser(0.,0.11);//1150);
    tmp.back()->Draw("hist");
  }
  else tmp.back()->Draw("same hist");
}
TFile *fin = TFile::Open("LJets2015/2016/mtop/sPlot/sPlot/TopMass_Data_sPlot_d0.root");
tmp.push_back(convert((RooPlot*)fin->Get("ptfrac_signal"),false,0,1.1));
tmp.back()->Scale(1./tmp.back()->Integral());
//tmp.push_back(convert((RooPlot*)fin->Get("ptfrac_signal")));
//tmp.back()->Rebin(2);
tmp.back()->SetLineColor(1);
tmp.back()->SetMarkerStyle(20);
tmp.back()->SetTitle("Data");
tmp.back()->Draw("same");
gStyle->SetOptStat(0);

TLegend *leg = c1->BuildLegend();
tmp.front()->SetTitle("D^{0} p_{T} / #Sigma p_{T}^{ch}");
tmp.front()->GetXaxis()->SetTitle("D^{0} p_{T} / #Sigma p_{T}^{ch}");
tmp.front()->GetYaxis()->SetTitle("Jets / 0.05");
tdr(tmp.front());
}
