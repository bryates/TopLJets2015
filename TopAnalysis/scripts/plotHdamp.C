#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TPad.h>
#include <RooPlot.h>
#include <RooHist.h>
#include "convert.h"
#include "LJets2015/2016/mtop/tdr.C"

void plotHdamp() {
std::vector<TString> tune = {"LJets2015/2016/mtop/sPlot/sPlot/TopMass_Data_sPlot_d0.root", "LJets2015/2016/mtop/sPlot/sPlot/TopMass_hdampup_sPlot_d0.root", "LJets2015/2016/mtop/sPlot/sPlot/TopMass_hdampdown_sPlot_d0.root"};
std::vector<TString> name = {"Data", "hdamp-up", "hdamp-down"};
/*
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up", "_fit" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055, 0.80};
std::vector<int> color = {kBlue, kCyan-3, kOrange+3, kRed+1, kOrange, kYellow-7, kViolet-4, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet, kBlack};
*/
std::vector<int> color = {kBlack, kBlue, kRed};

//RooHist *tmp;
std::vector<TH1F*> tmp;
//TPad *p1;
TCanvas *c1 = setupCanvas();
TPad *p1 = setupPad();
p1->cd();
float iniy=0.83;
float dy=0.03;
float ndy=name.size();
TLegend *leg = new TLegend(0.2, iniy-dy*ndy, 0.29, iniy+0.05);
leg->SetBorderSize(1);
leg->SetFillStyle(0);      
leg->SetTextFont(43);
leg->SetTextSize(12);
for(size_t i=0; i<tune.size(); i++){
//TString name(TString::Format("LJets2015/2016/mtop/TopMass_fsr-up%s_sPlot.root",tune[i].Data()));
std::cout << tune[i] << std::endl;
TFile *fin = TFile::Open(tune[i]);

  //tmp = ((RooPlot*)fin->Get("ptfrac_mu_tag_signal"))->getHist();
  tmp.push_back(convert((RooPlot*)fin->Get("ptfrac_signal"),false,0,1.1));
  //tmp.push_back(convert((RooPlot*)fin->Get("ptfrac_mu_tag_signal"),false,0,1.1));
  std::cout << "converted" << std::endl;
  tmp.back()->SetLineColor(color[i]);
  std::cout << "line color" << std::endl;
  //if(i==0) tmp->SetTitle("J/#Psi p_{T} / #Sigma p_{T}^{ch}");
  tmp.back()->SetTitle(name[i].Data());
  std::cout << "Title" << std::endl;
  if(tune[i] == TString("_fit")) tmp.back()->SetLineStyle(10);
  tmp.back()->Scale(1./tmp.back()->Integral());
  if(i==0) {
    //tmp.front()->GetXaxis()->SetTitle("J/#Psi p_{T} / #Sigma p_{T}^{ch}");
    tmp.front()->GetXaxis()->SetTitle("D^{0} p_{T} / #Sigma p_{T}^{ch}");
    //tmp.front()->GetXaxis()->SetTitle("(D^{0} p_{T} + #mu p_{T}) / #Sigma p_{T}^{ch}");
    tmp.front()->GetYaxis()->SetTitle("Jets / 0.05");
    tmp.back()->GetXaxis()->SetRangeUser(0.,1.1);
    //tmp.back()->GetYaxis()->SetRangeUser(0.,0.11);//1150);
    tmp.back()->Draw();
    tdr(tmp.back());
    tmp.back()->SetTitle(name[i].Data());
  }
  else tmp.back()->Draw("same hist");
  leg->AddEntry( tmp.back(), tmp.back()->GetTitle(),"lp");
}


//tdr(tmp.front());
leg->Draw();
tmp.front()->SetTitle("");
}
