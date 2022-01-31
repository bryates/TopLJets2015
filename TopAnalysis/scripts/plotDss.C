#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/mtop/tdr_sim.C"

void plotDss(bool norm=true, bool ratio=false) {
TFile *fin = new TFile("LJets2015/2016/bfrag/xb_fit.root");
std::vector<TString> hadron = {"inc", "noDss", "Dss"};//,  "411","10411","10421","10413","10423","20413","20423","415","425","10431","433","10433","20433","435"};//, "DssDz", "DssDzb"};
std::vector<TString> name   = {"Inc.", "No D^{**}", "D^{**}", "D^{**} #rightarrow D^{0}", "D^{**} #rightarrow #bar{D^{0}}"};
std::map<TString, TString> pdgids;
pdgids["inc"] = "Inc.";
pdgids["noDss"] = "No D^{**}";
pdgids["Dss"] = "D^{**}";
pdgids["DssDz"] = "D^{**} #rightarrow D^{0}";
pdgids["DssDzb"] = "D^{**} #rightarrow #bar{D^{0}}";
pdgids["411"] = "411";
pdgids["10411"] = "10411";
pdgids["10421"] = "10421";
pdgids["10413"] = "10413";
pdgids["10423"] = "10423";
pdgids["20413"] = "20413";
pdgids["20423"] = "20423";
pdgids["415"] = "415";
pdgids["425"] = "425";
pdgids["10431"] = "10431";
pdgids["433"] = "433";
pdgids["10433"] = "10433";
pdgids["20433"] = "20433";
pdgids["435"] = "435";
std::vector<int> color = {kBlack, kRed, kCyan-3, kOrange+3, kBlue+1, kCyan}; //kOrange, kYellow-7, kViolet-4, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet, kBlack};
//std::vector<int> color = {kBlack, kBlue, kCyan-3, kOrange+3, kRed+1, kCyan}; //kOrange, kYellow-7, kViolet-4, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet, kBlack};
c1 = setupCanvas();
setupPad()->cd();
gStyle->SetOptStat(0);
TLegend *leg = new TLegend(0.15, 0.5, 0.65, 0.8);
leg->SetTextFont(43);
leg->SetTextSize(22);
leg->SetBorderSize(0);
float sum = 0;
for(size_t ibin = 0; ibin < hadron.size(); ibin++) {
  TH1F* h;
  if(hadron[ibin]=="inc")
    h = (TH1F*)fin->Get("bfragAnalysis/xb_inc")->Clone();
  else
    h = (TH1F*)fin->Get(TString::Format("bfragAnalysis/xb_%s", hadron[ibin].Data()))->Clone();
  std::cout << hadron[ibin] << std::endl;
  h->SetDirectory(0);
  h->GetXaxis()->SetRangeUser(0, 1.01);
  //if(norm) h->SetFillColor(color[ibin]);
  h->SetLineColor(color[ibin]);
  h->SetMarkerColor(color[ibin]);
  h->SetLineWidth(2);
  //if(ibin==0) h->SetLineWidth(6);
  if(ibin>=name.size()) name.push_back(hadron[ibin]);
      //h->SetTitle(hadron[ibin]);
  h->SetTitle(name[ibin]);
  if(ibin == 0) {
    h->GetYaxis()->SetRangeUser(0,h->GetMaximum()*1.1);
    if(norm) h->GetYaxis()->SetTitle("1/N dN/d#it{x}_{b}");
    else h->GetYaxis()->SetTitle("dN/d#it{x}_{b}");
    tdr(h);
    //h->DrawNormalized("hist c");
    if(norm) h->DrawNormalized("hist c");
    else h->Draw("hist c");
    tdr(h);
  }
  else if(norm) h->DrawNormalized("same hist c");
  else h->Draw("same hist c");
  if(ratio) {
  if(ibin==1) h->Divide(h, (TH1F*)fin->Get("bfragAnalysis/xb_inc")->Clone(), 1/h->Integral(), 1/((TH1F*)fin->Get("bfragAnalysis/xb_inc"))->Integral());
  tdr(h);
  h->GetYaxis()->SetRangeUser(0.8,1.2);
  h->Draw("hist c");
  tdr(h);
  }
  std::cout << pdgids[hadron[ibin]] << "=" << h->Integral() << std::endl;
  if(ibin != 0)
  sum += h->Integral();
  if(ibin==1 || !ratio) leg->AddEntry(h, TString::Format("%s", pdgids[hadron[ibin]].Data()), "f");
  //else h->DrawNormalized("same hist c");
}
std::cout << "Total=" << sum << std::endl;
leg->Draw();

}
