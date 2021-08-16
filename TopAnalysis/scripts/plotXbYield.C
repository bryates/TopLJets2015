#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/mtop/tdr_sim.C"

void plotXbYield(bool norm=false) {
TFile *fin = new TFile("LJets2015/2016/bfrag/xb_fit.root");
std::vector<int> hadron        = {521,511,531,5122};
std::vector<float> hadronUncDz = {1.5,0.7,0.5,1.1};
std::vector<TString> name = {"B^{#pm}", "B^{0}", "B^{0}_{s}", "#Lambda^{0}_{b}"};
std::map<int, TString> pdgids;
pdgids[0]   = "Inc.";
pdgids[511] = "B^{0}";
pdgids[521] = "B^{#pm}";
pdgids[531] = "B^{0}_{s}";
pdgids[5122] = "#Lambda_{b}";
std::vector<int> color = {kBlue, kCyan-3, kOrange+3, kRed+1, kCyan}; //kOrange, kYellow-7, kViolet-4, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet, kBlack};
c1 = setupCanvas();
setupPad()->cd();
gStyle->SetOptStat(0);
TLegend *leg = new TLegend(0.15, 0.7, 0.65, 0.8);
leg->SetTextFont(43);
leg->SetTextSize(22);
leg->SetBorderSize(0);
THStack *hs = new THStack("hs", "hs");
for(int ibin = hadron.size()-1; ibin >= 0; ibin--) {
  TH1F* h;
  /*
  if(hadron[ibin]==0)
    h = (TH1F*)fin->Get("bfragAnalysis/xb_semilepDzuinc")->Clone();
  else
  */
    auto up = (TH1F*)fin->Get(TString::Format("bfragAnalysis/xb_semilepDzu%d", hadron[ibin]))->Clone();
    auto down = (TH1F*)fin->Get(TString::Format("bfragAnalysis/xb_semilepDzd%d", hadron[ibin]))->Clone();
    h = (TH1F*)up->Clone();
  h->Clear();
  h->SetDirectory(0);
  h->GetXaxis()->SetRangeUser(0, 1.01);
  h->GetYaxis()->SetRangeUser(0, 1e4);
  h->SetFillColor(color[ibin]);
  h->SetLineColor(color[ibin]);
  h->SetMarkerColor(color[ibin]);
  h->SetTitle(name[ibin]);
  for(int i = 1; i <= h->GetNbinsX(); i++) {
    float vup = up->GetBinContent(i);
    float eup = up->GetBinError(i);
    float vdown = down->GetBinContent(i);
    float edown = down->GetBinContent(i);
    h->SetBinContent(i, (vup-vdown)/2);
    h->SetBinError(i, sqrt((vup*vup+vdown*vdown)/2));
  }
  auto htmp = new TH1F("htmp", "htmp", 50, 0, 1);
  for(int i = 1; i <= htmp->GetNbinsX(); i++) {
    htmp->SetBinContent(i, h->GetBinContent(i));
  }
  htmp->SetFillColor(color[ibin]);
  htmp->SetLineColor(color[ibin]);
  htmp->SetMarkerColor(color[ibin]);
  htmp->SetTitle(name[ibin]);
  hs->Add(htmp);
  if(ibin == 0) {
    h->GetYaxis()->SetRangeUser(0,h->GetMaximum()*1.1);
    if(norm) h->GetYaxis()->SetTitle("1/N dN/d#it{x}_{b}");
    else h->GetYaxis()->SetTitle("dN/d#it{x}_{b}");
    tdr(h);
    //hs->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
    h->DrawNormalized("hist");
    tdr(h);
  }
  else h->Draw("hist same");
  leg->AddEntry(h, TString::Format("%s", pdgids[hadron[ibin]].Data()));
  //else h->DrawNormalized("same hist c");
}
leg->SetNColumns(hadron.size());
//hs->GetYaxis()->SetTitle("dN/d#it{x}_{b}");
hs->Draw();
hs->SetMaximum(hs->GetMaximum()*10);
hs->SetMinimum(1e-1);
gPad->SetLogy();
tdr(hs);
hs->GetXaxis()->SetTitle("#it{x}_{b}=#it{p}_{T}(B)/#it{p}_{T}(jet)");
hs->GetYaxis()->SetTitle("dN/d#it{x}_{b}");
leg->Draw();

}
