#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "LJets2015/2016/mtop/tdr_sim.C"

void plotXb(bool norm=true) {
TFile *fin = new TFile("LJets2015/2016/bfrag/xb_fit.root");
std::vector<int> hadron        = {0,521,511,531,5122};
std::vector<float> hadronUncDz = {0,1.5,0.7,0.5,1.1};
std::vector<TString> name = {"Inc.", "B^{#pm}", "B^{0}", "B^{0}_{s}", "#Lambda^{0}_{b}"};
std::map<int, TString> pdgids;
pdgids[0]   = "Inc.";
pdgids[511] = "B^{0}";
pdgids[521] = "B^{#pm}";
pdgids[531] = "B^{0}_{s}";
pdgids[5122] = "#Lambda_{b}";
std::vector<int> color = {kBlack, kBlue, kCyan-3, kOrange+3, kRed+1, kCyan}; //kOrange, kYellow-7, kViolet-4, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet, kBlack};
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
  if(hadron[ibin]==0)
    h = (TH1F*)fin->Get("bfragAnalysis/xb_inc")->Clone();
  else
    h = (TH1F*)fin->Get(TString::Format("bfragAnalysis/xb_%d", hadron[ibin]))->Clone();
  h->SetDirectory(0);
  h->GetXaxis()->SetRangeUser(0, 1.01);
  //if(norm) h->SetFillColor(color[ibin]);
  h->SetLineColor(color[ibin]);
  h->SetMarkerColor(color[ibin]);
  h->SetLineWidth(2);
  if(ibin==1) h->SetLineWidth(6);
  h->SetTitle(name[ibin]);
  for(int i = 1; i <= h->GetNbinsX(); i++) {
    h->SetBinError(i, h->GetBinContent(i)*sqrt(1/h->GetBinContent(i) + hadronUncDz[ibin]/h->GetBinContent(i)));
  }
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
  std::cout << pdgids[hadron[ibin]] << "=" << h->Integral() << std::endl;
  if(ibin != 0)
  sum += h->Integral();
  leg->AddEntry(h, TString::Format("%s", pdgids[hadron[ibin]].Data()), "f");
  //else h->DrawNormalized("same hist c");
}
std::cout << "Total=" << sum << std::endl;
leg->Draw();

}
