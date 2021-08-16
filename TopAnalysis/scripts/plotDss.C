#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "LJets2015/2016/mtop/tdr_sim.C"

void plotDss(bool norm=true, bool ratio=false) {
TFile *fin = new TFile("LJets2015/2016/bfrag/xb_fit.root");
std::vector<TString> hadron = {"inc", "noDss", "Dss"};//, "DssDz", "DssDzb"};
std::vector<TString> name   = {"Inc.", "No D^{**}", "D^{**}", "D^{**} #rightarrow D^{0}", "D^{**} #rightarrow #bar{D^{0}}"};
std::map<TString, TString> pdgids;
pdgids["inc"] = "Inc.";
pdgids["noDss"] = "No D^{**}";
pdgids["Dss"] = "D^{**}";
pdgids["DssDz"] = "D^{**} #rightarrow D^{0}";
pdgids["DssDzb"] = "D^{**} #rightarrow #bar{D^{0}}";
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
  h->SetDirectory(0);
  h->GetXaxis()->SetRangeUser(0, 1.01);
  //if(norm) h->SetFillColor(color[ibin]);
  h->SetLineColor(color[ibin]);
  h->SetMarkerColor(color[ibin]);
  h->SetLineWidth(2);
  //if(ibin==0) h->SetLineWidth(6);
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
