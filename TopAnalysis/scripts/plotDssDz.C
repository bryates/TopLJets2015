#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/mtop/tdr_sim.C"

void plotDssDz(bool norm=true, bool ratio=false) {
TFile *fin = new TFile("LJets2015/2016/bfrag/xb_fit.root");
TFile *fup = new TFile("LJets2015/2016/bfrag/xb_fitup.root");
TFile *fdown = new TFile("LJets2015/2016/bfrag/xb_fitdown.root");
std::vector<TString> hadron = {"noDssDz", "DssDz"};//,  "411","10411","10421","10413","10423","20413","20423","415","425","10431","433","10433","20433","435"};//, "DssDz", "DssDzb"};
//std::vector<TString> hadron = {"inc", "noDss", "Dss"};//,  "411","10411","10421","10413","10423","20413","20423","415","425","10431","433","10433","20433","435"};//, "DssDz", "DssDzb"};
std::vector<TString> name   = {"D^{**}", "D^{**} #rightarrow D^{0}", "D^{**} #rightarrow #bar{D^{0}}"};
std::map<TString, TString> pdgids;
pdgids["inc"] = "B #rightarrow inc.";
pdgids["noDssDz"] = "B #rightarrow D^{0}";
pdgids["DssDz"] = "B #rightarrow D^{**} #rightarrow D^{0}";
pdgids["DssDzu"] = "B #rightarrow D^{**} #rightarrow D^{0} up";
pdgids["DssDzd"] = "B #rightarrow D^{**} #rightarrow D^{0} down";
//pdgids["DssDz"] = "D^{**} #rightarrow D^{0}";
pdgids["DssDzb"] = "B #rightarrow D^{**} #rightarrow #bar{D^{0}}";
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
std::vector<TH1*> hists;
std::vector<TH1*> histsu;
std::vector<TH1*> histsd;
float sum = 0;
for(size_t ibin = 0; ibin < hadron.size(); ibin++) {
  TH1F* h,*hup,*hdown,*hnodss;
  if(hadron[ibin]=="noDssDz") {
    h = (TH1F*)fin->Get("bfragAnalysis/xb_noDssDz")->Clone();
    hnodss = (TH1F*)h->Clone();
    hup = (TH1F*)fup->Get("bfragAnalysis/xb_noDssDz")->Clone();
    hdown = (TH1F*)fdown->Get("bfragAnalysis/xb_noDssDz")->Clone();
  }
  else {
    h = (TH1F*)fin->Get(TString::Format("bfragAnalysis/xb_%s", hadron[ibin].Data()))->Clone();
    hup = (TH1F*)fup->Get(TString::Format("bfragAnalysis/xb_%s", hadron[ibin].Data()))->Clone();
    hdown = (TH1F*)fdown->Get(TString::Format("bfragAnalysis/xb_%s", hadron[ibin].Data()))->Clone();
    std::cout << "Ratio = " << h->Integral() / hnodss->Integral() << std::endl;
  }
  h->SetDirectory(0);
  h->GetXaxis()->SetRangeUser(0, 1.01);
  hup->GetXaxis()->SetRangeUser(0, 1.01);
  hdown->GetXaxis()->SetRangeUser(0, 1.01);
  //if(norm) h->SetFillColor(color[ibin]);
  h->SetLineColor(color[ibin]);
  h->SetMarkerColor(color[ibin]);
  h->SetLineWidth(2);
  h->SetLineWidth(6);
  hup->SetLineWidth(4);
  hdown->SetLineWidth(2);
  //if(ibin==0) h->SetLineWidth(6);
  if(ibin>=name.size()) name.push_back(hadron[ibin]);
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
  if(ibin>0) h->Divide(h, (TH1F*)fin->Get("bfragAnalysis/xb_noDssDz")->Clone(), 1/h->Integral(), 1/((TH1F*)fin->Get("bfragAnalysis/xb_noDssDz"))->Integral());
  if(ibin>0) hup->Divide(hup, (TH1F*)fup->Get("bfragAnalysis/xb_noDssDz")->Clone(), 1/hup->Integral(), 1/((TH1F*)fup->Get("bfragAnalysis/xb_noDssDz"))->Integral());
  if(ibin>0) hdown->Divide(hdown, (TH1F*)fdown->Get("bfragAnalysis/xb_noDssDz")->Clone(), 1/hdown->Integral(), 1/((TH1F*)fdown->Get("bfragAnalysis/xb_noDssDz"))->Integral());
  if(ibin>0) {
    hists.push_back((TH1F*)h->Clone());
    histsu.push_back((TH1F*)hup->Clone());
    histsd.push_back((TH1F*)hdown->Clone());
  }
  hup->SetLineColor(color[ibin]+2);
  hdown->SetLineColor(color[ibin]-2);
  tdr(h);
  h->GetYaxis()->SetRangeUser(0.6,1.4);
  h->Draw("hist c");
  hup->Divide(h);
  hdown->Divide(h);
  hup->Draw("same hist c");
  hup->GetYaxis()->SetRangeUser(0.8,1.2);
  hup->Draw("hist c");
  hdown->Draw("same hist c");
  tdr(h);
  }
  std::cout << pdgids[hadron[ibin]] << "=" << h->Integral() << std::endl;
  if(ibin != 0)
  sum += h->Integral();
  if(ibin>0 || !ratio) leg->AddEntry(h, TString::Format("%s", pdgids[hadron[ibin]].Data()), "f");
  if(ibin>0 || !ratio) leg->AddEntry(hup, TString::Format("%s", pdgids[hadron[ibin]].Data()), "f");
  if(ibin>0 || !ratio) leg->AddEntry(hdown, TString::Format("%s", pdgids[hadron[ibin]].Data()), "f");
  //else h->DrawNormalized("same hist c");
}
std::cout << "Total=" << sum << std::endl;
leg->Draw();

}
