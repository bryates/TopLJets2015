#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TPad.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooRealVar.h>
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/convert.h"
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/tdr.C"

void plotD0(TString sample="172v5", int epoch=0, bool fin=false, float best=0.) {
/*
std::vector<TString> tune = {"_700", "_ddown", "", "_ccentral", "_fit"};
std::vector<float> param = {0.700, 0.800, 0.855,0.900, 0.863};
std::vector<int> color = {kBlue, kCyan-3, kOrange+3, kRed+1, kBlack}; //kOrange, kYellow-7, kViolet-4, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kSpring+9, kViolet, kBlack};
*/
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_central", "_uuup", "_up"};//, "_fit" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.955, 0.975, 1.055, 0};
std::vector<int> color = {kBlue, kCyan-3, kOrange+3, kRed+1, kOrange, kYellow-7, kViolet-4, kGreen+2, kCyan-7, kBlue+2, kMagenta, kOrange+10, kViolet, kBlack};

//RooHist *tmp;
std::vector<TH1F*> tmp;
TString cut("j_pt_ch<150");
//TString sample("172v5");
//int bin(22);
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/param.h"
TCanvas *c1 = setupCanvas();
//TPad *p1 = setupPad();
//p1->cd();
TPad *p1 = new TPad("p1","p1",0.0,0.2,1.0,1.0);
//p1->SetLogy();
p1->SetRightMargin(0.05);
p1->SetLeftMargin(0.12);
p1->SetTopMargin(0.1);
p1->SetBottomMargin(0.12);
p1->Draw();
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

std::vector<TH1F*> fullMC;
TH1F *mf;

for(size_t i=0; i<tune.size(); i++){
TString name("");
if(best>0. && param[i] != best) continue;
name = TString::Format("LJets2015/2016/mtop/sPlot/sPlot/TopMass_%s%s_sPlot_d0.root",sample.Data(),tune[i].Data());
if(name.Contains("fit")) name.ReplaceAll(TString::Format("_%s%s_", sample.Data(), tune[i].Data()), "_172v5_");
name.ReplaceAll("d0.root","d01_xb.root");
TFile *f = TFile::Open(name);
if(f==nullptr) continue;
std::cout << name << std::endl;
/*
RooPlot *ptfrac_mc = (RooPlot*)f->Get("ptfrac_signal");
std::cout << name << std::endl;
TH1F *h1 = (TH1F*)convert(ptfrac_mc,false,bin);
*/
TH1F *h1 = (TH1F*)f->Get("ptfrac_signal_hist");
h1->SetDirectory(0);
std::cout << name << std::endl;
f->Close();
name.ReplaceAll("d01","d02");
f = TFile::Open(name);
if(f==nullptr) continue;
/*
ptfrac_mc = (RooPlot*)f->Get("ptfrac_signal");
std::cout << name << std::endl;
TH1F *h2 = (TH1F*)convert(ptfrac_mc,false,bin);
*/
TH1F *h2 = (TH1F*)f->Get("ptfrac_signal_hist");
h1->SetDirectory(0);
h1->Add(h2);
/*
if(i==0) fullMC = (TH1F*)h1->Clone("fullMC");
else fullMC->Add((TH1F*)h1->Clone());
fullMC->SetDirectory(0);
*/
auto tmph1 = (TH1F*)h1->Clone(TString(h1->GetName()) + "full");
tmph1->SetDirectory(0);
tmph1->SetLineColor(color[i]);
fullMC.push_back(tmph1);
mf = (TH1F*)tmph1->Clone("MCFull");
mf->SetDirectory(0);
h1->Scale(1./h1->Integral());
f->Close();
//ptfrac_hist->plotOn(ptfrac_mc);

  //tmp = ((RooPlot*)f->Get("ptfrac_signal"))->getHist();
  //tmp.push_back((TH1F*)convert(ptfrac_mc,true,0,1.1));
  tmp.push_back(h1);
  //tmp.back()->Rebin(2);
  //tmp.push_back(convert((RooPlot*)f->Get("ptfrac_signal")));
  tmp.back()->SetLineColor(color[i]);
  if(best>0.) tmp.back()->SetLineColor(kRed);
  //if(i==0) tmp->SetTitle("J/#Psi p_{T} / #Sigma p_{T}^{ch}");
  tmp.back()->SetTitle(TString::Format("%0.3f",param[i]));
  if(param[i]==0)
    tmp.back()->SetTitle("Fit");
  if(tune[i] == TString("_fit")) tmp.back()->SetLineStyle(10);
  tmp.back()->Scale(1./tmp.back()->Integral());
  if(i==0) {
    tmp.front()->GetXaxis()->SetTitle("D^{0} p_{T} / #Sigma p_{T}^{ch}");
    tmp.front()->GetYaxis()->SetTitle("Jets / 0.05");
    TString title = tmp.back()->GetTitle();
    tmp.back()->SetTitle("");
    tmp.back()->Draw("hist");
    tmp.back()->SetTitle(title);
    tdr(tmp.back(),epoch,fin);
    tmp.back()->SetTitle(TString::Format("%0.3f",param[0]));
    tmp.back()->GetXaxis()->SetRangeUser(0,1.);
    tmp.back()->GetYaxis()->SetRangeUser(0,0.16);
    if(sample.Contains("Bfrag_genreco"))
         tmp.back()->GetYaxis()->SetRangeUser(0,0.16);
    /*
    if(jpT) tmp.back()->GetYaxis()->SetRangeUser(0,0.16);
    else tmp.back()->GetYaxis()->SetRangeUser(0,0.24);
    */
  }
  else tmp.back()->Draw("same hist");
  tmp.back()->Draw("same e");
  leg->AddEntry( tmp.back(), tmp.back()->GetTitle(),"lp");
}
TString name(TString::Format("LJets2015/2016/mtop/sPlot/sPlot/TopMass_Data_sPlot_d01_xb.root"));
if(sample.Contains("Bfrag_genreco")) name.ReplaceAll("Data", "Data_Bfrag_genreco");
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
RooPlot *ptfrac_data = (RooPlot*)f->Get("ptfrac_signal");
TH1F *d1 = (TH1F*)convert(ptfrac_data,false,bin);
*/
TH1F *d1 = (TH1F*)f->Get("ptfrac_signal_hist");
d1->SetDirectory(0);
f->Close();
name.ReplaceAll("d01","d02");
f = TFile::Open(name);
/*
ptfrac_data = (RooPlot*)f->Get("ptfrac_signal");
std::cout << name << std::endl;
TH1F *d2 = (TH1F*)convert(ptfrac_data,false,bin);
*/
TH1F *d2 = (TH1F*)f->Get("ptfrac_signal_hist");
d2->SetDirectory(0);
d1->Add(d2);
auto df = (TH1F*)d1->Clone("dataFull");
df->SetDirectory(0);
d1->Scale(1./d1->Integral());
f->Close();
//ptfrac_data_hist->plotOn(ptfrac_data);
//tmp.push_back((TH1F*)convert(ptfrac_data,true,0,1.1));
tmp.push_back(d1);
//tmp.push_back(convert((RooPlot*)f->Get("ptfrac_signal")));
tmp.back()->Scale(1./tmp.back()->Integral());
//tmp.back()->Rebin(2);
tmp.back()->SetMarkerStyle(20);
tmp.back()->SetLineColor(1);
tmp.back()->SetMarkerColor(1);
tmp.back()->SetTitle("Data");
leg->AddEntry( tmp.back(), tmp.back()->GetTitle(),"p");
tmp.back()->Draw("same");
gStyle->SetOptStat(0);

c1->cd();
auto p2 = new TPad("p2","p2",0.0,0.0,1.0,0.26);
gStyle->SetPaintTextFormat("%0.2f");
p2->SetBottomMargin(0.4);
p2->SetRightMargin(0.05);
p2->SetLeftMargin(0.12);
p2->SetTopMargin(0.1);
p2->SetGridx(false);
p2->SetGridy(false);
p2->Draw();
p2->cd();
auto ratio = (TH1F*)df->Clone("");
ratio->Divide(df, mf, 1. / df->Integral(), 1. / mf->Integral());
for(int ibin = 1; ibin <= ratio->GetNbinsX(); ibin++) {
  float val = df->GetBinContent(ibin);
  float unc = df->GetBinError(ibin);
  float best = mf->GetBinContent(ibin);
  std::cout << val << "\t" << best << std::endl;
  //ratio->SetBinContent(ibin, 1.);
  //ratio->SetBinError(ibin, unc / val);
}
ratio->SetTitle("");
ratio->GetXaxis()->SetTitle("D^{0} p_{T} / #Sigma p_{T}^{ch}");
ratio->GetYaxis()->SetTitle("Data/MC");
ratio->SetFillColor(920);
ratio->SetMarkerColor(920);
ratio->GetXaxis()->SetNdivisions(10);
ratio->GetYaxis()->SetNdivisions(3);
ratio->GetXaxis()->SetLabelSize(0.1);
ratio->GetXaxis()->SetTitleSize(0.15);
ratio->GetYaxis()->SetLabelSize(0.1);
ratio->GetYaxis()->SetTitleSize(0.15);
ratio->GetYaxis()->SetTitleOffset(0.25);
//ratio->GetXaxis()->SetLabelSize(0.15);
//ratio->GetXaxis()->SetTitleSize(0.2);
//ratio->GetXaxis()->SetTitleOffset(0.8);

/*
ratio->GetXaxis()->SetTitleSize(0.1);
ratio->GetXaxis()->SetLabelSize(0.1);
ratio->GetYaxis()->SetTitleSize(0.1);
ratio->GetYaxis()->SetLabelSize(0.1);
*/
ratio->GetYaxis()->SetRangeUser(.77,1.23);
ratio->GetYaxis()->SetNdivisions(4);
ratio->Draw("e2");
for(auto f : fullMC) {
auto ratio = (TH1F*)d1->Clone("");
ratio->Divide(d1, f, 1, 1. / f->Integral());
ratio->SetTitle("");
ratio->GetXaxis()->SetTitle(d1->GetXaxis()->GetTitle());
ratio->GetYaxis()->SetTitle(TString::Format("Data/MC"));

ratio->SetLineColor(f->GetLineColor());
ratio->SetLineStyle(f->GetLineStyle());
if(TString(f->GetName()).Contains("863")) ratio->SetLineStyle(10);
ratio->Draw("hist same");
}

p1->cd();
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
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_d0" + ext + b + "_d0_final.png");
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_d0" + ext + b + "_d0_final.pdf");
}
else {
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_d0" + ext + b + "_d0.png");
c1->SaveAs("LJets2015/2016/mtop/www/meson/morph/weights_d0" + ext + b + "_d0.pdf");
}
/*
*/

}
