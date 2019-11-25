#include "TCanvas.h"
#include "TString.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFit.h"
#include "TH1F.h"
#include "TF1.h"
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/convert.h"
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/tdr.C"

void comb() {
TCanvas *c1 = setupCanvas();
TPad *p1 = setupPad();
p1->cd();
TH1F *mc;
gStyle->SetOptStat(0);
//std::vector<int> type = {25, 50, 75, 100};
//std::vector<float> type = {0.04, 0.06, 0.1};
//std::vector<int> type = {40, 60,0, 100
//std::vector<int> type = {425, 575, 800, 100};
//std::vector<int> type = {425, 575, 650, 725, 800, 875, 1000};
//std::vector<float> type = {0., 0.425, 0.65, 0.8, 1.};
//std::vector<int> type = {375, 600, 825, 1000};
//std::vector<int> type = {30, 45, 80, 100};
//std::vector<int> type = {30, 45, 70, 100};
//std::vector<float> bin = {-0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.875, 0.95, 1.0};
//std::vector<float> bin = {0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.875, 0.95, 1.0};
//std::vector<float> bin = {0, 0.075, 0.15, 0.225, 0.3, 0.45, 0.525, 0.6, 0.675, 0.75, 0.825, 0.9, 0.975, 1.0};
//bin = {-0.03, 0.045, 0.12, 0.195, 0.27, 0.345, 0.42, 0.495, 0.57, 0.645, 0.72, 0.795, 0.87, 0.945, 1.0};
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/param.h"
//bin = {0.2, 0.425, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
RooBinning bins(0,1.1);
for(int i = 0; i < bin.size(); i++) {
 bins.addBoundary(bin[i]);
}
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_up" };
//tune = { "_fit", "_fitup", "_fitdown" };
for(size_t s = 0; s < tune.size(); s++) {
TH1F *tot = nullptr;
TH1F *totm = nullptr;
TH1F *ptep1 = nullptr;
TH1F *ptep2 = nullptr;
//RooDataSet *ds;
RooRealVar *nsig_sw;
RooRealVar *nsig;
RooRealVar *nbkg;
//size_t s(7);
//if(s != 7) continue; //FIXME uncomment for Data
//if(tune[s] != "_700") continue;
for(size_t i = 1; i < type.size(); i++) {
  TH1F *ptep = nullptr;
  TH1F *massp = nullptr;
  for(int epoch = 1; epoch <= 2; epoch++) {
    //if(epoch == 2 && i == 0) continue;
    //if(epoch ==2 && type[i] <= 0.575) continue;
    TString fUrl(TString::Format("sPlot/sPlot//TopMass_172v5%s_sPlot_d0%d_%d.root",tune[s].Data(), epoch, int(type[i]*1000)));
    //fUrl.ReplaceAll("d01", "d0");
   /*
    if(epoch==1 && type[i]<=0.5) 
      fUrl = TString::Format("sPlot/sPlot//TopMass_172v5%s_sPlot_d0_%d.root",tune[s].Data(), int(type[i]*1000));
   */
    std::cout << fUrl << std::endl;
    TFile *f = TFile::Open(fUrl);
    if(f == nullptr) return;
    auto u = (RooWorkspace*)f->Get("w");
    auto sig = (RooDataSet*)u->data("sigData");
    nsig = u->var("nsig");
    /*
    if(ds == nullptr)
      //ds = (RooDataSet*)((RooDataSet*)((RooWorkspace*)f->Get("w"))->data("sigData"))->Clone();
      ds = (RooDataSet*)sig->Clone();
    else
      //ds->append(*(RooDataSet*)((RooDataSet*)((RooWorkspace*)f->Get("w"))->data("sigData"))->Clone());
      ds->append(*(RooDataSet*)sig->Clone());
    */
    RooPlot *ptfrac = (RooPlot*)f->Get("ptfrac_signal");
    ((RooDataSet*)((RooWorkspace*)f->Get("w"))->data("sigData"))->plotOn(ptfrac, RooFit::Binning(bins), RooFit::DataError(RooAbsData::SumW2));
    //((RooDataSet*)((RooWorkspace*)f->Get("w"))->data("sigData"))->plotOn(ptfrac, RooFit::DataError(RooAbsData::SumW2));
    //((RooDataSet*)((RooWorkspace*)f->Get("w"))->data("sigData"))->plotOn(ptfrac, RooFit::DataError(RooAbsData::SumW2));
    ptfrac->Draw();
    if(ptfrac == nullptr) return;
    TH1F *h = (TH1F*)convert(ptfrac, false, bin);
    h->SetDirectory(0);
    //if(type[i]>0.5) h->Scale(((RooWorkspace*)f->Get("w"))->var("nsig")->getVal() / h->Integral());
    if(tot == nullptr) tot = (TH1F*)h->Clone("ptfrac_signal_hist");
    else tot->Add(h);
    tot->SetDirectory(0);
    if(epoch == 1 && ptep1 == nullptr) ptep1 = (TH1F*)h->Clone(TString("proxy_BCDEF"));
    else if(epoch == 1) ptep1->Add((TH1F*)h->Clone(TString("ptep%d",epoch)));
    if(epoch == 2 && ptep2 == nullptr) ptep2 = (TH1F*)h->Clone(TString("proxy_GH"));
    else if(epoch == 2) ptep2->Add((TH1F*)h->Clone(TString("ptep%d",epoch)));
    if(ptep == nullptr) ptep = (TH1F*)h->Clone(TString("ptep"));
    else ptep->Add(h);
    ptep->SetDirectory(0);
    if(epoch==1 && ptep1 != nullptr) ptep1->SetDirectory(0);
    if(epoch==2 && ptep2 != nullptr) ptep2->SetDirectory(0);
    if(epoch==1 && ptep1 != nullptr) ptep1->Draw();
    else if(epoch==2 && ptep2 != nullptr) ptep2->Draw();
    c1->SaveAs(TString::Format("www/meson/tdr/ptfrac_signal_172v5%s_%d-%d_d0%d_noSumCh.png", tune[s].Data(), int(type[i-1]*1000), int(type[i]*1000), epoch));
    c1->SaveAs(TString::Format("www/meson/tdr/ptfrac_signal_172v5%s_%d-%d_d0%d_noSumCh.pdf", tune[s].Data(), int(type[i-1]*1000), int(type[i]*1000), epoch));
    //if(epoch==1) std::cout << ptep1->Integral() << std::endl;
    RooPlot *mass = (RooPlot*)f->Get("massD0");
    mass->SetTitle("");
    tdr(mass, epoch);
    mass->Draw();
    tdr(mass, epoch);
    c1->SaveAs(TString::Format("www/meson/tdr/massD0_172v5%s_%d-%d_d0%d_noSumCh.png", tune[s].Data(), int(type[i-1]*1000), int(type[i]*1000), epoch));
    c1->SaveAs(TString::Format("www/meson/tdr/massD0_172v5%s_%d-%d_d0%d_noSumCh.pdf", tune[s].Data(), int(type[i-1]*1000), int(type[i]*1000), epoch));
    std::cout << "mass" << std::endl;
    RooWorkspace *w = (RooWorkspace*)f->Get("w");
    RooDataSet *dsn = (RooDataSet*)w->data("dsSWeights");
    RooDataSet *dsS = new RooDataSet("ds", "ds", dsn, *dsn->get(), 0, "weight");
    mass = w->var("d0_mass")->frame();
    nsig_sw = (RooRealVar*)w->var("nsig_sw");
    dsS->plotOn(mass, RooFit::Binning(60), RooFit::DataError(RooAbsData::SumW2));
    h = (TH1F*)convert(mass, false, 1.7, 2);
    std::cout << "n mass: " << h->Integral() << std::endl;
    h->SetDirectory(0);
    if(epoch == 1) massp = h;//FIXME
    else massp->Add(h);
    massp->SetDirectory(0);
    if(tune[s] == "") {
    if(totm == nullptr) totm = (TH1F*)h->Clone();
    else totm->Add(h);
    totm->SetDirectory(0);
    }
    f->Close();
    delete f;
  }
  ptep->Draw();
  c1->SaveAs(TString::Format("www/meson/tdr/ptfrac_signal_172v5%s_%d-%d_d0.png", tune[s].Data(), int(type[i-1]*1000), int( type[i]*1000)));
  c1->SaveAs(TString::Format("www/meson/tdr/ptfrac_signal_172v5%s_%d-%d_d0.pdf", tune[s].Data(), int(type[i-1]*1000), int( type[i]*1000)));
  char c;
  delete ptep;
  ptep = nullptr;
  massp->SetMinimum(0);
  massp->Draw();
  c1->SaveAs(TString::Format("www/meson/tdr/massD0_172v5%s_%d-%d_d0.png", tune[s].Data(), int(type[i-1]*1000), int( type[i]*1000)));
  c1->SaveAs(TString::Format("www/meson/tdr/massD0_172v5%s_%d-%d_d0.pdf", tune[s].Data(), int(type[i-1]*1000), int( type[i]*1000)));
  delete massp;
  massp = nullptr;
}

if(tune[s] == "") {
  totm->SetMinimum(0);
  tot->GetYaxis()->SetRangeUser(0,1400);
totm->Draw();
    std::cout << "n tot mass: " << totm->Integral() << std::endl;
//fnc->Draw("same");
    c1->SaveAs(TString::Format("www/meson/tdr/massD0_172v5_d0.png"));
    c1->SaveAs(TString::Format("www/meson/tdr/massD0_172v5_d0.pdf"));
}
tot->Draw();
    c1->SaveAs(TString::Format("www/meson/tdr/ptfrac_signal_172v5%s_d0.png", tune[s].Data()));
    c1->SaveAs(TString::Format("www/meson/tdr/ptfrac_signal_172v5%s_d0.pdf", tune[s].Data()));
  if(ptep1!=nullptr) {
  ptep1->SetTitle(ptep1->GetName());
  ptep1->SetMinimum(0);
  tdr(ptep1,1);
  ptep1->GetXaxis()->SetRangeUser(0,1.1);
  ptep1->Draw();
  tdr(ptep1,1);
  c1->SaveAs(TString::Format("www/meson/tdr/ptfrac_signal_172v5%s_BCDEF_d0.png", tune[s].Data()));
  c1->SaveAs(TString::Format("www/meson/tdr/ptfrac_signal_172v5%s_BCDEF_d0.pdf", tune[s].Data()));
  delete ptep1;
  ptep1 = nullptr;
  }
  if(ptep2!=nullptr) {
  ptep2->SetMinimum(0);
  ptep2->Draw();
  ptep2->SetTitle(ptep2->GetName());
  c1->SaveAs(TString::Format("www/meson/tdr/ptfrac_signal_172v5%s_GH_d0.png", tune[s].Data()));
  c1->SaveAs(TString::Format("www/meson/tdr/ptfrac_signal_172v5%s_GH_d0.pdf", tune[s].Data()));
  delete ptep2;
  ptep2 = nullptr;
  }
if(tune[s] == "") mc = (TH1F*)tot->Clone();
TFile *fout = TFile::Open(TString::Format("sPlot/sPlot//TopMass_172v5%s_sPlot_d0_xb.root",tune[s].Data()),"RECREATE");
tot->SetDirectory(fout);
tot->SetMaximum(2000);
tot->Write();
tot->SetDirectory(0);
/*
RooWorkspace w("w","w");
w.import(*ds);
//w.import(*nsig_sw);
//w.import(*nsig);
//w.import(*nbkg);
w.Write();
*/
fout->Close();
delete fout;
delete tot;
delete totm;
//delete ds;
tot = nullptr;
totm = nullptr;
}

/*
TH1F *h;
for(int epoch = 1; epoch < 2; epoch++) {
  TString fUrl(TString::Format("sPlot/sPlot//TopMass_172v5_sPlot_d0%d.root", epoch));
  std::cout << fUrl << std::endl;
  TFile *f = new TFile(fUrl);
  if(f == nullptr) return;
  RooPlot *ptfrac = (RooPlot*)f->Get("ptfrac_signal");
  if(epoch == 1) h = (TH1F*)convert(ptfrac, false, bin);
  else h->Add((TH1F*)convert(ptfrac, false, bin));
}
h->Divide(mc);
p1->SetLogy();
h->Draw();
tdr(h);
h->Draw();
tdr(h);
*/
}
