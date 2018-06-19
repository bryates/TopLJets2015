#include <RooFit.h>
#include "RooGlobalFunc.h"
#include "RooWorkspace.h"
#include "RooHistPdf.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooIntegralMorph.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooAddition.h"
#include "RooArgSet.h"
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/convert.h"

using namespace RooFit;

TString report("");
float chi2_sim_test(TString tune="", TString sample="d0_mu_tag_mu");
std::vector<RooChi2Var> chi;
RooRealVar ptfrac;
int bins(22);
float chi2vals[12];

void chi2_sim() {
//std::vector<TString> samples = { "_d0_mu_tag_mu", "_jpsi" };
std::vector<TString> samples = { "_d0_mu_tag_mu", "_d0", "_jpsi" };
/*
std::vector<TString> tunes = {"_down", "", "_cccentral", "_central", "_up" };
std::vector<float> param = {0.755, 0.855, 0.875, 0.955, 1.055};
*/
std::vector<TString> tunes = {"_sdown", "_down", "_scentral", "", "_cccentral", "_central", "_up" };
std::vector<float> param = {0.655, 0.755, 0.825, 0.855, 0.875, 0.955, 1.055};
/*
std::vector<TString> tunes = {"_down", "_ddown", "_dddown", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.755, 0.775, 0.800, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055};
*/

TCanvas *c1 = new TCanvas("c1","c1");
TH1F *chiTest = new TH1F("chiTest","#chi^{2} test",100,0,2);
for(size_t i = 0; i < tunes.size(); i++) {
  float chi(0.);
  for(size_t j = 0; j < samples.size(); j++) {
    chi2vals[i + j * samples.size()] = chi2_sim_test(tunes[i], samples[j]);
    chi += chi2vals[i + j * samples.size()];
  }
  chiTest->SetBinContent(chiTest->FindBin(param[i]),chi);
}

for(size_t i = 0; i < tunes.size(); i++) {
  for(size_t j = 0; j < samples.size(); j++) {
    std::cout << chi2vals[i + j * samples.size()] << "  ";
  }
  std::cout << std::endl;
}

/*
for(size_t i = 1; i < tunes.size(); i++) {
  float chi(chi2vals[i-1]);
  for(size_t j = 1; j < samples.size(); j++) {
    chi += chi2vals[i + j * samples.size()];
  }
  chiTest->SetBinContent(chiTest->FindBin(param[i]),chi);
}
for(size_t i = 0; i < samples.size(); i++) {
  float chi(chi2vals[i]);
  for(size_t j = 1; j < tunes.size(); j++) { //D^0
    for(size_t k = 1; k < tunes.size(); k++) { //D^0_mu
      for(size_t l = 1; l < tunes.size(); l++) { //J/Psi
        chi += chi2vals[];
      }
    }
  }
  chiTest->SetBinContent(chiTest->FindBin(param[i]),chi);
}
*/
      
/*
for(size_t i = 0; i < tunes.size(); i++) {
  float chi = chi2_sim_test(tunes[i], samples[0]);
  for(size_t j = 0; j < tunes.size(); j++) {
    chi += chi2_sim_test(tunes[j], samples[1]);
    for(size_t k = 0; k < tunes.size(); k++) {
      chi += chi2_sim_test(tunes[k], samples[2]);
    }
  }
  chiTest->SetBinContent(chiTest->FindBin(param[i]),chi);
}
*/
chiTest->GetXaxis()->SetRangeUser(0.6,1.1);
chiTest->GetYaxis()->SetRangeUser(75,375);
chiTest->SetMarkerStyle(20);
chiTest->Draw("p9");
TFitResultPtr fit = chiTest->Fit("pol2","FS");
std:cout << report << std::endl;
float min = (-1)*fit->Parameter(1)/(2*fit->Parameter(2));
float chimin = fit->Parameter(0) + fit->Parameter(1)*min + fit->Parameter(2) * pow(min,2);
float err = (-1)*fit->Parameter(1) / (2 * fit->Parameter(2)) - sqrt(pow(fit->Parameter(1),2)
            - 4 * fit->Parameter(2) * (fit->Parameter(0) - chimin - 1)) / (2 * fit->Parameter(2));
//std::cout << "Minimum at x= " << (-1)*fit->Parameter(1)/(2*fit->Parameter(2)) << " +/- " << err << std::endl;
report = Form("Minimum at x= %0.3g +/- %0.2g",min, abs(min-err));
std::cout << report << std::endl;
std::cout << "chi^2_min + 1 at x= " << err << std::endl;

TPaveText *pt = new TPaveText(0.12,0.85,0.3,0.65,"NDC"); //NB blNDC
pt->SetFillStyle(0);
pt->SetTextAlign(11);
pt->SetBorderSize(0);
pt->SetTextFont(42);
pt->SetTextSize(0.046);
TString text = TString::Format("r_{B}= %.3f +/- %.2g",min,abs(min-err));
pt->AddText(text);
pt->Draw();
gStyle->SetOptStat(0);

}

float chi2_sim_test(TString tune="", TString sample="_d0_mu_tag_mu") {
std::cout << tune << " " << sample << std::endl;
TFile *fdata = TFile::Open(TString::Format("TopMass_Data_sPlot%s.root",sample.Data()));
TFile *fmc = TFile::Open(TString::Format("TopMass_172v5%s_sPlot%s.root",tune.Data(),sample.Data()));
//TFile *fmc = TFile::Open("TopMass_172v5_matched.root");

//load RooWorkspace and set binning
RooWorkspace *wmc = (RooWorkspace*)fmc->Get("w");
RooWorkspace *wdata = (RooWorkspace*)fdata->Get("w");
wmc->var("ptfrac")->setBins(22);
wdata->var("ptfrac")->setBins(22);
if(tune == "") ptfrac=*wmc->var("ptfrac");
RooDataSet *sigData;

//load MC into RooDataHist
TString cut("j_pt_ch<50");// && j_pt_ch<100");
sigData = (RooDataSet*)wmc->data("sigData")->reduce(cut);
//RooDataHist *ptfrac_mc_hist = new RooDataHist("ptfrac_mc_hist", "ptfrac_mc_hist", *wmc->var("ptfrac"), *wmc->data("sigData"));//*mcData);
RooDataHist *ptfrac_mc_hist = new RooDataHist("ptfrac_mc_hist", "ptfrac_mc_hist", *wmc->var("ptfrac"), *sigData);

sigData = (RooDataSet*)wdata->data("sigData")->reduce(cut);
//load Data into RooDataHist
//RooDataHist *ptfrac_data_hist = new RooDataHist("ptfrac_data_hist", "ptfrac_data_hist", *wdata->var("ptfrac"), *wdata->data("sigData"));//*dataData);
RooDataHist *ptfrac_data_hist = new RooDataHist("ptfrac_data_hist", "ptfrac_data_hist", *wdata->var("ptfrac"), *sigData);

//Chi2 fit
/*
TString title(TString::Format("ptfrac_chi2%s",tune.Data()));
RooChi2Var chi2_test(title,"#chi^{2} p_{T}", *ptfrac_mc_pdf, *ptfrac_data_hist);
RooMinuit m2_var(chi2_test);
m2_var.migrad();
m2_var.hesse();
RooFitResult *result_final = m2_var.save();
result_final->Print("v");
std::cout << chi2_test.getVal() << std::endl;
chi.push_back(chi2_test);

wmc->factory("Gaussian::g(x[-10,10],mean[5,0,1],sigma[0.01,0.001,0.1])");
RooGaussian g = *(RooGaussian*)wmc->pdf("g");
RooRealVar c("c", "c", 1, 0, 1.);
RooAddPdf model("model","c*ptfrac_mc_pdf + (1 - c)*g", *ptfrac_mc_pdf, g, c);
model.fitTo(*ptfrac_data_hist);
*/

RooPlot *frame = wdata->var("ptfrac")->frame();
ptfrac_data_hist->plotOn(frame);
//ptfrac_mc_hist->plotOn(frame,MarkerColor(38));
//model.plotOn(frame,LineColor(38));
frame->Draw();
//std::cout << frame->chiSquare("model","ptfrac_data_pdf") << std::endl;
/*
RooPlot *ptfrac_data = (RooPlot*)fdata->Get("ptfrac_mu_tag_signal");
RooPlot *ptfrac_mc = (RooPlot*)fmc->Get("ptfrac_mu_tag_signal");
*/
RooPlot *ptfrac_data = ptfrac_data_hist->plotOn(frame, RooFit::Binning(bins));
TH1F *data = (TH1F*)convert(ptfrac_data);
//data->Rebin(2);
frame = wdata->var("ptfrac")->frame();
RooPlot *ptfrac_mc = ptfrac_mc_hist->plotOn(frame, RooFit::Binning(bins));
TH1F *mc = (TH1F*)convert(ptfrac_mc);
/*
if(tune.Length() > 0) {
  TH1F *tuneWgt = (TH1F*)fmc->Get("tuneWgt");
  report += "ptfrac";
  report += tune;
  report += " scaled by ";
  report += tuneWgt->GetBinContent(1)/tuneWgt->GetBinContent(2);
  report += '\n';
  mc->Scale(tuneWgt->GetBinContent(1)/tuneWgt->GetBinContent(2));
}
*/
//mc->Rebin(2);
//TH1F *mc = (TH1F*)convert(((RooPlot*)fmc->Get("ptfrac_mu_tag_signal"))->getHist());
//TH1F *mc = (TH1F*)fmc->Get("ptfrac_signal_mu_tag_mu");
float chi2 = data->Chi2Test(mc, "CHI2");
std::cout << tune << " Chi2/ndf= " << chi2 << std::endl;

/*
delete ptfrac_data;
delete ptfrac_mc;
delete data;
delete mc;
fdata->Close();
fmc->Close();
*/
/*
c1->cd();
c1->SetFillStyle(400);
c1->SetBorderSize(0);
c1->SetLineWidth(0);
r->Draw();
gStyle->SetOptStat(0);
TLine *line = new TLine(0,1,1.2,1);
line->SetLineColor(kGreen+2);
line->SetLineStyle(10);
line->Draw();
TLine *line2 = new TLine(1,0,1,5);
line2->SetLineColor(kGreen+2);
line2->SetLineStyle(10);
line2->Draw();
c1->SaveAs("ratio_mu.png");
c1->SaveAs("ratio_mu.pdf");
*/
return chi2;

}
