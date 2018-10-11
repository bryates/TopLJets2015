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
TString name("");
float low(50.), high(50.),nom(0.818905),nerr(0.05);
bool tdr(1);

TString report("");
TString json("\"d0\" :      [");
float chi2_d0_test(TString tune="");
void run_chi2_d0(TString);
std::vector<RooChi2Var> chi;
RooRealVar ptfrac;

void chi2_d0() {
  run_chi2_d0("");
  /*
  run_chi2_d0("fsr-up");
  run_chi2_d0("fsr-down");
  run_chi2_d0("ueup");
  run_chi2_d0("uedown");
  run_chi2_d0("erdON");
  std::vector<TString> syst = {"LEP", "TRIGGER", "TRK", "PU", "PI", "JER" };
  for(auto & it : syst) {
    run_chi2_d0("up_"+it);
    run_chi2_d0("down_"+it);
  }
  */

  json += ("],");
  std::cout << json << std::endl;

}

void run_chi2_d0(TString lname="") {
low=999.;
high=0;
name=lname;
//gROOT->ProcessLine(".L convert.C");
//std::vector<TString> tune = {"", "_up", "_central", "_down"};
//std::vector<float> param = {0.855, 1.079, 0.8949, 0.6981};
/*
std::vector<TString> tune = {"_down", "", "_cccentral", "_central", "_up" };
std::vector<float> param = {0.755, 0.855, 0.875, 0.955, 1.055};
std::vector<TString> tune = {"_sdown", "_down", "_scentral", "", "_cccentral", "_central", "_up" };
std::vector<float> param = {0.655, 0.755, 0.825, 0.855, 0.875, 0.955, 1.055};
std::vector<TString> tune = {"_sdown", "_down", "_scentral", "", "_cccentral", "_925", "_central", "_up" };
std::vector<float> param = {0.655, 0.755, 0.825, 0.855, 0.875, 0.925, 0.955, 1.055};
*/
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_up" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975, 1.055, 0.802};
/*
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975, 1.000, 1.055};
std::vector<TString> tune = {"_down", "_ddown", "_dddown", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.755, 0.775, 0.800, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055};
*/
TCanvas *c1 = new TCanvas("c1","c1");
TH1F *chiTest = new TH1F("chiTest","",400,0,2);
chiTest->SetDirectory(0);
for(auto & it : tune) {
  int pos = &it - &tune[0];
  if(param[pos]>1) continue;
  std::cout << "Running on tune: " << it << std::endl;
  float chi = chi2_d0_test(it);
  chiTest->SetBinContent(chiTest->FindBin(param[pos]),chi);
}
/*
chi2("");
*/

/*
RooArgSet args;
for(auto &it : chi)
  args.add(it);
RooAddition sumchi2("sumchi2", "#Sigma #chi^{2}", args);

RooMinuit c2_var(sumchi2);
c2_var.migrad();
c2_var.hesse();
RooFitResult *result_final = c2_var.save();
result_final->Print("v");

RooRealVar alpha("alpha", "alpha", 0, 1.0);
RooIntegralMorph lmorph("lmorph", "lmorph", hists[0], hists[1], ptfrac, alpha);
RooPlot *frame = ptfrac.frame();
hists[0].plotOn(frame);
hists[1].plotOn(frame,LineColor(kRed+2));
alpha.setVal(0);
lmorph.plotOn(frame);
alpha.setVal(1);
lmorph.plotOn(frame,LineColor(kRed));
frame->Draw();
*/


//chiTest->GetXaxis()->SetRangeUser(0.65,1.055);
chiTest->GetXaxis()->SetRangeUser(0.65,0.976);//1.055);
//chiTest->GetYaxis()->SetRangeUser(55,90);
chiTest->GetYaxis()->SetRangeUser(int(low)-1,int(high)+2);
//chiTest->GetYaxis()->SetRangeUser(200,220);
chiTest->SetMarkerStyle(20);
chiTest->Draw("p9");
TLatex txt;
txt.SetNDC(true);
txt.SetTextFont(43);
txt.SetTextSize(16);
txt.SetTextAlign(12);
float iniy=0.95;// if self.wideCanvas else 0.95
float inix=0.12;// if noStack else 0.12
float lumi(35859.038);
if(lumi<100)
    txt.DrawLatex(inix,iniy,TString::Format("#bf{CMS} #it{Preliminary} %3.1f pb^{-1} (13 TeV)", (lumi) ));
else
    txt.DrawLatex(inix,iniy,TString::Format("#bf{CMS} #it{Preliminary} %3.1f fb^{-1} (13 TeV)", (lumi/1000.) ));
//TFitResultPtr fit = chiTest->Fit("pol3","FSEMQ","",0.6,1.055);
TFitResultPtr fit = chiTest->Fit("pol3","FSEMQ","",0.6,0.975);
//TFitResultPtr fit = chiTest->Fit("pol2","FSMEQ");
//TFitResultPtr fit = chiTest->Fit("pol2","FSMEQ","",0.8,1.0);
std:cout << report << std::endl;
/*
float min = (-1)*fit->Parameter(1)/(2*fit->Parameter(2));
float chimin = fit->Parameter(0) + fit->Parameter(1)*min + fit->Parameter(2) * pow(min,2);
float err = (-1)*fit->Parameter(1) / (2 * fit->Parameter(2)) - sqrt(pow(fit->Parameter(1),2)
            - 4 * fit->Parameter(2) * (fit->Parameter(0) - chimin - 1)) / (2 * fit->Parameter(2));
*/
float min = chiTest->GetFunction("pol3")->GetMinimumX(0.7,1.0);
float chimin = fit->Parameter(0) + fit->Parameter(1)*min + fit->Parameter(2) * pow(min,2) + fit->Parameter(3) * pow(min,3);
float err = chiTest->GetFunction("pol3")->GetX(chimin+1,0.7,1.0);
if(lname=="") { nom=min; nerr=err; }
report = Form("Minimum at x= %g +/- %0.6g",min, abs(min-err));
json += Form(" %.4f, %.4f,",min,abs(min-err));
//std::cout << "Minimum at x= " << min << " +/- " << abs(min - err) << std::endl;
std::cout << report << std::endl;
std::cout << "chi^2_min + 1 at x= " << err << std::endl;

TPaveText *pt = new TPaveText(0.12,0.85,0.3,0.65,"NDC"); //NB blNDC
pt->SetFillStyle(0);
pt->SetTextAlign(11);
pt->SetBorderSize(0);
pt->SetTextFont(42);
pt->SetTextSize(0.046);
TString text = TString::Format("r_{B}= %.4f +/- %.4f (stat)",min,abs(min-err));
if(name.Length() > 0)
  text += TString::Format(" %c %.4f (syst) +/- %.4f",(min<nom ? '-' : '+'), abs(nom-min), sqrt(abs(pow(nerr,2)-pow(abs(min-err),2))));
  //text += TString::Format(" %c %.4f (syst) +/- %.4f",(min<nom ? '-' : '+'), abs(nom-min), sqrt(abs(pow(0.0507584,2)-pow(abs(min-err),2))));
  //text += TString::Format(" %c %.4f (syst) +/- %.4f",(min<0.818905 ? '-' : '+'), abs(0.818905-min), sqrt(abs(pow(0.0507584,2)-pow(abs(min-err),2))));
pt->AddText(text);
if(!tdr) pt->Draw();
gStyle->SetOptStat(0);

if(name.Length()>0) name = "_" + name;
c1->SaveAs("chi2_d0"+name+".pdf");
c1->SaveAs("chi2_d0"+name+".png");

}

float chi2_d0_test(TString tune="") {
TFile *fdata = TFile::Open("sPlot/sPlot/TopMass_Data_sPlot_d0.root");//,"UPDATE");
TFile *fmc;
if(name.Length()==0)
fmc = TFile::Open(TString::Format("sPlot/sPlot/TopMass_172v5%s_sPlot_d0.root",tune.Data()),"UPDATE");
else
fmc = TFile::Open(TString::Format("sPlot/sPlot/TopMass_%s%s_sPlot_d0.root",name.Data(),tune.Data()),"UPDATE");
//TFile *fmc = TFile::Open(TString::Format("TopMass_ueup%s_sPlot_d0.root",tune.Data()));
//TFile *fmc = TFile::Open(TString::Format("TopMass_erdOn%s_sPlot_d0.root",tune.Data()));
//TFile *fmc = TFile::Open(TString::Format("TopMass_fsr-down%s_sPlot_d0.root",tune.Data()));
//TFile *fmc = TFile::Open(TString::Format("TopMass_ueup%s_sPlot_d0.root",tune.Data()));
//TFile *fmc = TFile::Open("TopMass_172v5_matched.root");

//load RooWorkspace and set binning
TString cut("j_pt_ch<75");

TH1F *mc,*data;
if(fmc->GetListOfKeys()->Contains("h_ptfrac_hist")) mc = (TH1F*)fmc->Get("h_ptfrac_hist");
else {
RooWorkspace *wmc = (RooWorkspace*)fmc->Get("w");
if(tune == "") ptfrac=*wmc->var("ptfrac");
wmc->var("ptfrac")->setBins(22);
RooDataSet *sigData = (RooDataSet*)wmc->data("sigData")->reduce(cut);
//load MC into RooDataHist
RooDataHist *ptfrac_mc_hist = new RooDataHist("ptfrac_hist", "ptfrac_hist", *wmc->var("ptfrac"), *sigData);//*mcData);
RooHistPdf *ptfrac_mc_pdf = new RooHistPdf("ptfrac_mc_pdf", "ptfrac_mc_pdf", RooArgList(*wmc->var("ptfrac")), *ptfrac_mc_hist);
RooPlot *ptfrac_mc = wmc->var("ptfrac")->frame();
ptfrac_mc_hist->plotOn(ptfrac_mc);
mc = (TH1F*)convert(ptfrac_mc,true,0,1.1);
mc->SetDirectory(fmc);
mc->SetTitle("ptfrac_sig");
mc->Write();
delete ptfrac_mc;
}

if(fdata->GetListOfKeys()->Contains("h_ptfrac_hist")) data = (TH1F*)fdata->Get("h_ptfrac_hist");
else {
RooWorkspace *wdata = (RooWorkspace*)fdata->Get("w");
if(tune == "") ptfrac=*wdata->var("ptfrac");
wdata->var("ptfrac")->setBins(22);
RooDataSet *sigData = (RooDataSet*)wdata->data("sigData")->reduce(cut);
//load Data into RooDataHist
RooDataHist *ptfrac_data_hist = new RooDataHist("ptfrac_hist", "ptfrac_hist", *wdata->var("ptfrac"), *sigData);//*dataData);
RooHistPdf *ptfrac_data_pdf = new RooHistPdf("ptfrac_data_pdf", "ptfrac_data_pdf", RooArgList(*wdata->var("ptfrac")), *ptfrac_data_hist);
RooPlot *ptfrac_data = wdata->var("ptfrac")->frame();
ptfrac_data_hist->plotOn(ptfrac_data);
data = (TH1F*)convert(ptfrac_data,true,0,1.1);
data->SetDirectory(fdata);
data->SetTitle("ptfrac_sig");
data->Write();
delete ptfrac_data;
delete sigData;
}

if(tune.Length() > 0) {
  TH1F *tuneWgt = (TH1F*)fmc->Get("tuneWgt");
  report += "ptfrac";
  report += tune;
  report += " scaled by ";
  report += tuneWgt->GetBinContent(1)/tuneWgt->GetBinContent(2);
  report += '\n';
  //mc->Scale(tuneWgt->GetBinContent(1)/tuneWgt->GetBinContent(2));
}
float chi2 = data->Chi2Test(mc, "CHI2 WW");
std::cout << tune << " Chi2/ndf= " << chi2 << std::endl;
if(chi2<low) low = chi2;
if(chi2>high) high = chi2;

delete data;
delete mc;
fdata->Close();
fmc->Close();
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
