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
float low(999.), high(0.),nom(0.818905),nerr(0.05);
bool tdr(0);

TString report("");
TString json("\"d0_mu\" :      [");
float chi2_d0_mu_tag_test(TString tune="");
void run_chi2_d0_mu_tag(TString);
std::vector<RooChi2Var> chi;
RooRealVar ptfrac;

void chi2_d0_mu_tag() {
  run_chi2_d0_mu_tag("");
  run_chi2_d0_mu_tag("fsr-down");
  run_chi2_d0_mu_tag("fsr-up");
  run_chi2_d0_mu_tag("uedown");
  run_chi2_d0_mu_tag("ueup");
  run_chi2_d0_mu_tag("erdON");
  std::vector<TString> syst = {"LEP", "TRIGGER", "TRK", "PU", "PI", "JER" };
  for(auto & it : syst) {
    run_chi2_d0_mu_tag("down_"+it);
    run_chi2_d0_mu_tag("up_"+it);
  }

  json += ("],");
  std::cout << json << std::endl;

}

void run_chi2_d0_mu_tag(TString lname="") {
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
TH1F *chiTest = new TH1F("chiTest","#chi^{2} test",400,0,2);
chiTest->SetDirectory(0);
for(auto & it : tune) {
  int pos = &it - &tune[0];
  //if(param[pos]>1) continue;
  std::cout << "Running on tune: " << it << std::endl;
  float chi = chi2_d0_mu_tag_test(it);
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
chiTest->GetXaxis()->SetRangeUser(0.65,0.975);
//chiTest->GetYaxis()->SetRangeUser(55,90);
chiTest->GetYaxis()->SetRangeUser(int(low)-1,int(high)+2);
//chiTest->GetYaxis()->SetRangeUser(200,220);
chiTest->SetMarkerStyle(20);
chiTest->Draw("p9");
TFitResultPtr fit = chiTest->Fit("pol3","FSEMQ","",0.6,1.055);
//TFitResultPtr fit = chiTest->Fit("pol3","FSEMQ","",0.6,0.975);//1.055);
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
float err = chiTest->GetFunction("pol3")->GetX((chimin+1),0.7,1.0);
//float err = chiTest->GetFunction("pol3")->GetX((chimin+1),min,min+1);
float errl = chiTest->GetFunction("pol3")->GetX((chimin+1),min-1,min);
if(lname=="") { nom=min; nerr=err; }
report = Form("Minimum at x= %g +/- %0.6g",min, abs(min-err));
//report = Form("Minimum at x= %g + %0.6g - %0.6g",min, abs(min-err), abs(min-errl));
json += Form(" %.4f, %.4f,",min,abs(min-err));
//json += Form(" %.4f, %.4f, %.4f,",min,abs(min-err),abs(min-errl));
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
c1->SaveAs("chi2_d0_mu_tag"+name+".pdf");
c1->SaveAs("chi2_d0_mu_tag"+name+".png");

}

float chi2_d0_mu_tag_test(TString tune="") {
TFile *fdata = TFile::Open("TopMass_Data_sPlot_d0_mu_tag_mu.root");
TFile *fmc;
if(name.Length()==0)
fmc = TFile::Open(TString::Format("TopMass_172v5%s_sPlot_d0_mu_tag_mu.root",tune.Data()));
else
fmc = TFile::Open(TString::Format("TopMass_%s%s_sPlot_d0_mu_tag_mu.root",name.Data(),tune.Data()));
//TFile *fmc = TFile::Open(TString::Format("TopMass_ueup%s_sPlot_d0.root",tune.Data()));
//TFile *fmc = TFile::Open(TString::Format("TopMass_erdOn%s_sPlot_d0_mu_tag.root",tune.Data()));
//TFile *fmc = TFile::Open(TString::Format("TopMass_fsr-down%s_sPlot_d0_mu_tag.root",tune.Data()));
//TFile *fmc = TFile::Open(TString::Format("TopMass_ueup%s_sPlot_d0_mu_tag.root",tune.Data()));
//TFile *fmc = TFile::Open("TopMass_172v5_matched.root");

//load RooWorkspace and set binning
RooWorkspace *wmc = (RooWorkspace*)fmc->Get("w");
RooWorkspace *wdata = (RooWorkspace*)fdata->Get("w");
wmc->var("ptfrac")->setBins(22);
wdata->var("ptfrac")->setBins(22);
if(tune == "") ptfrac=*wmc->var("ptfrac");

TString cut("j_pt_ch<75");
RooDataSet *sigData = (RooDataSet*)wmc->data("sigData");//->reduce(cut);
//load MC into RooDataHist
RooDataHist *ptfrac_mc_hist = new RooDataHist("ptfrac_hist", "ptfrac_hist", *wmc->var("ptfrac"), *sigData);//*mcData);
RooHistPdf *ptfrac_mc_pdf = new RooHistPdf("ptfrac_mc_pdf", "ptfrac_mc_pdf", RooArgList(*wmc->var("ptfrac")), *ptfrac_mc_hist);

sigData = (RooDataSet*)wdata->data("sigData");//->reduce(cut);
//load Data into RooDataHist
RooDataHist *ptfrac_data_hist = new RooDataHist("ptfrac_hist", "ptfrac_hist", *wdata->var("ptfrac"), *sigData);//*dataData);
RooHistPdf *ptfrac_data_pdf = new RooHistPdf("ptfrac_data_pdf", "ptfrac_data_pdf", RooArgList(*wdata->var("ptfrac")), *ptfrac_data_hist);

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

RooPlot *ptfrac_data = wdata->var("ptfrac")->frame();
ptfrac_data_hist->plotOn(ptfrac_data);
RooPlot *ptfrac_mc = wdata->var("ptfrac")->frame();
ptfrac_mc_hist->plotOn(ptfrac_mc);
TH1F *data = (TH1F*)convert(ptfrac_data,true,0,1.1);
//data->Scale(1./data->Integral());
TH1F *mc = (TH1F*)convert(ptfrac_mc,true,0,1.1);
//mc->Scale(1./mc->Integral());
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
