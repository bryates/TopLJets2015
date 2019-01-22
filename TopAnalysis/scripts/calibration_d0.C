#include <RooFit.h>
#include "RooPlot.h"
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

#include "TFile.h"
#include "TH1F.h"
#include "TPaveText.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TAxis.h"
#include "TString.h"

using namespace RooFit;

TString report("");
float chi2_d0_test(TString, TString);
void control(bool full=false);
void buildVector(bool full=false);
std::vector<RooChi2Var> chi;
RooRealVar ptfrac;

//std::vector<TString> tune = {"", "_up", "_central", "_down"};
//std::vector<float> param = {0.855, 1.079, 0.8949, 0.6981};
/*
std::vector<TString> tune = {"_down", "", "_cccentral", "_central", "_up" };
std::vector<float> param = {0.755, 0.855, 0.875, 0.955, 1.055};
*/
/*
std::vector<TString> tune = {"_sdown", "_down", "_scentral", "_cccentral", "_central", "_up" };
std::vector<float> param = {0.655, 0.755, 0.825, 0.875, 0.955, 1.055};
std::vector<TString> tune = {"_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup" };
std::vector<float> param = {0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000};
std::vector<TString> tune = {"_700", "_725", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup" };
std::vector<float> param = {0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000};
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975, 1.000, 1.055};
<<<<<<< HEAD
*/
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055};
=======
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055};
*/
std::vector<TString> tune = {"_555", "_600", "_625", "_sdown", "_700", "_725", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.555, 0.600, 0.625, 0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975, 1.000, 1.055};
>>>>>>> fit
/*
std::vector<TString> tune = {"_sdown", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.655, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055};
std::vector<TString> tune = {"_down", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.755, 0.775, 0.825, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055};
*/
TCanvas *c1 = new TCanvas("c1","c1");
std::vector<TH1F*> chiTest;

void calibration_d0(TString in="", bool full=false) {
  if(in == "") {
    buildVector(full);
    control(full);
    std::cout << "not making calibration plots" << std::endl;
  }
  else {
    if(in == "nominal") in = TString("");
    for(auto & mc : tune) {
      if(in != mc) continue;
      chiTest.push_back(new TH1F("chiTest",Form("#chi^{2} test%s",mc.Data()),100,0,2));
      std::cout << "Closure test with " << mc << std::endl;
      TFile *fout = TFile::Open(TString::Format("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/calibration_chi2_d0%s.root",mc.Data()), "RECREATE");
      for(auto & it : tune) {
        if(mc==it) continue;
        //if(it == "_sdown" || it=="_uuup" || it == "_uup" || it == "_up") continue;
        std::cout << "Running on tune: " << it << std::endl;
        float chi = chi2_d0_test(mc, it);
        int pos = &it - &tune[0];
        chiTest.back()->SetBinContent(chiTest.back()->FindBin(param[pos]),chi);
      }
    
<<<<<<< HEAD
      chiTest.back()->GetXaxis()->SetRangeUser(0.6,1.1);
=======
      chiTest.back()->GetXaxis()->SetRangeUser(0.5,1.1);
>>>>>>> fit
      chiTest.back()->SetMarkerStyle(20);
      chiTest.back()->Draw("p9");
      //TFitResultPtr fit = chiTest.back()->Fit("pol2","FS");
      chiTest.back()->SetDirectory(fout);
      chiTest.back()->Write();
      fout->Write();
      fout->Close();
    }
    //control();
  }
}

void buildVector(bool full=false) {
  std::cout << "loading external plots" << std::endl;
  for(auto & mc : tune) {
    TFile *fin = TFile::Open(TString::Format("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/calibration_chi2_d0%s.root",mc.Data()));
    std::cout << TString::Format("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/calibration_chi2_d0%s.root",mc.Data()) << std::endl;
    if(!fin) continue;
    TH1F *tmp = (TH1F*)(fin->Get("chiTest"));
    if(!tmp) {
      std::cout << TString::Format("chiTest%d",(int)(param[&mc - &tune[0]]*1000)) << " not found!" << std::endl;
      break;
    }
    tmp->SetName(TString::Format("chiTest%d",(int)(param[&mc - &tune[0]]*1000)));
    tmp->SetTitle(TString::Format("chiTest%d",(int)(param[&mc - &tune[0]]*1000)));
    TH1F *htmp = (TH1F*)tmp->Clone();
    htmp->Reset();
    //remove points more that +/- 0.05 from test value
    //tmp->SetBinContent(tmp->FindBin(0.655), 0); //manually remove 0.655
    //tmp->SetBinError(tmp->FindBin(0.655), 0); //manually remove 0.655
    //tmp->SetBinContent(tmp->FindBin(1.055), 0); //manually remove 1.055
    //tmp->SetBinError(tmp->FindBin(1.055), 0); //manually remove 1.055
    /*
    for(size_t i = 0; i < param.size(); i++) {
      if(param[i] == param[&mc - &tune[0]]) continue;
      if(abs(param[i] - param[&mc - &tune[0]]) < 0.1) continue;
      int dist(4);
      if(fabs((&mc - &tune[0])-(float)i) < dist) continue; //up to dist away from point
      if((&mc - &tune[0])-dist > 0 || tune.size() > (&mc - &tune[0])-dist) continue; //check boundaries
      tmp->SetBinContent(tmp->FindBin(param[i]), 0); //delete leftover points
    }
    */
    for(size_t i = 0; i < param.size(); i++) {
      if(param[i] == param[&mc - &tune[0]]) continue;
      //if(param[i]==0.655) continue; //manually remove 0.655
      //if(param[i]==1.055) continue; //manually remove 1.055
      if(!full) {
        /*
        if(abs(param[i] - param[&mc - &tune[0]]) > 0.1) continue;
        */
        int dist(4);
        if(fabs((&mc - &tune[0])-(float)i) > dist) continue; //up to dist away from point
        //if((&mc - &tune[0] - i)-dist < 0 || tune.size() > (&mc - &tune[0] - i)-dist) continue; //check boundaries
      }
      htmp->SetBinContent(htmp->FindBin(param[i]), tmp->GetBinContent(tmp->FindBin(param[i])));
      //htmp->SetBinError(htmp->FindBin(param[i]), tmp->GetBinError(tmp->FindBin(param[i])));
    }
    //if(mc == "_up") continue;
    chiTest.push_back(htmp);
    std::cout << chiTest.back()->GetName() << std::endl;
  }
  std::cout << "loading external plots DONE!" << std::endl;
}

void control(bool full=false) {
TH1F *curve = new TH1F("curve","#chi^{2} curve", 400, -2, 2);
TH1F *curveg = new TH1F("curveg","#chi^{2} curve", 40, 0.6, 1.1);
for(size_t i = 0; i < tune.size(); i++) {
  if(!chiTest[i]) continue;
  std::cout << chiTest[i]->GetTitle() << std::endl;
<<<<<<< HEAD
  chiTest[i]->GetXaxis()->SetRangeUser(0.6,1.1);
=======
  chiTest[i]->GetXaxis()->SetRangeUser(0.5,1.1);
>>>>>>> fit
  //chiTest->GetYaxis()->SetRangeUser(0,5);
  chiTest[i]->SetMarkerStyle(20);
  chiTest[i]->Draw("p9");
  std::cout << "fitting with parabola" << std::endl;
  TFitResultPtr fit = chiTest[i]->Fit("pol2","FSMEQ");
  std:cout << report << std::endl;
  float min = (-1)*fit->Parameter(1)/(2*fit->Parameter(2));
  float chimin = fit->Parameter(0) + fit->Parameter(1)*min + fit->Parameter(2) * pow(min,2);
  std::cout << min << " " << chimin << std::endl;
  float err = (-1)*fit->Parameter(1) / (2 * fit->Parameter(2)) - sqrt(pow(fit->Parameter(1),2)
              - 4 * fit->Parameter(2) * (fit->Parameter(0) - chimin - 1)) / (2 * fit->Parameter(2));
  if(err != err) err = (-1)*fit->Parameter(1) / (2 * fit->Parameter(2)) - sqrt(pow(fit->Parameter(1),2)
                       - 4 * fit->Parameter(2) * (fit->Parameter(0) - chimin + 1)) / (2 * fit->Parameter(2));
  std::cout << err << std::endl;
  //report = Form("Minimum at x= %0.3g +/- %0.2g",min, abs(min-err));
  //std::cout << "Minimum at x= " << min << " +/- " << abs(min - err) << std::endl;
  //std::cout << report << std::endl;
  std::cout << Form("Minimum at x= %0.3g +/- %0.2g",min, abs(min-err)) << std::endl;
  std::cout << "chi^2_min + 1 at x= " << err << std::endl;
  
  TPaveText *pt = new TPaveText(0.12,0.85,0.3,0.65,"NDC"); //NB blNDC
  pt->SetFillStyle(0);
  pt->SetTextAlign(11);
  pt->SetBorderSize(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.046);
  TString text = TString(Form("r_{B}^{GEN}= %.3f", param[i]));
  pt->AddText(text);
  text = TString::Format("r_{B}^{fit}= %.3f +/- %.2g",min,abs(min-err));
  pt->AddText(text);
  pt->Draw();
  gStyle->SetOptStat(0);
  c1->SaveAs(TString::Format("chi2_d0-calibration%s.png",tune[i].Data()));
  c1->SaveAs(TString::Format("chi2_d0-calibration%s.pdf",tune[i].Data()));

  if(err != err) continue; //check for NaN
  if(fit->Parameter(2) < 0) continue; //ignore maxima
  std::cout << "filling calibration" << std::endl;
  curve->SetBinContent(curve->FindBin(param[i]-0.855), min);
  curve->SetBinError(curve->FindBin(param[i]-0.855), abs(min-err));
  curveg->SetBinContent(curveg->FindBin(param[i]), min);
  curveg->SetBinError(curveg->FindBin(param[i]), abs(min-err));

}
c1->cd();
//curve->GetXaxis()->SetRangeUser(0.7,1.1);
curve->GetXaxis()->SetRangeUser(-0.1,0.2);
<<<<<<< HEAD
curve->GetYaxis()->SetRangeUser(0.7,1.2);
=======
curve->GetYaxis()->SetRangeUser(0.5,1.05);
>>>>>>> fit
curve->Draw();
TFitResultPtr f = curve->Fit("pol1","FSMEQ");
float avg_deltag(0.);
for(auto & it : param) {
  float itg(it-0.855);
  float calib = f->Parameter(0) + f->Parameter(1)*itg;// + f->Parameter(2)*itg*itg;
  float deltag = calib - it;
  avg_deltag += deltag;
}
avg_deltag /= param.size();
std::cout << "Average difference = " << avg_deltag << std::endl;
std::cout << "\tGEN\tfit\tfit-avg_diff" << std::endl;
for(auto & it : param) {
  float itg(it-0.855);
  float calib = f->Parameter(0) + f->Parameter(1)*itg;// + f->Parameter(2)*itg*itg;
  //std::cout << "r_B: " << it << " = " << calib << " -> " << calib - avg_deltag << " +/- " << sqrt(pow(f->ParError(0),2)+pow(f->ParError(1),2)) << std::endl;
  std::cout << TString::Format("r_B: %.3f  = %.3f -> %.3f +/- %.3f",it,calib,calib-avg_deltag,sqrt(pow(f->ParError(0),2)+pow(f->ParError(1),2))) << std::endl;
  //std::cout << "r_B: " << it << " = " << calib << " -> " << calib - avg_deltag << " +/- " << sqrt(pow(f->ParError(0),2)+pow(f->ParError(1),2)+pow(f->ParError(2),2)) << std::endl;
}
char sign = '+';
if(f->Parameter(1)<0) sign = '-';
char sign2 = '+';
if(f->Parameter(2)<0) sign2 = '-';
float yint = f->Parameter(0) - 0.855 * f->Parameter(1);
//TString text = TString::Format("Calibration curve : %0.2f (#pm%0.2f) %c r_{B} %0.2f (#pm%0.2f)",f->Parameter(0),abs(f->ParError(0)),sign,f->Parameter(1),abs(f->ParError(1)));
TString text = TString::Format("Calibration curve : %0.2f (#pm%0.2f) %c r_{B} %0.2f (#pm%0.2f)",yint,abs(f->ParError(0)),sign,f->Parameter(1),abs(f->ParError(1)));
curveg->Draw();
<<<<<<< HEAD
curveg->GetXaxis()->SetRangeUser(0.6,1.05);
curveg->GetYaxis()->SetRangeUser(0.6,1.05);
if(full) {
  curveg->GetXaxis()->SetRangeUser(0.6,1.05);
  curveg->GetYaxis()->SetRangeUser(0.6,1.45);
=======
curveg->GetXaxis()->SetRangeUser(0.5,1.05);
curveg->GetYaxis()->SetRangeUser(0.5,1.05);
if(full) {
  curveg->GetXaxis()->SetRangeUser(0.5,1.05);
  curveg->GetYaxis()->SetRangeUser(0.5,1.45);
>>>>>>> fit
}
TF1 *func = new TF1("pol1","pol1",0.5,1.1);
func->SetParameter(0, yint);
func->SetParameter(1, f->Parameter(1));
//curveg->Fit("pol1", "FS");
func->Draw("same");
TPaveText *pt = new TPaveText(0.12,0.85,0.3,0.65,"NDC"); //NB blNDC
pt->SetFillStyle(0);
pt->SetTextAlign(11);
pt->SetBorderSize(0);
pt->SetTextFont(42);
pt->SetTextSize(0.046);
//TString text = TString::Format("Calibration curve : %0.2f (#pm%0.2f) %c r_{B} %0.2f (#pm%0.2f) %c r_{B}^{2} %0.2f (#pm%0.2f)",f->Parameter(0),abs(f->ParError(0)),sign,f->Parameter(1),abs(f->ParError(1)),sign2,f->Parameter(2),abs(f->ParError(2)));
pt->AddText(text);
pt->SetFillColor(0);
pt->Draw();
c1->SaveAs("chi2_d0-calibration_curve.png");
c1->SaveAs("chi2_d0-calibration_curve.pdf");
}

float chi2_d0_test(TString mcTune="", TString tune="") {
<<<<<<< HEAD
TFile *fdata = TFile::Open(TString::Format("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/4dau/TopMass_172v5%s_sPlot_d0.root",mcTune.Data()));
TFile *fmc = TFile::Open(TString::Format("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/4dau/TopMass_172v5%s_sPlot_d0.root",tune.Data()));
=======
TFile *fdata = TFile::Open(TString::Format("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/TopMass_172v5%s_sPlot_d0.root",mcTune.Data()));
TFile *fmc = TFile::Open(TString::Format("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/TopMass_172v5%s_sPlot_d0.root",tune.Data()));
>>>>>>> fit

//load RooWorkspace and set binning
RooWorkspace *wmc = (RooWorkspace*)fmc->Get("w");
RooWorkspace *wdata = (RooWorkspace*)fdata->Get("w");
wmc->var("ptfrac")->setBins(22);
wdata->var("ptfrac")->setBins(22);
if(tune == "") ptfrac=*wmc->var("ptfrac");

//load MC into RooDataHist
RooDataSet *sigData = (RooDataSet*)wmc->data("sigData")->reduce("j_pt_ch<75");
RooDataHist *ptfrac_mc_hist = new RooDataHist("ptfrac_hist", "ptfrac_hist", *wmc->var("ptfrac"), *sigData);//*mcData);
RooHistPdf *ptfrac_mc_pdf = new RooHistPdf("ptfrac_mc_pdf", "ptfrac_mc_pdf", RooArgList(*wmc->var("ptfrac")), *ptfrac_mc_hist);

//load Data into RooDataHist
sigData = (RooDataSet*)wdata->data("sigData")->reduce("j_pt_ch<75");
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

RooPlot *frame = wdata->var("ptfrac")->frame();
ptfrac_data_hist->plotOn(frame);
//ptfrac_mc_hist->plotOn(frame,MarkerColor(38));
//model.plotOn(frame,LineColor(38));
frame->Draw();
//std::cout << frame->chiSquare("model","ptfrac_data_pdf") << std::endl;
RooPlot *ptfrac_data = (RooPlot*)fdata->Get("ptfrac_signal");
RooPlot *ptfrac_mc = (RooPlot*)fmc->Get("ptfrac_signal");
TH1F *data = (TH1F*)convert(ptfrac_data);
//data->Rebin(2);
TH1F *mc = (TH1F*)convert(ptfrac_mc);
if(tune.Length() > 0) {
  TH1F *tuneWgt = (TH1F*)fmc->Get("tuneWgt");
  report += "ptfrac";
  report += tune;
  report += " scaled by ";
  report += tuneWgt->GetBinContent(1)/tuneWgt->GetBinContent(2);
  report += '\n';
  mc->Scale(tuneWgt->GetBinContent(1)/tuneWgt->GetBinContent(2));
}
float chi2 = data->Chi2Test(mc, "CHI2");
std::cout << tune << " Chi2= " << chi2 << std::endl;

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
