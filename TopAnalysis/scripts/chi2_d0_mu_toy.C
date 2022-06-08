#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <RooFit.h>
#include "RooGlobalFunc.h"
#include "RooWorkspace.h"
#include "RooHistPdf.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "RooBinning.h"
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
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/tdr.C"

using namespace RooFit;
//TString name("");
float low(50.), high(50.),nom(0.8103),nerr(0.05);
bool fin(true);
bool TDR(1);
int epoch(0);
bool fullpt(0);
TString epoch_name[3] = {"_BCDEFGH", "_BCDEF", "_GH"};

TCanvas *c1 = setupCanvas();
TString report("");
TString json("\"d0_mu\" :   [");
float chi2_d0_mu_tag_test(TH1F *&shiftData, TString tune="", TString name="", float num=0.855);
void run_chi2_d0_mu_tag(TString);
RooRealVar ptfrac;
std::vector<TH1F*> mcHists;
TH1F *dataHist;

void getHist(TString name, TString tune, TH1F *&data, TH1F *&mc, int epoch, bool norm=true) {
TString fname = TString::Format("sPlot/sPlot/TopMass_Data_sPlot_d0_mu_tag_mu.root");
/*
if(name.Contains("_toys") && !name.Contains("FSR"))
  fname = TString::Format("sPlot/sPlot/TopMass_%s%s_sPlot_d0_mu_tag_mu.root",name.Data(),tune.Data());
*/
if(name.Contains("_toys")) fname.ReplaceAll("_toys","");
if(name.Contains("_toyData")) fname.ReplaceAll("_toyData","");
if(epoch>0) fname.ReplaceAll(".root",TString::Format("%d.root",epoch));
else if(epoch<0) fname.ReplaceAll(".root","_xb.root");
//if(epoch<0)fname.ReplaceAll("sPlot/sPlot/","");
if(fullpt) fname.ReplaceAll(".root","_jpT.root");
if(name.Contains("Bfrag_genreco")) fname.ReplaceAll("Data","Data_Bfrag_genreco");
if(name.Contains("FSR") && name.Contains("_toyData")) fname.ReplaceAll("Data","FSR_toyData");
else if(name.Contains("FSR")) fname.ReplaceAll("Data","FSR");
std::cout << fname << std::endl;
TFile *fdata = TFile::Open(fname);
if(name.Length()==0)
fname = TString::Format("sPlot/sPlot/TopMass_172v5%s_sPlot_d0_mu_tag_mu.root",tune.Data());
//fname = TString::Format("sPlot/sPlot/TopMass_172v5%s_toppTOFF_sPlot_d0_mu_tag_mu.root",tune.Data());
else
fname = TString::Format("sPlot/sPlot/TopMass_%s%s_sPlot_d0_mu_tag_mu.root",name.Data(),tune.Data());
if(name.Contains("_toys")) fname.ReplaceAll("_toys","");
if(name.Contains("_toyData")) fname.ReplaceAll("_toyData","");
if(epoch>0) fname.ReplaceAll(".root",TString::Format("%d.root",epoch));
else if(epoch>0) fname.ReplaceAll(".root","_inc.root");
if(epoch<0)fname.ReplaceAll("sPlot/sPlot/","");
if(fullpt) fname.ReplaceAll(".root","_jpT.root");
std::cout << fname << std::endl;
TFile *fmc = TFile::Open(fname);

RooBinning bins(0,1.1);
/*
std::vector<float> bin;
if(!fullpt) bin = {0, 0.2, 0.4, 0.6, 0.7, 0.75, 0.8, 0.82, 0.84,0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0};
bin = {0, 0.2, 0.4, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
//float bin[] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.1};
if(fullpt) bin = {0, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 1.0};
  bin = {0, 0.2, 0.4, 0.6, 0.7, 0.75, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0};
 bin = {0, 0.2, 0.4, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0}; 
*/
std::vector<float> bin = {0, 0.2, 0.4, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0}; 
for(int i = 0; i < bin.size(); i++) {
 bins.addBoundary(bin[i]);
}

RooPlot *tmp = nullptr;
if(epoch<0) mc = (TH1F*)fmc->Get("ptfrac_signal_smooth")->Clone();
//if(epoch<0) mc = (TH1F*)fmc->Get("ptfrac_signal_hist")->Clone();
else {
tmp = (RooPlot*)fmc->Get("ptfrac_mu_tag_signal")->Clone(TString::Format("ptfrac_signal_mc%s%s",name.Data(),tune.Data()));
if(tmp==nullptr) {std::cout << fname << std::endl; return;}
//delete tmp;
//float bin[] = {0, 0.2, 0.6, 0.7, 0.75, 0.8, 0.82, 0.84,0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0};
//float *bin;
delete tmp;
tmp = ((RooWorkspace*)fmc->Get("w"))->var("ptfrac")->frame();
//((RooDataSet*)((RooWorkspace*)fmc->Get("w"))->data("sigData"))->reduce("ptfrac<0.99")->plotOn(tmp, RooFit::Binning(bins));
((RooDataSet*)((RooWorkspace*)fmc->Get("w"))->data("sigData"))->plotOn(tmp, RooFit::Binning(bins));
//else ((RooDataSet*)((RooWorkspace*)fmc->Get("w"))->data("sigData"))->reduce("j_pt<200")->plotOn(tmp, RooFit::Binning(22));
//mc = (TH1F*)convert(tmp, norm, 0, 1.1);
mc = (TH1F*)convert(tmp, norm, bin);
}
//mc = (TH1F*)convert(tmp, norm, 13, bin);
mc->SetDirectory(0);
mc->SetTitle(mc->GetName());
delete tmp;
if(epoch<0) data = (TH1F*)fdata->Get("ptfrac_signal_hist")->Clone();
else {
tmp = (RooPlot*)fdata->Get("ptfrac_mu_tag_signal")->Clone(TString::Format("ptfrac_signal_data%s%s",name.Data(),tune.Data()));
delete tmp;
tmp = ((RooWorkspace*)fdata->Get("w"))->var("ptfrac")->frame();
((RooDataSet*)((RooWorkspace*)fdata->Get("w"))->data("sigData"))->plotOn(tmp, RooFit::Binning(bins));
//else ((RooDataSet*)((RooWorkspace*)fdata->Get("w"))->data("sigData"))->reduce("j_pt<200")->plotOn(tmp, RooFit::Binning(22));
//data = (TH1F*)convert(tmp, norm, 0, 1.1);
data = (TH1F*)convert(tmp, norm, bin);
}
//data = (TH1F*)convert(tmp, norm, 8, bin);
data->SetDirectory(0);
data->SetTitle(data->GetName());
delete tmp;

fdata->Close();
fmc->Close();
delete fdata;
delete fmc;
}

void chi2_d0_mu_toy(TString samp="") {
    if(samp != "")
      run_chi2_d0_mu_tag(samp);
    else {
      run_chi2_d0_mu_tag("");
      run_chi2_d0_mu_tag("isr-down");
      run_chi2_d0_mu_tag("isr-up");
      run_chi2_d0_mu_tag("fsr-down");
      run_chi2_d0_mu_tag("fsr-up");
      run_chi2_d0_mu_tag("uedown");
      run_chi2_d0_mu_tag("ueup");
      run_chi2_d0_mu_tag("erdON");
      //run_chi2_d0_mu_tag("GluonMove_erdON");
      //run_chi2_d0_mu_tag("QCD_erdON");
      std::vector<TString> syst = {"LEP", "PU", "PI", "TRIGGER", "JER" }; //no lepton tracker efficiencies are used!
      //std::vector<TString> syst = {"TRK", "LEP", "PU", "PI", "TRIGGER", "JER" };
      for(auto & it : syst) {
        run_chi2_d0_mu_tag("down_"+it);
        run_chi2_d0_mu_tag("up_"+it);
      }
      run_chi2_d0_mu_tag("hdampdown");
      run_chi2_d0_mu_tag("hdampup");
    }
    for(auto mc : mcHists)
      std::cout << mc->Integral() << std::endl;
  std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_up"};//, "_1075", "_1125" };
  std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975, 1.055};//, 1.075, 1.125, 0.802};
  TString json("\"d0_mu\" :   [");
  for(int i = 0; i < 1000; i++) {
    std::cout << "iteration: " << i << std::endl;
    if(!samp.Contains("_toy") && i > 0) continue;
    TRandom3 *rand = new TRandom3(0);
    dataHist = (TH1F*)mcHists[0]->Clone();
    std::cout << "Shifting toy data " << dataHist->Integral() << std::endl;
    for(int i = 1; i <= dataHist->GetNbinsX(); i++) {
      float y = dataHist->GetBinContent(i);
      float e(0.);
      e = dataHist->GetBinError(i);
      std::cout << y << " +/- " << e << " ";
      //shift each bin by random amount samples by the appropriate bin error
      if(!(samp == "" || samp == "172v5")) //MC errors otherwise (e.g. FSR)
        y += rand->Gaus(0, e);
      dataHist->SetBinContent(i,y);
    }
    delete rand;
    TH1F *chiTest = new TH1F("chiTest_"+samp,TString::Format("chiTest_%s",samp.Data()),400,0,2);
    chiTest->Sumw2();
    chiTest->GetXaxis()->SetRangeUser(0.65,0.976);//1.055);
    std::cout << "Fitting to toy data" << std::endl;
    low=999.;
    high=0;
    for(int pos = 1; pos < mcHists.size(); pos++) {
      TH1F *mc = mcHists[pos];
      if(pos==0) continue;
      //std::cout << "Fitting " << tune[pos-1] << std::endl;
      float chi2 = dataHist->Chi2Test(mc, "CHI2 P WW");
      std::cout << chi2 << std::endl;
      if(chi2<low) low = chi2;
      if(chi2>high) high = chi2;
      chiTest->GetYaxis()->SetRangeUser(max(int(low)-1,0),int(high)+2);
      chiTest->SetBinContent(chiTest->FindBin(param[pos-1]),chi2);
    }
    chiTest->SetMarkerStyle(20);
    chiTest->Draw("p9");
    std::cout << "Fitting pol3" << std::endl;
    ((TF1*)(gROOT->GetFunction("pol3")))->SetParameters(1., 1., 1., 1.);
    chiTest->Fit("pol3","FSMEQRW","",0.6,0.976);
    c1->SaveAs(TString::Format("chi2_d0_mu_tag_%s_%d.pdf", samp.Data(), i));
    c1->SaveAs(TString::Format("chi2_d0_mu_tag_%s_%d.png", samp.Data(), i));
    if(samp.Contains("fsr-down")) chiTest->Fit("pol3","FSMEQRW","",0.6,1.126);//0.976);
    std::cout << "Getting pol3 min" << std::endl;
    float min = chiTest->GetFunction("pol3")->GetMinimumX(0.6,1.126);
    chiTest->GetFunction("pol3")->Print("v");
    std::cout << "Getting pol3 err" << std::endl;
    float chimin = chiTest->GetFunction("pol3")->Eval(min);
    float err = chiTest->GetFunction("pol3")->GetX(chimin+1,0.6,1.126);
    if(samp=="") { nom=min; nerr=err; }
    report = Form("Minimum at x= %0.3f +/- %0.3f",min, abs(min-err));
    std::cout << i << ":" << report << std::endl;
    json += Form("%.4f, %.4f, ",min,abs(min-err));
/*  
    run_chi2_d0_mu_tag("tpt");
    run_chi2_d0_mu_tag("bkg");
*/  
/*  
    run_chi2_d0_mu_tag("as117");
    run_chi2_d0_mu_tag("as119");
    run_chi2_d0_mu_tag("m166v5");
    run_chi2_d0_mu_tag("m169v5");
    run_chi2_d0_mu_tag("m171v5");
    run_chi2_d0_mu_tag("m173v5");
    run_chi2_d0_mu_tag("m175v5");
    run_chi2_d0_mu_tag("m178v5");
*/  

  }
  json += ("]");
  std::cout << json << std::endl;

}

void run_chi2_d0_mu_tag(TString name="") {
low=999.;
high=0;
//name=lname;
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
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_up"};//, "_1075", "_1125" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975, 1.055};//, 1.075, 1.125, 0.802};
/*
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975, 1.000, 1.055};
std::vector<TString> tune = {"_down", "_ddown", "_dddown", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.755, 0.775, 0.800, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055};
*/
//TCanvas *c1 = new TCanvas("c1","c1");
TH1F *chiTest = new TH1F("chiTest_"+name,TString::Format("chiTest_%s",name.Data()),400,0,2);
chiTest->Sumw2();
//chiTest->SetDirectory(0);
for(auto & it : tune) {
  int pos = &it - &tune[0];
  if(param[pos]<0.65 && !name.Contains("fsr-up")) continue;
  if(!fullpt && param[pos]>1 && !name.Contains("fsr-down")) continue;
  if(epoch==1 && param[pos]==0.875) continue; //FIXME
  if(name.Contains("m173v5") and it == "") continue; //FIXME
  std::cout << "Running on tune: " << it << std::endl;
  TH1F *shiftData;
  float chi = chi2_d0_mu_tag_test(shiftData, it, name, param[pos]);
  if(isnan(chi)) continue;
  if(chi<low) low = chi;
  if(chi>high) high = chi;
  chiTest->GetYaxis()->SetRangeUser(max(int(low)-1,0),int(high)+2);
  chiTest->SetBinContent(chiTest->FindBin(param[pos]),chi);
  //chiTest->SetBinError(chiTest->FindBin(param[pos]),1);
}

if(fullpt)
chiTest->GetXaxis()->SetRangeUser(0.65,1.125);
else
chiTest->GetXaxis()->SetRangeUser(0.65,0.976);//1.055);
if(name.Contains("fsr-down")) chiTest->GetXaxis()->SetRangeUser(0.65,1.126);
//chiTest->GetYaxis()->SetRangeUser(55,90);
chiTest->GetYaxis()->SetRangeUser(int(low)-1,int(high)+2);
//chiTest->GetYaxis()->SetRangeUser(200,220);
chiTest->SetMarkerStyle(20);
chiTest->Draw("p9");
std::cout << chiTest->GetName() << std::endl;
std::cout << chiTest->GetTitle() << std::endl;
tdr(chiTest,epoch,fin);
/*
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
*/
((TF1*)(gROOT->GetFunction("pol3")))->SetParameters(1., 1., 1., 1.);
//TFitResultPtr fit = chiTest->Fit("pol3","FSEMQ","",0.6,1.055);
if(fullpt)
chiTest->Fit("pol3","FSMEQW");//,"",0.6,0.976);
else
chiTest->Fit("pol3","FSMEQRW","",0.6,0.976);
if(name.Contains("fsr-down")) chiTest->Fit("pol3","FSMEQRW","",0.6,1.126);//0.976);
//TFitResultPtr fit = chiTest->Fit("pol3","FSEMQ","",0.6,0.975);
//TFitResultPtr fit = chiTest->Fit("pol2","FSMEQ");
//TFitResultPtr fit = chiTest->Fit("pol2","FSMEQ","",0.8,1.0);
/*
float min = (-1)*fit->Parameter(1)/(2*fit->Parameter(2));
float chimin = fit->Parameter(0) + fit->Parameter(1)*min + fit->Parameter(2) * pow(min,2);
float err = (-1)*fit->Parameter(1) / (2 * fit->Parameter(2)) - sqrt(pow(fit->Parameter(1),2)
            - 4 * fit->Parameter(2) * (fit->Parameter(0) - chimin - 1)) / (2 * fit->Parameter(2));
*/
float min = chiTest->GetFunction("pol3")->GetMinimumX(0.6,1.126);
chiTest->GetFunction("pol3")->Print("v");
//float chimin = fit->Parameter(0) + fit->Parameter(1)*min + fit->Parameter(2) * pow(min,2) + fit->Parameter(3) * pow(min,3);
float chimin = chiTest->GetFunction("pol3")->Eval(min);
float err = chiTest->GetFunction("pol3")->GetX(chimin+1,0.6,1.126);
if(name=="") { nom=min; nerr=err; }
report = Form("Minimum at x= %0.3f +/- %0.3f",min, abs(min-err));
json += Form("%.4f, %.4f, ",min,abs(min-err));
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
if(!TDR) pt->Draw();
gStyle->SetOptStat(0);

if(name.Length()>0) name = "_" + name;
name += epoch_name[epoch];
if(fullpt) name += "_jpT";
if(fin) {
c1->SaveAs("chi2_d0_mu_tag"+name+"_final.pdf");
c1->SaveAs("chi2_d0_mu_tag"+name+"_final.png");
}
else {
c1->SaveAs("chi2_d0_mu_tag"+name+".pdf");
c1->SaveAs("chi2_d0_mu_tag"+name+".png");
}

delete pt;
chiTest->Delete();
//delete chiTest;
//delete c1;
}

float chi2_d0_mu_tag_test(TH1F *&shiftData, TString tune="", TString name="", float num=0.855) {
TH1F *data, *mc;
TH1F *tmpData, *tmpData2;
if(epoch!=0) {
getHist(name, tune, data, mc, epoch);
}
else {
getHist(name, tune, data, mc, 1, false);
TH1F *data2, *mc2;
getHist(name, tune, data2, mc2, 2, false);
data->Add(data2);
mc->Add(mc2);
mc->GetXaxis()->SetTitle("(D^{0} #it{p}_{T} + #mu #it{p}_{T}) / #Sigma_{ch} #it{p}_{T}");
data->GetXaxis()->SetTitle("(D^{0} #it{p}_{T} + #mu #it{p}_{T}) / #Sigma_{ch} #it{p}_{T}");
if(name.Contains("_toy") && tune == "_sdown") { //only compute modified hist once
    std::cout << "here1" << std::endl;
    if((name == "" || name == "172v5")) { //MC errors otherwise (e.g. FSR)
    std::cout << "here2" << std::endl;
      std::cout << "loading MC for toy data!" << std::endl;
      TH1F *data2;
      std::cout << "creating toy data" << std::endl;
    std::cout << "here3" << std::endl;
      getHist(name, tune, data, mc, 1, false);
    std::cout << "here4" << std::endl;
      getHist(name, tune, data2, mc2, 2, false);
    std::cout << "here5" << std::endl;
      mc->Add(mc2);
      std::cout << "loading MC for toy data DONE!" << std::endl;
      data->Add(data2);
      delete data2;
    }
    /*
    else
      data = (TH1F*)tmpData->Clone("");
    */
    data->SetTitle("frame_ptfrac_toyData_hist");
    std::cout << data->GetTitle() << std::endl;
    std::cout << "creating toy data DONE!" << std::endl;
  std::cout << data->GetTitle() << std::endl;
  //Clone toy data for shifting in each for loop (one shift per for-loop iteration)
    std::cout << "here6" << std::endl;
  shiftData = (TH1F*)data->Clone("");
    std::cout << "here7" << std::endl;

  //vary MC bin content by data uncertainty
  TRandom3 *rand = new TRandom3(0);
  for(int i = 1; i <= data->GetNbinsX(); i++) {
    float y = shiftData->GetBinContent(i);
    float e(0.);
    e = shiftData->GetBinError(i);
    /*
    if(name == "") //Data errors for nominal MC only
      e = data->GetBinError(i);
    else //MC errors otherwise (e.g. FSR)
    */
      //toyData sample has same number of events as data, similar statistics
      //e = shiftData->GetBinError(i);
    std::cout << y << " +/- " << e << " ";
    //shift each bin by random amount samples by the appropriate bin error
    if(!(name == "" || name == "172v5")) //MC errors otherwise (e.g. FSR)
      y += rand->Gaus(0, e);
    std::cout << y << std::endl;
    shiftData->SetBinContent(i,y);
    //if(name.Length()==0) shiftData->SetBinError(i,shiftData->GetBinError(i) * sqrt(float(nmc)/float(nentries)));
    //shiftData->SetBinError(i,e);
  }
  mcHists.push_back((TH1F*)shiftData->Clone());
  mcHists.back()->SetDirectory(0);
  delete rand;
  shiftData->GetXaxis()->SetRangeUser(0,1.1);
  std::cout << data->GetEntries() << " " << shiftData->GetEntries() << std::endl;
  //shiftData->Scale(data->Integral()/shiftData->Integral());
}
else if(name.Contains("_toy")) { //only compute modified hist once
  data = (TH1F*)shiftData->Clone();
}
mcHists.push_back((TH1F*)mc->Clone());
mcHists.back()->SetDirectory(0);
delete data2;
delete mc2;
}
mc->Scale(1./mc->Integral());
data->Scale(1./data->Integral());
data->GetXaxis()->SetRangeUser(0.,0.95);
mc->GetXaxis()->SetRangeUser(0.,0.95);
mc->GetYaxis()->SetRangeUser(0.,.26);
data->GetYaxis()->SetRangeUser(0.,.26);
setupPad()->cd();
mc->Draw();
tdr(mc, epoch, fin);
if(epoch==0) {
gStyle->SetOptStat(0);
mc->SetMarkerStyle(20);
mc->SetMarkerColor(kBlack);
mc->SetLineColor(kBlack);
mc->SetLineWidth(2);
TString namet(name);
if(num==0) num=0.855;
if(namet == "") namet = "172v5";
//if(tunet == "") tunet = "855";
TString mcvdname(TString::Format("www/meson/morph/ptfrac/ptfrac_signal_%s_%d_BCDEFGH_d0mu",namet.Data(),int(num*1000)));
if(fullpt) mcvdname += "_jpT";
if(fin) mcvdname += "_final";
c1->SaveAs(mcvdname + ".pdf");
c1->SaveAs(mcvdname + ".png");

if(namet=="172v5" && num > 0.825 && num < 0.875) {
data->SetTitle("");
//data->GetXaxis()->SetRangeUser(0,1.1);
mc->SetMarkerStyle(20);
data->SetMarkerStyle(20);
data->SetMarkerColor(kBlack);
data->SetLineColor(kBlack);
data->SetLineWidth(2);
data->Draw();
tdr(data, epoch, fin);
TString ptname("www/meson/morph/ptfrac/ptfrac_signal_Data_BCDEFGH_d0mu");
if(fullpt) ptname += "_jpT";
if(fin) ptname += "_final";
c1->SaveAs(ptname + ".pdf");
c1->SaveAs(ptname + ".png");
}
}

/*
if(tune=="" && name=="") {
TCanvas *c1 = setupCanvas();
TPad *p1 = setupPad();
p1->cd();
data->Draw();
gStyle->SetOptStat(0);
tdr(data,0);
data->SetMarkerStyle(20);
data->SetMarkerColor(kBlack);
data->SetLineColor(kBlack);
data->SetLineWidth(2);
c1->SaveAs("ptfrac_signal_Data_"+name+"d0.pdf");
c1->SaveAs("ptfrac_signal_Data_"+name+"d0.png");

}
*/

//data->Scale(1./data->Integral());
//mc->Scale(1./mc->Integral());
data->SetLineColor(kBlack);
data->SetMarkerColor(kBlack);
data->SetMarkerStyle(20);
data->SetLineWidth(2);
mc->SetLineColor(kRed);
mc->SetMarkerColor(kRed);
mc->SetMarkerStyle(1);
mc->GetYaxis()->SetRangeUser(0.,.4);
data->GetYaxis()->SetRangeUser(0.,.4);
data->GetXaxis()->SetRangeUser(0.2,1.);
mc->GetXaxis()->SetRangeUser(0.2,1.);
if(fullpt) {
data->GetXaxis()->SetRangeUser(0.21,0.89);
mc->GetXaxis()->SetRangeUser(0.21,0.89);
}
gStyle->SetOptStat(0);
mc->SetTitle("");
data->SetTitle("");
mc->Draw("hist");
mc->Draw("same e");
data->Draw("same");
tdr(mc, epoch, fin);
if(num==0) num=0.855;
if(name=="") name="172v5";
mc->SetTitle(TString::Format("mcVdata_%s_%d_d0mu_tag",name.Data(),(int)(num*1000)) + epoch_name[epoch] + "_d0mu");
TString title(TString::Format("mcVdata_%s_%d_d0mu_tag",name.Data(),(int)(num*1000)) + epoch_name[epoch] + "_d0mu");
if(fullpt) title += "_jpT";
if(fin) {
title += "_final";
}
/*
else {
c1->SaveAs(TString::Format("mcVdata_%s_%d_d0mu_tag",name.Data(),(int)(num*1000)) + epoch_name[epoch] + "_d0mu.pdf");
c1->SaveAs(TString::Format("mcVdata_%s_%d_d0mu_tag",name.Data(),(int)(num*1000)) + epoch_name[epoch] + "_d0mu.png");
}
*/
c1->SaveAs(title + ".pdf");
c1->SaveAs(title + ".png");
data->GetXaxis()->SetRangeUser(0.2,1.);//0.985);
mc->GetXaxis()->SetRangeUser(0.2,1.);//0.985);
if(fullpt) {
data->GetXaxis()->SetRangeUser(0.,0.89);
mc->GetXaxis()->SetRangeUser(0.,0.89);
}
float chi2 = data->Chi2Test(mc, "CHI2 P WW");
std::cout << tune << " Chi2= " << chi2 << std::endl;
if(chi2<low) low = chi2;
if(chi2>high) high = chi2;

delete data;
delete mc;

return chi2;

}
