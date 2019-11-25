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
//#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/convert.h"
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/tdr.C"
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/splot_d0.C"

using namespace RooFit;
//TString name("");
float low(50.), high(50.),nom(0.8103),nerr(0.05);
bool TDR(1);
int epoch(0);
bool fullpt(0);
TString epoch_name[3] = {"_BCDEFGH", "_BCDEF", "_GH"};

float N(1.);
TCanvas *c1 = setupCanvas();
TString report("");
TString json("\"d0\" :      [");
float chi2_d0_toy_test(TH1F *&data, TString sample="d0", TString tune="", TString name="", float num=0.855, int iteration=0);
//float chi2_d0_toy_test(TH1F *&data, TH1F *&mc, TString sample="d0", TString tune="", TString name="", float num=0.855);
void run_chi2_d0_toy(TString, int iteration);
//void run_chi2_d0_toy(std::vector<TH1F*>&, std::vector<TH1F*>&, TString);
RooRealVar ptfrac;

//holds MC for shifting toy models, updates ones per for-loop iteration
TH1F *shiftData;
  //std::vector<int> ndata = {0, 0, 0};
  //std::vector<int> nmc = {0, 0, 0};


//Holds MC hists for repeated runs in for loop
std::vector<TH1F*> hists;
//std::vector<TH1F*> hists;// = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

//TH1F *hists[] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
//std::vector<TH1F*> hists;
//polulates initial with nullptr for each rb template (saves time by not opening MC files in each for-loop iteration)
/*
for(size_t i = 0; i < 14; i++) {
  hists.push_back(nullptr);
}
*/

//Holds toy data before shifts (MC rb=0.855)
TH1F* data[] = {nullptr,nullptr,nullptr};

void getHist(int &nentries, int &nmc, TString sample, TString name, TString tune, TH1F *&pdata, TH1F *&mc, int epoch, bool norm=true, bool toyData=false, int iteration=0) {
if(name == "MC") name = "172v5";
//std::vector<float> bin;
RooBinning bins(0,1.1);
/*
if(sample.Contains("mu_tag"))
  bin = {0, 0.2, 0.4, 0.6, 0.7, 0.75, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0};
else if(sample.Contains("d0"))
bin = {0, 0.2, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
else if(sample.Contains("d0"))
  bin = {-0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.875, 0.95, 1.0};
*/
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/param.h"
TString fname = TString::Format("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_%s_sPlot_%s.root",name.Data(),sample.Data());
if(toyData) {
  splot_d0(pdata, TString::Format("toyData%d",iteration), false, "", epoch, false);
  //fname = TString::Format("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_toyData%d_sPlot_%s.root",iteration,sample.Data());
}
if(fullpt) fname.ReplaceAll(".root","_jpT.root");
if(epoch>0) fname.ReplaceAll(".root",TString::Format("%d.root",epoch));
std::cout << fname << std::endl;
TFile *fdata = TFile::Open(fname);
if(name.Length()==0)
fname = TString::Format("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_172v5%s_sPlot_%s.root",tune.Data(),sample.Data());
//fname = TString::Format("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/sPlot/sPlot/morph/TopMass_172v5%s_sPlot_d0.root",tune.Data());
else
fname = TString::Format("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_%s%s_sPlot_%s.root",name.Data(),tune.Data(),sample.Data());
//fname = TString::Format("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/sPlot/sPlot/morph/TopMass_%s%s_sPlot_d0.root",name.Data(),tune.Data());
if(epoch>0) fname.ReplaceAll(".root",TString::Format("%d.root",epoch));
if(fullpt) fname.ReplaceAll(".root","_jpT.root");
TFile *fmc = TFile::Open(fname);

RooPlot *tmp = nullptr;
TString Tptfrac("ptfrac_signal");
if(sample.Contains("mu_tag")) Tptfrac = "ptfrac_mu_tag_signal";
//bin = {0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675, 0.75, 0.825, 0.9, 0.975, 1.0};
//bin = {-0.025, 0.025, 0.075, 0.125,  0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.0};
//bin = {0, 0.05, 0.1, 0.15,  0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 1.0};
for(int i = 0; i < bin.size(); i++) {
 bins.addBoundary(bin[i]);
}
tmp = (RooPlot*)fdata->Get(Tptfrac)->Clone(TString::Format("ptfrac_signal_Data"));
if(!toyData) pdata = (TH1F*)convert(tmp, norm, bin);
if(sample.Contains("mu_tag") || sample.Contains("d0")) {
  //std::cout << "variable binning" << std::endl;
  delete tmp;
  tmp = ((RooWorkspace*)fdata->Get("w"))->var("ptfrac")->frame();
  ((RooDataSet*)((RooWorkspace*)fdata->Get("w"))->data("sigData"))->plotOn(tmp, RooFit::Binning(bins), DataError(RooAbsData::SumW2));
  if(!toyData) pdata = (TH1F*)convert(tmp, norm, bin);
}
std::cout << fname << std::endl;
pdata->SetTitle(pdata->GetName());
delete tmp;
/*
if(fmc->Get("ptfrac")->GetBinContent(0) != 0) {
  mc = (TH1F*)fmc->Get("ptfrac");
  return;
}
*/
tmp = (RooPlot*)fmc->Get(Tptfrac)->Clone(TString::Format("ptfrac_signal_mc%s%s",name.Data(),tune.Data()));
//if(nmc == 0) ((RooWorkspace*)fmc->Get("w"))->data("sigData")->numEntries();
//tmp = ((RooWorkspace*)fmc->Get("w"))->var("ptfrac")->frame();
//((RooDataSet*)((RooWorkspace*)fmc->Get("w"))->data("sigData"))->plotOn(tmp, RooFit::Binning(bins), DataError(RooAbsData::SumW2));
if(tmp==nullptr) {std::cout << fname << "prfrac" << std::endl; return;}
mc = (TH1F*)convert(tmp, norm, bin);
//mc = (TH1F*)convert(tmp, norm, bin);
if(sample.Contains("mu_tag") || sample.Contains("d0")) {
  //std::cout << "variable binning" << std::endl;
  delete tmp;
  tmp = ((RooWorkspace*)fmc->Get("w"))->var("ptfrac")->frame();
  ((RooDataSet*)((RooWorkspace*)fmc->Get("w"))->data("sigData"))->plotOn(tmp, RooFit::Binning(bins), DataError(RooAbsData::SumW2));
  mc = (TH1F*)convert(tmp, norm, bin);
}
pdata->SetDirectory(0);
mc->SetDirectory(0);
mc->SetTitle(mc->GetName());
delete tmp;

fdata->Close();
fmc->Close();
delete fdata;
delete fmc;
}

void chi2_d0_toy(TString set="", int queue=0) {
  int max(1);
  for(int i = 0; i < max; i++) {
    std::cout << std::endl << "iteration " << i+1 << "/" << max << std::endl << std::endl;
    run_chi2_d0_toy(set, i + max*queue);
  }
/*
  run_chi2_d0_toy("isr-down" ndata, nmc);
  run_chi2_d0_toy("isr-up" ndata, nmc);
  run_chi2_d0_toy("fsr-down" ndata, nmc);
  run_chi2_d0_toy("fsr-up" ndata, nmc);
  run_chi2_d0_toy("uedown" ndata, nmc);
  run_chi2_d0_toy("ueup" ndata, nmc);
  //run_chi2_d0_toy("erdON" ndata, nmc);
  run_chi2_d0_toy("GluonMove_erdON" ndata, nmc);
  //run_chi2_d0_toy("QCD_erdON" ndata, nmc);
  std::vector<TString> syst = {"LEP", "PU", "PI", "TRIGGER", "JER", "JSF" }; //no lepton tracker efficiencies are used!
  //std::vector<TString> syst = {"TRK", "LEP", "PU", "PI", "TRIGGER", "JER" };
  for(auto & it : syst) {
    run_chi2_d0_toy("down_"+it ndata, nmc);
    run_chi2_d0_toy("up_"+it ndata, nmc);
  }
  run_chi2_d0_toy("hdampdown" ndata, nmc);
  run_chi2_d0_toy("hdampup" ndata, nmc);
  run_chi2_d0_toy("tpt" ndata, nmc);
  run_chi2_d0_toy("bkg" ndata, nmc);
  run_chi2_d0_toy("as117" ndata, nmc);
  run_chi2_d0_toy("as119" ndata, nmc);
  run_chi2_d0_toy("m166v5" ndata, nmc);
  run_chi2_d0_toy("m169v5" ndata, nmc);
  run_chi2_d0_toy("m171v5" ndata, nmc);
  run_chi2_d0_toy("m173v5" ndata, nmc);
  run_chi2_d0_toy("m175v5" ndata, nmc);
  run_chi2_d0_toy("m178v5" ndata, nmc);
*/

  json += ("],");
  std::cout << json << std::endl;

}

void run_chi2_d0_toy(TString name="", int iteration=0) {
//void run_chi2_d0_toy(std::vector<int> &nentries, std::vector<int> &nmc, TString name="") {
gROOT->Reset();
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
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_up" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975, 1.055, 0.802};
/*
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975, 1.000, 1.055};
std::vector<TString> tune = {"_down", "_ddown", "_dddown", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.755, 0.775, 0.800, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055};
*/
//TCanvas *c1 = new TCanvas("c1","c1");
//TCanvas *c1 = setupCanvas();
TH1F *chiTest = new TH1F("chiTest_"+name,TString::Format("chiTest_%s",name.Data()),400,0,2);
chiTest->Sumw2();
//chiTest->SetDirectory(0);

//for(size_t pos = 0; pos < tune.size(); pos++) {
for(auto & it : tune) {
  //TString it = tune[pos];
  int pos = &it - &tune[0];
  if(param[pos]>1) continue;
  if(name == "up_PI" && param[pos]>0.975 && fullpt) continue; //remove up_PI 0.975 with large chi^2
  if(name == "GluonMove_erdON" && param[pos]==0.900 && epoch==0) continue; //remove up_PI 0.975 with large chi^2
  if(name == "GluonMove_erdON" && param[pos]<0.700 && epoch==1) continue; //remove up_PI 0.975 with large chi^2
  //if(name == "down_PU" && it=="_down" && epoch==2) continue; //remove down_UP 0.755 FIXME
  
  std::cout << "Running on tune: " << it << std::endl;
  /*
  float chi = chi2_d0_toy_test(data[0], "d0", it, name, param[pos]);
  chi += chi2_d0_toy_test(data[1], "d0", it, name, param[pos]);
  chi += chi2_d0_toy_test(data[2], "d0", it, name, param[pos]);
  */
  //if(hists[pos] != nullptr) std::cout << hists[pos]->GetTitle() << std::endl;
  float chi = chi2_d0_toy_test(data[1], "d0", it, name, param[pos], iteration);
  //float chi = chi2_d0_toy_test(data[1], hists[pos], "d0", it, name, param[pos]);
  std::cout << data[1]->GetTitle() << std::endl;
  //if(hists[pos] != nullptr) std::cout << hists[pos]->GetTitle() << std::endl;
  if(chi<low) low = chi;
  if(chi>high) high = chi;
  chiTest->GetYaxis()->SetRangeUser(int(low)-1,int(high)+2);
  chiTest->SetBinContent(chiTest->FindBin(param[pos]),chi);
  //chiTest->SetBinError(chiTest->FindBin(param[pos]),sqrt(1./N));
}

//chiTest->GetXaxis()->SetRangeUser(0.65,1.055);
chiTest->GetXaxis()->SetRangeUser(0.65,0.976);//1.055);
//chiTest->GetYaxis()->SetRangeUser(55,90);
chiTest->GetYaxis()->SetRangeUser(int(low)-1,int(high)+2);
//chiTest->GetYaxis()->SetRangeUser(200,220);
chiTest->SetMarkerStyle(20);
chiTest->Draw("p9");
std::cout << chiTest->GetName() << std::endl;
std::cout << chiTest->GetTitle() << std::endl;
tdr(chiTest,epoch);
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
chiTest->Fit("pol3","FSMEQRW","",0.6,0.976);
//TFitResultPtr fit = chiTest->Fit("pol3","FSEMQ","",0.6,0.975);
//TFitResultPtr fit = chiTest->Fit("pol2","FSMEQ");
//TFitResultPtr fit = chiTest->Fit("pol2","FSMEQ","",0.8,1.0);
/*
float min = (-1)*fit->Parameter(1)/(2*fit->Parameter(2));
float chimin = fit->Parameter(0) + fit->Parameter(1)*min + fit->Parameter(2) * pow(min,2);
float err = (-1)*fit->Parameter(1) / (2 * fit->Parameter(2)) - sqrt(pow(fit->Parameter(1),2)
            - 4 * fit->Parameter(2) * (fit->Parameter(0) - chimin - 1)) / (2 * fit->Parameter(2));
*/
float min = chiTest->GetFunction("pol3")->GetMinimumX(0.6,1.075);
//float chimin = fit->Parameter(0) + fit->Parameter(1)*min + fit->Parameter(2) * pow(min,2) + fit->Parameter(3) * pow(min,3);
float chimin = chiTest->GetFunction("pol3")->Eval(min);
float err = chiTest->GetFunction("pol3")->GetX(chimin+1,0.6,1.075);
if(name=="") { nom=min; nerr=err; }
report = Form("Minimum at x= %g +/- %0.6g",min, abs(min-err));
json += Form("%.4f, %.4f, ",min,abs(min-err));
//json += Form("%.4f, %.4f, ",min,abs(min-err));
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
//c1->SaveAs("chi2_d0"+name+"_toy.pdf");
//c1->SaveAs("chi2_sim"+name+"_toy.png");

delete pt;
//chiTest->Delete();
//delete chiTest;
//delete c1;
}

float chi2_d0_toy_test(TH1F *&data, TString sample="d0", TString tune="", TString name="", float num=0.855, int iteration=0) {
//float chi2_d0_toy_test(TH1F *&data, TH1F *&mc, TString sample="d0", TString tune="", TString name="", float num=0.855) {
if((sample.Contains("d0_mu") || sample.Contains("d0")) && name.Contains("GluonMove_erdON"))
  name = "erdON";
TH1F *mc, *mc2;
TH1F *tmpData, *tmpData2;
int nentries(0),nmc(0);
//TH1F *data, *data2, *mc, *mc2;
if(epoch>0) {
getHist(nentries, nmc, sample, name, tune, data, mc, epoch, 1, iteration);
}
else {
//only load each MC hist once (instead of each time in for loop)
//if(mc == nullptr) {
  std::cout << "loading MC" << std::endl;
  getHist(nentries, nmc, sample, name, tune, tmpData, mc, 1, false, false, iteration);
  getHist(nentries, nmc, sample, name, tune, tmpData2, mc2, 2, false, false, iteration);
  mc->Add(mc2);
  mc = (TH1F*)mc2->Clone();
  std::cout << "loading MC DONE!" << std::endl;
  tmpData->Add(tmpData2);
  tmpData = (TH1F*)tmpData2->Clone();
//}
if(tune == "_sdown") { //only compute modified hist once
    if((name == "" || name == "172v5")) { //MC errors otherwise (e.g. FSR)
  //only load each nominal MC hist once (instead of each time in for loop)
  //if(mc == nullptr || data == nullptr) {
    std::cout << "loading MC for toy data!" << std::endl;
    TH1F *data2;
    std::cout << "creating toy data" << std::endl;
    //delete mc;
    //delete mc2;
    getHist(nentries,nmc,sample, name, tune, data, mc, 1, false, true, iteration);
    getHist(nentries, nmc, sample, name, tune, data2, mc2, 2, false, true, iteration);
    mc->Add(mc2);
  mc = (TH1F*)mc2->Clone();
    std::cout << "loading MC for toy data DONE!" << std::endl;
    data->Add(data2);
  data = (TH1F*)data2->Clone();
    delete data2;
    }
    else
    data = (TH1F*)tmpData->Clone("");
    data->SetTitle("frame_ptfrac_toyData_hist");
    std::cout << data->GetTitle() << std::endl;
    std::cout << "creating toy data DONE!" << std::endl;
  //}
  std::cout << data->GetTitle() << std::endl;
  //scale hist to ensure the number of events hasn't changed
  //Clone toy data for shifting in each for loop (one shift per for-loop iteration)
  shiftData = (TH1F*)data->Clone("");

  //vary MC bin content by data uncertainty
  TRandom3 *rand = new TRandom3(0);
  for(int i = 1; i <= tmpData->GetNbinsX(); i++) {
    float y = shiftData->GetBinContent(i);
    float e(0.);
    e = shiftData->GetBinError(i);
    /*
    if(name == "") //Data errors for nominal MC only
      e = tmpData->GetBinError(i);
    else //MC errors otherwise (e.g. FSR)
    */
      //toyData sample has same number of events as data, similar statistics
      //e = shiftData->GetBinError(i);
    //shift each bin by random amount samples by the appropriate bin error
    if(!(name == "" || name == "172v5")) { //MC errors otherwise (e.g. FSR)
    std::cout << y << " +/- " << e << " ";
      y += rand->Gaus(0, e);
    std::cout << y << std::endl;
    shiftData->SetBinContent(i,y);
    }
    //if(name.Length()==0) shiftData->SetBinError(i,shiftData->GetBinError(i) * sqrt(float(nmc)/float(nentries)));
    //shiftData->SetBinError(i,e);
  }
  delete rand;
  shiftData->GetXaxis()->SetRangeUser(0,1.1);
  std::cout << data->GetEntries() << " " << shiftData->GetEntries() << std::endl;
  //shiftData->Scale(data->Integral()/shiftData->Integral());
}

delete tmpData;
delete tmpData2;
TCanvas *c1 = setupCanvas();
bool toy = tune =="_sdown" ? true : false;
std::cout << shiftData->GetTitle() << std::endl;
std::cout << "using toy data" << std::endl;
data->SetDirectory(0);
mc->SetDirectory(0);
mc->GetXaxis()->SetTitle("D^{0} #it{p}_{T} / #Sigma #it{p}_{T}^{ch}");
shiftData->GetXaxis()->SetTitle("D^{0} #it{p}_{T} / #Sigma #it{p}_{T}^{ch}");
delete mc2;
mc->Scale(1./mc->Integral());
if(toy) shiftData->Scale(1./shiftData->Integral());
setupPad()->cd();
tdr(mc, epoch);
if(fullpt) mc->GetXaxis()->SetTitle("D^{0} #it{p}_{T} / jet #it{p}_{T}");
mc->Draw();
tdr(mc, epoch);
if(epoch==0) {
gStyle->SetOptStat(0);
TString namet(name);
shiftData->SetMarkerStyle(20);
shiftData->SetMarkerColor(kBlack);
shiftData->SetLineColor(kBlack);
shiftData->SetLineWidth(2);
if(num==0) num=0.855;
if(namet == "") namet = "172v5";
//if(tunet == "") tunet = "855";
//c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_signal_%s_%d%s_sim%s.pdf",namet.Data(),int(num*1000), epoch_name[epoch].Data(), (fullpt ? "_jpT" : "")));
//c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_signal_%s_%d%s_sim%s.png",namet.Data(),int(num*1000), epoch_name[epoch].Data(), (fullpt ? "_jpT" : "")));

if(namet=="172v5" && num < 0.875) {
//if(namet=="172v5" && num > 0.825 && num < 0.875) {
shiftData->SetTitle("");
shiftData->GetXaxis()->SetRangeUser(0,1.1);
mc->SetMarkerStyle(20);
shiftData->SetMarkerStyle(20);
shiftData->SetMarkerColor(kBlack);
shiftData->SetLineColor(kBlack);
shiftData->SetLineWidth(2);
shiftData->GetYaxis()->SetRangeUser(0,.16);
tdr(data, epoch);
shiftData->Draw();
tdr(shiftData, epoch);
//c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_signal_toy_BCDEFGH_d0%d.pdf",iteration));
//c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_signal_toy_BCDEFGH_d0%d.png",iteration));
}
}

/*
if(tune=="" && name=="") {
TCanvas *c1 = setupCanvas();
TPad *p1 = setupPad();
p1->cd();
shiftData->Draw();
gStyle->SetOptStat(0);
tdr(data,0);
shiftData->SetMarkerStyle(20);
shiftData->SetMarkerColor(kBlack);
shiftData->SetLineColor(kBlack);
shiftData->SetLineWidth(2);
c1->SaveAs("ptfrac_signal_Data_"+name+"d0.pdf");
c1->SaveAs("ptfrac_signal_Data_"+name+"d0.png");

}
*/

N = mc->Integral();
}

shiftData->GetXaxis()->SetRangeUser(0.2,0.975);
mc->GetXaxis()->SetRangeUser(0.2,0.975);
if(fullpt) {
shiftData->GetXaxis()->SetRangeUser(0.0,0.7);
mc->GetXaxis()->SetRangeUser(0.0,0.7);
}
shiftData->SetLineColor(kBlack);
shiftData->SetMarkerColor(kBlack);
shiftData->SetMarkerStyle(20);
shiftData->SetLineWidth(2);
mc->SetLineColor(kRed);
mc->SetMarkerColor(kRed);
mc->SetMarkerStyle(1);
mc->SetLineWidth(1);
mc->GetYaxis()->SetRangeUser(0.,.16);
shiftData->GetYaxis()->SetRangeUser(0.,.16);
if(fullpt) {
mc->GetYaxis()->SetRangeUser(0.,.16);
shiftData->GetYaxis()->SetRangeUser(0.,.16);
}
/*
TCanvas *c1 = setupCanvas();
setupPad()->cd();
*/
mc->Draw("hist");
tdr(mc, epoch);
mc->Draw("same e");
shiftData->Draw("same");
if(num==0) num=0.855;
if(name=="") name="172v5";
TString mcvname(TString::Format("mcVdata_%s_%d_d0",name.Data(),(int)(num*1000)) + epoch_name[epoch]);
if(fullpt) mcvname += "_jpT";
//c1->SaveAs(mcvname + "_toy.pdf");
//c1->SaveAs(mcvname + "_toy.png");
float chi2 = shiftData->Chi2Test(mc, "CHI2 P WW");
/*
chi2 = 0.;
float sum1(0.);
float sum2(0.);
for(int i = 0; i < shiftData->GetNbinsX(); i++) {
  float exp = mc->GetBinContent(i);
  float obs = shiftData->GetBinContent(i);
  sum1 += obs;
  sum2 += exp;
  
}
std::cout << sum1 << " " << sum2 << std::endl;
float ndf = shiftData->GetXaxis()->GetLast() - shiftData->GetXaxis()->GetFirst();
for(int i = 0; i < shiftData->GetNbinsX(); i++) {
  float cnt1 = shiftData->GetBinContent(i);
  float cnt2 = mc->GetBinContent(i);
  float e1sq = pow(shiftData->GetBinError(i),2);
  float e2sq = pow(mc->GetBinError(i),2);
  if (cnt1 * cnt1 == 0 && cnt2 * cnt2 == 0){  continue; }
  //if (cnt1 * cnt1 == 0 && cnt2 * cnt2 == 0){ ndf--;  continue; }
  float delta = sum1 * cnt2 - sum2 * cnt1;
  float sigma = sum1 * sum1 * e2sq + sum2 * sum2 *e1sq;
  chi2 += delta * delta / sigma;
}
std::cout << ndf << std::endl;
*/
std::cout << tune << " Chi2= " << chi2 << std::endl;
if(chi2<low) low = chi2;
if(chi2>high) high = chi2;

/*
delete data;
delete mc;
*/

return chi2;

}
