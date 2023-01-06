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
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/tdr.C"
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/param.h"

using namespace RooFit;
//TString name("");
float low(50.), high(50.),nom(0.8103),nerr(0.05),pre(0);
TString errtxt;
bool TDR(1);
int epoch(0);
bool fullpt(0);
bool fin(false);
TString epoch_name[4] = {"_xb", "_BCDEFGH", "_BCDEF", "_GH"};

float N(1.);
TCanvas *c1 = setupCanvas();
TString report("");
TString json("\"j\" :      [");
TString hepdata("");
TString rbtest[3] = {"", "", ""};
float chi2_jpsi_(TH1F *& shiftData, TString tune="", TString name="", float num=0.855, int syst = 0, TString postfix="");
void run_chi2_jpsi(TString, int syst = 0, TString postfix = "");
RooRealVar ptfrac;

void getHist(TString name, TString tune, TH1F *&data, TH1F *&mc, TH1F *&bkg, int epoch, bool norm=true, int syst = 0) {
TString fname = TString::Format("/eos/cms/store/group/phys_top/byates/sPlot//TopMass_Data_sPlot_jpsi.root");
//TString fname = TString::Format("/eos/cms/store/group/phys_top/byates/sPlot//TopMass_Data_var_sPlot_jpsi.root");
if(name.Contains("_bkg")) fname = TString::Format("/eos/cms/store/group/phys_top/byates/sPlot//TopMass_Data_bkg_sPlot_jpsi.root");
if(name.Contains("_toys")) fname.ReplaceAll("_toys","");
if(epoch>0) fname.ReplaceAll(".root",TString::Format("%d.root",epoch));
else if(epoch<0) fname.ReplaceAll(".root","_xb.root");
if(name.Contains("_xb")) {
fname.ReplaceAll("_xb","");
fname.ReplaceAll("_noSmooth","");
fname.ReplaceAll(".root","_xb.root");
}
if(name.Contains("noDup_shift")) fname.ReplaceAll("Data", "Data_noDup_shift");
else if(name.Contains("noDup")) fname.ReplaceAll("Data", "Data_noDup");
if(name.Contains("noTrkSF")) fname.ReplaceAll("Data", "Data_noTrkSF");
if(name.Contains("ep")) fname.ReplaceAll("Data", "Data_ep");
if(name.Contains("noHT")) fname.ReplaceAll("Data", "Data_noHT");
if(name.Contains("jpsikk")) fname.ReplaceAll("Data", "Data_jpsikk");
//if(name.Contains("FSR")) fname.ReplaceAll("Data","FSR");
if(name.Contains("FSR")) fname.ReplaceAll("Data","FSR");
if(name.Contains("dupBest")) fname.ReplaceAll("Data","Data_dupBest");
if(name.Contains("BDzb")) fname.ReplaceAll("Data","Data_dupBest");
if(name.Contains("_ds")) fname.ReplaceAll("Data","Data_ds");
if(name.Contains("_dupNest")) fname.ReplaceAll("Data","Data_dupNest");
if(name.Contains("172v5_no")) fname.ReplaceAll("Data","Data_no");
if(name.Contains("172v5_var_up")) fname.ReplaceAll("Data","Data_var_up");
else if(name.Contains("172v5_var_down")) fname.ReplaceAll("Data","Data_var_down");
else if(name.Contains("172v5_var")) fname.ReplaceAll("Data","Data_var");
if(name.Contains("172v5_mass")) fname.ReplaceAll("Data","Data_mass");
if(name.Contains("172v5_width")) fname.ReplaceAll("Data","Data_width");
if(name.Contains("172v5_bestmass")) fname.ReplaceAll("Data","Data_bestmass");
if(name.Contains("172v5_randmass_var_hand")) fname.ReplaceAll("Data","Data_randmass_var_hand");
else if(name.Contains("172v5_randmass_var")) fname.ReplaceAll("Data","Data_randmass_var");
else if(name.Contains("172v5_randmass")) fname.ReplaceAll("Data","Data_randmass");
else if(name.Contains("172v5_Dz")) fname.ReplaceAll("Data","Data_randmass");
else if(name.Contains("172v5_Bfrag_rebin")) fname.ReplaceAll("Data","Data_rebin");
else if(name.Contains("172v5_Bfrag_test")) fname.ReplaceAll("Data","Data_randmass");
else if(name.Contains("172v5_Bfrag_genreco")) fname.ReplaceAll("Data","Data_Bfrag_genreco");
else if(name.Contains("Bfrag_genreco")) fname.ReplaceAll("Data","Data_Bfrag_genreco");
//else if(name.Contains("172v5_Bfrag_genreco")) fname.ReplaceAll("Data","Data_randmass");
else if(name.Contains("172v5_Bfrag_j5")) fname.ReplaceAll("Data","Data_Bfrag_genreco");
else if(name.Contains("172v5_Bfrag_csv")) fname.ReplaceAll("Data","Data_Bfrag_csv");
//lse if(name.Contains("172v5_Bfrag_genreco")) fname.ReplaceAll("Data","Data_randmass");//Bfrag_genreco");
else if(name.Contains("172v5_Bfrag")) fname.ReplaceAll("Data","Data_randmass");
if(name.Contains("172v5_test")) fname.ReplaceAll("Data","Data_randmass");
if(name.Contains("172v5_sumchptrand")) fname.ReplaceAll("Data","Data_randmass");
if(name.Contains("172v5_ctauup")) fname.ReplaceAll("Data","Data_ctauup");
if(name.Contains("172v5_ctaudown")) fname.ReplaceAll("Data","Data_ctaudown");
if(name.Contains("isr") || name.Contains("fsr") || name.Contains("ue") || name.Contains("erdON") || name.Contains("hdamp"))
  //fname.ReplaceAll("Data","172v5");
std::cout << name << std::endl;
if(fullpt) fname.ReplaceAll(".root","_jpT.root");
std::cout << fname << std::endl;
TFile *fdata = TFile::Open(fname);
if(name.Length()==0)
fname = TString::Format("/eos/cms/store/group/phys_top/byates/sPlot//morph/TopMass_172v5%s_sPlot_jpsi.root",tune.Data());
else {
/*
if(name.Contains("_no")) {
  name.ReplaceAll("_no", "");
  fname = TString::Format("/eos/cms/store/group/phys_top/byates/sPlot//TopMass_%s%s_no_sPlot_jpsi.root",name.Data(),tune.Data());
  std::cout << fname << std::endl;
}
else
*/
  fname = TString::Format("/eos/cms/store/group/phys_top/byates/sPlot//TopMass_%s%s_sPlot_jpsi.root",name.Data(),tune.Data());
}
if(name.Contains("_toys")) fname.ReplaceAll("_toys","");
if(epoch>0) fname.ReplaceAll(".root",TString::Format("%d.root",epoch));
else if(epoch<0) fname.ReplaceAll(".root","_xb.root");
if(name.Contains("_xb")) {
fname.ReplaceAll("_xb_","_");
fname.ReplaceAll("_noSmooth","");
fname.ReplaceAll(".root","_xb.root");
}
if(name.Contains("_charm")) {
fname.ReplaceAll("_charmup","");
fname.ReplaceAll("_charmdown","");
}
if(fullpt) fname.ReplaceAll(".root","_jpT.root");
std::cout << fname << std::endl;
TFile *fmc = TFile::Open(fname);

RooPlot *tmp = nullptr;
RooBinning bins(0,1.1);
/*
std::vector<float> bin;
bin = {-0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.875, 0.95, 1.0};
bin = {0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675, 0.75, 0.825, 0.9, 0.975, 1.0};
*/
//bin = {0.025, 0.1, 0.175, 0.25, 0.325, 0.4, 0.475, 0.55, 0.625, 0.7, 0.775, 0.85, 0.925, 0.975, 1.0};
//bin = {0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675, 0.75, 0.825, 0.9, 0.975, 1.0};
//bin = {-0.025, 0.025, 0.075, 0.125,  0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.0};
//bin = {0, 0.05, 0.1, 0.15,  0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 1.0};
bin = {0-0.025, 0.2-0.025, 0.4-0.025, 0.55-0.025, 0.6-0.025, 0.65-0.025, 0.7-0.025, 0.75-0.025, 0.8-0.025, 0.85-0.025, 0.9-0.025, 0.95-0.025, 1.0};
bin = {0-0.025, 0.2-0.025, 0.4-0.025, 0.55-0.025, 0.65-0.025, 0.75-0.025, 0.85-0.025, 0.95-0.025, 1.0};
if(fullpt)
bin = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0};
for(int i = 0; i < bin.size(); i++) {
 std::cout << "bin " << i << "= " << bin[i] << std::endl;
 bins.addBoundary(bin[i]);
}
//if(epoch<0) std::cout << "getting hist" << std::endl;
//if(epoch<0 || fname.Contains("_xb")) mc = (TH1F*)fmc->Get("ptfrac_signal_hist")->Clone();
if(epoch<0 || fname.Contains("_xb")) bkg = (TH1F*)fmc->Get("ptfrac_bkg_hist")->Clone();
bkg->SetDirectory(0);
if((epoch<0 || fname.Contains("_xb")) && syst > 0) { std::cout << "Loading smoothUp" << syst << std::endl; mc = (TH1F*)fmc->Get(TString::Format("ptfrac_signal_smoothUp%d", syst))->Clone(); }
//else if(name.Contains("_noSmooth")) mc = (TH1F*)fmc->Get("ptfrac_signal_hist")->Clone();
//else if((epoch<0 || fname.Contains("_xb"))) mc = (TH1F*)fmc->Get("ptfrac_signal_smooth")->Clone();
else if(epoch<0 || fname.Contains("_xb")) mc = (TH1F*)fmc->Get("ptfrac_signal_hist")->Clone();
//if(epoch<0 || fname.Contains("_xb")) mc->Add((TH1F*)fmc->Get("ptfrac_signal_bkgW")->Clone());
//else if(fname.Contains("no_sPlot"))  mc = (TH1F*)fmc->Get("ptfrac_signal")->Clone(TString::Format("ptfrac_signal_data%s%s",name.Data(),tune.Data()));
else if(fname.Contains("no_sPlot"))  mc->Scale(1/10.);
else { tmp = (RooPlot*)fmc->Get("ptfrac_signal")->Clone(TString::Format("ptfrac_signal_mc%s%s",name.Data(),tune.Data()));
//tmp = ((RooWorkspace*)fmc->Get("w"))->var("ptfrac")->frame();
((RooDataSet*)((RooWorkspace*)fmc->Get("w"))->data("sigData"))->plotOn(tmp, RooFit::Binning(bins), DataError(RooAbsData::SumW2));
if(tmp==nullptr) {std::cout << fname << std::endl; return;}
//mc = (TH1F*)convert(tmp, norm, 0, 1.1);
for(int i = 0; i < bin.size(); i++)
  std::cout << bin[i] << std::endl;
mc = (TH1F*)convert(tmp, norm, bin);
}

/*
if(tune.Length()>0 && !tune.Contains("Dz")){
std::cout << "here0" << std::endl;
auto bfragfin = TFile::Open("/afs/cern.ch/user/b/byates/TopAnalysis/data/era2016/bfragweights.root");
std::cout << "here1" << std::endl;
TH1F *hg = (TH1F*)bfragfin->Get(TString::Format("%sgen_jpsiFrag", tune.ReplaceAll("_","").Data()))->Clone();
//TH1F *hg = (TH1F*)bfragfin->Get(TString::Format("%sgen_jpsiFrag", tune.ReplaceAll("_","").Data()))->Clone();
std::cout << "here2" << std::endl;
hg->SetDirectory(0);
std::cout << "here3" << std::endl;
std::cout << hg->Integral() << "\t" << hg->GetName() << std::endl;
std::cout << "here4" << std::endl;
for(int ibin = 1; ibin < mc->GetNbinsX(); ibin++) {
std::cout << "here5" << std::endl;
  mc->SetBinContent(ibin, mc->GetBinContent(ibin) * hg->GetBinContent(ibin));
}
std::cout << "here6" << std::endl;
bfragfin->Close();
delete bfragfin;
delete hg;
}
*/
if(tune.Length()>0 && tune.Contains("Dz")){
auto bfragfin = TFile::Open("/afs/cern.ch/user/b/byates/TopAnalysis/data/era2016/bfragweights.root");
TGraph *hg = (TGraph*)bfragfin->Get(TString::Format("%sFrag", tune.ReplaceAll("_","").Data()))->Clone();
//TH1F *hg = (TH1F*)bfragfin->Get(TString::Format("%sgen_jpsiFrag", tune.ReplaceAll("_","").Data()))->Clone();
for(int ibin = 1; ibin < mc->GetNbinsX(); ibin++) {
  mc->SetBinContent(ibin, mc->GetBinContent(ibin) * hg->Eval(mc->GetBinCenter(ibin)));
}
bfragfin->Close();
delete bfragfin;
delete hg;
}
if(name.Contains("_charmup")) {
  std::cout << name << std::endl;
  auto fnomname = TString::Format("/eos/cms/store/group/phys_top/byates/sPlot//TopMass_172v5%s_sPlot_jpsi%d_xb.root",tune.Data(), epoch);
  fnomname.ReplaceAll("charmup","");
  TFile *fnom = TFile::Open(fnomname);
  TFile *fup = TFile::Open(TString::Format("/eos/cms/store/group/phys_top/byates/sPlot//TopMass_172v5_charmup_sPlot_jpsi%d_xb.root", epoch));
  std::cout << fnom->GetName() << std::endl;
  std::cout << fup->GetName() << std::endl;
  auto tmpnom = (TH1F*)fnom->Get("ptfrac_bkgc_hist")->Clone("nom");
  auto tmpup = (TH1F*)fup->Get("ptfrac_bkgc_hist")->Clone("charmup");
  for(int i = 0; i < tmpnom->GetNbinsX(); i++) {
    std::cout << i << ":\t" << tmpnom->GetBinContent(i) << "\t" << tmpup->GetBinContent(i) << std::endl;
  }
  mc->Add(tmpnom);
  mc->Add(tmpup,-1);
  mc->SetDirectory(0);
  delete tmpup;
  fnom->Close();
  fup->Close();
  delete fnom;
  delete fup;
}
if(name.Contains("_charmdown")) {
  auto fnom = TFile::Open(TString::Format("/eos/cms/store/group/phys_top/byates/sPlot//TopMass_172v5_sPlot_jpsi%d_xb.root",epoch));
  auto fdown = TFile::Open(TString::Format("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_172v5_charmdown_sPlot_jpsi%d_xb.root",epoch));
  std::cout << fnom->GetName() << std::endl;
  std::cout << fdown->GetName() << std::endl;
  auto tmpnom = (TH1F*)fnom->Get("ptfrac_bkgc_hist")->Clone("nom");
  auto tmpdown = (TH1F*)fdown->Get("ptfrac_bkgc_hist")->Clone("charmdown");
  mc->Add(tmpnom);
  mc->Add(tmpdown,-1);
  mc->SetDirectory(0);
  delete tmpdown;
  fnom->Close();
  fdown->Close();
  delete fnom;
  delete fdown;
}
mc->SetDirectory(0);
mc->SetTitle(mc->GetName());
std::cout << mc->GetTitle() << std::endl;
//mc->Rebin();
delete tmp;
//if(epoch<0) std::cout << "getting hist" << std::endl;
//if(epoch<0 || fname.Contains("_xb")) std::cout << "getting hist" << std::endl;
//if(epoch<0) std::cout << "getting hist" << std::endl;
if(epoch<0 || fname.Contains("_xb")) data = (TH1F*)fdata->Get("ptfrac_signal_hist")->Clone();
//if(epoch<0 || fname.Contains("_xb")) data = (TH1F*)fdata->Get("ptfrac_tot_hist")->Clone();
//if(epoch<0 || fname.Contains("_xb")) data->Add((TH1F*)fdata->Get("ptfrac_bkg_hist")->Clone());
else if(fname.Contains("no_sPlot"))  data = (TH1F*)fdata->Get("ptfrac_signal")->Clone(TString::Format("ptfrac_signal_data%s%s",name.Data(),tune.Data()));
else { tmp = (RooPlot*)fdata->Get("ptfrac_signal")->Clone(TString::Format("ptfrac_signal_data%s%s",name.Data(),tune.Data()));
//tmp = ((RooWorkspace*)fdata->Get("w"))->var("ptfrac")->frame();
((RooDataSet*)((RooWorkspace*)fdata->Get("w"))->data("sigData"))->plotOn(tmp, RooFit::Binning(bins), DataError(RooAbsData::SumW2));
//data = (TH1F*)convert(tmp, norm, 0, 1.1);
data = (TH1F*)convert(tmp, norm, bin);
}
data->SetDirectory(0);
data->SetTitle(data->GetName());
//data->Rebin();
std::cout << data->GetTitle() << std::endl;
delete tmp;

fdata->Close();
fmc->Close();
delete fdata;
delete fmc;
}

void chi2_jpsi(int ep=epoch, TString samp="", bool isFinal=false, int syst = 0, TString postfix="") {
  fin = isFinal;
  epoch = ep;
  if(samp != "")
  run_chi2_jpsi(samp, syst, postfix);
  else {
  run_chi2_jpsi("");
  run_chi2_jpsi("isr-down");
  run_chi2_jpsi("isr-up");
  run_chi2_jpsi("fsr-down");
  run_chi2_jpsi("fsr-up");
  run_chi2_jpsi("uedown");
  run_chi2_jpsi("ueup");
  //run_chi2_jpsi("erdON");
  run_chi2_jpsi("GluonMove_erdON");
  //run_chi2_jpsi("QCD_erdON");
  std::vector<TString> syst = {"LEP", "PU", "PI", "TRIGGER", "JER"};//, "JSF" }; //no lepton tracker efficiencies are used!
  //std::vector<TString> syst = {"TRK", "LEP", "PU", "PI", "TRIGGER", "JER" };
  for(auto & it : syst) {
    run_chi2_jpsi("down_"+it);
    run_chi2_jpsi("up_"+it);
  }
  run_chi2_jpsi("hdampdown");
  run_chi2_jpsi("hdampup");
  run_chi2_jpsi("172v5_Wdown");
  run_chi2_jpsi("172v5_Wup");
  run_chi2_jpsi("172v5_cdown");
  run_chi2_jpsi("172v5_cup");
  run_chi2_jpsi("fit-fcn");
  }
/*
  run_chi2_jpsi("tpt");
  run_chi2_jpsi("bkg");
  run_chi2_jpsi("bkg");
  run_chi2_jpsi("as117");
  run_chi2_jpsi("as119");
  */
  /*
  run_chi2_jpsi("m166v5");
  run_chi2_jpsi("m169v5");
  run_chi2_jpsi("m171v5");
  run_chi2_jpsi("m173v5");
  run_chi2_jpsi("m175v5");
  run_chi2_jpsi("m178v5");
  */

  json += ("],");
  std::cout << errtxt << std::endl;
  std::cout << json << std::endl;
  std::cout << hepdata << std::endl;
  /*
  std::cout << rbtest[0] << std::endl;
  std::cout << rbtest[1] << std::endl;
  std::cout << rbtest[2] << std::endl;
  */

}

void run_chi2_jpsi(TString name="", int syst = 0, TString postfix="") {
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
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_925", "_central", "_uuup" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.925, 0.955, 0.975};
//std::vector<TString> tune = {"_sdown", "_700", "_725", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup"};
//std::vector<float> param = {0.655, 0.700, 0.725, 0.775, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975};
//std::vector<TString> tune = {"_sdown", "_700", "_725", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup"};
//std::vector<float> param = {0.655, 0.700, 0.725, 0.775, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975};
//std::vector<TString> tune = {"_sdown", "_700", "_725", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup"};
//std::vector<float> param = {0.655, 0.700, 0.725, 0.775, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975};
/* FIXME these were for the B->D^0 GEN level weight cross checks
tune = {"_Dzdown", "_Dzdchargedinc", "", "_Dzuchargedinc", "_Dzup"};
param = {0.755, 0.810, 0.855, 0.906, 1.055};
tune = {"_Dzdown", "", "_Dzup"};
param = {0.755, 0.855, 1.055};
*/
//std::vector<TString> tune = {"_sdown", "_700", "_725", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup"};
//std::vector<float> param = {0.655, 0.700, 0.725, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975};
//std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup" };
//std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975};
//std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup" };
//std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975};
//std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup" };
//std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975};
//std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_up", "_1125" };
//std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975, 1.055, 1.125, 0.802};
/*
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_ddown", "_dddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.655, 0.700, 0.725, 0.755, 0.775, 0.800, 0.825, 0.855, 0.875, 0.900, 0.925, 0.955, 0.975, 1.000, 1.055};
std::vector<TString> tune = {"_down", "_ddown", "_dddown", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up" };
std::vector<float> param = {0.755, 0.775, 0.800, 0.855, 0.875, 0.900, 0.955, 0.975, 1.000, 1.055};
*/
//TCanvas *c1 = new TCanvas("c1","c1");
//TCanvas *c1 = setupCanvas();
TH1F *chiTest = new TH1F("chiTest_"+name,TString::Format("chiTest_%s;#it{r}_{b};Fit #chi^{2}",name.Data()),400,0,2);
//TH1F *chiTest = new TH1F("chiTest_"+name,TString::Format("chiTest_%s",name.Data()),400,0,2);
chiTest->Sumw2();
//chiTest->SetDirectory(0);
for(auto & it : tune) {
  int pos = &it - &tune[0];
  if(param[pos]>1.1) continue;
  //if(param[pos]<0.675) continue;
  //if(param[pos]>1 && !name.Contains("fsr-down")) continue;
  //if(name == "up_PI" && param[pos]>0.975 && fullpt) continue; //remove up_PI 0.975 with large chi^2
  if(name == "GluonMove_erdON" && param[pos]==0.900 && epoch==0) continue; //remove up_PI 0.975 with large chi^2
  if(name == "GluonMove_erdON" && param[pos]<0.700 && epoch==1) continue; //remove up_PI 0.975 with large chi^2
  //if(name == "down_PU" && it=="_down" && epoch==2) continue; //remove down_UP 0.755 FIXME
  //if(name.Contains("172v5_W") && it=="_cccentral") continue; //remove down_UP 0.755 FIXME
  //if(name.Contains("172v5_W") && it=="_ccentral") continue; //remove down_UP 0.755 FIXME
  //if(name.Contains("172v5_W") && it=="_central") continue; //remove down_UP 0.755 FIXME
  std::cout << "Running on tune: " << it << std::endl;
  TH1F *shiftData;
  float chi = chi2_jpsi_(shiftData, it, name, param[pos], syst, postfix);
  std:cout << "chi=" << chi << std::endl;
  if(abs(pre - chi) > 5 && abs(chiTest->FindBin(param[pos-2]) - chi) > 5 && pos>0) errtxt += Form("submit(\"%s\", false, \"%s\", 1);, ", name.Data(), tune[pos-1].Data());
  //if(abs(pre - chi) > chi) errtxt += Form("%s %0.3f, ", name.Data(), param[pos-1]);
  pre = chi;
  if(chi<low) low = chi;
  if(chi>high) high = chi;
  chiTest->GetYaxis()->SetRangeUser(int(low)-1,int(high)+2);
  chiTest->SetBinContent(chiTest->FindBin(param[pos]),chi);
  std::cout << chiTest->GetBinContent(chiTest->FindBin(param[pos])) << std::endl;
  hepdata += TString::Format("%.3f %.1f\n", param[pos], chi);
  //chiTest->SetBinError(chiTest->FindBin(param[pos]),sqrt(1./N));
}

//chiTest->GetXaxis()->SetRangeUser(0.65,1.255);
chiTest->GetXaxis()->SetRangeUser(0.65,0.976);//1.055);
//chiTest->GetYaxis()->SetRangeUser(55,90);
chiTest->GetYaxis()->SetRangeUser(int(low)-0.5,int(high)+1.5);
//chiTest->GetYaxis()->SetRangeUser(200,220);
chiTest->SetMarkerStyle(20);
chiTest->Draw("p9");
std::cout << chiTest->GetName() << std::endl;
std::cout << chiTest->GetTitle() << std::endl;
tdr(chiTest,epoch, fin);
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
((TF1*)(gROOT->GetFunction("pol3")))->SetParameters(-88.3245, 1045.87, -2049.8, 1137.63);
TString pol = "pol3";
 //TFitResultPtr fit = chiTest->Fit("pol3","FSEMQ","",0.6,1.055);
chiTest->Fit(pol,"FSMEQRW","",0.6,1.055);//0.976);
//TFitResultPtr fit = chiTest->Fit("pol3","FSEMQ","",0.6,0.975);
//TFitResultPtr fit = chiTest->Fit("pol2","FSMEQ");
//TFitResultPtr fit = chiTest->Fit("pol2","FSMEQ","",0.8,1.0);
/*
float min = (-1)*fit->Parameter(1)/(2*fit->Parameter(2));
float chimin = fit->Parameter(0) + fit->Parameter(1)*min + fit->Parameter(2) * pow(min,2);
float err = (-1)*fit->Parameter(1) / (2 * fit->Parameter(2)) - sqrt(pow(fit->Parameter(1),2)
            - 4 * fit->Parameter(2) * (fit->Parameter(0) - chimin - 1)) / (2 * fit->Parameter(2));
*/
float min = chiTest->GetFunction(pol)->GetMinimumX(0.6,1.075);
//float chimin = fit->Parameter(0) + fit->Parameter(1)*min + fit->Parameter(2) * pow(min,2) + fit->Parameter(3) * pow(min,3);
//std::cout << std::endl << chiTest->GetFunction("pol3")->GetParameter(0) << std::endl << chiTest->GetFunction("pol3")->GetParameter(1) << std::endl << chiTest->GetFunction("pol3")->GetParameter(2) << std::endl << chiTest->GetFunction("pol3")->GetParameter(3) << std::endl << std::endl; 
float chimin = chiTest->GetFunction(pol)->Eval(min);
float err = chiTest->GetFunction(pol)->GetX(chimin+1,0.6,1.075);
chiTest->GetYaxis()->SetRangeUser(int(std::min(low,chimin))-0.5,int(high)+1.5);
if(name=="") { nom=min; nerr=err; }
report = Form("Minimum at x= %g +/- %0.6g",min, abs(min-err));
json += Form("%.3f +/- %.3f ",min,abs(min-err));
//std::cout << "Minimum at x= " << min << " +/- " << abs(min - err) << std::endl;
std::cout << report << std::endl;
std::cout << "chi^2 at min= " << chimin << " " << TMath::Prob(chimin,10) << std::endl;
std::cout << "chi^2_min + 1 at x= " << err << std::endl;

TPaveText *pt = new TPaveText(0.12,0.85,0.3,0.65,"NDC"); //NB blNDC
pt->SetFillStyle(0);
pt->SetTextAlign(11);
pt->SetBorderSize(0);
pt->SetTextFont(42);
pt->SetTextSize(0.046);
TString text = TString::Format("r_{b}= %.4f +/- %.4f (stat)",min,abs(min-err));
if(name.Length() > 0)
  text += TString::Format(" %c %.4f (syst) +/- %.4f",(min<nom ? '-' : '+'), abs(nom-min), sqrt(abs(pow(nerr,2)-pow(abs(min-err),2))));
  //text += TString::Format(" %c %.4f (syst) +/- %.4f",(min<nom ? '-' : '+'), abs(nom-min), sqrt(abs(pow(0.0507584,2)-pow(abs(min-err),2))));
  //text += TString::Format(" %c %.4f (syst) +/- %.4f",(min<0.818905 ? '-' : '+'), abs(0.818905-min), sqrt(abs(pow(0.0507584,2)-pow(abs(min-err),2))));
pt->AddText(text);
if(!TDR) pt->Draw();
gStyle->SetOptStat(0);
TPaveText *var = new TPaveText(0.44,0.53,0.59,0.65,"NDC"); //NB blNDC
var->SetFillStyle(0);
var->SetTextAlign(21);
var->SetBorderSize(0);
var->SetTextFont(42);
var->SetTextSize(0.046);
TString tvar = TString::Format("J/#psi channel");
var->AddText(tvar);
tvar = TString::Format("(#it{ndf}=%d)", int(bin.size()/2) + 1);
var->AddText(tvar);
var->Draw();

if(name.Length()>0) name = "_" + name;
name += epoch_name[epoch+1];
if(fullpt) name += "_jpT";
if(postfix.Length()>0) postfix = "_" + postfix;
c1->SaveAs("chi2_jpsi"+name+(fin ? "_final" : "")+(postfix.Length()>0 ? postfix : "")+".pdf");
c1->SaveAs("chi2_jpsi"+name+(fin ? "_final" : "")+(postfix.Length()>0 ? postfix : "")+".png");

delete pt;
//chiTest->Delete();
//delete chiTest;
//delete c1;
}

float chi2_jpsi_(TH1F *&shiftData, TString tune="", TString name="", float num=0.855, int syst = 0, TString postfix="") {
if(postfix.Length()>0) postfix = "_" + postfix;
TH1F *data, *data2, *mc, *mc2;
TH1F *bkg1, *bkg2;
if(epoch!=0) { //FIXME !=
getHist(name, tune, data, mc, bkg1, epoch, true, syst);
}
else {
getHist(name, tune, data, mc, bkg1, 1, false, syst);
getHist(name, tune, data2, mc2, bkg2, 2, false, syst);
data->Add(data2);
mc->Add(mc2);
/*
bkg1->Add(bkg2);
data->Add(bkg1, -1);
*/
//std::vector<float> bin_jpsi = {0-0.025, 0.2-0.025, 0.4-0.025, 0.55-0.025, 0.6-0.025, 0.65-0.025, 0.7-0.025, 0.75-0.025, 0.8125-0.025, 0.88-0.025, 0.9425-0.025, 1.0};
//std::vector<double> rebin = {0-0.025, 0.2-0.025, 0.4-0.025, 0.6-0.025, 0.65-0.025, 0.7-0.025, 0.75-0.025, 0.9425-0.025, 1.0};
/*
std::vector<double> rebin = {0-0.025, 0.4-0.025, 0.6-0.025, 0.65-0.025, 0.7-0.025, 0.8125-0.025, 0.9425-0.025, 1.0};
data = (TH1F*)data->Rebin(rebin.size()-1, "", rebin.data());
mc = (TH1F*)mc->Rebin(rebin.size()-1, "", rebin.data());
*/
/*
for(int i = 1; i <= data->GetNbinsX(); i++) {
  data->SetBinContent(i, data->GetBinContent(i) / data->GetBinWidth(i));
  data->SetBinError(i, data->GetBinError(i) / data->GetBinWidth(i));
  mc->SetBinContent(i, mc->GetBinContent(i) / mc->GetBinWidth(i));
  mc->SetBinError(i, mc->GetBinError(i) / mc->GetBinWidth(i));
}
*/
if(name.Contains("toys") && tune == "_sdown") { //only compute modified hist once
  TH1F *dtmp, *dtmp2, *mtmp, *mtmp2;
  getHist(name, "", dtmp, mtmp, bkg1, 1, false, syst);
  getHist(name, "", dtmp2, mtmp2, bkg2, 2, false, syst);
  mtmp->Add(mtmp2);
  data = (TH1F*)mtmp->Clone();
  delete dtmp, dtmp2, mtmp, mtmp2;
  data->SetTitle("frame_ptfrac_toyData_hist");
  std::cout << data->GetTitle() << std::endl;
  std::cout << "creating toy data DONE!" << std::endl;
  std::cout << data->GetTitle() << std::endl;
  //scale hist to ensure the number of events hasn't changed
  //Clone toy data for shifting in each for loop (one shift per for-loop iteration)
  shiftData = (TH1F*)data->Clone("");
  //TH1F *shiftData = (TH1F*)data->Clone("");
  shiftData->SetDirectory(0);

  //vary MC bin content by data uncertainty
  TRandom3 *rand = new TRandom3(0);
  for(int i = 1; i <= data->GetNbinsX(); i++) {
    float y = shiftData->GetBinContent(i);
    float e(0.);
    e = shiftData->GetBinError(i);
    std::cout << y << " +/- " << e << " ";
    //shift each bin by random amount samples by the appropriate bin error
    if(!(name == "" || name == "172v5")) //MC errors otherwise (e.g. FSR)
      y += rand->Gaus(0, e);
    std::cout << y << std::endl;
    shiftData->SetBinContent(i,y);
    //if(name.Length()==0) shiftData->SetBinError(i,shiftData->GetBinError(i) * sqrt(float(nmc)/float(nentries)));
    shiftData->SetBinError(i,e);
  }
  delete rand;
  std::cout << data->GetEntries() << " " << shiftData->GetEntries() << std::endl;
  //shiftData->Scale(data->Integral()/shiftData->Integral());
  data = (TH1F*)shiftData->Clone();
  data->Print();
  shiftData->Print();
}
else if(name.Contains("toys"))
  data = (TH1F*)shiftData->Clone();
rbtest[2] += "{";
for(int ibin = 1; ibin <= mc->GetNbinsX(); ibin++) {
  rbtest[0] += TString::Format("%0.3f, ", num);
  rbtest[1] += TString::Format("%d, ", ibin);
  rbtest[2] += TString::Format("%0.2f, ", mc->GetBinContent(ibin));
  //rbtest[2] += TString::Format("{%0.2f, %0.3f}, ", mc->GetBinContent(ibin), mc->GetBinError(ibin));
}
rbtest[2] += "}";
delete data2;
delete mc2;
}
mc->GetXaxis()->SetTitle("J/#psi #it{p}_{T} / #Sigma #it{p}_{T}^{ch}");
data->GetXaxis()->SetTitle("J/#psi #it{p}_{T} / #Sigma #it{p}_{T}^{ch}");
setupPad()->cd();
tdr(mc, epoch, fin);
if(fullpt) mc->GetXaxis()->SetTitle("J/#psi #it{p}_{T}/ jet #it{p_}{T}");
mc->GetYaxis()->SetRangeUser(0,2000);
mc->Draw();
tdr(mc, epoch, fin);
//if(epoch>=0) {
gStyle->SetOptStat(0);
TString namet(name);
data->SetMarkerStyle(20);
data->SetMarkerColor(kBlack);
data->SetLineColor(kBlack);
data->SetLineWidth(2);
if(num==0) num=0.855;
if(namet == "") namet = "172v5";
//if(tunet == "") tunet = "855";
c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_signal_%s_%d%s_jpsi%s%s.pdf",namet.Data(),int(num*1000), epoch_name[epoch+1].Data(), (fullpt ? "_jpT" : ""), postfix.Data()));
c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_signal_%s_%d%s_jpsi%s%s.png",namet.Data(),int(num*1000), epoch_name[epoch+1].Data(), (fullpt ? "_jpT" : ""), postfix.Data()));

std::cout << "" << std::endl;
if(namet=="172v5" && num > 0.825 && num < 0.875) {
data->SetTitle("");
data->GetXaxis()->SetRangeUser(0,1.);
//data->GetYaxis()->SetRangeUser(0,0.145);
mc->SetMarkerStyle(20);
data->SetMarkerStyle(20);
data->SetMarkerColor(kBlack);
data->SetLineColor(kBlack);
data->SetLineWidth(2);
tdr(data, epoch, fin);
data->Draw();
tdr(data, epoch, fin);
//c1->SaveAs("www/meson/morph/ptfrac/ptfrac_signal_Data_BCDEFGH_jpsi.pdf");
//c1->SaveAs("www/meson/morph/ptfrac/ptfrac_signal_Data_BCDEFGH_jpsi.png");
c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_signal_Data%s_jpsi%s.pdf", epoch_name[epoch+1].Data(), (fullpt ? "_jpT" : "")));
c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_signal_Data%s_jpsi%s.png", epoch_name[epoch+1].Data(), (fullpt ? "_jpT" : "")));
}
//}

/*
if(tune=="" && name=="") {
TCanvas *c1 = setupCanvas();
TPad *p1 = setupPad();
p1->cd();
data->Draw();
gStyle->SetOptStat(0);
tdr(data,0, fin);
data->SetMarkerStyle(20);
data->SetMarkerColor(kBlack);
data->SetLineColor(kBlack);
data->SetLineWidth(2);
c1->SaveAs("ptfrac_signal_Data_"+name+"jpsi.pdf");
c1->SaveAs("ptfrac_signal_Data_"+name+"jpsi.png");

}
*/

N = mc->Integral();
mc->Scale(1./mc->Integral());
data->Scale(1./data->Integral());

data->GetXaxis()->SetRangeUser(0.6,1.);//0.975);
mc->GetXaxis()->SetRangeUser(0.6,1.);//0.975);
if(epoch<0 && 0) {
data->GetXaxis()->SetRangeUser(0.5,1.0);
mc->GetXaxis()->SetRangeUser(0.5,1.0);
}
if(fullpt) {
data->GetXaxis()->SetRangeUser(0.0,0.7);
mc->GetXaxis()->SetRangeUser(0.0,0.7);
}
data->SetLineColor(kBlack);
data->SetMarkerColor(kBlack);
data->SetMarkerStyle(20);
data->SetLineWidth(2);
mc->SetLineColor(kRed);
mc->SetMarkerColor(kRed);
mc->SetMarkerStyle(1);
mc->SetLineWidth(1);
mc->GetYaxis()->SetRangeUser(0.,.2);
data->GetYaxis()->SetRangeUser(0.,.2);
if(fullpt) {
mc->GetYaxis()->SetRangeUser(0.,.17);
data->GetYaxis()->SetRangeUser(0.,.17);
}
/*
int bin(mc->FindBin(0.425));
mc->SetBinContent(bin, 0);
data->SetBinContent(bin, 0);
*/
mc->Draw("hist");
tdr(mc, epoch, fin);
mc->Draw("same e");
data->Draw("same");
if(num==0) num=0.855;
if(name=="") name="172v5";
TPaveText *txt = new TPaveText(0.3,0.90,0.4,0.70,"NDC"); //NB blNDC
txt->SetFillStyle(0);
txt->SetTextAlign(11);
txt->SetBorderSize(0);
txt->SetTextFont(42);
txt->SetTextSize(0.046);
TString text = TString::Format("#it{r}_{B} = %.3f", num);
txt->AddText(text);
txt->Draw();
TString mcvname(TString::Format("mcVdata_%s_%d_jpsi",name.Data(),(int)(num*1000)) + epoch_name[epoch+1]);
if(fullpt) mcvname += "_jpT";
c1->SaveAs(mcvname + (postfix.Length()>0 ? postfix : "") + ".pdf");
c1->SaveAs(mcvname + (postfix.Length()>0 ? postfix : "") + ".png");
//data->SetBinContent(data->FindBin(0.4), 0);
//mc->SetBinContent(data->FindBin(0.4), 0);
float chi2 = data->Chi2Test(mc, "CHI2 P WW");
/*
chi2 = 0.;
float sum1(0.);
float sum2(0.);
for(int i = 0; i < data->GetNbinsX(); i++) {
  float exp = mc->GetBinContent(i);
  float obs = data->GetBinContent(i);
  sum1 += obs;
  sum2 += exp;
  
}
std::cout << sum1 << " " << sum2 << std::endl;
float ndf = data->GetXaxis()->GetLast() - data->GetXaxis()->GetFirst();
for(int i = 0; i < data->GetNbinsX(); i++) {
  float cnt1 = data->GetBinContent(i);
  float cnt2 = mc->GetBinContent(i);
  float e1sq = pow(data->GetBinError(i),2);
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

delete data;
delete mc;

return chi2;

}
