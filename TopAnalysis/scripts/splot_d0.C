//
// Execute the macro with
// root -l roofit_d0.C++
//

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLegend.h"
#include "TMath.h"
#include <RooFit.h>
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooStats/SPlot.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooGamma.h"
#include "RooGlobalFunc.h"
#include <RooFitResult.h>
#include "RooWorkspace.h"
#include "RooMinuit.h"
#include <vector>
#include "/afs/cern.ch/user/b/byates/TopAnalysis/interface/CharmEvent.h"
#include "convert.h"
#ifndef CHARM
#define CHARM
#include "/afs/cern.ch/user/b/byates/TopAnalysis/src/CharmEvent.cc"
#endif
//#include "TopAnalysis/interface/CharmEvent.h"
using namespace RooFit;
using namespace RooStats;

bool GET_BIT(short x, short b) { return (x & (1<<b) ); }

//void roofit_mtop_BCDEFGH(TString mass="166v5", TString file="d0_fit.root") {
//void mtop_norm(std::vector<pair<float,float>> &p, TString mass="171.5", short flags=0b00) {
void splot_d0(TH1F *&ptfrac_signal, TString mass="172.5", bool isData=false, TString fragWeight="", int ep=0, bool jpT=false, int xb=0) {
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/param.h"
  RooWorkspace w("w",mass);
  float wind(0.045);
  //std::vector<float> bin;
  //bin = {0.025, 0.1, 0.175, 0.25, 0.325, 0.4, 0.475, 0.55, 0.625, 0.7, 0.775, 0.85, 0.925, 0.975, 1.0};
  //bin = {0, 0.075, 0.15, 0.225, 0.3, 0.375, 0.45, 0.525, 0.6, 0.675, 0.75, 0.825, 0.9, 0.975, 1.0};
  //bin = {0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.875, 0.95, 1.0};
  //bin = {-0.025, 0.05, 0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.875, 0.95, 1.0};
  //std::vector<float> type = {0., 0.45, 0.525, 0.6, 0.675, 0.75, 0.825, 0.9, 1.};
  //std::vector<float> type = {0., 0.425, 0.575, 0.65, 0.725, 0.8, 0.875, 0.95, 1.};
  //std::vector<float> type = {0., 0.425, 0.65, 0.8, 1.};
  //std::vector<float> type = {0., 0.375, 0.6, 0.825, 1.};
  //std::vector<float> type = {0., 0.4, 0.6, 0.8, 1.};
  //std::vector<float> type = {0., 0.4, 0.6, 1.};
  //std::vector<float> type = {0., 0.3, 0.45, 0.8, 1.};
  //std::vector<float> type = {0., 0.25, 0.5, 0.7, 1.};
  TString cut("");
  if(xb != 0)
    cut = TString::Format("ptfrac > %f && ptfrac < %f", type[xb-1], type[xb]);
  /*
  bool tpt(true);
  int pdfs[4] = {0, 1, 110, 111}; //nominal, a_s=0.118, 0.117, 0.119
  int pdf(0);
  */
/*
void splot_d0_mu(RooWorkspace &w, TString mass="172.5", bool isData=false) {
*/
  //TFile *f = new TFile("../BatchJobs/merged.root"); 
  //TFile *f = new TFile("plots/plotter_mtop_BCDEFGH.root");
  TString syst("");
  //TString dir("/eos/cms/store/user/byates/top18012/Chunks/");
  TString dir("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/test/Chunks/");
  //TString dir("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/etaPiK/Chunks/");
  mass.ReplaceAll(".","v");
  if(mass.Contains("v5") && !mass.Contains("172v5")) mass = "m" + mass;
  if(mass.Contains("sr") || mass.Contains("erdON") || mass.Contains("Move") || mass.Contains("ue") || mass.Contains("hdamp")) { //ISR,FSR,CR,UE
    syst = mass;
    std::cout << "Processing systematics " << TString::Format("MC13TeV_TTJets_%s",syst.Data()) << std::endl;
  }
  else if(mass.Contains("m1")) { //top mass tunes
    syst = mass;
    //mass = "172.5";
    std::cout << "Processing systematics " << TString::Format("MC13TeV_%s",syst.Data()) << std::endl;
  }
  /*
  else if(mass.Contains("tpt")) { //no top pt reweighting
    tpt = false;
  }
  else if(mass.Contains("as118")) {
    pdf = 1;
  }
  else if(mass.Contains("as117")) {
    pdf = 2;
  }
  else if(mass.Contains("as119")) {
    pdf = 3;
  }
  */
  TFile *fin;
  TGraph *g=nullptr;
  //std::vector<TGraph*> fwgt;
  TH1F *tuneWgt = new TH1F("tuneWgt","tuneWgt",2,0,2);
  //TString dir("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/ndau/Chunks/");
  TFile *f = new TFile(dir+"../MC13TeV_TTJets_powheg.root"); //open a randome file to get correction histograms
  //TFile *f = new TFile(dir+"../MC13TeV_TTJets_powheg.root"); //open a randome file to get correction histograms
  TChain *data = new TChain("data");
  /*
  if(isData) {
    data->Add(dir+"Data13TeV_*");
    mass="Data";
  }
  */
  if(isData) {
    std::vector<TString> dataSamples;
    if(ep == 1)
      dataSamples = {"Data13TeV_DoubleEG_2016B","Data13TeV_DoubleEG_2016C","Data13TeV_DoubleEG_2016D","Data13TeV_DoubleEG_2016E","Data13TeV_DoubleEG_2016F","Data13TeV_DoubleMuon_2016B","Data13TeV_DoubleMuon_2016C","Data13TeV_DoubleMuon_2016D","Data13TeV_DoubleMuon_2016E","Data13TeV_DoubleMuon_2016F","Data13TeV_MuonEG_2016B","Data13TeV_MuonEG_2016C","Data13TeV_MuonEG_2016D","Data13TeV_MuonEG_2016E","Data13TeV_MuonEG_2016F","Data13TeV_SingleElectron_2016B","Data13TeV_SingleElectron_2016C","Data13TeV_SingleElectron_2016D","Data13TeV_SingleElectron_2016E","Data13TeV_SingleElectron_2016F","Data13TeV_SingleMuon_2016B","Data13TeV_SingleMuon_2016C","Data13TeV_SingleMuon_2016D","Data13TeV_SingleMuon_2016E","Data13TeV_SingleMuon_2016F"};
    else if(ep == 2)
      dataSamples = {"Data13TeV_SingleMuon_2016G","Data13TeV_SingleMuon_2016H_v2","Data13TeV_SingleMuon_2016H_v3","Data13TeV_DoubleEG_2016G","Data13TeV_DoubleEG_2016H_v2","Data13TeV_DoubleEG_2016H_v3","Data13TeV_DoubleMuon_2016G","Data13TeV_DoubleMuon_2016H_v2","Data13TeV_DoubleMuon_2016H_v3","Data13TeV_MuonEG_2016G","Data13TeV_MuonEG_2016H_v2","Data13TeV_MuonEG_2016H_v3","Data13TeV_SingleElectron_2016G","Data13TeV_SingleElectron_2016H_v2","Data13TeV_SingleElectron_2016H_v3"};
    else
      dataSamples = {"Data13TeV_DoubleEG_2016B","Data13TeV_DoubleEG_2016C","Data13TeV_DoubleEG_2016D","Data13TeV_DoubleEG_2016E","Data13TeV_DoubleEG_2016F","Data13TeV_DoubleMuon_2016B","Data13TeV_DoubleMuon_2016C","Data13TeV_DoubleMuon_2016D","Data13TeV_DoubleMuon_2016E","Data13TeV_DoubleMuon_2016F","Data13TeV_MuonEG_2016B","Data13TeV_MuonEG_2016C","Data13TeV_MuonEG_2016D","Data13TeV_MuonEG_2016E","Data13TeV_MuonEG_2016F","Data13TeV_SingleElectron_2016B","Data13TeV_SingleElectron_2016C","Data13TeV_SingleElectron_2016D","Data13TeV_SingleElectron_2016E","Data13TeV_SingleElectron_2016F","Data13TeV_SingleMuon_2016B","Data13TeV_SingleMuon_2016C","Data13TeV_SingleMuon_2016D","Data13TeV_SingleMuon_2016E","Data13TeV_SingleMuon_2016F","Data13TeV_SingleMuon_2016G","Data13TeV_SingleMuon_2016H_v2","Data13TeV_SingleMuon_2016H_v3","Data13TeV_DoubleEG_2016G","Data13TeV_DoubleEG_2016H_v2","Data13TeV_DoubleEG_2016H_v3","Data13TeV_DoubleMuon_2016G","Data13TeV_DoubleMuon_2016H_v2","Data13TeV_DoubleMuon_2016H_v3","Data13TeV_MuonEG_2016G","Data13TeV_MuonEG_2016H_v2","Data13TeV_MuonEG_2016H_v3","Data13TeV_SingleElectron_2016G","Data13TeV_SingleElectron_2016H_v2","Data13TeV_SingleElectron_2016H_v3"};
    for(auto & it : dataSamples) {
      TString dataname(dir);
      dataname += it;
      dataname += "_*.root";
      std::cout << "Adding " << it;
      int num = data->Add(dataname);
      std::cout << " (" << num << " files)" << std::endl;
    }
  }
  //else data->Add("Chunks/MC13TeV_TTJets_m"+mass+"_*.root");
  else {
    if(mass.Contains("TRK") || mass.Contains("TRIGGER") || mass.Contains("LEP") ||mass.Contains("PU") || mass.Contains("PI") || mass.Contains("JER") || mass.Contains("JSF")) {  // SFs
      TString tmp(mass);
      tmp.ReplaceAll("_","/");
      //dir = TString("/eos/cms/store/user/byates/top18012/" + tmp + "/Chunks/");
      dir = TString("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/test/" + tmp + "/");
      //dir = TString("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/test/" + tmp + "/Chunks/");
      //if(mass.Contains("JSF")) dir = TString("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/test/" + tmp + "/Chunks/");
      std::cout << dir << std::endl;
    }
    else if(mass.Contains("_noHT")) {
      dir = TString("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/noHT/Chunks/");
      std::cout << dir << std::endl;

    }
    //std::vector<TString> mcSamples = { "MC13TeV_TTJets_powheg" };
    std::vector<TString> mcSamples;
    //PS weights instead of deficated FSR samples
    if(syst.Length()>0 && !mass.Contains("FSR")) mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW",TString::Format("MC13TeV_TTJets_%s*",syst.Data()),"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    else if(mass.Contains("FSR-up")) mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW",TString::Format("MC13TeV_TTJets_fsr-up"),"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    else if(mass.Contains("FSR-down")) mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW",TString::Format("MC13TeV_TTJets_fsr-down"),"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    else mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_TTJets_powheg","MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    //mcSamples = { "MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets" }; // For W background only
    if(mass.Contains("QCD") || mass.Contains("GluonMove")) mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW",TString::Format("MC13TeV_%s*",syst.Data()),"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    for(auto & it : mcSamples) {
      //TString mcname = "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/Chunks/";
      TString mcname(dir);
      mcname += it;
      mcname += "_*.root";
      std::cout << "Adding " << it;
      int num = data->Add(mcname);
      std::cout << " (" << num << " files)" << std::endl;
    }
    if(fragWeight.Length() > 0) {
      if(fragWeight.Contains("lep")) {
        fin = TFile::Open("/eos/cms/store/user/byates/top18012/bfragweights_Markus.root");
      }
      else fin = TFile::Open("/eos/cms/store/user/byates/top18012/bfragweights.root");
      mass += "_" + fragWeight;
      fragWeight.ReplaceAll("lep","");
      g = (TGraph*)fin->Get(fragWeight+"Frag");
      std::cout << g->GetName() << std::endl;
    }
  }
  TString fUrl("/eos/cms/store/group/phys_top/byates/sPlot//TopMass_"+mass+"_sPlot_d0.root");
  if(ep>0) fUrl.ReplaceAll(".root",TString::Format("%d.root",ep));
  if(jpT) fUrl.ReplaceAll(".root","_jpT.root");
  if(xb != 0) fUrl.ReplaceAll(".root",TString::Format("_%d.root",int(type[xb]*1000)));
  std::cout << "creating file: "  << fUrl<< std::endl;
  TFile *fout;// = new TFile(fUrl,"RECREATE");
  if(!mass.Contains("toyData"))
  fout = new TFile(fUrl,"RECREATE");
  //data->Add("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/Chunks/MC13TeV_W*Jets_*.root");
  //disabled top pT reweighting
  //data->Draw("d0_mass>>h(30,1.7,2.0)","norm*sfs*puwgt*(meson_id==42113 && d0_l3d/d0_sigmal3d>10 && HT>180 && d0_sigmal3d>2E-4)");// && j_hadflav[d0_j]==5 && d0_pi_mother==421 && d0_k_mother==421)","goff");
  //data->Draw("d0_mass>>h(30,1.7,2.0)","norm*sfs*puwgt*topptwgt*(meson_id==42113 && d0_l3d/d0_sigmal3d>10 && HT>180 && d0_sigmal3d>2E-4)");// && j_hadflav[d0_j]==5 && d0_pi_mother==421 && d0_k_mother==421)","goff");
  
  CharmEvent_t ev;
  attachToCharmEventTree(data,ev);
  TH1F *h = (TH1F*)gDirectory->Get("h");
  //if(!isData) h->Scale(35864);
  //if(!isData) h->Scale(832*35864);
  //f->ls(); 
  TString name = "massD0_l_all_d0";
  //TString name = "massJPsi_l_all_d0_BCDEF";
  TH1F *pu1 = (TH1F*) f->Get("puwgtctr_BCDEF");
  TH1F *pu2 = (TH1F*) f->Get("puwgtctr_GH");
  TH1F *top1 = (TH1F*) f->Get("topptwgt_BCDEF");
  TH1F *top2 = (TH1F*) f->Get("topptwgt_GH");
  float puSF1 = pu1->GetBinContent(1)/pu1->GetBinContent(2);
  float puSF2 = pu2->GetBinContent(1)/pu2->GetBinContent(2);
  float topSF1 = top1->GetBinContent(2)/top1->GetBinContent(1);
  float topSF2 = top2->GetBinContent(2)/top2->GetBinContent(1);
  float jsfSF1(1.);
  float jsfSF2(1.);
  if(mass.Contains("JSF") && 0) {
    TH1F *jsf1 = (TH1F*) f->Get("jsfwgt_BCDEF");
    TH1F *jsf2 = (TH1F*) f->Get("jsfwgt_GH");
    jsfSF1 = jsf1->GetBinContent(1)/jsf1->GetBinContent(2);
    jsfSF2 = jsf2->GetBinContent(1)/jsf2->GetBinContent(2);
  }
  pu1->SetDirectory(0);
  pu2->SetDirectory(0);
  top1->SetDirectory(0);
  top2->SetDirectory(0);
  f->Close();
  if(!isData) {
  cout << "PU normalization " << puSF1 << endl;
  cout << "top pT normalization " << topSF1 << endl;
  cout << "PU normalization " << puSF2 << endl;
  cout << "top pT normalization " << topSF2 << endl;
  }
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  if(xb!=0) std::cout << "slecting " << type[xb-1] << " < ptfrac < " << type[xb] << std::endl;
  cout << "loaded!" << endl;

  // Declare observable x
  RooRealVar d0_mass("d0_mass","D^{0} mass", 1.7, 2, "GeV") ;
  RooRealVar d0_mass2("d0_mass2","D^{0} mass2", 1.7, 2, "GeV") ;
  RooRealVar weight("weight","weight",1.,0.,2.);
  //RooRealVar weight("weight","weight",1.,0.,36000.);
  RooRealVar meson_l_mass("d0_l_mass","D^{0}+l mass", 0, 250, "GeV") ;
  RooRealVar ptfrac("ptfrac","D^{0} #it{p}_{T} / #Sigma #it{p}_{T}^{ch}", 0, 1.1, "") ;
  RooRealVar d0_pt("d0_pt","D^{0} p_{T}", 0, 250, "GeV");
  RooRealVar j_pt_ch("j_pt_ch","j p_{T} charged", 0, 400, "GeV");
  RooRealVar j_pt("j_pt","j p_{T}", 0, 400, "GeV");
  RooRealVar epoch("epoch","epoch",1,2);
  RooRealVar tuneW = RooRealVar("tuneW", "tuneW", 1., 0, 2.);
  
  //cout << "creating dataset" << endl;
  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'

  //RooDataSet dsn("dsn", "dsn", RooArgSet(meson_id,d0_mass,ptfrac,meson_l_mass,weight), Import(*data), Cut("meson_id==42113"));
  RooDataSet dsn("dsn", "dsn", RooArgSet(RooArgList(d0_mass,ptfrac,meson_l_mass,weight,d0_pt,epoch,j_pt_ch,j_pt,tuneW,""),RooArgList(d0_mass2,"")));
  //RooDataSet dsn("dsn", "dsn", args);
  int nmc = ep == 1 ? 430267 : 555621;
  int nentries = ep == 1 ? 97791 : 92135;
  int nset(0);
  std::cout << "Total events: " << nmc  << std::endl;
  /*
  int start = nmc;
  if(mass.Contains("toyData")) {
    std::cout << "Creating toy dataset" << std::endl;
    TRandom3 *rand = new TRandom3(0);
    //pick a random starting point that will encompass a subet of MC containing N points where N=num of points in the data
    std::cout << nmc << " " << nentries << " " <<  int(float(nmc)/float(nentries)) << std::endl;
    while(start + nentries > nmc) {
      start = rand->Uniform(0, nentries);
      std::cout << start << " " << start + nentries << " " << nentries << std::endl;
    }
    delete rand;
  }
  */
  for(int i=0; i < data->GetEntries(); i++) {
    ev = {};
    data->GetEntry(i);
    //scale *= float(nmc)/float(nentries);
    for(int j=0; j<ev.nmeson; j++) {
      float scale = 1.;
      //int j(0);
      float pi = 0.;
      float k = 0.;
      /* FIXME
      if(ev.nmeson>1) {
        //if(pi == ev.d0_k_pt[j] && k == ev.d0_pi_pt[j]) continue;
        if(ev.nmeson > j+2) //entries are 2 apart (j=0 epoch 1 j=1 epoch 2, j=2 epoch 1 j=3 epoch 2
          if(ev.d0_pi_pt[j] == ev.d0_k_pt[j+2] && ev.d0_k_pt[j] == ev.d0_pi_pt[j+2] && ev.meson_id[j] == ev.meson_id[j+2]) { scale *= 0.5; }
        if(ev.nmeson > j+1) //entries are 2 apart (j=0 epoch 1 j=1 epoch 2, j=2 epoch 1 j=3 epoch 2
          if(ev.d0_pi_pt[j] == ev.d0_k_pt[j+1] && ev.d0_k_pt[j] == ev.d0_pi_pt[j+1] && ev.meson_id[j] == ev.meson_id[j+1]) { scale *= 0.5; }
        pi = ev.d0_pi_pt[j];
        k = ev.d0_k_pt[j];
      }
      */
      //if(ev.j_pt_charged[j] == 0) continue; //some events are bad
      if(ev.meson_id[j] != 421) continue;
      if(ev.d0_mass[j] < 1.7) continue;
      if(ev.d0_mass[j] > 2.0) continue;
      if(ev.ht < 180) continue;
      if(ep>0 && ev.epoch[j]!=ep) continue;
      /*
      if(!jpT) { //missing in raw plots
        //if(ev.j_pt[j]>150) continue;
        if(ev.j_pt_charged[j]>100) continue;
      }
      */
      if(mass.Contains("toyData")) {
        TRandom3 *rand = new TRandom3(0);
        if(rand->Uniform(0, 1) > float(nentries)/float(nmc)) continue;
        if(nset > nentries) continue;
        nset++;
      }
      /*
      else
        if(ev.j_pt[j]>200) continue;
      */
      /*
      if(ev.d0_l3d[j]/ev.d0_sigmal3d[j]<10) continue;
      if(ev.d0_sigmal3d[j]<2E-4) continue;
      if(ev.ht<180) continue;
      if(ev.j_pt[j]>150) continue;
      */
      //if(ev.j_pt_charged[j]>75) continue;
      //if(ev.epoch[j]!=1 && ev.epoch[j]!=2) continue;
      //remove events below xB=0.2
      if(ev.d0_pt[j]/ev.j_pt_charged[j] < 0.2) continue;
      if(xb != 0) {
         /*
	if(type[xb] < 1) { cut = TString::Format("ptfrac > %f && ptfrac < %f",type[xb-1], type[xb]); }
	else { cut = TString::Format("ptfrac > %f",type[xb-1]); }
         */
	if(type[xb] < 1) { if(ev.d0_pt[j]/ev.j_pt_charged[j] > type[xb] || ev.d0_pt[j]/ev.j_pt_charged[j] < type[xb-1]) continue; }
        else if(ev.d0_pt[j]/ev.j_pt_charged[j] < type[xb-1]) continue; 
      }
      /*
      if(!ev.isData) {
        if(ev.j_hadflav[ev.d0_j[j]]!=5) continue;
        if(ev.d0_pi_mother[j]!=421) continue;
        if(ev.d0_k_mother[j]!=421) continue;
      }
      */
      //if(j_pt[(int)d0_j[j]]<30) continue;
      /*
      if(!(mesonlm[j] > 0)) continue;
      if(mesonlm[j] > 250) continue;
      */
      tuneW.setVal(1.);
      scale *= ev.sfs[j]; //Data has sfs=0.5 if duplicates (piK and Kpi)

      //Observed small differences between historams and treess
      if(ev.epoch[j]==1)
        scale *= 1.06; //BCDEF normalization const
      if(ev.epoch[j]==2)
        scale *= 1.13; //GH normalization const

      if(!isData) {
        if(mass.Contains("toyData"))
        scale *= ev.puwgt[j] * ev.topptwgt;// * topSF * puSF;
        else
        scale *= ev.norm * ev.xsec * ev.puwgt[j] * ev.topptwgt;// * topSF * puSF;
        //if(pdf>0) scale *= pdfs[pdf];//alternate PDF event weight
        //if(tpt) scale *= ev.topptwgt;// * topSF * puSF;
        //scale *= ev.norm * ev.xsec * ev.sfs[j] * ev.puwgt[j] * ev.topptwgt;// * topSF * puSF;
        //PS weights
        /*
        if(mass.Contains("FSR-up")) scale *= ev.gfsr[0];
        else if(mass.Contains("FSR-down")) scale *= ev.gfsr[1];
        */
        //scale *= 1.11; //GH normalization const
        //if(!jpT) scale *= ev.pitrk[j];
        //scale = norm * sfs[j] * puwgt[j] * topptwgt * topSF * puSF;
        if(ev.epoch[j]==1) {
          if(mass.Contains("toyData"))
          scale =  scale * puSF1 * topSF1;
          else
          scale =  scale * 19712.86 * puSF1 * topSF1;
          //if(tpt) scale =  scale * topSF1;
          //scale *= 1.06; //GH normalization const
          //h1->Fill(mesonlm[j], scale);
        }
        else if(ev.epoch[j]==2) {
          if(mass.Contains("toyData"))
          scale = scale * puSF2 * topSF2;
          else
          scale = scale * 16146.178 * puSF2 * topSF2;
          //if(tpt) scale =  scale * topSF1;
          //if(tpt) scale =  scale * topSF2;
          //scale *= 1.13; //GH normalization const
          //h2->Fill(mesonlm[j], scale);
        }
        else
          continue;
        tuneWgt->Fill(0.,1.0);
        if(fragWeight.Length() > 0) {
          double frac(ev.d0_pt[j]/ev.j_pt_charged[j]);
          if(frac>1.) frac = 1.;
          if(frac<0.1) frac = 0.1;
          scale *= g->Eval(frac);
          tuneWgt->Fill(1.,g->Eval(frac));
        }
        else tuneWgt->Fill(1., 1.0); //unit weights to make re-weighting easier to loop over
      }
      epoch.setVal(ev.epoch[j]);
      //d0_mass = ev.d0_mass[j];
      d0_mass.setVal(ev.d0_mass[j]);
      d0_mass2.setVal(ev.d0_mass[j]);
      //if(ev.d0_mass[j]>(1.864-wind) && ev.d0_mass[j]<(1.864+wind) && ev.d0_pt[j]>0) { //Symmetric around PDG mass
      //if(ev.d0_mass[j]>(1.864-0.034) && ev.d0_mass[j]<(1.864+0.034) && ev.d0_pt[j]>0) { //Symmetric around PDG mass
      //if(ev.d0_mass[j]>1.83 && ev.d0_mass[j]<1.9 && ev.d0_pt[j]>0) {
      if(!jpT) ptfrac.setVal(ev.d0_pt[j]/ev.j_pt_charged[j]);
      else ptfrac.setVal(ev.d0_pt[j]/ev.j_pt[j]);
      d0_pt.setVal(ev.d0_pt[j]);
      j_pt_ch.setVal(ev.j_pt_charged[j]);
      j_pt.setVal(ev.j_pt[j]);
      //meson_l_mass = ev.d0_l_mass[j];
      meson_l_mass.setVal(ev.d0_l_mass[j]);
      //if(ev.d0_pt[j]/ev.j_pt_charged[j] * zwgt > 0.79 && ev.d0_pt[j]/ev.j_pt_charged[j] * zwgt < 0.81) std::cout << i << " " << ev.epoch[j] << " " << ev.d0_pt[j]/ev.j_pt_charged[j] * zwgt << " " << zwgt << " " << ev.d0_pt[j]/ev.j_pt_charged[j] <<  " " << ev.d0_pi_pt[j] << " " << ev.d0_pi_eta[j] << " " << ev.d0_k_pt[j] << " " << ev.d0_k_eta[j] << " " << scale << std::endl;
      //}
      weight.setVal(scale);
      tuneW.setVal(scale);// * tuneW.getVal());
      dsn.add(RooArgSet(RooArgList(d0_mass,meson_l_mass,ptfrac,weight,d0_pt,j_pt_ch,j_pt,epoch,tuneW,""),RooArgList(d0_mass2,"")));
      //dsn.add(RooArgSet(j_pt));
      //dsn.add(RooArgSet(d0_mass,meson_l_mass,ptfrac,weight), scale);
      //within D^0 mass peak (1.864 +/- 0.05)
      //if(ev.d0_l_mass[j] > 1.8 && ev.d0_l_mass[j] < 1.93)
        //dsn.add(RooArgSet(meson_l_mass,ptfrac), scale);

    }
  }
  delete data;
  std::cout << "Done loading tree" << std::endl;

  tuneW.setConstant();
  if(tuneWgt->GetBinContent(-2) > 0) tuneW.setVal(tuneWgt->GetBinContent(2)/tuneWgt->GetBinContent(1));
  /*
  RooDataSet tmp("tmp", "tmp", &dsn, *dsn.get(), 0, "weight");
  RooDataSet ds("ds", "ds", &tmp, *tmp.get(), 0, "tuneW");
  */
  double tuneShape = 1.;
  if(tuneWgt->GetBinContent(2) > 0) tuneShape = tuneWgt->GetBinContent(1)/tuneWgt->GetBinContent(2);
  std::cout << "Norm weight: " << tuneShape << std::endl;
  //very inefficient second loop to adjust normalization
  std::cout << "Re-weighting based on norm weight" << std::endl;
  for(int i = 0; i < dsn.numEntries(); i++) {
    double tmpW = dsn.get(i)->getRealValue("tuneW");
    tuneW.setVal(tmpW * tuneShape);
    weight.setVal(tuneW.getVal());
    //std::cout << tuneW.getVal() << " " << weight.getVal() << std::endl;;
  }
  std::cout << "Re-weighting done!" << std::endl;
  /*
  tmp.Print();
  tmp.addColumn(tuneW, tuneShape);
  tmp.Print();
  //RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, "weight");
  RooDataSet ds("ds", "ds", &tmp, *tmp.get(), 0, "tuneW");
  RooDataSet tmp("tmp", "tmp", &dsn, *dsn.get(), 0, "weight");
  tmp.add(tuneW, tuneShape);
  */
  //RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, "weight");
  //RooRealVar *wgt = (RooRealVar*)dsn.addColumn(tuneW, tuneShape);
  RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, "tuneW");
  //RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, "weight");
  //RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, wgt->GetName());
  /*
  for(int i = 0; i < dsn.numEntries(); i++) {
    double tmpW = dsn.get(i)->getRealValue("tuneW");
    double frac = dsn.get(i)->getRealValue("ptfrac");
    tuneW.setVal(tmpW * tuneShape);
    weight.setVal(tuneW.getVal());
    //std::cout << tuneW.getVal() << " " << weight.getVal() << std::endl;;
        tuneWgt->Fill(0.,1.0);
        if(fragWeight.Length() > 0) {
          if(frac>1.) frac = 1.;
          if(frac<0.1) frac = 0.1;
          double scale = g->Eval(frac);
          ptfrac.setVal(frac * scale);
          tuneWgt->Fill(1.,g->Eval(frac));
        }
  }
  */
  dsn.Print("v");
  ds.Print("v");
  ds.Print();
  std::cout << ds.weight() << std::endl;
  cout << "plotting dataset" << endl;
  // Make plot of binned dataset showing Poisson error bars (RooFit default)

  cout << "defining variables" << endl;

  // Mean of the J/psi mass peak
  /*
  RooRealVar mean("mean","mean", 1.864, 1.8, 1.93);

  // Construct Gaussian1 PDF for signal
  RooRealVar sigma("sigma","sigma", 0.02, 0.005, 0.06);
  RooRealVar ngsig("ngsig","ngsignal", 5000, 100, 10000000);
  RooGaussian gauss("gauss","gauss", d0_mass, mean, sigma);

  // Construct exponential PDF to fit the bkg component
  RooRealVar lambda("lambda", "slope", -3., -10, 0.);
  //RooRealVar lambda("lambda", "slope", -2.9, -10, 0.);
  RooExponential expo("expo", "exponential PDF", d0_mass, lambda);
  
  // Construct signal + bkg PDF
  RooRealVar nsig("nsig","#signal events", 1300, 0, 10000) ;
  RooRealVar nbkg("nbkg","#background events", 5000, 0, 10000) ;
  RooAddPdf model("model","g+exp", RooArgList(gauss, expo), RooArgList(nsig,nbkg)) ;
  */
  RooRealVar mean("mean","mean", 1.864, 1.864-wind, 1.864+wind);

  // Construct Crystal Ball PDF for signal
  RooRealVar cbmean("cbmean", "cbmean" , 1.864, 1.85, 1.87); 
  RooRealVar cbsigma("csigma", "cbsigma" , 0.02, 0.001, 0.25); 
  RooRealVar ncbsig("ncbsig", "ncbsignal", 4000, 0, 10000); 
  RooRealVar n("n","cbn", 5, 0, 10);
  RooRealVar alpha("#alpha","cbalpha", 1, 0, 5);
  RooCBShape cball("cball", "crystal ball", d0_mass, mean, cbsigma, alpha, n);

  //Gamma
  RooRealVar gam("gam","gamma", 2.5, 0, 10);
  RooRealVar beta("beta","beta", 1, 0, 40);
  RooRealVar mu("mu","mu", 1.864, 0, 50);

  RooGamma gamma("gamma","gamma",d0_mass,gam,beta,mu);

  //Polynomial
  RooRealVar p0("p0", "p0", 1, 0, 1000);
  RooRealVar p1("p1", "p1", 1, 0, 10);
  RooRealVar p2("p2", "p2", 1, 0, 10);
  RooPolynomial poly("poly", "poly", d0_mass, RooArgList(p0,p1,p2));


  // Construct Gaussian1 PDF for signal
  RooRealVar sigma("sigma","sigma", 0.0195, 0.005, 0.5);
  RooRealVar ngsig("ngsig","ngsignal", 1000, 100, 10000000);
  RooGaussian gauss("gauss","gauss", d0_mass, mean, sigma);
  RooRealVar sigma1("sigma1","sigma1", 0.06, 0.005, 0.5);
  RooRealVar ngsig1("ngsig1","ngsignal1", 1000, 100, 10000000);
  RooGaussian gauss1("gauss1","gauss1", d0_mass, mean, sigma1);
  RooPoisson pois("pois","pois", d0_mass, mean);

  // Gaussian for D0->KK
  float fNgsigkk[2] = {1.11557e+02, 8.20522e+02};
  float fSigmakk[2] = {1.94549e-02, 1.99998e-02};
  RooRealVar meankk("meankk","meankk", 1.793);//, 1.793-wind, 1.793+0.02);
  RooRealVar sigmakk("sigmakk","sigmakk", 2e-2, 0.0, 0.02);
  RooRealVar ngsigkk("ngsigkk","ngsignalkk", 100, 1, 10000);
  RooGaussian gausskk("gausskk","gausskk", d0_mass, meankk, sigmakk);

  // Gaussian for D0 from W
  // parameters split by xB
  std::vector< std::vector< std::vector<float> > > param_xb = {{{1.86574, 3.84928e+02, 1.71942e-02}, {1.86956, 3.15513e+01, 5.00000e-03}, {1.84927, 6.91239e+01, 2.39324e-02}, {1.86667, 1.14463e+02, 1.43698e-02}, {1.86901, 1.72984e+02, 1.83646e-02}}, {{1.86583, 3.30189e+02, 1.67071e-02}, {1.87131, 2.56132e+01, 5.00016e-03}, {1.85611, 5.57578e+01, 1.79126e-02}, {1.86782, 1.16715e+02, 1.45450e-02}, {1.86888, 1.31075e+02, 2.04610e-02}}};
  float fMeanW = param_xb[ep-1][xb][0];
  float fNgsigW = param_xb[ep-1][xb][1];
  float fSigmaW = param_xb[ep-1][xb][2];
  std::cout << fMeanW << "\t" << fNgsigW << "\t" << fSigmaW << std::endl;
  RooRealVar meanW("meanW","meanW", fMeanW);//1.86583, 1.793-wind, 1.793+wind);
  RooRealVar sigmaW("sigmaW","sigmaW", fSigmaW);//1.67071e-02, 0, 0.02);
  RooRealVar ngsigW("ngsigW","ngsignalW", fNgsigW);//*(1-.25));//3.30189e+02, 0, 1000);
  // uncomment for W background samples only
  /*
  RooRealVar meanW("meanW","meanW", 1.86583, 1.793-wind, 1.793+wind);
  RooRealVar sigmaW("sigmaW","sigmaW", 1.67071e-02, 0, 0.02);
  RooRealVar ngsigW("ngsigW","ngsignalW", 3.30189e+02, 0, 1000);
  */
  RooGaussian gaussW("gaussW","gaussW", d0_mass, meanW, sigmaW);

  // Construct exponential PDF to fit the bkg component
  RooRealVar lambda("lambda", "slope", -2.76, -10, 10);
  RooExponential expo("expo", "exponential PDF", d0_mass, lambda);
  RooRealVar blambda("blambda", "slope", -0.1, -10, 10);
  RooRealVar bsigma("bsigma","bsigma", 1, 0.1, 10);
  RooGaussian bgauss("bgauss","bgauss", d0_mass, blambda, bsigma);
  RooProdPdf prod("bkg", "expo*bgauss", RooArgList(expo, bgauss));
  
  // Construct signal + bkg PDF
  RooRealVar nsig("nsig","#signal events", 5994, 400, 70000) ;
  RooRealVar nbkge("nbkge","#background events", 84745, 8000, 800000) ;
  RooRealVar nbkg("nbkg","#background events", 84745, 8000, 800000) ;
  //RooAddPdf model("model","g+exp", RooArgList(cball, expo), RooArgList(nsig,nbkg)) ;
  //RooAddPdf model("model","g+exp", RooArgList(gauss, prod), RooArgList(nsig,nbkg)) ;
  //RooAddPdf model("model","g+exp", RooArgList(poly, expo), RooArgList(nsig,nbkg)) ;
  //RooAddPdf signalModel("signalModel","gauss1+gauss2",RooArgList(gauss,gauss1),RooArgList(ngsig,ngsig1));
  //RooAddPdf model("model","g+exp", RooArgList(signalModel, expo), RooArgList(nsig,nbkg)) ;
  RooAddPdf bkgModel("bkgModel","expo+gausskk",RooArgList(expo,gausskk),RooArgList(nbkge,ngsigkk));
  //RooAddPdf bkgModel("bkgModel","expo+gausskk",RooArgList(expo,gausskk,gaussW),RooArgList(nbkge,ngsigkk,ngsigW));
  RooAddPdf model("model","g+exp", RooArgList(gauss, bkgModel), RooArgList(nsig,nbkg)) ;
  //RooAddPdf model("model","g+exp", RooArgList(gauss, expo), RooArgList(nsig,nbkg)) ;

  cout << "fitting model" << endl;
  RooPlot* frame = d0_mass.frame() ;
  RooFitResult *r = model.fitTo(ds, SumW2Error(kTRUE), RooFit::Save(kTRUE), Name("data"));
  //RooFitResult *r = model.fitTo(ds, Extended(), SumW2Error(kTRUE), RooFit::Save(kTRUE), Name("data"));
  //ds.plotOn(frame,Binning(bins));
  //ds.plotOn(frame,Binning(60));
  ds.plotOn(frame,Binning(60));
  /*
  RooAbsReal *nll = model.createNLL(ds, NumCPU(8), SumW2Error(kTRUE));
  RooMinuit m(*nll);
  m.setPrintLevel(-1); 
  m.setPrintEvalErrors(-1);
  m.migrad();
  m.hesse();
  RooFitResult *r = m.save();
  */

  model.plotOn(frame, Name("model"));
  float chi2 = frame->chiSquare();//"model", "data", 3);
  //model.plotOn(frame, Components(expo),LineStyle(kDashed));
  model.plotOn(frame, Components(bkgModel),LineStyle(kDashed));
  frame->Draw();
  frame->SetName("massD0");
  frame->Write();
  std::cout << std::endl << "chi2 " << chi2 << std::endl << std::endl;
  //r->Print();

  if(!mass.Contains("toyData")) {
  c1->SaveAs("massD0_"+mass+"_d0.pdf");
  c1->SaveAs("massD0_"+mass+"_d0.png");
  }

  frame = ptfrac.frame();
  RooBinning bins(0,1.);
  for(int i = 0; i < bin.size(); i++)
    bins.addBoundary(bin[i]);
  ds.plotOn(frame,Binning(bins));
 
  //model.paramOn(frame);

  // Crystal ball
  /*
  cbmean.setConstant();
  cbsigma.setConstant();
  ncbsig.setConstant();
  n.setConstant();
  alpha.setConstant();
  */

  mean.setConstant();
  sigma.setConstant();
  lambda.setConstant();
  meson_l_mass.setConstant();
  ptfrac.setConstant();
  weight.setConstant();
  RooMsgService::instance().setSilentMode(true);
  //RooDataSet ds1 = RooDataSet("ds1", "ds1", &ds, *ds.get(), TString::Format("d0_mass>%f && d0_mass<%f",(1.864-wind),(1.864+wind)).Data(), "weight");
  SPlot sData("sData","An SPlot from mass", ds, &model, RooArgList(nsig,nbkg));
  std::cout << std::endl <<  "Yield of nsig is "
            << nsig.getVal() << ".  From sWeights it is "
            << sData.GetYieldFromSWeight("nsig") << std::endl;
  for(Int_t i=0; i < 10; i++) {
      std::cout << "nsig Weight   " << sData.GetSWeight(i,"nsig")
                << "   nbkg Weight   " << sData.GetSWeight(i,"nbkg")
                << "  Total Weight   " << sData.GetSumOfEventSWeight(i)
                << std::endl;
  }
  RooDataSet sigData = RooDataSet(ds.GetName(), ds.GetTitle(), &ds, *ds.get(), "", "nsig_sw");
  RooDataSet bkgData = RooDataSet(ds.GetName(), ds.GetTitle(), &ds, *ds.get(), "", "nbkg_sw");
  std::cout << sigData.weight() << std::endl;

  std::cout << sData.GetYieldFromSWeight("nsig") << std::endl;
  std::cout << sData.GetYieldFromSWeight("nbkg") << std::endl;
  std::cout << sData.GetYieldFromSWeight("nsig")/sData.GetYieldFromSWeight("nbkg") << std::endl;
  /*
  for(int i = 0; i < 10; i++) {
    std::cout << "weight " << sData.GetSWeight(i,"nsig) << std::endl;
    std::cout << "Total Weight   " << sData.GetSumOfEventSWeight(i) << std::endl;
  }
  */
  /*
  sigData.SetName("sSignal");
  sigData.SetTitle("sSignal");
  bkgData.SetName("sBackground");
  bkgData.SetTitle("sBackground");
  */
  w.import(d0_mass);
  w.import(d0_mass2);
  w.import(meson_l_mass);
  w.import(ptfrac);
  w.import(d0_pt);
  w.import(j_pt_ch);
  w.import(j_pt);
  w.import(epoch);
  w.import(weight);
  w.import(tuneW);
  w.import(model);
  //w.import(bkgData, Rename("bkgData")); //Redundant, just clone dataset and change weight to nsig_bkg
  w.import(sigData, Rename("sigData"));
  w.import(ds, Rename("dsSWeights"));
  /*
  for(auto &it : frgWgt)
    w.import(it);
  */
  /*
  for(int i=0; i < 10; i++) {
    std::cout << "Weight" << sData.GetSWeight(i,"nsig") << std::endl;
    std::cout << "Weight" << sData.GetSWeight(i,"nbkg") << std::endl;
    std::cout << "Total weight" << sData.GetSumOfEventSWeight(i) << std::endl;
  }
  */

  frame->Draw();

  int binning(22);

  /*
  frame = d0_mass.frame();
  sigData.plotOn(frame, DataError(RooAbsData::SumW2),
                 RooFit::Name("massD0_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(30));
  frame->SetName("massD0_signal");
  frame->Draw();
  frame->Write();
  c1->SaveAs("massD0_signal_"+mass+".pdf");
  c1->SaveAs("massD0_signal_"+mass+".png");

  frame = d0_mass.frame();
  bkgData.plotOn(frame, DataError(RooAbsData::SumW2),
                 RooFit::Name("massD0_background"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(30));
  frame->SetName("massD0_background");
  frame->Draw();
  frame->Write();
  c1->SaveAs("massD0_background_"+mass+".pdf");
  c1->SaveAs("massD0_background_"+mass+".png");

  RooPlot *frame2 = d0_mass.frame();
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2));
  frame2->Draw();
  */

  RooPlot *frame2 = ptfrac.frame();
  bkgData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("ptfrac_bkg"), RooFit::MarkerColor(419),
                 RooFit::MarkerStyle(24), RooFit::MarkerStyle(24),
                 RooFit::LineWidth(2), RooFit::LineColor(419), Binning(binning));
  frame2->Draw();
  frame2->SetName("ptfrac_bkg");
  TGraph *h_frac = (TGraph*)c1->GetPrimitive("ptfrac_bkg");
  if(!mass.Contains("toyData")) {
  c1->SaveAs("ptfrac_bkg_"+mass+"_d0.pdf");
  c1->SaveAs("ptfrac_bkg_"+mass+"_d0.png");
  std::cout << "Background: " << frame2->getHist()->Integral() << std::endl;
  }
  frame2 = ptfrac.frame();
  /*
  ds.plotOn(frame, RooFit::Binning(110));
  TH1F *mc = (TH1F*)convert(frame, false, 0, 1.1);
  TH1F *mc_new = (TH1F*)mc->Clone("morph");
  mc_new->Reset();
  mc->Rebin(5);
  float prev(-1);
  for(int i = 0; i < mc->GetNbinsX(); i++) {
    double frac = mc->GetBinCenter(i);
    if(frac == prev) continue;
    prev = frac;
    if(fragWeight.Length() > 0) {
      if(frac>1.) frac = 1.;
      if(frac<0.1) frac = 0.1;
      float scale = g->Eval(frac);
      int binnew = mc_new->FindBin(scale*frac);
      float n = mc->GetBinContent(i);
      float ne = mc->GetBinError(i);
      mc_new->SetBinContent(binnew, n + mc_new->GetBinContent(binnew));
      mc_new->SetBinError(binnew, sqrt(pow(ne, 2) + pow(mc->GetBinError(i), 2));
    }
    //else tuneWgt->Fill(1., 1.0); //unit weights to make re-weighting easier to loop over
  }
  mc_new->Rebin(5);
  mc_new->Write();
  */
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::Name("ptfrac_signal"), RooFit::MarkerColor(1),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(bins));
  frame2->Draw();
  std::cout << "Signal: " << frame2->getHist()->Integral() << std::endl;
  frame2->SetName("ptfrac_signal");
  frame2->SetTitle("");
  frame2->Write();

  if(jpT) mass = mass + "_jpT";
  if(!mass.Contains("toyData")) {
  c1->SaveAs("ptfrac_signal_"+mass+"_d0.pdf");
  c1->SaveAs("ptfrac_signal_"+mass+"_d0.png");
  }

  frame2 = ptfrac.frame();
  bkgData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("ptfrac_bkg"), RooFit::MarkerColor(419),
                 RooFit::MarkerStyle(24), RooFit::MarkerStyle(24),
                 RooFit::LineWidth(2), RooFit::LineColor(419), Binning(binning));
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("ptfrac_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(bins));
  frame2->SetName("ptfrac");
  frame2->SetTitle("");
  frame2->Draw();

  if(!mass.Contains("toyData")) {
  c1->SaveAs("ptfrac_"+mass+"_d0.pdf");
  c1->SaveAs("ptfrac_"+mass+"_d0.png");
  }

  frame2 = ptfrac.frame();
  bkgData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("ptfrac_bkg"), RooFit::MarkerColor(419),
                 RooFit::MarkerStyle(24), RooFit::MarkerStyle(24),
                 RooFit::LineWidth(2), RooFit::LineColor(419), Binning(bins));
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("ptfrac_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(bins));
  frame2->Draw();
  ptfrac_signal = (TH1F*)convert(frame2, false, bin);
  ptfrac_signal->SetDirectory(0);

  //show a legend
  c1->cd();
  TLegend *leg = new TLegend(0.6,0.75,0.9,0.9,"","brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry("ptfrac_signal", "Signal","p");
  leg->AddEntry("ptfrac_bkg","Background","p");
  leg->Draw();

  if(!mass.Contains("toyData")) {
  c1->SaveAs("ptfrac_"+mass+"_d0.pdf");
  c1->SaveAs("ptfrac_"+mass+"_d0.png");
  }

  binning = 25;

  frame2 = meson_l_mass.frame();
  bkgData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_bkg"), RooFit::MarkerColor(419),
                 RooFit::MarkerStyle(24), RooFit::MarkerStyle(24),
                 RooFit::LineWidth(2), RooFit::LineColor(419), Binning(binning));
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(binning));
  frame2->Draw();
  delete leg;

  //show a legend
  c1->cd();
  leg = new TLegend(0.6,0.75,0.9,0.9,"","brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry("meson_l_mass_signal", "Signal","p");
  leg->AddEntry("meson_l_mass_bkg","Background","p");
  leg->Draw();

  if(!mass.Contains("toyData")) {
  c1->SaveAs("meson_l_mass_"+mass+"_d0.pdf");
  c1->SaveAs("meson_l_mass_"+mass+"_d0.png");
  }

  /*
  frame2 = meson_l_mass.frame();
  bkgData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_bkg"), RooFit::MarkerColor(419),
                 RooFit::MarkerStyle(24), RooFit::MarkerStyle(24),
                 RooFit::LineWidth(2), RooFit::LineColor(419), Binning(binning));
  frame2->Draw();
  c1->SaveAs(""+mass+".pdf");
  c1->SaveAs(""+mass+".png");
  std::cout << "Background: " << frame2->getHist()->Integral() << std::endl;
  frame2 = meson_l_mass.frame();
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(binning));
  */
  frame2->Draw();
  std::cout << "Signal: " << frame2->getHist()->Integral() << std::endl;
  frame2->SetName("meson_l_mass_signal");
  frame2->Write();

  if(!mass.Contains("toyData")) {
  c1->SaveAs("meson_l_mass_signal_"+mass+"_d0.pdf");
  c1->SaveAs("meson_l_mass_signal_"+mass+"_d0.png");
  }

  /*
  TH1F *signalGr=(TH1F*)c1->GetPrimitive("ptfrac_signal");
  signalGr->SaveAs("/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/mtop/ptfrac_signal.pdf");
  signalGr->SaveAs("/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/mtop/ptfrac_signal.png");
  */

  TH1F *num = new TH1F("num", "num", 1, 0, 2);
  num->Fill(1,ds.numEntries());
  if(!mass.Contains("toyData")) {
  num->SetDirectory(fout);
  num->Write();
  }

  std::cout << "writting workspace to file" << std::endl;
  w.Write();
  std::cout << "writting tuneWgt to file" << std::endl;
  if(!mass.Contains("toyData")) {
  tuneWgt->Write();
  std::cout << "closing file" << std::endl;
  fout->Close();
  }
  std::cout << "DONE!" << std::endl;

  return;
}

void splot_d0(TString mass="172.5", bool isData=false, TString fragWeight="", int ep=0, bool jpT=false, int xb = 0) {
  TH1F *pdata;
  splot_d0(pdata, mass, isData, fragWeight, ep, jpT, xb);
  delete pdata;
}
