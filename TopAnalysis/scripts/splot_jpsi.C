//
// Execute the macro with
// root -l roofit_jpsi.C++
//

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLegend.h"
#include <RooFit.h>
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooStats/SPlot.h"
#include "RooPlot.h"
#include "RooGaussian.h"
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

//bool GET_BIT(short x, short b) { return (x & (1<<b) ); }

//void roofit_mtop_BCDEFGH(TString mass="166v5", TString file="jpsi_fit.root") {
//void mtop_norm(std::vector<pair<float,float>> &p, TString mass="171.5", short flags=0b00) {
void splot_jpsi(TH1F *&ptfrac_signal, TString mass="172.5", bool isData=false, TString fragWeight="", int ep=0, bool jpT=false, int xb=0) {
  RooWorkspace w("w",mass);
  float wind(0.11);
  std::vector<float> type = {0., 0.7, 1.};
  int pdf(0); //nominal, a_s=0.118, 0.117, 0.119
  mass.ReplaceAll(".","v");
  if(mass.Contains("v5") && mass != "172v5") mass = "m" + mass;
/*
void splot_jpsi_mu(RooWorkspace &w, TString mass="172.5", bool isData=false) {
*/
  //TFile *f = new TFile("../BatchJobs/merged.root"); 
  //TFile *f = new TFile("plots/plotter_mtop_BCDEFGH.root");
  TString syst("");
  //TString dir("/eos/cms/store/user/byates/top18012/Chunks/");
  TString dir("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/test/Chunks/");
  if(mass.Contains("sr") || mass.Contains("erdON") || mass.Contains("Move") || mass.Contains("ue") || mass.Contains("hdamp") || mass.Contains("m1")) { //ISR,FSR,CR,UE
    syst = mass;
    std::cout << "Processing systematics " << TString::Format("MC13TeV_TTJets_%s",syst.Data()) << std::endl;
    syst.ReplaceAll("_shift","");
  }
  else if(mass.Contains("m1")) { //top mass tunes
    syst = mass;
    //mass = "172.5";
    std::cout << "Processing systematics " << TString::Format("MC13TeV_%s",syst.Data()) << std::endl;
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
  TFile *fin;
  TGraph *g;
  TH1F *tuneWgt = new TH1F("tuneWgt","tuneWgt",2,0,2);
  //TFile *f = new TFile("MC13TeV_TTJets_m"+mass+".root");
  TFile *f = new TFile(dir+"../MC13TeV_TTJets_powheg.root"); //open a randome file to get correction histograms
  //TFile *f = new TFile(dir+"MC13TeV_TTJets_m"+mass+"_0.root"); //open a randome file to get correction histograms
  TChain *data = new TChain("data");
  if(isData) {
    data->Add(dir+"Data13TeV_*");
    //mass="Data";
  }
  //else data->Add("Chunks/MC13TeV_TTJets_m"+mass+"_*.root");
  else {
    if(mass.Contains("TRK") || mass.Contains("TRIGGER") || mass.Contains("LEP") ||mass.Contains("PU") || mass.Contains("PI") || mass.Contains("JER") || mass.Contains("JSF")) {  // SFs
      TString tmp(mass);
      tmp.ReplaceAll("_","/");
      //dir = TString("/eos/cms/store/user/byates/top18012/" + tmp + "/Chunks/");
      dir = TString("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/test/" + tmp + "/");//Chunks/");
      std::cout << dir << std::endl;
    }
    //std::vector<TString> mcSamples = { "MC13TeV_TTJets_powheg" };
    std::vector<TString> mcSamples;
    if(syst.Length()>0) mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW",TString::Format("MC13TeV_TTJets_%s",syst.Data()),"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    //GluonMove and QCD only!
    //if(syst.Length()>0) mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW",TString::Format("MC13TeV_%s",syst.Data()),"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    //else mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_TTJets_m"+mass,"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    else if(mass.Contains("psweights") || mass.Contains("FSR")) mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_2016_TTJets_psweights","MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    else mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_TTJets_powheg","MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    if(mass.Contains("QCD") || mass.Contains("GluonMove")) mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW",TString::Format("MC13TeV_%s*",syst.Data()),"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    for(auto & it : mcSamples) {
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
  TString fUrl("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_"+mass+"_sPlot_jpsi.root");
  if(ep>0) fUrl.ReplaceAll(".root",TString::Format("%d.root",ep));
  if(jpT) fUrl.ReplaceAll(".root", "_jpT.root");
  if(xb != 0) fUrl.ReplaceAll(".root",TString::Format("_%d.root",int(type[xb]*100)));
  std::cout << "creating file: "  << fUrl<< std::endl;
  TFile *fout;// = new TFile(fUrl,"RECREATE");
  if(!mass.Contains("toyData"))
  fout = new TFile(fUrl,"RECREATE");
  //data->Add("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/Chunks/MC13TeV_W*Jets_*.root");
  //disabled top pT reweighting
  //data->Draw("jpsi_mass>>h(30,1.7,2.0)","norm*sfs*puwgt*(meson_id==443 && jpsi_l3d/jpsi_sigmal3d>10 && HT>100 && jpsi_sigmal3d>2E-4)");// && j_hadflav[jpsi_j]==5 && jpsi_pi_mother==421 && jpsi_k_mother==421)","goff");
  //data->Draw("jpsi_mass>>h(30,1.7,2.0)","norm*sfs*puwgt*topptwgt*(meson_id==42113 && jpsi_l3d/jpsi_sigmal3d>10 && HT>180 && jpsi_sigmal3d>2E-4)");// && j_hadflav[jpsi_j]==5 && jpsi_pi_mother==421 && jpsi_k_mother==421)","goff");
  
  CharmEvent_t ev;
  attachToCharmEventTree(data,ev);
  TH1F *h = (TH1F*)gDirectory->Get("h");
  //if(!isData) h->Scale(35864);
  //if(!isData) h->Scale(832*35864);
  //f->ls(); 
  TString name = "massJPsi_l_all_d0";
  //TString name = "massJPsi_l_all_d0_BCDEF";
  TH1F *pu1 = (TH1F*) f->Get("puwgtctr_BCDEF");
  TH1F *pu2 = (TH1F*) f->Get("puwgtctr_GH");
  TH1F *top1 = (TH1F*) f->Get("topptwgt_BCDEF");
  TH1F *top2 = (TH1F*) f->Get("topptwgt_GH");
  float puSF1 = pu1->GetBinContent(1)/pu1->GetBinContent(2);
  float puSF2 = pu2->GetBinContent(1)/pu2->GetBinContent(2);
  float topSF1 = top1->GetBinContent(2)/top1->GetBinContent(1);
  float topSF2 = top2->GetBinContent(2)/top2->GetBinContent(1);
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
  cout << "loaded!" << endl;

std::vector<TH1F*> hCorrBtop;
std::vector<TH1F*> hCorrBpi;
std::vector<TH1F*> hCorrBpu;
std::vector<TH1F*> hCorrGtop;
std::vector<TH1F*> hCorrGpi;
std::vector<TH1F*> hCorrGpu;
/*
{
std::vector<TString> mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_TTJets_powheg","MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
for(auto it : mcSamples) {
  TFile *file = new TFile("LJets2015/2016/test/" + it + ".root");
  hCorrBtop.push_back((TH1F*)file->Get("topptwgt_BCDEF")->Clone());
  hCorrBpi.push_back((TH1F*)file->Get("piwgt_BCDEF")->Clone());
  hCorrBpu.push_back((TH1F*)file->Get("puwgtctr_BCDEF")->Clone());

  hCorrGtop.push_back((TH1F*)file->Get("topptwgt_GH")->Clone());
  //hCorrGpi.push_back((TH1F*)file->Get("piwgt_GH")->Clone());
  hCorrGpu.push_back((TH1F*)file->Get("puwgtctr_GH")->Clone());
}
}
*/

  // Declare observable x
  RooRealVar jpsi_mass("jpsi_mass","J/#Psi mass", 2.8, 3.4, "GeV") ;
  RooRealVar weight("weight","weight",1.,0.,36000.);
  RooRealVar meson_l_mass("jpsi_l_mass","J/#Psi+l mass", 0, 250, "GeV") ;
  RooRealVar ptfrac("ptfrac","J/#Psi #it{p}_{T} / #Sigma #it{p}_{T}^{ch}", 0, 1.1, "") ;
  RooRealVar jpsi_pt("jpsi_pt","J/#Psi #it{p}_{T}", 0, 250, "GeV");
  RooRealVar j_pt_ch("j_pt_ch","j #it{p}_{T} charged", 0, 400, "GeV");
  RooRealVar j_pt("j_pt","j #it{p}_{T}", 0, 400, "GeV");
  RooRealVar j_hadflav("j_hadflav","j hadflav", 0, 6, "");
  /*
  RooCategory mcFile("mcFile","mcFile");
  {
    std::vector<TString> mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_TTJets_powheg","MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
  for(auto & it : mcSamples) {
    mcFile.defineType(it.Data());
  }
  }
  */
  //RooRealVar mcFile("mcFile","mcFile", 0, 6, "");
  RooRealVar mcFile("mcFile","mcFile", 0, 20);
  RooRealVar epoch("epoch","epoch",1,2);
  RooRealVar tuneW = RooRealVar("tuneW", "tuneW", 1., 0, 2.);
  
  //cout << "creating dataset" << endl;
  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'

  //RooDataSet dsn("dsn", "dsn", RooArgSet(meson_id,jpsi_mass,ptfrac,meson_l_mass,weight), Import(*data), Cut("meson_id==42113"));
  //RooDataSet dsn("dsn", "dsn", RooArgSet(jpsi_mass,ptfrac,meson_l_mass,weight,jpsi_pt,j_pt_ch,j_pt,tuneW));
  RooDataSet dsn("dsn", "dsn", RooArgSet(RooArgList(jpsi_mass,ptfrac,meson_l_mass,weight,jpsi_pt,epoch,j_pt_ch,j_pt,tuneW,""),RooArgList(j_hadflav,mcFile,"")));
  //RooDataSet dsn("dsn", "dsn", RooArgSet(jpsi_mass,ptfrac,meson_l_mass,weight));
  int nmc = ep == 1 ? 7327 : 9539;
  int nentries = ep == 1 ? 1339 : 1329;
  int nset(0);
  float sumFSRwgt(0);
  std::cout << "Total events: " << data->GetEntries() << std::endl;
  for(int i=0; i < data->GetEntries(); i++) {
    ev = {};
    data->GetEntry(i);
    for(int j=0; j<ev.nmeson; j++) {
      //int j(0);
      if(mass.Contains("_c") && 0) {
        if(ev.j_hadflav[j] == 5 && TString(data->GetFile()->GetName()).Contains("powheg")) continue;
        else if(mass.Contains("ttbar_c") && !TString(data->GetFile()->GetName()).Contains("powheg")) continue;
        else if(mass.Contains("W_c") && TString(data->GetFile()->GetName()).Contains("powheg")) continue;
      }
      if(ev.meson_id[j] != 443) continue;
      if(ev.jpsi_mass[j] < 2.8) continue;
      if(ev.jpsi_mass[j] > 3.4) continue;
      if(ep>0 && ev.epoch[j]!=ep) continue;
      /*
      if(ev.jpsi_l3d[j]/ev.jpsi_sigmal3d[j]<10) continue;
      if(ev.jpsi_sigmal3d[j]<2E-4) continue;
      */
      if(ev.ht<80) continue;
      if(mass.Contains("toyData")) {
        /*
        if(i < start || i > start + nentries) continue;
        */
        TRandom3 *rand = new TRandom3(0);
        /*
        if(ep==1 && rand->Uniform(0, 1) > 0.2) continue;
        if(ep==2 && rand->Uniform(0, 1) > 0.14) continue;
        */
        if(rand->Uniform(0, 1) > float(nentries)/float(nmc)) continue;
        if(nset > nentries) continue;
        nset++;
      }
      /*
      if(!ev.isData) {
        if(ev.j_hadflav[ev.jpsi_j[j]]!=5) continue;
        if(ev.jpsi_pi_mother[j]!=421) continue;
        if(ev.jpsi_k_mother[j]!=421) continue;
      }
      */
      //if(j_pt[(int)jpsi_j[j]]<30) continue;
      /*
      if(!(mesonlm[j] > 0)) continue;
      if(mesonlm[j] > 250) continue;
      */
      if(xb != 0) {
	if(type[xb] < 1) { if(ev.jpsi_pt[j]/ev.j_pt_charged[j] > type[xb] || ev.jpsi_pt[j]/ev.j_pt_charged[j] < type[xb-1]) continue; }
        else if(ev.jpsi_pt[j]/ev.j_pt_charged[j] < type[xb-1]) continue; 
      }
      float scale = 1.;
      tuneW.setVal(1.);
      scale *= ev.sfs[j];
      if(!isData) {
        if(mass.Contains("toyData"))
        scale *= ev.puwgt[j] * ev.topptwgt;// * topSF * puSF;
        else
        scale *= ev.norm * ev.xsec * ev.puwgt[j] * ev.topptwgt;// * topSF * puSF;
        if(pdf>0 && ev.ttbar_nw > pdf) scale *= ev.ttbar_w[pdf];//alternate PDF event weight
        //scale = ev.norm * ev.xsec * ev.sfs[j] * ev.puwgt[j] * ev.topptwgt;// * topSF * puSF;
        //scale = norm * sfs[j] * puwgt[j] * topptwgt * topSF * puSF;
        //PS weights
        if(mass.Contains("FSR-up")) { scale *= ev.gfsr[0]; sumFSRwgt += ev.gfsr[0]; }
        else if(mass.Contains("FSR-down")) { scale *= ev.gfsr[1]; sumFSRwgt += ev.gfsr[1]; }
        if(ev.epoch[j]==1) {
          //scale *= 1.11;
          if(mass.Contains("toyData"))
          scale =  scale * puSF1 * topSF1;
          else
          scale =  scale * 19712.86 * puSF1 * topSF1;
          //h1->Fill(mesonlm[j], scale);
        }
        else if(ev.epoch[j]==2) {
          //scale *= 1.01;
          if(mass.Contains("toyData"))
          scale = scale * puSF2 * topSF2;
          else
          scale = scale * 16146.178 * puSF2 * topSF2;
          //h2->Fill(mesonlm[j], scale);
        }
        else
          continue;
        tuneWgt->Fill(0.,1.0);
        if(fragWeight.Length() > 0) {
          double frac(ev.jpsi_pt[j]/ev.j_pt_charged[j]);
          if(frac>1.) frac = 1.;
          if(frac<0.1) frac = 0.1;
          scale *= g->Eval(frac);
          tuneWgt->Fill(1.,g->Eval(frac));
        }
        else tuneWgt->Fill(1., 1.0); //unit weights to make re-weighting easier to loop over
      }
      if(ev.epoch[j]!=1 && ev.epoch[j]!=2) continue;
      epoch.setVal(ev.epoch[j]);
      //jpsi_mass = ev.jpsi_mass[j];
      jpsi_mass.setVal(ev.jpsi_mass[j]);
      //if(ev.jpsi_mass[j]>(2.95) && ev.jpsi_mass[j]<(3.097+wind) && ev.jpsi_pt[j]>0) { //Symmetric around PDG mass
      //if(ev.jpsi_mass[j]>(3.097-wind) && ev.jpsi_mass[j]<(3.097+wind) && ev.jpsi_pt[j]>0) { //Symmetric around PDG mass
      //if(ev.jpsi_mass[j]>3.0 && ev.jpsi_mass[j]<3.2 && ev.jpsi_pt[j]>0) {
      //if(ev.jpsi_mass[j]>1.81 && ev.jpsi_mass[j]<1.91 && ev.jpsi_pt[j]>0) {
      if(!jpT) ptfrac.setVal(ev.jpsi_pt[j]/ev.j_pt_charged[j]);
      else ptfrac.setVal(ev.jpsi_pt[j]/ev.j_pt[j]);
      jpsi_pt.setVal(ev.jpsi_pt[j]);
      j_pt_ch.setVal(ev.j_pt_charged[j]);
      j_pt.setVal(ev.j_pt[j]);
      if(!isData) j_hadflav.setVal(ev.j_hadflav[j]);
      if(!isData) {
/*
      if(TString(data->GetFile()->GetName()).Contains("powheg")) mcFile.setVal(5);
      else if(TString(data->GetFile()->GetName()).Contains("W1Jets")) mcFile.setVal(4);
      else if(TString(data->GetFile()->GetName()).Contains("W2Jets")) mcFile.setVal(4);
      else if(TString(data->GetFile()->GetName()).Contains("W3Jets")) mcFile.setVal(4);
      else if(TString(data->GetFile()->GetName()).Contains("W4Jets")) mcFile.setVal(4);
      else mcFile.setVal(0);
*/
       TString tmp = TString(data->GetFile()->GetName());
       tmp.ReplaceAll("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/test/Chunks/", "");
       int ind(-1);
       if(tmp.Contains("ext"))
       ind = tmp.Index(TRegexp("ext_[0-9]*.root"));
       else
       ind = tmp.Index(TRegexp("_[0-9]*.root"));
       tmp.Remove(ind, tmp.Length()-ind);
       //std::cout << tmp << std::endl;
    std::vector<TString> mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_TTJets_powheg","MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
       ind = distance(mcSamples.begin(), find(mcSamples.begin(), mcSamples.end(), tmp));
       if(ev.epoch[j]==1 && 0) {
         scale /= (puSF1 * topSF1);
         scale *= hCorrBtop[ind]->GetBinContent(2)/hCorrBtop[ind]->GetBinContent(1);
         scale *= hCorrBpi[ind]->GetBinContent(1)/hCorrBpi[ind]->GetBinContent(2);
         scale *= hCorrBpu[ind]->GetBinContent(1)/hCorrBpu[ind]->GetBinContent(2);
       }
       else if(ev.epoch[j]==2 && 0) {
         scale /= (puSF1 * topSF1);
         scale *= hCorrGtop[ind]->GetBinContent(2)/hCorrGtop[ind]->GetBinContent(1);
         scale *= hCorrGpu[ind]->GetBinContent(1)/hCorrGpu[ind]->GetBinContent(2);
       }
       //std::cout << ind << std::endl;
       //mcFile.setLabel(tmp);
       mcFile.setVal(ind);
       //mcFile.setVal(data->GetFile()->GetName());
      }
      //meson_l_mass = ev.jpsi_l_mass[j];
      meson_l_mass.setVal(ev.jpsi_l_mass[j]);
      //}
      weight.setVal(scale);
      tuneW.setVal(scale * tuneW.getVal());
      //dsn.add(RooArgSet(jpsi_mass,meson_l_mass,ptfrac,weight,jpsi_pt,j_pt_ch,j_pt,epoch,tuneW));
      dsn.add(RooArgSet(RooArgList(jpsi_mass,meson_l_mass,ptfrac,weight,jpsi_pt,j_pt_ch,j_pt,epoch,tuneW,""),RooArgList(j_hadflav,mcFile,"")));
      //dsn.add(RooArgSet(j_pt));
      //dsn.add(RooArgSet(jpsi_mass,meson_l_mass,ptfrac,weight));
      //dsn.add(RooArgSet(jpsi_mass,meson_l_mass,ptfrac,weight), scale);
      //within D^0 mass peak (1.864 +/- 0.05)
      //if(ev.jpsi_l_mass[j] > 1.8 && ev.jpsi_l_mass[j] < 1.93)
        //dsn.add(RooArgSet(meson_l_mass,ptfrac), scale);

    }
  }

  /*
  double tuneShape = 1.;
  if(tuneWgt->GetBinContent(2) > 0) tuneShape = tuneWgt->GetBinContent(1)/tuneWgt->GetBinContent(2);
  std::cout << "Norm weight: " << tuneShape << std::endl;
  //very inefficient second loop to adjust normalization
  std::cout << "Re-weighting based on norm weight" << std::endl;
  for(int i = 0; i < dsn.numEntries(); i++) {
    double tmpW = dsn.get(i)->getRealValue("tuneW");
    //tuneW.setVal(tmpW * tuneShape);
    weight.setVal(tuneW.getVal() / sumFSRwgt);
    //std::cout << tuneW.getVal() << " " << weight.getVal() << std::endl;;
  }
  std::cout << "Re-weighting done!" << std::endl;
  std::cout << dsn.numEntries() / sumFSRwgt << std::endl;
  */
  //RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, "weight");
  RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, "tuneW");
  dsn.Print("v");
  std::cout << dsn.weight() << std::endl;
  ds.Print("v");
  ds.Print();
  std::cout << ds.weight() << std::endl;
  cout << "plotting dataset" << endl;
  // Make plot of binned dataset showing Poisson error bars (RooFit default)

  cout << "defining variables" << endl;

  // Mean of the J/psi mass peak
  RooRealVar mean("mean","mean", 3.097, 3.097-wind, 3.097+wind);
  RooRealVar mean3("mean3","mean3", 3.097, 3.097-wind, 3.097+wind);
  //mean.setVal(3.0969);
  //mean.setConstant();

  // Construct Crystal Ball PDF for signal
  RooRealVar cbsigma("csigma", "cbsigma" , 0.02, 0.001, 0.25); 
  RooRealVar ncbsig("ncbsig", "ncbsignal", 4000, 0, 10000); 
  RooRealVar n("n","cbn", 5, 0, 10);
  RooRealVar alpha("alpha","cbalpha", 1, 0, 5);
  RooCBShape cball("cball", "crystal ball", jpsi_mass, mean, cbsigma, alpha, n);

  // Construct Gaussian1 PDF for signal
  RooRealVar sigma1("sigma1","sigma1", 0.02, 0.001, 0.3);
  RooRealVar ngsig1("ngsig1","ngsignal1", 100, 0, 10000);
  RooGaussian gauss1("gauss1","gauss1", jpsi_mass, mean, sigma1);

  // Construct Gaussian2 PDF for signal
  RooRealVar sigma2("sigma2","sigma2", 0.02, 0.01, 0.5);
  RooRealVar ngsig2("ngsig2","ngsignal2", 100, 0, 10000);
  RooGaussian gauss2("gauss2","gaus2s", jpsi_mass, mean, sigma2);

  RooRealVar sigma3("sigma3","sigma3", 2, 1, 5);
  RooRealVar ngsig3("ngsig3","ngsignal3", 100, 0, 10000);
  RooGaussian gauss3("gauss3","gaus3s", jpsi_mass, mean3, sigma3);

  // Construct a double Gaussian function to fit the signal component
  RooAddPdf signalModel("signal model","gauss1+gauss2",RooArgList(gauss1,gauss2),RooArgList(ngsig1,ngsig2));
  //RooAddPdf signalModel("signal model","gauss1+gauss2",RooArgList(gauss1,gauss2,gauss3),RooArgList(ngsig1,ngsig2,ngsig3));

  // Construct exponential PDF to fit the bkg component
  RooRealVar lambda("lambda", "slope", -2, -5, 5.);
  RooExponential expo("expo", "exponential PDF", jpsi_mass, lambda);

  RooRealVar cbmean("cbmean", "cbmean" , 3.093);
  RooRealVar bkgcbsigma("csigma", "cbsigma" , 0.034);
  RooRealVar bkgncbsig("ncbsig", "ncbsignal", 194);
  RooRealVar bkgn("n","cbn", 10);
  RooRealVar bkgalpha("alpha","cbalpha", 1.1);
  RooCBShape bkgcball("bkgcball", "crystal ball", jpsi_mass, cbmean, bkgcbsigma, bkgalpha, bkgn);

  //std::vector< std::vector<float> > param_c = {{1.86644,1.36245e+03,2.04736e-02,-2.66732,2.40061e+04}, {1.86562,1.27583e+03,2.01833e-02,-2.70638,2.71743e+04}};
  std::vector< std::vector<float> > param_c = {{3.09803e+00,1.17511e+00,9.99884e+00,1.93583e+02,3.49406e-02,-1.36613e+00,4.28522e+01}, {3.09766e+00,1.05342e+00,9.99590e+00,1.93656e+02,3.27095e-02,-1.01739e+00,4.88592e+01}};

  float fMeanc = param_c[ep-1][0];
  float fNc = param_c[ep-1][1];
  float fNalphac = param_c[ep-1][2];
  float fNcbsigc = param_c[ep-1][3];
  float fNsigmac = param_c[ep-1][4];
  float fLambdac = param_c[ep-1][5];
  float fNbkgec = param_c[ep-1][6];
  //if(name.Contains("Wup")) fNcbsigc *= 1.25; 
  //if(name.Contains("Wdown")) fNcbsigc *= 0.75; 
  RooRealVar meanc("meanc", "meanc" , fMeanc);
  RooRealVar ncbsigc("ncbsigc","ncbsignalc", fNcbsigc*(1+.25));//3.30189e+02, 0, 1000);
  RooRealVar cbsigmac("csigmac", "cbsigmac", fNsigmac);
  RooRealVar nc("nc","cbnc", fNc);
  RooRealVar alphac("alphac","cbalphac", fNalphac);
  RooCBShape cballc("cballc", "crystal ballc", jpsi_mass, meanc, cbsigmac, alphac, nc);
  RooRealVar lambdac("lambdac", "slope", fLambdac);
  RooRealVar nbkgec("nbkgec","nbkgenalc", fNbkgec*(1-.25));//3.30189e+02, 0, 1000);
  RooExponential expoc("expoc", "exponential PDF", jpsi_mass, lambdac);
  
  // Construct signal + bkg PDF
  RooRealVar nsig("nsig","#signal events", 4000, 0, 10000) ;
  RooRealVar nbkg("nbkg","#background events", 4000, 0, 10000) ;
  //RooAddPdf bkgModel("bkgModel","expo+gausskk",RooArgList(expo,expoc,gaussc),RooArgList(nbkge,nbkgec,ngsigc));
  //RooAddPdf model("model","g+a", RooArgList(cball, expo), RooArgList(nsig,nbkg)) ;
  RooAddPdf model("model","g+a", RooArgList(cball, expo, cballc), RooArgList(nsig,nbkg,ncbsigc)) ;
  //RooAddPdf model("model","g+a", RooArgList(cball, expo, bkgcball), RooArgList(nsig,nbkg,bkgncbsig)) ;
  //RooAddPdf model("model","g+a", RooArgList(signalModel, expo, cballc), RooArgList(nsig,nbkg,ncbsigc)) ;

  cout << "fitting model" << endl;
  meson_l_mass.setConstant();
  //ptfrac.setConstant();
  weight.setConstant();
  RooPlot* frame = jpsi_mass.frame() ;
  std::vector<float> bin;
  RooBinning bins(0,1.1);
  //bin = {0, 0.2, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
  bin = {0, 0.2, 0.4, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
  for(int i = 0; i < bin.size(); i++)
    bins.addBoundary(bin[i]);
  model.fitTo(ds, Extended(), SumW2Error(kTRUE));
  ds.plotOn(frame,Binning(30));
  /*
  RooAbsReal *nll = model.createNLL(ds, NumCPU(8), SumW2Error(kTRUE));
  RooMinuit m(*nll);
  m.setPrintLevel(-1); 
  m.setPrintEvalErrors(-1);
  m.migrad();
  m.hesse();
  RooFitResult *r = m.save();
  */

  model.plotOn(frame);
  float chi2 = frame->chiSquare();//"model", "data", 3);
  //model.plotOn(frame, Components(RooArgSet(expo)),LineStyle(kDashed));
  model.plotOn(frame, Components(RooArgSet(expo,cballc)),LineStyle(kDashed));
  frame->Draw();
  frame->SetName("massJPsi");
  frame->Write();
  std::cout << std::endl << "chi2 " << chi2 << std::endl << std::endl;

  if(!mass.Contains("toyData")) {
  c1->SaveAs("massJPsi_"+mass+".pdf");
  c1->SaveAs("massJPsi_"+mass+".png");
  }

  frame = ptfrac.frame();
  ds.plotOn(frame,Binning(24));
 
  //model.paramOn(frame);

  SPlot sData("sData","An SPlot from mass", ds, &model, RooArgList(nsig,nbkg));
  RooDataSet sigData = RooDataSet(ds.GetName(), ds.GetTitle(), &ds, *ds.get(), "", "nsig_sw");
  //sigData = *(RooDataSet*)sigData.reduce(Cut("jpsi_mass>1.8 && jpsi_mass<1.93"));
  std::cout << sigData.weight() << std::endl;
  RooDataSet bkgData = RooDataSet(ds.GetName(), ds.GetTitle(), &ds, *ds.get(), "", "nbkg_sw");
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
  w.import(jpsi_mass, Silence());
  w.import(meson_l_mass, Silence());
  w.import(ptfrac, Silence());
  w.import(jpsi_pt, Silence());
  w.import(j_pt_ch, Silence());
  w.import(j_pt, Silence());
  w.import(j_hadflav);
  w.import(mcFile);
  w.import(epoch, Silence());
  w.import(weight, Silence());
  w.import(model);
  w.import(sigData, Rename("sigData"), Silence());
  //w.import(bkgData, Rename("bkgData"), Silence());
  w.import(ds, Rename("dsSWeights"), Silence());
  /*
  for(int i=0; i < 10; i++) {
    std::cout << "Weight" << sData.GetSWeight(i,"nsig") << std::endl;
    std::cout << "Weight" << sData.GetSWeight(i,"nbkg") << std::endl;
    std::cout << "Total weight" << sData.GetSumOfEventSWeight(i) << std::endl;
  }
  */

  frame->Draw();

  int binning(22);

  frame = jpsi_mass.frame();
  sigData.plotOn(frame, DataError(RooAbsData::SumW2),
                 RooFit::Name("massJPsi_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(30));
  frame->SetName("massJPsi_signal");
  frame->Draw();
  frame->Write();
  /*
  c1->SaveAs("massJPsi_signal_"+mass+".pdf");
  c1->SaveAs("massJPsi_signal_"+mass+".png");
  */

  frame = jpsi_mass.frame();
  bkgData.plotOn(frame, DataError(RooAbsData::SumW2),
                 RooFit::Name("massJPsi_background"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(30));
  frame->SetName("massJPsi_background");
  frame->Draw();
  frame->Write();
  /*
  c1->SaveAs("massJPsi_background_"+mass+".pdf");
  c1->SaveAs("massJPsi_background_"+mass+".png");
  */

  /*
  RooPlot *frame2 = jpsi_mass.frame();
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
  frame2->Write();
  TGraph *h_frac = (TGraph*)c1->GetPrimitive("ptfrac_bkg");
  mass = mass + "_jpsi";
  /*
  c1->SaveAs("ptfrac_bkg_"+mass+".pdf");
  c1->SaveAs("ptfrac_bkg_"+mass+".png");
  */
  std::cout << "Background: " << frame2->getHist()->Integral() << std::endl;
  frame2 = ptfrac.frame();
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("ptfrac_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(bins));
  frame2->Draw();
  std::cout << "Signal: " << frame2->getHist()->Integral() << std::endl;
  frame2->SetName("ptfrac_signal");
  frame2->Write();

  if(jpT) mass = mass + "_jpT";
  if(!mass.Contains("toyData")) {
  c1->SaveAs("ptfrac_signal_"+mass+".pdf");
  c1->SaveAs("ptfrac_signal_"+mass+".png");
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
  frame2->Draw();
  //ptfrac_graph = (TGraph*)frame2->getHist()->Clone();
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

  /*
  c1->SaveAs("ptfrac_"+mass+".pdf");
  c1->SaveAs("ptfrac_"+mass+".png");
  */

  binning = 25;

  frame2 = meson_l_mass.frame();
  bkgData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_jpsi_bkg"), RooFit::MarkerColor(419),
                 RooFit::MarkerStyle(24), RooFit::MarkerStyle(24),
                 RooFit::LineWidth(2), RooFit::LineColor(419), Binning(binning));
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_jpsi_signal"), RooFit::MarkerColor(1),
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
  leg->AddEntry("meson_l_mass_jpsi_signal", "Signal","p");
  leg->AddEntry("meson_l_mass_jpsi_bkg","Background","p");
  leg->Draw();

  /*
  c1->SaveAs("meson_l_mass_jpsi_"+mass+".pdf");
  c1->SaveAs("meson_l_mass_jpsi_"+mass+".png");
  */

  frame2 = meson_l_mass.frame();
  bkgData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_jpsi_bkg"), RooFit::MarkerColor(419),
                 RooFit::MarkerStyle(24), RooFit::MarkerStyle(24),
                 RooFit::LineWidth(2), RooFit::LineColor(419), Binning(binning));
  frame2->Draw();
  /*
  c1->SaveAs("meson_l_mass_jpsi_bkg_"+mass+".pdf");
  c1->SaveAs("meson_l_mass_jpsi_bkg_"+mass+".png");
  */
  std::cout << "Background: " << frame2->getHist()->Integral() << std::endl;
  frame2 = meson_l_mass.frame();
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_jpsi_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(binning));
  frame2->Draw();
  std::cout << "Signal: " << frame2->getHist()->Integral() << std::endl;
  frame2->SetName("meson_l_mass_jpsi_signal");
  frame2->Write();

  /*
  c1->SaveAs("meson_l_mass_jpsi_signal_"+mass+".pdf");
  c1->SaveAs("meson_l_mass_jpsi_signal_"+mass+".png");
  */

  /*
  TH1F *signalGr=(TH1F*)c1->GetPrimitive("ptfrac_jpsi_signal");
  signalGr->SaveAs("ptfrac_signal.pdf");
  signalGr->SaveAs("ptfrac_signal.png");
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

void splot_jpsi(TString mass="172.5", bool isData=false, TString fragWeight="", int ep=0, bool jpT=false, int xb=0) {
  TH1F *pdata;
  splot_jpsi(pdata, mass, isData, fragWeight, ep, jpT, xb);
  delete pdata;
}

/*
void splot_jpsi(TString mass="172.5", bool isData=false, TString fragWeight="", int ep=0, bool jpT=false) {
  TH1F *pdata;
  TGraph *pgraph;
  splot_jpsi(pdata, pgraph, mass, isData, fragWeight, ep, jpT);
  delete pdata;
  delete pgraph;
}

void splot_jpsi(TGraph *&pgraph, TString mass="172.5", bool isData=false, TString fragWeight="", int ep=0, bool jpT=false) {
  TH1F *pdata;
  splot_jpsi(pdata, pgraph, mass, isData, fragWeight, ep, jpT);
  delete pdata;
}
*/
