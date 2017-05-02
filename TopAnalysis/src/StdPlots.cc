#include <iostream>
#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"
#include "TopLJets2015/TopAnalysis/interface/StdPlots.h"
#include "TopLJets2015/TopAnalysis/interface/Jet.h"

StdPlots::StdPlots(TString runPeriod, bool debug) {  
  runPeriod_ = "_"+runPeriod;
  debug_ = debug;
  norm_ = 1.;
  sfs_ = 1.;
  puWgt_ = 1.;

  std::vector<TString> lfsVec = { "_all", "_e", "_ee", "_em", "_mm", "_m" }; 
  std::vector<TString> cutVec = { "", "_lep", "_lepjets", "_jpsi", "_csv", "_meson" };
  std::vector<TString> wgtVec = { "", "_no_weight" };
  for(int i = 0; i < (int)lfsVec.size(); i++) {
  for(int j = 0; j < (int)cutVec.size(); j++) {
  for(int k = 0; k < (int)wgtVec.size(); k++) {
    TString tag(lfsVec[i]);
    TString cut(cutVec[j]);
    TString weight(wgtVec[k]);
    allPlots["lp_pt_iso"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt_iso"+tag+cut+weight+runPeriod_,";Lepton P_{T} [GeV] after cleaning;Events / 10 GeV", 20, 0,200);
    allPlots["lp_pt_veto"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt_veto"+tag+cut+weight+runPeriod_,";Lepton P_{T} [GeV] after veto;Events / 10 GeV", 20, 0,200);
    allPlots["lp_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt_low"+tag+cut+weight+runPeriod_,";Leading Lepton P_{T} [GeV];Events / 1 GeV", 20, 20,40);
    allPlots["lp_pt"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt"+tag+cut+weight+runPeriod_,";Leading lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["l2p_pt"+tag+cut+weight+runPeriod_] = new TH1F("l2p_pt"+tag+cut+weight+runPeriod_,";Sub-leading lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_pt"+tag+cut+weight+runPeriod_] = new TH1F("dilp_pt"+tag+cut+weight+runPeriod_,";Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_m"+tag+cut+weight+runPeriod_] = new TH1F("dilp_m"+tag+cut+weight+runPeriod_,";M_{ll} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["j_pt"+tag+cut+weight+runPeriod_] = new TH1F("j_pt"+tag+cut+weight+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["lj_pt"+tag+cut+weight+runPeriod_] = new TH1F("lj_pt"+tag+cut+weight+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["bj_pt"+tag+cut+weight+runPeriod_] = new TH1F("bj_pt"+tag+cut+weight+runPeriod_,";Leading b Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["j_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("j_pt_low"+tag+cut+weight+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["lj_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("lj_pt_low"+tag+cut+weight+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["bj_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("bj_pt_low"+tag+cut+weight+runPeriod_,";Leading b Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["nlp"+tag+cut+weight+runPeriod_]     = new TH1F("nlp"+tag+cut+weight+runPeriod_,";N_{l};Events" ,3,0.,3.);
    allPlots["ndilp"+tag+cut+weight+runPeriod_]     = new TH1F("ndilp"+tag+cut+weight+runPeriod_,";N_{ll};Events" ,3,0.,3.);
    allPlots["nj"+tag+cut+weight+runPeriod_]     = new TH1F("nj"+tag+cut+weight+runPeriod_,";N_{jets} (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["nlj"+tag+cut+weight+runPeriod_]     = new TH1F("nlj"+tag+cut+weight+runPeriod_,";N_{jets} (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["nbj"+tag+cut+weight+runPeriod_]     = new TH1F("nbj"+tag+cut+weight+runPeriod_,";N_{b-jets} (CSV > 0.8);Events" ,4,1.,5.);
    allPlots["npf"+tag+cut+weight+runPeriod_]     = new TH1F("npf"+tag+cut+weight+runPeriod_,";N_{pf};Events / 10" ,5,0.,5.);
    allPlots["lp_eta"+tag+cut+weight+runPeriod_]  = new TH1F("lp_eta"+tag+cut+weight+runPeriod_,";Leading lepton #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["l2p_eta"+tag+cut+weight+runPeriod_]  = new TH1F("l2p_eta"+tag+cut+weight+runPeriod_,";Sub-Leading lepton #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["lp_phi"+tag+cut+weight+runPeriod_]  = new TH1F("lp_phi"+tag+cut+weight+runPeriod_,";Leading lepton #phi; Events", 50, -3.14,3.14);
    allPlots["l2p_phi"+tag+cut+weight+runPeriod_]  = new TH1F("l2p_phi"+tag+cut+weight+runPeriod_,";Sub-Leading lepton #phi; Events", 50, -3.14,3.14);
    allPlots["nstart"+tag+cut+weight+runPeriod_]     = new TH1F("jetindex"+tag+cut+weight+runPeriod_,";N_{jetindex};Events" ,5,0.,5.);
    allPlots["pfid"+tag+cut+weight+runPeriod_]     = new TH1F("pfid"+tag+cut+weight+runPeriod_,";PFID;Events" ,440,-220.,220.);
    allPlots["massJPsi"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi"+tag+cut+weight+runPeriod_,";M_{ll};Events / 18 MeV" ,50,2.5,3.4);
    allPlots["massJPsiK"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsiK"+tag+cut+weight+runPeriod_,";M_{llk};Events / 15 MeV" ,100,4.5,6);
    allPlots["massD0"+tag+cut+weight+runPeriod_]     = new TH1F("massD0"+tag+cut+weight+runPeriod_,";M_{D^{0}};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_lep"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_lep"+tag+cut+weight+runPeriod_,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_mu"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_mu"+tag+cut+weight+runPeriod_,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_e"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_ele"+tag+cut+weight+runPeriod_,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massDsmD0loose"+tag+cut+weight+runPeriod_]     = new TH1F("massDsmD0loose"+tag+cut+weight+runPeriod_,";M_{K#pi#pi} - M_{K#pi};Events / 0.5 MeV" ,20,0.14,0.16);
    allPlots["massDsmD0"+tag+cut+weight+runPeriod_]     = new TH1F("massDsmD0"+tag+cut+weight+runPeriod_,";M_{K#pi#pi} - M_{K#pi};Events / 0.5 MeV" ,20,0.14,0.16);
    allPlots["massDs"+tag+cut+weight+runPeriod_]     = new TH1F("massDs"+tag+cut+weight+runPeriod_,";M_{D^{*}};Events / 10 MeV" ,200,0.,2.0);
    allPlots["pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("pi_pt"+tag+cut+weight+runPeriod_,";#pi^{#pm} P_{T} [GeV];Events / 5 GeV", 10, 0,50);
    allPlots["MET"+tag+cut+weight+runPeriod_] = new TH1F("MET"+tag+cut+weight+runPeriod_,";MET [GeV];Events / 20 GeV", 10,0,200);
    allPlots["HT"+tag+cut+weight+runPeriod_] = new TH1F("HT"+tag+cut+weight+runPeriod_,";HT [GeV];Events / 20 GeV", 10,0,200);
    allPlots["ST"+tag+cut+weight+runPeriod_] = new TH1F("ST"+tag+cut+weight+runPeriod_,";ST [GeV];Events / 20 GeV", 10,0,200);
    allPlots["MET2oST"+tag+cut+weight+runPeriod_] = new TH1F("MET2oST"+tag+cut+weight+runPeriod_,";MET2oST [GeV];Events / 20 GeV", 10,0,200);
    allPlots["charge"+tag+cut+weight+runPeriod_] = new TH1F("charge"+tag+cut+weight+runPeriod_,";Charge(l_{1}*l_{2});Events", 5,-2,2);
    allPlots["csv"+tag+cut+weight+runPeriod_] = new TH1F("CSV"+tag+cut+weight+runPeriod_,";Jet CSV;Events / 0.1", 10,0,1);
    allPlots["dR"+tag+cut+weight+runPeriod_] = new TH1F("dR"+tag+cut+weight+runPeriod_,";dR;Events / 0.05", 20,0.0,1.);
    allPlots["pflp_pt"+tag+cut+weight+runPeriod_] = new TH1F("pflp_pt"+tag+cut+weight+runPeriod_,";PF lepton P_{T} [GeV];Events / 0.2 GeV", 15, 0,3);
    allPlots["massZ"+tag+cut+weight+runPeriod_]     = new TH1F("massZ_control"+tag+cut+weight+runPeriod_,";M_{ll};Events / 1.0 GeV" ,30,81,111);
    allPlots["chargeZ"+tag+cut+weight+runPeriod_]     = new TH1F("chargeZ_control"+tag+cut+weight+runPeriod_,";M_{ll};Events / 1.0 GeV" ,5,-2,2);
    allPlots["nevt"+tag+cut+weight+runPeriod_]     = new TH1F("nevt"+tag+cut+weight+runPeriod_,";N_{events};Events" ,1,1.,2.);
    allPlots["weight"+tag+cut+weight+runPeriod_]     = new TH1F("weight"+tag+cut+weight+runPeriod_,";N_{events};Events/ 1.0" ,20,0.,2.);
    allPlots["norm"+tag+cut+weight+runPeriod_]     = new TH1F("norm"+tag+cut+weight+runPeriod_,";N_{events};Events / 1.0" ,2,0.,2.);
    allPlots["relIso"+tag+cut+weight+runPeriod_] = new TH1F("relIso"+tag+cut+weight+runPeriod_,";relIso;Events / 0.01", 25,0,0.25);
    allPlots["nvtx"+tag+cut+weight+runPeriod_]     = new TH1F("nvtx"+tag+cut+weight+runPeriod_,";N_{PV};Events / 1.0" ,50,0.,50.);
  }
  }
  }
  /*
  h_particles = new TH1F("h_"+name+"_N"  , "Number of "+name,   20,0,20);
  h_pt        = new TH1F("h_"+name+"_pt" , name+" p_{t} (GeV)", 100,0,1000);
  h_eta       = new TH1F("h_"+name+"_eta", name+" #eta",        100,-5,5);
  h_phi       = new TH1F("h_"+name+"_phi", name+" #phi",        100,-3.1415926535,3.1415926535);
  */

}

/*
StdPlots::~StdPlots() {
  delete h_particles;
  delete h_pt;
  delete h_eta;
  delete h_phi;
}
*/

StdPlots::~StdPlots() {};

void StdPlots::SetNorm(float norm) {
  if(debug_) std::cout << "Setting norm= " << norm << std::endl;
  norm_ = norm;
}

void StdPlots::SetSFs(float sfs) {
  if(debug_) std::cout << "Setting SFs= " << sfs << std::endl;
  sfs_ = sfs;
}

void StdPlots::SetPuWgt(float puWgt) {
  if(debug_) std::cout << "Setting puWgt= " << puWgt << std::endl;
  puWgt_ = puWgt;
}

void StdPlots::Fill(double nevt, double nvtx, double HT, double ST, double MET, TString chTag, TString name) {
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = chTag+"_"+name;
  else name = chTag;
  name += runPeriod_;
  if(debug_) std::cout << "Filling " << name << " with wgt=" << wgt << std::endl;

  allPlots["nevt"+name]->Fill(nevt,norm_);
  //allPlots["weight"+name]->Fill(wgt,norm_);
  //allPlots["norm_"+name]->Fill(norm_,norm_);
  allPlots["nvtx"+name]->Fill(nvtx,wgt);

  allPlots["HT"+name]->Fill(HT,wgt);
  allPlots["ST"+name]->Fill(ST,wgt);
  allPlots["MET2oST"+name]->Fill(pow(MET,2)/ST,wgt);
}

void StdPlots::Fill(Leptons leptons, TString chTag, TString name) {
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = chTag+"_"+name;
  else name = chTag;
  name += runPeriod_;
  if(debug_) std::cout << "Filling " << name << " with wgt=" << wgt << std::endl;

  allPlots["nlp"+name]->Fill(leptons.size(),wgt);
  //allPlots["nlp"+name+"_no_weight"]->Fill(leptons.size(),norm_);
  if(leptons.size() > 0) {
    allPlots["lp_pt_low"+name]->Fill(leptons[0].Pt(),wgt);
    allPlots["lp_pt"+name]->Fill(leptons[0].Pt(),wgt);
    //allPlots["lp_pt"+name+"_no_weight"]->Fill(leptons[0].Pt(),norm_);
    //allPlots["relIso"+name]->Fill(leptons[0].getRelIso(),wgt);
    allPlots["lp_eta"+name]->Fill(leptons[0].Eta(),wgt);
    //allPlots["lp_eta"+name+"_no_weight"]->Fill(leptons[0].Eta(),norm_);
    allPlots["lp_phi"+name]->Fill(leptons[0].Phi(),wgt);
    //allPlots["lp_phi"+name+"_no_weight"]->Fill(leptons[0].Phi(),norm_);
  }
}

void StdPlots::Fill(std::vector<Jet> lightJetsVec, std::vector<Jet> bJetsVec, std::vector<Jet> allJetsVec, TString name) {
  float wgt = norm_ * sfs_ * puWgt_;

  allPlots["nj"+name]->Fill(allJetsVec.size(),wgt);
  allPlots["nlj"+name]->Fill(lightJetsVec.size(),wgt);
  allPlots["nbj"+name]->Fill(bJetsVec.size(),wgt);
  allPlots["nlj"+name+"_no_weight"]->Fill(lightJetsVec.size(),norm_);
  allPlots["nbj"+name+"_no_weight"]->Fill(bJetsVec.size(),norm_);
  allPlots["nlj_all"]->Fill(lightJetsVec.size(),wgt);
  allPlots["nbj_all"]->Fill(bJetsVec.size(),wgt);

  if(allJetsVec.size() > 0) {
    allPlots["j_pt"+name]->Fill(allJetsVec[0].getVec().Pt(),wgt);
    allPlots["j_pt_low"+name]->Fill(allJetsVec[0].getVec().Pt(),wgt);
    allPlots["j_pt"+name+"_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm_);
    allPlots["j_pt_all"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
  }
  if(lightJetsVec.size() > 0) {
    allPlots["lj_pt_low"+name]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
    allPlots["lj_pt"+name]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
    allPlots["lj_pt"+name+"_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm_);
    allPlots["lj_pt_all"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
  }
  if(bJetsVec.size() > 0) {
    allPlots["bj_pt_low"+name]->Fill(bJetsVec[0].getVec().Pt(),wgt);
    allPlots["bj_pt"+name]->Fill(bJetsVec[0].getVec().Pt(),wgt);
    allPlots["bj_pt"+name+"_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm_);
    allPlots["bj_pt_all"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
  }

}

void StdPlots::Fill(Double_t weight, Int_t N, Double_t pt, Double_t eta, Double_t phi) {
  h_particles->Fill(N,weight/N);
  h_pt->Fill(pt,weight);
  h_eta->Fill(eta,weight);
  h_phi->Fill(phi,weight);
}

void StdPlots::Write() {
  if(debug_) std::cout << "writing histograms" << std::endl;

  for (auto& it : allPlots)  { 
    if(debug_) std::cout << it.second->GetName() << std::endl;
    if(debug_) std::cout << it.second->GetEntries() << std::endl;

    //fOut->cd( dir );
    //it.second->SetDirectory(fOut); 
    it.second->Write(); 
  }
 /*
 h_particles->Write();
 h_pt->Write();
 h_eta->Write();
 h_phi->Write();
 */
}

