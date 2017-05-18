#include <iostream>
#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"
#include "TopLJets2015/TopAnalysis/interface/StdPlots.h"
#include "TopLJets2015/TopAnalysis/interface/Jet.h"

StdPlots::StdPlots(TString runPeriod, TString name, bool debug) {  
  runPeriod_ = "_"+runPeriod;
  name_ = name;
  //Check for multiple run periods
  if(runPeriod.Length() > 1) {
    for(int i=0; i < runPeriod.Length(); i++) {
      TString tmp(runPeriod[i]);
      if(name.Contains("MC") || (name.Contains("Data") && name.Contains("2016"+tmp))) {
        isGood_ = true;
        break;
      }
    }
  }
  //Only use MC or Data with correct run period
  else if(name.Contains("MC") || (name.Contains("Data") && name.Contains("2016"+runPeriod))) isGood_ = true;
  else isGood_ = false;
  debug_ = debug;
  norm_ = 1.;
  sfs_ = 1.;
  puWgt_ = 1.;

  if(debug_ && isGood_)
    std::cout << "Initializing run" << runPeriod_ << std::endl;

  //PU plot
  allPlots["puwgtctr"+runPeriod_] = new TH1F("puwgtctr"+runPeriod_,"Weight sums",4,0,4);
  std::vector<TString> lfsVec = { "_all", "_e", "_ee", "_em", "_mm", "_m" }; 
  std::vector<TString> cutVec = { "", "_lep", "_lepjets", "_jpsi", "_csv", "_meson" };
  std::vector<TString> wgtVec = { "", "_no_weight" };
  for(int i = 0; i < (int)lfsVec.size(); i++) {
  for(int j = 0; j < (int)cutVec.size(); j++) {
  for(int k = 0; k < (int)wgtVec.size(); k++) {
    TString tag(lfsVec[i]);
    TString cut(cutVec[j]);
    TString weight(wgtVec[k]);

    // Lepton plots
    allPlots["lp_pt_iso"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt_iso"+tag+cut+weight+runPeriod_,";Lepton P_{T} [GeV] after cleaning;Events / 10 GeV", 20, 0,200);
    allPlots["lp_pt_veto"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt_veto"+tag+cut+weight+runPeriod_,";Lepton P_{T} [GeV] after veto;Events / 10 GeV", 20, 0,200);
    allPlots["lp_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt_low"+tag+cut+weight+runPeriod_,";Leading Lepton P_{T} [GeV];Events / 0.5 GeV", 20, 20,40);
    allPlots["lp_pt"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt"+tag+cut+weight+runPeriod_,";Leading lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["l2p_pt"+tag+cut+weight+runPeriod_] = new TH1F("l2p_pt"+tag+cut+weight+runPeriod_,";Sub-leading lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["lp_eta"+tag+cut+weight+runPeriod_]  = new TH1F("lp_eta"+tag+cut+weight+runPeriod_,";Leading lepton #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["l2p_eta"+tag+cut+weight+runPeriod_]  = new TH1F("l2p_eta"+tag+cut+weight+runPeriod_,";Sub-Leading lepton #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["lp_phi"+tag+cut+weight+runPeriod_]  = new TH1F("lp_phi"+tag+cut+weight+runPeriod_,";Leading lepton #phi; Events", 50, -3.14,3.14);
    allPlots["l2p_phi"+tag+cut+weight+runPeriod_]  = new TH1F("l2p_phi"+tag+cut+weight+runPeriod_,";Sub-Leading lepton #phi; Events", 50, -3.14,3.14);
    allPlots["dilp_pt"+tag+cut+weight+runPeriod_] = new TH1F("dilp_pt"+tag+cut+weight+runPeriod_,";Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_m"+tag+cut+weight+runPeriod_] = new TH1F("dilp_m"+tag+cut+weight+runPeriod_,";M_{ll} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["ndilp"+tag+cut+weight+runPeriod_]     = new TH1F("ndilp"+tag+cut+weight+runPeriod_,";Di-Lepton Multiplicity;Events" ,3,0.,3.);

    //Jet plots
    allPlots["j_pt"+tag+cut+weight+runPeriod_] = new TH1F("j_pt"+tag+cut+weight+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["lj_pt"+tag+cut+weight+runPeriod_] = new TH1F("lj_pt"+tag+cut+weight+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["bj_pt"+tag+cut+weight+runPeriod_] = new TH1F("bj_pt"+tag+cut+weight+runPeriod_,";Leading b-Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["j_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("j_pt_low"+tag+cut+weight+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 0.5 GeV", 20, 30,50);
    allPlots["lj_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("lj_pt_low"+tag+cut+weight+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 0.5 GeV", 20, 30,50);
    allPlots["bj_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("bj_pt_low"+tag+cut+weight+runPeriod_,";Leading b-Jet P_{T} [GeV];Events / 0.5 GeV", 20, 30,50);
    allPlots["nlp"+tag+cut+weight+runPeriod_]     = new TH1F("nlp"+tag+cut+weight+runPeriod_,";Lepton Multiplicity;Events" ,3,0.,3.);
    allPlots["nj"+tag+cut+weight+runPeriod_]     = new TH1F("nj"+tag+cut+weight+runPeriod_,";Jet Multiplicity (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["nlj"+tag+cut+weight+runPeriod_]     = new TH1F("nlj"+tag+cut+weight+runPeriod_,";Light-Jet Multiplicity (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["nbj"+tag+cut+weight+runPeriod_]     = new TH1F("nbj"+tag+cut+weight+runPeriod_,";b-Jet Multiplicity (CSV > 0.8);Events" ,4,1.,5.);

    
    //J/Psi plots
    allPlots["massJPsi"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi"+tag+cut+weight+runPeriod_,";M_{#mu^{#pm}#mu^{#mp}} [GeV];Events / 18 MeV" ,50,2.5,3.4);
    allPlots["massJPsi_mu"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_mu"+tag+cut+weight+runPeriod_,";M_{J/#Psi+#mu} [GeV];Events / 10 GeV" ,25,0,250);
    allPlots["massJPsi_e"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_e"+tag+cut+weight+runPeriod_,";M_{J/#Psi+e} [GeV];Events / 10 GeV" ,25,0,250);
    allPlots["massJPsi_l"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_l"+tag+cut+weight+runPeriod_,";M_{J/#Psi+l} [GeV];Events / 10 GeV" ,25,0,250);
    allPlots["massJPsiK"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsiK"+tag+cut+weight+runPeriod_,";M_{llk} [GeV];Events / 15 MeV" ,100,4.5,6);
    allPlots["JPsi_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_pt"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi) [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["JPsi_eta"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_eta"+tag+cut+weight+runPeriod_,";J/#Psi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["JPsi_phi"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_phi"+tag+cut+weight+runPeriod_,";J/#Psi #phi; Events", 50, -3.14,3.14);
    allPlots["JPsi_mu_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu_pt"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+#mu) [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["JPsi_e_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_e_pt"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+e) [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["JPsi_l_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_l_pt"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+l) [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["JPsi_mu1_eta"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu1_eta"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{1} #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["JPsi_mu2_eta"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu2_eta"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{2} #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["JPsi_mu1_phi"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu1_phi"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{1} #phi; Events", 50, -3.14,3.14);
    allPlots["JPsi_mu2_phi"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu2_phi"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{2} #phi; Events", 50, -3.14,3.14);
    allPlots["JPsioJet_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi)/P_{T}(jet);Events / 0.02", 10, 0,1);
    allPlots["dR_JPsi_mu"+tag+cut+weight+runPeriod_] = new TH1F("dR_JPsi_mu"+tag+cut+weight+runPeriod_,";#DeltaR(J/#Psi-leading #mu);Events / 0.01", 10, 0,5);
    allPlots["dR_JPsi_e"+tag+cut+weight+runPeriod_] = new TH1F("dR_JPsi_e"+tag+cut+weight+runPeriod_,";#DeltaR(J/#Psi-leading e);Events / 0.01", 10, 0,5);
    allPlots["dR_JPsi_l"+tag+cut+weight+runPeriod_] = new TH1F("dR_JPsi_l"+tag+cut+weight+runPeriod_,";#DeltaR(J/#Psi-leading l);Events / 0.01", 10, 0,5);
    allPlots["dR_JPsi_mumu"+tag+cut+weight+runPeriod_] = new TH1F("dR_JPsi_mumu"+tag+cut+weight+runPeriod_,";#DeltaR(J/#Psi_{#mu1},J/#Psi_{#mu2});Events / 0.01", 10, 0,0.5);
    allPlots["JPsi_mumu_pt_ratio"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mumu_pt_ratio"+tag+cut+weight+runPeriod_,";(#mu_{1 P_{T}}-#mu_{2 P_{T}})/(#mu_{1 P_{T}}+#mu_{2 P_{T}});Events / 0.01", 10, 0,1);
 
    //D meson plots
    allPlots["massD0"+tag+cut+weight+runPeriod_]     = new TH1F("massD0"+tag+cut+weight+runPeriod_,";M_{D^{0}};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_lep"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_lep"+tag+cut+weight+runPeriod_,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_mu"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_mu"+tag+cut+weight+runPeriod_,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_e"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_ele"+tag+cut+weight+runPeriod_,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massDsmD0loose"+tag+cut+weight+runPeriod_]     = new TH1F("massDsmD0loose"+tag+cut+weight+runPeriod_,";M_{K#pi#pi} - M_{K#pi};Events / 0.5 MeV" ,20,0.14,0.16);
    allPlots["massDsmD0"+tag+cut+weight+runPeriod_]     = new TH1F("massDsmD0"+tag+cut+weight+runPeriod_,";M_{K#pi#pi} - M_{K#pi};Events / 0.5 MeV" ,20,0.14,0.16);
    allPlots["massDs"+tag+cut+weight+runPeriod_]     = new TH1F("massDs"+tag+cut+weight+runPeriod_,";M_{D^{*}};Events / 10 MeV" ,200,0.,2.0);

    allPlots["pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("pi_pt"+tag+cut+weight+runPeriod_,";#pi^{#pm} P_{T} [GeV];Events / 5 GeV", 10, 0,50);

    allPlots["MET"+tag+cut+weight+runPeriod_] = new TH1F("MET"+tag+cut+weight+runPeriod_,";E^{miss}_{T} [GeV];Events / 20 GeV", 10,0,200);
    allPlots["HT"+tag+cut+weight+runPeriod_] = new TH1F("HT"+tag+cut+weight+runPeriod_,";H_{T} [GeV];Events / 20 GeV", 50,0,1000);
    allPlots["ST"+tag+cut+weight+runPeriod_] = new TH1F("ST"+tag+cut+weight+runPeriod_,";S_{T} [GeV];Events / 20 GeV", 50,0,1000);
    allPlots["MET2oST"+tag+cut+weight+runPeriod_] = new TH1F("MET2oST"+tag+cut+weight+runPeriod_,";MET^{2}/ST [GeV];Events / 20 GeV", 10,0,200);
    allPlots["charge"+tag+cut+weight+runPeriod_] = new TH1F("charge"+tag+cut+weight+runPeriod_,";Charge(l_{1}*l_{2});Events", 5,-2,2);
    allPlots["csv"+tag+cut+weight+runPeriod_] = new TH1F("CSV"+tag+cut+weight+runPeriod_,";Jet CSV;Events / 0.1", 10,0,1);
    allPlots["dR"+tag+cut+weight+runPeriod_] = new TH1F("dR"+tag+cut+weight+runPeriod_,";dR;Events / 0.05", 20,0.0,1.);
    allPlots["pflp_pt"+tag+cut+weight+runPeriod_] = new TH1F("pflp_pt"+tag+cut+weight+runPeriod_,";PF lepton P_{T} [GeV];Events / 0.2 GeV", 15, 0,3);

    //Z control plots
    allPlots["massZ"+tag+cut+weight+runPeriod_]     = new TH1F("massZ_control"+tag+cut+weight+runPeriod_,";M_{l^{#pm}l^{#mp}};Events / 1.0 GeV" ,30,81,111);
    allPlots["chargeZ"+tag+cut+weight+runPeriod_]     = new TH1F("chargeZ_control"+tag+cut+weight+runPeriod_,";Charage (l^#pm) * Charge(l^#mp);Events / 1.0 GeV" ,5,-2,2);

    //Event plots
allPlots["nevt"+tag+cut+weight+runPeriod_]     = new TH1F("nevt"+tag+cut+weight+runPeriod_,";Event Multiplicity;Events" ,1,1.,2.);
allPlots["weight"+tag+cut+weight+runPeriod_]     = new TH1F("weight"+tag+cut+weight+runPeriod_,";weights;Events/ 1.0" ,20,0.,2.);
allPlots["norm"+tag+cut+weight+runPeriod_]     = new TH1F("norm"+tag+cut+weight+runPeriod_,";norm;Events / 1.0" ,2,0.,2.);
allPlots["relIso"+tag+cut+weight+runPeriod_] = new TH1F("relIso"+tag+cut+weight+runPeriod_,";relIso;Events / 0.01", 25,0,0.25);
allPlots["nvtx"+tag+cut+weight+runPeriod_]     = new TH1F("nvtx"+tag+cut+weight+runPeriod_,";Vertex Multiplicity;Events / 1.0" ,50,0.,50.);
  }
  }
  }

}

StdPlots::~StdPlots() {
  //for (auto it : allPlots) delete it.second;
}

void StdPlots::SetNorm(float norm) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting norm= " << norm << std::endl;
  norm_ = norm;
}

void StdPlots::SetSFs(float sfs) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting SFs= " << sfs << std::endl;
  sfs_ = sfs;
}

void StdPlots::SetPuWgt(float puWgt) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting puWgt= " << puWgt << std::endl;
  puWgt_ = puWgt;
  allPlots["puwgtctr"+runPeriod_]->Fill(0.,1.0);
  allPlots["puwgtctr"+runPeriod_]->Fill(1,puWgt_);
}

void StdPlots::Fill(double nevt, double nvtx, double HT, double ST, double MET, TString chTag, TString name) {
  if(!isGood_) return;
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = "_" + name;
  if(debug_) std::cout << "Filling nvtx" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;

  allPlots["nevt"+chTag+name+runPeriod_]->Fill(nevt,wgt);
  //allPlots["weight"+chTag+name+runPeriod_]->Fill(wgt,norm_);
  //allPlots["norm"+chTag+name+runPeriod_]->Fill(norm_,norm_);
  allPlots["nvtx"+chTag+name+runPeriod_]->Fill(nvtx,wgt);

  allPlots["HT"+chTag+name+runPeriod_]->Fill(HT,wgt);
  allPlots["ST"+chTag+name+runPeriod_]->Fill(ST,wgt);
  allPlots["MET2oST"+chTag+name+runPeriod_]->Fill(pow(MET,2)/ST,wgt);
}

void StdPlots::Fill(Leptons leptons, TString chTag, TString name) {
  if(!isGood_) return;
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = "_" + name;
  if(debug_) std::cout << "Filling lep" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;

  allPlots["nlp"+chTag+name+runPeriod_]->Fill(leptons.size(),wgt);
  //allPlots["nlp"+chTag+name+"_no_weight"+"_"+runPeiod_]->Fill(leptons.size(),norm_);
  if(leptons.size() > 0) {
    allPlots["lp_pt_low"+chTag+name+runPeriod_]->Fill(leptons[0].Pt(),wgt);
    allPlots["lp_pt"+chTag+name+runPeriod_]->Fill(leptons[0].Pt(),wgt);
    //allPlots["lp_pt"+chTag+name+"_"+"_no_weight"+runPeiod_]->Fill(leptons[0].Pt(),norm_);
    //allPlots["relIso"+chTag+name+runPeriod_]->Fill(leptons[0].getRelIso(),wgt);
    allPlots["lp_eta"+chTag+name+runPeriod_]->Fill(leptons[0].Eta(),wgt);
    //allPlots["lp_eta"+chTag+name+"_"+"_no_weight"+"_"+runPeiod_]->Fill(leptons[0].Eta(),norm_);
    allPlots["lp_phi"+chTag+name+runPeriod_]->Fill(leptons[0].Phi(),wgt);
    //allPlots["lp_phi"+chTag+name+"_"+"_no_weight"+"_"+runPeiod_]->Fill(leptons[0].Phi(),norm_);
  }
}

void StdPlots::Fill(std::vector<Jet> lightJetsVec, std::vector<Jet> bJetsVec, std::vector<Jet> allJetsVec, TString chTag, TString name) {
  if(!isGood_) return;
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = "_" + name;
  if(debug_) std::cout << "Filling jet" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;

  allPlots["nj"+chTag+name+runPeriod_]->Fill(allJetsVec.size(),wgt);
  allPlots["nlj"+chTag+name+runPeriod_]->Fill(lightJetsVec.size(),wgt);
  allPlots["nbj"+chTag+name+runPeriod_]->Fill(bJetsVec.size(),wgt);
  //allPlots["nlj"+chTag+name+"_no_weight"+runPeriod_]->Fill(lightJetsVec.size(),norm_);
  //allPlots["nbj"+chTag+name+"_no_weight"+runPeriod_]->Fill(bJetsVec.size(),norm_);
  //allPlots["nlj_all"+name+runPeriod_]->Fill(lightJetsVec.size(),wgt);
  //allPlots["nbj_all"+name+runPeriod_]->Fill(bJetsVec.size(),wgt);

  if(allJetsVec.size() > 0) {
    allPlots["j_pt"+chTag+name+runPeriod_]->Fill(allJetsVec[0].getVec().Pt(),wgt);
    allPlots["j_pt_low"+chTag+name+runPeriod_]->Fill(allJetsVec[0].getVec().Pt(),wgt);
    //allPlots["j_pt"+chTag+name+"_no_weight"+runPeriod_]->Fill(allJetsVec[0].getVec().Pt(), norm_);
    //allPlots["j_pt_all"+name+runPeriod_]->Fill(allJetsVec[0].getVec().Pt(),wgt);
    for(size_t ij = 0; ij < allJetsVec.size(); ij++) {
      float csv = allJetsVec.at(ij).getCSV();
      allPlots["csv"+chTag+name+runPeriod_]->Fill(csv,wgt);
    }
    
    for(auto it : allJetsVec)
      allPlots["csv_all"+name+runPeriod_]->Fill(it.getCSV(),wgt);
  }
  if(lightJetsVec.size() > 0) {
    allPlots["lj_pt_low"+chTag+name+runPeriod_]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
    allPlots["lj_pt"+chTag+name+runPeriod_]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
    //allPlots["lj_pt"+chTag+name+"_no_weight"+runPeriod_]->Fill(lightJetsVec[0].getVec().Pt(),norm_);
    //allPlots["lj_pt_all"+name+runPeriod_]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
  }
  if(bJetsVec.size() > 0) {
    allPlots["bj_pt_low"+chTag+name+runPeriod_]->Fill(bJetsVec[0].getVec().Pt(),wgt);
    allPlots["bj_pt"+chTag+name+runPeriod_]->Fill(bJetsVec[0].getVec().Pt(),wgt);
    //allPlots["bj_pt"+chTag+name+"_no_weight"+runPeriod_]->Fill(bJetsVec[0].getVec().Pt(),norm_);
    //allPlots["bj_pt_all"+name+runPeriod_]->Fill(bJetsVec[0].getVec().Pt(),wgt);
  }

}
//Called by Fill(std::vector<pfTrack> pfmuCands, Leptons lep, TString chTag, TString name)
void StdPlots::Fill(std::vector<pfTrack> pfmuCands, TString chTag, TString name) {
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = "_" + name;
  if(debug_) std::cout << "Filling J/Psi" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;
  if(name.Contains("jpsi")) {
    if(debug_) std::cout << "is " << name << std::endl;
    if(pfmuCands[0].getPfid() != -pfmuCands[1].getPfid()) return;
    TLorentzVector jpsi = pfmuCands[0].getVec() + pfmuCands[1].getVec();
    //float mass12((pfmuCands[0].getVec() + pfmuCands[1].getVec()).M());
    //J/Psi mass in slightly wide window
    if(jpsi.M()<2.5 || jpsi.M()>3.4) return; //Loose window for mass only
    if(debug_) std::cout << "is J/Psi" << name << std::endl;
    allPlots["massJPsi"+chTag+name+runPeriod_]->Fill(jpsi.M(),wgt);
    allPlots["massJPsi_all"+name+runPeriod_]->Fill(jpsi.M(),wgt);
    //float pt12((pfmuCands[0].getVec() + pfmuCands[1].getVec()).Pt());

    if(jpsi.M()<3.0 || jpsi.M()>3.2) return; //Window in Elvire's AN
    allPlots["JPsi_pt"+chTag+name+runPeriod_]->Fill(jpsi.Pt(),wgt);
    allPlots["JPsi_pt_all"+name+runPeriod_]->Fill(jpsi.Pt(),wgt);
    allPlots["JPsi_eta"+chTag+name+runPeriod_]->Fill(jpsi.Eta(),wgt);
    allPlots["JPsi_eta_all"+name+runPeriod_]->Fill(jpsi.Eta(),wgt);
    allPlots["JPsi_phi"+chTag+name+runPeriod_]->Fill(jpsi.Phi(),wgt);
    allPlots["JPsi_phi_all"+name+runPeriod_]->Fill(jpsi.Phi(),wgt);
  }
}

//Called by Fill(std::vector<pfTrack> pfmuCands, Leptons lep, Jet jet, TString chTag, TString name)
void StdPlots::Fill(std::vector<pfTrack> pfmuCands, Leptons lep, TString chTag, TString name) {
  Fill(pfmuCands, chTag, name); //Fill meson only plots
  Fill(lep, chTag, name);       //Fill lepton only plots
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = "_" + name;
  if(debug_) std::cout << "Filling J/Psi+lep" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;
  //J/Psi events
  if(name.Contains("jpsi")) {
    if(debug_) std::cout << "is " << name << std::endl;
    if(pfmuCands[0].getPfid() != -pfmuCands[1].getPfid()) return;
    float mass12((pfmuCands[0].getVec() + pfmuCands[1].getVec()).M());
    //J/Psi mass in slightly wide window
    if(mass12<3.0 || mass12>3.2) return; //Window in Elvire's AN
    if(debug_) std::cout << "is J/Psi" << name << std::endl;
    float pt123((pfmuCands[0].getVec() + pfmuCands[1].getVec()+lep[0].getVec()).Pt());
    float mass123((pfmuCands[0].getVec() + pfmuCands[1].getVec()+lep[0].getVec()).M());
    float dRJPsil((pfmuCands[0].getVec() + pfmuCands[1].getVec()).DeltaR(lep[0].getVec()));
    float dRmumu(pfmuCands[0].getVec().DeltaR(pfmuCands[1].getVec()));
    //J/Psi mass in slightly wide window
    allPlots["massJPsi_l"+chTag+name+runPeriod_]->Fill(mass123,wgt);
    if(lep[0].isMuon())
      allPlots["massJPsi_mu"+chTag+name+runPeriod_]->Fill(mass123,wgt);
    else if(lep[0].isElectron())
      allPlots["massJPsi_e"+chTag+name+runPeriod_]->Fill(mass123,wgt);

    if(lep[0].isMuon()) {
      allPlots["JPsi_mu_pt"+chTag+name+runPeriod_]->Fill(pt123,wgt);
      allPlots["dR_JPsi_mu"+chTag+name+runPeriod_]->Fill(dRJPsil,wgt);
      allPlots["JPsi_mu_pt_all"+name+runPeriod_]->Fill(pt123,wgt);
      allPlots["massJPsi_mu_all"+name+runPeriod_]->Fill(mass123,wgt);
      allPlots["dR_JPsi_mu_all"+name+runPeriod_]->Fill(dRJPsil,wgt);
    }
    else if(lep[0].isElectron()) {
      allPlots["JPsi_e_pt"+chTag+name+runPeriod_]->Fill(pt123,wgt);
      allPlots["dR_JPsi_e"+chTag+name+runPeriod_]->Fill(dRJPsil,wgt);
      allPlots["JPsi_e_pt_all"+name+runPeriod_]->Fill(pt123,wgt);
      allPlots["massJPsi_e_all"+name+runPeriod_]->Fill(mass123,wgt);
      allPlots["dR_JPsi_e_all"+name+runPeriod_]->Fill(dRJPsil,wgt);
    }
    allPlots["JPsi_l_pt"+chTag+name+runPeriod_]->Fill(pt123,wgt);
    allPlots["dR_JPsi_l"+chTag+name+runPeriod_]->Fill(dRJPsil,wgt);
    allPlots["JPsi_l_pt_all"+name+runPeriod_]->Fill(pt123,wgt);
    allPlots["massJPsi_l_all"+name+runPeriod_]->Fill(mass123,wgt);
    allPlots["dR_JPsi_l_all"+name+runPeriod_]->Fill(dRJPsil,wgt);

    allPlots["dR_JPsi_mumu"+chTag+name+runPeriod_]->Fill(dRmumu,wgt);
    allPlots["JPsi_mu1_eta"+chTag+name+runPeriod_]->Fill(pfmuCands[0].Eta(),wgt);
    allPlots["JPsi_mu2_eta"+chTag+name+runPeriod_]->Fill(pfmuCands[1].Eta(),wgt);
    allPlots["JPsi_mu1_phi"+chTag+name+runPeriod_]->Fill(pfmuCands[0].Eta(),wgt);
    allPlots["JPsi_mu2_phi"+chTag+name+runPeriod_]->Fill(pfmuCands[1].Eta(),wgt);
    float pt_ratio = (pfmuCands[0].Pt() - pfmuCands[1].Pt())/(pfmuCands[0].Pt()+ pfmuCands[1].Pt());
    allPlots["JPsi_mumu_pt_ratio"+chTag+name+runPeriod_]->Fill(pt_ratio,wgt);
    allPlots["JPsi_mumu_pt_ratio_all"+name+runPeriod_]->Fill(pt_ratio,wgt);
 
  }
}

//Called by Fill(std::vector<pfTrack> pfmuCands, Leptons lep, Jet jet, TString chTag, TString name)
void StdPlots::Fill(std::vector<pfTrack> pfmuCands, Jet jet, TString chTag, TString name) {
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = "_" + name;
  if(debug_) std::cout << "Filling J/Psi+jet" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;
  if(name.Contains("jpsi")) {
    if(debug_) std::cout << "is " << name << std::endl;
    if(pfmuCands[0].getPfid() != -pfmuCands[1].getPfid()) return;
    float mass12((pfmuCands[0].getVec() + pfmuCands[1].getVec()).M());
    if(mass12<3.0 || mass12>3.2) return; //Window in Elvire's AN
    if(debug_) std::cout << "is J/Psi" << name << std::endl;
    float pt12((pfmuCands[0].getVec() + pfmuCands[1].getVec()).Pt());
    float jpt(jet.getPt());
    allPlots["JPsioJet_pt"+chTag+name+runPeriod_]->Fill(pt12/jpt,wgt);
    allPlots["JPsioJet_pt_all"+name+runPeriod_]->Fill(pt12/jpt,wgt);
    float csv = jet.getCSV();
    allPlots["csv"+chTag+name+runPeriod_]->Fill(csv,wgt);
    allPlots["csv_all"+name+runPeriod_]->Fill(csv,wgt);
  }
}

//Fill for mesons (currently only J/Psi
void StdPlots::Fill(std::vector<pfTrack> pfmuCands, Leptons lep, Jet jet, TString chTag, TString name) {
  Fill(pfmuCands, lep, chTag, name); //Fill meson+lep plots
  Fill(pfmuCands, jet, chTag, name); //Fill meson+jet run2
}

StdPlots Combine(const StdPlots &plots1, const StdPlots &plots2) {
  std::map<TString, float> lumi;
  lumi["all"] = 35740.161;
  lumi["BCDEF"] = 19593.811;
  lumi["GH"] = 16146.17;
  //Load run period from plots
  //e.g. "_GH"
  TString runPeriod(plots2.runPeriod_);
  //Strip the '_'
  //e.g. "GH"
  runPeriod.Replace(0,1,"");
  //if(debug_) std::cout << "Adding plots from runPeriod " << runPeriod << std::endl;
  //Append to current run period
  //e.g. "_BCDEF"->"_BCDEFGH"
  runPeriod = plots1.runPeriod_ + runPeriod;
  runPeriod.Replace(0,1,"");
  if(plots1.debug_) std::cout << "New runPeriod " << runPeriod << std::endl;
  //Add plot to new run
  //e.g. CSV_m_BCDEFGH->Add(CSV_m_BCDEFGH)
  StdPlots run(runPeriod, plots1.name_, plots1.debug_);
  for(auto& it : plots1.allPlots) {
    if(!plots1.isGood_) break;
    //Strip run period from plot being added and append new run period
    if(it.second->Integral() == 0) continue;
    TString plotName(it.first);
    plotName.Resize((it.first.Length() - plots1.runPeriod_.Length()));
    TString tmp(plotName + plots2.runPeriod_);
    plotName += "_"+runPeriod;
    TString tmpl(plots1.runPeriod_);
    tmpl.Replace(0,1,"");
    //1 - wgt due to 0 < Integral < 1
    //e.g. BCDEF has larger lumi, should get larger weight
    //float wgt = it.second->Integral()/(it.second->Integral() + plots2.allPlots.at(tmp)->Integral());
    float wgt = lumi[tmpl]/lumi["all"];
    if(plots2.debug_) std::cout << tmp << " " << wgt << std::endl;
    if(run.name_.Contains("Data")) wgt=1.0;
    it.second->Scale(wgt);
    run.allPlots[plotName]->Add(it.second);
  }
  //Add plot from plots to new
  //e.g. CSV_m_GH->Add(CSV_m_BCDEFGH)
  for(auto& it : plots2.allPlots) {
    if(!plots2.isGood_) break;
    //Strip run period from plot being added and append new run period
    if(it.second->Integral() == 0) continue;
    TString plotName(it.first);
    plotName.Resize((it.first.Length() - plots2.runPeriod_.Length()));
    TString tmp(plotName + plots1.runPeriod_);
    plotName += "_"+runPeriod;
    TString tmpl(plots2.runPeriod_);
    tmpl.Replace(0,1,"");
    //float wgt = it.second->Integral()/(it.second->Integral() + plots1.allPlots.at(tmp)->Integral());
    float wgt = lumi[tmpl]/lumi["all"];
    if(plots1.debug_) std::cout << tmp << " " << wgt << std::endl;
    if(run.name_.Contains("Data")) wgt=1.0;
    if(plots1.debug_) std::cout << it.second->Integral() << std::endl;
    it.second->Scale(wgt);
    run.allPlots[plotName]->Add(it.second);
  }
  return run;
}

//StdPlots operator+(const StdPlots &lhs, const StdPlots &rhs) { return Combine(lhs,rhs); };

void StdPlots::Write() {
  if(!isGood_) return;
  if(debug_) std::cout << "writing histograms" << std::endl;

  for (auto& it : allPlots)  { 
    if(debug_) std::cout << it.second->GetName() << std::endl;
    if(debug_) std::cout << it.second->GetEntries() << std::endl;

    it.second->Write(); 
  }
}

