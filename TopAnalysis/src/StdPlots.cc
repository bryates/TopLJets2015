#include <iostream>
#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"
#include "TopLJets2015/TopAnalysis/interface/StdPlots.h"
#include "TopLJets2015/TopAnalysis/interface/Jet.h"

StdPlots::StdPlots(TString runPeriod, TString name, bool debug) {  
  runPeriod_ = runPeriod;
  if(name.Contains("MC") || (name.Contains("Data") && name.Contains("2016"+runPeriod_))) isGood_ = true;
  else isGood_ = false;
  debug_ = debug;
  norm_ = 1.;
  sfs_ = 1.;
  puWgt_ = 1.;

  if(debug_ && isGood_)
    std::cout << "Initializing run" << runPeriod_ << std::endl;

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
    allPlots["lp_pt_iso"+tag+cut+weight+"_"+runPeriod_] = new TH1F("lp_pt_iso"+tag+cut+weight+"_"+runPeriod_,";Lepton P_{T} [GeV] after cleaning;Events / 10 GeV", 20, 0,200);
    allPlots["lp_pt_veto"+tag+cut+weight+"_"+runPeriod_] = new TH1F("lp_pt_veto"+tag+cut+weight+"_"+runPeriod_,";Lepton P_{T} [GeV] after veto;Events / 10 GeV", 20, 0,200);
    allPlots["lp_pt_low"+tag+cut+weight+"_"+runPeriod_] = new TH1F("lp_pt_low"+tag+cut+weight+"_"+runPeriod_,";Leading Lepton P_{T} [GeV];Events / 1 GeV", 20, 20,40);
    allPlots["lp_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("lp_pt"+tag+cut+weight+"_"+runPeriod_,";Leading lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["l2p_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("l2p_pt"+tag+cut+weight+"_"+runPeriod_,";Sub-leading lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["lp_eta"+tag+cut+weight+"_"+runPeriod_]  = new TH1F("lp_eta"+tag+cut+weight+"_"+runPeriod_,";Leading lepton #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["l2p_eta"+tag+cut+weight+"_"+runPeriod_]  = new TH1F("l2p_eta"+tag+cut+weight+"_"+runPeriod_,";Sub-Leading lepton #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["lp_phi"+tag+cut+weight+"_"+runPeriod_]  = new TH1F("lp_phi"+tag+cut+weight+"_"+runPeriod_,";Leading lepton #phi; Events", 50, -3.14,3.14);
    allPlots["l2p_phi"+tag+cut+weight+"_"+runPeriod_]  = new TH1F("l2p_phi"+tag+cut+weight+"_"+runPeriod_,";Sub-Leading lepton #phi; Events", 50, -3.14,3.14);
    allPlots["dilp_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("dilp_pt"+tag+cut+weight+"_"+runPeriod_,";Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_m"+tag+cut+weight+"_"+runPeriod_] = new TH1F("dilp_m"+tag+cut+weight+"_"+runPeriod_,";M_{ll} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["ndilp"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("ndilp"+tag+cut+weight+"_"+runPeriod_,";N_{ll};Events" ,3,0.,3.);

    //Jet plots
    allPlots["j_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("j_pt"+tag+cut+weight+"_"+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["lj_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("lj_pt"+tag+cut+weight+"_"+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["bj_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("bj_pt"+tag+cut+weight+"_"+runPeriod_,";Leading b Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["j_pt_low"+tag+cut+weight+"_"+runPeriod_] = new TH1F("j_pt_low"+tag+cut+weight+"_"+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["lj_pt_low"+tag+cut+weight+"_"+runPeriod_] = new TH1F("lj_pt_low"+tag+cut+weight+"_"+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["bj_pt_low"+tag+cut+weight+"_"+runPeriod_] = new TH1F("bj_pt_low"+tag+cut+weight+"_"+runPeriod_,";Leading b Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["nlp"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("nlp"+tag+cut+weight+"_"+runPeriod_,";N_{l};Events" ,3,0.,3.);
    allPlots["nj"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("nj"+tag+cut+weight+"_"+runPeriod_,";N_{jets} (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["nlj"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("nlj"+tag+cut+weight+"_"+runPeriod_,";N_{jets} (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["nbj"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("nbj"+tag+cut+weight+"_"+runPeriod_,";N_{b-jets} (CSV > 0.8);Events" ,4,1.,5.);

    allPlots["npf"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("npf"+tag+cut+weight+"_"+runPeriod_,";N_{pf};Events / 10" ,5,0.,5.);
    allPlots["nstart"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("jetindex"+tag+cut+weight+"_"+runPeriod_,";N_{jetindex};Events" ,5,0.,5.);
    allPlots["pfid"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("pfid"+tag+cut+weight+"_"+runPeriod_,";PFID;Events" ,440,-220.,220.);
    
    //J/Psi plots
    allPlots["massJPsi"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massJPsi"+tag+cut+weight+"_"+runPeriod_,";M_{ll};Events / 18 MeV" ,50,2.5,3.4);
    allPlots["massJPsi_mu"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massJPsi_mu"+tag+cut+weight+"_"+runPeriod_,";M_{J/#Psi+#mu};Events / 10 GeV" ,25,0,250);
    allPlots["massJPsi_e"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massJPsi_e"+tag+cut+weight+"_"+runPeriod_,";M_{J/#Psi+e};Events / 10 GeV" ,25,0,250);
    allPlots["massJPsi_l"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massJPsi_l"+tag+cut+weight+"_"+runPeriod_,";M_{J/#Psi+l};Events / 10 GeV" ,25,0,250);
    allPlots["massJPsiK"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massJPsiK"+tag+cut+weight+"_"+runPeriod_,";M_{llk};Events / 15 MeV" ,100,4.5,6);
    allPlots["JPsi_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("JPsi_pt"+tag+cut+weight+"_"+runPeriod_,";P_{T}(J/#Psi) [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["JPsi_mu_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("JPsi_mu_pt"+tag+cut+weight+"_"+runPeriod_,";P_{T}(J/#Psi+#mu) [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["JPsi_e_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("JPsi_e_pt"+tag+cut+weight+"_"+runPeriod_,";P_{T}(J/#Psi+e) [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["JPsi_l_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("JPsi_l_pt"+tag+cut+weight+"_"+runPeriod_,";P_{T}(J/#Psi+l) [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["JPsioJet_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("JPsioJet_pt"+tag+cut+weight+"_"+runPeriod_,";P_{T}(J/#Psi)/P_{T}(jet) [GeV];Events / 0.05 GeV", 10, 0,1);
    allPlots["dR_JPsi_mu"+tag+cut+weight+"_"+runPeriod_] = new TH1F("dR_JPsi_mu"+tag+cut+weight+"_"+runPeriod_,";#DeltaR(J/#Psi-leading #mu) [GeV];Events / 0.1 GeV", 10, 0,5);
    allPlots["dR_JPsi_e"+tag+cut+weight+"_"+runPeriod_] = new TH1F("dR_JPsi_e"+tag+cut+weight+"_"+runPeriod_,";#DeltaR(J/#Psi-leading e) [GeV];Events / 0.1 GeV", 10, 0,5);
    allPlots["dR_JPsi_l"+tag+cut+weight+"_"+runPeriod_] = new TH1F("dR_JPsi_l"+tag+cut+weight+"_"+runPeriod_,";#DeltaR(J/#Psi-leading l) [GeV];Events / 0.1 GeV", 10, 0,5);
 
    //D meson plots
    allPlots["massD0"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massD0"+tag+cut+weight+"_"+runPeriod_,";M_{D^{0}};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_lep"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massD0_lep"+tag+cut+weight+"_"+runPeriod_,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_mu"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massD0_mu"+tag+cut+weight+"_"+runPeriod_,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_e"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massD0_ele"+tag+cut+weight+"_"+runPeriod_,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massDsmD0loose"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massDsmD0loose"+tag+cut+weight+"_"+runPeriod_,";M_{K#pi#pi} - M_{K#pi};Events / 0.5 MeV" ,20,0.14,0.16);
    allPlots["massDsmD0"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massDsmD0"+tag+cut+weight+"_"+runPeriod_,";M_{K#pi#pi} - M_{K#pi};Events / 0.5 MeV" ,20,0.14,0.16);
    allPlots["massDs"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massDs"+tag+cut+weight+"_"+runPeriod_,";M_{D^{*}};Events / 10 MeV" ,200,0.,2.0);

    allPlots["pi_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("pi_pt"+tag+cut+weight+"_"+runPeriod_,";#pi^{#pm} P_{T} [GeV];Events / 5 GeV", 10, 0,50);

    allPlots["MET"+tag+cut+weight+"_"+runPeriod_] = new TH1F("MET"+tag+cut+weight+"_"+runPeriod_,";MET [GeV];Events / 20 GeV", 10,0,200);
    allPlots["HT"+tag+cut+weight+"_"+runPeriod_] = new TH1F("HT"+tag+cut+weight+"_"+runPeriod_,";HT [GeV];Events / 20 GeV", 10,0,200);
    allPlots["ST"+tag+cut+weight+"_"+runPeriod_] = new TH1F("ST"+tag+cut+weight+"_"+runPeriod_,";ST [GeV];Events / 20 GeV", 10,0,200);
    allPlots["MET2oST"+tag+cut+weight+"_"+runPeriod_] = new TH1F("MET2oST"+tag+cut+weight+"_"+runPeriod_,";MET2oST [GeV];Events / 20 GeV", 10,0,200);
    allPlots["charge"+tag+cut+weight+"_"+runPeriod_] = new TH1F("charge"+tag+cut+weight+"_"+runPeriod_,";Charge(l_{1}*l_{2});Events", 5,-2,2);
    allPlots["csv"+tag+cut+weight+"_"+runPeriod_] = new TH1F("CSV"+tag+cut+weight+"_"+runPeriod_,";Jet CSV;Events / 0.1", 10,0,1);
    allPlots["dR"+tag+cut+weight+"_"+runPeriod_] = new TH1F("dR"+tag+cut+weight+"_"+runPeriod_,";dR;Events / 0.05", 20,0.0,1.);
    allPlots["pflp_pt"+tag+cut+weight+"_"+runPeriod_] = new TH1F("pflp_pt"+tag+cut+weight+"_"+runPeriod_,";PF lepton P_{T} [GeV];Events / 0.2 GeV", 15, 0,3);

    //Z control plots
    allPlots["massZ"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("massZ_control"+tag+cut+weight+"_"+runPeriod_,";M_{ll};Events / 1.0 GeV" ,30,81,111);
    allPlots["chargeZ"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("chargeZ_control"+tag+cut+weight+"_"+runPeriod_,";M_{ll};Events / 1.0 GeV" ,5,-2,2);

    //Event plots
    allPlots["nevt"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("nevt"+tag+cut+weight+"_"+runPeriod_,";N_{events};Events" ,1,1.,2.);
    allPlots["weight"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("weight"+tag+cut+weight+"_"+runPeriod_,";N_{events};Events/ 1.0" ,20,0.,2.);
    allPlots["norm"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("norm"+tag+cut+weight+"_"+runPeriod_,";N_{events};Events / 1.0" ,2,0.,2.);
    allPlots["relIso"+tag+cut+weight+"_"+runPeriod_] = new TH1F("relIso"+tag+cut+weight+"_"+runPeriod_,";relIso;Events / 0.01", 25,0,0.25);
    allPlots["nvtx"+tag+cut+weight+"_"+runPeriod_]     = new TH1F("nvtx"+tag+cut+weight+"_"+runPeriod_,";N_{PV};Events / 1.0" ,50,0.,50.);
  }
  }
  }

}

StdPlots::~StdPlots() {
  for (auto& it : allPlots) delete it.second;
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
}

void StdPlots::Fill(double nevt, double nvtx, double HT, double ST, double MET, TString chTag, TString name) {
  if(!isGood_) return;
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = chTag+"_"+name;
  else name = chTag;
  if(debug_) std::cout << "Filling " << name << "_" << runPeriod_ << " with wgt=" << wgt << std::endl;

  allPlots["nevt"+name+"_"+runPeriod_]->Fill(nevt,norm_);
  //allPlots["weight"+name+"_"+runPeriod_]->Fill(wgt,norm_);
  //allPlots["norm"+name+"_"+runPeriod_]->Fill(norm_,norm_);
  allPlots["nvtx"+name+"_"+runPeriod_]->Fill(nvtx,wgt);

  allPlots["HT"+name+"_"+runPeriod_]->Fill(HT,wgt);
  allPlots["ST"+name+"_"+runPeriod_]->Fill(ST,wgt);
  allPlots["MET2oST"+name+"_"+runPeriod_]->Fill(pow(MET,2)/ST,wgt);
}

void StdPlots::Fill(Leptons leptons, TString chTag, TString name) {
  if(!isGood_) return;
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = chTag+"_"+name;
  else name = chTag;
  if(debug_) std::cout << "Filling " << name << "_" << runPeriod_ << " with wgt=" << wgt << std::endl;

  allPlots["nlp"+name+"_"+runPeriod_]->Fill(leptons.size(),wgt);
  //allPlots["nlp"+name+"_no_weight"+"_"+runPeiod_]->Fill(leptons.size(),norm_);
  if(leptons.size() > 0) {
    allPlots["lp_pt_low"+name+"_"+runPeriod_]->Fill(leptons[0].Pt(),wgt);
    allPlots["lp_pt"+name+"_"+runPeriod_]->Fill(leptons[0].Pt(),wgt);
    //allPlots["lp_pt"+name+"_"+"_no_weight"+runPeiod_]->Fill(leptons[0].Pt(),norm_);
    //allPlots["relIso"+name+"_"+runPeriod_]->Fill(leptons[0].getRelIso(),wgt);
    allPlots["lp_eta"+name+"_"+runPeriod_]->Fill(leptons[0].Eta(),wgt);
    //allPlots["lp_eta"+name+"_"+"_no_weight"+"_"+runPeiod_]->Fill(leptons[0].Eta(),norm_);
    allPlots["lp_phi"+name+"_"+runPeriod_]->Fill(leptons[0].Phi(),wgt);
    //allPlots["lp_phi"+name+"_"+"_no_weight"+"_"+runPeiod_]->Fill(leptons[0].Phi(),norm_);
  }
}

void StdPlots::Fill(std::vector<Jet> lightJetsVec, std::vector<Jet> bJetsVec, std::vector<Jet> allJetsVec, TString chTag, TString name) {
  if(!isGood_) return;
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = chTag+"_"+name;
  else name = chTag;
  if(debug_) std::cout << "Filling " << name << "_" << runPeriod_ << " with wgt=" << wgt << std::endl;

  allPlots["nj"+name+"_"+runPeriod_]->Fill(allJetsVec.size(),wgt);
  allPlots["nlj"+name+"_"+runPeriod_]->Fill(lightJetsVec.size(),wgt);
  allPlots["nbj"+name+"_"+runPeriod_]->Fill(bJetsVec.size(),wgt);
  //allPlots["nlj"+name+"_no_weight"+"_"+runPeriod_]->Fill(lightJetsVec.size(),norm_);
  //allPlots["nbj"+name+"_no_weight"+"_"+runPeriod_]->Fill(bJetsVec.size(),norm_);
  //allPlots["nlj_all"]->Fill(lightJetsVec.size(),wgt);
  //allPlots["nbj_all"]->Fill(bJetsVec.size(),wgt);

  if(allJetsVec.size() > 0) {
    allPlots["j_pt"+name+"_"+runPeriod_]->Fill(allJetsVec[0].getVec().Pt(),wgt);
    allPlots["j_pt_low"+name+"_"+runPeriod_]->Fill(allJetsVec[0].getVec().Pt(),wgt);
    //allPlots["j_pt"+name+"_no_weight"+"_"+runPeriod_]->Fill(allJetsVec[0].getVec().Pt(), norm_);
    //allPlots["j_pt_all"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
    for(size_t ij = 0; ij < allJetsVec.size(); ij++) {
      float csv = allJetsVec.at(ij).getCSV();
      allPlots["csv"+name+"_"+runPeriod_]->Fill(csv,wgt);
    }
    /*
    for(auto it : allJetsVec)
      allPlots["csv"+name+"_"+runPeriod_]->Fill(it.getCSV(),wgt);
    */
  }
  if(lightJetsVec.size() > 0) {
    allPlots["lj_pt_low"+name+"_"+runPeriod_]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
    allPlots["lj_pt"+name+"_"+runPeriod_]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
    //allPlots["lj_pt"+name+"_no_weight"+"_"+runPeriod_]->Fill(lightJetsVec[0].getVec().Pt(),norm_);
    //allPlots["lj_pt_all"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
  }
  if(bJetsVec.size() > 0) {
    allPlots["bj_pt_low"+name+"_"+runPeriod_]->Fill(bJetsVec[0].getVec().Pt(),wgt);
    allPlots["bj_pt"+name+"_"+runPeriod_]->Fill(bJetsVec[0].getVec().Pt(),wgt);
    //allPlots["bj_pt"+name+"_no_weight"+"_"+runPeriod_]->Fill(bJetsVec[0].getVec().Pt(),norm_);
    //allPlots["bj_pt_all"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
  }

}

void StdPlots::Fill(std::vector<pfTrack> pfmuCands, TString chTag, TString name) {
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = chTag+"_"+name;
  else name = chTag;
  if(debug_) std::cout << "Filling " << name << "_" << runPeriod_ << " with wgt=" << wgt << std::endl;
  if(name.Contains("jpsi") && pfmuCands[0].getPfid() == -pfmuCands[1].getPfid()) {
    float mass12((pfmuCands[0].getVec() + pfmuCands[1].getVec()).M());
    float pt12((pfmuCands[0].getVec() + pfmuCands[1].getVec()).Pt());
    //float mass123( kaonCands.size()>0 ? (pfmuCands[0].getVec()+pfmuCands[1].getVec()+kaonCands[0].getVec()).M() : -1);

    allPlots["massJPsi"+name+"_"+runPeriod_]->Fill(mass12,wgt);
    allPlots["JPsi_pt"+name+"_"+runPeriod_]->Fill(pt12,wgt);
    //allPlots["massJPsi_all"]->Fill(mass12,wgt);
  }
}

void StdPlots::Fill(std::vector<pfTrack> pfmuCands, Leptons lep, TString chTag, TString name) {
  Fill(pfmuCands, chTag, name);
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = chTag+"_"+name;
  else name = chTag;
  if(debug_) std::cout << "Filling " << name << "_" << runPeriod_ << " with wgt=" << wgt << std::endl;
  if(name.Contains("jpsi") && pfmuCands[0].getPfid() == -pfmuCands[1].getPfid()) {
    float pt123((pfmuCands[0].getVec() + pfmuCands[1].getVec()+lep[0].getVec()).Pt());
    float mass123((pfmuCands[0].getVec() + pfmuCands[1].getVec()+lep[0].getVec()).M());
    float dRJPsil((pfmuCands[0].getVec() + pfmuCands[1].getVec()).DeltaR(lep[0].getVec()));
    if(lep[0].isMuon()) {
      allPlots["JPsi_mu_pt"+name+"_"+runPeriod_]->Fill(pt123,wgt);
      allPlots["massJPsi_mu"+name+"_"+runPeriod_]->Fill(mass123,wgt);
      allPlots["dR_JPsi_mu"+name+"_"+runPeriod_]->Fill(dRJPsil,wgt);
    }
    else if(lep[0].isElectron()) {
      allPlots["JPsi_e_pt"+name+"_"+runPeriod_]->Fill(pt123,wgt);
      allPlots["massJPsi_e"+name+"_"+runPeriod_]->Fill(mass123,wgt);
      allPlots["dR_JPsi_mu"+name+"_"+runPeriod_]->Fill(dRJPsil,wgt);
      allPlots["dR_JPsi_e"+name+"_"+runPeriod_]->Fill(dRJPsil,wgt);
    }
    allPlots["JPsi_l_pt"+name+"_"+runPeriod_]->Fill(pt123,wgt);
    allPlots["massJPsi_l"+name+"_"+runPeriod_]->Fill(mass123,wgt);
    allPlots["dR_JPsi_l"+name+"_"+runPeriod_]->Fill(dRJPsil,wgt);
  }
}

void StdPlots::Fill(std::vector<pfTrack> pfmuCands, Jet jet, TString chTag, TString name) {
  float wgt = norm_ * sfs_ * puWgt_;
  if(!name.EqualTo("")) name = chTag+"_"+name;
  else name = chTag;
  if(debug_) std::cout << "Filling " << name << "_" << runPeriod_ << " with wgt=" << wgt << std::endl;
  if(name.Contains("jpsi") && pfmuCands[0].getPfid() == -pfmuCands[1].getPfid()) {
    float mass12((pfmuCands[0].getVec() + pfmuCands[1].getVec()).M());
    if(mass12<3.0 || mass12>3.2) return;
    float pt12((pfmuCands[0].getVec() + pfmuCands[1].getVec()).Pt());
    float jpt(jet.getPt());
    allPlots["JPsioJet_pt"+name+"_"+runPeriod_]->Fill(pt12/jpt,wgt);
  }
}

void StdPlots::Write() {
  if(!isGood_) return;
  if(debug_) std::cout << "writing histograms" << std::endl;

  for (auto& it : allPlots)  { 
    if(debug_) std::cout << it.second->GetName() << std::endl;
    if(debug_) std::cout << it.second->GetEntries() << std::endl;

    it.second->Write(); 
  }
}

