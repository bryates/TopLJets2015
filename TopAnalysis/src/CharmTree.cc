#include <iostream>
#include "TopLJets2015/TopAnalysis/interface/CharmEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CharmTree.h"
#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"
#include "TopLJets2015/TopAnalysis/interface/Jet.h"

#include "Math/VectorUtil.h"

CharmTree::CharmTree(TTree *t, TString runPeriod, TString name, bool debug) {  
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
  top_pt_wgt_ = 1.;
  //t_ = t;
  //attachToCharmEventTree(t_,ev_);
  if(debug_ && isGood_)
    std::cout << "Initializing run" << runPeriod_ << std::endl;
}

CharmTree::~CharmTree() {
  //for (auto it : allPlots) delete it.second;
}

void CharmTree::SetNorm(float norm) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting norm= " << norm << std::endl;
  norm_ = norm;
  //ev_.norm = norm_;
}

void CharmTree::SetSFs(float sfs) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting SFs= " << sfs << std::endl;
  sfs_ = sfs;
  //ev_.sfs = sfs_;
}

void CharmTree::SetPuWgt(float puWgt) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting puWgt= " << puWgt << std::endl;
  puWgt_ = puWgt;
  //ev_.puwgt = puWgt_;
}

void CharmTree::SetTopPtWgt(float top_pt_wgt) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting top pT weight= " << top_pt_wgt << std::endl;
  top_pt_wgt_ = top_pt_wgt;
  //ev_.topptwgt = top_pt_wgt_;
}

void CharmTree::Fill(CharmEvent_t &ev_, std::vector<pfTrack> &pfCands, Leptons lep, Jet jet, TString chTag, TString name) {
  //attachToCharmEventTree(t_,ev_);
  //Fill(pfCands, lep, chTag, name); //Fill meson+lep plots
  //Fill(pfCands, jet, chTag, name); //Fill meson+jet run2
  if(!name.EqualTo("")) name = "_" + name;
  if(name.Contains("jpsi")) {
    TLorentzVector jpsi = pfCands[0].getVec() + pfCands[1].getVec();
    if(jpsi.M()<2.5 || jpsi.M()>3.4) return; //Loose window for mass resonance
    ev_.jpsi_mass[ev_.njpsi] = jpsi.M();

    int epoch(0);
    if(runPeriod_.Contains("BCDEF"))
      epoch = 1;
    else if(runPeriod_.Contains("GH"))
      epoch = 2;

    ev_.epoch[ev_.njpsi] = epoch;
    ev_.norm = norm_;
    ev_.puwgt[ev_.njpsi] = puWgt_;
    ev_.topptwgt = top_pt_wgt_;
    ev_.sfs[ev_.njpsi] = sfs_;

    //if(jpsi.M()<3.0 || jpsi.M()>3.2) ev_.njpsi++;
    //if(jpsi.M()<3.0 || jpsi.M()>3.2) return; //Window in Elvire's AN
    //float jpt(jet.getPt());
    //float jpt_charged(jet.getChargedPt());
    //float jpt_pf(jet.getPFPt());
    //float dRJPsil(jpsi.DeltaR(lep[0].getVec()));
    float mass123((jpsi + lep[0].getVec()).M());
    int ilep = 0;
    float dRJPsil(jpsi.DeltaR(lep[ilep].getVec()));
    //Find closest isolated lepton in di-lepton event
    if(lep.size()>1) {
      for(ilep = 0; ilep < (int)lep.size(); ilep++) {
        float tmpdRJPsil(jpsi.DeltaR(lep[ilep].getVec()));
        if(tmpdRJPsil < dRJPsil) {
          dRJPsil = tmpdRJPsil;
          mass123 = (jpsi + lep[ilep].getVec()).M();
        }
      }
    }


    ev_.jpsi_pt[ev_.njpsi] = jpsi.Pt();
    ev_.jpsi_eta[ev_.njpsi] = jpsi.Eta();
    ev_.jpsi_phi[ev_.njpsi] = jpsi.Phi();
    ev_.jpsi_p[ev_.njpsi] = jpsi.P();
    ev_.jpsi_pz[ev_.njpsi] = jpsi.Pz();
    ev_.jpsi_j[ev_.njpsi] = ev_.nj;
    ev_.jpsi_ptrel[ev_.njpsi] = ROOT::Math::VectorUtil::Perp(jpsi.Vect(),jet.getVec().Vect());

    ev_.jpsi_l[ev_.njpsi] = ilep;
    ev_.jpsi_l_mass[ev_.njpsi] = mass123;
    ev_.jpsi_l_dR[ev_.njpsi] = dRJPsil;

    ev_.jpsi_mu1_pt[ev_.njpsi] = pfCands[0].Pt();
    ev_.jpsi_mu1_eta[ev_.njpsi] = pfCands[0].Eta();
    ev_.jpsi_mu1_phi[ev_.njpsi] = pfCands[0].Phi();
    ev_.jpsi_mu2_pt[ev_.njpsi] = pfCands[1].Pt();
    ev_.jpsi_mu2_eta[ev_.njpsi] = pfCands[1].Eta();
    ev_.jpsi_mu2_phi[ev_.njpsi] = pfCands[1].Phi();

    ev_.j_pt[ev_.njpsi] = jet.getPt();
    ev_.j_pt_charged[ev_.njpsi] = jet.getChargedPt();
    ev_.j_pt_pf[ev_.njpsi] = jet.getPFPt();
    //ev_.j_p[ev_.njpsi] = jet.getP();
    //ev_.j_p_charged[ev_.njpsi] = jet.getChargedP();
    //ev_.j_p_pf[ev_.njpsi] = jet.getPFP();
    //ev_.j_pz[ev_.njpsi] = jet.getPz();
    //ev_.j_pz_charged[ev_.njpsi] = jet.getChargedPz();
    //ev_.j_pz_pf[ev_.njpsi] = jet.getPFPz();
    ev_.njpsi++;
    ev_.nj++;
  }
}

void CharmTree::Write() {
  if(!isGood_) return;
  if(debug_) std::cout << "writing tree" << std::endl;

  //t_->Fill();
}

