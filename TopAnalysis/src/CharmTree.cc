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
  sfs_ = std::pair<float,float>(1.,0.);
  puWgt_ = 1.;
  top_pt_wgt_ = 1.;
  tracker_wgt_ = 1.;
  pi_wgt_ = std::pair<float,float>(1.,0.);
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

void CharmTree::SetSFs(float sfs, float unc) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting SFs= " << sfs << std::endl;
  sfs_.first = sfs;
  sfs_.second= unc;
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

void CharmTree::SetTrackerWgt(float tracker_wgt) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting SFs= " << tracker_wgt << std::endl;
  tracker_wgt_ = tracker_wgt;
}

void CharmTree::SetPiWgt(float pi_wgt, float unc) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting SFs= " << pi_wgt << std::endl;
  pi_wgt_.first = pi_wgt;
  pi_wgt_.second = unc;
}

void CharmTree::SetLumi(float lumi) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting lumi= " << lumi << std::endl;
  lumi_ = lumi;
  //ev_.lumi = lumi_;
}

void CharmTree::SetXsec(float xsec) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting xsec= " << xsec << std::endl;
  xsec_ = xsec;
  //ev_.xsec = xsec_;
}

void CharmTree::Fill(CharmEvent_t &ev_, double nvtx, double HT, double ST, double MET, std::vector<Jet> lightJets) {
  if(!isGood_) return;
  ev_.ht = HT;
  ev_.met = MET;
  ev_.st = ST;
  ev_.nlj = lightJets.size();
  ev_.nvtx = nvtx;
}

void CharmTree::CheckRunPeriod(TString runPeriod) {
  if(name_.Contains("Data")) return;
  runPeriod_ = runPeriod;
}

void CharmTree::Fill(CharmEvent_t &ev_, std::vector<pfTrack>& pfCands, Leptons lep, Jet jet, TString, TString name, int event, std::vector<pfTrack> genMatch, std::vector<float> frag, Jet genJet) {
  if(!isGood_) return;
  //attachToCharmEventTree(t_,ev_);
  //Fill(pfCands, lep, chTag, name); //Fill meson+lep plots
  //Fill(pfCands, jet, chTag, name); //Fill meson+jet run2
  if(!name.EqualTo("")) name = "_" + name;
  if(name.Contains("jpsi")) {
    TLorentzVector jpsi = pfCands[0].getVec() + pfCands[1].getVec();
    if(jpsi.M()<2.5 || jpsi.M()>3.4) return; //Loose window for mass resonance
    ev_.jpsi_mass[ev_.nmeson] = jpsi.M();
    ev_.meson_id[ev_.nmeson] = abs(pfCands[0].getMotherId());

    int epoch(0);
    if(runPeriod_.Contains("BCDEF"))
      epoch = 1;
    else if(runPeriod_.Contains("GH"))
      epoch = 2;

    if(event>0) ev_.event = event;
    ev_.epoch[ev_.nmeson] = epoch;
    ev_.norm = norm_;
    ev_.xsec = xsec_;
    ev_.puwgt[ev_.nmeson] = puWgt_;
    ev_.topptwgt = top_pt_wgt_;
    ev_.sfs[ev_.nmeson] = sfs_.first;
    ev_.sfsu[ev_.nmeson] = sfs_.second;

    /*
    ev_.peterson[ev_.nj] = frag[0];
    ev_.up[ev_.nj] = frag[1];
    ev_.central[ev_.nj] = frag[2];
    ev_.down[ev_.nj] = frag[3];
    */
    ev_.up[ev_.nj] = frag[0];
    ev_.central[ev_.nj] = frag[1];
    ev_.down[ev_.nj] = frag[2];

    //if(jpsi.M()<3.0 || jpsi.M()>3.2) ev_.nmeson++;
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


    ev_.jpsi_pt[ev_.nmeson] = jpsi.Pt();
    ev_.jpsi_eta[ev_.nmeson] = jpsi.Eta();
    ev_.jpsi_phi[ev_.nmeson] = jpsi.Phi();
    ev_.jpsi_p[ev_.nmeson] = jpsi.P();
    ev_.jpsi_pz[ev_.nmeson] = jpsi.Pz();
    ev_.jpsi_j[ev_.nmeson] = ev_.nj;
    ev_.jpsi_ptrel[ev_.nmeson] = ROOT::Math::VectorUtil::Perp(jpsi.Vect(),jet.getVec().Vect());

    ev_.jpsi_l[ev_.nmeson] = ilep;
    ev_.jpsi_l_mass[ev_.nmeson] = mass123;
    ev_.jpsi_l_dR[ev_.nmeson] = dRJPsil;
    ev_.jpsi_l3d[ev_.nmeson] = pfCands[0].getL3D();
    ev_.jpsi_sigmal3d[ev_.nmeson] = pfCands[0].getSigmaL3D();

    ev_.jpsi_mu1_pt[ev_.nmeson] = pfCands[0].Pt();
    ev_.jpsi_mu1_eta[ev_.nmeson] = pfCands[0].Eta();
    ev_.jpsi_mu1_phi[ev_.nmeson] = pfCands[0].Phi();
    ev_.jpsi_mu2_pt[ev_.nmeson] = pfCands[1].Pt();
    ev_.jpsi_mu2_eta[ev_.nmeson] = pfCands[1].Eta();
    ev_.jpsi_mu2_phi[ev_.nmeson] = pfCands[1].Phi();

    ev_.jpsi_chi2[ev_.nmeson] = pfCands[0].chi2();

    ev_.j_pt[ev_.nj] = jet.getPt();
    ev_.j_pt_charged[ev_.nj] = jet.getChargedPt();
    ev_.j_pt_pf[ev_.nj] = jet.getPFPt();
    ev_.j_p[ev_.nj] = jet.getP();
    ev_.j_p_charged[ev_.nj] = jet.getChargedP();
    ev_.j_p_pf[ev_.nj] = jet.getPFP();
    ev_.j_pz[ev_.nj] = jet.getPz();
    ev_.j_pz_charged[ev_.nj] = jet.getChargedPz();
    ev_.j_pz_pf[ev_.nj] = jet.getPFPz();
    ev_.j_csv[ev_.nj] = jet.getCSV();
    ev_.j_hadflav[ev_.nj] = jet.getHadFlav();
    ev_.nmeson++;
    ev_.nj++;
  }

  else if(name.Contains("meson")) {
    if(pfCands.size()<2) return;
    TLorentzVector D0 = pfCands[0].getVec() + pfCands[1].getVec();
    if(D0.M()<1.7 || D0.M()>2.0) return; //Loose window for mass resonance
    ev_.d0_mass[ev_.nmeson] = D0.M();
    ev_.meson_id[ev_.nmeson] = abs(pfCands[0].getMotherId());
    if(pfCands.size()>2 && abs(pfCands[2].getPdgId())==13) {
      ev_.meson_id[ev_.nmeson] = 42113;
      //if(D0.M()>1.8 && D0.M()<1.93) {
        ev_.d0_mu_pt[ev_.nmeson] = pfCands[2].Pt();
        ev_.d0_mu_eta[ev_.nmeson] = pfCands[2].Eta();
        ev_.d0_mu_phi[ev_.nmeson] = pfCands[2].Phi();
        ev_.d0_mu_tag_mu_pt[ev_.nmeson] = (D0+pfCands[2].getVec()).Pt();
      //}
    }
    if(pfCands.size()>2 && abs(pfCands[2].getPdgId())==211) {
      float mass123 = (D0+pfCands[2].getVec()).M();
      float deltam = mass123 - D0.M();
      if(deltam<0.14 || deltam>0.16) return;
      ev_.ds_mass[ev_.nmeson] = mass123;
      ev_.ds_pi2_pt[ev_.nmeson] = pfCands[2].Pt();
      ev_.ds_pi2_eta[ev_.nmeson] = pfCands[2].Eta();
      ev_.ds_pi2_phi[ev_.nmeson] = pfCands[2].Phi();
    }
    /*
    if(abs(pfCands[0].getMotherId())==42113) {
      ev_.meson_id[ev_.nmeson] = 42113;
    }
    */

    int epoch(0);
    if(runPeriod_.Contains("BCDEF"))
      epoch = 1;
    else if(runPeriod_.Contains("GH"))
      epoch = 2;

    if(event>0) ev_.event = event;
    ev_.epoch[ev_.nmeson] = epoch;
    ev_.norm = norm_;
    ev_.lumi = lumi_;
    ev_.puwgt[ev_.nmeson] = puWgt_;
    ev_.topptwgt = top_pt_wgt_;
    ev_.sfs[ev_.nmeson] = sfs_.first*tracker_wgt_*pi_wgt_.first;
    ev_.sfsu[ev_.nmeson] = sqrt(pow(sfs_.second,2)+pow(pi_wgt_.second,2));

    /*
    ev_.peterson[ev_.nj] = frag[0];
    ev_.central[ev_.nj] = frag[2];
    */
    ev_.up[ev_.nj] = frag[0];
    ev_.down[ev_.nj] = frag[1];
    /*
    ev_.gj_pt[ev_.nj] = genJet.Pt();
    ev_.gj_eta[ev_.nj] = genJet.getVec().Eta();
    ev_.gj_phi[ev_.nj] = genJet.getVec().Phi();
    */

    float mass123((D0 + lep[0].getVec()).M());
    int ilep = 0;
    float dRJPsil(D0.DeltaR(lep[ilep].getVec()));
    //Find closest isolated lepton in di-lepton event
    if(lep.size()>1) {
      for(ilep = 0; ilep < (int)lep.size(); ilep++) {
        float tmpdRJPsil(D0.DeltaR(lep[ilep].getVec()));
        if(tmpdRJPsil < dRJPsil) {
          dRJPsil = tmpdRJPsil;
          mass123 = (D0 + lep[ilep].getVec()).M();
        }
      }
    }


    ev_.d0_pt[ev_.nmeson] = D0.Pt();
    ev_.d0_eta[ev_.nmeson] = D0.Eta();
    ev_.d0_phi[ev_.nmeson] = D0.Phi();
    ev_.d0_p[ev_.nmeson] = D0.P();
    ev_.d0_pz[ev_.nmeson] = D0.Pz();
    ev_.d0_j[ev_.nmeson] = ev_.nj;
    ev_.d0_ptrel[ev_.nmeson] = ROOT::Math::VectorUtil::Perp(D0.Vect(),jet.getVec().Vect());

    ev_.d0_l[ev_.nmeson] = ilep;
    ev_.d0_l_mass[ev_.nmeson] = mass123;
    ev_.d0_l_dR[ev_.nmeson] = dRJPsil;
    ev_.d0_l3d[ev_.nmeson] = pfCands[0].getL3D();
    ev_.d0_lx[ev_.nmeson] = pfCands[0].getLX();
    ev_.d0_ly[ev_.nmeson] = pfCands[0].getLY();
    ev_.d0_lz[ev_.nmeson] = pfCands[0].getLZ();
    ev_.d0_sigmal3d[ev_.nmeson] = pfCands[0].getSigmaL3D();
    ev_.d0_sigmax[ev_.nmeson] = pfCands[0].getSigmaX();
    ev_.d0_sigmay[ev_.nmeson] = pfCands[0].getSigmaY();
    ev_.d0_sigmaz[ev_.nmeson] = pfCands[0].getSigmaZ();

    ev_.d0_pi_pt[ev_.nmeson] = pfCands[0].Pt();
    ev_.d0_pi_eta[ev_.nmeson] = pfCands[0].Eta();
    ev_.d0_pi_phi[ev_.nmeson] = pfCands[0].Phi();
    ev_.d0_pi_dxy[ev_.nmeson] = pfCands[0].getDxy();
    ev_.d0_pi_dxyE[ev_.nmeson] = pfCands[0].getDxyE();
    ev_.d0_pi_dz[ev_.nmeson] = pfCands[0].getDz();
    ev_.d0_k_pt[ev_.nmeson] = pfCands[1].Pt();
    ev_.d0_k_eta[ev_.nmeson] = pfCands[1].Eta();
    ev_.d0_k_phi[ev_.nmeson] = pfCands[1].Phi();
    ev_.d0_k_dxy[ev_.nmeson] = pfCands[1].getDxy();
    ev_.d0_k_dxyE[ev_.nmeson] = pfCands[1].getDxyE();
    ev_.d0_k_dz[ev_.nmeson] = pfCands[1].getDz();

    ev_.d0_chi2[ev_.nmeson] = pfCands[0].chi2();

    if(genMatch.size()>0) {
      ev_.d0_pi_mother[ev_.nmeson] = genMatch[0].getMotherId();
      ev_.d0_k_mother[ev_.nmeson]  = genMatch[1].getMotherId();
    }

    ev_.j_pt[ev_.nj] = jet.getPt();
    ev_.j_eta[ev_.nj] = genJet.getVec().Eta();
    ev_.j_phi[ev_.nj] = genJet.getVec().Phi();
    ev_.j_pt_charged[ev_.nj] = jet.getChargedPt();
    ev_.j_pt_pf[ev_.nj] = jet.getPFPt();
    ev_.j_p[ev_.nj] = jet.getP();
    ev_.j_p_charged[ev_.nj] = jet.getChargedP();
    ev_.j_p_pf[ev_.nj] = jet.getPFP();
    ev_.j_pz[ev_.nj] = jet.getPz();
    ev_.j_pz_charged[ev_.nj] = jet.getChargedPz();
    ev_.j_pz_pf[ev_.nj] = jet.getPFPz();
    ev_.j_csv[ev_.nj] = jet.getCSV();
    ev_.j_hadflav[ev_.nj] = jet.getHadFlav();
    ev_.j_ntk[ev_.nj] = jet.getTracks().size();
    ev_.nmeson++;
    ev_.nj++;
    //std::cout << "tree done" << std::endl;
  }

}

void CharmTree::Write() {
  if(!isGood_) return;
  if(debug_) std::cout << "writing tree" << std::endl;

  //t_->Fill();
}

