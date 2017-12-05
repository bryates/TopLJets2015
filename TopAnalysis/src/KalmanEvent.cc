#include <TTree.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <iostream>
#include "TopLJets2015/TopAnalysis/interface/Jet.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/KalmanEvent.h"

KalmanEvent::KalmanEvent(bool debug) {
  debug_ = debug;
  row_ = 0;
  njpsi_ = 0;
  nmeson_ = 0;
}

KalmanEvent::~KalmanEvent() { ev_ = {}; jets_.clear(); }

void KalmanEvent::loadEvent(const MiniEvent_t &ev) {
  ev_ = ev;
  if(debug_ && ev_.nmeson>0) std::cout << "KalmanEvent loading event: " << ev.event << std::endl;
  njpsi_ = ev_.njpsi;
  nmeson_ = ev_.nmeson;
  vtxProb_ = 0.02;
  chi2_ = 5.; //Same as Elvire's, chi2=5.365 at vtxProb>0.02
  l3dsig_ = 10.; //Elvire used 20 but prompt becomes ~1% at L3D=0.01 in https://byates.web.cern.ch/byates/Top2016/2016/test/JPsi/L3D/l3d_ratio_B.png
  csv_ = 0.5426;
  //csv_ = 0.8484;
  //buildJets();
  if(nmeson_) buildJets();
}

void KalmanEvent::buildJets() {
  jets_.clear(); //start off fresh
  for(int ij=0; ij<ev_.nj; ij++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev_.j_pt[ij],ev_.j_eta[ij],ev_.j_phi[ij],ev_.j_mass[ij]);
    Jet tmpj(jp4, ev_.j_csv[ij], ij, ev_.j_pt_charged[ij], ev_.j_pz_charged[ij], ev_.j_p_charged[ij], ev_.j_pt_pf[ij], ev_.j_pz_pf[ij], ev_.j_p_pf[ij], ev_.j_g[ij]); //Store pt of charged and total PF tracijs and gen matched index
    if(debug_) std::cout << "jet pT=" << tmpj.getPt() << std::endl;
    for(int ipf = 0; ipf < ev_.nkpf; ipf++) {
      if(ev_.k_j[ipf] != ij) continue; //skip if PF track doesn't belong to current jet
      //if(ev_.k_vtxProb[ipf]<vtxProb_) continue;
      if(ev_.k_chi2[ipf]>chi2_) continue;
      if(debug_) std::cout << "passed chi^2 < " << chi2_ << std::endl; 
      //if(ev_.k_mass[ipf]<2.5 || ev_.k_mass[ipf]>3.4) continue;
      if(!ev_.k_mass[ipf]) continue;
      if(debug_) std::cout << "passed mass window" << std::endl; 
      //if(ev_.k_l3d[ipf]/ev_.k_sigmal3d[ipf] < l3dsig_) continue; //proper decay length significance > 20
      //if(debug_) std::cout << "passed l3d/sigmal3d < " << l3dsig_ << std::endl; 
      //testing CSV
      //if(ev_.j_csv[ev_.k_j[ipf]]<csv_) continue;
      TLorentzVector tkP4(0,0,0,0);
      tkP4.SetPtEtaPhiM(ev_.k_pf_pt[ipf],ev_.k_pf_eta[ipf],ev_.k_pf_phi[ipf],ev_.k_pf_m[ipf]);
      pfTrack pftk(tkP4, ev_.k_mass[ipf], ev_.k_l3d[ipf], ev_.k_sigmal3d[ipf], ev_.k_chi2[ipf], ev_.k_vtxProb[ipf], ev_.k_pf_id[ipf], ev_.k_id[ipf]);
      if(debug_) { std::cout << "pfTrack "; pftk.print(); }
      if(debug_) std::cout << "Kalman jet " << ev_.k_j[ipf] << " with pT=" << ev_.j_pt[ev_.k_j[ipf]] << std::endl;
      tmpj.addTrack(pftk);
      //std::cout << std::endl << ev_.event << ": " << pftk.Pt() << " " << pftk.Eta() << " " << pftk.Phi() << " " << ev_.k_mass[ipf] << std::endl;
      //std::cout << ev_.event << ": " << ev_.k_pf_pt[ipf] << " " << ev_.k_pf_eta[ipf] << " " << ev_.k_pf_phi[ipf] << " " << ev_.k_mass[ipf] << std::endl;
    }
    if(!tmpj.getTracks().size()) continue; //skip empty jets
    tmpj.sortTracksByPt();
    jets_.push_back( tmpj );
  }
}

std::vector<Jet> KalmanEvent::getJets()
{
  return jets_;
}
  
bool KalmanEvent::isGoodJet(int idx) {
  //buildJets();
  if(getNMeson()==0) return false;
  for(auto &jet : jets_) {
    if(jet.getJetIndex()==idx) return true;
  }
  return false;
}

