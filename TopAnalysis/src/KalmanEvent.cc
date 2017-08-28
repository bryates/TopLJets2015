#include <TTree.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <iostream>
#include "TopLJets2015/TopAnalysis/interface/Jet.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/KalmanEvent.h"

KalmanEvent::KalmanEvent(TTree *tree) {
  loadTree(tree);
  row_ = 0;
  njpsi_ = 0;
  nmeson_ = 0;
}

KalmanEvent::KalmanEvent(bool debug) {
  debug_ = debug;
  row_ = 0;
  njpsi_ = 0;
  nmeson_ = 0;
}

KalmanEvent::~KalmanEvent() { }
//pass it "kalman/data"
void KalmanEvent::loadTree(TTree *tree) {
  tree_ = tree;
  attachToMiniEventTree(tree_,ev_,true);
}

void KalmanEvent::loadEvent(int event) {
  char selection[500];
  sprintf(selection,"event==%d",event);
  int tmp = tree_->Draw("Entry$",selection,"goff");
  if(!tmp) { njpsi_=0;nmeson_=0; return; }
  event_ = event;
  TH1F *h = (TH1F*)tree_->GetHistogram();
  row_ = h->GetBinLowEdge(h->FindFirstBinAbove(0));
  tree_->GetEntry(row_);
  njpsi_ = ev_.njpsi;
  nmeson_ = ev_.nmeson;
}

void KalmanEvent::loadEvent(const MiniEvent_t &ev) {
  ev_ = ev;
  njpsi_ = ev_.njpsi;
  nmeson_ = ev_.nmeson;
  buildJets();
  vtxProb_ = 0.02;
  chi2_ = 6;
}

void KalmanEvent::buildJets() {
  for(int ij=0; ij<ev_.nj; ij++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev_.j_pt[ij],ev_.j_eta[ij],ev_.j_phi[ij],ev_.j_mass[ij]);
    //jp4.SetPtEtaPhiM(ev_.k_j_pt[ij],ev_.k_j_eta[ij],ev_.k_j_phi[ij],ev_.k_j_mass[ij]);
    Jet tmpj(jp4, 1, ij);
    if(debug_) std::cout << "jet pT=" << tmpj.getPt() << std::endl;
    for(int ipf = 0; ipf < ev_.npf; ipf++) {
      if(ev_.k_j[ipf] != ij) continue; //skip if PF track doesn't belong to current jet
      //if(ev_.k_vtxProb[ipf]<0.02) continue;
      if(ev_.k_vtxProb[ipf]<vtxProb_) continue;
      //if(ev_.k_chi2[ipf]<chi2_) continue;
      if(!ev_.k_mass[ipf]) continue;
      TLorentzVector tkP4(0,0,0,0);
      tkP4.SetPtEtaPhiM(ev_.k_pf_pt[ipf],ev_.k_pf_eta[ipf],ev_.k_pf_phi[ipf],ev_.k_pf_m[ipf]);
      pfTrack pftk(tkP4, ev_.k_chi2[ipf], ev_.k_vtxProb[ipf], ev_.k_pf_id[ipf]);
      if(debug_) { std::cout << "pfTrack "; pftk.print(); }
      tmpj.addTrack(pftk);
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
  buildJets();
  for(auto &jet : jets_) {
    if(jet.getJetIndex()==idx) return true;
  }
  return false;
}

/*
std::vector<Jet> KalmanEvent::getJets() {
  std::vector<Jet> jets;
  for(int ij=0; ij<ev_.nj; ij++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev_.j_pt[ij],ev_.j_eta[ij],ev_.j_phi[ij],ev_.j_mass[ij]);
    //jp4.SetPtEtaPhiM(ev_.k_j_pt[ij],ev_.k_j_eta[ij],ev_.k_j_phi[ij],ev_.k_j_mass[ij]);
    Jet tmpj(jp4, 0, ij);
    for(int ipf = 0; ipf < ev_.npf; ipf++) {
      if(ev_.k_j[ipf] != ij) continue; //skip if PF track doesn't belong to current jet
      if(ev_.k_vtxProb[ipf]<0.02) continue;
      if(!ev_.k_mass[ipf]) continue;
      if(debug_) std::cout << "Kalman jet: " << ev_.k_j[ipf] << " " << ev_.k_pf_pt[ipf] << std::endl;
      //float mass(0);
      //if(abs(ev_.k_id[ipf])==13) mass=0.1057;
      TLorentzVector tkP4(0,0,0,0);
      //tkP4.SetPtEtaPhiM(ev_.pf_pt[ipf],ev_.pf_eta[ipf],ev_.pf_phi[ipf],mass);
      tkP4.SetPtEtaPhiM(ev_.k_pf_pt[ipf],ev_.k_pf_eta[ipf],ev_.k_pf_phi[ipf],ev_.k_pf_m[ipf]);
      //pfTrack pftk(tkP4, ev_.k_chi2[ipf], ev_.k_vtxProb[ipf], ev_.pf_id[ipf]);
      pfTrack pftk(tkP4, ev_.k_chi2[ipf], ev_.k_vtxProb[ipf], ev_.k_pf_id[ipf]);
      if(debug_) pftk.print();
      tmpj.addTrack(pftk);
    }
    tmpj.sortTracksByPt();
    //if(debug_) tmpj.getTracks().back().print();
    jets.push_back( tmpj );
  }
  return jets;
}
*/

/*
          Jet tmpj(jp4, csv, k, ev.j_pt_charged[k], ev.j_pt_pf[k], ev.j_g[k]); //Store pt of charged and total PF tracks and gen matched index
	  for(int ipf = 0; ipf < ev.npf; ipf++) {
	    if(ev.pf_j[ipf] != k) continue; //skip if PF track doesn't belong to current jet
	    if(ev.pf_c[ipf]==0) continue;   //skip if PF track is neutral
	    TLorentzVector tkP4(0,0,0,0);
	    tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);
            pfTrack pftk(tkP4, ev.pf_dxy[ipf], ev.pf_dxyE[ipf], ev.pf_dz[ipf], ev.pf_dzE[ipf], ev.pf_id[ipf],ev.pf_quality[ipf],ev.pf_highPurity[ipf]);
            if(abs(pftk.getPdgId())==13) {
              pftk.setGlobalMuon(ev.pf_globalMuon[ipf]);
              pftk.setTrackerMuon(ev.pf_trackerMuon[ipf]);
            }
	    tmpj.addTrack(pftk); //,ev.pf_id[ipf]);
	  }
          tmpj.sortTracksByPt();
*/
