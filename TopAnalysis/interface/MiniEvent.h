#ifndef _minievent_h_
#define _minievent_h_

#include "TTree.h"

struct MiniEvent_t
{
  Int_t run,event,lumi; 
  Int_t nvtx,pu;
  Bool_t isFiducial;
  Int_t l_id,l_charge, ttbar_nw, ttbar_allmepartons, ttbar_matchmepartons;
  Float_t l_pt,l_eta,l_phi,l_chargedHadronIso,l_neutralHadronIso,l_photonIso,l_puChargedHadronIso, l_tmass;
  Float_t rho;
  Float_t ttbar_w[500];

  Int_t nj;
  Float_t j_pt[1000],j_eta[1000],j_phi[1000],j_csv[1000],j_vtxmass[1000],j_vtx3DVal[1000],j_vtx3DSig[1000],j_puid[1000], j_chpt[1000],j_cheta[1000],j_chphi[1000], j_mass[1000];
  Float_t genj_pt[1000],genj_eta[1000],genj_phi[1000],genj_mass[1000];
  
  Int_t j_vtxNtracks[1000],j_flav[1000],j_pid[1000];
  Float_t met_pt,met_phi,mt;
  Float_t chmet_pt,chmet_phi,chmt;
  Bool_t isData, isMC;

};

void createMiniEventTree(TTree *t,MiniEvent_t &ev);
void attachToMiniEventTree(TTree *t, MiniEvent_t &ev);

#endif
