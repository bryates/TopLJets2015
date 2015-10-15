#ifndef _minievent_h_
#define _minievent_h_

#include "TTree.h"

struct MiniEvent_t
{
  Bool_t isData;
  Int_t run,event,lumi;

  Int_t ttbar_nw, ttbar_allmepartons, ttbar_matchmepartons,ttbar_genId;
  Float_t ttbar_w[500];
  Bool_t isFiducial;

  Int_t nvtx,pu,putrue;
  Float_t rho;
  
  Int_t muTrigger,elTrigger;
  
  Bool_t isPromptFinalState, isDirectPromptTauDecayProductFinalState;
  Int_t l_id,l_charge;
  Float_t l_pt,l_eta,l_phi, l_mass, l_chargedHadronIso, l_neutralHadronIso, l_photonIso, l_puChargedHadronIso,l_ip3d,l_ip3dsig;

  Int_t nj,ngenj;
  Float_t j_pt[200],j_eta[200],j_phi[200],j_mass[200],j_area[200];
  Float_t genj_pt[200],genj_eta[200],genj_phi[200],genj_mass[200];
  Float_t j_csv[200],j_vtxmass[200],j_vtx3DVal[200],j_vtx3DSig[200],j_puid[200],j_vtxpx[200],j_vtxpy[200],j_vtxpz[200];
  Int_t j_vtxNtracks[200],j_flav[200],j_pid[200],j_hadflav[200];

  Float_t met_pt,met_phi,mt;

  Int_t me_id,me_np,me_pid[25];
  Float_t me_px[25],me_py[25],me_pz[25],me_mass[25];
};

void createMiniEventTree(TTree *t,MiniEvent_t &ev);
void attachToMiniEventTree(TTree *t, MiniEvent_t &ev);

#endif
