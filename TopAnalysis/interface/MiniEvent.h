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

  Int_t nl;
  Bool_t isPromptFinalState[50], isDirectPromptTauDecayProductFinalState[50];
  Int_t l_id[50],l_charge[50],l_pid[50];
  Float_t l_pt[50],l_eta[50],l_phi[50], l_mass[50], l_miniIso[50], l_chargedHadronIso[50], l_relIso[50], l_ip3d[50], l_ip3dsig[50];

  Int_t nj,ngenj;
  Float_t j_pt[200],j_eta[200],j_phi[200],j_mass[200],j_area[200],j_rawsf[200];
  Float_t genj_pt[200],genj_eta[200],genj_phi[200],genj_mass[200];
  Float_t j_csv[200],j_vtxmass[200],j_vtx3DVal[200],j_vtx3DSig[200],j_puid[200],j_qg[200],j_vtxpx[200],j_vtxpy[200],j_vtxpz[200];
  Int_t j_vtxNtracks[200],j_flav[200],j_pid[200],j_hadflav[200];

  Int_t npf,pf_j[5000];
  Int_t pf_id[5000],pf_charge[5000];
  Float_t pf_px[5000],pf_py[5000],pf_pz[5000];

  Int_t ngen,g_j[5000];
  Int_t g_id[5000],g_charge[5000];
  Float_t g_px[5000],g_py[5000],g_pz[5000];

  Int_t ngenHardProc;
  Int_t ghp_id[100];
  Float_t ghp_pt[100],ghp_eta[100],ghp_phi[100],ghp_m[100];

  Int_t nmet;
  Float_t met_pt[10],met_phi[10];
};

void createMiniEventTree(TTree *t,MiniEvent_t &ev);
void attachToMiniEventTree(TTree *t, MiniEvent_t &ev);

#endif
