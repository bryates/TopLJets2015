#ifndef _charmevent_h_
#define _charmevent_h_

#include "TTree.h"

struct CharmEvent_t
{
  CharmEvent_t()
  {
    ttbar_nw=0;
    nl=0; nj=0; npf=0;
    njpsi=0; nmeson=0;
  }

  Bool_t isData;
  Int_t run,event,lumi,epoch[5];
  Float_t norm,puwgt[5],topptwgt;
  Float_t sfs[5]; //need to split into separate parts


  //gen level event
  Int_t ttbar_nw;
  Float_t ttbar_w[500];

  //reco level event
  Int_t nvtx;
  Int_t nl,nleptons;
  Int_t l_id[50],l_charge[50],l_pid[50],l_g[200];
  Float_t l_pt[50],l_eta[50],l_phi[50], l_mass[50];
  Float_t l_chi2norm[50], l_dxy[50], l_dxyE[50], l_dz[50], l_dzE[50];

  Int_t nj;
  Float_t j_pt[200],j_pt_pf[200],j_pt_charged[200],j_eta[200],j_phi[200],j_mass[200];
  Float_t j_p[200],j_p_pf[200],j_p_charged[200];
  Float_t j_pz[200],j_pz_pf[200],j_pz_charged[200];
  Float_t j_csv[200];

  //PF candidates
  Int_t npf,nmu,pf_j[5000];
  Int_t pf_id[5000],pf_c[5000];
  Float_t pf_pt[5000],pf_eta[5000],pf_phi[5000],pf_m[5000],pf_dxy[5000],pf_dxyE[5000],pf_dz[5000],pf_dzE[5000];

  //JPsi candidates
  Int_t njpsi,nmeson;
  Float_t jpsi_mass[5],jpsi_pt[5],jpsi_eta[5],jpsi_phi[5],jpsi_p[5],jpsi_pz[5];
  Float_t jpsi_j[5],jpsi_ptrel[5],jpsi_l[5];//,jpsi_j_dR[5];
  Float_t jpsi_mu1_pt[5],jpsi_mu1_eta[5],jpsi_mu1_phi[5];
  Float_t jpsi_mu2_pt[5],jpsi_mu2_eta[5],jpsi_mu2_phi[5];
  Float_t jpsi_l_mass[5],jpsi_l_dR[5];
  Float_t jpsi_l3d[5],jpsi_sigmal3d[5];

};

void createCharmEventTree(TTree *t, CharmEvent_t &ev);
void attachToCharmEventTree(TTree *t, CharmEvent_t &ev);

#endif
