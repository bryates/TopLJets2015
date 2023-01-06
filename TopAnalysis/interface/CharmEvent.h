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
    ht=0; met=0; st=0;
    nlj=0; event=0; nvtx=0;
  }

  Bool_t isData;
  Int_t run,event,epoch[50];
  Float_t norm,xsec,lumi,lum,puwgt[50],topptwgt;
  Float_t sfs[50],sfsu[50],pitrk[50],piptsf[50],kptsf[50]; //need to split into separate parts
  Float_t ht,met,st;
  Int_t nlj;


  //gen level event
  Int_t ttbar_nw, ngfsr;
  Float_t ttbar_w[500], gfsr[500];

  //reco level event
  Int_t nvtx;
  Int_t nl,nleptons;
  Int_t l_id[50],l_charge[50],l_pid[50],l_g[500];
  Float_t l_pt[50],l_eta[50],l_phi[50], l_mass[50];
  Float_t l_chi2norm[50], l_dxy[50], l_dxyE[50], l_dz[50], l_dzE[50];

  Int_t nj;
  Float_t j_pt[500],j_pt_pf[500],j_pt_charged[500],j_eta[500],j_phi[500],j_mass[500],j_ntk[50];
  Float_t j_p[500],j_p_pf[500],j_p_charged[500];
  Float_t j_pz[500],j_pz_pf[500],j_pz_charged[500];
  Float_t j_csv[500];
  Int_t j_hadflav[500], j_flav[500], j_pid[500];
  //Float_t gj_pt[500],gj_eta[500],gj_phi[500];

  //PF candidates
  Int_t npf,nmu,pf_j[5000];
  Int_t pf_id[5000],pf_c[5000];
  Float_t pf_pt[5000],pf_eta[5000],pf_phi[5000],pf_m[5000],pf_dxy[5000],pf_dxyE[5000],pf_dz[5000],pf_dzE[5000];

  //JPsi candidates
  Int_t njpsi;
  Float_t jpsi_mass[50],jpsi_pt[50],jpsi_eta[50],jpsi_phi[50],jpsi_p[50],jpsi_pz[50];
  Float_t jpsi_j[50],jpsi_ptrel[50],jpsi_l[50];//,jpsi_j_dR[50];
  Float_t jpsi_mu1_pt[50],jpsi_mu1_eta[50],jpsi_mu1_phi[50];
  Float_t jpsi_mu2_pt[50],jpsi_mu2_eta[50],jpsi_mu2_phi[50];
  Float_t jpsi_k_pt[50],jpsi_k_eta[50],jpsi_k_phi[50],jpsi_k_dxy[50],jpsi_k_dz[50],jpsi_chi2[50];
  Float_t jpsi_l_mass[50],jpsi_l_dR[50];
  Float_t jpsi_l3d[50],jpsi_sigmal3d[50];
  Float_t jpsi_lx[5000],jpsi_ly[5000],jpsi_lz[5000];
  Float_t jpsi_sigmax[5000],jpsi_sigmay[5000],jpsi_sigmaz[5000];

  //D meson candidates
  Int_t nmeson, meson_id[50],d0_j[50],d0_l[50],d0_ndau[50];
  Float_t d0_mass[50],d0_pt[50],d0_eta[50],d0_phi[50],d0_p[50],d0_pz[50];
  Float_t ds_mass[50],ds_pt[50],ds_eta[50],ds_phi[50],ds_p[50],ds_pz[50];
  Float_t d0_ptrel[50];//,d0_j_dR[50];
  Float_t d0_pi_pt[50],d0_pi_eta[50],d0_pi_phi[50],d0_pi_dxy[50],d0_pi_dxyE[50],d0_pi_dz[50];
  Float_t d0_k_pt[50],d0_k_eta[50],d0_k_phi[50],d0_k_dxy[50],d0_k_dxyE[50],d0_k_dz[50],d0_chi2[50],d0_opang[50];
  Float_t d0_mu_pt[50],d0_mu_eta[50],d0_mu_phi[50],d0_mu_tag_mu_pt[50];
  Float_t ds_pi2_pt[50],ds_pi2_eta[50],ds_pi2_phi[50];
  Float_t d0_l_mass[50],d0_l_dR[50];
  Float_t d0_l3d[50],d0_sigmal3d[50];
  Float_t d0_lx[5000],d0_ly[5000],d0_lz[5000];
  Float_t d0_sigmax[5000],d0_sigmay[5000],d0_sigmaz[5000];
  Float_t d0_pi_mother[50],d0_k_mother[50];

  //Fragmentation
  Float_t peterson[50], up[50], down[50], central[50];

};

void createCharmEventTree(TTree *t, CharmEvent_t &ev);
void attachToCharmEventTree(TTree *t, CharmEvent_t &ev);

#endif
