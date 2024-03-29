#ifndef _minievent_h_
#define _minievent_h_

#include "TTree.h"

struct MiniEvent_t
{
  MiniEvent_t()
  {
    ttbar_nw=0;
    ng=0; ngtop=0; ngpf=0; ngjpsi=0; ngmeson=0; ngmeson_daug=0; ngmeson_daug_daug=0;
    nl=0; nj=0; nmet=0; npf=0;
    ngpsw=0;
  }

  Bool_t isData;
  Int_t run,event,lumi;


  //gen level event
  Int_t pu,putrue;
  Int_t ttbar_nw, ttbar_allmepartons, ttbar_matchmepartons;
  Float_t ttbar_w[500];
  Int_t ng,ngj,ngjets,ngbjets,ngleptons,ngtop,ngpf,ngjpsi,ngmeson,ngmeson_daug,ngmeson_daug_daug;
  Float_t gtop_pt_wgt;
  Int_t g_id[500];
  Float_t g_pt[500],g_eta[500],g_phi[500],g_m[500]; 
  Int_t g_B[500];
  Int_t gtop_id[15];
  Float_t gtop_pt[15],gtop_eta[15],gtop_phi[15],gtop_m[15]; 
  Int_t gpf_id[5000],gpf_c[5000],gpf_g[5000],gpf_mother[5000];
  Float_t gpf_pt[5000],gpf_eta[5000],gpf_phi[5000],gpf_m[5000];
  Int_t gmeson_id[5000],gmeson_daug_id[5000],gmeson_mother_id[5000];
  Int_t gmeson_daug_daug_id[5000];
  Float_t gmeson_pt[5000],gmeson_eta[5000],gmeson_phi[5000],gmeson_m[5000],gmeson_daug_dR[5000],gmeson_index[5000];
  Float_t gmeson_daug_pt[5000],gmeson_daug_eta[5000],gmeson_daug_phi[5000],gmeson_daug_m[5000],gmeson_daug_meson_index[5000],gmeson_daug_dxy[5000],gmeson_daug_dxyE[5000],gmeson_daug_dz[5000],gmeson_daug_dzE[5000];
  Float_t gmeson_daug_daug_pt[5000],gmeson_daug_daug_eta[5000],gmeson_daug_daug_phi[5000],gmeson_daug_daug_m[5000],gmeson_daug_daug_meson_index[5000],gmeson_daug_daug_dxy[5000],gmeson_daug_daug_dxyE[5000],gmeson_daug_daug_dz[5000],gmeson_daug_daug_dzE[5000];
  Int_t ngpsw;
  Float_t gpsw[500];

  //reco level event
  Int_t nvtx;
  Int_t muTrigger,elTrigger;
  Float_t rho;
  Int_t nl,nleptons;
  Bool_t isPromptFinalState[50], isDirectPromptTauDecayProductFinalState[50];
  Int_t l_id[50],l_charge[50],l_pid[50],l_g[200];
  Float_t l_pt[50],l_eta[50],l_phi[50], l_mass[50], l_miniIso[50], l_chargedHadronIso[50], l_relIso[50], l_ip3d[50], l_ip3dsig[50];
  Float_t l_chi2norm[50], l_dxy[50], l_dxyE[50], l_dz[50];
  Bool_t l_global[50], l_pf[50];
  Float_t l_nValTrackerHits[50], l_globalTrackNumberOfValidHits[50], l_nValPixelHits[50], l_nMatchedStations[50], l_pixelLayerWithMeasurement[50], l_trackerLayersWithMeasurement[50], l_validFraction[50], l_chi2LocalPosition[50], l_trkKink[50];

  Int_t nj;
  Float_t j_pt[200],j_pt_pf[200],j_pt_charged[200],j_eta[200],j_phi[200],j_mass[200],j_area[200],j_rawsf[200];
  Float_t j_p[200],j_p_pf[200],j_p_charged[200],j_pz[200],j_pz_pf[200],j_pz_charged[200];
  Float_t j_csv[200],j_cvsl[200],j_cvsb[200],j_vtxmass[200],j_vtx3DVal[200],j_vtx3DSig[200],j_puid[200],j_vtxpx[200],j_vtxpy[200],j_vtxpz[200];
  Int_t j_vtxNtracks[200],j_flav[200],j_pid[200],j_hadflav[200],j_g[200];

  //met
  Int_t nmet;
  Float_t met_pt[10],met_phi[10];

  //PF candidates
  Int_t npf,nmu,pf_j[5000],pf_jnpf[5000],pf_jnhppf[5000];
  Int_t pf_id[5000],pf_c[5000],pf_fromPV[5000],pf_quality[5000],pf_mother[5000];
  Float_t pf_pt[5000],pf_eta[5000],pf_phi[5000],pf_m[5000],pf_puppiWgt[5000],pf_dxy[5000],pf_dxyE[5000],pf_dz[5000],pf_dzE[5000], pf_chi2ndof[5000], pf_vtxchi2ndof[5000],pf_relIso[5000];
  Bool_t pf_highPurity[5000],pf_muon[5000],pf_standAloneMuon[5000],pf_globalMuon[5000],pf_trackerMuon[5000];

  //Kalman Filter
  Int_t njpsi,nmeson,nkj,nkpf;
  Int_t pf_idx[5000];
  Int_t k_j[5000],k_ndof[5000],k_id[5000],k_pf_ndau[5000];
  Float_t k_j_pt[5000],k_j_eta[5000],k_j_phi[5000],k_j_mass[5000];
  Float_t k_pf_id[5000],k_pf_pt[5000],k_pf_ptE[5000],k_pf_eta[5000],k_pf_phi[5000],k_pf_m[5000];
  Float_t k_pf_dxy[5000],k_pf_dxyE[5000],k_pf_dz[5000],k_pf_dzE[5000];
  Bool_t k_pf_tracker[5000],k_pf_global[5000];
  Float_t k_mass[5000],k_chi2[5000],k_vtxProb[5000];
  Float_t k_dxy[5000],k_dxyE[5000],k_opang[5000];
  Float_t k_l3d[5000],k_sigmal3d[5000];
  Float_t k_lx[5000],k_ly[5000],k_lz[5000];
  Float_t k_sigmax[5000],k_sigmay[5000],k_sigmaz[5000];

  //Fragmentation
  Float_t peterson[500], up[500], down[500], central[500];
  Float_t xb[500];
};

void createMiniEventTree(TTree *t,MiniEvent_t &ev);
void attachToMiniEventTree(TTree *t, MiniEvent_t &ev,bool full=false);

#endif
