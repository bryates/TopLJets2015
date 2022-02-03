#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

//
void createMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  //event header
  t->Branch("isData",     &ev.isData,     "isData/O");
  t->Branch("run",       &ev.run,       "run/I");
  t->Branch("event",     &ev.event,     "event/I");
  t->Branch("lumi",      &ev.lumi,      "lumi/I");

  //generator level event
  t->Branch("pu",      &ev.pu,      "pu/I");
  t->Branch("putrue",      &ev.putrue,      "putrue/I");
  t->Branch("ttbar_nw",        &ev.ttbar_nw,        "ttbar_nw/I");
  t->Branch("ttbar_allmepartons",        &ev.ttbar_allmepartons,        "ttbar_allmepartons/I");
  t->Branch("ttbar_matchmepartons",        &ev.ttbar_matchmepartons,        "ttbar_matchmepartons/I");
  t->Branch("ttbar_w",        ev.ttbar_w,        "ttbar_w[ttbar_nw]/F");
  t->Branch("ngpsw",    &ev.ngpsw,   "ngpsw/I");
  t->Branch("gpsw",      ev.gpsw,    "gpsw[ngpsw]/F");

  //gen event (jets and dressed leptons)
  t->Branch("ngjets",       &ev.ngjets,       "ngjets/I");
  t->Branch("ng",           &ev.ng,           "ng/I");
  t->Branch("ngbjets",       &ev.ngbjets,       "ngbjets/I");
  t->Branch("ngleptons",       &ev.ngleptons,       "ngleptons/I");
  t->Branch("ngj",      &ev.ngj,      "ngj/I");
  t->Branch("g_id",      ev.g_id,     "g_id[ng]/I");
  t->Branch("g_pt",      ev.g_pt,     "g_pt[ng]/F");
  t->Branch("g_eta",     ev.g_eta,    "g_eta[ng]/F");
  t->Branch("g_phi",     ev.g_phi,    "g_phi[ng]/F");
  t->Branch("g_m",       ev.g_m,      "g_m[ng]/F");
  t->Branch("g_B",       ev.g_B,      "g_B[ng]/I");

  //gen level J/Psi
  t->Branch("ngjpsi",       &ev.ngjpsi, "ngjpsi/I");
  t->Branch("ngmeson",      &ev.ngmeson, "ngmeson/I");
  t->Branch("ngmeson_daug", &ev.ngmeson_daug, "ngmeson_daug/I");
  t->Branch("gmeson_mother_id",  ev.gmeson_mother_id, "gmeson_mother_id[ngmeson]/I");
  t->Branch("gmeson_daug_id",  ev.gmeson_daug_id, "gmeson_daug_id[ngmeson_daug]/I");
  t->Branch("gmeson_daug_pt",  ev.gmeson_daug_pt, "gmeson_daug_pt[ngmeson_daug]/F");
  t->Branch("gmeson_daug_eta", ev.gmeson_daug_eta, "gmeson_daug_eta[ngmeson_daug]/F");
  t->Branch("gmeson_daug_phi", ev.gmeson_daug_phi, "gmeson_daug_phi[ngmeson_daug]/F");
  t->Branch("ngmeson_daug_daug", &ev.ngmeson_daug_daug, "ngmeson_daug_daug/I");
  t->Branch("gmeson_mother_id",  ev.gmeson_mother_id, "gmeson_mother_id[ngmeson]/I");
  t->Branch("gmeson_daug_daug_id",  ev.gmeson_daug_daug_id, "gmeson_daug_daug_id[ngmeson_daug_daug]/I");
  t->Branch("gmeson_daug_daug_pt",  ev.gmeson_daug_daug_pt, "gmeson_daug_daug_pt[ngmeson_daug_daug]/F");
  t->Branch("gmeson_daug_daug_eta", ev.gmeson_daug_daug_eta, "gmeson_daug_daug_eta[ngmeson_daug_daug]/F");
  t->Branch("gmeson_daug_daug_phi", ev.gmeson_daug_daug_phi, "gmeson_daug_daug_phi[ngmeson_daug_daug]/F");
  t->Branch("gmeson_id",     ev.gmeson_id, "gmeson_id[ngmeson]/I");
  t->Branch("gmeson_pt",     ev.gmeson_pt, "gmeson_pt[ngmeson]/F");
  t->Branch("gmeson_eta",    ev.gmeson_eta, "gmeson_eta[ngmeson]/F");
  t->Branch("gmeson_phi",    ev.gmeson_phi, "gmeson_phi[ngmeson]/F");
  t->Branch("gmeson_m",      ev.gmeson_m, "gmeson_m[ngmeson]/F");
  t->Branch("gmeson_daug_dR",  ev.gmeson_daug_dR, "gmeson_daug_dR[ngmeson]/F");
  t->Branch("gmeson_index",  ev.gmeson_index, "gmeson_index[ngmeson]/F");
  t->Branch("gmeson_daug_meson_index",  ev.gmeson_daug_meson_index, "gmeson_daug_meson_index[ngmeson_daug]/F");
  t->Branch("gmeson_daug_daug_meson_index",  ev.gmeson_daug_daug_meson_index, "gmeson_daug_daug_meson_index[ngmeson_daug_daug]/F");

  //top (lastCopy and pseudo-top)
  t->Branch("ngtop",     &ev.ngtop,      "ngtop/I");
  t->Branch("gtop_id",    ev.gtop_id,    "gtop_id[ngtop]/I");
  t->Branch("gtop_pt",    ev.gtop_pt,    "gtop_pt[ngtop]/F");
  t->Branch("gtop_eta",   ev.gtop_eta,   "gtop_eta[ngtop]/F");
  t->Branch("gtop_phi",   ev.gtop_phi,   "gtop_phi[ngtop]/F");
  t->Branch("gtop_m",     ev.gtop_m,     "gtop_m[ngtop]/F");

  //final state
  t->Branch("ngpf",       &ev.ngpf,       "ngpf/I");
  t->Branch("gpf_id",      ev.gpf_id,     "gpf_id[ngpf]/I");
  t->Branch("gpf_c",       ev.gpf_c,      "gpf_c[ngpf]/I");
  t->Branch("gpf_g",       ev.gpf_g,      "gpf_g[ngpf]/I");
  t->Branch("gpf_mother",  ev.gpf_mother, "gpf_mother[ngpf]/I");
  t->Branch("gpf_pt",      ev.gpf_pt,     "gpf_pt[ngpf]/F");
  t->Branch("gpf_eta",     ev.gpf_eta,    "gpf_eta[ngpf]/F");
  t->Branch("gpf_phi",     ev.gpf_phi,    "gpf_phi[ngpf]/F");
  t->Branch("gpf_m",       ev.gpf_m,      "gpf_m[ngpf]/F");

  //reco level event
  t->Branch("nvtx",      &ev.nvtx,      "nvtx/I");
  t->Branch("muTrigger",        &ev.muTrigger,        "muTrigger/I");
  t->Branch("elTrigger",        &ev.elTrigger,        "elTrigger/I");

  t->Branch("nleptons", &ev.nleptons, "nleptons/I");
  t->Branch("nl", &ev.nl, "nl/I");
  t->Branch("isPromptFinalState",                         ev.isPromptFinalState,        "isPromptFinalState[nl]/O");
  t->Branch("isDirectPromptTauDecayProductFinalState",    ev.isDirectPromptTauDecayProductFinalState,        "isDirectPromptTauDecayProductFinalState[nl]/O");
  t->Branch("l_id",       ev.l_id,      "l_id[nl]/I");
  t->Branch("l_pid",      ev.l_pid,     "l_pid[nl]/I");
  t->Branch("l_g",        ev.l_g,       "l_g[nl]/I");
  t->Branch("l_charge",   ev.l_charge,  "l_charge[nl]/I");
  t->Branch("l_pt",       ev.l_pt,      "l_pt[nl]/F");
  t->Branch("l_eta",      ev.l_eta,     "l_eta[nl]/F");
  t->Branch("l_phi",      ev.l_phi,     "l_phi[nl]/F");
  t->Branch("l_mass",     ev.l_mass,    "l_mass[nl]/F");
  t->Branch("l_chargedHadronIso", ev.l_chargedHadronIso, "l_chargedHadronIso[nl]/F");
  t->Branch("l_miniIso",          ev.l_miniIso,          "l_miniIso[nl]/F");
  t->Branch("l_relIso",           ev.l_relIso,           "l_relIso[nl]/F");
  t->Branch("l_ip3d",             ev.l_ip3d,             "l_ip3d[nl]/F");
  t->Branch("l_ip3dsig",          ev.l_ip3dsig,          "l_ip3dsig[nl]/F");
  t->Branch("l_chi2norm",         ev.l_chi2norm,         "l_chi2norm[nl]/F");
  t->Branch("l_dxy",              ev.l_dxy,              "dxy[nl]/F");
  t->Branch("l_dxyE",             ev.l_dxyE,             "dxyE[nl]/F");
  t->Branch("l_dz",               ev.l_dz,               "dz[nl]/F");
  t->Branch("l_global",           ev.l_global,           "global[nl]/O");
  t->Branch("l_pf",               ev.l_pf,               "pf[nl]/O");
  t->Branch("l_nValTrackerHits",  ev.l_nValTrackerHits,  "nValTrackerHits[nl]/F");
  t->Branch("l_globalTrackNumberOfValidHits",                 ev.l_globalTrackNumberOfValidHits,                 "globalTrackNumberOfValidHits[nl]/F");
  t->Branch("l_nValPixelHits",    ev.l_nValPixelHits,    "nValPixelHits[nl]/F");
  t->Branch("l_pixelLayerWithMeasurement",    ev.l_pixelLayerWithMeasurement,    "pixelLayerWithMeasurement[nl]/F");
  t->Branch("l_nMatchedStations",    ev.l_nMatchedStations,    "nMatchedStations[nl]/F");
  t->Branch("l_trackerLayersWithMeasurement",    ev.l_trackerLayersWithMeasurement,    "trackerLayersWithMeasurement[nl]/F");
  t->Branch("l_validFraction",               ev.l_validFraction,               "pf[nl]/F");
  t->Branch("l_chi2LocalPosition",               ev.l_chi2LocalPosition,               "pf[nl]/F");
  t->Branch("l_trkKink",               ev.l_trkKink,               "pf[nl]/F");


  //jet info
  t->Branch("nj",        &ev.nj,        "nj/I");
  t->Branch("j_g",        ev.j_g,   "j_g[nj]/I");
  t->Branch("j_area",     ev.j_area,      "j_area[nj]/F");
  t->Branch("j_rawsf",    ev.j_rawsf,      "j_rawsf[nj]/F");
  t->Branch("j_pt",       ev.j_pt,      "j_pt[nj]/F");
  t->Branch("j_pt_pf",    ev.j_pt_pf,   "j_pt_pf[nj]/F");
  t->Branch("j_pt_charged", ev.j_pt_charged, "j_pt_charged[nj]/F");
  t->Branch("j_p",       ev.j_p,      "j_p[nj]/F");
  t->Branch("j_p_pf",    ev.j_p_pf,   "j_p_pf[nj]/F");
  t->Branch("j_p_charged", ev.j_p_charged, "j_p_charged[nj]/F");
  t->Branch("j_pz",       ev.j_pz,      "j_pz[nj]/F");
  t->Branch("j_pz_pf",    ev.j_pz_pf,   "j_pz_pf[nj]/F");
  t->Branch("j_pz_charged", ev.j_pz_charged, "j_pz_charged[nj]/F");
  t->Branch("j_eta",      ev.j_eta,     "j_eta[nj]/F");
  t->Branch("j_phi",      ev.j_phi,     "j_phi[nj]/F");
  t->Branch("j_mass",     ev.j_mass,     "j_mass[nj]/F");
  t->Branch("j_csv",      ev.j_csv,     "j_csv[nj]/F");
  t->Branch("j_csvl",     ev.j_cvsl,     "j_cvsl[nj]/F");
  t->Branch("j_cvsb",     ev.j_cvsb,     "j_cvsb[nj]/F");
  t->Branch("j_vtxpx",    ev.j_vtxpx, "j_vtxpx[nj]/F");
  t->Branch("j_vtxpy",    ev.j_vtxpy, "j_vtxpy[nj]/F");
  t->Branch("j_vtxpz",    ev.j_vtxpz, "j_vtxpz[nj]/F");
  t->Branch("j_vtxmass",  ev.j_vtxmass, "j_vtxmass[nj]/F");
  t->Branch("j_vtxNtracks",  ev.j_vtxNtracks, "j_vtxNtracks[nj]/I");
  t->Branch("j_vtx3DVal",    ev.j_vtx3DVal, "j_vtx3DVal[nj]/F");
  t->Branch("j_vtx3DSig",    ev.j_vtx3DSig, "j_vtx3DSig[nj]/F");
  t->Branch("j_puid",        ev.j_puid,    "j_puid[nj]/F");  
  t->Branch("j_flav",        ev.j_flav,    "j_flav[nj]/I");
  t->Branch("j_hadflav",     ev.j_hadflav,    "j_hadflav[nj]/I");
  t->Branch("j_pid",         ev.j_pid,     "j_pid[nj]/I");

  //pf candidates (only charged if outside jets)
  t->Branch("npf",        &ev.npf,         "npf/I");
  t->Branch("pf_j",        ev.pf_j,        "pf_j[npf]/I");
  t->Branch("pf_jnpf",     ev.pf_jnpf,     "pf_jnpf[npf]/I");
  t->Branch("pf_jnhppf",   ev.pf_jnhppf,   "pf_jnhppf[npf]/I");
  t->Branch("pf_id",       ev.pf_id,       "pf_id[npf]/I");
  t->Branch("pf_mother",   ev.pf_mother,   "pf_mother[npf]/I");
  t->Branch("pf_fromPV",   ev.pf_fromPV,   "pf_fromPV[npf]/I");
  t->Branch("pf_c",        ev.pf_c,        "pf_c[npf]/I");
  t->Branch("pf_pt",       ev.pf_pt,       "pf_pt[npf]/F");
  t->Branch("pf_eta",      ev.pf_eta,      "pf_eta[npf]/F");
  t->Branch("pf_phi",      ev.pf_phi,      "pf_phi[npf]/F");
  t->Branch("pf_m",        ev.pf_m,        "pf_m[npf]/F");
  t->Branch("pf_puppiWgt", ev.pf_puppiWgt, "pf_puppiWgt[npf]/F");
  t->Branch("pf_dxy",      ev.pf_dxy,      "pf_dxy[npf]/F");
  t->Branch("pf_dxyE",     ev.pf_dxyE,     "pf_dxyE[npf]/F");
  t->Branch("pf_dz",       ev.pf_dz,       "pf_dz[npf]/F");
  t->Branch("pf_dzE",      ev.pf_dzE,      "pf_dzE[npf]/F");
  t->Branch("pf_highPurity",      ev.pf_highPurity,      "pf_highPurity[npf]/O");
  t->Branch("pf_quality",      ev.pf_quality,      "pf_quality[npf]/I");
  t->Branch("pf_muon",   ev.pf_muon,   "pf_muon[npf]/O");
  t->Branch("pf_standAloneMuon",   ev.pf_standAloneMuon,   "pf_standAloneMuon[npf]/O");
  t->Branch("pf_globalMuon",   ev.pf_globalMuon,   "pf_globalMuon[npf]/O");
  t->Branch("pf_trackerMuon",   ev.pf_trackerMuon,   "pf_trackerMuon[npf]/O");
  t->Branch("pf_chi2ndof",   ev.pf_chi2ndof,   "pf_chi2ndof[npf]/F");
  t->Branch("pf_vtxchi2ndof",   ev.pf_vtxchi2ndof,   "pf_vtxchi2ndof[npf]/F");
  t->Branch("pf_relIso",   ev.pf_relIso,   "pf_relIso[npf]/F");

  //MET
  t->Branch("nmet",      &ev.nmet,     "nmet/I");
  t->Branch("met_pt",     ev.met_pt,   "met_pt[nmet]/F");
  t->Branch("met_phi",    ev.met_phi,  "met_phi[nmet]/F");

  //Kalman Filter
  t->Branch("njpsi",      &ev.njpsi,        "njpsi/I");
  t->Branch("nmeson",     &ev.nmeson,       "nmeson/I");
  t->Branch("nkj",        &ev.nkj,          "nkj/I");
  t->Branch("nkpf",       &ev.nkpf,         "nkpf/I");
  t->Branch("pf_idx",      ev.pf_idx,       "pf_idx[nkpf]/I");
  t->Branch("k_j",         ev.k_j,          "k_j[nj]/I");
  t->Branch("k_j_pt",      ev.k_j_pt,       "k_j_pt[nj]/F");
  t->Branch("k_j_eta",     ev.k_j_eta,      "k_j_eta[nj]/F");
  t->Branch("k_j_phi",     ev.k_j_phi,      "k_j_phi[nj]/F");
  t->Branch("k_j_mass",    ev.k_j_mass,     "k_j_mass[nj]/F");
  t->Branch("k_pf_ndau",   ev.k_pf_ndau,    "k_pf_ndau[nkpf]/I");
  t->Branch("k_pf_id",     ev.k_pf_id,      "k_pf_id[nkpf]/F");
  t->Branch("k_pf_pt",     ev.k_pf_pt,      "k_pf_pt[nkpf]/F");
  t->Branch("k_pf_ptE",    ev.k_pf_ptE,     "k_pf_ptE[nkpf]/F");
  t->Branch("k_pf_eta",    ev.k_pf_eta,     "k_pf_eta[nkpf]/F");
  t->Branch("k_pf_phi",    ev.k_pf_phi,     "k_pf_phi[nkpf]/F");
  t->Branch("k_pf_m",      ev.k_pf_m,       "k_pf_m[nkpf]/F");
  t->Branch("k_pf_dxy",    ev.k_pf_dxy,     "k_pf_dxy[nkpf]/F");
  t->Branch("k_pf_dxyE",   ev.k_pf_dxyE,    "k_pf_dxyE[nkpf]/F");
  t->Branch("k_pf_dz",     ev.k_pf_dz,      "k_pf_dz[nkpf]/F");
  t->Branch("k_pf_dzE",    ev.k_pf_dzE,     "k_pf_dzE[nkpf]/F");
  t->Branch("k_pf_tracker",      ev.k_pf_tracker,       "k_pf_tracker[nkpf]/O");
  t->Branch("k_pf_global",      ev.k_pf_global,       "k_pf_global[nkpf]/O");
  t->Branch("k_mass",      ev.k_mass,       "k_mass[nkpf]/F");
  t->Branch("k_chi2",      ev.k_chi2,       "k_chi2[nkpf]/F");
  t->Branch("k_ndof",      ev.k_ndof,       "k_ndof[nkpf]/I");
  t->Branch("k_vtxProb",   ev.k_vtxProb,    "k_vtxProb[nkpf]/F");
  t->Branch("k_id",        ev.k_id,         "k_id[nkpf]/I");
  t->Branch("k_dxy",       ev.k_dxy,        "k_dxy[nkpf]/F");
  t->Branch("k_dxyE",      ev.k_dxyE,       "k_dxyE[nkpf]/F");
  t->Branch("k_opang",     ev.k_opang,      "k_opang[nkpf]/F");
  t->Branch("k_l3d",       ev.k_l3d,        "k_l3d[nkpf]/F");
  t->Branch("k_lx",        ev.k_lx,         "k_lx[nkpf]/F");
  t->Branch("k_ly",        ev.k_ly,         "k_ly[nkpf]/F");
  t->Branch("k_lz",        ev.k_lz,         "k_lz[nkpf]/F");
  t->Branch("k_sigmal3d",  ev.k_sigmal3d,   "k_sigmal3d[nkpf]/F");
  t->Branch("k_sigmax",    ev.k_sigmax,     "k_sigmax[nkpf]/F");
  t->Branch("k_sigmay",    ev.k_sigmay,     "k_sigmay[nkpf]/F");
  t->Branch("k_sigmaz",    ev.k_sigmaz,     "k_sigmaz[nkpf]/F");

  //Fragmentation
  t->Branch("peterson",    ev.peterson,     "peterson[ng]/F");
  t->Branch("up",          ev.up,           "up[ng]/F");
  t->Branch("central",     ev.central,      "central[ng]/F");
  t->Branch("down",        ev.down,         "down[ng]/F");
  t->Branch("xb",          ev.xb,           "xb[ngjets]/F");
}

//
void attachToMiniEventTree(TTree *t,MiniEvent_t &ev,bool full)
{
  //event header
  t->SetBranchAddress("isData",     &ev.isData);
  t->SetBranchAddress("run",       &ev.run);
  t->SetBranchAddress("event",     &ev.event);
  t->SetBranchAddress("lumi",      &ev.lumi);

  //generator level event
  t->SetBranchAddress("pu",      &ev.pu);
  t->SetBranchAddress("putrue",      &ev.putrue);
  t->SetBranchAddress("ttbar_nw",        &ev.ttbar_nw);
  t->SetBranchAddress("ttbar_allmepartons",        &ev.ttbar_allmepartons);
  t->SetBranchAddress("ttbar_matchmepartons",        &ev.ttbar_matchmepartons);
  t->SetBranchAddress("ttbar_w",        ev.ttbar_w);
  if(t->GetListOfBranches()->FindObject("ngpsw")) {
  t->SetBranchAddress("ngpsw",      &ev.ngpsw);
  t->SetBranchAddress("gpsw",        ev.gpsw);
  }

  //gen event (jets and dressed leptons)
  t->SetBranchAddress("ngjets",       &ev.ngjets);
  t->SetBranchAddress("ng",           &ev.ng);
  t->SetBranchAddress("ngbjets",       &ev.ngbjets);
  t->SetBranchAddress("ngleptons",       &ev.ngleptons);
  t->SetBranchAddress("ng",       &ev.ng);
  t->SetBranchAddress("g_id",      ev.g_id);
  t->SetBranchAddress("g_pt",      ev.g_pt);
  t->SetBranchAddress("g_eta",     ev.g_eta);
  t->SetBranchAddress("g_phi",     ev.g_phi);
  t->SetBranchAddress("g_m",       ev.g_m);
  t->SetBranchAddress("g_B",       ev.g_B);

  //top (lastCopy and pseudo-top)
  t->SetBranchAddress("ngtop",     &ev.ngtop);
  t->SetBranchAddress("gtop_id",    ev.gtop_id);
  t->SetBranchAddress("gtop_pt",    ev.gtop_pt);
  t->SetBranchAddress("gtop_eta",   ev.gtop_eta);
  t->SetBranchAddress("gtop_phi",   ev.gtop_phi);
  t->SetBranchAddress("gtop_m",     ev.gtop_m);

  //final state
  if(full)
    {
      t->SetBranchAddress("ngpf",       &ev.ngpf);
      t->SetBranchAddress("gpf_id",      ev.gpf_id);
      t->SetBranchAddress("gpf_c",       ev.gpf_c);
      t->SetBranchAddress("gpf_g",       ev.gpf_g); 
      t->SetBranchAddress("gpf_mother",  ev.gpf_mother); 
      t->SetBranchAddress("gpf_pt",      ev.gpf_pt);
      t->SetBranchAddress("gpf_eta",     ev.gpf_eta);
      t->SetBranchAddress("gpf_phi",     ev.gpf_phi);
      t->SetBranchAddress("gpf_m",       ev.gpf_m);
    }

  //gen level J/Psi
  t->SetBranchAddress("ngjpsi",       &ev.ngjpsi);
  t->SetBranchAddress("ngmeson",      &ev.ngmeson);
  t->SetBranchAddress("ngmeson_daug", &ev.ngmeson_daug);
  t->SetBranchAddress("gmeson_mother_id",  ev.gmeson_mother_id);
  t->SetBranchAddress("gmeson_daug_id",  ev.gmeson_daug_id);
  t->SetBranchAddress("gmeson_daug_pt",  ev.gmeson_daug_pt);
  t->SetBranchAddress("gmeson_daug_eta", ev.gmeson_daug_eta);
  t->SetBranchAddress("gmeson_daug_phi", ev.gmeson_daug_phi);
  t->SetBranchAddress("gmeson_id",     ev.gmeson_id);
  t->SetBranchAddress("gmeson_pt",     ev.gmeson_pt);
  t->SetBranchAddress("gmeson_eta",    ev.gmeson_eta);
  t->SetBranchAddress("gmeson_phi",    ev.gmeson_phi);
  t->SetBranchAddress("gmeson_m",      ev.gmeson_m);
  t->SetBranchAddress("gmeson_daug_dR",  ev.gmeson_daug_dR);
  t->SetBranchAddress("gmeson_index",  ev.gmeson_index);
  t->SetBranchAddress("gmeson_daug_meson_index",  ev.gmeson_daug_meson_index);

  //reco level event
  t->SetBranchAddress("nvtx",      &ev.nvtx);
  t->SetBranchAddress("muTrigger",        &ev.muTrigger);
  t->SetBranchAddress("elTrigger",        &ev.elTrigger);

  t->SetBranchAddress("nleptons", &ev.nleptons);
  t->SetBranchAddress("nl", &ev.nl);
  t->SetBranchAddress("isPromptFinalState",                         ev.isPromptFinalState);
  t->SetBranchAddress("isDirectPromptTauDecayProductFinalState",    ev.isDirectPromptTauDecayProductFinalState);
  t->SetBranchAddress("l_id",       ev.l_id);
  t->SetBranchAddress("l_pid",      ev.l_pid);
  t->SetBranchAddress("l_g",        ev.l_g);
  t->SetBranchAddress("l_charge",   ev.l_charge);
  t->SetBranchAddress("l_pt",       ev.l_pt);
  t->SetBranchAddress("l_eta",      ev.l_eta);
  t->SetBranchAddress("l_phi",      ev.l_phi);
  t->SetBranchAddress("l_mass",     ev.l_mass);
  t->SetBranchAddress("l_chargedHadronIso", ev.l_chargedHadronIso);
  t->SetBranchAddress("l_miniIso",          ev.l_miniIso);
  t->SetBranchAddress("l_relIso",           ev.l_relIso);
  t->SetBranchAddress("l_ip3d",             ev.l_ip3d);
  t->SetBranchAddress("l_ip3dsig",          ev.l_ip3dsig);
  t->SetBranchAddress("l_chi2norm",         ev.l_chi2norm);
  t->SetBranchAddress("l_dxy",              ev.l_dxy);
  t->SetBranchAddress("l_dz",               ev.l_dz);
  t->SetBranchAddress("l_global",           ev.l_global);
  t->SetBranchAddress("l_pf",               ev.l_pf);
  t->SetBranchAddress("l_nValTrackerHits",  ev.l_nValTrackerHits);
  t->SetBranchAddress("l_globalTrackNumberOfValidHits",               ev.l_globalTrackNumberOfValidHits);
  t->SetBranchAddress("l_nValPixelHits",    ev.l_nValPixelHits);
  t->SetBranchAddress("l_pixelLayerWithMeasurement",    ev.l_pixelLayerWithMeasurement);
  t->SetBranchAddress("l_nMatchedStations",    ev.l_nMatchedStations);
  t->SetBranchAddress("l_trackerLayersWithMeasurement",    ev.l_trackerLayersWithMeasurement);
  t->SetBranchAddress("l_validFraction",               ev.l_validFraction);
  t->SetBranchAddress("l_chi2LocalPosition",               ev.l_chi2LocalPosition);
  t->SetBranchAddress("l_trkKink",               ev.l_trkKink);

  //jet info
  t->SetBranchAddress("nj",        &ev.nj);
  t->SetBranchAddress("j_g",        ev.j_g);
  t->SetBranchAddress("j_area",     ev.j_area);
  t->SetBranchAddress("j_rawsf",    ev.j_rawsf);
  t->SetBranchAddress("j_pt",       ev.j_pt);
  t->SetBranchAddress("j_pt_pf", ev.j_pt_pf);
  t->SetBranchAddress("j_pt_charged", ev.j_pt_charged);
  t->SetBranchAddress("j_p",       ev.j_p);
  t->SetBranchAddress("j_p_pf", ev.j_p_pf);
  t->SetBranchAddress("j_p_charged", ev.j_p_charged);
  t->SetBranchAddress("j_pz",       ev.j_pz);
  t->SetBranchAddress("j_pz_pf", ev.j_pz_pf);
  t->SetBranchAddress("j_pz_charged", ev.j_pz_charged);
  t->SetBranchAddress("j_eta",      ev.j_eta);
  t->SetBranchAddress("j_phi",      ev.j_phi);
  t->SetBranchAddress("j_mass",     ev.j_mass);
  t->SetBranchAddress("j_csv",      ev.j_csv);
  t->SetBranchAddress("j_csvl",     ev.j_cvsl);
  t->SetBranchAddress("j_cvsb",     ev.j_cvsb);
  t->SetBranchAddress("j_vtxpx",    ev.j_vtxpx);
  t->SetBranchAddress("j_vtxpy",    ev.j_vtxpy);
  t->SetBranchAddress("j_vtxpz",    ev.j_vtxpz);
  t->SetBranchAddress("j_vtxmass",  ev.j_vtxmass);
  t->SetBranchAddress("j_vtxNtracks",  ev.j_vtxNtracks);
  t->SetBranchAddress("j_vtx3DVal",    ev.j_vtx3DVal);
  t->SetBranchAddress("j_vtx3DSig",    ev.j_vtx3DSig);
  t->SetBranchAddress("j_puid",        ev.j_puid);
  t->SetBranchAddress("j_flav",        ev.j_flav);
  t->SetBranchAddress("j_hadflav",     ev.j_hadflav);
  t->SetBranchAddress("j_pid",         ev.j_pid);

  //pf candidates (only charged if outside jets)
  if(full)
    {
      t->SetBranchAddress("npf",        &ev.npf);
      t->SetBranchAddress("pf_j",        ev.pf_j);
      t->SetBranchAddress("pf_jnpf",     ev.pf_jnpf);
      t->SetBranchAddress("pf_jnhppf",   ev.pf_jnhppf);
      t->SetBranchAddress("pf_id",       ev.pf_id);
      t->SetBranchAddress("pf_mother",   ev.pf_mother); 
      t->SetBranchAddress("pf_fromPV",   ev.pf_fromPV);
      t->SetBranchAddress("pf_c",        ev.pf_c);
      t->SetBranchAddress("pf_pt",       ev.pf_pt);
      t->SetBranchAddress("pf_eta",      ev.pf_eta);
      t->SetBranchAddress("pf_phi",      ev.pf_phi);
      t->SetBranchAddress("pf_m",        ev.pf_m);
      t->SetBranchAddress("pf_puppiWgt", ev.pf_puppiWgt);

      t->SetBranchAddress("pf_dxy",      ev.pf_dxy);
      t->SetBranchAddress("pf_dxyE",     ev.pf_dxyE);
      t->SetBranchAddress("pf_dz",       ev.pf_dz);
      t->SetBranchAddress("pf_dzE",      ev.pf_dzE);
      t->SetBranchAddress("pf_highPurity",      ev.pf_highPurity);
      t->SetBranchAddress("pf_quality",      ev.pf_quality);
      t->SetBranchAddress("pf_muon",       ev.pf_muon);
      t->SetBranchAddress("pf_standAloneMuon",   ev.pf_standAloneMuon);
      t->SetBranchAddress("pf_globalMuon",       ev.pf_globalMuon);
      t->SetBranchAddress("pf_trackerMuon",      ev.pf_trackerMuon);
      t->SetBranchAddress("pf_chi2ndof",         ev.pf_chi2ndof);
      t->SetBranchAddress("pf_vtxchi2ndof",      ev.pf_vtxchi2ndof);
      t->SetBranchAddress("pf_relIso",      ev.pf_relIso);
    }

  //MET
  t->SetBranchAddress("nmet",      &ev.nmet);
  t->SetBranchAddress("met_pt",    ev.met_pt);
  t->SetBranchAddress("met_phi",   ev.met_phi);

  //Kalman Filter
  t->SetBranchAddress("njpsi",      &ev.njpsi);
  t->SetBranchAddress("nmeson",     &ev.nmeson);
  t->SetBranchAddress("nkpf",       &ev.nkpf);
  t->SetBranchAddress("pf_idx",      ev.pf_idx);
  t->SetBranchAddress("k_j",         ev.k_j);
  t->SetBranchAddress("k_j_pt",      ev.k_j_pt);
  t->SetBranchAddress("k_j_eta",     ev.k_j_eta);
  t->SetBranchAddress("k_j_phi",     ev.k_j_phi);
  t->SetBranchAddress("k_j_mass",    ev.k_j_mass);
  t->SetBranchAddress("k_pf_ndau",   ev.k_pf_ndau);
  t->SetBranchAddress("k_pf_id",     ev.k_pf_id);
  t->SetBranchAddress("k_pf_pt",     ev.k_pf_pt);
  t->SetBranchAddress("k_pf_ptE",    ev.k_pf_ptE);
  t->SetBranchAddress("k_pf_eta",    ev.k_pf_eta);
  t->SetBranchAddress("k_pf_phi",    ev.k_pf_phi);
  t->SetBranchAddress("k_pf_m",      ev.k_pf_m);
  t->SetBranchAddress("k_pf_dxy",    ev.k_pf_dxy);
  t->SetBranchAddress("k_pf_dxyE",   ev.k_pf_dxyE);
  t->SetBranchAddress("k_pf_dz",     ev.k_pf_dz);
  t->SetBranchAddress("k_pf_dzE",    ev.k_pf_dzE);
  t->SetBranchAddress("k_pf_tracker",      ev.k_pf_tracker);
  t->SetBranchAddress("k_pf_global",      ev.k_pf_global);
  t->SetBranchAddress("k_mass",      ev.k_mass);
  t->SetBranchAddress("k_chi2",      ev.k_chi2);
  t->SetBranchAddress("k_ndof",      ev.k_ndof);
  t->SetBranchAddress("k_vtxProb",   ev.k_vtxProb);
  t->SetBranchAddress("k_id",        ev.k_id);
  t->SetBranchAddress("k_dxy",       ev.k_dxy);
  t->SetBranchAddress("k_dxyE",      ev.k_dxyE);
  t->SetBranchAddress("k_opang",  ev.k_opang);
  t->SetBranchAddress("k_l3d",       ev.k_l3d);
  t->SetBranchAddress("k_lx",        ev.k_lx);
  t->SetBranchAddress("k_ly",        ev.k_ly);
  t->SetBranchAddress("k_lz",        ev.k_lz);
  t->SetBranchAddress("k_sigmal3d",  ev.k_sigmal3d);
  t->SetBranchAddress("k_sigmax",  ev.k_sigmax);
  t->SetBranchAddress("k_sigmay",  ev.k_sigmay);
  t->SetBranchAddress("k_sigmaz",  ev.k_sigmaz);

  //Fragmentation
  t->SetBranchAddress("peterson",    ev.peterson);
  t->SetBranchAddress("up",          ev.up);
  t->SetBranchAddress("central",     ev.central);
  t->SetBranchAddress("down",        ev.down);
  t->SetBranchAddress("xb",          ev.xb);
}
