#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

//
void createMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  //event header
  t->Branch("isData",     &ev.isData,     "isData/O");
  t->Branch("run",       &ev.run,       "run/I");
  t->Branch("event",     &ev.event,     "event/I");
  t->Branch("lumi",      &ev.lumi,      "lumi/I");

  //generator level weights
  t->Branch("ttbar_nw",        &ev.ttbar_nw,        "ttbar_nw/I");
  t->Branch("ttbar_allmepartons",        &ev.ttbar_allmepartons,        "ttbar_allmepartons/I");
  t->Branch("ttbar_matchmepartons",        &ev.ttbar_matchmepartons,        "ttbar_matchmepartons/I");
  t->Branch("ttbar_w",        ev.ttbar_w,        "ttbar_w[ttbar_nw]/F");
  t->Branch("ttbar_genId",    &ev.ttbar_genId,    "ttbar_genId/I");

  //generator level flag for events in the fiducial region (1l+1j)
  t->Branch("isFiducial",       &ev.isFiducial,       "isFiducial/O");

  //trigger information
  t->Branch("muTrigger",        &ev.muTrigger,        "muTrigger/I");
  t->Branch("elTrigger",        &ev.elTrigger,        "elTrigger/I");

  //pileup information
  t->Branch("nvtx",      &ev.nvtx,      "nvtx/I");
  t->Branch("pu",      &ev.pu,      "pu/I");
  t->Branch("putrue",      &ev.putrue,      "putrue/I");

  //lepton information
  t->Branch("nl", &ev.nl, "nl/I");
  t->Branch("isPromptFinalState",        ev.isPromptFinalState,        "isPromptFinalState[nl]/O");
  t->Branch("isDirectPromptTauDecayProductFinalState",       ev.isDirectPromptTauDecayProductFinalState,        "isDirectPromptTauDecayProductFinalState[nl]/O");
  t->Branch("l_id",      ev.l_id,          "l_id[nl]/I");
  t->Branch("l_pid",      ev.l_pid,          "l_pid[nl]/I");
  t->Branch("l_charge",  ev.l_charge,  "l_charge[nl]/I");
  t->Branch("l_pt",      ev.l_pt,      "l_pt[nl]/F");
  t->Branch("l_eta",     ev.l_eta,     "l_eta[nl]/F");
  t->Branch("l_phi",     ev.l_phi,     "l_phi[nl]/F");
  t->Branch("l_mass",    ev.l_mass,    "l_mass[nl]/F");
  t->Branch("l_chargedHadronIso",       ev.l_chargedHadronIso,     "l_miniIso[nl]/F");
  t->Branch("l_miniIso",        ev.l_miniIso,     "l_miniIso[nl]/F");
  t->Branch("l_relIso",               ev.l_relIso,     "l_relIso[nl]/F");
  t->Branch("l_ip3d",      ev.l_ip3d,      "l_ip3d[nl]/F");
  t->Branch("l_ip3dsig",      ev.l_ip3dsig,      "l_ip3dsig[nl]/F");

  //jet info
  t->Branch("nj",        &ev.nj,        "nj/I");
  t->Branch("ngenj",        &ev.ngenj,        "ngenj/I");
  t->Branch("j_area",       ev.j_area,      "j_area[nj]/F");
  t->Branch("j_pt",       ev.j_pt,      "j_pt[nj]/F");
  t->Branch("j_eta",      ev.j_eta,     "j_eta[nj]/F");
  t->Branch("j_phi",      ev.j_phi,     "j_phi[nj]/F");
  t->Branch("j_mass",      ev.j_mass,     "j_mass[nj]/F");
  t->Branch("genj_pt",       ev.genj_pt,      "genj_pt[nj]/F");
  t->Branch("genj_eta",      ev.genj_eta,     "genj_eta[nj]/F");
  t->Branch("genj_phi",      ev.genj_phi,     "genj_phi[nj]/F");
  t->Branch("genj_mass",       ev.genj_mass,      "genj_mass[nj]/F");
  t->Branch("j_csv",      ev.j_csv,     "j_csv[nj]/F");
  t->Branch("j_vtxpx",  ev.j_vtxpx, "j_vtxpx[nj]/F");
  t->Branch("j_vtxpy",  ev.j_vtxpy, "j_vtxpy[nj]/F");
  t->Branch("j_vtxpz",  ev.j_vtxpz, "j_vtxpz[nj]/F");
  t->Branch("j_vtxmass",  ev.j_vtxmass, "j_vtxmass[nj]/F");
  t->Branch("j_vtxNtracks",  ev.j_vtxNtracks, "j_vtxNtracks[nj]/I");
  t->Branch("j_vtx3DVal",  ev.j_vtx3DVal, "j_vtx3DVal[nj]/F");
  t->Branch("j_vtx3DSig",  ev.j_vtx3DSig, "j_vtx3DSig[nj]/F");
  t->Branch("j_puid",     ev.j_puid,    "j_puid[nj]/F");
  t->Branch("j_qg",     ev.j_qg,    "j_qg[nj]/F");
  t->Branch("j_flav",     ev.j_flav,    "j_flav[nj]/I");
  t->Branch("j_hadflav",     ev.j_hadflav,    "j_hadflav[nj]/I");
  t->Branch("j_pid",      ev.j_pid,     "j_pid[nj]/I");

  //pf candidates in jets
  t->Branch("npf",        &ev.npf,      "npf/I");
  t->Branch("pf_j",       ev.pf_j,     "pf_j[npf]/I");
  t->Branch("pf_id",      ev.pf_id,     "pf_id[npf]/I");
  t->Branch("pf_charge",  ev.pf_charge, "pf_charge[npf]/I");
  t->Branch("pf_px",      ev.pf_px,     "pf_px[npf]/F");
  t->Branch("pf_py",      ev.pf_py,     "pf_py[npf]/F");
  t->Branch("pf_pz",      ev.pf_pz,     "pf_pz[npf]/F");

  //gen particles in gen jets
  t->Branch("ngen",     &ev.ngen,     "ngen/I");
  t->Branch("g_j",       ev.g_j,      "g_j[ngen]/I");
  t->Branch("g_id",      ev.g_id,     "g_id[ngen]/I");
  t->Branch("g_charge",  ev.g_charge, "g_charge[ngen]/I");
  t->Branch("g_px",      ev.g_px,     "g_px[ngen]/F");
  t->Branch("g_py",      ev.g_py,     "g_py[ngen]/F");
  t->Branch("g_pz",      ev.g_pz,     "g_pz[ngen]/F");

  //gen particles from hard process
  t->Branch("ngenHardProc",     &ev.ngenHardProc,     "ngenHardProc/I");
  t->Branch("ghp_id",      ev.ghp_id,     "ghp_id[ngenHardProc]/I");
  t->Branch("ghp_pt",      ev.ghp_pt,     "ghp_pt[ngenHardProc]/F");
  t->Branch("ghp_eta",      ev.ghp_eta,     "ghp_eta[ngenHardProc]/F");
  t->Branch("ghp_phi",      ev.ghp_phi,     "ghp_phi[ngenHardProc]/F");
  t->Branch("ghp_m",      ev.ghp_m,     "ghp_m[ngenHardProc]/F");

  t->Branch("nmet",      &ev.nmet,     "nmet/I");
  t->Branch("met_pt",    ev.met_pt,    "met_pt[nmet]/F");
  t->Branch("met_phi",   ev.met_phi,   "met_phi[nmet]/F");
}

//
void attachToMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  t->SetBranchAddress("isData",     &ev.isData);

  t->SetBranchAddress("ttbar_nw",              &ev.ttbar_nw);
  t->SetBranchAddress("ttbar_w",                ev.ttbar_w);
  t->SetBranchAddress("ttbar_allmepartons",    &ev.ttbar_allmepartons);
  t->SetBranchAddress("ttbar_matchmepartons",  &ev.ttbar_matchmepartons);
  t->SetBranchAddress("ttbar_genId",           &ev.ttbar_genId);

  t->SetBranchAddress("run",       &ev.run);
  t->SetBranchAddress("event",     &ev.event);
  t->SetBranchAddress("lumi",      &ev.lumi);

  t->SetBranchAddress("isFiducial",       &ev.isFiducial);
  t->SetBranchAddress("muTrigger",        &ev.muTrigger);
  t->SetBranchAddress("elTrigger",        &ev.elTrigger);

  t->SetBranchAddress("nvtx",    &ev.nvtx);
  t->SetBranchAddress("pu",      &ev.pu);
  t->SetBranchAddress("putrue",  &ev.putrue);

  t->SetBranchAddress("nl",      &ev.nl);
  t->SetBranchAddress("isPromptFinalState",                      ev.isPromptFinalState);
  t->SetBranchAddress("isDirectPromptTauDecayProductFinalState", ev.isDirectPromptTauDecayProductFinalState);
  t->SetBranchAddress("l_id",      ev.l_id);
  t->SetBranchAddress("l_pid",     ev.l_pid);
  t->SetBranchAddress("l_charge",  ev.l_charge);
  t->SetBranchAddress("l_pt",      ev.l_pt);
  t->SetBranchAddress("l_eta",     ev.l_eta);
  t->SetBranchAddress("l_phi",     ev.l_phi);
  t->SetBranchAddress("l_mass",    ev.l_mass);
  t->SetBranchAddress("l_chargedHadronIso",   ev.l_chargedHadronIso);
  t->SetBranchAddress("l_miniIso",          ev.l_miniIso);
  t->SetBranchAddress("l_relIso", ev.l_relIso);
  t->SetBranchAddress("l_ip3d",               ev.l_ip3d);
  t->SetBranchAddress("l_ip3dsig",            ev.l_ip3dsig);

  t->SetBranchAddress("ngenj",       &ev.ngenj);
  t->SetBranchAddress("nj",          &ev.nj);
  t->SetBranchAddress("j_area",       ev.j_area);
  t->SetBranchAddress("j_pt",         ev.j_pt);
  t->SetBranchAddress("j_eta",        ev.j_eta);
  t->SetBranchAddress("j_phi",        ev.j_phi);
  t->SetBranchAddress("j_mass",       ev.j_mass);
  t->SetBranchAddress("genj_pt",      ev.genj_pt);
  t->SetBranchAddress("genj_eta",     ev.genj_eta);
  t->SetBranchAddress("genj_phi",     ev.genj_phi);
  t->SetBranchAddress("genj_mass",    ev.genj_mass);
  t->SetBranchAddress("j_csv",        ev.j_csv);
  t->SetBranchAddress("j_vtxpy",      ev.j_vtxpx);
  t->SetBranchAddress("j_vtxpx",      ev.j_vtxpy);
  t->SetBranchAddress("j_vtxpz",      ev.j_vtxpz);
  t->SetBranchAddress("j_vtxmass",    ev.j_vtxmass);
  t->SetBranchAddress("j_vtxNtracks", ev.j_vtxNtracks);
  t->SetBranchAddress("j_vtx3DVal",   ev.j_vtx3DVal);
  t->SetBranchAddress("j_vtx3DSig",   ev.j_vtx3DSig);
  t->SetBranchAddress("j_puid",       ev.j_puid);
  t->SetBranchAddress("j_qg",         ev.j_qg);
  t->SetBranchAddress("j_flav",       ev.j_flav);
  t->SetBranchAddress("j_hadflav",    ev.j_hadflav);
  t->SetBranchAddress("j_pid",         ev.j_pid);

  t->SetBranchAddress("nmet",    &ev.nmet);
  t->SetBranchAddress("met_pt",    ev.met_pt);
  t->SetBranchAddress("met_phi",   ev.met_phi);
  
  if(t->GetBranch("npf"))
    {
      t->SetBranchAddress("npf",        &ev.npf);
      t->SetBranchAddress("pf_j",       ev.pf_j);
      t->SetBranchAddress("pf_id",      ev.pf_id);
      t->SetBranchAddress("pf_charge",  ev.pf_charge);
      t->SetBranchAddress("pf_px",      ev.pf_px);
      t->SetBranchAddress("pf_py",      ev.pf_py);
      t->SetBranchAddress("pf_pz",      ev.pf_pz);
    }
  if(t->GetBranch("ngen"))
    {
      t->SetBranchAddress("ngen",     &ev.ngen);
      t->SetBranchAddress("g_j",       ev.g_j);
      t->SetBranchAddress("g_id",      ev.g_id);
      t->SetBranchAddress("g_charge",  ev.g_charge);
      t->SetBranchAddress("g_px",      ev.g_px);
      t->SetBranchAddress("g_py",      ev.g_py);
      t->SetBranchAddress("g_pz",      ev.g_pz);
    }
}
