#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

//
void createMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  t->Branch("isData",     &ev.isData,     "isData/O");

  t->Branch("ttbar_nw",        &ev.ttbar_nw,        "ttbar_nw/I");
  t->Branch("ttbar_allmepartons",        &ev.ttbar_allmepartons,        "ttbar_allmepartons/I");
  t->Branch("ttbar_matchmepartons",        &ev.ttbar_matchmepartons,        "ttbar_matchmepartons/I");
  t->Branch("ttbar_w",        ev.ttbar_w,        "ttbar_w[ttbar_nw]/F");
  t->Branch("ttbar_genId",    &ev.ttbar_genId,    "ttbar_genId/I");

  t->Branch("run",       &ev.run,       "run/I");
  t->Branch("event",     &ev.event,     "event/I");
  t->Branch("lumi",      &ev.lumi,      "lumi/I");

  t->Branch("isFiducial",       &ev.isFiducial,       "isFiducial/O");
  t->Branch("muTrigger",        &ev.muTrigger,        "muTrigger/I");
  t->Branch("elTrigger",        &ev.elTrigger,        "elTrigger/I");

  t->Branch("nvtx",      &ev.nvtx,      "nvtx/I");
  t->Branch("pu",      &ev.pu,      "pu/I");
  t->Branch("putrue",      &ev.putrue,      "putrue/I");

  t->Branch("isPromptFinalState",        &ev.isPromptFinalState,        "isPromptFinalState/O");
  t->Branch("isDirectPromptTauDecayProductFinalState",        &ev.isDirectPromptTauDecayProductFinalState,        "isDirectPromptTauDecayProductFinalState/O");
  t->Branch("l_id",      &ev.l_id,      "l_id/I");
  t->Branch("l_charge",  &ev.l_charge,  "l_charge/I");
  t->Branch("l_pt",      &ev.l_pt,      "l_pt/F");
  t->Branch("l_eta",     &ev.l_eta,     "l_eta/F");
  t->Branch("l_phi",     &ev.l_phi,     "l_phi/F");
  t->Branch("l_mass",    &ev.l_mass,    "l_mass/F");
  t->Branch("l_chargedHadronIso",        &ev.l_chargedHadronIso,     "l_chargedHadronIso/F");
  t->Branch("l_neutralHadronIso",        &ev.l_neutralHadronIso,     "l_neutralHadronIso/F");
  t->Branch("l_photonIso",               &ev.l_photonIso,     "l_photonIso/F");
  t->Branch("l_puChargedHadronIso",      &ev.l_puChargedHadronIso,     "l_puChargedHadronIso/F");
  t->Branch("l_ip3d",      &ev.l_ip3d,      "l_ip3d/F");
  t->Branch("l_ip3dsig",      &ev.l_ip3dsig,      "l_ip3dsig/F");

  t->Branch("me_id",        &ev.me_id,        "me_id/I");
  t->Branch("me_np",        &ev.me_np,        "me_np/I");
  t->Branch("me_pid",       ev.me_pid,        "me_pid/I");
  t->Branch("me_px",       ev.me_px,        "me_px/F");
  t->Branch("me_py",       ev.me_py,        "me_py/F");
  t->Branch("me_pz",       ev.me_pz,        "me_pz/F");
  t->Branch("me_mass",       ev.me_mass,        "me_mass/F");

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

  t->Branch("met_pt",    &ev.met_pt,    "met_pt/F");
  t->Branch("met_phi",   &ev.met_phi,   "met_phi/F");
  t->Branch("mt",        &ev.mt,        "mt/F");
}

//
void attachToMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  t->SetBranchAddress("isData",     &ev.isData);

  t->SetBranchAddress("ttbar_nw",        &ev.ttbar_nw);
  t->SetBranchAddress("ttbar_allmepartons",        &ev.ttbar_allmepartons);
  t->SetBranchAddress("ttbar_matchmepartons",        &ev.ttbar_matchmepartons);
  t->SetBranchAddress("ttbar_w",        ev.ttbar_w);
  t->SetBranchAddress("ttbar_genId",    &ev.ttbar_genId);

  t->SetBranchAddress("me_id",        &ev.me_id);
  t->SetBranchAddress("me_np",        &ev.me_np);
  t->SetBranchAddress("me_pid",       ev.me_pid);
  t->SetBranchAddress("me_px",       ev.me_px);
  t->SetBranchAddress("me_py",       ev.me_py);
  t->SetBranchAddress("me_pz",       ev.me_pz);
  t->SetBranchAddress("me_mass",       ev.me_mass);

  t->SetBranchAddress("run",       &ev.run);
  t->SetBranchAddress("event",     &ev.event);
  t->SetBranchAddress("lumi",      &ev.lumi);

  t->SetBranchAddress("isFiducial",       &ev.isFiducial);
  t->SetBranchAddress("muTrigger",        &ev.muTrigger);
  t->SetBranchAddress("elTrigger",        &ev.elTrigger);

  t->SetBranchAddress("nvtx",      &ev.nvtx);
  t->SetBranchAddress("pu",      &ev.pu);
  t->SetBranchAddress("putrue",      &ev.putrue);

  t->SetBranchAddress("isPromptFinalState",        &ev.isPromptFinalState);
  t->SetBranchAddress("isDirectPromptTauDecayProductFinalState",        &ev.isDirectPromptTauDecayProductFinalState);
  t->SetBranchAddress("l_id",      &ev.l_id);
  t->SetBranchAddress("l_charge",  &ev.l_charge);
  t->SetBranchAddress("l_pt",      &ev.l_pt);
  t->SetBranchAddress("l_eta",     &ev.l_eta);
  t->SetBranchAddress("l_phi",     &ev.l_phi);
  t->SetBranchAddress("l_mass",    &ev.l_mass);
  t->SetBranchAddress("l_chargedHadronIso",        &ev.l_chargedHadronIso);
  t->SetBranchAddress("l_neutralHadronIso",        &ev.l_neutralHadronIso);
  t->SetBranchAddress("l_photonIso",               &ev.l_photonIso);
  t->SetBranchAddress("l_puChargedHadronIso",      &ev.l_puChargedHadronIso);
  t->SetBranchAddress("l_ip3d",      &ev.l_ip3d);
  t->SetBranchAddress("l_ip3dsig",      &ev.l_ip3dsig);

  t->SetBranchAddress("ngenj",        &ev.ngenj);
  t->SetBranchAddress("nj",        &ev.nj);
  t->SetBranchAddress("j_area",       ev.j_area);
  t->SetBranchAddress("j_pt",       ev.j_pt);
  t->SetBranchAddress("j_eta",      ev.j_eta);
  t->SetBranchAddress("j_phi",      ev.j_phi);
  t->SetBranchAddress("j_mass",      ev.j_mass);
  t->SetBranchAddress("genj_pt",       ev.genj_pt);
  t->SetBranchAddress("genj_eta",      ev.genj_eta);
  t->SetBranchAddress("genj_phi",      ev.genj_phi);
  t->SetBranchAddress("genj_mass",       ev.genj_mass);
  t->SetBranchAddress("j_csv",      ev.j_csv);
  t->SetBranchAddress("j_vtxpy",  ev.j_vtxpx);
  t->SetBranchAddress("j_vtxpx",  ev.j_vtxpy);
  t->SetBranchAddress("j_vtxpz",  ev.j_vtxpz);
  t->SetBranchAddress("j_vtxmass",  ev.j_vtxmass);
  t->SetBranchAddress("j_vtxNtracks",  ev.j_vtxNtracks);
  t->SetBranchAddress("j_vtx3DVal",  ev.j_vtx3DVal);
  t->SetBranchAddress("j_vtx3DSig",  ev.j_vtx3DSig);
  t->SetBranchAddress("j_puid",     ev.j_puid);
  t->SetBranchAddress("j_qg",     ev.j_qg);
  t->SetBranchAddress("j_flav",     ev.j_flav);
  t->SetBranchAddress("j_hadflav",     ev.j_hadflav);
  t->SetBranchAddress("j_pid",      ev.j_pid);

  t->SetBranchAddress("met_pt",    &ev.met_pt);
  t->SetBranchAddress("met_phi",   &ev.met_phi);
  t->SetBranchAddress("mt",        &ev.mt);
}
