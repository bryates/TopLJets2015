#include "UserCode/TopAnalysis/interface/MiniEvent.h"

//
void createMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  t->Branch("run",       &ev.run,       "run/I");
  t->Branch("event",     &ev.event,     "event/I");
  t->Branch("lumi",      &ev.lumi,      "lumi/I");
  t->Branch("isData",    &ev.isData,    "isData/O");
  t->Branch("isMC",      &ev.isMC,      "isMC/O");
  t->Branch("elTrigger", &ev.elTrigger, "elTrigger/O");
  t->Branch("muTrigger", &ev.muTrigger, "muTrigger/O");
  
  t->Branch("nvtx",      &ev.nvtx,      "nvtx/I");
  t->Branch("l_id",      &ev.l_id,      "l_id/I");
  t->Branch("l_charge",  &ev.l_charge,  "l_charge/I");
  t->Branch("l_pt",      &ev.l_pt,      "l_pt/F");
  t->Branch("l_eta",     &ev.l_eta,     "l_eta/F");
  t->Branch("l_phi",     &ev.l_phi,     "l_phi/F");
  t->Branch("l_chargedHadronIso",        &ev.l_chargedHadronIso,     "l_chargedHadronIso/F");
  t->Branch("l_neutralHadronIso",        &ev.l_neutralHadronIso,     "l_neutralHadronIso/F");
  t->Branch("l_photonIso",               &ev.l_photonIso,            "l_photonIso/F");
  t->Branch("l_puChargedHadronIso",      &ev.l_puChargedHadronIso,   "l_puChargedHadronIso/F");
  t->Branch("l_passVetoId",              &ev.l_passVetoId,           "l_passVetoId/O");
  t->Branch("l_passLooseId",             &ev.l_passLooseId,          "l_passLooseId/O");
  t->Branch("l_passMediumId",            &ev.l_passMediumId,         "l_passMediumId/O");
  t->Branch("l_passTightId",             &ev.l_passTightId,          "l_passTightId/O");
  t->Branch("l_passMediumMVAId",            &ev.l_passMediumMVAId,   "l_passMediumMVAId/O");
  t->Branch("l_passTightMVAId",             &ev.l_passTightMVAId,    "l_passTightMVAId/O");
  t->Branch("nj",           &ev.nj,        "nj/I");
  t->Branch("j_pt",         ev.j_pt,       "j_pt[nj]/F");
  t->Branch("j_energy",     ev.j_energy,   "j_energy[nj]/F");
  t->Branch("j_eta",        ev.j_eta,      "j_eta[nj]/F");
  t->Branch("j_phi",        ev.j_phi,      "j_phi[nj]/F");
  t->Branch("j_nch",        ev.j_nch,      "j_nch[nj]/I");
  t->Branch("j_chpt",       ev.j_chpt,     "j_chpt[nj]/F");
  t->Branch("j_chsumpt",    ev.j_chsumpt,  "j_chsumpt[nj]/F");
  t->Branch("j_cheta",      ev.j_cheta,    "j_cheta[nj]/F");
  t->Branch("j_chphi",      ev.j_chphi,    "j_chphi[nj]/F");
  t->Branch("j_csv",        ev.j_csv,      "j_csv[nj]/F");
  t->Branch("j_vtxmass",    ev.j_vtxmass,  "j_vtxmass[nj]/F");
  t->Branch("j_vtxNtracks", ev.j_vtxNtracks, "j_vtxNtracks[nj]/I");
  t->Branch("j_vtx3DVal",   ev.j_vtx3DVal, "j_vtx3DVal[nj]/F");
  t->Branch("j_vtx3DSig",   ev.j_vtx3DSig, "j_vtx3DSig[nj]/F");
  t->Branch("j_puid",       ev.j_puid,     "j_puid[nj]/F");
  t->Branch("j_flav",       ev.j_flav,     "j_flav[nj]/I");
  t->Branch("j_pid",        ev.j_pid,      "j_pid[nj]/I");

  t->Branch("met_pt",       &ev.met_pt,    "met_pt/F");
  t->Branch("met_phi",      &ev.met_phi,   "met_phi/F");
  t->Branch("mt",           &ev.mt,        "mt/F");
  t->Branch("chmet_pt",     &ev.chmet_pt,    "chmet_pt/F");
  t->Branch("chmet_phi",    &ev.chmet_phi,   "chmet_phi/F");
  t->Branch("chmt",         &ev.chmt,        "chmt/F");
}

//
void attachToMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  t->SetBranchAddress("run",       &ev.run);
  t->SetBranchAddress("event",     &ev.event);
  t->SetBranchAddress("lumi",      &ev.lumi);
  t->SetBranchAddress("isData",    &ev.isData);
  t->SetBranchAddress("isMC",      &ev.isMC);
  t->SetBranchAddress("elTrigger", &ev.elTrigger);
  t->SetBranchAddress("muTrigger", &ev.muTrigger);
  
  t->SetBranchAddress("nvtx",      &ev.nvtx);
  t->SetBranchAddress("l_id",      &ev.l_id);
  t->SetBranchAddress("l_charge",  &ev.l_charge);
  t->SetBranchAddress("l_pt",      &ev.l_pt);
  t->SetBranchAddress("l_eta",     &ev.l_eta);
  t->SetBranchAddress("l_phi",     &ev.l_phi);
  t->SetBranchAddress("l_chargedHadronIso",     &ev.l_chargedHadronIso);
  t->SetBranchAddress("l_neutralHadronIso",     &ev.l_neutralHadronIso);
  t->SetBranchAddress("l_photonIso",            &ev.l_photonIso);
  t->SetBranchAddress("l_puChargedHadronIso",   &ev.l_puChargedHadronIso);
  t->SetBranchAddress("l_passVetoId",           &ev.l_passVetoId);
  t->SetBranchAddress("l_passLooseId",          &ev.l_passLooseId);
  t->SetBranchAddress("l_passMediumId",         &ev.l_passMediumId);
  t->SetBranchAddress("l_passTightId",          &ev.l_passTightId);
  t->SetBranchAddress("l_passMediumMVAId",      &ev.l_passMediumMVAId);
  t->SetBranchAddress("l_passTightMVAId",       &ev.l_passTightMVAId);
  t->SetBranchAddress("nj",        &ev.nj);
  t->SetBranchAddress("j_pt",       ev.j_pt);
  t->SetBranchAddress("j_energy",   ev.j_energy);
  t->SetBranchAddress("j_eta",      ev.j_eta);
  t->SetBranchAddress("j_phi",      ev.j_phi);
  t->SetBranchAddress("j_nch",      ev.j_nch);
  t->SetBranchAddress("j_chpt",     ev.j_chpt);
  t->SetBranchAddress("j_chsumpt",  ev.j_chsumpt);
  t->SetBranchAddress("j_cheta",    ev.j_cheta);
  t->SetBranchAddress("j_chphi",    ev.j_chphi);
  t->SetBranchAddress("j_csv",      ev.j_csv);
  t->SetBranchAddress("j_vtxmass",  ev.j_vtxmass);
  t->SetBranchAddress("j_vtxNtracks",  ev.j_vtxNtracks);
  t->SetBranchAddress("j_vtx3DVal",  ev.j_vtx3DVal);
  t->SetBranchAddress("j_vtx3DSig",  ev.j_vtx3DSig);
  t->SetBranchAddress("j_puid",     ev.j_puid);
  t->SetBranchAddress("j_flav",     ev.j_flav);
  t->SetBranchAddress("j_pid",      ev.j_pid);
  t->SetBranchAddress("met_pt",    &ev.met_pt);
  t->SetBranchAddress("met_phi",   &ev.met_phi);
  t->SetBranchAddress("mt",        &ev.mt);
  t->SetBranchAddress("chmet_pt",    &ev.chmet_pt);
  t->SetBranchAddress("chmet_phi",   &ev.chmet_phi);
  t->SetBranchAddress("chmt",        &ev.chmt);
}
