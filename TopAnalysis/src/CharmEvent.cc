#include "TopLJets2015/TopAnalysis/interface/CharmEvent.h"

//
void createCharmEventTree(TTree *t, CharmEvent_t &ev)
{
  //event header
  t->Branch("isData",       &ev.isData,       "isData/O");
  t->Branch("run",          &ev.run,          "run/I");
  t->Branch("event",        &ev.event,        "event/I");
  t->Branch("norm",         &ev.norm,         "norm/F");
  t->Branch("xsec",         &ev.xsec,         "xsec/F");
  t->Branch("lumi",         &ev.lumi,         "lumi/F");
  t->Branch("lum",          &ev.lum,          "lum/F");
  t->Branch("topptwgt",     &ev.topptwgt,     "topptwgt/F");


  //generator level event
  t->Branch("ttbar_nw",      &ev.ttbar_nw,       "ttbar_nw/I");
  t->Branch("ttbar_w",        ev.ttbar_w,        "ttbar_w[ttbar_nw]/F");
  t->Branch("ngfsr",         &ev.ngfsr,          "ngfsr/I");
  t->Branch("gfsr",           ev.gfsr,           "gfsr[ngfsr]/F");

  //reco level event
  t->Branch("nvtx",      &ev.nvtx,      "nvtx/I");
  t->Branch("HT",           &ev.ht,           "H_{T}/F");
  t->Branch("MET",          &ev.met,          "MET/F");
  t->Branch("ST",           &ev.st,           "S_{T}/F");
  t->Branch("nlj",          &ev.nlj,          "nlj/I");

  t->Branch("nleptons",  &ev.nleptons, "nleptons/I");
  t->Branch("nl",        &ev.nl,       "nl/I");
  t->Branch("l_id",       ev.l_id,      "l_id[nl]/I");
  t->Branch("l_pid",      ev.l_pid,     "l_pid[nl]/I");
  t->Branch("l_g",        ev.l_g,       "l_g[nl]/I");
  t->Branch("l_charge",   ev.l_charge,  "l_charge[nl]/I");
  t->Branch("l_pt",       ev.l_pt,      "l_pt[nl]/F");
  t->Branch("l_eta",      ev.l_eta,     "l_eta[nl]/F");
  t->Branch("l_phi",      ev.l_phi,     "l_phi[nl]/F");
  t->Branch("l_mass",     ev.l_mass,    "l_mass[nl]/F");
  t->Branch("l_chi2norm",         ev.l_chi2norm,         "l_chi2norm[nl]/F");
  t->Branch("l_dxy",              ev.l_dxy,              "dxy[nl]/F");
  t->Branch("l_dxyE",             ev.l_dxyE,             "dxyE[nl]/F");
  t->Branch("l_dz",               ev.l_dz,               "dz[nl]/F");
  t->Branch("l_dzE",              ev.l_dzE,              "dzE[nl]/F");

  //jet info
  t->Branch("nj",        &ev.nj,        "nj/I");
  t->Branch("j_pt",       ev.j_pt,      "j_pt[nj]/F");
  t->Branch("j_pt_pf",    ev.j_pt_pf,   "j_pt_pf[nj]/F");
  t->Branch("j_pt_charged", ev.j_pt_charged, "j_pt_charged[nj]/F");
  t->Branch("j_eta",      ev.j_eta,     "j_eta[nj]/F");
  t->Branch("j_phi",      ev.j_phi,     "j_phi[nj]/F");
  t->Branch("j_p",        ev.j_p,       "j_p[nj]/F");
  t->Branch("j_p_pf",     ev.j_p_pf,    "j_p_pf[nj]/F");
  t->Branch("j_p_charged",  ev.j_p_charged,   "j_p_charged[nj]/F");
  t->Branch("j_pz",       ev.j_pz,      "j_pz[nj]/F");
  t->Branch("j_pz_pf",     ev.j_pz_pf,    "j_pz_pf[nj]/F");
  t->Branch("j_pz_charged",  ev.j_pz_charged,   "j_pz_charged[nj]/F");
  t->Branch("j_mass",     ev.j_mass,     "j_mass[nj]/F");
  t->Branch("j_csv",      ev.j_csv,     "j_csv[nj]/F");
  t->Branch("j_flav",     ev.j_flav,    "j_flav[nj]/I");
  t->Branch("j_hadflav",     ev.j_hadflav,    "j_hadflav[nj]/I");
  t->Branch("j_pid",     ev.j_pid,    "j_pid[nj]/I");
  /*
  t->Branch("gj_pt",       ev.gj_pt,      "gj_pt[nj]/F");
  t->Branch("gj_eta",      ev.gj_eta,     "gj_eta[nj]/F");
  t->Branch("gj_phi",      ev.gj_phi,     "gj_phi[nj]/F");
  */

  //pf candidates (only charged if outside jets)
  t->Branch("npf",        &ev.npf,         "npf/I");
  t->Branch("pf_j",        ev.pf_j,        "pf_j[npf]/I");
  t->Branch("pf_id",       ev.pf_id,       "pf_id[npf]/I");
  t->Branch("pf_c",        ev.pf_c,        "pf_c[npf]/I");
  t->Branch("pf_pt",       ev.pf_pt,       "pf_pt[npf]/F");
  t->Branch("pf_eta",      ev.pf_eta,      "pf_eta[npf]/F");
  t->Branch("pf_phi",      ev.pf_phi,      "pf_phi[npf]/F");
  t->Branch("pf_m",        ev.pf_m,        "pf_m[npf]/F");
  t->Branch("pf_dxy",      ev.pf_dxy,      "pf_dxy[npf]/F");
  t->Branch("pf_dxyE",     ev.pf_dxyE,     "pf_dxyE[npf]/F");
  t->Branch("pf_dz",       ev.pf_dz,       "pf_dz[npf]/F");
  t->Branch("pf_dzE",      ev.pf_dzE,      "pf_dzE[npf]/F");

  //JPsi candidates
  t->Branch("njpsi",        &ev.njpsi,        "njpsi/I");
  t->Branch("nmeson",       &ev.nmeson,       "nmeson/I");
  t->Branch("nk",           &ev.nk,           "nk/I");
  t->Branch("meson_id",      ev.meson_id,     "meson_id[nmeson]/I");
  t->Branch("epoch" ,        ev.epoch,        "epoch[nmeson]/I");
  t->Branch("puwgt",         ev.puwgt,        "puwgt[nmeson]/F");
  t->Branch("sfs",           ev.sfs,          "sfs[nmeson]/F");
  t->Branch("sfsu",          ev.sfsu,         "sfsu[nmeson]/F");
  t->Branch("pitrk",         ev.pitrk,        "pitrk[nmeson]/F");
  t->Branch("piptsf",        ev.piptsf,       "piptsf[nmeson]/F");
  t->Branch("kptsf",         ev.kptsf,        "kptsf[nmeson]/F");

  t->Branch("jpsi_mass",    ev.jpsi_mass,     "jpsi_mass[nmeson]/F");
  t->Branch("jpsi_b_mass",  ev.jpsi_b_mass,   "jpsi_b_mass[nk]/F");
  t->Branch("jpsi_pt",      ev.jpsi_pt,       "jpsi_pt[nmeson]/F");
  t->Branch("jpsi_eta",     ev.jpsi_eta,      "jpsi_eta[nmeson]/F");
  t->Branch("jpsi_phi",     ev.jpsi_phi,      "jpsi_phi[nmeson]/F");
  t->Branch("jpsi_p",       ev.jpsi_p,        "jpsi_p[nmeson]/F");
  t->Branch("jpsi_pz",      ev.jpsi_pz,       "jpsi_pz[nmeson]/F");
  t->Branch("jpsi_mu1_pt",  ev.jpsi_mu1_pt,   "jpsi_mu1_pt[nmeson]/F");
  t->Branch("jpsi_mu1_eta", ev.jpsi_mu1_eta,  "jpsi_mu1_eta[nmeson]/F");
  t->Branch("jpsi_mu1_phi", ev.jpsi_mu1_phi,  "jpsi_mu1_phi[nmeson]/F");
  t->Branch("jpsi_mu2_pt",  ev.jpsi_mu2_pt,   "jpsi_mu2_pt[nmeson]/F");
  t->Branch("jpsi_mu2_eta", ev.jpsi_mu2_eta,  "jpsi_mu2_eta[nmeson]/F");
  t->Branch("jpsi_mu2_phi", ev.jpsi_mu2_phi,  "jpsi_mu2_phi[nmeson]/F");
  t->Branch("jpsi_k_pt",   ev.jpsi_k_pt,      "jpsi_k_pt[nmeson]/F");
  t->Branch("jpsi_k_eta",  ev.jpsi_k_eta,     "jpsi_k_eta[nmeson]/F");
  t->Branch("jpsi_k_phi",  ev.jpsi_k_phi,     "jpsi_k_phi[nmeson]/F");
  t->Branch("jpsi_k_dxy",  ev.jpsi_k_dxy,     "jpsi_k_dxy[nmeson]/F");
  t->Branch("jpsi_k_dz",   ev.jpsi_k_dz,      "jpsi_k_dz[nmeson]/F");
  t->Branch("jpsi_chi2",   ev.jpsi_chi2,      "jpsi_chi2[nmeson]/F");
  t->Branch("jpsi_j",       ev.jpsi_j,        "jpsi_j[nmeson]/F");
  t->Branch("jpsi_ptrel",   ev.jpsi_ptrel,    "jpsi_ptrel[nmeson]/F");
  t->Branch("jpsi_l",       ev.jpsi_l,        "jpsi_l[nmeson]/F");
  t->Branch("jpsi_l_mass",  ev.jpsi_l_mass,   "jpsi_l_mass[nmeson]/F");
  t->Branch("jpsi_l_dR",    ev.jpsi_l_dR,     "jpsi_l_dR[nmeson]/F");
  t->Branch("jpsi_l3d",     ev.jpsi_l3d,      "jpsi_l3d[nmeson]/F");
  t->Branch("jpsi_lx",        ev.jpsi_lx,         "jpsi_lx[nmeson]/F");
  t->Branch("jpsi_ly",        ev.jpsi_ly,         "jpsi_ly[nmeson]/F");
  t->Branch("jpsi_lz",        ev.jpsi_lz,         "jpsi_lz[nmeson]/F");
  t->Branch("jpsi_sigmal3d",  ev.jpsi_sigmal3d,   "jpsi_sigmal3d[nmeson]/F");
  t->Branch("jpsi_sigmax",    ev.jpsi_sigmax,     "jpsi_sigmax[nmeson]/F");
  t->Branch("jpsi_sigmay",    ev.jpsi_sigmay,     "jpsi_sigmay[nmeson]/F");
  t->Branch("jpsi_sigmaz",    ev.jpsi_sigmaz,     "jpsi_sigmaz[nmeson]/F");

  //D meson candidates
  t->Branch("d0_mass",       ev.d0_mass,        "d0_mass[nmeson]/F");
  t->Branch("d0_pt",         ev.d0_pt,          "d0_pt[nmeson]/F");
  t->Branch("d0_eta",        ev.d0_eta,         "d0_eta[nmeson]/F");
  t->Branch("d0_phi",        ev.d0_phi,         "d0_phi[nmeson]/F");
  t->Branch("d0_p",          ev.d0_p,           "d0_p[nmeson]/F");
  t->Branch("d0_pz",         ev.d0_pz,          "d0_pz[nmeson]/F");
  t->Branch("ds_mass",       ev.ds_mass,        "ds_mass[nmeson]/F");
  t->Branch("ds_pt",         ev.ds_pt,          "ds_pt[nmeson]/F");
  t->Branch("ds_eta",        ev.ds_eta,         "ds_eta[nmeson]/F");
  t->Branch("ds_phi",        ev.ds_phi,         "ds_phi[nmeson]/F");
  t->Branch("ds_p",          ev.ds_p,           "ds_p[nmeson]/F");
  t->Branch("ds_pz",         ev.ds_pz,          "ds_pz[nmeson]/F");
  t->Branch("d0_pi_pt",      ev.d0_pi_pt,       "d0_pi_pt[nmeson]/F");
  t->Branch("d0_pi_eta",     ev.d0_pi_eta,      "d0_pi_eta[nmeson]/F");
  t->Branch("d0_pi_phi",     ev.d0_pi_phi,      "d0_pi_phi[nmeson]/F");
  t->Branch("d0_pi_dxy",     ev.d0_pi_dxy,      "d0_pi_dxy[nmeson]/F");
  t->Branch("d0_pi_dxyE",    ev.d0_pi_dxyE,     "d0_pi_dxyE[nmeson]/F");
  t->Branch("d0_pi_dz",      ev.d0_pi_dz,       "d0_pi_dz[nmeson]/F");
  t->Branch("d0_k_pt",       ev.d0_k_pt,        "d0_k_pt[nmeson]/F");
  t->Branch("d0_k_eta",      ev.d0_k_eta,       "d0_k_eta[nmeson]/F");
  t->Branch("d0_k_phi",      ev.d0_k_phi,       "d0_k_phi[nmeson]/F");
  t->Branch("d0_k_dxy",      ev.d0_k_dxy,       "d0_k_dxy[nmeson]/F");
  t->Branch("d0_k_dxyE",     ev.d0_k_dxyE,      "d0_k_dxyE[nmeson]/F");
  t->Branch("d0_k_dz",       ev.d0_k_dz,        "d0_k_dz[nmeson]/F");
  t->Branch("d0_chi2",       ev.d0_chi2,        "d0_chi2[nmeson]/F");
  t->Branch("d0_opang",      ev.d0_opang,       "d0_opang[nmeson]/F");
  t->Branch("d0_mu_pt",      ev.d0_mu_pt,       "d0_mu_pt[nmeson]/F");
  t->Branch("d0_mu_eta",     ev.d0_mu_eta,      "d0_mu_eta[nmeson]/F");
  t->Branch("d0_mu_phi",     ev.d0_mu_phi,      "d0_mu_phi[nmeson]/F");
  t->Branch("d0_mu_tag_mu_pt",      ev.d0_mu_tag_mu_pt,       "d0_mu_tag_mu_pt[nmeson]/F");
  t->Branch("ds_pi2_pt",     ev.ds_pi2_pt,      "ds_pi2_pt[nmeson]/F");
  t->Branch("ds_pi2_eta",    ev.ds_pi2_eta,     "ds_pi2_eta[nmeson]/F");
  t->Branch("ds_pi2_phi",    ev.ds_pi2_phi,     "ds_pi2_phi[nmeson]/F");
  t->Branch("d0_j",          ev.d0_j,           "d0_j[nmeson]/I");
  t->Branch("d0_ptrel",      ev.d0_ptrel,       "d0_ptrel[nmeson]/F");
  t->Branch("d0_ndau",       ev.d0_ndau,        "d0_ndau[nmeson]/I");
  t->Branch("d0_l",          ev.d0_l,           "d0_l[nmeson]/I");
  t->Branch("d0_l_mass",     ev.d0_l_mass,      "d0_l_mass[nmeson]/F");
  t->Branch("d0_l_dR",       ev.d0_l_dR,        "d0_l_dR[nmeson]/F");
  t->Branch("d0_l3d",        ev.d0_l3d,         "d0_l3d[nmeson]/F");
  t->Branch("d0_lx",         ev.d0_lx,          "d0_lx[nmeson]/F");
  t->Branch("d0_ly",         ev.d0_ly,          "d0_ly[nmeson]/F");
  t->Branch("d0_lz",         ev.d0_lz,          "d0_lz[nmeson]/F");
  t->Branch("d0_sigmal3d",   ev.d0_sigmal3d,    "d0_sigmal3d[nmeson]/F");
  t->Branch("d0_sigmax",     ev.d0_sigmax,      "d0_sigmax[nmeson]/F");
  t->Branch("d0_sigmay",     ev.d0_sigmay,      "d0_sigmay[nmeson]/F");
  t->Branch("d0_sigmaz",     ev.d0_sigmaz,      "d0_sigmaz[nmeson]/F");
  t->Branch("d0_pi_mother",  ev.d0_pi_mother,   "d0_pi_mother[nmeson]/F");
  t->Branch("d0_k_mother",   ev.d0_k_mother,    "d0_k_mother[nmeson]/F");
  t->Branch("d0_pi_gpt",     ev.d0_pi_gpt,      "d0_pi_gpt[nmeson]/F");
  t->Branch("d0_k_gpt",      ev.d0_k_gpt,       "d0_k_gpt[nmeson]/F");
  t->Branch("d0_pi_gpid",    ev.d0_pi_gpid,     "d0_pi_gpid[nmeson]/F");
  t->Branch("d0_k_gpid",     ev.d0_k_gpid,      "d0_k_gpid[nmeson]/F");
  t->Branch("d0_gen",        ev.d0_gen,         "d0_gen[nmeson]/O");

  //Fragmentation
  t->Branch("peterson",    ev.peterson,     "peterson[nj]/F");
  t->Branch("up",          ev.up,           "up[nj]/F");
  t->Branch("central",     ev.central,      "central[nj]/F");
  t->Branch("down",        ev.down,         "down[nj]/F");
}

//
void attachToCharmEventTree(TTree *t, CharmEvent_t &ev)
{
  //event header
  t->SetBranchAddress("isData",    &ev.isData);
  t->SetBranchAddress("run",       &ev.run);
  t->SetBranchAddress("event",     &ev.event);
  t->SetBranchAddress("lumi",      &ev.lumi);
  if(t->GetListOfBranches()->FindObject("lum"))
  t->SetBranchAddress("lum",       &ev.lum);
  t->SetBranchAddress("epoch",      ev.epoch);
  t->SetBranchAddress("norm",      &ev.norm);
  t->SetBranchAddress("xsec",      &ev.xsec);
  t->SetBranchAddress("puwgt",      ev.puwgt);
  t->SetBranchAddress("topptwgt",  &ev.topptwgt);
  t->SetBranchAddress("sfs",        ev.sfs);
  t->SetBranchAddress("sfsu",       ev.sfsu);
  t->SetBranchAddress("pitrk",        ev.pitrk);

  //generator level event
  t->SetBranchAddress("ttbar_nw",      &ev.ttbar_nw);
  t->SetBranchAddress("ttbar_w",        ev.ttbar_w);
  if(t->GetListOfBranches()->FindObject("ngfsr")) {
  t->SetBranchAddress("ngfsr",         &ev.ngfsr);
  t->SetBranchAddress("gfsr",           ev.gfsr);
  }

  //reco level event
  t->SetBranchAddress("nvtx",      &ev.nvtx);
  t->SetBranchAddress("HT",        &ev.ht);
  t->SetBranchAddress("MET",       &ev.met);
  t->SetBranchAddress("ST",        &ev.st);
  t->SetBranchAddress("nlj",       &ev.nlj);

  t->SetBranchAddress("nleptons",  &ev.nleptons);
  t->SetBranchAddress("nl",        &ev.nl);
  t->SetBranchAddress("l_id",       ev.l_id);
  t->SetBranchAddress("l_pid",      ev.l_pid);
  t->SetBranchAddress("l_g",        ev.l_g);
  t->SetBranchAddress("l_charge",   ev.l_charge);
  t->SetBranchAddress("l_pt",       ev.l_pt);
  t->SetBranchAddress("l_eta",      ev.l_eta);
  t->SetBranchAddress("l_phi",      ev.l_phi);
  t->SetBranchAddress("l_mass",     ev.l_mass);
  t->SetBranchAddress("l_dxy",              ev.l_dxy);
  t->SetBranchAddress("l_dz",               ev.l_dz);
  t->SetBranchAddress("l_dzE",              ev.l_dzE);

  //jet info
  t->SetBranchAddress("nj",        &ev.nj);
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
  if(t->GetListOfBranches()->FindObject("j_ntk"))
    t->SetBranchAddress("j_ntk",      ev.j_ntk);
  t->SetBranchAddress("j_csv",      ev.j_csv);
  t->SetBranchAddress("j_hadflav",     ev.j_hadflav);
  /*
  t->SetBranchAddress("gj_pt",       ev.gj_pt);
  t->SetBranchAddress("gj_eta",      ev.gj_eta);
  t->SetBranchAddress("gj_phi",      ev.gj_phi);
  */

  //pf candidates (only charged if outside jets)
  t->SetBranchAddress("npf",        &ev.npf);
  t->SetBranchAddress("pf_j",        ev.pf_j);
  t->SetBranchAddress("pf_id",       ev.pf_id);
  t->SetBranchAddress("pf_c",        ev.pf_c);
  t->SetBranchAddress("pf_pt",       ev.pf_pt);
  t->SetBranchAddress("pf_eta",      ev.pf_eta);
  t->SetBranchAddress("pf_phi",      ev.pf_phi);
  t->SetBranchAddress("pf_m",        ev.pf_m);

  t->SetBranchAddress("pf_dxy",      ev.pf_dxy);
  t->SetBranchAddress("pf_dxyE",     ev.pf_dxyE);
  t->SetBranchAddress("pf_dz",       ev.pf_dz);
  t->SetBranchAddress("pf_dzE",      ev.pf_dzE);

  //JPsi candidates
  t->SetBranchAddress("njpsi",       &ev.njpsi);
  t->SetBranchAddress("jpsi_mass",    ev.jpsi_mass);
  t->SetBranchAddress("jpsi_b_mass",  ev.jpsi_b_mass);
  t->SetBranchAddress("jpsi_pt",      ev.jpsi_pt);
  t->SetBranchAddress("jpsi_eta",     ev.jpsi_eta);
  t->SetBranchAddress("jpsi_phi",     ev.jpsi_phi);
  t->SetBranchAddress("jpsi_p",       ev.jpsi_p);
  t->SetBranchAddress("jpsi_pz",      ev.jpsi_pz);
  t->SetBranchAddress("jpsi_eta",     ev.jpsi_eta);
  t->SetBranchAddress("jpsi_eta",     ev.jpsi_eta);
  t->SetBranchAddress("jpsi_mu1_pt",  ev.jpsi_mu1_pt);
  t->SetBranchAddress("jpsi_mu1_eta", ev.jpsi_mu1_eta);
  t->SetBranchAddress("jpsi_mu1_phi", ev.jpsi_mu1_phi);
  t->SetBranchAddress("jpsi_mu2_pt",  ev.jpsi_mu2_pt);
  t->SetBranchAddress("jpsi_mu2_eta", ev.jpsi_mu2_eta);
  t->SetBranchAddress("jpsi_mu2_phi", ev.jpsi_mu2_phi);
  t->SetBranchAddress("jpsi_k_pt",    ev.jpsi_k_pt);
  t->SetBranchAddress("jpsi_k_eta",   ev.jpsi_k_eta);
  t->SetBranchAddress("jpsi_k_phi",   ev.jpsi_k_phi);
  t->SetBranchAddress("jpsi_k_dxy",   ev.jpsi_k_dxy);
  t->SetBranchAddress("jpsi_k_dz",    ev.jpsi_k_dz);
  t->SetBranchAddress("jpsi_chi2",    ev.jpsi_chi2);
  t->SetBranchAddress("jpsi_j",       ev.jpsi_j);
  t->SetBranchAddress("jpsi_ptrel",   ev.jpsi_ptrel);
  t->SetBranchAddress("jpsi_l",       ev.jpsi_l);
  t->SetBranchAddress("jpsi_l_mass",  ev.jpsi_l_mass);
  t->SetBranchAddress("jpsi_l_dR",    ev.jpsi_l_dR);
  t->SetBranchAddress("jpsi_l3d",     ev.jpsi_l3d);
  t->SetBranchAddress("jpsi_lx",        ev.jpsi_lx);
  t->SetBranchAddress("jpsi_ly",        ev.jpsi_ly);
  t->SetBranchAddress("jpsi_lz",        ev.jpsi_lz);
  t->SetBranchAddress("jpsi_sigmal3d",  ev.jpsi_sigmal3d);
  t->SetBranchAddress("jpsi_sigmal3d",  ev.jpsi_sigmal3d);
  t->SetBranchAddress("jpsi_sigmax",    ev.jpsi_sigmax);
  t->SetBranchAddress("jpsi_sigmay",    ev.jpsi_sigmay);
  t->SetBranchAddress("jpsi_sigmaz",    ev.jpsi_sigmaz);

  //D meson candidates
  t->SetBranchAddress("nmeson",       &ev.nmeson);
  t->SetBranchAddress("nk",           &ev.nk);
  t->SetBranchAddress("meson_id",      ev.meson_id);
  t->SetBranchAddress("d0_mass",       ev.d0_mass);
  t->SetBranchAddress("d0_pt",         ev.d0_pt);
  t->SetBranchAddress("d0_eta",        ev.d0_eta);
  t->SetBranchAddress("d0_phi",        ev.d0_phi);
  t->SetBranchAddress("d0_p",          ev.d0_p);
  t->SetBranchAddress("d0_pz",         ev.d0_pz);
  t->SetBranchAddress("ds_mass",       ev.ds_mass);
  t->SetBranchAddress("ds_pt",         ev.ds_pt);
  t->SetBranchAddress("ds_eta",        ev.ds_eta);
  t->SetBranchAddress("ds_phi",        ev.ds_phi);
  t->SetBranchAddress("ds_p",          ev.ds_p);
  t->SetBranchAddress("ds_pz",         ev.ds_pz);
  t->SetBranchAddress("d0_pi_pt",      ev.d0_pi_pt);
  t->SetBranchAddress("d0_pi_eta",     ev.d0_pi_eta);
  t->SetBranchAddress("d0_pi_phi",     ev.d0_pi_phi);
  t->SetBranchAddress("d0_pi_dxy",     ev.d0_pi_dxy);
  t->SetBranchAddress("d0_pi_dxyE",    ev.d0_pi_dxyE);
  t->SetBranchAddress("d0_pi_dz",      ev.d0_pi_dz);
  t->SetBranchAddress("d0_k_pt",       ev.d0_k_pt);
  t->SetBranchAddress("d0_k_eta",      ev.d0_k_eta);
  t->SetBranchAddress("d0_k_phi",      ev.d0_k_phi);
  t->SetBranchAddress("d0_k_dxy",      ev.d0_k_dxy);
  t->SetBranchAddress("d0_k_dxyE",     ev.d0_k_dxyE);
  t->SetBranchAddress("d0_k_dz",       ev.d0_k_dz);
  t->SetBranchAddress("d0_chi2",       ev.d0_chi2);
  t->SetBranchAddress("d0_opang",      ev.d0_opang);
  t->SetBranchAddress("d0_mu_pt",      ev.d0_mu_pt);
  t->SetBranchAddress("d0_mu_eta",     ev.d0_mu_eta);
  t->SetBranchAddress("d0_mu_phi",     ev.d0_mu_phi);
  t->SetBranchAddress("d0_mu_tag_mu_pt",  ev.d0_mu_tag_mu_pt);
  t->SetBranchAddress("ds_pi2_pt",     ev.ds_pi2_pt);
  t->SetBranchAddress("ds_pi2_eta",    ev.ds_pi2_eta);
  t->SetBranchAddress("ds_pi2_phi",    ev.ds_pi2_phi);
  t->SetBranchAddress("d0_j",          ev.d0_j);
  t->SetBranchAddress("d0_ptrel",      ev.d0_ptrel);
  t->SetBranchAddress("d0_ndau",       ev.d0_ndau);
  t->SetBranchAddress("d0_l",          ev.d0_l);
  t->SetBranchAddress("d0_l_mass",     ev.d0_l_mass);
  t->SetBranchAddress("d0_l_dR",       ev.d0_l_dR);
  t->SetBranchAddress("d0_l3d",        ev.d0_l3d);
  t->SetBranchAddress("d0_lx",         ev.d0_lx);
  t->SetBranchAddress("d0_ly",         ev.d0_ly);
  t->SetBranchAddress("d0_lz",         ev.d0_lz);
  t->SetBranchAddress("d0_sigmal3d",   ev.d0_sigmal3d);
  t->SetBranchAddress("d0_sigmax",     ev.d0_sigmax);
  t->SetBranchAddress("d0_sigmay",     ev.d0_sigmay);
  t->SetBranchAddress("d0_sigmaz",     ev.d0_sigmaz);
  t->SetBranchAddress("d0_pi_mother",  ev.d0_pi_mother);
  t->SetBranchAddress("d0_k_mother",   ev.d0_k_mother);
  t->SetBranchAddress("d0_pi_gpt",  ev.d0_pi_gpt);
  t->SetBranchAddress("d0_k_gpt",   ev.d0_k_gpt);
  t->SetBranchAddress("d0_pi_gpid",  ev.d0_pi_gpid);
  t->SetBranchAddress("d0_k_gpid",   ev.d0_k_gpid);
  t->SetBranchAddress("d0_gen",      ev.d0_gen);

  //Fragmentation
  t->SetBranchAddress("peterson",    ev.peterson);
  t->SetBranchAddress("up",          ev.up);
  t->SetBranchAddress("central",     ev.central);
  t->SetBranchAddress("down",        ev.down);
}
