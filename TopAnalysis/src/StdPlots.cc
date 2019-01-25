#include <iostream>
#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"
#include "TopLJets2015/TopAnalysis/interface/StdPlots.h"
#include "TopLJets2015/TopAnalysis/interface/Jet.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"

StdPlots::StdPlots(TString runPeriod, TString name, bool debug) {
  runPeriod_ = "_"+runPeriod;
  name_ = name;
  //Check for multiple run periods
  if(runPeriod.Length() > 1) {
    for(int i=0; i < runPeriod.Length(); i++) {
      TString tmp(runPeriod[i]);
      if(name.Contains("MC") || (name.Contains("Data") && name.Contains("2016"+tmp))) {
        isGood_ = true;
        break;
      }
    }
  }
  //Only use MC or Data with correct run period
  else if(name.Contains("MC") || (name.Contains("Data") && name.Contains("2016"+runPeriod))) isGood_ = true;
  else isGood_ = false;
  debug_ = debug;
  norm_ = 1.;
  sfs_.first = 1.; sfs_.second.first = 0.; sfs_.second.second = 0.;
  puWgt_ = 1.;
  top_pt_wgt_ = 1.;
  tracker_wgt_ = 1.;
  pi_wgt_.first = 1.; pi_wgt_.second.first = 0.; pi_wgt_.second.second = 0.;
  rbWgt_ = 1.;

  if(debug_ && isGood_)
    std::cout << "Initializing run" << runPeriod_ << std::endl;
  if(!isGood_) return;

  //PU plot
  allPlots["puwgtctr"+runPeriod_] = new TH1F("puwgtctr"+runPeriod_,"Weight sums",4,0,4);
  allPlots["topptwgt"+runPeriod_] = new TH1F("topptwgt"+runPeriod_,"Top p_{T} weights", 2, 0, 2);
  allPlots["rbwgt_jpsi"+runPeriod_] = new TH1F("rbwgt_jpsi"+runPeriod_,"r_{B} event weights", 2, 0, 2);
  allPlots["rbwgt_d0"+runPeriod_] = new TH1F("rbwgt_d0"+runPeriod_,"r_{B} event weights", 2, 0, 2);
  allPlots["rbwgt_d0mu"+runPeriod_] = new TH1F("rbwgt_d0mu"+runPeriod_,"r_{B} event weights", 2, 0, 2);
  std::vector<TString> lfsVec = { "_all", "_e", "_ee", "_em", "_mm", "_m" }; 
  std::vector<TString> cutVec = { "", "_lep", "_lepjets", "_csv", "_jpsi", "_gjpsi", "_rgjpsi", "_kjpsi", "_meson", "_gmeson", "_rgmeson" };
  if(name_.Contains("Data")) //Gen plots only for MC
    cutVec = { "", "_lep", "_lepjets", "_csv", "_jpsi", "_kjpsi", "_meson" };
  std::vector<TString> wgtVec = { "" }; //, "_no_weight" };
  for(int i = 0; i < (int)lfsVec.size(); i++) {
  for(int j = 0; j < (int)cutVec.size(); j++) {
  for(int k = 0; k < (int)wgtVec.size(); k++) {
    TString tag(lfsVec[i]);
    TString cut(cutVec[j]);
    TString weight(wgtVec[k]);

    //Plots to hold up/down uncertainties
    //allPlots["unc"+tag+cut+weight+runPeriod_] = new TH1F("unc"+tag+cut+weight+runPeriod_,";Uncertainties (up,down);Events", 2, 0,2);

    // Lepton plots
    allPlots["lp_pt_iso"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt_iso"+tag+cut+weight+runPeriod_,";Lepton P_{T} [GeV] after cleaning;Events / 10 GeV", 50, 0,500);
    allPlots["lp_pt_veto"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt_veto"+tag+cut+weight+runPeriod_,";Lepton P_{T} [GeV] after veto;Events / 10 GeV", 60, 0,600);
    allPlots["lp_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt_low"+tag+cut+weight+runPeriod_,";Leading Lepton P_{T} [GeV];Events / 0.5 GeV", 30, 10,40);
    allPlots["lp_pt"+tag+cut+weight+runPeriod_] = new TH1F("lp_pt"+tag+cut+weight+runPeriod_,";Leading lepton P_{T} [GeV];Events / 10 GeV", 30, 0,300);
    allPlots["l2p_pt"+tag+cut+weight+runPeriod_] = new TH1F("l2p_pt"+tag+cut+weight+runPeriod_,";Sub-leading lepton P_{T} [GeV];Events / 10 GeV", 60, 0,600);
    allPlots["lp_eta"+tag+cut+weight+runPeriod_]  = new TH1F("lp_eta"+tag+cut+weight+runPeriod_,";Leading lepton #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["l2p_eta"+tag+cut+weight+runPeriod_]  = new TH1F("l2p_eta"+tag+cut+weight+runPeriod_,";Sub-Leading lepton #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["lp_phi"+tag+cut+weight+runPeriod_]  = new TH1F("lp_phi"+tag+cut+weight+runPeriod_,";Leading lepton #phi; Events", 50, -3.14,3.14);
    allPlots["l2p_phi"+tag+cut+weight+runPeriod_]  = new TH1F("l2p_phi"+tag+cut+weight+runPeriod_,";Sub-Leading lepton #phi; Events", 50, -3.14,3.14);
    allPlots["dilp_pt"+tag+cut+weight+runPeriod_] = new TH1F("dilp_pt"+tag+cut+weight+runPeriod_,";Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_m"+tag+cut+weight+runPeriod_] = new TH1F("dilp_m"+tag+cut+weight+runPeriod_,";M_{ll} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["ndilp"+tag+cut+weight+runPeriod_]     = new TH1F("ndilp"+tag+cut+weight+runPeriod_,";Di-Lepton Multiplicity;Events" ,3,0.,3.);

    //Jet plots
    allPlots["j_pt"+tag+cut+weight+runPeriod_] = new TH1F("j_pt"+tag+cut+weight+runPeriod_,";Leading Jet P_{T} [GeV];Events / 10 GeV", 40, 0,400);
    allPlots["j_pt_ch"+tag+cut+weight+runPeriod_] = new TH1F("j_pt_ch"+tag+cut+weight+runPeriod_,";PF #Sigma P_{T}^{ch} [GeV];Events / 10 GeV", 40, 0,400);
    allPlots["lj_pt"+tag+cut+weight+runPeriod_] = new TH1F("lj_pt"+tag+cut+weight+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 10 GeV", 40, 0, 400);
    allPlots["lj_pt_ch"+tag+cut+weight+runPeriod_] = new TH1F("lj_pt_ch"+tag+cut+weight+runPeriod_,";Light PF #Sigma P_{T}^{ch} [GeV];Events / 10 GeV", 40, 0,400);
    allPlots["kj_pt"+tag+cut+weight+runPeriod_] = new TH1F("kj_pt"+tag+cut+weight+runPeriod_,";Leading Kalman Jet P_{T} [GeV];Events / 10 GeV", 40, 0,400);
    allPlots["kj_pt_ch"+tag+cut+weight+runPeriod_] = new TH1F("kj_pt_ch"+tag+cut+weight+runPeriod_,";Kalman PF #Sigma P_{T}^{ch} [GeV];Events / 10 GeV", 40, 0,400);
    allPlots["j_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("j_pt_low"+tag+cut+weight+runPeriod_,";Leading Jet P_{T} [GeV];Events / 0.5 GeV", 20, 30,50);
    allPlots["lj_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("lj_pt_low"+tag+cut+weight+runPeriod_,";Leading light Jet P_{T} [GeV];Events / 0.5 GeV", 20, 30,50);
    allPlots["kj_pt_low"+tag+cut+weight+runPeriod_] = new TH1F("kj_pt_low"+tag+cut+weight+runPeriod_,";Leading Kalman Jet P_{T} [GeV];Events / 0.5 GeV", 20, 30,50);
    allPlots["nlp"+tag+cut+weight+runPeriod_]     = new TH1F("nlp"+tag+cut+weight+runPeriod_,";Lepton Multiplicity;Events" ,3,0.,3.);
    allPlots["nj"+tag+cut+weight+runPeriod_]     = new TH1F("nj"+tag+cut+weight+runPeriod_,";Jet Multiplicity (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["njtk"+tag+cut+weight+runPeriod_]     = new TH1F("njtk"+tag+cut+weight+runPeriod_,";Jet Track Multiplicity (Jet P_{T} > 30 GeV);Events" ,20,0,20.);
    allPlots["j_tk_pt"+tag+cut+weight+runPeriod_] = new TH1F("j_tk_pt"+tag+cut+weight+runPeriod_,";Jet Tracks P_{T} [GeV];Events / 10 GeV", 40, 0,400);
    allPlots["nlj"+tag+cut+weight+runPeriod_]     = new TH1F("nlj"+tag+cut+weight+runPeriod_,";Light-Jet Multiplicity (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["nkj"+tag+cut+weight+runPeriod_]     = new TH1F("nkj"+tag+cut+weight+runPeriod_,";Kalman Jet Multiplicity (#chi^{2} < 5);Events" ,4,1.,5.);
    allPlots["j_sphericity"+tag+cut+weight+runPeriod_] = new TH1F("j_sphericity"+tag+cut+weight+runPeriod_,";Jet Sphericity;Events / 0.01", 100, 0,1);
    allPlots["j_planarity"+tag+cut+weight+runPeriod_] = new TH1F("j_planarity"+tag+cut+weight+runPeriod_,";Jet Planarity;Events / 0.01", 50, 0,0.5);

    allPlots["pf_dxy"+tag+cut+weight+runPeriod_] = new TH1F("pf_dxy"+tag+cut+weight+runPeriod_,";d_{xy} [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dz"+tag+cut+weight+runPeriod_] = new TH1F("pf_dz"+tag+cut+weight+runPeriod_,";d_{z} [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dxyE"+tag+cut+weight+runPeriod_] = new TH1F("pf_dxyE"+tag+cut+weight+runPeriod_,";#sigma(d_{xy}) [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dzE"+tag+cut+weight+runPeriod_] = new TH1F("pf_dzE"+tag+cut+weight+runPeriod_,";#sigma(d_{z}) [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dxy_sig"+tag+cut+weight+runPeriod_] = new TH1F("pf_dxy_significance"+tag+cut+weight+runPeriod_,";d_{xy};Events / 1", 30, 0, 30);
    allPlots["pf_dz_sig"+tag+cut+weight+runPeriod_] = new TH1F("pf_dz_significance"+tag+cut+weight+runPeriod_,";d_{z};Events / 1", 30, 0, 30);

    allPlots["mt"+tag+cut+weight+runPeriod_]         = new TH1F("mt"+tag+cut+weight+runPeriod_,";m_{T};Events / 10 GeV" ,200,0,2000);
    
    //J/Psi plots
    if(cut.Contains("jpsi")) { // reduce number of empty histograms
    allPlots["massJPsi"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi"+tag+cut+weight+runPeriod_,";M_{#mu^{#pm}#mu^{#mp}} [GeV];Events / 20 MeV" ,30,2.8,3.4);
    allPlots["massJPsi_mu"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_mu"+tag+cut+weight+runPeriod_,";M_{J/#Psi+#mu} [GeV];Events / 10 GeV" ,30,0,300);
    allPlots["massJPsi_mu_low_dR"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_mu_low_dR"+tag+cut+weight+runPeriod_,";M_{J/#Psi+#mu} [GeV] (#DeltaR(J/#Psi,l)<2.0);Events / 10 GeV" ,30,0,300);
    allPlots["massJPsi_mu_high_dR"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_mu_high_dR"+tag+cut+weight+runPeriod_,";M_{J/#Psi+#mu} [GeV] (#DeltaR(J/#Psi,l)>2.0);Events / 10 GeV" ,30,0,300);
    allPlots["massJPsi_e"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_e"+tag+cut+weight+runPeriod_,";M_{J/#Psi+e} [GeV];Events / 10 GeV" ,30,0,300);
    allPlots["massJPsi_e_low_dR"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_e_low_dR"+tag+cut+weight+runPeriod_,";M_{J/#Psi+#e} [GeV] (#DeltaR(J/#Psi,l)<2.0);Events / 10 GeV" ,30,0,300);
    allPlots["massJPsi_e_high_dR"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_e_high_dR"+tag+cut+weight+runPeriod_,";M_{J/#Psi+#e} [GeV] (#DeltaR(J/#Psi,l)>2.0);Events / 10 GeV" ,30,0,300);
    allPlots["massJPsi_l"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_l"+tag+cut+weight+runPeriod_,";M_{J/#Psi+l} [GeV];Events / 10 GeV" ,30,0,300);
    allPlots["massJPsi_l_low_dR"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_l_low_dR"+tag+cut+weight+runPeriod_,";M_{J/#Psi+#l} [GeV] (#DeltaR(J/#Psi,l)<2.0);Events / 10 GeV" ,30,0,300);
    allPlots["massJPsi_l_high_dR"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_l_high_dR"+tag+cut+weight+runPeriod_,";M_{J/#Psi+#l} [GeV] (#DeltaR(J/#Psi,l)>2.0);Events / 10 GeV" ,30,0,300);
    allPlots["massJPsiK"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsiK"+tag+cut+weight+runPeriod_,";M_{llk} [GeV];Events / 15 MeV" ,100,4.5,6);
    allPlots["JPsi_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_pt"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi) [GeV];Events / 10 GeV", 20, 0, 200);
    allPlots["JPsi_eta"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_eta"+tag+cut+weight+runPeriod_,";J/#Psi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["JPsi_phi"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_phi"+tag+cut+weight+runPeriod_,";J/#Psi #phi; Events", 50, -3.14,3.14);
    allPlots["JPsi_mu_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu_pt"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+#mu) [GeV];Events / 10 GeV", 40, 0, 400);
    allPlots["JPsi_mu_pt_low_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu_pt_low_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+#mu) [GeV] (#DeltaR(J/#Psi,l)<2.0);Events / 10 GeV", 40, 0, 400);
    allPlots["JPsi_mu_pt_high_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu_pt_high_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+#mu) [GeV] (#DeltaR(J/#Psi,l)>2.0);Events / 10 GeV", 40, 0, 400);
    allPlots["JPsi_e_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_e_pt"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+e) [GeV];Events / 10 GeV", 40, 0, 400);
    allPlots["JPsi_e_pt_low_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_e_pt_low_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+#mu) [GeV] (#DeltaR(J/#Psi,l)<2.0);Events / 10 GeV", 40, 0, 400);
    allPlots["JPsi_e_pt_high_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_e_pt_high_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+#mu) [GeV] (#DeltaR(J/#Psi,l)>2.0);Events / 10 GeV", 40, 0, 400);
    allPlots["JPsi_l_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_l_pt"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+l) [GeV];Events / 10 GeV", 40, 0,400);
    allPlots["JPsi_l_pt_low_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_l_pt_low_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+#mu) [GeV] (#DeltaR(J/#Psi,l)<2.0);Events / 10 GeV", 40, 0, 400);
    allPlots["JPsi_l_pt_high_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_l_pt_high_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi+#mu) [GeV] (#DeltaR(J/#Psi,l)>2.0);Events / 10 GeV", 40, 0, 400);
    allPlots["JPsi_mu1_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu1_pt"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{1} P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["JPsi_mu2_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu2_pt"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{2} P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["JPsi_mu1_eta"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu1_eta"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{1} #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["JPsi_mu2_eta"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu2_eta"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{2} #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["JPsi_mu1_phi"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu1_phi"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{1} #phi; Events", 50, -3.14,3.14);
    allPlots["JPsi_mu2_phi"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu2_phi"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{2} #phi; Events", 50, -3.14,3.14);

    /*
    allPlots["JPsi_mu1_quality"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu1_quality"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{1} track quality; Events", 8, 0, 8);
    allPlots["JPsi_mu2_quality"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu2_quality"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{2} track quality; Events", 8, 0, 8);
    allPlots["JPsi_mu1_highPurity"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu1_highPurity"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{1} track highPurity; Events", 2, 0, 2);
    allPlots["JPsi_mu2_highPurity"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu2_highPurity"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{2} track highPurity; Events", 2, 0, 2);
    */
    allPlots["JPsioJet_pt"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi)/P_{T}(jet);Events / 0.025", 40, 0,1);
    allPlots["JPsioJet_pt_mu"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt_mu"+tag+cut+weight+runPeriod_,";P_{T}(#mu)/P_{T}(jet);Events / 0.025", 40, 0,1);
    allPlots["JPsioJet_pt_hard"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt_hard"+tag+cut+weight+runPeriod_,";P_{T}(hardest)/P_{T}(jet);Events / 0.025", 40, 0,1);
    allPlots["JPsioJet_pt_charged"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt_charged"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi)/P_{T}#Sigma p_{T}^{ch};Events / 0.05", 20, 0,1);
    allPlots["JPsioJet_pt_pf"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt_pf"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi)/P_{T}(jet PF tracks);Events / 0.025", 40, 0,1);

    allPlots["JPsioJet_pt_low_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt_low_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi)/P_{T}(jet) (#DeltaR(J/#Psi,l)<2.0);Events / 0.02", 10, 0,1);
    allPlots["JPsioJet_pt_high_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt_high_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi)/P_{T}(jet) (#DeltaR(J/#Psi,l)>2.0);Events / 0.02", 10, 0,1);
    allPlots["JPsioJet_pt_charged_low_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt_charged_low_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi)/P_{T}(jet charged PF tracks) (#DeltaR(J/#Psi,l)<2.0);Events / 0.02", 10, 0,1);
    allPlots["JPsioJet_pt_charged_high_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt_charged_high_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi)/P_{T}(jet charged PF tracks) (#DeltaR(J/#Psi,l)>2.0);Events / 0.02", 10, 0,1);

    allPlots["JPsioJet_pt_pf_low_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt_pf_low_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi)/P_{T}(jet all PF tracks) (#DeltaR(J/#Psi,l)<2.0);Events / 0.02", 10, 0,1);
    allPlots["JPsioJet_pt_pf_high_dR"+tag+cut+weight+runPeriod_] = new TH1F("JPsioJet_pt_pf_high_dR"+tag+cut+weight+runPeriod_,";P_{T}(J/#Psi)/P_{T}(jet all PF tracks) (#DeltaR(J/#Psi,l)>2.0);Events / 0.02", 10, 0,1);
    allPlots["dR_JPsi_mu"+tag+cut+weight+runPeriod_] = new TH1F("dR_JPsi_mu"+tag+cut+weight+runPeriod_,";#DeltaR(J/#Psi,leading #mu);Events / 0.01", 10, 0,5);
    allPlots["dR_JPsi_e"+tag+cut+weight+runPeriod_] = new TH1F("dR_JPsi_e"+tag+cut+weight+runPeriod_,";#DeltaR(J/#Psi,leading e);Events / 0.01", 10, 0,5);
    allPlots["dR_JPsi_l"+tag+cut+weight+runPeriod_] = new TH1F("dR_JPsi_l"+tag+cut+weight+runPeriod_,";#DeltaR(J/#Psi,leading l);Events / 0.01", 10, 0,5);
    allPlots["dR_JPsi_mumu"+tag+cut+weight+runPeriod_] = new TH1F("dR_JPsi_mumu"+tag+cut+weight+runPeriod_,";#DeltaR(J/#Psi_{#mu1},J/#Psi_{#mu2});Events / 0.01", 10, 0,0.5);
    allPlots["JPsi_mumu_pt_ratio"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mumu_pt_ratio"+tag+cut+weight+runPeriod_,";(#mu_{1 P_{T}}-#mu_{2 P_{T}})/(#mu_{1 P_{T}}+#mu_{2 P_{T}});Events / 0.01", 10, 0,1);
    allPlots["JPsi_l3d"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_l3d"+tag+cut+weight+runPeriod_,";J/#Psi c#tau;Events / 0.002",50, 0,0.1);
    allPlots["JPsi_sigmal3d"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_sigmal3d"+tag+cut+weight+runPeriod_,";J/#Psi #sigma_{c#tau};Events / 0.002", 50, 0,0.01);
    }
    if(cut.Contains("gjpsi")) {
      allPlots["dR_JPsi_l_good"+tag+cut+weight+runPeriod_] = new TH1F("dR_JPsi_l_good"+tag+cut+weight+runPeriod_,";#DeltaR(gen J/#Psi,leading l);Events / 0.01", 10, 0,5);
      allPlots["dR_JPsi_l_bad"+tag+cut+weight+runPeriod_] = new TH1F("dR_JPsi_l_bad"+tag+cut+weight+runPeriod_,";#DeltaR(gen J/#Psi,sub-leading l);Events / 0.01", 10, 0,5);
      allPlots["JPsi_mu1_pt_good"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu1_pt_good"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{1} P_{T} (good gen #DeltaR) [GeV];Events / 1 GeV", 50, 0,50);
      allPlots["JPsi_mu2_pt_good"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu2_pt_good"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{2} P_{T} (good gen #DeltaR(J/#Psi,closest l)) [GeV];Events / 1 GeV", 50, 0,50);
      allPlots["JPsi_mu1_pt_bad"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu1_pt_bad"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{1} P_{T} (bad gen #DeltaR(J/#Psi,closest l)) [GeV];Events / 1 GeV", 50, 0,50);
      allPlots["JPsi_mu2_pt_bad"+tag+cut+weight+runPeriod_] = new TH1F("JPsi_mu2_pt_bad"+tag+cut+weight+runPeriod_,";J/#Psi #mu_{2} P_{T} (bad gen #DeltaR(J/#Psi,closest l)) [GeV];Events / 1 GeV", 50, 0,50);
      allPlots["massJPsi_l_good_dR"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_l_good_dR"+tag+cut+weight+runPeriod_,";M_{J/#Psi+#l} [GeV] (good gen #DeltaR(J/#Psi,l));Events / 10 GeV" ,30,0,300);
      allPlots["massJPsi_l_bad_dR"+tag+cut+weight+runPeriod_]     = new TH1F("massJPsi_l_bad_dR"+tag+cut+weight+runPeriod_,";M_{J/#Psi+#l} [GeV] (bad gen #DeltaR(J/#Psi,l));Events / 10 GeV" ,30,0,300);
    }
 
    //D meson plots
    if(cut.Contains("meson")) {
    allPlots["massD0"+tag+cut+weight+runPeriod_]     = new TH1F("massD0"+tag+cut+weight+runPeriod_,";M_{K#pi};Events / 5 MeV" ,60,1.7,2.0);
    allPlots["massD0_fromDs"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_fromDs"+tag+cut+weight+runPeriod_,";M_{K#pi}M_{D^{0} from D*-D^{0} peak};Events / 5 MeV" ,60,1.7,2.0);
    allPlots["massD0_fromDsloose"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_fromDsloose"+tag+cut+weight+runPeriod_,"M_{K#pi};M_{D^{0} from D*-D^{0} peak};Events / 5 MeV" ,60,1.7,2.0);
    allPlots["massD0_notDs"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_notDs"+tag+cut+weight+runPeriod_,";M_{K#pi}M_{D^{0} not from D*-D^{0} peak};Events / 5 MeV" ,60,1.7,2.0);
    allPlots["massD0_notDsloose"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_notDsloose"+tag+cut+weight+runPeriod_,";M_{K#pi}M_{D^{0} not from D*-D^{0} peak};Events / 5 MeV" ,60,1.7,2.0);
    allPlots["massD0_lep_tag"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_lep_tag"+tag+cut+weight+runPeriod_,";M_{K#pi+l};Events / 5 MeV" ,60,1.7,2.0);
    allPlots["massD0_mu_tag"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_mu_tag"+tag+cut+weight+runPeriod_,";M_{K#pi}+#mu;Events / 10 MeV" ,30,1.7,2.0);
    /*
    allPlots["massD0_mu_tagB"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_mu_tagB"+tag+cut+weight+runPeriod_,";M_{K#pi+#mu+#mu}(B);Events / 5 MeV" ,200,5,7);
    allPlots["massD0_e_tag"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_e_tag"+tag+cut+weight+runPeriod_,";M_{K#pi+e}(D^{0})+e;Events / 5 MeV" ,60,1.7,2.0);
    allPlots["massD0_mu"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_mu"+tag+cut+weight+runPeriod_,";M_{K#pi+#mu};Events / 10 GeV" ,30,0,300);
    allPlots["massD0_e"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_e"+tag+cut+weight+runPeriod_,";M_{K#pi+e};Events / 10 GeV" ,30,0,300);
    */
    allPlots["massD0_l"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_l"+tag+cut+weight+runPeriod_,";M_{K#pi+l};Events / 10 GeV" ,30,0,300);
    allPlots["massD0_mu_tag_l"+tag+cut+weight+runPeriod_]     = new TH1F("massD0_mu_tag_l"+tag+cut+weight+runPeriod_,";M_{K#pi+l}+#mu;Events / 10 GeV" ,30,0,300);
    allPlots["D0_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_pt"+tag+cut+weight+runPeriod_,";D^{0} P_{T} [GeV];Events / 1 GeV", 150, 0,150);
    allPlots["D0_pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_pi_pt"+tag+cut+weight+runPeriod_,";D^{0} #pi P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["D0_K_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_K_pt"+tag+cut+weight+runPeriod_,";D^{0} K P_{T} [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["D0_p"+tag+cut+weight+runPeriod_] = new TH1F("D0_p"+tag+cut+weight+runPeriod_,";D^{0} P [GeV];Events / 1 GeV", 150, 0,150);
    allPlots["D0_pi_p"+tag+cut+weight+runPeriod_] = new TH1F("D0_pi_p"+tag+cut+weight+runPeriod_,";D^{0} #pi P [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["D0_K_p"+tag+cut+weight+runPeriod_] = new TH1F("D0_K_p"+tag+cut+weight+runPeriod_,";D^{0} K P [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["D0_dEtapik"+tag+cut+weight+runPeriod_] = new TH1F("D0_dEtapik"+tag+cut+weight+runPeriod_,";D^{0} #Delta#eta(pi,K) [GeV];Events / 0.05", 40, -1., 1.);
    allPlots["D0_dPhipik"+tag+cut+weight+runPeriod_] = new TH1F("D0_dPhipik"+tag+cut+weight+runPeriod_,";D^{0} #Delta#phi(pi,K) [GeV];Events / 0.05", 40, -1., 1.);
    allPlots["D0_dRpik"+tag+cut+weight+runPeriod_] = new TH1F("D0_dRpik"+tag+cut+weight+runPeriod_,";D^{0} #DeltaR(pi,K) [GeV];Events / 0.01", 100, 0, 1.);
    allPlots["D0_mu_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_pt"+tag+cut+weight+runPeriod_,";D^{0} #mu P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["D0_mu_tag_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_pt"+tag+cut+weight+runPeriod_,";D^{0}_{#mu} P_{T} [GeV];Events / 1 GeV", 150, 0,150);
    allPlots["D0_mu_tag_pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_pi_pt"+tag+cut+weight+runPeriod_,";D^{0}_{#mu} #pi P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["D0_mu_tag_K_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_K_pt"+tag+cut+weight+runPeriod_,";D^{0}_{#mu} K P_{T} [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["D0_pi_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_pi_eta"+tag+cut+weight+runPeriod_,";D^{0} #pi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["D0_K_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_K_eta"+tag+cut+weight+runPeriod_,";D^{0} K #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["D0_mu_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_eta"+tag+cut+weight+runPeriod_,";D^{0} #mu #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["D0_pi_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_pi_phi"+tag+cut+weight+runPeriod_,";D^{0} #pi #phi; Events", 50, -3.14,3.14);
    allPlots["D0_K_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_K_phi"+tag+cut+weight+runPeriod_,";D^{0} K #phi; Events", 50, -3.14,3.14);
    allPlots["D0_mu_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_phi"+tag+cut+weight+runPeriod_,";D^{0} #mu #phi; Events", 50, -3.14,3.14);

    allPlots["D0_mu_tag_pi_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_pi_eta"+tag+cut+weight+runPeriod_,";D^{0}_{#mu} #pi #eta [GeV];Events / 0.1", 30, -2.5, 2.5);
    allPlots["D0_mu_tag_K_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_K_eta"+tag+cut+weight+runPeriod_,";D^{0}_{#mu} Kpi #eta [GeV];Events / 0.1", 30, -2.5, 2.5);
    allPlots["D0_mu_tag_pi_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_pi_phi"+tag+cut+weight+runPeriod_,";D^{0}_{#mu} #pi #phi [GeV];Events", 50, -3.14,3.14);
    allPlots["D0_mu_tag_K_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_K_phi"+tag+cut+weight+runPeriod_,";D^{0}_{#mu} K #phi [GeV];Events", 50, -3.14,3.14);
    /*
    allPlots["D0_pi_quality"+tag+cut+weight+runPeriod_] = new TH1F("D0_pi_quality"+tag+cut+weight+runPeriod_,";D^{0} #pi track quality; Events", 8, 0, 8);
    allPlots["D0_K_quality"+tag+cut+weight+runPeriod_] = new TH1F("D0_K_quality"+tag+cut+weight+runPeriod_,";D^{0} K track quality; Events", 8, 0, 8);
    allPlots["D0_pi_highPurity"+tag+cut+weight+runPeriod_] = new TH1F("D0_pi_highPurity"+tag+cut+weight+runPeriod_,";D^{0} #pi track highPurity; Events", 2, 0, 2);
    allPlots["D0_K_highPurity"+tag+cut+weight+runPeriod_] = new TH1F("D0_K_highPurity"+tag+cut+weight+runPeriod_,";D^{0} K track highPurity; Events", 2, 0, 2);
    allPlots["D0_mu_highPurity"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_highPurity"+tag+cut+weight+runPeriod_,";D^{0} #mu track highPurity; Events", 2, 0, 2);
    */
    //loose Ds
    allPlots["D0_fromDsloose_pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDsloose_pi_pt"+tag+cut+weight+runPeriod_,";D^{0} #pi P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["D0_fromDsloose_K_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDsloose_K_pt"+tag+cut+weight+runPeriod_,";D^{0} K P_{T} [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["D0_fromDsloose_pi_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDsloose_pi_eta"+tag+cut+weight+runPeriod_,";D^{0} #pi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["D0_fromDsloose_K_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDsloose_K_eta"+tag+cut+weight+runPeriod_,";D^{0} K #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["D0_fromDsloose_pi_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDsloose_pi_phi"+tag+cut+weight+runPeriod_,";D^{0} #pi #phi; Events", 50, -3.14,3.14);
    allPlots["D0_fromDsloose_K_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDsloose_K_phi"+tag+cut+weight+runPeriod_,";D^{0} K #phi; Events", 50, -3.14,3.14);
    allPlots["Ds_fromDsloose_pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDsloose_pi_pt"+tag+cut+weight+runPeriod_,";D* #pi P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["Ds_fromDsloose_K_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDsloose_K_pt"+tag+cut+weight+runPeriod_,";D* K P_{T} [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["Ds_fromDsloose_pi_eta"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDsloose_pi_eta"+tag+cut+weight+runPeriod_,";D* #pi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["Ds_fromDsloose_K_eta"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDsloose_K_eta"+tag+cut+weight+runPeriod_,";D* K #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["Ds_fromDsloose_pi_phi"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDsloose_pi_phi"+tag+cut+weight+runPeriod_,";D* #pi #phi; Events", 50, -3.14,3.14);
    allPlots["Ds_fromDsloose_K_phi"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDsloose_K_phi"+tag+cut+weight+runPeriod_,";D* K #phi; Events", 50, -3.14,3.14);
    //tight D^0
    allPlots["D0_fromDs_pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDs_pi_pt"+tag+cut+weight+runPeriod_,";D^{0} #pi P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["D0_fromDs_K_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDs_K_pt"+tag+cut+weight+runPeriod_,";D^{0} K P_{T} [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["D0_fromDs_pi_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDs_pi_eta"+tag+cut+weight+runPeriod_,";D^{0} #pi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["D0_fromDs_K_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDs_K_eta"+tag+cut+weight+runPeriod_,";D^{0} K #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["D0_fromDs_pi_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDs_pi_phi"+tag+cut+weight+runPeriod_,";D^{0} #pi #phi; Events", 50, -3.14,3.14);
    allPlots["D0_fromDs_K_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_fromDs_K_phi"+tag+cut+weight+runPeriod_,";D^{0} K #phi; Events", 50, -3.14,3.14);
    //tight D*
    allPlots["Ds_fromDs_pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDs_pi_pt"+tag+cut+weight+runPeriod_,";D* #pi P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["Ds_fromDs_pi2_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDs_pi2_pt"+tag+cut+weight+runPeriod_,";D* #pi_{2} P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["Ds_fromDs_K_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDs_K_pt"+tag+cut+weight+runPeriod_,";D* K P_{T} [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["Ds_fromDs_pi_eta"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDs_pi_eta"+tag+cut+weight+runPeriod_,";D* #pi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["Ds_fromDs_K_eta"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDs_K_eta"+tag+cut+weight+runPeriod_,";D* K #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["Ds_fromDs_pi_phi"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDs_pi_phi"+tag+cut+weight+runPeriod_,";D* #pi #phi; Events", 50, -3.14,3.14);
    allPlots["Ds_fromDs_K_phi"+tag+cut+weight+runPeriod_] = new TH1F("Ds_fromDs_K_phi"+tag+cut+weight+runPeriod_,";D* K #phi; Events", 50, -3.14,3.14);

    allPlots["D0oJet_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0oJet_pt"+tag+cut+weight+runPeriod_,";P_{T}(D^{0})/P_{T}(jet);Events / 0.05", 20, 0,1);
    allPlots["D0oJet_pt_mu"+tag+cut+weight+runPeriod_] = new TH1F("D0oJet_pt_mu"+tag+cut+weight+runPeriod_,";P_{T}(#mu)/P_{T}(jet);Events / 0.02", 10, 0,1);
    allPlots["D0oJet_pt_hard"+tag+cut+weight+runPeriod_] = new TH1F("D0oJet_pt_hard"+tag+cut+weight+runPeriod_,";P_{T}(hardest)/P_{T}(jet);Events / 0.02", 10, 0,1);
    allPlots["D0oJet_pt_charged"+tag+cut+weight+runPeriod_] = new TH1F("D0oJet_pt_charged"+tag+cut+weight+runPeriod_,";P_{T}(D^{0})/#Sigma p_{T}^{ch};Events / 0.05", 20, 0,1);
    allPlots["D0oJet_pt_pf"+tag+cut+weight+runPeriod_] = new TH1F("D0oJet_pt_pf"+tag+cut+weight+runPeriod_,";P_{T}(D^{0})/P_{T}(jet PF tracks);Events / 0.05", 20, 0,1);

    allPlots["D0_mu_tag_oJet_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_oJet_pt"+tag+cut+weight+runPeriod_,";P_{T}(D^{0}_{#mu})/P_{T}(jet);Events / 0.05", 20, 0,1);
    allPlots["D0_mu_tag_oJet_pt_mu"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_oJet_pt_mu"+tag+cut+weight+runPeriod_,";P_{T}(#mu)/P_{T}(jet);Events / 0.02", 10, 0,1);
    allPlots["D0_mu_tag_oJet_pt_hard"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_oJet_pt_hard"+tag+cut+weight+runPeriod_,";P_{T}(hardest)/P_{T}(jet);Events / 0.02", 10, 0,1);
    allPlots["D0_mu_tag_oJet_pt_charged"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_oJet_pt_charged"+tag+cut+weight+runPeriod_,";P_{T}(D^{0}_{#mu})/#Sigma p_{T}^{ch};Events / 0.05", 20, 0,1);
    allPlots["D0_mu_tag_mu_oJet_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_mu_oJet_pt"+tag+cut+weight+runPeriod_,";P_{T}(D^{0}_{#mu}+#mu)/P_{T} (jet);Events / 0.05", 20, 0,1);
    allPlots["D0_mu_tag_mu_oJet_pt_charged"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_mu_oJet_pt_charged"+tag+cut+weight+runPeriod_,";P_{T}(D^{0}_{#mu}+#mu)/#Sigma p_{T}^{ch};Events / 0.05", 20, 0,1);
    allPlots["D0_mu_tag_oJet_pt_pf"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_oJet_pt_pf"+tag+cut+weight+runPeriod_,";P_{T}(D^{0}_{#mu})/P_{T}(jet PF tracks);Events / 0.05", 20, 0,1);

    allPlots["D0dotJet"+tag+cut+weight+runPeriod_] = new TH1F("D0dotJet"+tag+cut+weight+runPeriod_,";P(D^{0})#upointP(jet)/|P(jet)|;Events / 1", 200, 0, 200);
    allPlots["D0cosJet"+tag+cut+weight+runPeriod_] = new TH1F("D0cosJet"+tag+cut+weight+runPeriod_,";cos(P(D^{0}),P(jet));Events / 0.0002", 200, 0.98, 1.0);
    allPlots["D0oJet_pt_fromDs"+tag+cut+weight+runPeriod_]     = new TH1F("D0oJet_pt_fromDs"+tag+cut+weight+runPeriod_,";P_{T}(D^{0})/P_{T}(jet) from D*-D^{0} peak};Events / 0.02", 10, 0, 1);
    allPlots["D0oJet_pt_mu_fromDs"+tag+cut+weight+runPeriod_]     = new TH1F("D0oJet_pt_mu_fromDs"+tag+cut+weight+runPeriod_,";P_{T}(#mu)/P_{T}(jet) from D*-D^{0} peak};Events / 0.02", 10, 0, 1);
    allPlots["D0oJet_pt_hard_fromDs"+tag+cut+weight+runPeriod_]     = new TH1F("D0oJet_pt_hard_fromDs"+tag+cut+weight+runPeriod_,";P_{T}(hardest)/P_{T}(jet) from D*-D^{0} peak};Events / 0.02", 10, 0, 1);
    allPlots["D0oJet_pt_charged_fromDs"+tag+cut+weight+runPeriod_]     = new TH1F("D0oJet_pt_charged_fromDs"+tag+cut+weight+runPeriod_,";P_{T}(D^{0})/P_{T}(jet charged PF tracks) from D*-D^{0} peak};Events / 0.02", 10, 0, 1);
    allPlots["D0oJet_pt_pf_fromDs"+tag+cut+weight+runPeriod_]     = new TH1F("D0oJet_pt_pf_fromDs"+tag+cut+weight+runPeriod_,";P_{T}(#pf)/P_{T}(jet PF tracks) from D*-D^{0} peak};Events / 0.02", 10, 0, 1);
    allPlots["D0_l3d"+tag+cut+weight+runPeriod_] = new TH1F("D0_l3d"+tag+cut+weight+runPeriod_,";D^{0} c#tau;Events / 0.002",100, -0.1,0.1);
    allPlots["D0_sigmal3d"+tag+cut+weight+runPeriod_] = new TH1F("D0_sigmal3d"+tag+cut+weight+runPeriod_,";D^{0} #sigma_{c#tau};Events / 0.0001", 50, 0,0.01);
    allPlots["D0_mu_tag_l3d"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_l3d"+tag+cut+weight+runPeriod_,";D^{0}_{#mu} c#tau;Events / 0.002",100, -0.1,0.1);
    allPlots["D0_mu_tag_sigmal3d"+tag+cut+weight+runPeriod_] = new TH1F("D0_mu_tag_sigmal3d"+tag+cut+weight+runPeriod_,";D^{0}_{#mu} #sigma_{c#tau};Events / 0.002", 50, 0,0.01);
    if(cut.Contains("gmeson")) {
    allPlots["D0_l3d_good"+tag+cut+weight+runPeriod_] = new TH1F("D0_l3d_good"+tag+cut+weight+runPeriod_,";D^{0} c#tau;Events / 0.001",200, -0.1,0.1);
    allPlots["D0_sigmal3d_good"+tag+cut+weight+runPeriod_] = new TH1F("D0_sigmal3d_good"+tag+cut+weight+runPeriod_,";D^{0} #sigma_{c#tau};Events / 0.0001", 50, 0,0.01);
    allPlots["D0_l3d_bad"+tag+cut+weight+runPeriod_] = new TH1F("D0_l3d_bad"+tag+cut+weight+runPeriod_,";D^{0} c#tau;Events / 0.001",200, -0.1,0.1);
    allPlots["D0_sigmal3d_bad"+tag+cut+weight+runPeriod_] = new TH1F("D0_sigmal3d_bad"+tag+cut+weight+runPeriod_,";D^{0} #sigma_{c#tau};Events / 0.0001", 50, 0,0.01);
    }

    //D*
    allPlots["massDsmD0loose"+tag+cut+weight+runPeriod_]     = new TH1F("massDsmD0loose"+tag+cut+weight+runPeriod_,";M_{K#pi#pi} - M_{K#pi};Events / 0.25 MeV" ,40,0.14,0.16);
    allPlots["massDsmD0"+tag+cut+weight+runPeriod_]     = new TH1F("massDsmD0"+tag+cut+weight+runPeriod_,";M_{K#pi#pi} - M_{K#pi};Events / 0.25 MeV" ,40,0.14,0.16);
    allPlots["massDs"+tag+cut+weight+runPeriod_]     = new TH1F("massDs"+tag+cut+weight+runPeriod_,";M_{D*};Events / 6 MeV" ,200,1.6,2.2);

    allPlots["Ds_pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_pi_pt"+tag+cut+weight+runPeriod_,";D* #pi P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["Ds_pi2_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_pi2_pt"+tag+cut+weight+runPeriod_,";D* #pi_{2} P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["Ds_K_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_K_pt"+tag+cut+weight+runPeriod_,";D* K P_{T} [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["Ds_pi_eta"+tag+cut+weight+runPeriod_] = new TH1F("Ds_pi_eta"+tag+cut+weight+runPeriod_,";D* #pi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["Ds_K_eta"+tag+cut+weight+runPeriod_] = new TH1F("Ds_K_eta"+tag+cut+weight+runPeriod_,";D* K #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["Ds_pi_phi"+tag+cut+weight+runPeriod_] = new TH1F("Ds_pi_phi"+tag+cut+weight+runPeriod_,";D* #pi #phi; Events", 50, -3.14,3.14);
    allPlots["Ds_K_phi"+tag+cut+weight+runPeriod_] = new TH1F("Ds_K_phi"+tag+cut+weight+runPeriod_,";D* K #phi; Events", 50, -3.14,3.14);
    /*
    allPlots["Ds_pi_quality"+tag+cut+weight+runPeriod_] = new TH1F("Ds_pi_quality"+tag+cut+weight+runPeriod_,";D* #pi track quality; Events", 8, 0, 8);
    allPlots["Ds_K_quality"+tag+cut+weight+runPeriod_] = new TH1F("Ds_K_quality"+tag+cut+weight+runPeriod_,";D* K track quality; Events", 8, 0, 8);
    allPlots["Ds_pi_highPurity"+tag+cut+weight+runPeriod_] = new TH1F("Ds_pi_highPurity"+tag+cut+weight+runPeriod_,";D* #pi track highPurity; Events", 2, 0, 2);
    allPlots["Ds_K_highPurity"+tag+cut+weight+runPeriod_] = new TH1F("Ds_K_highPurity"+tag+cut+weight+runPeriod_,";D* K track highPurity; Events", 2, 0, 2);
    */

    allPlots["massDs_l"+tag+cut+weight+runPeriod_]     = new TH1F("massDs_l"+tag+cut+weight+runPeriod_,";M_{K#pi+l};Events / 10 GeV" ,30,0,300);
    allPlots["massDs_mu"+tag+cut+weight+runPeriod_]     = new TH1F("massDs_mu"+tag+cut+weight+runPeriod_,";M_{K#pi+#mu};Events / 10 GeV" ,30,0,300);
    allPlots["massDs_e"+tag+cut+weight+runPeriod_]     = new TH1F("massDs_ele"+tag+cut+weight+runPeriod_,";M_{K#pi+e};Events / 10 GeV" ,30,0,300);
    allPlots["massDs_fromDs"+tag+cut+weight+runPeriod_] = new TH1F("massDs_fromDs"+tag+cut+weight+runPeriod_,";M_{D* from D*-D^{0} peak};Events / 5 MeV" ,40,1.9,2.1);
    allPlots["massDs_fromDsloose"+tag+cut+weight+runPeriod_] = new TH1F("massDs_fromDsloose"+tag+cut+weight+runPeriod_,";M_{D* from D*-D^{0} peak};Events / 5 MeV" ,40,1.9,2.1);
    allPlots["massDs_notDs"+tag+cut+weight+runPeriod_]  = new TH1F("massDs_notDs"+tag+cut+weight+runPeriod_,";M_{D* not from D*-D^{0} peak};Events / 5 MeV" ,40,1.9,2.1);
    allPlots["massDs_notDsloose"+tag+cut+weight+runPeriod_]  = new TH1F("massDs_notDsloose"+tag+cut+weight+runPeriod_,";M_{D* not from D*-D^{0} peak};Events / 5 MeV" ,40,1.9,2.1);
    allPlots["D0_notDsloose_pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDsloose_pi_pt"+tag+cut+weight+runPeriod_,";D^{0} #pi P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["D0_notDsloose_K_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDsloose_K_pt"+tag+cut+weight+runPeriod_,";D^{0} K P_{T} [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["D0_notDsloose_pi_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDsloose_pi_eta"+tag+cut+weight+runPeriod_,";D^{0} #pi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["D0_notDsloose_K_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDsloose_K_eta"+tag+cut+weight+runPeriod_,";D^{0} K #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["D0_notDsloose_pi_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDsloose_pi_phi"+tag+cut+weight+runPeriod_,";D^{0} #pi #phi; Events", 50, -3.14,3.14);
    allPlots["D0_notDsloose_K_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDsloose_K_phi"+tag+cut+weight+runPeriod_,";D^{0} K #phi; Events", 50, -3.14,3.14);
    allPlots["D0_notDs_pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDs_pi_pt"+tag+cut+weight+runPeriod_,";D^{0} #pi P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["D0_notDs_pi2_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDs_pi2_pt"+tag+cut+weight+runPeriod_,";D^{0} #pi_{2} P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["D0_notDs_K_pt"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDs_K_pt"+tag+cut+weight+runPeriod_,";D^{0} K P_{T} [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["D0_notDs_pi_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDs_pi_eta"+tag+cut+weight+runPeriod_,";D^{0} #pi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["D0_notDs_K_eta"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDs_K_eta"+tag+cut+weight+runPeriod_,";D^{0} K #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["D0_notDs_pi_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDs_pi_phi"+tag+cut+weight+runPeriod_,";D^{0} #pi #phi; Events", 50, -3.14,3.14);
    allPlots["D0_notDs_K_phi"+tag+cut+weight+runPeriod_] = new TH1F("D0_notDs_K_phi"+tag+cut+weight+runPeriod_,";D^{0} K #phi; Events", 50, -3.14,3.14);

    allPlots["Ds_notDsloose_pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDsloose_pi_pt"+tag+cut+weight+runPeriod_,";D* #pi P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["Ds_notDsloose_K_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDsloose_K_pt"+tag+cut+weight+runPeriod_,";D* K P_{T} [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["Ds_notDsloose_pi_eta"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDsloose_pi_eta"+tag+cut+weight+runPeriod_,";D* #pi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["Ds_notDsloose_K_eta"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDsloose_K_eta"+tag+cut+weight+runPeriod_,";D* K #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["Ds_notDsloose_pi_phi"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDsloose_pi_phi"+tag+cut+weight+runPeriod_,";D* #pi #phi; Events", 50, -3.14,3.14);
    allPlots["Ds_notDsloose_K_phi"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDsloose_K_phi"+tag+cut+weight+runPeriod_,";D* K #phi; Events", 50, -3.14,3.14);
    allPlots["Ds_notDs_pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDs_pi_pt"+tag+cut+weight+runPeriod_,";D* #pi P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["Ds_notDs_pi2_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDs_pi2_pt"+tag+cut+weight+runPeriod_,";D* #pi_{2} P_{T} [GeV];Events / 1 GeV", 50, 0,50);
    allPlots["Ds_notDs_K_pt"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDs_K_pt"+tag+cut+weight+runPeriod_,";D* K P_{T} [GeV];Events / 1 GeV", 25, 0,25);
    allPlots["Ds_notDs_pi_eta"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDs_pi_eta"+tag+cut+weight+runPeriod_,";D* #pi #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["Ds_notDs_K_eta"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDs_K_eta"+tag+cut+weight+runPeriod_,";D* K #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["Ds_notDs_pi_phi"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDs_pi_phi"+tag+cut+weight+runPeriod_,";D* #pi #phi; Events", 50, -3.14,3.14);
    allPlots["Ds_notDs_K_phi"+tag+cut+weight+runPeriod_] = new TH1F("Ds_notDs_K_phi"+tag+cut+weight+runPeriod_,";D* K #phi; Events", 50, -3.14,3.14);

    allPlots["DsoJet_pt"+tag+cut+weight+runPeriod_] = new TH1F("DsoJet_pt"+tag+cut+weight+runPeriod_,";P_{T}(D*)/P_{T}(jet);Events / 0.05", 20, 0,1);
    allPlots["DsoJet_pt_mu"+tag+cut+weight+runPeriod_] = new TH1F("DsoJet_pt_mu"+tag+cut+weight+runPeriod_,";P_{T}(#mu)/P_{T}(jet);Events / 0.02", 10, 0,1);
    allPlots["DsoJet_pt_hard"+tag+cut+weight+runPeriod_] = new TH1F("DsoJet_pt_hard"+tag+cut+weight+runPeriod_,";P_{T}(hardest)/P_{T}(jet);Events / 0.02", 10, 0,1);
    allPlots["DsoJet_pt_charged"+tag+cut+weight+runPeriod_] = new TH1F("DsoJet_pt_charged"+tag+cut+weight+runPeriod_,";P_{T}(D*)/P_{T}(jet charged PF tracks);Events / 0.05", 20, 0,1);
    allPlots["DsoJet_pt_pf"+tag+cut+weight+runPeriod_] = new TH1F("DsoJet_pt_pf"+tag+cut+weight+runPeriod_,";P_{T}(D*)/P_{T}(jet PF tracks);Events / 0.05", 20, 0,1);
    }

    allPlots["pi_pt"+tag+cut+weight+runPeriod_] = new TH1F("pi_pt"+tag+cut+weight+runPeriod_,";#pi^{#pm} P_{T} [GeV];Events / 5 GeV", 10, 0,50);

    allPlots["MET"+tag+cut+weight+runPeriod_] = new TH1F("MET"+tag+cut+weight+runPeriod_,";E^{miss}_{T} [GeV];Events / 20 GeV", 11,0,220);
    allPlots["HT"+tag+cut+weight+runPeriod_] = new TH1F("HT"+tag+cut+weight+runPeriod_,";H_{T} [GeV];Events / 20 GeV", 55,0,1100);
    allPlots["H"+tag+cut+weight+runPeriod_] = new TH1F("H"+tag+cut+weight+runPeriod_,";H [GeV];Events / 20 GeV", 55,0,1100);
    allPlots["ST"+tag+cut+weight+runPeriod_] = new TH1F("ST"+tag+cut+weight+runPeriod_,";S_{T} [GeV];Events / 20 GeV", 55,0,1100);
    allPlots["MET2oST"+tag+cut+weight+runPeriod_] = new TH1F("MET2oST"+tag+cut+weight+runPeriod_,";MET^{2}/ST [GeV];Events / 20 GeV", 10,0,200);
    allPlots["charge"+tag+cut+weight+runPeriod_] = new TH1F("charge"+tag+cut+weight+runPeriod_,";Charge(l_{1}*l_{2});Events", 5,-2,2);
    allPlots["csv"+tag+cut+weight+runPeriod_] = new TH1F("CSV"+tag+cut+weight+runPeriod_,";Jet CSV;Events / 0.1", 10,0,1);
    allPlots["dR"+tag+cut+weight+runPeriod_] = new TH1F("dR"+tag+cut+weight+runPeriod_,";dR;Events / 0.05", 20,0.0,1.);
    allPlots["pflp_pt"+tag+cut+weight+runPeriod_] = new TH1F("pflp_pt"+tag+cut+weight+runPeriod_,";PF lepton P_{T} [GeV];Events / 0.2 GeV", 15, 0,3);

    //gen-level plots
    allPlots["gtop_pt"+tag+cut+weight+runPeriod_] = new TH1F("gtop_pt"+tag+cut+weight+runPeriod_,";Generator top P_{T} [GeV];Events / 10 GeV", 40, 0,400);

    //Z control plots
    allPlots["massZ"+tag+cut+weight+runPeriod_]     = new TH1F("massZ_control"+tag+cut+weight+runPeriod_,";M_{l^{#pm}l^{#mp}};Events / 1.0 GeV" ,30,81,111);
    //allPlots["chargeZ"+tag+cut+weight+runPeriod_]     = new TH1F("chargeZ_control"+tag+cut+weight+runPeriod_,";Charage (l^#pm) * Charge(l^#mp);Events / 1.0 GeV" ,5,-2,2);

    //Event plots
    allPlots["nevt"+tag+cut+weight+runPeriod_]     = new TH1F("nevt"+tag+cut+weight+runPeriod_,";Event Multiplicity;Events" ,1,1.,2.);
    allPlots["weight"+tag+cut+weight+runPeriod_]     = new TH1F("weight"+tag+cut+weight+runPeriod_,";weights;Events/ 1.0" ,20,0.,2.);
    allPlots["norm"+tag+cut+weight+runPeriod_]     = new TH1F("norm"+tag+cut+weight+runPeriod_,";norm;Events / 1.0" ,2,0.,2.);
    allPlots["relIso"+tag+cut+weight+runPeriod_] = new TH1F("relIso"+tag+cut+weight+runPeriod_,";relIso;Events / 0.01", 25,0,0.25);
    allPlots["nvtx"+tag+cut+weight+runPeriod_]     = new TH1F("nvtx"+tag+cut+weight+runPeriod_,";Vertex Multiplicity;Events / 1.0" ,50,0.,50.);
  }
  }
  }
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }

}

StdPlots::~StdPlots() {
  //for (auto it : allPlots) delete it.second;
}

void StdPlots::SetNorm(float norm) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting norm= " << norm << std::endl;
  norm_ = norm;
}

void StdPlots::SetSFs(float sfs) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting SFs= " << sfs << std::endl;
  sfs_.first = sfs;
  //sfs_.second = 0.;
  sfs_.second.first = 0.;
  sfs_.second.second = 0.;
}

void StdPlots::SetSFs(float sfs, float unc) {
  SetSFs(sfs, unc/2, unc/2);
}

void StdPlots::SetSFs(float sfs, float unc_u, float unc_d) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting SFs= " << sfs << std::endl;
  sfs_.first = sfs;
  //sfs_.second = unc_u;
  sfs_.second.first = unc_u;
  sfs_.second.second = unc_d;
}

void StdPlots::SetPuWgt(float puWgt) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting puWgt= " << puWgt << std::endl;
  puWgt_ = puWgt;
  allPlots["puwgtctr"+runPeriod_]->Fill(0.,1.0);
  allPlots["puwgtctr"+runPeriod_]->Fill(1,puWgt_);
}

void StdPlots::SetTopPtWgt(float top_pt_wgt) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting top pT weight= " << top_pt_wgt << std::endl;
  allPlots["topptwgt"+runPeriod_]->Fill(0.,1.0);
  allPlots["topptwgt"+runPeriod_]->Fill(1.,top_pt_wgt);
  top_pt_wgt_ = top_pt_wgt;
  /*
  top_pt_wgt_vec.push_back(top_pt_wgt);
  */
}

void StdPlots::SetTrackerWgt(float tracker_wgt) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting SFs= " << tracker_wgt << std::endl;
  tracker_wgt_ = tracker_wgt;
}

void StdPlots::SetPiWgt(float pi_wgt, float unc) {
  SetPiWgt(pi_wgt, unc/2, unc/2);
}

void StdPlots::SetPiWgt(float pi_wgt, float unc_u, float unc_d) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting SFs= " << pi_wgt << std::endl;
  pi_wgt_.first = pi_wgt;
  //pi_wgt_.second = unc_u;
  pi_wgt_.second.first = unc_u;
  pi_wgt_.second.second = unc_d;
}

void StdPlots::SetRbWgt(float rbWgt, CharmMeson meson) {
  if(!isGood_) return;
  if(debug_) std::cout << "Setting rbWgt= " << rbWgt << std::endl;
  rbWgt_ = rbWgt;
  switch(meson) {
    case JPsi: allPlots["rbwgt_jpsi"+runPeriod_]->Fill(0.,1.0);
               allPlots["rbwgt_jpsi"+runPeriod_]->Fill(1.,rbWgt);
               break;
    case D0:   allPlots["rbwgt_d0"+runPeriod_]->Fill(0.,1.0);
               allPlots["rbwgt_d0"+runPeriod_]->Fill(1.,rbWgt);
               break;
    case D0mu: allPlots["rbwgt_d0mu"+runPeriod_]->Fill(0.,1.0);
               allPlots["rbwgt_d0mu"+runPeriod_]->Fill(1.,rbWgt);
               break;
  }
}

void StdPlots::CheckRunPeriod(TString runPeriod) {
  if(name_.Contains("Data")) return;
  if(runPeriod_.Contains(runPeriod))
    isGood_ = true;
  else
    isGood_ = false;
}

void StdPlots::Fill(double nevt, double nvtx, double HT, double ST, double MET, TString chTag, TString name) {
  if(!isGood_) return;
  float wgt = norm_ * sfs_.first * puWgt_ * top_pt_wgt_;
  if(!name.EqualTo("")) name = "_" + name;
  if(debug_) std::cout << "Filling nvtx" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;

  allPlots["nevt"+chTag+name+runPeriod_]->Fill(nevt,getWgt());
  //allPlots["weight"+chTag+name+runPeriod_]->Fill(wgt,norm_);
  //allPlots["norm"+chTag+name+runPeriod_]->Fill(norm_,norm_);
  allPlots["nvtx"+chTag+name+runPeriod_]->Fill(nvtx,getWgt());

  //allPlots["HT"+chTag+name+runPeriod_]->Fill(HT,getWgt());
  allPlots["ST"+chTag+name+runPeriod_]->Fill(ST,getWgt());
  allPlots["MET2oST"+chTag+name+runPeriod_]->Fill(pow(MET,2)/ST,getWgt());
  allPlots["MET"+chTag+name+runPeriod_]->Fill(MET,getWgt());

  allPlots["nevt_all"+name+runPeriod_]->Fill(nevt,getWgt());
  allPlots["nvtx_all"+name+runPeriod_]->Fill(nvtx,getWgt());

  //allPlots["HT_all"+name+runPeriod_]->Fill(HT,getWgt());
  allPlots["ST_all"+name+runPeriod_]->Fill(ST,getWgt());
  allPlots["MET2oST_all"+name+runPeriod_]->Fill(pow(MET,2)/ST,getWgt());
  allPlots["MET_all"+name+runPeriod_]->Fill(MET,getWgt());
}

void StdPlots::FillGen(std::vector<Particle> tops, TString chTag, TString name) {
  if(!isGood_) return;
  float wgt = norm_ * sfs_.first * puWgt_; //No top_pt_wgt_
  if(!name.EqualTo("")) name = "_" + name;
  if(debug_) std::cout << "Filling gen-level top" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;
  for(auto& it : tops) {
    allPlots["gtop_pt"+chTag+name+runPeriod_]->Fill(it.Pt(),wgt);
    allPlots["gtop_pt_all"+name+runPeriod_]->Fill(it.Pt(),wgt);
  }
}

void StdPlots::Fill(Leptons leptons, TString chTag, TString name) {
  if(!isGood_) return;
  if(debug_) std::cout << "Filling leptons only" << std::endl;
  float wgt = norm_ * sfs_.first * puWgt_ * top_pt_wgt_;
  if(!name.EqualTo("")) name = "_" + name;
  if(debug_) std::cout << "Filling lep" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;

  allPlots["nlp"+chTag+name+runPeriod_]->Fill(leptons.size(),getWgt());
  //allPlots["nlp"+chTag+name+"_no_weight"+"_"+runPeiod_]->Fill(leptons.size(),norm_);
  if(leptons.size() > 0) {
    allPlots["lp_pt_low"+chTag+name+runPeriod_]->Fill(leptons[0].Pt(),getWgt());
    allPlots["lp_pt"+chTag+name+runPeriod_]->Fill(leptons[0].Pt(),getWgt());
    allPlots["lp_pt_all"+name+runPeriod_]->Fill(leptons[0].Pt(),getWgt());
    //allPlots["lp_pt"+chTag+name+"_"+"_no_weight"+runPeiod_]->Fill(leptons[0].Pt(),norm_);
    allPlots["relIso"+chTag+name+runPeriod_]->Fill(leptons[0].getRelIso(),getWgt());
    allPlots["lp_eta"+chTag+name+runPeriod_]->Fill(leptons[0].Eta(),getWgt());
    //allPlots["lp_eta"+chTag+name+"_"+"_no_weight"+"_"+runPeiod_]->Fill(leptons[0].Eta(),norm_);
    allPlots["lp_phi"+chTag+name+runPeriod_]->Fill(leptons[0].Phi(),getWgt());
    //allPlots["lp_phi"+chTag+name+"_"+"_no_weight"+"_"+runPeiod_]->Fill(leptons[0].Phi(),norm_);
    if(leptons.size()>1) {
      bool isZ(false);
      TLorentzVector dilp4 = leptons[0].getVec()+leptons[1].getVec();
      if(leptons[0].getPdgId() == -leptons[1].getPdgId() &&
         fabs(dilp4.M()-91)<15) { isZ=true; }
      if(isZ) {
        allPlots["massZ"+chTag+name+runPeriod_]->Fill(dilp4.M(),getWgt());
      }
      else {
        allPlots["dilp_m"+chTag+name+runPeriod_]->Fill(dilp4.M(),getWgt());
        allPlots["dilp_m_all"+name+runPeriod_]->Fill(dilp4.M(),getWgt());
        allPlots["dilp_pt"+chTag+name+runPeriod_]->Fill(dilp4.Pt(),getWgt());
        allPlots["dilp_pt_all"+name+runPeriod_]->Fill(dilp4.Pt(),getWgt());
      }
    }
  }
}

void StdPlots::Fill(Leptons &leptons, std::vector<Jet> &lightJetsVec, std::vector<Jet> &kJetsVec, std::vector<Jet> &allJetsVec, TString chTag, TString name) {
  if(!isGood_) return;
  Fill(lightJetsVec, kJetsVec, allJetsVec, chTag, name);
  float wgt = norm_ * sfs_.first * puWgt_ * top_pt_wgt_;
  if(!name.EqualTo("")) name = "_" + name;
  if(debug_) std::cout << "Filling jet" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;

  TLorentzVector pmt;
  for(auto &it : allJetsVec)
    pmt += it.getVec();
  for(size_t il = 0; il < leptons.size(); il++)
    pmt += leptons[il].getVec();
  allPlots["mt"+chTag+name+runPeriod_]->Fill(pmt.Mt(),wgt);
  allPlots["mt_all"+name+runPeriod_]->Fill(pmt.Mt(),wgt);
}

void StdPlots::Fill(std::vector<Jet> &lightJetsVec, std::vector<Jet> &kJetsVec, std::vector<Jet> &allJetsVec, TString chTag, TString name) {
  if(!isGood_) return;
  float wgt = norm_ * sfs_.first * puWgt_ * top_pt_wgt_;
  if(!name.EqualTo("")) name = "_" + name;
  if(debug_) std::cout << "Filling jet" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;

  allPlots["nj"+chTag+name+runPeriod_]->Fill(allJetsVec.size(),getWgt());
  allPlots["nj_all"+name+runPeriod_]->Fill(allJetsVec.size(),getWgt());
  allPlots["nlj"+chTag+name+runPeriod_]->Fill(lightJetsVec.size(),getWgt());
  allPlots["nlj_all"+name+runPeriod_]->Fill(lightJetsVec.size(),getWgt());
  allPlots["nkj"+chTag+name+runPeriod_]->Fill(kJetsVec.size(),getWgt());
  allPlots["nkj_all"+name+runPeriod_]->Fill(kJetsVec.size(),getWgt());
  //allPlots["nlj"+chTag+name+"_no_weight"+runPeriod_]->Fill(lightJetsVec.size(),norm_);
  //allPlots["nkj"+chTag+name+"_no_weight"+runPeriod_]->Fill(kJetsVec.size(),norm_);
  //allPlots["nlj_all"+name+runPeriod_]->Fill(lightJetsVec.size(),getWgt());
  //allPlots["nkj_all"+name+runPeriod_]->Fill(kJetsVec.size(),getWgt());

  Float_t htsum(0),hsum(0);
  for(auto & it : allJetsVec) {
    htsum += it.Pt();
    hsum  += it.P();
  }
  if(allJetsVec.size() > 0) {
    allPlots["HT"+chTag+name+runPeriod_]->Fill(htsum,getWgt());
    allPlots["HT_all"+name+runPeriod_]->Fill(htsum,getWgt());
    allPlots["H"+chTag+name+runPeriod_]->Fill(hsum,getWgt());
    allPlots["H_all"+name+runPeriod_]->Fill(hsum,getWgt());
    allPlots["j_pt"+chTag+name+runPeriod_]->Fill(allJetsVec[0].getVec().Pt(),getWgt());
    allPlots["j_pt_low"+chTag+name+runPeriod_]->Fill(allJetsVec[0].getVec().Pt(),getWgt());
    //allPlots["j_pt"+chTag+name+"_no_weight"+runPeriod_]->Fill(allJetsVec[0].getVec().Pt(), norm_);
    allPlots["j_pt_all"+name+runPeriod_]->Fill(allJetsVec[0].getVec().Pt(),getWgt());

    std::pair<float,float> sphericity = Sphericity(allJetsVec);
    if(sphericity.first>=0) {
      allPlots["j_sphericity"+chTag+name+runPeriod_]->Fill(sphericity.first,getWgt());
      allPlots["j_sphericity_all"+name+runPeriod_]->Fill(sphericity.first,getWgt());
    }
    if(sphericity.second>=0) {
      allPlots["j_planarity"+chTag+name+runPeriod_]->Fill(sphericity.second,getWgt());
      allPlots["j_planarity_all"+name+runPeriod_]->Fill(sphericity.second,getWgt());
    }
  }
  if(lightJetsVec.size() > 0) {
    allPlots["lj_pt_low"+chTag+name+runPeriod_]->Fill(lightJetsVec[0].getVec().Pt(),getWgt());
    allPlots["lj_pt"+chTag+name+runPeriod_]->Fill(lightJetsVec[0].getVec().Pt(),getWgt());
    //allPlots["lj_pt"+chTag+name+"_no_weight"+runPeriod_]->Fill(lightJetsVec[0].getVec().Pt(),norm_);
    allPlots["lj_pt_all"+name+runPeriod_]->Fill(lightJetsVec[0].getVec().Pt(),getWgt());
    allPlots["lj_pt_ch"+chTag+name+runPeriod_]->Fill(lightJetsVec[0].getChargedPt(),getWgt());
    allPlots["lj_pt_ch_all"+name+runPeriod_]->Fill(lightJetsVec[0].getChargedPt(),getWgt());
  }
  if(kJetsVec.size() > 0) {
    allPlots["kj_pt_low"+chTag+name+runPeriod_]->Fill(kJetsVec[0].getVec().Pt(),getWgt());
    allPlots["kj_pt"+chTag+name+runPeriod_]->Fill(kJetsVec[0].getVec().Pt(),getWgt());
    //allPlots["kj_pt"+chTag+name+"_no_weight"+runPeriod_]->Fill(kJetsVec[0].getVec().Pt(),norm_);
    allPlots["kj_pt_all"+name+runPeriod_]->Fill(kJetsVec[0].getVec().Pt(),getWgt());
    allPlots["kj_pt_ch"+chTag+name+runPeriod_]->Fill(kJetsVec[0].getChargedPt(),getWgt());
    allPlots["kj_pt_ch_all"+name+runPeriod_]->Fill(kJetsVec[0].getChargedPt(),getWgt());
    for(size_t ij = 0; ij < kJetsVec.size(); ij++) {
      float csv = kJetsVec.at(ij).getCSV();
      allPlots["csv"+chTag+name+runPeriod_]->Fill(csv,getWgt());
    }
    for(auto it : kJetsVec)
      allPlots["csv_all"+name+runPeriod_]->Fill(it.getCSV(),getWgt());
  }
}

//Called by Fill(std::vector<pfTrack> puCands, Leptons lep, TString chTag, TString name)
void StdPlots::Fill(std::vector<pfTrack> &pfCands, TString chTag, TString name) {
  if(!isGood_) return;
  if(debug_) std::cout << "Filling meson only" << std::endl;
  float wgt = norm_ * sfs_.first * puWgt_ * top_pt_wgt_;
  if(!name.EqualTo("")) name = "_" + name;
  if(name.Contains("jpsi")) {
    if(debug_) std::cout << "Filling J/Psi" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;
    if(debug_) std::cout << "is " << name << std::endl;
    if(pfCands[0].getPfid() != -pfCands[1].getPfid()) return;
    TLorentzVector jpsi = pfCands[0].getVec() + pfCands[1].getVec();
    //float mass12((pfCands[0].getVec() + pfCands[1].getVec()).M());
    //J/Psi mass in slightly wide window
    if(jpsi.M()<2.8 || jpsi.M()>3.4) return; //Loose window for mass only
    if(debug_) std::cout << "is J/Psi" << name << std::endl;
    allPlots["massJPsi"+chTag+name+runPeriod_]->Fill(jpsi.M(),getWgt());
    allPlots["massJPsi_all"+name+runPeriod_]->Fill(jpsi.M(),getWgt());
    //float pt12((pfCands[0].getVec() + pfCands[1].getVec()).Pt());

    //if(jpsi.M()<3.0 || jpsi.M()>3.2) return; //Window in Elvire's AN
    if(abs(jpsi.M()-3.097) > 0.1) return;
    allPlots["JPsi_pt"+chTag+name+runPeriod_]->Fill(jpsi.Pt(),getWgt());
    allPlots["JPsi_pt_all"+name+runPeriod_]->Fill(jpsi.Pt(),getWgt());
    allPlots["JPsi_eta"+chTag+name+runPeriod_]->Fill(jpsi.Eta(),getWgt());
    allPlots["JPsi_eta_all"+name+runPeriod_]->Fill(jpsi.Eta(),getWgt());
    allPlots["JPsi_phi"+chTag+name+runPeriod_]->Fill(jpsi.Phi(),getWgt());
    allPlots["JPsi_phi_all"+name+runPeriod_]->Fill(jpsi.Phi(),getWgt());

    allPlots["JPsi_mu1_pt"+chTag+name+runPeriod_]->Fill(pfCands[0].Pt(),getWgt());
    allPlots["JPsi_mu1_pt_all"+name+runPeriod_]->Fill(pfCands[0].Pt(),getWgt());
    allPlots["JPsi_mu2_pt"+chTag+name+runPeriod_]->Fill(pfCands[1].Pt(),getWgt());
    allPlots["JPsi_mu2_pt_all"+name+runPeriod_]->Fill(pfCands[1].Pt(),getWgt());
    allPlots["JPsi_mu1_eta"+chTag+name+runPeriod_]->Fill(pfCands[0].Eta(),getWgt());
    allPlots["JPsi_mu1_eta_all"+name+runPeriod_]->Fill(pfCands[0].Eta(),getWgt());
    allPlots["JPsi_mu2_eta"+chTag+name+runPeriod_]->Fill(pfCands[1].Eta(),getWgt());
    allPlots["JPsi_mu2_eta_all"+name+runPeriod_]->Fill(pfCands[1].Eta(),getWgt());
    allPlots["JPsi_mu1_phi"+chTag+name+runPeriod_]->Fill(pfCands[0].Phi(),getWgt());
    allPlots["JPsi_mu1_phi_all"+name+runPeriod_]->Fill(pfCands[0].Phi(),getWgt());
    allPlots["JPsi_mu2_phi"+chTag+name+runPeriod_]->Fill(pfCands[1].Phi(),getWgt());
    allPlots["JPsi_mu2_phi_all"+name+runPeriod_]->Fill(pfCands[1].Phi(),getWgt());

    /*
    allPlots["JPsi_mu1_quality"+chTag+name+runPeriod_]->Fill(pfCands[0].getQuality(),getWgt());
    allPlots["JPsi_mu1_quality_all"+name+runPeriod_]->Fill(pfCands[0].getQuality(),getWgt());
    allPlots["JPsi_mu2_quality"+chTag+name+runPeriod_]->Fill(pfCands[1].getQuality(),getWgt());
    allPlots["JPsi_mu2_quality_all"+name+runPeriod_]->Fill(pfCands[1].getQuality(),getWgt());
    allPlots["JPsi_mu1_highPurity"+chTag+name+runPeriod_]->Fill(pfCands[0].highPurity(),getWgt());
    allPlots["JPsi_mu1_highPurity_all"+name+runPeriod_]->Fill(pfCands[0].highPurity(),getWgt());
    allPlots["JPsi_mu2_highPurity"+chTag+name+runPeriod_]->Fill(pfCands[1].highPurity(),getWgt());
    allPlots["JPsi_mu2_highPurity_all"+name+runPeriod_]->Fill(pfCands[1].highPurity(),getWgt());
    */

    allPlots["JPsi_l3d"+chTag+name+runPeriod_]->Fill(pfCands[0].getL3D(),getWgt());
    allPlots["JPsi_l3d_all"+name+runPeriod_]->Fill(pfCands[0].getL3D(),getWgt());
    allPlots["JPsi_sigmal3d"+chTag+name+runPeriod_]->Fill(pfCands[0].getSigmaL3D(),getWgt());
    allPlots["JPsi_sigmal3d_all"+name+runPeriod_]->Fill(pfCands[0].getSigmaL3D(),getWgt());

    /*
    allPlots["pf_dxy"+chTag+name+runPeriod_]->Fill(abs(pfmuCands[i].getDxy()),wgt);
    allPlots["pf_dxy_all"+name+runPeriod_]->Fill(abs(pfmuCands[i].getDxy()),wgt);
    allPlots["pf_dz"+chTag+name+runPeriod_]->Fill(abs(pfmuCands[i].getDz()),wgt);
    allPlots["pf_dz_all"+name+runPeriod_]->Fill(abs(pfmuCands[i].getDz()),wgt);
    allPlots["pf_dxyE"+chTag+name+runPeriod_]->Fill(abs(pfmuCands[i].getDxyE()),wgt);
    allPlots["pf_dxyE_all"+name+runPeriod_]->Fill(abs(pfmuCands[i].getDxyE()),wgt);
    allPlots["pf_dzE"+chTag+name+runPeriod_]->Fill(abs(pfmuCands[i].getDzE()),wgt);
    allPlots["pf_dzE_all"+name+runPeriod_]->Fill(abs(pfmuCands[i].getDzE()),wgt);
    allPlots["pf_dz_sig"+chTag+name+runPeriod_]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
    allPlots["pf_dz_sig_all"+name+runPeriod_]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
    allPlots["pf_dxy_sig"+chTag+name+runPeriod_]->Fill(abs(pfmuCands[i].getDxy())/abs(pfmuCands[i].getDxyE()),wgt);
    allPlots["pf_dxy_sig_all"+name+runPeriod_]->Fill(abs(pfmuCands[i].getDxy())/abs(pfmuCands[i].getDxyE()),wgt);
    allPlots["pf_dz_sig"+chTag+name+runPeriod_]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
    allPlots["pf_dz_sig_all"+name+runPeriod_]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
    */
  }
  else if(name.Contains("meson")) {
    if(debug_) std::cout << "Filling D_meson" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;
    if(pfCands.size()<2) return;
    if(pfCands[0].getMotherId()!=421) return;
    if(pfCands[1].getMotherId()!=421) return;
    TLorentzVector D0 = pfCands[0].getVec() + pfCands[1].getVec();
    float mass12 = D0.M();
    pfTrack *pi, *k;
    if(pfCands[0].M() < 0.15) { //m_pi = 0.1396 and m_K = 0.4937
      pi = &pfCands[0];
      k = &pfCands[1];
    }
    else {
      k = &pfCands[0];
      pi = &pfCands[1];
    }
    //if(pi->Pt() < 6) return; //cut from ovelaying D0_pi and D0_formDs_pi
    //if(!pi->highPurity()) return;
    //if(!k->highPurity()) return;

    //if(mass12<1.7 || mass12>2.0) return; //D^0 mass window
    if(pfCands.size()==2 && mass12>1.7 && mass12<2.0) { //Plot D0 only
      if(debug_) std::cout << "Filling D0" << chTag << name << runPeriod_ << " with wgt=" << getWgt() << std::endl;
      allPlots["massD0"+chTag+name+runPeriod_]->Fill(mass12,getWgt());
      allPlots["massD0_all"+name+runPeriod_]->Fill(mass12,getWgt());

      allPlots["D0_pt"+chTag+name+runPeriod_]->Fill(D0.Pt(),getWgt());
      allPlots["D0_pt_all"+name+runPeriod_]->Fill(D0.Pt(),getWgt());
      allPlots["D0_pi_pt"+chTag+name+runPeriod_]->Fill(pi->Pt(),getWgt());
      allPlots["D0_pi_pt_all"+name+runPeriod_]->Fill(pi->Pt(),getWgt());
      allPlots["D0_K_pt"+chTag+name+runPeriod_]->Fill(k->Pt(),getWgt());
      allPlots["D0_K_pt_all"+name+runPeriod_]->Fill(k->Pt(),getWgt());
      allPlots["D0_p"+chTag+name+runPeriod_]->Fill(D0.P(),getWgt());
      allPlots["D0_p_all"+name+runPeriod_]->Fill(D0.P(),getWgt());
      allPlots["D0_pi_p"+chTag+name+runPeriod_]->Fill(pi->P(),getWgt());
      allPlots["D0_pi_p_all"+name+runPeriod_]->Fill(pi->P(),getWgt());
      allPlots["D0_K_p"+chTag+name+runPeriod_]->Fill(k->P(),getWgt());
      allPlots["D0_K_p_all"+name+runPeriod_]->Fill(k->P(),getWgt());
      allPlots["D0_pi_eta"+chTag+name+runPeriod_]->Fill(pi->Eta(),getWgt());
      allPlots["D0_pi_eta_all"+name+runPeriod_]->Fill(pi->Eta(),getWgt());
      allPlots["D0_K_eta"+chTag+name+runPeriod_]->Fill(k->Eta(),getWgt());
      allPlots["D0_K_eta_all"+name+runPeriod_]->Fill(k->Eta(),getWgt());
      allPlots["D0_pi_phi"+chTag+name+runPeriod_]->Fill(pi->Phi(),getWgt());
      allPlots["D0_pi_phi_all"+name+runPeriod_]->Fill(pi->Phi(),getWgt());
      allPlots["D0_K_phi"+chTag+name+runPeriod_]->Fill(k->Phi(),getWgt());
      allPlots["D0_K_phi_all"+name+runPeriod_]->Fill(k->Phi(),getWgt());
      allPlots["D0_dRpik"+chTag+name+runPeriod_]->Fill(pi->DeltaR(*k),getWgt());
      allPlots["D0_dRpik_all"+name+runPeriod_]->Fill(pi->DeltaR(*k),getWgt());
      allPlots["D0_dEtapik"+chTag+name+runPeriod_]->Fill((pi->Eta())-(k->Eta()),getWgt());
      allPlots["D0_dEtapik_all"+name+runPeriod_]->Fill((pi->Eta())-(k->Eta()),getWgt());
      allPlots["D0_dPhipik"+chTag+name+runPeriod_]->Fill(pi->getVec().DeltaPhi(k->getVec()),getWgt());
      allPlots["D0_dPhipik_all"+name+runPeriod_]->Fill(pi->getVec().DeltaPhi(k->getVec()),getWgt());

      /*
      allPlots["D0_pi_quality"+chTag+name+runPeriod_]->Fill(pi->getQuality(),getWgt());
      allPlots["D0_pi_quality_all"+name+runPeriod_]->Fill(pi->getQuality(),getWgt());
      allPlots["D0_K_quality"+chTag+name+runPeriod_]->Fill(k->getQuality(),getWgt());
      allPlots["D0_K_quality_all"+name+runPeriod_]->Fill(k->getQuality(),getWgt());
      allPlots["D0_pi_highPurity"+chTag+name+runPeriod_]->Fill(pi->highPurity(),getWgt());
      allPlots["D0_pi_highPurity_all"+name+runPeriod_]->Fill(pi->highPurity(),getWgt());
      allPlots["D0_K_highPurity_all"+name+runPeriod_]->Fill(k->highPurity(),getWgt());
      */

      if(mass12>1.8 && mass12<1.9) {
        allPlots["D0_l3d"+chTag+name+runPeriod_]->Fill(pi->getL3D(),getWgt());
        allPlots["D0_l3d_all"+name+runPeriod_]->Fill(pi->getL3D(),getWgt());
        allPlots["D0_sigmal3d"+chTag+name+runPeriod_]->Fill(pi->getSigmaL3D(),getWgt());
        allPlots["D0_sigmal3d_all"+name+runPeriod_]->Fill(pi->getSigmaL3D(),getWgt());
        if(name.Contains("gmeson") && abs(pfCands[0].getMotherId())==421) {
          allPlots["D0_l3d_good_all"+name+runPeriod_]->Fill(pi->getL3D(),getWgt());
          allPlots["D0_sigmal3d_good_all"+name+runPeriod_]->Fill(pi->getSigmaL3D(),getWgt());
        }
        else if(name.Contains("gmeson")) {
          allPlots["D0_l3d_bad_all"+name+runPeriod_]->Fill(pi->getL3D(),getWgt());
          allPlots["D0_sigmal3d_bad_all"+name+runPeriod_]->Fill(pi->getSigmaL3D(),getWgt());
        }
      }
    }
    if(pfCands.size()<3) return; // D^* and flavor tagging D^0 
    if(abs(pfCands[2].getPdgId())==13 && pfCands[2].getMotherId()==42113) {
    //if(abs(pfCands[0].getMotherId())==42113) {
      //if(pfCands[1].charge() == pfCands[2].charge()) { //reinforce kaon and lepton have same sign
        if(mass12>1.7 && mass12<2.0) { //Plot D0 only
          allPlots["massD0_mu_tag"+chTag+name+runPeriod_]->Fill(mass12,getWgt());
          allPlots["massD0_lep_tag"+chTag+name+runPeriod_]->Fill(mass12,getWgt());
          allPlots["massD0_mu_tag_all"+name+runPeriod_]->Fill(mass12,getWgt());
          allPlots["massD0_lep_tag_all"+name+runPeriod_]->Fill(mass12,getWgt());
        }
        if(mass12>1.8 && mass12<1.93) { //Plot D0 kinematics
          allPlots["D0_mu_tag_pt"+chTag+name+runPeriod_]->Fill(D0.Pt(),getWgt());
          allPlots["D0_mu_tag_pt_all"+name+runPeriod_]->Fill(D0.Pt(),getWgt());
          allPlots["D0_mu_tag_pi_pt"+chTag+name+runPeriod_]->Fill(pi->Pt(),getWgt());
          allPlots["D0_mu_tag_pi_pt_all"+name+runPeriod_]->Fill(pi->Pt(),getWgt());
          allPlots["D0_mu_tag_K_pt"+chTag+name+runPeriod_]->Fill(k->Pt(),getWgt());
          allPlots["D0_mu_tag_K_pt_all"+name+runPeriod_]->Fill(k->Pt(),getWgt());
          /*
          allPlots["massD0_mu_tagB"+chTag+name+runPeriod_]->Fill((D0+pfCands[2].getVec()).M(),getWgt());
          allPlots["massD0_mu_tagB_all"+name+runPeriod_]->Fill((D0+pfCands[2].getVec()).M(),getWgt());
          allPlots["D0_mu_highPurity"+chTag+name+runPeriod_]->Fill(pfCands[2].highPurity(),getWgt());
          allPlots["D0_mu_highPurity_all"+name+runPeriod_]->Fill(pfCands[2].highPurity(),getWgt());
          */
          allPlots["D0_mu_pt"+chTag+name+runPeriod_]->Fill(pfCands[2].Pt(),getWgt());
          allPlots["D0_mu_pt_all"+name+runPeriod_]->Fill(pfCands[2].Pt(),getWgt());
          //allPlots["D0_mu_eta"+chTag+name+runPeriod_]->Fill(pfCands[2].Eta(),getWgt());
          //allPlots["D0_mu_eta_all"+name+runPeriod_]->Fill(pfCands[2].Eta(),getWgt());
          allPlots["D0_mu_phi"+chTag+name+runPeriod_]->Fill(pfCands[2].Phi(),getWgt());
          allPlots["D0_mu_phi_all"+name+runPeriod_]->Fill(pfCands[2].Phi(),getWgt());
          allPlots["D0_mu_tag_pi_eta"+chTag+name+runPeriod_]->Fill(pi->Eta(),getWgt());
          allPlots["D0_mu_tag_pi_eta_all"+name+runPeriod_]->Fill(pi->Eta(),getWgt());
          allPlots["D0_mu_tag_pi_phi"+chTag+name+runPeriod_]->Fill(pi->Phi(),getWgt());
          allPlots["D0_mu_tag_pi_phi_all"+name+runPeriod_]->Fill(pi->Phi(),getWgt());
          allPlots["D0_mu_tag_K_eta"+chTag+name+runPeriod_]->Fill(k->Eta(),getWgt());
          allPlots["D0_mu_tag_K_eta_all"+name+runPeriod_]->Fill(k->Eta(),getWgt());
          allPlots["D0_mu_tag_K_phi"+chTag+name+runPeriod_]->Fill(k->Phi(),getWgt());
          allPlots["D0_mu_tag_K_phi_all"+name+runPeriod_]->Fill(k->Phi(),getWgt());
        }
      //}
      if(mass12>1.8 && mass12<1.93) {
        allPlots["D0_mu_tag_l3d"+chTag+name+runPeriod_]->Fill(pi->getL3D(),getWgt());
        allPlots["D0_mu_tag_l3d_all"+name+runPeriod_]->Fill(pi->getL3D(),getWgt());
        allPlots["D0_mu_tag_sigmal3d"+chTag+name+runPeriod_]->Fill(pi->getSigmaL3D(),getWgt());
        allPlots["D0_mu_tag_sigmal3d_all"+name+runPeriod_]->Fill(pi->getSigmaL3D(),getWgt());
      }
    }
    /*
    if(abs(pfCands[2].getPdgId())==11) {
      //if(pfCands[1].charge() == pfCands[2].charge()) { //reinforce kaon and lepton have same sign
        if(pfCands.size()==2 && mass12>1.7 && mass12<2.0) { //Plot D0 only
          allPlots["massD0_e_tag"+chTag+name+runPeriod_]->Fill(mass12,getWgt());
          allPlots["massD0_lep_tag"+chTag+name+runPeriod_]->Fill(mass12,getWgt());
          allPlots["massD0_e_tag_all"+name+runPeriod_]->Fill(mass12,getWgt());
          allPlots["massD0_lep_tag_all"+name+runPeriod_]->Fill(mass12,getWgt());
        }
      //}
    }
    */
    //recheck first pi and K
    if(pfCands[0].getPdgId()*pfCands[1].getPdgId() != -211*211) return;
    //if (mass12<1.65 && mass12>2.0) return;

    if(fabs(pfCands[2].getPdgId()) != 211) return; //reinforce pion
    if(pfCands[2].getMotherId()!=413) return;
    // Kaon and pion have opposite charges
    // I.e. correct mass assumption
    //if(pfCands[1].charge() != -pfCands[2].charge()) return;
    if(k->charge() != -pfCands[2].charge()) return;

    pfTrack *pi2;
    pi2 = &pfCands[2];
    //if(!pi2->highPurity()) return;
    //if(pi2->Pt() > 10) return; //cut from ovelaying Ds_pi2 and Ds_formDs_pi2

    if(debug_) std::cout << "Filling D*" << chTag << name << runPeriod_ << " with wgt=" << getWgt() << std::endl;
    TLorentzVector p_cand = pfCands[0].getVec()+pfCands[1].getVec()+pfCands[2].getVec();
    allPlots["massDs"+chTag+name+runPeriod_]->Fill(p_cand.M(), getWgt());
    allPlots["massDs_all"+name+runPeriod_]->Fill(p_cand.M(), getWgt());

    allPlots["Ds_pi_pt"+chTag+name+runPeriod_]->Fill(pi->Pt(),getWgt());
    allPlots["Ds_pi_pt_all"+name+runPeriod_]->Fill(pi->Pt(),getWgt());
    allPlots["Ds_K_pt"+chTag+name+runPeriod_]->Fill(k->Pt(),getWgt());
    allPlots["Ds_K_pt_all"+name+runPeriod_]->Fill(k->Pt(),getWgt());
    allPlots["Ds_pi2_pt"+chTag+name+runPeriod_]->Fill(pi2->Pt(),getWgt());
    allPlots["Ds_pi2_pt_all"+name+runPeriod_]->Fill(pi2->Pt(),getWgt());
    allPlots["Ds_pi_eta"+chTag+name+runPeriod_]->Fill(pi->Eta(),getWgt());
    allPlots["Ds_pi_eta_all"+name+runPeriod_]->Fill(pi->Eta(),getWgt());
    allPlots["Ds_K_eta"+chTag+name+runPeriod_]->Fill(k->Eta(),getWgt());
    allPlots["Ds_K_eta_all"+name+runPeriod_]->Fill(k->Eta(),getWgt());
    allPlots["Ds_pi_phi"+chTag+name+runPeriod_]->Fill(pi->Phi(),getWgt());
    allPlots["Ds_pi_phi_all"+name+runPeriod_]->Fill(pi->Phi(),getWgt());
    allPlots["Ds_K_phi"+chTag+name+runPeriod_]->Fill(k->Phi(),getWgt());
    allPlots["Ds_K_phi_all"+name+runPeriod_]->Fill(k->Phi(),getWgt());

    /*
    allPlots["Ds_pi_quality"+chTag+name+runPeriod_]->Fill(pi->getQuality(),getWgt());
    allPlots["Ds_pi_quality_all"+name+runPeriod_]->Fill(pi->getQuality(),getWgt());
    allPlots["Ds_K_quality"+chTag+name+runPeriod_]->Fill(k->getQuality(),getWgt());
    allPlots["Ds_K_quality_all"+name+runPeriod_]->Fill(k->getQuality(),getWgt());
    allPlots["Ds_pi_highPurity"+chTag+name+runPeriod_]->Fill(pi->highPurity(),getWgt());
    allPlots["Ds_pi_highPurity_all"+name+runPeriod_]->Fill(pi->highPurity(),getWgt());
    allPlots["Ds_K_highPurity_all"+name+runPeriod_]->Fill(k->highPurity(),getWgt());
    */

    if(fabs(mass12-1.864) > 0.10) return; // mass window cut
    float deltam = p_cand.M() - mass12;

    if(debug_) std::cout << "Filling D*-D0" << chTag << name << runPeriod_ << " with wgt=" << getWgt() << std::endl;
    allPlots["massDsmD0loose"+chTag+name+runPeriod_]->Fill(deltam, getWgt());
    allPlots["massDsmD0loose_all"+name+runPeriod_]->Fill(deltam, getWgt());

    if(deltam>0.14 && deltam<0.15) {
      allPlots["massD0_fromDsloose"+chTag+name+runPeriod_]->Fill(D0.M(), getWgt());
      allPlots["massD0_fromDsloose_all"+name+runPeriod_]->Fill(D0.M(), getWgt());

      allPlots["massDs_fromDsloose"+chTag+name+runPeriod_]->Fill(p_cand.M(), getWgt());
      allPlots["massDs_fromDsloose_all"+name+runPeriod_]->Fill(p_cand.M(), getWgt());

      //D^0 from D*-D^0 loose
      allPlots["D0_fromDsloose_pi_pt"+chTag+name+runPeriod_]->Fill(pi->Pt(),getWgt());
      allPlots["D0_fromDsloose_pi_pt_all"+name+runPeriod_]->Fill(pi->Pt(),getWgt());
      allPlots["D0_fromDsloose_K_pt"+chTag+name+runPeriod_]->Fill(k->Pt(),getWgt());
      allPlots["D0_fromDsloose_K_pt_all"+name+runPeriod_]->Fill(k->Pt(),getWgt());
      allPlots["D0_fromDsloose_pi_eta"+chTag+name+runPeriod_]->Fill(pi->Eta(),getWgt());
      allPlots["D0_fromDsloose_pi_eta_all"+name+runPeriod_]->Fill(pi->Eta(),getWgt());
      allPlots["D0_fromDsloose_K_eta"+chTag+name+runPeriod_]->Fill(k->Eta(),getWgt());
      allPlots["D0_fromDsloose_K_eta_all"+name+runPeriod_]->Fill(k->Eta(),getWgt());
      allPlots["D0_fromDsloose_pi_phi"+chTag+name+runPeriod_]->Fill(pi->Phi(),getWgt());
      allPlots["D0_fromDsloose_pi_phi_all"+name+runPeriod_]->Fill(pi->Phi(),getWgt());
      allPlots["D0_fromDsloose_K_phi"+chTag+name+runPeriod_]->Fill(k->Phi(),getWgt());
      allPlots["D0_fromDsloose_K_phi_all"+name+runPeriod_]->Fill(k->Phi(),getWgt());

      //D* from D*-D^0 loose
      allPlots["Ds_fromDsloose_pi_pt"+chTag+name+runPeriod_]->Fill(pi->Pt(),getWgt());
      allPlots["Ds_fromDsloose_pi_pt_all"+name+runPeriod_]->Fill(pi->Pt(),getWgt());
      allPlots["Ds_fromDsloose_K_pt"+chTag+name+runPeriod_]->Fill(k->Pt(),getWgt());
      allPlots["Ds_fromDsloose_K_pt_all"+name+runPeriod_]->Fill(k->Pt(),getWgt());
      allPlots["Ds_fromDsloose_pi_eta"+chTag+name+runPeriod_]->Fill(pi->Eta(),getWgt());
      allPlots["Ds_fromDsloose_pi_eta_all"+name+runPeriod_]->Fill(pi->Eta(),getWgt());
      allPlots["Ds_fromDsloose_K_eta"+chTag+name+runPeriod_]->Fill(k->Eta(),getWgt());
      allPlots["Ds_fromDsloose_K_eta_all"+name+runPeriod_]->Fill(k->Eta(),getWgt());
      allPlots["Ds_fromDsloose_pi_phi"+chTag+name+runPeriod_]->Fill(pi->Phi(),getWgt());
      allPlots["Ds_fromDsloose_pi_phi_all"+name+runPeriod_]->Fill(pi->Phi(),getWgt());
      allPlots["Ds_fromDsloose_K_phi"+chTag+name+runPeriod_]->Fill(k->Phi(),getWgt());
      allPlots["Ds_fromDsloose_K_phi_all"+name+runPeriod_]->Fill(k->Phi(),getWgt());
    }

    if(deltam>0.15 && deltam<0.16) { //window of 0.01 to match window of good peak
      allPlots["massD0_notDsloose"+chTag+name+runPeriod_]->Fill(mass12, getWgt());
      allPlots["massD0_notDsloose_all"+name+runPeriod_]->Fill(mass12, getWgt());
      allPlots["massDs_notDsloose"+chTag+name+runPeriod_]->Fill(p_cand.M(), getWgt());
      allPlots["massDs_notDsloose_all"+name+runPeriod_]->Fill(p_cand.M(), getWgt());
      if(fabs(mass12-1.864) < 0.05) {
        allPlots["massD0_notDs"+chTag+name+runPeriod_]->Fill(mass12, getWgt());
        allPlots["massD0_notDs_all"+name+runPeriod_]->Fill(mass12, getWgt());
        allPlots["massDs_notDs"+chTag+name+runPeriod_]->Fill(p_cand.M(), getWgt());
        allPlots["massDs_notDs_all"+name+runPeriod_]->Fill(p_cand.M(), getWgt());
        allPlots["massD0_notDs"+chTag+name+runPeriod_]->Fill(mass12, getWgt());
        allPlots["massD0_notDs_all"+name+runPeriod_]->Fill(mass12, getWgt());
        allPlots["massDs_notDs"+chTag+name+runPeriod_]->Fill(p_cand.M(), getWgt());
        allPlots["massDs_notDs_all"+name+runPeriod_]->Fill(p_cand.M(), getWgt());

        allPlots["D0_notDs_pi_pt"+chTag+name+runPeriod_]->Fill(pi->Pt(),getWgt());
        allPlots["D0_notDs_pi_pt_all"+name+runPeriod_]->Fill(pi->Pt(),getWgt());
        allPlots["D0_notDs_K_pt"+chTag+name+runPeriod_]->Fill(k->Pt(),getWgt());
        allPlots["D0_notDs_K_pt_all"+name+runPeriod_]->Fill(k->Pt(),getWgt());
        allPlots["D0_notDs_pi_eta"+chTag+name+runPeriod_]->Fill(pi->Eta(),getWgt());
        allPlots["D0_notDs_pi_eta_all"+name+runPeriod_]->Fill(pi->Eta(),getWgt());
        allPlots["D0_notDs_K_eta"+chTag+name+runPeriod_]->Fill(k->Eta(),getWgt());
        allPlots["D0_notDs_K_eta_all"+name+runPeriod_]->Fill(k->Eta(),getWgt());
        allPlots["D0_notDs_pi_phi"+chTag+name+runPeriod_]->Fill(pi->Phi(),getWgt());
        allPlots["D0_notDs_pi_phi_all"+name+runPeriod_]->Fill(pi->Phi(),getWgt());
        allPlots["D0_notDs_K_phi"+chTag+name+runPeriod_]->Fill(k->Phi(),getWgt());
        allPlots["D0_notDs_K_phi_all"+name+runPeriod_]->Fill(k->Phi(),getWgt());

        //D* not D*-D^0
        allPlots["Ds_notDs_pi_pt"+chTag+name+runPeriod_]->Fill(pi->Pt(),getWgt());
        allPlots["Ds_notDs_pi_pt_all"+name+runPeriod_]->Fill(pi->Pt(),getWgt());
        allPlots["Ds_notDs_pi2_pt"+chTag+name+runPeriod_]->Fill(pi2->Pt(),getWgt());
        allPlots["Ds_notDs_pi_pt_all"+name+runPeriod_]->Fill(pi->Pt(),getWgt());
        allPlots["Ds_notDs_K_pt"+chTag+name+runPeriod_]->Fill(k->Pt(),getWgt());
        allPlots["Ds_notDs_K_pt_all"+name+runPeriod_]->Fill(k->Pt(),getWgt());
        allPlots["Ds_notDs_pi_eta"+chTag+name+runPeriod_]->Fill(pi->Eta(),getWgt());
        allPlots["Ds_notDs_pi_eta_all"+name+runPeriod_]->Fill(pi->Eta(),getWgt());
        allPlots["Ds_notDs_K_eta"+chTag+name+runPeriod_]->Fill(k->Eta(),getWgt());
        allPlots["Ds_notDs_K_eta_all"+name+runPeriod_]->Fill(k->Eta(),getWgt());
        allPlots["Ds_notDs_pi_phi"+chTag+name+runPeriod_]->Fill(pi->Phi(),getWgt());
        allPlots["Ds_notDs_pi_phi_all"+name+runPeriod_]->Fill(pi->Phi(),getWgt());
        allPlots["Ds_notDs_K_phi"+chTag+name+runPeriod_]->Fill(k->Phi(),getWgt());
        allPlots["Ds_notDs_K_phi_all"+name+runPeriod_]->Fill(k->Phi(),getWgt());
      }
    }
    if(fabs(mass12-1.864) > 0.05) return; // tighter mass window cut
    if(debug_) std::cout << "Masses: " << pfCands[0].getVec().M() << " ";
    if(debug_) std::cout << pfCands[1].getVec().M() << " ";
    if(debug_) std::cout << pfCands[2].getVec().M() << std::endl;

    allPlots["massDsmD0"+chTag+name+runPeriod_]->Fill(deltam, getWgt());
    allPlots["massDsmD0_all"+name+runPeriod_]->Fill(deltam, getWgt());

    if(deltam<0.14 || deltam>0.15) return;
    allPlots["massD0_fromDs"+chTag+name+runPeriod_]->Fill(mass12, getWgt());
    allPlots["massD0_fromDs_all"+name+runPeriod_]->Fill(mass12, getWgt());
    allPlots["massDs_fromDs"+chTag+name+runPeriod_]->Fill(p_cand.M(), getWgt());
    allPlots["massDs_fromDs_all"+name+runPeriod_]->Fill(p_cand.M(), getWgt());

    allPlots["D0_fromDs_pi_pt"+chTag+name+runPeriod_]->Fill(pi->Pt(),getWgt());
    allPlots["D0_fromDs_pi_pt_all"+name+runPeriod_]->Fill(pi->Pt(),getWgt());
    allPlots["D0_fromDs_K_pt"+chTag+name+runPeriod_]->Fill(k->Pt(),getWgt());
    allPlots["D0_fromDs_K_pt_all"+name+runPeriod_]->Fill(k->Pt(),getWgt());
    allPlots["D0_fromDs_pi_eta"+chTag+name+runPeriod_]->Fill(pi->Eta(),getWgt());
    allPlots["D0_fromDs_pi_eta_all"+name+runPeriod_]->Fill(pi->Eta(),getWgt());
    allPlots["D0_fromDs_K_eta"+chTag+name+runPeriod_]->Fill(k->Eta(),getWgt());
    allPlots["D0_fromDs_K_eta_all"+name+runPeriod_]->Fill(k->Eta(),getWgt());
    allPlots["D0_fromDs_pi_phi"+chTag+name+runPeriod_]->Fill(pi->Phi(),getWgt());
    allPlots["D0_fromDs_pi_phi_all"+name+runPeriod_]->Fill(pi->Phi(),getWgt());
    allPlots["D0_fromDs_K_phi"+chTag+name+runPeriod_]->Fill(k->Phi(),getWgt());
    allPlots["D0_fromDs_K_phi_all"+name+runPeriod_]->Fill(k->Phi(),getWgt());

    //D* from D*-D^0
    allPlots["Ds_fromDs_pi_pt"+chTag+name+runPeriod_]->Fill(pi->Pt(),getWgt());
    allPlots["Ds_fromDs_pi_pt_all"+name+runPeriod_]->Fill(pi->Pt(),getWgt());
    allPlots["Ds_fromDs_K_pt"+chTag+name+runPeriod_]->Fill(k->Pt(),getWgt());
    allPlots["Ds_fromDs_K_pt_all"+name+runPeriod_]->Fill(k->Pt(),getWgt());
    allPlots["Ds_fromDs_pi2_pt"+chTag+name+runPeriod_]->Fill(pi2->Pt(),getWgt());
    allPlots["Ds_fromDs_pi2_pt_all"+name+runPeriod_]->Fill(pi2->Pt(),getWgt());
    allPlots["Ds_fromDs_pi_eta"+chTag+name+runPeriod_]->Fill(pi->Eta(),getWgt());
    allPlots["Ds_fromDs_pi_eta_all"+name+runPeriod_]->Fill(pi->Eta(),getWgt());
    allPlots["Ds_fromDs_K_eta"+chTag+name+runPeriod_]->Fill(k->Eta(),getWgt());
    allPlots["Ds_fromDs_K_eta_all"+name+runPeriod_]->Fill(k->Eta(),getWgt());
    allPlots["Ds_fromDs_pi_phi"+chTag+name+runPeriod_]->Fill(pi->Phi(),getWgt());
    allPlots["Ds_fromDs_pi_phi_all"+name+runPeriod_]->Fill(pi->Phi(),getWgt());
    allPlots["Ds_fromDs_K_phi"+chTag+name+runPeriod_]->Fill(k->Phi(),getWgt());
    allPlots["Ds_fromDs_K_phi_all"+name+runPeriod_]->Fill(k->Phi(),getWgt());
  }
}

//Called by Fill(std::vector<pfTrack> pfCands, Leptons lep, Jet jet, TString chTag, TString name)
void StdPlots::Fill(std::vector<pfTrack> &pfCands, Leptons lep, TString chTag, TString name) {
  if(!isGood_) return;
  Fill(pfCands, chTag, name); //Fill meson only plots
  Fill(lep, chTag, name);       //Fill lepton only plots
  float wgt = norm_ * sfs_.first * puWgt_ * top_pt_wgt_;
  if(!name.EqualTo("")) name = "_" + name;
  //J/Psi events
  if(name.Contains("jpsi")) {
    if(debug_) std::cout << "Filling J/Psi+lep" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;
    if(debug_) std::cout << "is " << name << std::endl;
    if(pfCands[0].getPfid() != -pfCands[1].getPfid()) return;
    TLorentzVector jpsi = pfCands[0].getVec() + pfCands[1].getVec();
    if(jpsi.M()<3.0 || jpsi.M()>3.2) return; //Window in Elvire's AN
    if(debug_) std::cout << "is J/Psi" << name << std::endl;
    float pt123((jpsi + lep[0].getVec()).Pt());
    float mass123((jpsi + lep[0].getVec()).M());
    int ilep = 0;
    float dRJPsil(jpsi.DeltaR(lep[ilep].getVec()));
    //Find closest isolated lepton in di-lepton event
    if(lep.size()>1) {
      for(ilep = 0; ilep < (int)lep.size(); ilep++) {
        float tmpdRJPsil(jpsi.DeltaR(lep[ilep].getVec()));
        if(tmpdRJPsil < dRJPsil) {
          dRJPsil = tmpdRJPsil;
          mass123 = (jpsi + lep[ilep].getVec()).M();
        }
      }
    }
    float dRmumu(pfCands[0].getVec().DeltaR(pfCands[1].getVec()));
    //J/Psi mass in slightly wide window
    allPlots["massJPsi_l"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
    allPlots["massJPsi_l_all"+name+runPeriod_]->Fill(mass123,getWgt());
    if(lep[ilep].isMuon())
      allPlots["massJPsi_mu"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
    else if(lep[ilep].isElectron())
      allPlots["massJPsi_e"+chTag+name+runPeriod_]->Fill(mass123,getWgt());

    if(lep[ilep].isMuon()) {
      allPlots["JPsi_mu_pt"+chTag+name+runPeriod_]->Fill(pt123,getWgt());
      allPlots["dR_JPsi_mu"+chTag+name+runPeriod_]->Fill(dRJPsil,getWgt());
      allPlots["JPsi_mu_pt_all"+name+runPeriod_]->Fill(pt123,getWgt());
      allPlots["massJPsi_mu_all"+name+runPeriod_]->Fill(mass123,getWgt());
      allPlots["dR_JPsi_mu_all"+name+runPeriod_]->Fill(dRJPsil,getWgt());
    }
    else if(lep[ilep].isElectron()) {
      allPlots["JPsi_e_pt"+chTag+name+runPeriod_]->Fill(pt123,getWgt());
      allPlots["dR_JPsi_e"+chTag+name+runPeriod_]->Fill(dRJPsil,getWgt());
      allPlots["JPsi_e_pt_all"+name+runPeriod_]->Fill(pt123,getWgt());
      allPlots["massJPsi_e_all"+name+runPeriod_]->Fill(mass123,getWgt());
      allPlots["dR_JPsi_e_all"+name+runPeriod_]->Fill(dRJPsil,getWgt());
    }
    if(name.Contains("gjpsi") && abs(pfCands[0].getGenT())==6) {
      if(lep[ilep].getPdgId()*pfCands[0].getGenT() > 0) {
        allPlots["dR_JPsi_l_good"+chTag+name+runPeriod_]->Fill(dRJPsil,getWgt());
        allPlots["dR_JPsi_l_good_all"+name+runPeriod_]->Fill(dRJPsil,getWgt());
        allPlots["JPsi_mu1_pt_good"+chTag+name+runPeriod_]->Fill(pfCands[0].Pt(),getWgt());
        allPlots["JPsi_mu1_pt_good_all"+name+runPeriod_]->Fill(pfCands[0].Pt(),getWgt());
        allPlots["JPsi_mu2_pt_good"+chTag+name+runPeriod_]->Fill(pfCands[1].Pt(),getWgt());
        allPlots["JPsi_mu2_pt_good_all"+name+runPeriod_]->Fill(pfCands[1].Pt(),getWgt());
        allPlots["massJPsi_l_good_dR"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
        allPlots["massJPsi_l_good_dR_all"+name+runPeriod_]->Fill(mass123,getWgt());
      }
      if(lep[ilep].getPdgId()*pfCands[0].getGenT() < 0) {
        allPlots["dR_JPsi_l_bad"+chTag+name+runPeriod_]->Fill(dRJPsil,getWgt());
        allPlots["dR_JPsi_l_bad_all"+name+runPeriod_]->Fill(dRJPsil,getWgt());
        allPlots["JPsi_mu1_pt_bad"+chTag+name+runPeriod_]->Fill(pfCands[0].Pt(),getWgt());
        allPlots["JPsi_mu1_pt_bad_all"+name+runPeriod_]->Fill(pfCands[0].Pt(),getWgt());
        allPlots["JPsi_mu2_pt_bad"+chTag+name+runPeriod_]->Fill(pfCands[1].Pt(),getWgt());
        allPlots["JPsi_mu2_pt_bad_all"+name+runPeriod_]->Fill(pfCands[1].Pt(),getWgt());
        allPlots["massJPsi_l_bad_dR"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
        allPlots["massJPsi_l_bad_dR_all"+name+runPeriod_]->Fill(mass123,getWgt());
      }
    }
    allPlots["JPsi_l_pt"+chTag+name+runPeriod_]->Fill(pt123,getWgt());
    allPlots["dR_JPsi_l"+chTag+name+runPeriod_]->Fill(dRJPsil,getWgt());
    allPlots["JPsi_l_pt_all"+name+runPeriod_]->Fill(pt123,getWgt());
    allPlots["massJPsi_l_all"+name+runPeriod_]->Fill(mass123,getWgt());
    allPlots["dR_JPsi_l_all"+name+runPeriod_]->Fill(dRJPsil,getWgt());

    allPlots["dR_JPsi_mumu"+chTag+name+runPeriod_]->Fill(dRmumu,getWgt());
    //allPlots["JPsi_mu1_eta"+chTag+name+runPeriod_]->Fill(pfCands[0].Eta(),getWgt());
    //allPlots["JPsi_mu2_eta"+chTag+name+runPeriod_]->Fill(pfCands[1].Eta(),getWgt());
    //allPlots["JPsi_mu1_phi"+chTag+name+runPeriod_]->Fill(pfCands[0].Eta(),getWgt());
    //allPlots["JPsi_mu2_phi"+chTag+name+runPeriod_]->Fill(pfCands[1].Eta(),getWgt());
    float pt_ratio = (pfCands[0].Pt() - pfCands[1].Pt())/(pfCands[0].Pt()+ pfCands[1].Pt());
    allPlots["JPsi_mumu_pt_ratio"+chTag+name+runPeriod_]->Fill(pt_ratio,getWgt());
    allPlots["JPsi_mumu_pt_ratio_all"+name+runPeriod_]->Fill(pt_ratio,getWgt());

    //dR test
    if(dRJPsil<2.0) {
      allPlots["massJPsi_l_low_dR"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
      allPlots["massJPsi_l_low_dR_all"+name+runPeriod_]->Fill(mass123,getWgt());
      if(lep[ilep].isMuon())
        allPlots["massJPsi_mu_low_dR"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
      if(lep[ilep].isElectron())
        allPlots["massJPsi_e_low_dR"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
    }
    else {
      allPlots["massJPsi_l_high_dR"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
      allPlots["massJPsi_l_high_dR_all"+name+runPeriod_]->Fill(mass123,getWgt());
      if(lep[ilep].isMuon())
        allPlots["massJPsi_mu_high_dR"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
      if(lep[ilep].isElectron())
        allPlots["massJPsi_e_high_dR"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
    }
  }
  else if(name.Contains("meson")) {
    if(debug_) std::cout << "Filling D_meson+lep" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;
    if(pfCands.size()<2) return;
    if(pfCands[0].getPdgId()*pfCands[1].getPdgId() != -211*211) return; // D^0 cut
    TLorentzVector D0 = pfCands[0].getVec() + pfCands[1].getVec();
    float mass12 = D0.M();
    if(mass12<1.7 || mass12>2.0) return; //D^0 mass window
    float mass123 = (D0 + lep[0].getVec()).M();
    int ilep = 0;
    float dRD0l(D0.DeltaR(lep[ilep].getVec()));
    //Find closest isolated lepton in di-lepton event
    if(lep.size()>1) {
      for(ilep = 0; ilep < (int)lep.size(); ilep++) {
        float tmpdRD0l(D0.DeltaR(lep[ilep].getVec()));
        if(tmpdRD0l < dRD0l) {
          ilep = ilep; 
          dRD0l = tmpdRD0l;
          mass123 = (D0 + lep[ilep].getVec()).M();
        }
      }
    }
    if(pfCands.size()==2) { //Plot D0 only
      allPlots["massD0_l"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
      allPlots["massD0_l_all"+name+runPeriod_]->Fill(mass123,getWgt());
      /*
      if(lep[ilep].isMuon())
        allPlots["massD0_mu"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
      else if(lep[ilep].isElectron())
        allPlots["massD0_e"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
      */
    }
    if(pfCands.size()<3) return;
    if(abs(pfCands[2].getPdgId())==13) {
      //if(mass12>1.7 && mass12<2.0) { //Plot D0 only
      if(mass12>1.8 && mass12<1.93) { //Plot D0 only
        allPlots["massD0_mu_tag_l"+chTag+name+runPeriod_]->Fill(mass123,getWgt());
        allPlots["massD0_mu_tag_l_all"+name+runPeriod_]->Fill(mass123,getWgt());
      }
    }
    if(fabs(mass12-1.864) > 0.05) return; // tighter mass window cut on D^0 for D^*
    // Kaon and pion have opposite charges
    // I.e. correct mass assumption
    if(pfCands[1].charge() == pfCands[2].charge()) return;
    TLorentzVector p_cand = pfCands[0].getVec()+pfCands[1].getVec()+pfCands[2].getVec();
    float deltam = p_cand.M() - mass12;
    if(deltam<0.14 || deltam>0.15) return; // cut on D^* mass window
    float mass1234 = (p_cand + lep[ilep].getVec()).M();
    allPlots["massDs_l"+chTag+name+runPeriod_]->Fill(mass1234,getWgt());
    allPlots["massDs_l_all"+name+runPeriod_]->Fill(mass1234,getWgt());
    if(lep[ilep].isMuon())
      allPlots["massDs_mu"+chTag+name+runPeriod_]->Fill(mass1234,getWgt());
    else if(lep[ilep].isElectron())
      allPlots["massDs_e"+chTag+name+runPeriod_]->Fill(mass1234,getWgt());
  }
}

//Called by Fill(std::vector<pfTrack> pfCands, Leptons lep, Jet jet, TString chTag, TString name)
void StdPlots::Fill(std::vector<pfTrack> &pfCands, Jet jet, TString chTag, TString name) {
  if(!isGood_) return;
  float wgt = norm_ * sfs_.first * puWgt_ * top_pt_wgt_;
  if(!name.EqualTo("")) name = "_" + name;
  //set unc^2
  //std::pair<float, float> unc = std::pair<float,float>(pow(sfs_.second.first,2), pow(pi_wgt_.second.first,2),
                                                       //pow(sfs_.second.second,2), pow(pi_wgt_.second.second,2));
  /*
  std::pair<float, float> unc = getUnc();
  //add unc^2 from hist and take sqrt
  unc.first = sqrt(unc.first + pow(allPlots["unc"+chTag+name+runPeriod_]->GetBinContent(1),2));
  unc.second = sqrt(unc.second + pow(allPlots["unc"+chTag+name+runPeriod_]->GetBinContent(2),2));
  allPlots["unc"+chTag+name+runPeriod_]->SetBinContent(1,unc.first);
  allPlots["unc_all"+name+runPeriod_]->SetBinContent(1,unc.first);
  allPlots["unc"+chTag+name+runPeriod_]->SetBinContent(2,unc.second);
  allPlots["unc_all"+name+runPeriod_]->SetBinContent(2,unc.second);
  float unc = sqrt(getUnc() + pow(allPlots["unc"+chTag+name+runPeriod_]->GetBinContent(1),2));
  allPlots["unc"+chTag+name+runPeriod_]->SetBinContent(1,unc);
  allPlots["unc_all"+name+runPeriod_]->SetBinContent(1,unc);
  //allPlots["unc"+chTag+name+runPeriod_]->SetBinError(1,0);
  //allPlots["unc"+chTag+name+runPeriod_]->SetBinError(2,0);
  */
  if(name.Contains("jpsi")) {
    if(debug_) std::cout << "Filling J/Psi+jet" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;
    if(debug_) std::cout << "is " << name << std::endl;
    if(pfCands[0].getPfid() != -pfCands[1].getPfid()) return;
    TLorentzVector jpsi = pfCands[0].getVec() + pfCands[1].getVec();
    //if(jpsi.M()<3.0 || jpsi.M()>3.2) return; //Window in Elvire's AN
    if(abs(jpsi.M()-3.097) > 0.1) return;
    if(debug_) std::cout << "is J/Psi" << name << std::endl;
    float jpt(jet.getPt());
    float jpt_charged(jet.getChargedPt());
    float jpt_pf(jet.getPFPt());
    allPlots["JPsioJet_pt"+chTag+name+runPeriod_]->Fill(jpsi.Pt()/jpt,getWgt());
    allPlots["JPsioJet_pt_all"+name+runPeriod_]->Fill(jpsi.Pt()/jpt,getWgt());
    allPlots["JPsioJet_pt_charged"+chTag+name+runPeriod_]->Fill(jpsi.Pt()/jpt_charged,getWgt());
    allPlots["JPsioJet_pt_charged_all"+name+runPeriod_]->Fill(jpsi.Pt()/jpt_charged,getWgt());
    allPlots["JPsioJet_pt_pf"+chTag+name+runPeriod_]->Fill(jpsi.Pt()/jpt_pf,getWgt());
    allPlots["JPsioJet_pt_pf_all"+name+runPeriod_]->Fill(jpsi.Pt()/jpt_pf,getWgt());

    std::vector<pfTrack> &tracks = jet.getTracks();
    allPlots["JPsioJet_pt_hard"+chTag+name+runPeriod_]->Fill(tracks[0].Pt()/jpt_charged,getWgt());
    allPlots["JPsioJet_pt_hard_all"+name+runPeriod_]->Fill(tracks[0].Pt()/jpt_charged,getWgt());
    allPlots["j_pt_ch"+chTag+name+runPeriod_]->Fill(jpt_charged,getWgt());
    allPlots["j_pt_ch_all"+name+runPeriod_]->Fill(jpt_charged,getWgt());
    int mu=-1;
    for(size_t itk=0;itk<tracks.size();itk++) {
      if(abs(tracks[itk].getPdgId())==13) {
        mu = itk;
        break;
      }
    }
    if(mu>=0) {
      allPlots["JPsioJet_pt_mu"+chTag+name+runPeriod_]->Fill(tracks[mu].Pt()/jpt_charged,getWgt());
      allPlots["JPsioJet_pt_mu_all"+name+runPeriod_]->Fill(tracks[mu].Pt()/jpt_charged,getWgt());
    }
    /*
    float csv = jet.getCSV();
    allPlots["csv"+chTag+name+runPeriod_]->Fill(csv,getWgt());
    allPlots["csv_all"+name+runPeriod_]->Fill(csv,getWgt());
    */
  }
  else if(name.Contains("meson")) {
    if(debug_) std::cout << "Filling D_meson+jet" << chTag << name << runPeriod_ << " with wgt=" << wgt << std::endl;
    if(pfCands.size()<2) return;
    //recheck first pi and K
    if(pfCands[0].getPdgId()*pfCands[1].getPdgId() != -211*211) return;
    /*
    if(pfCands[0].getMotherId()!=421) return;
    if(pfCands[1].getMotherId()!=421) return;
    */
    TLorentzVector D0 = pfCands[0].getVec() + pfCands[1].getVec();
    float mass12 = D0.M();
    float jpt(jet.getPt());
    float jpt_charged(jet.getChargedPt());
    float jpt_pf(jet.getPFPt());
    if (mass12<1.7 && mass12>2.0) return;
    if(pfCands.size()==2) { //Plot D0 only
      if(abs(mass12-1.864) > 0.034) return;  // tighter mass window cut
      if(debug_) std::cout << "Filling D0" << chTag << name << runPeriod_ << " with wgt=" << getWgt() << std::endl;
      allPlots["nkj"+chTag+name+runPeriod_]->Fill(1,getWgt());
      allPlots["D0oJet_pt"+chTag+name+runPeriod_]->Fill(D0.Pt()/jpt,getWgt());
      allPlots["D0oJet_pt_all"+name+runPeriod_]->Fill(D0.Pt()/jpt,getWgt());
      allPlots["D0oJet_pt_charged"+chTag+name+runPeriod_]->Fill(D0.Pt()/jpt_charged,getWgt());
      allPlots["D0oJet_pt_charged_all"+name+runPeriod_]->Fill(D0.Pt()/jpt_charged,getWgt());
      allPlots["D0oJet_pt_pf"+chTag+name+runPeriod_]->Fill(D0.Pt()/jpt_pf,getWgt());
      allPlots["D0oJet_pt_pf_all"+name+runPeriod_]->Fill(D0.Pt()/jpt_pf,getWgt());

      std::vector<pfTrack> &tracks = jet.getTracks();
      allPlots["D0oJet_pt_hard"+chTag+name+runPeriod_]->Fill(tracks[0].Pt()/jpt_charged,getWgt());
      allPlots["D0oJet_pt_hard_all"+name+runPeriod_]->Fill(tracks[0].Pt()/jpt_charged,getWgt());
      allPlots["njtk"+chTag+name+runPeriod_]->Fill(tracks.size(),getWgt());
      allPlots["njtk_all"+name+runPeriod_]->Fill(tracks.size(),getWgt());
      allPlots["j_tk_pt"+chTag+name+runPeriod_]->Fill(jpt,getWgt());
      allPlots["j_tk_pt_all"+name+runPeriod_]->Fill(jpt,getWgt());
      allPlots["j_pt_ch"+chTag+name+runPeriod_]->Fill(jpt_charged,getWgt());
      allPlots["j_pt_ch_all"+name+runPeriod_]->Fill(jpt_charged,getWgt());
      int mu=-1;
      for(size_t itk=0;itk<tracks.size();itk++) {
        if(abs(tracks[itk].getPdgId())==13) {
          mu = itk;
          break;
        }
      }
      if(mu>=0) {
        allPlots["D0oJet_pt_mu"+chTag+name+runPeriod_]->Fill(tracks[mu].Pt()/jpt_charged,getWgt());
        allPlots["D0oJet_pt_mu_all"+name+runPeriod_]->Fill(tracks[mu].Pt()/jpt_charged,getWgt());
      }

      float cosDjet = (D0.Vect()).Dot(jet.getVec().Vect())/(jet.getVec().P());
      allPlots["D0dotJet"+chTag+name+runPeriod_]->Fill(cosDjet,getWgt());
      allPlots["D0dotJet_all"+name+runPeriod_]->Fill(cosDjet,getWgt());
      cosDjet /= D0.P();
      allPlots["D0cosJet"+chTag+name+runPeriod_]->Fill(cosDjet,getWgt());
      allPlots["D0cosJet_all"+name+runPeriod_]->Fill(cosDjet,getWgt());
    }
    if(pfCands.size()<3) return; // D^0+mu and D^*

    //if(abs(pfCands[2].getPdgId()) == 13 && abs(mass12-1.864) > 0.05) {  // Plot D0+mu tighter mass window cut
    if(fabs(pfCands[2].getPdgId()) == 13 && abs(mass12-1.864)<0.036) {
      TLorentzVector D0mu = D0 + pfCands[2].getVec();
      allPlots["D0_mu_tag_oJet_pt"+chTag+name+runPeriod_]->Fill(D0.Pt()/jpt,getWgt());
      allPlots["D0_mu_tag_oJet_pt_all"+name+runPeriod_]->Fill(D0.Pt()/jpt,getWgt());
      allPlots["D0_mu_tag_oJet_pt_charged"+chTag+name+runPeriod_]->Fill(D0.Pt()/jpt_charged,getWgt());
      allPlots["D0_mu_tag_oJet_pt_charged_all"+name+runPeriod_]->Fill(D0.Pt()/jpt_charged,getWgt());
      allPlots["D0_mu_tag_mu_oJet_pt"+chTag+name+runPeriod_]->Fill(D0mu.Pt()/jpt_charged,getWgt());
      allPlots["D0_mu_tag_mu_oJet_pt_all"+name+runPeriod_]->Fill(D0mu.Pt()/jpt_charged,getWgt());
      allPlots["D0_mu_tag_mu_oJet_pt_charged"+chTag+name+runPeriod_]->Fill(D0mu.Pt()/jpt_charged,getWgt());
      allPlots["D0_mu_tag_mu_oJet_pt_charged_all"+name+runPeriod_]->Fill(D0mu.Pt()/jpt_charged,getWgt());
      allPlots["D0_mu_tag_oJet_pt_pf"+chTag+name+runPeriod_]->Fill(D0.Pt()/jpt_pf,getWgt());
      allPlots["D0_mu_tag_oJet_pt_pf_all"+name+runPeriod_]->Fill(D0.Pt()/jpt_pf,getWgt());
    }

    if(fabs(pfCands[2].getPdgId()) != 211) return; //reinforce pion
    // Kaon and pion have opposite charges
    // I.e. correct mass assumption
    if(pfCands[1].charge() != -pfCands[2].charge()) return;

    if(debug_) std::cout << "Filling D*" << chTag << name << runPeriod_ << " with wgt=" << getWgt() << std::endl;
    TLorentzVector p_cand = pfCands[0].getVec()+pfCands[1].getVec()+pfCands[2].getVec();
    float deltam = p_cand.M() - mass12;
    if(debug_) std::cout << "Filling D*-D0" << chTag << name << runPeriod_ << " with wgt=" << getWgt() << std::endl;
    if(fabs(mass12-1.864) > 0.05) return; // tighter mass window cut
    if(deltam<0.14 || deltam>0.15) return;
    allPlots["D0oJet_pt_fromDs"+chTag+name+runPeriod_]->Fill(D0.Pt()/jpt,getWgt());
    allPlots["D0oJet_pt_fromDs_all"+name+runPeriod_]->Fill(D0.Pt()/jpt,getWgt());
    allPlots["D0oJet_pt_charged_fromDs"+chTag+name+runPeriod_]->Fill(D0.Pt()/jpt_charged,getWgt());
    allPlots["D0oJet_pt_charged_fromDs_all"+name+runPeriod_]->Fill(D0.Pt()/jpt_charged,getWgt());
    allPlots["D0oJet_pt_pf_fromDs"+chTag+name+runPeriod_]->Fill(D0.Pt()/jpt_pf,getWgt());
    allPlots["D0oJet_pt_pf_fromDs_all"+name+runPeriod_]->Fill(D0.Pt()/jpt_pf,getWgt());

    std::vector<pfTrack> &tracks = jet.getTracks();
    allPlots["D0oJet_pt_hard_fromDs"+chTag+name+runPeriod_]->Fill(tracks[0].Pt()/jpt_charged,getWgt());
    allPlots["D0oJet_pt_hard_fromDs_all"+name+runPeriod_]->Fill(tracks[0].Pt()/jpt_charged,getWgt());
    int mu=-1;
    for(size_t itk=0;itk<tracks.size();itk++) {
      if(abs(tracks[itk].getPdgId())==13) {
        mu = itk;
        break;
      }
    }
    if(mu>=0) {
      allPlots["D0oJet_pt_mu_fromDs"+chTag+name+runPeriod_]->Fill(tracks[mu].Pt()/jpt_charged,getWgt());
      allPlots["D0oJet_pt_mu_fromDs_all"+name+runPeriod_]->Fill(tracks[mu].Pt()/jpt_charged,getWgt());
    }
    //D*
    allPlots["DsoJet_pt"+chTag+name+runPeriod_]->Fill(p_cand.Pt()/jpt,getWgt());
    allPlots["DsoJet_pt_all"+name+runPeriod_]->Fill(p_cand.Pt()/jpt,getWgt());
    allPlots["DsoJet_pt_charged"+chTag+name+runPeriod_]->Fill(p_cand.Pt()/jpt_charged,getWgt());
    allPlots["DsoJet_pt_charged_all"+name+runPeriod_]->Fill(p_cand.Pt()/jpt_charged,getWgt());
    allPlots["DsoJet_pt_pf"+chTag+name+runPeriod_]->Fill(p_cand.Pt()/jpt_pf,getWgt());
    allPlots["DsoJet_pt_pf_all"+name+runPeriod_]->Fill(p_cand.Pt()/jpt_pf,getWgt());

    //std::vector<pfTrack> &tracks = jet.getTracks();
    allPlots["DsoJet_pt_hard"+chTag+name+runPeriod_]->Fill(tracks[0].Pt()/jpt_charged,getWgt());
    allPlots["DsoJet_pt_hard_all"+name+runPeriod_]->Fill(tracks[0].Pt()/jpt_charged,getWgt());
    mu=-1;
    for(size_t itk=0;itk<tracks.size();itk++) {
      if(abs(tracks[itk].getPdgId())==13) {
        mu = itk;
        break;
      }
    }
    if(mu>=0) {
      allPlots["DsoJet_pt_mu"+chTag+name+runPeriod_]->Fill(tracks[mu].Pt()/jpt_charged,getWgt());
      allPlots["DsoJet_pt_mu_all"+name+runPeriod_]->Fill(tracks[mu].Pt()/jpt_charged,getWgt());
    }
  }
}

//Fill for mesons (currently only J/Psi
void StdPlots::Fill(std::vector<pfTrack> &pfCands, Leptons lep, Jet jet, TString chTag, TString name) {
  if(!isGood_) return;
  Fill(pfCands, lep, chTag, name); //Fill meson+lep plots
  Fill(pfCands, jet, chTag, name); //Fill meson+jet run2
  if(!name.EqualTo("")) name = "_" + name;
  if(name.Contains("jpsi")) {
    TLorentzVector jpsi = pfCands[0].getVec() + pfCands[1].getVec();
    //if(jpsi.M()<3.0 || jpsi.M()>3.2) return; //Window in Elvire's AN
    if(abs(jpsi.M()-3.097) > 0.1) return;
    float jpt(jet.getPt());
    float jpt_charged(jet.getChargedPt());
    float jpt_pf(jet.getPFPt());
    float dRJPsil(jpsi.DeltaR(lep[0].getVec()));
    if(dRJPsil < 2.0) {
      allPlots["JPsioJet_pt_low_dR"+chTag+name+runPeriod_]->Fill(jpsi.Pt()/jpt,getWgt());
      allPlots["JPsioJet_pt_low_dR_all"+name+runPeriod_]->Fill(jpsi.Pt()/jpt,getWgt());
      allPlots["JPsioJet_pt_charged_low_dR"+chTag+name+runPeriod_]->Fill(jpsi.Pt()/jpt_charged,getWgt());
      allPlots["JPsioJet_pt_charged_low_dR_all"+name+runPeriod_]->Fill(jpsi.Pt()/jpt_charged,getWgt());
      allPlots["JPsioJet_pt_pf_low_dR"+chTag+name+runPeriod_]->Fill(jpsi.Pt()/jpt_pf,getWgt());
      allPlots["JPsioJet_pt_pf_low_dR_all"+name+runPeriod_]->Fill(jpsi.Pt()/jpt_pf,getWgt());
    }
    else {
      allPlots["JPsioJet_pt_high_dR"+chTag+name+runPeriod_]->Fill(jpsi.Pt()/jpt,getWgt());
      allPlots["JPsioJet_pt_high_dR_all"+name+runPeriod_]->Fill(jpsi.Pt()/jpt,getWgt());
      allPlots["JPsioJet_pt_charged_high_dR"+chTag+name+runPeriod_]->Fill(jpsi.Pt()/jpt_charged,getWgt());
      allPlots["JPsioJet_pt_charged_high_dR_all"+name+runPeriod_]->Fill(jpsi.Pt()/jpt_charged,getWgt());
      allPlots["JPsioJet_pt_pf_high_dR"+chTag+name+runPeriod_]->Fill(jpsi.Pt()/jpt_pf,getWgt());
      allPlots["JPsioJet_pt_pf_high_dR_all"+name+runPeriod_]->Fill(jpsi.Pt()/jpt_pf,getWgt());
    }
  }
}

StdPlots Combine(const StdPlots &plots1, const StdPlots &plots2) {
  std::map<TString, float> lumi;
  lumi["all"] = 35740.161;
  lumi["BCDEF"] = 19593.811;
  lumi["GH"] = 16146.17;
  //Load run period from plots
  //e.g. "_GH"
  TString runPeriod(plots2.runPeriod_);
  //Strip the '_'
  //e.g. "GH"
  runPeriod.Replace(0,1,"");
  //if(debug_) std::cout << "Adding plots from runPeriod " << runPeriod << std::endl;
  //Append to current run period
  //e.g. "_BCDEF"->"_BCDEFGH"
  runPeriod = plots1.runPeriod_ + runPeriod;
  runPeriod.Replace(0,1,"");
  if(plots1.debug_) std::cout << "New runPeriod " << runPeriod << std::endl;
  //Add plot to new run
  //e.g. CSV_m_BCDEFGH->Add(CSV_m_BCDEFGH)
  StdPlots run(runPeriod, plots1.name_, plots1.debug_);
  for(auto& it : plots1.allPlots) {
    if(!plots1.isGood_) break;
    //Strip run period from plot being added and append new run period
    if(it.second->Integral() == 0) continue;
    TString plotName(it.first);
    plotName.Resize((it.first.Length() - plots1.runPeriod_.Length()));
    TString tmp(plotName + plots2.runPeriod_);
    plotName += "_"+runPeriod;
    TString tmpl(plots1.runPeriod_);
    tmpl.Replace(0,1,"");
    //1 - wgt due to 0 < Integral < 1
    //e.g. BCDEF has larger lumi, should get larger weight
    //float wgt = it.second->Integral()/(it.second->Integral() + plots2.allPlots.at(tmp)->Integral());
    float wgt = lumi[tmpl]/lumi["all"];
    if(plots2.debug_) std::cout << tmp << " " << wgt << std::endl;
    if(run.name_.Contains("Data")) wgt=1.0;
    it.second->Scale(wgt);
    run.allPlots[plotName]->Add(it.second);
  }
  //Add plot from plots to new
  //e.g. CSV_m_GH->Add(CSV_m_BCDEFGH)
  for(auto& it : plots2.allPlots) {
    if(!plots2.isGood_) break;
    //Strip run period from plot being added and append new run period
    if(it.second->Integral() == 0) continue;
    TString plotName(it.first);
    plotName.Resize((it.first.Length() - plots2.runPeriod_.Length()));
    TString tmp(plotName + plots1.runPeriod_);
    plotName += "_"+runPeriod;
    TString tmpl(plots2.runPeriod_);
    tmpl.Replace(0,1,"");
    //float wgt = it.second->Integral()/(it.second->Integral() + plots1.allPlots.at(tmp)->Integral());
    float wgt = lumi[tmpl]/lumi["all"];
    if(plots1.debug_) std::cout << tmp << " " << wgt << std::endl;
    if(run.name_.Contains("Data")) wgt=1.0;
    if(plots1.debug_) std::cout << it.second->Integral() << std::endl;
    it.second->Scale(wgt);
    run.allPlots[plotName]->Add(it.second);
  }
  return run;
}

//StdPlots operator+(const StdPlots &lhs, const StdPlots &rhs) { return Combine(lhs,rhs); };

void StdPlots::Write() {
  if(!isGood_) return;
  if(debug_) std::cout << "writing histograms" << std::endl;

  for (auto& it : allPlots)  { 
    //if(debug_) std::cout << it.second->GetName() << std::endl;
    //if(debug_) std::cout << it.second->GetEntries() << std::endl;

    it.second->Write(); 
  }

  /*
  TTree *t = new TTree("top_pt_wgt","Top pT");
  t->Branch("top_pt_wgt",&top_pt_wgt_vec);
  t->Fill();
  t->Write();
  */
}

