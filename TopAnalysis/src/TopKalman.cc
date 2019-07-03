#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CharmEvent.h"
#include "TopLJets2015/TopAnalysis/interface/FragEvent.h"
#include "TopLJets2015/TopAnalysis/interface/TOPWidth.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"
#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"

//#include "TopLJets2015/TopAnalysis/interface/OtherFunctions.h"
#include "TopLJets2015/TopAnalysis/interface/Trigger.h"
#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"
#include "TopLJets2015/TopAnalysis/interface/Jet.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h" //FIXME
#include "TopLJets2015/TopAnalysis/interface/StdPlots.h"
#include "TopLJets2015/TopAnalysis/interface/CharmTree.h"
#include "TopLJets2015/TopAnalysis/interface/KalmanEvent.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <map>
#include <unordered_map>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

using namespace std;

enum systBit { TRIGGER_BIT=1, LEP_BIT, PU_BIT, PI_BIT, JER_BIT, JSF_BIT };
//enum systBit { TRIGGER_BIT=1, LEP_BIT, TRK_BIT, PU_BIT, PI_BIT, JER_BIT };
int passBit(int syst, int BIT) {
  if(syst==0) return 0;
  return (syst/abs(syst) * ((abs(syst)>>BIT)&0x1));
}

float customSF(pfTrack &pi, TString runPeriod, float &ptsf, float &etasf, TH2 *pf_eff_ = nullptr) {
  float sf(1.); //normalization derived from GH epoch, may account for MC track discrepancy
  if(pf_eff_ != nullptr) {
  float minEtaForEff( pf_eff_->GetXaxis()->GetXmin() ), maxEtaForEff( pf_eff_->GetXaxis()->GetXmax()-0.01 );
  float etaForEff=TMath::Max(TMath::Min(pi.Eta(),maxEtaForEff),minEtaForEff);
  Int_t etaBinForEff=pf_eff_->GetXaxis()->FindBin(etaForEff);
  
  float minPtForEff( pf_eff_->GetYaxis()->GetXmin() ), maxPtForEff( pf_eff_->GetYaxis()->GetXmax()-0.01 );
  float ptForEff=TMath::Max(TMath::Min(pi.Pt(),maxPtForEff),minPtForEff);
  //float ptForEff=TMath::Max(TMath::Min(TMath::Max(pi.Pt(),5.),maxPtForEff),minPtForEff);
  Int_t ptBinForEff=pf_eff_->GetYaxis()->FindBin(ptForEff);

  sf=pf_eff_->GetBinContent(etaBinForEff,ptBinForEff);
  pi.setPtCorrection(sf);
  }
  //pi.setEtaCorrection(1.11);
  return sf;
  /*
  // pol1
  if(runPeriod=="BCDEF") {
    if(pi.Eta() < -1.1)
      sf *= 1.33e+00 + pi.Eta() * 3.35e-01;
    else if(pi.Eta() > -1.1 && pi.Eta() < -0.8)
      sf *= 5.490e-01 + pi.Eta() * -3.934e-01;
    else if(pi.Eta() > -0.8 && pi.Eta() < -0.4)
      sf *= 8.12e-01 + pi.Eta() * -8.78e-02;
    else if(pi.Eta() > -0.4 && pi.Eta() < 0.0)
      sf *= 8.47e-01 + pi.Eta() * -1.37e-01;
    else if(pi.Eta() > 0.0 && pi.Eta() < 0.5)
      sf *= 8.30e-01 + pi.Eta() * 1.92e-01;
    else if(pi.Eta() > 0.5 && pi.Eta() < 0.8)
      sf *= 6.82e-01 + pi.Eta() * 2.87e-01;
    else if(pi.Eta() > 0.8 && pi.Eta() < 1.1)
      sf *= 8.50e-01 + pi.Eta() * 6.62e-02;
    else
      sf *= 1.02e+00 + pi.Eta() * -6.01e-02;
  }
  else if(runPeriod=="GH") {
    sf = 1.11; //const normalization
  }
  if(etasf > 0) etasf = sf;
  pi.setEtaCorrection(sf);
  if(runPeriod=="BCDEF") {
    float ptcor = 8.65e-01 + pi.Pt() * 7.63e-03;
  ptcor = 1.; //FIXME
    sf *= ptcor;
    if(ptsf > 0) ptsf = ptcor;
    if(ptsf > 0) pi.setPtCorrection(ptcor);
  }
  //return 1.0; //FIXME
  float v0(1.),v1(1.),t(0.),pt(pi.Pt());
  if(runPeriod=="BCDEF") {
    //if(pt<25) return 0.8199 + 0.0108 * pt;
    //else return 0.9242 + 0.00463 * pt;
    //tails
    v0 = 0.8199 + 0.0108 * pt;
    v1 = 0.9212 + 0.00465 * pt;
  
  }
  else if(runPeriod=="GH") {
    //if(pt<25) return 1.070 + 0.00254 * pt;
    //else return 1.073 + 0.00197 * pt;
    //tails
    v0 = 1.070 + 0.00256 * pt;
    v1 = 1.070 + 0.00209 * pt;
  }
  //else return 1.;
  if(pt>25) t = 1.;
  else if(pt>24)
    //map 24-25 to 0-1
    t = 0 + (1 - 0) * ((pt - 24) / (25 - 24));
  //lerp
  return (1 - t) * v0 + t * v1;
  */
}
float customSF(pfTrack &pi, TString runPeriod) {
  float ptsf(0),etasf(0);
  return customSF(pi, runPeriod, ptsf, etasf); 
}

void RunTopKalman(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
                 TString era,
                 TString runPeriod,
                 Bool_t debug=false,
                 short runSysts, //0 for nominal, 1 for up, -1 for down
                 short rbFit)
{
  if(debug) cout << "in RunTopKalman" << endl;

  bool isTTbar( filename.Contains("_TTJets") );
  bool isData( filename.Contains("Data13TeV") );

  //CREATE CHARM TREE IN FILE
  TTree *cht = new TTree("data","Charm tree");
  cht->SetDirectory(0);
  CharmEvent_t evch;
  createCharmEventTree(cht,evch);

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << (runSysts ? " syst variations will be considered" : "") << endl;

  TString updown = runSysts > 0 ? "up" : "down";
  if(abs(passBit(runSysts,TRIGGER_BIT))) std::cout << TString::Format("running %s trigger systematics", updown.Data()) << std::endl;
  if(abs(passBit(runSysts,LEP_BIT))) std::cout << TString::Format("running %s lepton selection systematics", updown.Data()) << std::endl;
  //if(abs(passBit(runSysts,TRK_BIT))) std::cout << TString::Format("running tracker %s efficiency systematics", updown.Data()) << std::endl;
  if(abs(passBit(runSysts,PU_BIT))) std::cout << TString::Format("running %s PU systematics", updown.Data()) << std::endl;
  if(abs(passBit(runSysts,PI_BIT))) std::cout << TString::Format("running %s pi systematics", updown.Data()) << std::endl;
  if(abs(passBit(runSysts,JER_BIT))) std::cout << TString::Format("running %s JER systematics", updown.Data()) << std::endl;
  if(abs(passBit(runSysts,JSF_BIT))) std::cout << TString::Format("running %s JSF systematics", updown.Data()) << std::endl;
  
  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  //std::vector<TGraph *>puWgt;
  std::vector<TH1D *>puWgt;
  TString tmpRun ("BCDEFGH");
  std::vector<TH1D*> puWgtsRun;
  if(!isData)
    {
      if(debug) cout << "loading pileup weight" << endl;
      /*
      TString puWgtUrl(era+"/pileupWgts"+runPeriod+".root");
      gSystem->ExpandPathName(puWgtUrl);
      TFile *fIn=TFile::Open(puWgtUrl);
      for(size_t i=0; i<3; i++)
	{
	  TString grName("pu_nom");
          if(i==1) grName="pu_down";
          if(i==2) grName="pu_up";
          TH1D *puWgtData=(TH1D *)fIn->Get("puwgts_nom");
          puWgt.push_back(puWgtData);
	}
      */

      //Load PU plots for all run periods
      /*
      for(int ic = 0; ic < tmpRun.Length(); ic++) {
        TString puWgtRunUrl(era+"/pileupWgts"+tmpRun[ic]+".root");
        gSystem->ExpandPathName(puWgtRunUrl);
        fIn=TFile::Open(puWgtRunUrl);
        TH1D *puWgtDataRun=(TH1D *)fIn->Get("puwgts_nom");
        puWgtDataRun->SetDirectory(0);
	puWgtsRun.push_back(puWgtDataRun);
        fIn->Close();
        
      }
      */
      TString puWgtRunUrl(era+"/pileupWgtsBCDEF.root");
      gSystem->ExpandPathName(puWgtRunUrl);
      TFile *fIn=TFile::Open(puWgtRunUrl);
      for(size_t i = 0; i < 3; i++) {
        TString  grName("puwgts_down");
        if(i==1) grName="puwgts_nom";
        if(i==2) grName="puwgts_up";
        TH1D *puWgtDataRun=(TH1D *)fIn->Get(grName);
        puWgtDataRun->SetDirectory(0);
        puWgtsRun.push_back(puWgtDataRun);
      }
      fIn->Close();
      puWgtRunUrl = era+"/pileupWgtsGH.root";
      gSystem->ExpandPathName(puWgtRunUrl);
      fIn=TFile::Open(puWgtRunUrl);
      /*
      puWgtDataRun=(TH1D *)fIn->Get("puwgts_nom");
      puWgtDataRun->SetDirectory(0);
      puWgtsRun.push_back(puWgtDataRun);
      */
      for(size_t i = 0; i < 3; i++) { // [BCDEFup, BCDEF, BCDEFdown, GHup, GH, GHdown]
        TString  grName("puwgts_down");
        if(i==1) grName="puwgts_nom";
        if(i==2) grName="puwgts_up";
        TH1D *puWgtDataRun=(TH1D *)fIn->Get(grName);
        puWgtDataRun->SetDirectory(0);
        puWgtsRun.push_back(puWgtDataRun);
      }
      fIn->Close();
      if(debug) cout << "loading pileup weight DONE" << endl;
    }

  //r_B EVENT WEIGHT
  TString rbWgtUrl(era+"/bfragweights.root");
  gSystem->ExpandPathName(rbWgtUrl);
  TFile *fIn=TFile::Open(rbWgtUrl);
  TGraph *rbWgt=(TGraph*)fIn->Get("fitFrag");
  TGraph *rbWgts[3];
  rbWgts[JPsi]=(TGraph*)fIn->Get("jpsiFrag");
  rbWgts[D0]=(TGraph*)fIn->Get("d0Frag");
  rbWgts[D0mu]=(TGraph*)fIn->Get("d0muFrag");
  fIn->Close();

  //LEPTON EFFICIENCIES
  //LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era,runPeriod,debug);
  LeptonEfficiencyWrapper lepEffH_BCDEF(filename.Contains("Data13TeV"),era,"BCDEF",debug);
  LeptonEfficiencyWrapper lepEffH_GH(filename.Contains("Data13TeV"),era,"GH",debug);


  //jet energy uncertainties
  TString jecUncUrl(era+"/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  JetCorrectorParameters *jecParam = new JetCorrectorParameters(jecUncUrl.Data(), "Total");
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty( *jecParam );

  //jet energy resolution
  TString jerUrl(era+"/Summer16_25nsV1_MC_EtaResolution_AK4PF.txt");
  gSystem->ExpandPathName(jerUrl);
  JME::JetResolution *jres = new JME::JetResolution(jerUrl.Data());
  TString jerSFUrl(era+"/Summer16_25nsV1_MC_SF_AK4PFchs.txt");
  gSystem->ExpandPathName(jerSFUrl);
  JME::JetResolutionScaleFactor *jres_sf = new JME::JetResolutionScaleFactor(jerSFUrl.Data());
  
  
  //PF TRACK EFFICIENCIES (B-F / GH)
  TString pfUrl(era+"/pf_tracks.root");
  gSystem->ExpandPathName(pfUrl);
  fIn=TFile::Open(pfUrl);
  TH2 *pf_eff=(TH2 *)fIn->Get("eta_pt")->Clone();
  TH2 *pf_eff_pi=(TH2 *)fIn->Get("pi_eta_pt")->Clone("eta_pt_pi");
  TH2 *pf_eff_k=(TH2 *)fIn->Get("k_eta_pt")->Clone("eta_pt_k");
  pf_eff->SetDirectory(0);
  pf_eff_pi->SetDirectory(0);
  pf_eff_k->SetDirectory(0);
  fIn->Close();
  //new blank hist for plottinf pi
  TH2F *pi_pf_eff = (TH2F*)pf_eff->Clone("pi_eta_pt");
  pfUrl=era+"/pf_tracks_jpsi.root";
  fIn=TFile::Open(pfUrl);
  gSystem->ExpandPathName(pfUrl);
  TH2 *pf_eff_jpsi=(TH2 *)fIn->Get("eta_pt")->Clone();
  pf_eff_jpsi->SetDirectory(0);
  fIn->Close();
  TString epoch("BCDEF");
  bool isGood_(false);
  //Dirty hack to enforce data epoch
  if(isData) {
    for(int i=0; i < epoch.Length(); i++) {
      TString tmp(epoch[i]);
      if(filename.Contains("MC") || (filename.Contains("Data") && filename.Contains("2016"+tmp))) {
        isGood_ = true;
        break;
      }
    }
    if(!isGood_) epoch = "GH";
  }
  pi_pf_eff->SetName("pi_eta_pt_" + epoch);
  pi_pf_eff->SetTitle("pi_eta_pt_" + epoch);
  pi_pf_eff->Reset();
  pi_pf_eff->SetDirectory(0);

  //LIST OF SYSTEMATICS
  
  //HISTOGRAMS BY RUNPERIOD
  StdPlots runBCDEF("BCDEF", outname, debug);
  StdPlots runGH("GH", outname, debug);
  
  CharmTree treeBCDEF(t, "BCDEF", outname, debug);
  CharmTree treeGH(t, "GH", outname, debug);

  std::vector<float> dpt;
  std::vector<float> deta;
  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  std::map<TString, TH2 *> allPlots2D; //hash map insted of map (hash faster than red-black tree: O(1) vs O(log(N)))
  std::map<TString, TH3 *> allPlots3D;
  //std::unordered_map<TString, TH2 *> allPlots2D; //has map insted of map (hash faster than red-black tree: O(1) vs O(log(N)))
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",4,0,4);
  allPlots["puwgtgr"] = new TH1F("puwgtgr","PU Weights (calc)",75,0,75);
  allPlots["puwgt"] = new TH1F("puwgt","PU Weights (data)",75,0,75);
  allPlots["topptwgt"] = new TH1F("topptwgt","Top p_{T} weights", 2, 0, 2);
  allPlots2D["toppteta"] = new TH2F("toppteta","Top p_{T} weights vs. #eta", 100, -2.4, 2.4, 100, 0, 2);
  allPlots["jerwgt"] = new TH1F("jerwgt","JER weights", 2, 0, 2);
  allPlots["jsfwgt"] = new TH1F("jsfwgt","JSF weights", 2, 0, 2);
  allPlots["piwgt_BCDEF"] = new TH1F("piwgt_BCDEF","pi weights", 2, 0, 2);
  allPlots["piwgt_GH"] = new TH1F("piwgt_GH","pi weights", 2, 0, 2);
  allPlots2D["jereta"] = new TH2F("jereta","JER weights vs. #eta", 100, -2.4, 2.4, 200, 0, 2);
  allPlots3D["jeceta"] = new TH3F("jeceta","JEC weights (uncorrectd p_{T} / corrected p_{T}) vs. p{T} vs. #eta", 500, 0, 500, 100, -2.4, 2.4, 200, 0, 2);
  std::vector<TString> lfsVec = { "_all", "_e", "_ee", "_em", "_mm", "_m" }; 
  std::vector<TString> cutVec = { "", "_lep", "_lepjets", "_jpsi", "_csv", "_meson" };
  std::vector<TString> wgtVec = { "", "_no_weight" };

  for(int i = 0; i < (int)lfsVec.size(); i++) {
    TString tag(lfsVec[i]);
    allPlots["pid"+tag] = new TH1F("pid"+tag,";pid;Events triggered",3,0,3);
  for(int j = 0; j < (int)cutVec.size(); j++) {
  for(int k = 0; k < (int)wgtVec.size(); k++) {
    TString cut(cutVec[j]);
    TString weight(wgtVec[k]);
    allPlots["lp_pt_iso"+tag+cut+weight] = new TH1F("lp_pt_iso"+tag+cut+weight,";Lepton P_{T} [GeV] after cleaning;Events / 10 GeV", 20, 0,200);
    allPlots["lp_pt_veto"+tag+cut+weight] = new TH1F("lp_pt_veto"+tag+cut+weight,";Lepton P_{T} [GeV] after veto;Events / 10 GeV", 20, 0,200);
    allPlots["lp_pt_low"+tag+cut+weight] = new TH1F("lp_pt_low"+tag+cut+weight,";Leading Lepton P_{T} [GeV];Events / 1 GeV", 20, 20,40);
    allPlots["lp_pt"+tag+cut+weight] = new TH1F("lp_pt"+tag+cut+weight,";Leading lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["l2p_pt"+tag+cut+weight] = new TH1F("l2p_pt"+tag+cut+weight,";Sub-leading lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_pt"+tag+cut+weight] = new TH1F("dilp_pt"+tag+cut+weight,";Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_m"+tag+cut+weight] = new TH1F("dilp_m"+tag+cut+weight,";M_{ll} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["j_pt"+tag+cut+weight] = new TH1F("j_pt"+tag+cut+weight,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["lj_pt"+tag+cut+weight] = new TH1F("lj_pt"+tag+cut+weight,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["kj_pt"+tag+cut+weight] = new TH1F("kj_pt"+tag+cut+weight,";Leading b Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["j_pt_low"+tag+cut+weight] = new TH1F("j_pt_low"+tag+cut+weight,";Leading light Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["lj_pt_low"+tag+cut+weight] = new TH1F("lj_pt_low"+tag+cut+weight,";Leading light Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["kj_pt_low"+tag+cut+weight] = new TH1F("kj_pt_low"+tag+cut+weight,";Leading b Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["nlp"+tag+cut+weight]     = new TH1F("nlp"+tag+cut+weight,";N_{l};Events" ,3,0.,3.);
    allPlots["ndilp"+tag+cut+weight]     = new TH1F("ndilp"+tag+cut+weight,";N_{ll};Events" ,3,0.,3.);
    allPlots["nj"+tag+cut+weight]     = new TH1F("nj"+tag+cut+weight,";N_{jets} (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["nlj"+tag+cut+weight]     = new TH1F("nlj"+tag+cut+weight,";N_{jets} (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["nkj"+tag+cut+weight]     = new TH1F("nkj"+tag+cut+weight,";N_{b-jets} (Klaman jet);Events" ,4,1.,5.);
    allPlots["nbj"+tag+cut+weight]     = new TH1F("nbj"+tag+cut+weight,";N_{b-jets};Events" ,4,1.,5.);
    allPlots["npf"+tag+cut+weight]     = new TH1F("npf"+tag+cut+weight,";N_{pf};Events / 10" ,5,0.,5.);
    allPlots["lp_eta"+tag+cut+weight]  = new TH1F("lp_eta"+tag+cut+weight,";Leading lepton #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["l2p_eta"+tag+cut+weight]  = new TH1F("l2p_eta"+tag+cut+weight,";Sub-Leading lepton #eta; Events / 0.1", 30, -2.5,2.5);
    allPlots["lp_phi"+tag+cut+weight]  = new TH1F("lp_phi"+tag+cut+weight,";Leading lepton #phi; Events", 50, -3.14,3.14);
    allPlots["l2p_phi"+tag+cut+weight]  = new TH1F("l2p_phi"+tag+cut+weight,";Sub-Leading lepton #phi; Events", 50, -3.14,3.14);
    allPlots["nstart"+tag+cut+weight]     = new TH1F("jetindex"+tag+cut+weight,";N_{jetindex};Events" ,5,0.,5.);
    allPlots["pfid"+tag+cut+weight]     = new TH1F("pfid"+tag+cut+weight,";PFID;Events" ,440,-220.,220.);
    allPlots["massJPsi"+tag+cut+weight]     = new TH1F("massJPsi"+tag+cut+weight,";M_{ll};Events / 18 MeV" ,50,2.5,3.4);
    allPlots["massJPsiK"+tag+cut+weight]     = new TH1F("massJPsiK"+tag+cut+weight,";M_{llk};Events / 15 MeV" ,100,4.5,6);
    allPlots["mt"+tag+cut+weight]         = new TH1F("mt"+tag+cut+weight,";m_{T};Events / 1 GeV" ,200,0,200);
    allPlots["massD0"+tag+cut+weight]     = new TH1F("massD0"+tag+cut+weight,";M_{D^{0}};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_1j"+tag+cut+weight]     = new TH1F("massD0_1j"+tag+cut+weight,";M_{D^{0}} from 1 light jet;Events / 5 MeV" ,60,1.7,2.0);
    allPlots["massD0_2j"+tag+cut+weight]     = new TH1F("massD0_2j"+tag+cut+weight,";M_{D^{0}} from 2 light jets;Events / 5 MeV" ,60,1.7,2.0);
    allPlots["massD0_3j"+tag+cut+weight]     = new TH1F("massD0_3j"+tag+cut+weight,";M_{D^{0}} from 3 light jets;Events / 5 MeV" ,60,1.7,2.0);
    allPlots["massD0_4j"+tag+cut+weight]     = new TH1F("massD0_4j"+tag+cut+weight,";M_{D^{0}} from 4 light jets;Events / 5 MeV" ,60,1.7,2.0);
    allPlots["massD0_lep"+tag+cut+weight]     = new TH1F("massD0_lep"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_mu"+tag+cut+weight]     = new TH1F("massD0_mu"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_e"+tag+cut+weight]     = new TH1F("massD0_ele"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massDsmD0loose"+tag+cut+weight]     = new TH1F("massDsmD0loose"+tag+cut+weight,";M_{K#pi#pi} - M_{K#pi};Events / 0.5 MeV" ,20,0.14,0.16);
    allPlots["massDsmD0"+tag+cut+weight]     = new TH1F("massDsmD0"+tag+cut+weight,";M_{K#pi#pi} - M_{K#pi};Events / 0.5 MeV" ,20,0.14,0.16);
    allPlots["massDs"+tag+cut+weight]     = new TH1F("massDs"+tag+cut+weight,";M_{D*};Events / 6 MeV" ,200,1.6,2.2);
    allPlots["pi_pt"+tag+cut+weight] = new TH1F("pi_pt"+tag+cut+weight,";#pi^{#pm} P_{T} [GeV];Events / 5 GeV", 10, 0,50);
    allPlots["MET"+tag+cut+weight] = new TH1F("MET"+tag+cut+weight,";MET [GeV];Events / 20 GeV", 10,0,200);
    allPlots["HT"+tag+cut+weight] = new TH1F("HT"+tag+cut+weight,";HT [GeV];Events / 20 GeV", 55,0,1100);
    allPlots["H"+tag+cut+weight] = new TH1F("H"+tag+cut+weight,";H [GeV];Events / 20 GeV", 55,0,1100);
    allPlots["HTb"+tag+cut+weight] = new TH1F("HTb"+tag+cut+weight,";HT (b-tagged) [GeV];Events / 20 GeV", 55,0,1100);
    allPlots["ST"+tag+cut+weight] = new TH1F("ST"+tag+cut+weight,";ST [GeV];Events / 20 GeV", 10,0,200);
    allPlots["MET2oST"+tag+cut+weight] = new TH1F("MET2oST"+tag+cut+weight,";MET2oST [GeV];Events / 20 GeV", 10,0,200);
    allPlots["charge"+tag+cut+weight] = new TH1F("charge"+tag+cut+weight,";Charge(l_{1}*l_{2});Events", 5,-2,2);
    allPlots["dR"+tag+cut+weight] = new TH1F("dR"+tag+cut+weight,";dR;Events / 0.05", 20,0.0,1.);
    allPlots["pflp_pt"+tag+cut+weight] = new TH1F("pflp_pt"+tag+cut+weight,";PF lepton P_{T} [GeV];Events / 0.2 GeV", 15, 0,3);
    allPlots["massZ"+tag+cut+weight]     = new TH1F("massZ_control"+tag+cut+weight,";M_{ll};Events / 1.0 GeV" ,30,81,111);
    allPlots["chargeZ"+tag+cut+weight]     = new TH1F("chargeZ_control"+tag+cut+weight,";M_{ll};Events / 1.0 GeV" ,5,-2,2);
    allPlots["nevt"+tag+cut+weight]     = new TH1F("nevt"+tag+cut+weight,";N_{events};Events" ,1,1.,2.);
    allPlots["weight"+tag+cut+weight]     = new TH1F("weight"+tag+cut+weight,";N_{events};Events/ 1.0" ,20,0.,2.);
    allPlots["norm"+tag+cut+weight]     = new TH1F("norm"+tag+cut+weight,";N_{events};Events / 1.0" ,2,0.,2.);
    allPlots["relIso"+tag+cut+weight] = new TH1F("relIso"+tag+cut+weight,";relIso;Events / 0.01", 25,0,0.25);
    allPlots["nvtx"+tag+cut+weight]     = new TH1F("nvtx"+tag+cut+weight,";N_{PV};Events / 1.0" ,50,0.,50.);
    allPlots["chi2"+tag+cut+weight] = new TH1F("normchi2"+tag+cut+weight,";#chi^2/n.d.o.f.;Events", 10,0.,10.);
    allPlots["lp_dxy"+tag+cut+weight] = new TH1F("lp_dxy"+tag+cut+weight,";d_{xy} [cm];Events / 0.01 #mum", 20, 0, 0.2);
    allPlots["lp_dz"+tag+cut+weight] = new TH1F("lp_dz"+tag+cut+weight,";d_{z} [cm];Events / 0.01 #mum", 50, 0, 0.5);
    allPlots["pf_dxy"+tag+cut+weight] = new TH1F("pf_dxy"+tag+cut+weight,";d_{xy} [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dz"+tag+cut+weight] = new TH1F("pf_dz"+tag+cut+weight,";d_{z} [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dxyE"+tag+cut+weight] = new TH1F("pf_dxyE"+tag+cut+weight,";#sigma(d_{xy}) [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dzE"+tag+cut+weight] = new TH1F("pf_dzE"+tag+cut+weight,";#sigma(d_{z}) [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dxy_sig"+tag+cut+weight] = new TH1F("pf_dxy_significance"+tag+cut+weight,";d_{xy};Events / 1", 30, 0, 30);
    allPlots["pf_dz_sig"+tag+cut+weight] = new TH1F("pf_dz_significance"+tag+cut+weight,";d_{z};Events / 1", 30, 0, 30);

    //gen-level plots
    allPlots["gtop_pt"+tag+cut+weight+runPeriod] = new TH1F("gtop_pt"+tag+cut+weight+runPeriod,";Generator top P_{T} [GeV];Events / 10 GeV", 40, 0,400);

  }
  }
  }
    allPlots["nevt_iso"] = new TH1F("nevt_iso",";After Isolation;Events", 1,1.,2.);
    allPlots["nevt_veto"] = new TH1F("nevt_veto",";After Veto;Events", 1,1.,2.);
    allPlots["nevt_trigger"] = new TH1F("nevt_trigger",";After Trigger;Events", 1,1.,2.);
    allPlots["norm_iso"] = new TH1F("norm_iso",";After Isolation;Events / 1.0", 20,0,2.);
    allPlots["norm_veto"] = new TH1F("norm_veto",";After Veto;Events / 1.0", 20,0.,2.);
    allPlots["nvtx_iso"]     = new TH1F("nvtx_iso",";N_{PV};Events / 1.0" ,50,0.,50.);
    allPlots["nvtx_veto"]     = new TH1F("nvtx_veto",";N_{PV};Events / 1.0" ,50,0.,50.);
    allPlots["chi2_iso"] = new TH1F("normchi2_iso",";#chi^2/n.d.o.f.;Events", 10,0.,10.);
    allPlots["chi2_veto"] = new TH1F("normchi2_veto",";#chi^2/n.d.o.f.;Events", 10,0.,10.);
    allPlots["lp_dxy_iso"] = new TH1F("lp_dxy_iso",";d_{xy} [cm];Events / 0.01 #mum", 20, 0, 0.2);
    allPlots["lp_dxy_veto"] = new TH1F("lp_dxy_veto",";d_{xy} [cm];Events / 0.01 #mum", 20, 0, 0.2);
    allPlots["lp_dz_iso"] = new TH1F("lp_dz_iso",";d_{z} [cm];Events / 0.01 #mum", 50, 0, 0.5);
    allPlots["lp_dz_veto"] = new TH1F("lp_dz_veto",";d_{z} [cm];Events / 0.01 #mum", 50, 0, 0.5);


  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  //for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  //GET LUMI INFO
  std::vector<RunPeriod_t> runPeriods = getRunPeriods(era);
  float totalLumi(0.);
  for(auto &it : runPeriods)
    totalLumi += it.second;

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      //evch = {};
      t->GetEntry(iev);
      //if(ev.nvtx>20) continue;
      if(iev%5000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));
      //Normalize to XSec and lumi
      float norm(1.0);
      float xsec(1.0);
      if(!isData) {
        norm =  normH ? normH->GetBinContent(1) : 1.0;
        xsec = normH ? normH->GetBinContent(2) : 0.;
        //if(xsec) norm*=xsec;
	//update nominal event weight
	if(ev.ttbar_nw>0) norm*=ev.ttbar_w[0];

        //Random run period based on lumi
        /*
        TString period = assignRunPeriod(runPeriods);
        runBCDEF.CheckRunPeriod(period);
        runGH.CheckRunPeriod(period);
        */
        //Scale wgt to run period's portion of total lumi
        /*
        for(auto &it : runPeriods) {
          if(it.first != period) continue;
          norm *= totalLumi/it.second;
        }
        */
      }
      runBCDEF.SetNorm(norm);
      runGH.SetNorm(norm);
      treeBCDEF.SetNorm(norm);
      treeGH.SetNorm(norm);
      treeBCDEF.SetLumi(19712.86);
      treeGH.SetLumi(16146.178);
      treeBCDEF.SetXsec(xsec);
      treeGH.SetXsec(xsec);
      int *piSFB = new int[ev.npf+50]();
      int *sumChBidx = new int[ev.npf]();
      int *keep = new int[ev.nj]();
      int *piSFG = new int[ev.npf+50]();

      //Apply top pT weight to ttbar events
      //https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Run_2_strategy
      //Particle(float pt, float eta, float phi, float mass, int pdgId, float relIso, int pid);
      std::vector<Particle> tops;
      if(isTTbar) {
        float top_pt_wgt(1.0);
        vector<float> pt;
        for(int i = 0; i < ev.ngtop; i++) {
          //Aviod stops stored by ntupelizer
          if(abs(ev.gtop_id[i]) != 6) continue;
          float tpt = ev.gtop_pt[i];
          if(tpt > 800) tpt = 800;
          pt.push_back(tpt);
          if(debug) std::cout << "Top pT= " << tpt << std::endl;
          tops.push_back(Particle(ev.gtop_pt[i], ev.gtop_eta[i], ev.gtop_phi[i], ev.gtop_m[i], ev.gtop_id[i], 0, 0));
        }
        //Save onlt hardest two tops (ttbar)
        //FIXME try max heap (priority_queue) to reduce time sorting
        //Might not need after imopsing |PdgId|==6
        std::sort(pt.begin(), pt.end(),
                  [] (float a, float b) { return a < b; } );
        //std::reverse(pt.begin(), pt.end());
        std::sort(tops.begin(), tops.end(), 
                 [] (Particle a, Particle b) { return a.Pt() < b.Pt(); } );
/*
            sort(genTracks.begin(), genTracks.end(),
                 [] (pfTrack a, pfTrack b) { return a.Pt() > b.Pt(); } );
*/
        //Calculate SFs based on expontial
        //https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting#Eventweight
        top_pt_wgt *= TMath::Exp(0.0615 - 0.0005*pt[0]);
        top_pt_wgt *= TMath::Exp(0.0615 - 0.0005*pt[1]);
        top_pt_wgt = TMath::Sqrt(top_pt_wgt);
        allPlots["topptwgt"]->Fill(0.,1.0);
        allPlots["topptwgt"]->Fill(1.,top_pt_wgt);
        allPlots2D["toppteta"]->Fill(tops[0].Eta(),top_pt_wgt);
        allPlots2D["toppteta"]->Fill(tops[1].Eta(),top_pt_wgt);
        runBCDEF.SetTopPtWgt(top_pt_wgt);
        runGH.SetTopPtWgt(top_pt_wgt);
        treeBCDEF.SetTopPtWgt(top_pt_wgt);
        treeGH.SetTopPtWgt(top_pt_wgt);
      }

      allPlots["nevt_all"]->Fill(1,norm);
      allPlots["nevt_all_no_weight"]->Fill(1);
      allPlots["norm_all"]->Fill(norm,norm);
      allPlots["nvtx_all"]->Fill(ev.nvtx,norm);

      //Basic lepton kinematics
      std::vector<int> tightLeptons,vetoLeptons;
      Leptons Muons(Tight,debug);
      Leptons Electrons(TightNoIso,debug);
      Leptons VetoLeptons(Veto,Loose,debug); // Designate Veto, only veto on Loose

      /* muon
      Muons.setMinPt(26);
      Muons.setMaxEta(2.1);
      Muons.setMaxRelIso(0.15);
      */
      Muons.setMinPt(20);
      Muons.setMaxEta(2.4);
      Muons.setMaxRelIso(0.15);

      //Electrons.setMinPt(12);
      Electrons.setMinPt(30);
      Electrons.setMaxEta(2.1); //trigger is eta2p1
      Electrons.setMaxRelIso(0.15);

      VetoLeptons.setMinPt(15);
      VetoLeptons.setMaxEta(2.4);
      VetoLeptons.setMaxRelIso(0.24);
      VetoLeptons.setMaxType(Loose);

      for(int il=0; il<ev.nl; il++)
	{
          //cout << "in lepton selection" << endl;
          Particle p(ev.l_pt[il], ev.l_eta[il], ev.l_phi[il], ev.l_mass[il], ev.l_id[il]*ev.l_charge[il], ev.l_relIso[il], ev.l_pid[il]);
          if(p.isMuon()) {
            Muons.addParticle(p); //only accepts tight
            VetoLeptons.addParticle(p); //only accepts loose FIXME think of more elegant way
          }
          else if(p.isElectron()) {
            Electrons.addParticle(p);
            VetoLeptons.addParticle(p);
            if(p.getType()==TightNoIso)
              allPlots["pid_all"]->Fill(1);
            if(p.getType()==Tight)
              allPlots["pid_all"]->Fill(2);
            else
              allPlots["pid_all"]->Fill(0);
          }
	}

      //Single Muon has tighter constraints
      if(Muons.size() == 1) {
        Muons.changeMinPt(26);
        Muons.changeMaxEta(2.1);
      }
      if(Muons.size() == 1) {
	allPlots["lp_pt_iso_m"]->Fill(Muons.getElement(0).Pt(),norm);
        if(VetoLeptons.size()==0)
	  allPlots["lp_pt_veto_m"]->Fill(Muons.getElement(0).Pt(),norm);
      }
      if(Electrons.size() == 1) {
        Electrons.changeMinPt(35);
        Electrons.changeParticleType(Tight); //TightNoIso -> Tight
      }
      if(Electrons.size() == 1) {
	allPlots["lp_pt_iso_e"]->Fill(Electrons.getElement(0).Pt(),norm);
        if(VetoLeptons.size()==0)
	  allPlots["lp_pt_veto_e"]->Fill(Electrons.getElement(0).Pt(),norm);
      }
      if(debug) cout << "lepton selection DONE" << endl;
      Leptons leptons(Tight,debug);
      leptons.combineLeptons(Muons);
      leptons.combineLeptons(Electrons);
      if(debug) cout << "sorting leptons" << endl;
      //FIXME try max heap (priority_queue) to reduce time sorting
      leptons.sortLeptonsByPt();

      allPlots["nevt_iso"]->Fill(1,norm);
      allPlots["norm_iso"]->Fill(norm,norm);
      allPlots["nvtx_iso"]->Fill(ev.nvtx,norm);
      
      //USE VETO HERE
      if(VetoLeptons.size()>0) continue; //veto only on lep+jets
      allPlots["nevt_veto"]->Fill(1,norm);
      allPlots["norm_veto"]->Fill(norm,norm);
      allPlots["nvtx_veto"]->Fill(ev.nvtx,norm);

      //check if triggers have fired
      //Trigger(muonTriggers, electronTriggers, debug=0)
      //Parse triggers
      Trigger trigger = Trigger(ev.muTrigger, ev.elTrigger, debug);
      //Check filetype (M/E/MM/EE/EM)
      trigger.setDataType(filename);

      //Dielectron
      //trigger.addRequiredDoubleElectronTrigger({"HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"});
      trigger.addRequiredDoubleElectronTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      //non-DZ
      //trigger.addRequiredDoubleElectronTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"); //FIXME not in Pedro's code, Carmen says prescaled
      //Prescaled

      //Dimuon
      trigger.addRequiredDoubleMuonTrigger({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"});
      //non-DZ
      //trigger.addRequiredDoubleMuonTrigger({"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"}); //Not in Pedro's code
      //Prescaled


      //Electron Muon (ME as well)
      trigger.addRequiredEMTrigger({"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"});
      //non-DZ
      trigger.addRequiredEMTrigger({"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"});

      //Single muon
      //Check only these triggers
      trigger.addRequiredMuonTrigger("HLT_IsoMu24_v");
      trigger.addRequiredMuonTrigger("HLT_IsoTkMu24_v");
      //trigger.addRequiredMuonTrigger{""HLT_IsoMu24_v","HLT_IsoTkMu24_v"}) also works


      //Single electron
      //trigger.addRequiredElectronTrigger("HLT_Ele27_WPTight_Gsf_v");
      trigger.addRequiredElectronTrigger("HLT_Ele32_eta2p1_WPTight_Gsf_v");
      //trigger.addRequiredElectronTrigger("HLT_Ele25_eta2p1_WPTight_Gsf_v");

      //Only triggers in Pedro's code
      //trigger.addRequiredElectronTrigger({"HLT_Ele32_eta2p1_WPTight_Gsf_v", "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v"});

      //decide the channel
      if(debug) cout << "decide channel" << endl;
      TString chTag("");

      if(debug) if(leptons.size()==0) cout << "NO LEPTONS!!" << endl;
      if(leptons.size()==0) continue;
      if(trigger.isSingleElectronEvent(leptons)) chTag="e";
      else if(trigger.isSingleMuonEvent(leptons)) chTag="m";
      else if(trigger.isDoubleElectronEvent(leptons)) chTag="ee";
      else if(trigger.isDoubleMuonEvent(leptons)) chTag="mm";
      else if(trigger.isEMEvent(leptons)) chTag="em";
      else continue;
      if(debug) cout << "check if Single Electron fired Single Muon trigger" << endl;
      if(leptons.size() == 2 && trigger.muonFired() && trigger.isElectronFile()) continue;
      chTag = "_"+chTag;
      if(debug) cout << "decide channel DONE" << endl;
      if(debug) cout << "Event: " << iev << endl;
      allPlots["nevt_trigger"]->Fill(1,norm);
      for(size_t i=0; i<leptons.size(); i++) {
        if(leptons[i].getType()==TightNoIso)
          allPlots["pid"+chTag]->Fill(1);
        else if(leptons[i].getType()==Tight)
          allPlots["pid"+chTag]->Fill(2);
        else
          allPlots["pid"+chTag]->Fill(0);
      }

      //one good lepton either isolated or in the non-isolated sideband or a Z candidate
      Bool_t isZ(false);//,isZPassingSIP3d(false);
      TLorentzVector l1p4,l2p4,dilp4;
      if(leptons.size()==2)
	{	  
          if(debug) cout << "di-lepton" << endl;
	  l1p4 = leptons[0].getVec();
	  l2p4 = leptons[1].getVec();
	  dilp4=l1p4+l2p4;
          if(leptons[0].getPdgId() == -leptons[1].getPdgId() &&
          //if(chTag=="mm" || chTag=="ee" &&
	     fabs(dilp4.M()-91)<15)
	    { 
	      isZ=true; 
	    }
          if(debug) cout << "di-lepton DONE" << endl;
	}

      //save lepton kinematics
      Float_t stsum(ev.met_pt[0]);
      for(size_t il=0; il<leptons.size(); il++)
	{
          stsum += leptons[il].Pt();
	}

      if(debug) cout << "Pion scale factors" << endl;
      //******************************
      //Pion tracker SFs
      if(!isData) {
        std::map<TString, std::map<TString, std::vector<double> > > trackEffMap =  getTrackingEfficiencyMap(era);
          applyTrackingEfficiencySF(ev, pf_eff, piSFB, passBit(runSysts,PI_BIT));
        /*
        if(runSysts==0) {
          applyTrackingEfficiencySF(ev, pf_eff, piSFB);
          //applyEtaDepTrackingEfficiencySF(ev, trackEffMap["BCDEF"]["nominal"], trackEffMap["BCDEF"]["binning"], piSFB);
          //applyEtaDepTrackingEfficiencySF(ev, trackEffMap["GH"]["nominal"], trackEffMap["GH"]["binning"], piSFG);
        }
        else if(passBit(runSysts,PI_BIT)<0) {
          applyEtaDepTrackingEfficiencySF(ev, trackEffMap["BCDEF"]["down"], trackEffMap["BCDEF"]["binning"], piSFB);
          applyEtaDepTrackingEfficiencySF(ev, trackEffMap["GH"]["down"], trackEffMap["GH"]["binning"], piSFG);
        }
        else if(passBit(runSysts,PI_BIT)>0) {
          applyEtaDepTrackingEfficiencySF(ev, trackEffMap["BCDEF"]["up"], trackEffMap["BCDEF"]["binning"], piSFB);
          applyEtaDepTrackingEfficiencySF(ev, trackEffMap["GH"]["up"], trackEffMap["GH"]["binning"], piSFG);
        }
        */
      }
      //******************************
      if(debug) cout << "Pion scale factors DONE!" << endl;

      //select jets
      float *pt_chargedB = new float[ev.nj]();
      float *pt_chargedG = new float[ev.nj]();
      float *j_pt_corrB = new float[ev.nj]();
      float pisfB(1.), pisfG(1.);
      for(int ij=0; ij<ev.nj; ij++)
        j_pt_corrB[ij]=1.;
      for(int ipf=0; ipf<ev.npf; ipf++) {
        if(ev.pf_fromPV[ipf]<2) continue;
        if(ev.pf_j[ipf]==-1) continue;
        //if(ev.pf_pt[ipf]<15) continue;
        if(ev.pf_c[ipf]==0) continue;
	TLorentzVector tkP4(0,0,0,0);
	tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);
        pfTrack pftk(tkP4, ev.pf_dxy[ipf], ev.pf_dxyE[ipf], ev.pf_dz[ipf], ev.pf_dzE[ipf], ev.pf_id[ipf],ev.pf_quality[ipf],ev.pf_highPurity[ipf]);
        float ptsf(1.),etasf(1.);
        //customSF(pftk, "BCDEF", ptsf, etasf, pf_eff);
        j_pt_corrB[ev.pf_j[ipf]] *= pftk.getPtCorrection();
        //pisfB *= customSF(pftk, "BCDEF", ptsf,etasf);
        //pisfB /= ptsf;
        if(ptsf > 0) {
          allPlots["piwgt_BCDEF"]->Fill(0.,1.0);
          allPlots["piwgt_BCDEF"]->Fill(1.,ptsf);
        }
        if(!isData && piSFB[ipf]>=0) { //FIXME
          sumChBidx[ipf]++;
          keep[ev.pf_j[ipf]]++;
          pt_chargedB[ev.pf_j[ipf]] += ev.pf_pt[ipf];
          //if(piSFB[ipf]>1 || (piSFB[ipf]==0 && piSFG[ipf]<1)) pt_chargedB[ev.pf_j[ipf]] += ev.pf_pt[ipf];
        }
        if(isData) pt_chargedB[ev.pf_j[ipf]] += ev.pf_pt[ipf];
        /*
        if(piSFG[ipf]>=0) {
          if(piSFG[ipf]>1 || (piSFG[ipf]==0 && piSFB[ipf]<1)) pt_chargedG[ev.pf_j[ipf]] += ev.pf_pt[ipf];
          pisfG *= customSF(pftk, "GH");
        }
        pt_chargedB[ev.pf_j[ipf]] += ev.pf_pt[ipf];
        pt_chargedG[ev.pf_j[ipf]] += ev.pf_pt[ipf];
        */
        pt_chargedG[ev.pf_j[ipf]] += ev.pf_pt[ipf];
        //Stor all PF pT and eta for data driven corrections
      }
      Float_t htsum(0),hsum(0),htbsum(0);
      int nbj(0);
      TLorentzVector jetDiff(0,0,0,0);
      std::vector<Jet> kJetsVec, lightJetsVec, allJetsVec, genJetsVec;
      KalmanEvent kalman(debug);
      kalman.loadEvent(ev);
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
	  //jp4=updateJES(jp4,ev.j_rawsf[k],ev.j_area[k],ev.rho,ev.nvtx,jetCorr);

	  //cross clean with respect to leptons deltaR<0.4
          bool overlapsWithLepton(false);
          for(size_t il=0; il<leptons.size(); il++) {
            if(jp4.DeltaR(leptons[il].getVec())>0.4) continue;  //Jet is fine
	    overlapsWithLepton=true;                   //Jet ovelaps with an "isolated" lepton, event is bad
          }
          if(overlapsWithLepton) continue;
          if(debug) cout << "Overlap with lepton DONE" << endl;

	  //smear jet energy resolution for MC
	  //jetDiff -= jp4;
	  float genJet_pt(0);
	  if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
	  if(!isData) { // && genJet_pt>0)
            float jerSmear(1.);
            double jer_sf = jres_sf->getScaleFactor({{JME::Binning::JetEta, jp4.Eta()}}, Variation::NOMINAL);
	    if(genJet_pt>0) {
              /* FIXME
              int smearIdx(0);
              if(passBit(runSysts,JER_BIT)<0) smearIdx=1;
              else if(passBit(runSysts,JER_BIT)>0) smearIdx=2;
	      jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt)[0];
              if(debug) std::cout << "JER " << jerSmear << std::endl;
	      if(passBit(runSysts,JER_BIT)<0) jerSmear=abs(passBit(runSysts,JER_BIT))*getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt)[1]; //systematics down
	      if(passBit(runSysts,JER_BIT)>0) jerSmear=abs(passBit(runSysts,JER_BIT))*getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt)[2]; //systematics up
              if(passBit(runSysts,JER_BIT)!=0 && debug) std::cout << "(+syst) JER " << jerSmear << std::endl;
              */
              double dPt = jp4.Pt() - genJet_pt;
              jerSmear = 1 + (jer_sf - 1.) * dPt / jp4.Pt();
	    }
            else { //stochastic smearing
              /*
              int smearIdx(0);
              if(passBit(runSysts,JER_BIT)<0) smearIdx=1;
              else if(passBit(runSysts,JER_BIT)>0) smearIdx=2;
              */
              double jet_resolution = jres->getResolution({{JME::Binning::JetPt, ev.j_pt[k]}, {JME::Binning::JetEta, ev.j_eta[k]}, {JME::Binning::Rho, ev.rho}});
              /*
              float jcor=getJetResolutionScales(jp4.Pt(),jp4.Eta())[0];
	      jcor+=passBit(runSysts,JER_BIT)*getJetResolutionScales(jp4.Pt(),jp4.Eta())[1];
              */
              //double jer_sf = jres_sf->getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::NOMINAL);
              TRandom *random = new TRandom3(0);
              jerSmear = 1 + random->Gaus(0, jet_resolution) * sqrt( std::max(jer_sf*jer_sf - 1., 0.) );
              delete random;
            }
            jp4 *= jerSmear;
            allPlots["jerwgt"]->Fill(0.,1.0);
            allPlots["jerwgt"]->Fill(1.,jerSmear);
            allPlots2D["jereta"]->Fill(jp4.Eta(),jer_sf);
            allPlots3D["jeceta"]->Fill(jp4.Pt(),jp4.Eta(),ev.j_rawsf[k]);
          }
	  //jetDiff += jp4;

	  // re-inforce kinematics cuts
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  htsum += jp4.Pt();
          hsum += jp4.P();
          if(ev.j_csv[k]>0.8484) { htbsum += jp4.Pt(); nbj++; }


	  //save jet
          //FIXME
          //Jet tmpj(jp4, ev.j_csv[k], k, pt_charged, ev.j_pz_charged[k], ev.j_p_charged[k], ev.j_pt_pf[k], ev.j_pz_pf[k], ev.j_p_pf[k], ev.j_g[k]); //Store pt of charged and total PF tracks and gen matched index
          Jet tmpj(jp4, ev.j_csv[k], k, ev.j_pt_charged[k], ev.j_pz_charged[k], ev.j_p_charged[k], ev.j_pt_pf[k], ev.j_pz_pf[k], ev.j_p_pf[k], ev.j_g[k]); //Store pt of charged and total PF tracks and gen matched index
          tmpj.updateChargedPt(pt_chargedB[k],pt_chargedG[k]);
          //delete[] pt_charged;
	  for(int ipf = 0; ipf < ev.npf; ipf++) {
	    if(ev.pf_j[ipf] != k) continue; //skip if PF track doesn't belong to current jet
	    if(ev.pf_c[ipf]==0) continue;   //skip if PF track is neutral
            if(ev.pf_eta[ipf]>2.4) continue; // |eta|<2.4
	    TLorentzVector tkP4(0,0,0,0);
	    tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);
            pfTrack pftk(tkP4, ev.pf_dxy[ipf], ev.pf_dxyE[ipf], ev.pf_dz[ipf], ev.pf_dzE[ipf], ev.pf_id[ipf],ev.pf_quality[ipf],ev.pf_highPurity[ipf]);
            if(abs(pftk.getPdgId())==13) {
              pftk.setGlobalMuon(ev.pf_globalMuon[ipf]);
              pftk.setTrackerMuon(ev.pf_trackerMuon[ipf]);
            }
            if(abs(pftk.getPdgId())==211 && ev.pf_mother[ipf]!=0) { //Handle promoted GEN particles from pion SFs
              pftk.setPromoted();
              pftk.setMother(ev.pf_mother[ipf]);
            }
	    tmpj.addTrack(pftk); //,ev.pf_id[ipf]);
	    if(passBit(runSysts,JSF_BIT)) {
              if(jecUnc) {
                jecUnc->setJetEta(jp4.Eta());
                jecUnc->setJetPt(jp4.Pt());
                                                   //true=UP; false=DOWN
                float jes = jecUnc->getUncertainty(passBit(runSysts,JSF_BIT)>0);
                jp4.Print();
	        jp4+=jp4*passBit(runSysts,JSF_BIT)*jes;
                jp4.Print();
                allPlots["jsfwgt"]->Fill(1.,jes);
                allPlots["jsfwgt"]->Fill(0.,1.);
                if(debug) std::cout << "jes= " << jes << std::endl;
              }
            }
	  }
          tmpj.sortTracksByPt();

          //if(debug && kalman.isGoodJet(k)) cout << "k=" << k << " jet pT=" << jp4.Pt() << endl;
          //if(kalman.isGoodJet(k)) kJetsVec.push_back(tmpj);
          if(!kalman.isGoodJet(k)) lightJetsVec.push_back(tmpj);
          allJetsVec.push_back(tmpj);
	}
      //after main jet b/c smearing MiniEvet_t is in there
      kalman.loadEvent(ev); //load again to get smeared versions
      for(auto &jet : kalman.getJets()) {
        jet.updateChargedPt(pt_chargedB[jet.getJetIndex()],pt_chargedG[jet.getJetIndex()]);
        if(jet.getTracks().size()<1) continue; //skip jets with no sub-structure
	TLorentzVector jp4 = jet.getVec();
	//cross clean with respect to leptons deltaR<0.4
        bool overlapsWithLepton(false);
        for(size_t il=0; il<leptons.size(); il++) {
          if(jp4.DeltaR(leptons[il].getVec())>0.4) continue;  //Jet is fine
	  overlapsWithLepton=true;                   //Jet ovelaps with an "isolated" lepton, event is bad
        }
        if(overlapsWithLepton) continue;
        if(debug) cout << "Overlap with lepton DONE" << endl;
        jet.sortTracksByPt();
        kJetsVec.push_back(jet);
        //allJetsVec.push_back(jet);
      }
      for (int k=0; k<ev.ng;k++) {
        //check kinematics
        if(abs(ev.g_id[k])==13 || abs(ev.g_id[k])==11) continue; //skip leptons since they are stored after genJets
	TLorentzVector jp4;
	jp4.SetPtEtaPhiM(ev.g_pt[k],ev.g_eta[k],ev.g_phi[k],ev.g_m[k]);

	// re-inforce kinematics cuts
	if(jp4.Pt()<30) continue;
	if(fabs(jp4.Eta()) > 2.4) continue;
	
	//save jet
        Jet tmpgj(jp4, ev.g_id[k], k, ev.xb[k]);
        genJetsVec.push_back(tmpgj);
      }
      if(debug) cout << kJetsVec.size() << " Kalman jets found" << endl;
      if(debug) cout << allJetsVec.size() << " jets found" << endl;
      //delete[] pt_chargedB;
      //delete[] pt_chargedG;

      if(htsum<80) continue;
      //if(htsum<160) continue;
      //if(allJetsVec[0].getVec().Pt()<50) continue;
      stsum += htsum;

      
      //event weight
      float wgt(1.0);

      std::vector<float> puWgts(3,1.0),topPtWgts(2,1.0);
      EffCorrection_t lepSelCorrWgt(1.0,0.0), triggerCorrWgt(1.0,0.0);
      if(debug) cout << "Lepton scale factors" << endl;
      if(!isData)
	{
	  //account for pu weights and effect on normalization
	  allPlots["puwgtctr"]->Fill(0.,1.0);
	  allPlots["puwgtgr"]->Fill(0.,1.0);
	  allPlots["puwgt"]->Fill(0.,1.0);
	  if(debug) cout << "getting puWgts" << endl;
          /*
          for(int xbin=1; xbin<=puWgt[0]->GetXaxis()->GetNbins(); xbin++) {
	    Double_t yobs;
	    yobs = puWgt[0]->GetBinContent(xbin);
            allPlots["puwgt"]->Fill(xbin,yobs);
          }
          */
          //Set PU weight for each run period [BCDEFup, BCDEF, BCDEFdown, GHup, GH, GHdown]
          runBCDEF.SetPuWgt(puWgtsRun[1+passBit(runSysts,PU_BIT)]->GetBinContent(ev.putrue));
          runGH.SetPuWgt(puWgtsRun[4+passBit(runSysts,PU_BIT)]->GetBinContent(ev.putrue));
          treeBCDEF.SetPuWgt(puWgtsRun[1+passBit(runSysts,PU_BIT)]->GetBinContent(ev.putrue));
          treeGH.SetPuWgt(puWgtsRun[4+passBit(runSysts,PU_BIT)]->GetBinContent(ev.putrue));

	  if(debug) cout << "getting puWgts DONE!" << endl;
	  //trigger/id+iso efficiency corrections
          if(debug) cout << "calling trigger function" << endl;
          std::vector<int> pdgIds; //vector of IDs for trigger correction function
	  //triggerCorrWgt=lepEffH.getTriggerCorrection(leptons);
          if(debug) cout << "calling trigger function DONE!" << endl;

          // ** SFs for BCDEF and GH separately
          EffCorrection_t lepSelCorrWgt_BCDEF(1.0,0.0), triggerCorrWgt_BCDEF(1.0,0.0);
          EffCorrection_t lepSelCorrWgt_GH(1.0,0.0), triggerCorrWgt_GH(1.0,0.0);
	  triggerCorrWgt_BCDEF=lepEffH_BCDEF.getTriggerCorrection(leptons);
	  triggerCorrWgt_GH=lepEffH_GH.getTriggerCorrection(leptons);
          // Systematics for Trigger
	  triggerCorrWgt_BCDEF.first+=passBit(runSysts,TRIGGER_BIT)*triggerCorrWgt_BCDEF.second;
	  triggerCorrWgt_GH.first+=passBit(runSysts,TRIGGER_BIT)*triggerCorrWgt_GH.second;
	  for(size_t il=0; il<leptons.size(); il++) {
	    EffCorrection_t selSF_BCDEF=lepEffH_BCDEF.getOfflineCorrection(leptons[il], ev.nvtx);
	    EffCorrection_t selSF_GH=lepEffH_GH.getOfflineCorrection(leptons[il], ev.nvtx);
	    lepSelCorrWgt_BCDEF.second = sqrt( pow(lepSelCorrWgt_BCDEF.first*selSF_BCDEF.second,2)+pow(lepSelCorrWgt_BCDEF.second*selSF_BCDEF.first,2));
	    lepSelCorrWgt_GH.second = sqrt( pow(lepSelCorrWgt_GH.first*selSF_GH.second,2)+pow(lepSelCorrWgt_GH.second*selSF_GH.first,2));
            if(debug) cout << "lepSelCorrWgt_BCDEF=" << lepSelCorrWgt_BCDEF.first << endl;
            if(debug) cout << "selSF=" << selSF_BCDEF.first << endl;
            if(debug) cout << "lepSelCorrWgt_GH=" << lepSelCorrWgt_GH.first << endl;
            if(debug) cout << "selSF=" << selSF_GH.first << endl;
                                                                // Systematics for lepton selection
	    lepSelCorrWgt_BCDEF.first *= selSF_BCDEF.first + passBit(runSysts,LEP_BIT)*selSF_BCDEF.second;
	    lepSelCorrWgt_GH.first *= selSF_GH.first + passBit(runSysts,LEP_BIT)*selSF_GH.second;
            if(passBit(runSysts,LEP_BIT)) {
              if(debug) cout << "(+syst) lepSelCorrWgt_BCDEF=" << lepSelCorrWgt_BCDEF.first << endl;
              if(debug) cout << "(+syst) lepSelCorrWgt_GH=" << lepSelCorrWgt_GH.first << endl;
             }
	  }
          runBCDEF.SetSFs(triggerCorrWgt_BCDEF.first*lepSelCorrWgt_BCDEF.first,lepSelCorrWgt_BCDEF.second);
          runGH.SetSFs(triggerCorrWgt_GH.first*lepSelCorrWgt_GH.first,lepSelCorrWgt_GH.second);
          treeBCDEF.SetSFs(triggerCorrWgt_BCDEF.first*lepSelCorrWgt_BCDEF.first);//,lepSelCorrWgt_BCDEF.second);
          treeGH.SetSFs(triggerCorrWgt_GH.first*lepSelCorrWgt_GH.first);//,lepSelCorrWgt_GH.second);
          // **
	  wgt=triggerCorrWgt_BCDEF.first*lepSelCorrWgt_BCDEF.first*(puWgtsRun[1]->GetBinContent(ev.putrue))*norm;

          if(debug) cout << "weight=" << wgt << endl;
          if(debug) cout << "Trigger=" << triggerCorrWgt.first << endl << "Lepton=" << lepSelCorrWgt.first << endl << "PU=" << puWgts[0] << endl << "norm=" << norm  << endl;
	}
      else
        if(debug) cout << "weight=" << wgt << " norm=" << norm << endl;
      if(debug) cout << "Lepton scale factors DONE!" << endl;


      //sort by Pt
      if(debug) cout << "sorting jets" << endl;
      //FIXME try max heap (priority_queue) to reduce time sorting
      sort(lightJetsVec.begin(),    lightJetsVec.end(),   sortJetsByPt);
      sort(kJetsVec.begin(),    kJetsVec.end(),   sortJetsByPt);
      sort(allJetsVec.begin(),  allJetsVec.end(), sortJetsByPt);

      for(size_t il=0; il<leptons.size(); il++) {
        for(size_t ij=0; ij<allJetsVec.size(); ij++) {
	  TLorentzVector jp4=allJetsVec[ij].getVec();
          allPlots["dR"+chTag]->Fill(jp4.DeltaR(leptons[il].getVec()),wgt);
          allPlots["dR"+chTag+"_no_weight"]->Fill(jp4.DeltaR(leptons[il].getVec()),norm);
        }
      }


      //MET and transverse mass
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
      met+=jetDiff;
      met.SetPz(0.); met.SetE(met.Pt());
      float mt( computeMT(isZ ? dilp4: leptons[0].getVec(),met) );

      //simple fill
      bool singleLep(false);
      bool doubleLep(false);
      bool minJets(false);

      if(debug) cout << "starting simple plots" << endl;
      allPlots["nevt"+chTag]->Fill(1,norm);
      allPlots["weight"+chTag]->Fill(wgt,norm);
      allPlots["norm"+chTag]->Fill(norm,norm);
      allPlots["nvtx"+chTag]->Fill(ev.nvtx,wgt);
      runBCDEF.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag);
      runGH.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag);

      runBCDEF.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag);
      runGH.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag);

      runBCDEF.Fill(leptons, chTag);
      runGH.Fill(leptons, chTag);

      allPlots["nj"+chTag]->Fill(allJetsVec.size(),wgt);
      allPlots["nlj"+chTag]->Fill(lightJetsVec.size(),wgt);
      allPlots["nkj"+chTag]->Fill(kJetsVec.size(),wgt);
      allPlots["nlj"+chTag+"_no_weight"]->Fill(lightJetsVec.size(),norm);
      allPlots["nkj"+chTag+"_no_weight"]->Fill(kJetsVec.size(),norm);
      allPlots["nlj_all"]->Fill(lightJetsVec.size(),wgt);
      allPlots["nkj_all"]->Fill(kJetsVec.size(),wgt);
      allPlots["nlp"+chTag]->Fill(leptons.size(),wgt);
      allPlots["nlp"+chTag+"_no_weight"]->Fill(leptons.size(),norm);
      allPlots["HT"+chTag]->Fill(htsum,wgt);
      allPlots["ST"+chTag]->Fill(stsum,wgt);
      allPlots["MET2oST"+chTag]->Fill(pow(ev.met_pt[0],2)/stsum,wgt);

      if(leptons.size() > 0) {
        allPlots["lp_pt_low"+chTag]->Fill(leptons[0].Pt(),wgt);
        allPlots["lp_pt"+chTag]->Fill(leptons[0].Pt(),wgt);
        allPlots["lp_pt"+chTag+"_no_weight"]->Fill(leptons[0].Pt(),norm);
        allPlots["lp_pt_all"]->Fill(leptons[0].Pt(),wgt);
        allPlots["relIso"+chTag]->Fill(leptons[0].getRelIso(),wgt);
        allPlots["lp_eta"+chTag]->Fill(leptons[0].Eta(),wgt);
        allPlots["lp_eta"+chTag+"_no_weight"]->Fill(leptons[0].Eta(),norm);
        allPlots["lp_phi"+chTag]->Fill(leptons[0].Phi(),wgt);
        allPlots["lp_phi"+chTag+"_no_weight"]->Fill(leptons[0].Phi(),norm);
      }

      if(isZ) {
        allPlots["massZ"+chTag]->Fill(dilp4.M(),wgt);
        allPlots["massZ"+chTag+"_no_weight"]->Fill(dilp4.M(),norm);
      }

      if(allJetsVec.size() > 0) {
        allPlots["j_pt"+chTag]->Fill(allJetsVec[0].getVec().Pt(),wgt);
        allPlots["j_pt_low"+chTag]->Fill(allJetsVec[0].getVec().Pt(),wgt);
        allPlots["j_pt"+chTag+"_no_weight"]->Fill(allJetsVec[0].getVec().Pt(),norm);
        allPlots["j_pt_all"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
      }
      if(lightJetsVec.size() > 0) {
        allPlots["lj_pt_low"+chTag]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
        allPlots["lj_pt"+chTag]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
        allPlots["lj_pt"+chTag+"_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
        allPlots["lj_pt_all"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
      }
      if(kJetsVec.size() > 0) {
        allPlots["kj_pt_low"+chTag]->Fill(kJetsVec[0].getVec().Pt(),wgt);
        allPlots["kj_pt"+chTag]->Fill(kJetsVec[0].getVec().Pt(),wgt);
        allPlots["kj_pt"+chTag+"_no_weight"]->Fill(kJetsVec[0].getVec().Pt(),norm);
        allPlots["kj_pt_all"]->Fill(kJetsVec[0].getVec().Pt(),wgt);
      }

      if(debug) cout << "starting lep/jets plots" << endl;
      //Require exactly 1 lepton
      if(leptons.size() == 1) {
        if(debug) cout << "single lepton" << endl;
        singleLep = true;

        allPlots["lp_pt"+chTag+"_lep"]->Fill(leptons[0].Pt(),wgt);
        allPlots["lp_pt"+chTag+"_lep_no_weight"]->Fill(leptons[0].Pt(),norm);

        allPlots["lp_pt_all_lep"]->Fill(leptons[0].Pt(),wgt);
        allPlots["lp_eta"+chTag+"_lep"]->Fill(leptons[0].Eta(),wgt);
        allPlots["lp_eta"+chTag+"_lep"+"_no_weight"]->Fill(leptons[0].Eta(),norm);
        allPlots["lp_phi"+chTag+"_lep"]->Fill(leptons[0].Phi(),wgt);
        allPlots["lp_phi"+chTag+"_lep"+"_no_weight"]->Fill(leptons[0].Phi(),norm);
        allPlots["MET"+chTag+"_lep"]->Fill(ev.met_pt[0],wgt);
        allPlots["HT"+chTag+"_lep"]->Fill(htsum,wgt);
        allPlots["ST"+chTag+"_lep"]->Fill(stsum,wgt);
        allPlots["MET2oST"+chTag+"_lep"]->Fill(pow(ev.met_pt[0],2)/stsum,wgt);
        allPlots["MET"+chTag+"_lep_no_weight"]->Fill(ev.met_pt[0],norm);
        allPlots["relIso"+chTag+"_lep"]->Fill(leptons[0].getRelIso(),wgt);
        if(debug) cout << "is" << (kalman.isGoodEvent() ? "" : " not") << " a Kalman event" << endl;
        //Require at least 1 Kalman jet and 1 at least light jet
        if(kalman.isGoodEvent() && kJetsVec.size()>0 && lightJetsVec.size() >= 1) {
          if(debug) cout << "jet requirements" << endl;
          minJets = true;

          allPlots["lp_pt"+chTag+"_lepjets"]->Fill(leptons[0].Pt(),wgt);
          //allPlots["lp_pt"+chTag+"_lepjets_no_weight"]->Fill(leptons[0].Pt(),norm);

          allPlots["lp_pt_all_lepjets"]->Fill(leptons[0].Pt(),wgt);
          allPlots["lp_eta"+chTag+"_lepjets"]->Fill(leptons[0].Eta(),wgt);
          allPlots["lp_eta"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons[0].Eta(),norm);
          allPlots["lp_phi"+chTag+"_lepjets"]->Fill(leptons[0].Phi(),wgt);
          allPlots["lp_phi"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons[0].Phi(),norm);
          allPlots["MET"+chTag+"_lepjets"]->Fill(ev.met_pt[0],wgt);
          allPlots["HT"+chTag+"_lepjets"]->Fill(htsum,wgt);
          allPlots["ST"+chTag+"_lepjets"]->Fill(stsum,wgt);
          allPlots["MET2oST"+chTag+"_lepjets"]->Fill(pow(ev.met_pt[0],2)/stsum,wgt);
          allPlots["relIso"+chTag+"_lepjets"]->Fill(leptons[0].getRelIso(),wgt);
          //allPlots["MET"+chTag+"_lepjets__no_weight"]->Fill(ev.met_pt[0],norm);
        }
      }
      //Require exactly 2 leptons
      else if(leptons.size() == 2) {
        if(debug) cout << "dilepton" << endl;
        doubleLep = true;
        //Z control plot
        if(isZ) {
          allPlots["massZ"+chTag+"_lep"]->Fill(dilp4.M(),wgt);
          allPlots["massZ"+chTag+"_lep_no_weight"]->Fill(dilp4.M(),norm);
        }
        if(abs(dilp4.M()-91)<15)
          allPlots["chargeZ"+chTag]->Fill(leptons[0].charge()*leptons[1].charge(),wgt); 
        //Exclude Z and low mass and require same falvor dilepton MET > 40 GeV
        if(!isZ && (dilp4.M() > 20 && leptons[0].getPdgId()==-leptons[1].getPdgId()) &&
          ((abs(leptons[0].getPdgId())!=abs(leptons[1].getPdgId())) || (leptons[0].getPdgId()==-leptons[1].getPdgId() && met.Pt() > 40))) {
          allPlots["ndilp"+chTag+"_lep"]->Fill(leptons.size(),wgt);
          allPlots["dilp_pt"+chTag+"_lep"]->Fill(dilp4.Pt(),wgt);
          allPlots["dilp_m"+chTag+"_lep"]->Fill(dilp4.M(),wgt);
          allPlots["ndilp"+chTag+"_lep"+"_no_weight"]->Fill(leptons.size(),norm);
          allPlots["dilp_pt"+chTag+"_lep"+"_no_weight"]->Fill(dilp4.Pt(),norm);
          allPlots["dilp_m"+chTag+"_lep"+"_no_weight"]->Fill(dilp4.M(),norm);
          allPlots["lp_pt"+chTag+"_lep"]->Fill(leptons[0].Pt(),wgt);
          allPlots["l2p_pt"+chTag+"_lep"]->Fill(leptons[1].Pt(),wgt);
          allPlots["lp_pt"+chTag+"_lep"+"_no_weight"]->Fill(leptons[0].Pt(),norm);
          allPlots["l2p_pt"+chTag+"_lep"+"_no_weight"]->Fill(leptons[1].Pt(),norm);
          allPlots["MET"+chTag+"_lep"]->Fill(ev.met_pt[0],wgt);
          allPlots["charge"+chTag+"_lep"]->Fill(leptons[0].charge()*leptons[1].charge(),wgt);
          allPlots["MET"+chTag+"_lep"+"_no_weight"]->Fill(ev.met_pt[0],norm);
          allPlots["charge"+chTag+"_lep"+"_no_weight"]->Fill(leptons[0].charge()*leptons[1].charge(),norm);
          allPlots["relIso"+chTag+"_lep"]->Fill(leptons[0].getRelIso(),wgt);
          allPlots["relIso"+chTag+"_lep"]->Fill(leptons[1].getRelIso(),wgt);
        }
        //Require at least 1 Kalman jet
        if(kalman.isGoodEvent() && kJetsVec.size()>0) {
          if(debug) cout << "jet requirements" << endl;
          //Z control plot
          minJets = true;
          if(isZ) {
            allPlots["massZ"+chTag+"_lepjets"]->Fill(dilp4.M(),wgt);
            allPlots["massZ"+chTag+"_lepjets_no_weight"]->Fill(dilp4.M(),norm);
          }
          //Exclude Z mass
          if(isZ) minJets = false;
          //Exclude low mass (M < 20 GeV)
          if(dilp4.M() < 20 && leptons[0].getPdgId()==-leptons[1].getPdgId()) minJets = false;
          //Require same falvor dilepton MET > 40 GeV
          if(leptons[0].getPdgId()==-leptons[1].getPdgId() && met.Pt() < 40) minJets = false;

          if(minJets) {
            allPlots["ndilp"+chTag+"_lepjets"]->Fill(leptons.size(),wgt);
            allPlots["dilp_pt"+chTag+"_lepjets"]->Fill(dilp4.Pt(),wgt);
            allPlots["dilp_m"+chTag+"_lepjets"]->Fill(dilp4.M(),wgt);
            allPlots["ndilp"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons.size(),norm);
            allPlots["dilp_pt"+chTag+"_lepjets"+"_no_weight"]->Fill(dilp4.Pt(),norm);
            allPlots["dilp_m"+chTag+"_lepjets"+"_no_weight"]->Fill(dilp4.M(),norm);
            allPlots["lp_pt"+chTag+"_lepjets"]->Fill(leptons[0].Pt(),wgt);
            allPlots["l2p_pt"+chTag+"_lepjets"]->Fill(leptons[1].Pt(),wgt);
            allPlots["lp_pt"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons[0].Pt(),norm);
            allPlots["l2p_pt"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons[1].Pt(),norm);
            allPlots["lp_eta"+chTag+"_lepjets"]->Fill(leptons[0].Eta(),wgt);
            allPlots["l2p_eta"+chTag+"_lepjets"]->Fill(leptons[1].Eta(),wgt);
            allPlots["lp_eta"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons[0].Eta(),norm);
            allPlots["l2p_eta"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons[1].Eta(),norm);
            allPlots["lp_phi"+chTag+"_lepjets"]->Fill(leptons[0].Phi(),wgt);
            allPlots["l2p_phi"+chTag+"_lepjets"]->Fill(leptons[1].Phi(),wgt);
            allPlots["lp_phi"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons[0].Phi(),norm);
            allPlots["l2p_phi"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons[1].Phi(),norm);
            allPlots["MET"+chTag+"_lepjets"]->Fill(ev.met_pt[0],wgt);
            allPlots["charge"+chTag+"_lepjets"]->Fill(leptons[0].charge()*leptons[1].charge(),wgt);
            allPlots["MET"+chTag+"_lepjets"+"_no_weight"]->Fill(ev.met_pt[0],norm);
            allPlots["charge"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons[0].charge()*leptons[1].charge(),norm);
            allPlots["relIso"+chTag+"_lepjets"]->Fill(leptons[0].getRelIso(),wgt);
            allPlots["relIso"+chTag+"_lepjets"]->Fill(leptons[1].getRelIso(),wgt);
          }
        }
      }
      if(debug) cout << "simple plots DONE" << endl;


      //Require lep+jets or dilepton
      if(!singleLep && !doubleLep) continue;
      if(debug) cout << "passed lep requirements" << endl;

      runBCDEF.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag, "lep");
      runGH.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag, "lep");

      runBCDEF.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag, "lep");
      runGH.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag, "lep");
      
      runBCDEF.Fill(leptons, chTag, "lep");
      runGH.Fill(leptons, chTag, "lep");

      allPlots["npf"+chTag+"_lep"]->Fill(ev.npf,wgt);
      //allPlots["npf"+chTag+"_lep"+"_no_weight"]->Fill(ev.npf,norm);
      allPlots["nevt"+chTag+"_lep"]->Fill(1,norm);
      allPlots["weight"+chTag+"_lep"]->Fill(wgt,norm);
      allPlots["norm"+chTag+"_lep"]->Fill(norm,norm);
      allPlots["nvtx"+chTag+"_lep"]->Fill(ev.nvtx,wgt);
      allPlots["nevt_all_lep"]->Fill(1,norm);

      allPlots["nj"+chTag+"_lep"]->Fill(allJetsVec.size(),wgt);
      allPlots["nlj"+chTag+"_lep"]->Fill(lightJetsVec.size(),wgt);
      allPlots["nkj"+chTag+"_lep"]->Fill(kJetsVec.size(),wgt);
      //allPlots["nj"+chTag+"_lep"+"_no_weight"]->Fill(lightJetsVec.size(),norm);
      //allPlots["nkj"+chTag+"_lep"+"_no_weight"]->Fill(kJetsVec.size(),norm);
      allPlots["nlp"+chTag+"_lep"]->Fill(leptons.size(),wgt);

      if(lightJetsVec.size() > 0 and kJetsVec.size() > 0) {
        allPlots["j_pt"+chTag+"_lep"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
        allPlots["lj_pt"+chTag+"_lep"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
        allPlots["kj_pt"+chTag+"_lep"]->Fill(kJetsVec[0].getVec().Pt(),wgt);
      }

      /*
      TLorentzVector pmt;
      for(auto &it : allJetsVec)
        pmt += it.getVec();
      for(size_t il = 0; il < leptons.size(); il++)
        pmt += leptons[il].getVec();
      allPlots["mt"+chTag+"_lep"]->Fill(pmt.Mt(),wgt);
      allPlots["mt_all_lep"]->Fill(pmt.Mt(),wgt);
      

      //Require b-tagged and light jets
      if(!minJets) continue;

      pmt = TLorentzVector();
      for(auto &it : allJetsVec)
        pmt += it.getVec();
      for(size_t il = 0; il < leptons.size(); il++)
        pmt += leptons[il].getVec();
      allPlots["mt"+chTag+"_lepjets"]->Fill(pmt.Mt(),wgt);
      allPlots["mt_all_lepjets"]->Fill(pmt.Mt(),wgt);
      */
      allPlots["mt"+chTag+"_lep"]->Fill(mt,wgt);
      allPlots["mt_all_lep"]->Fill(mt,wgt);

      if(debug) cout << "passed jet requirements" << endl;

      //Fill gen-level top plots
      if(isTTbar) {
        runBCDEF.FillGen(tops, chTag, "lepjets");
        runGH.FillGen(tops, chTag, "lepjets");
      }

      runBCDEF.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag, "lepjets");
      runGH.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag, "lepjets");
      
      runBCDEF.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag, "lepjets");
      runGH.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag, "lepjets");
      
      runBCDEF.Fill(leptons, chTag, "lepjets");
      runGH.Fill(leptons, chTag, "lepjets");

      allPlots["npf"+chTag+"_lepjets"]->Fill(ev.npf,wgt);
      allPlots["npf"+chTag+"_lepjets"+"_no_weight"]->Fill(ev.npf,norm);
      allPlots["nevt"+chTag+"_lepjets"]->Fill(1,norm);
      allPlots["weight"+chTag+"_lepjets"]->Fill(wgt,norm);
      allPlots["norm"+chTag+"_lepjets"]->Fill(norm,norm);
      allPlots["nvtx"+chTag+"_lepjets"]->Fill(ev.nvtx,wgt);
      allPlots["nevt_all_lepjets"]->Fill(1,norm);

      allPlots["nj"+chTag+"_lepjets"]->Fill(allJetsVec.size(),wgt);
      allPlots["nlj"+chTag+"_lepjets"]->Fill(lightJetsVec.size(),wgt);
      allPlots["nkj"+chTag+"_lepjets"]->Fill(kJetsVec.size(),wgt);
      allPlots["nlp"+chTag+"_lepjets"]->Fill(leptons.size(),wgt);

      if(allJetsVec.size()>0)   allPlots["j_pt"+chTag+"_lepjets"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
      if(lightJetsVec.size()>0) allPlots["lj_pt"+chTag+"_lepjets"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
      if(kJetsVec.size()>0)     allPlots["kj_pt"+chTag+"_lepjets"]->Fill(kJetsVec[0].getPt(),wgt);


      //charmed resonance analysis : use only jets with CSV>CSVL, up to two per event
      if(htsum<100) continue;
      evch.njpsi=0;
      evch.nmeson=0;
      evch.nj=0;
      //Better J/Psi (Just look how much shorter it is!)
      const float gMassMu(0.1057),gMassK(0.4937),gMassPi(0.1396);
      if(kalman.isJPsiEvent()) {
        for(auto &jet : kJetsVec) {
          vector<pfTrack> muTracks;
          for(auto &track : jet.getTracks()) {
            if(track.getMotherId()!=443) continue;
            if(abs(track.getPdgId())==13) { track.setMass(gMassMu); muTracks.push_back(track); }
            //if(abs(track.getPdgId())==13) { cout << endl << ev.event << ": " << track.Pt() << " " << track.Eta() << " " << track.Phi() <<  " " << ev.k_mass[0] << endl; }
            //if(abs(track.getPdgId())==13) { cout << endl << ev.event << ": " << ev.k_pf_pt[0] << " " << ev.k_pf_eta[0] << " " << ev.k_pf_phi[0] <<  " " << ev.k_mass[0] << endl; }
            //if(abs(track.getPdgId())==13) { cout << endl << ev.event << ": " << ev.k_pf_pt[1] << " " << ev.k_pf_eta[1] << " " << ev.k_pf_phi[1] <<  " " << ev.k_mass[1] << endl; }
          }
          if(muTracks.size()<2) continue;
          Jet genJet;
          for(auto & gjet : genJetsVec) {
            if(gjet.getVec().DeltaR(jet.getVec())>0.1) continue; //find good dR match
            if(gjet.getVec().DeltaR(jet.getVec()) > genJet.getVec().DeltaR(jet.getVec())) continue; //find tighter dR match
            genJet = gjet;
            break;
          }
          //std::vector<float> frag = {ev.xb[genJet.getJetIndex()],ev.peterson[genJet.getJetIndex()],ev.up[genJet.getJetIndex()],ev.central[genJet.getJetIndex()],ev.down[genJet.getJetIndex()],};
          std::vector<float> frag = {ev.up[genJet.getJetIndex()],ev.central[genJet.getJetIndex()],ev.down[genJet.getJetIndex()]};

          std::vector<pfTrack> pfmuMatched, pfmuReject;
          //Gen-matching
          if(!isData) {
            std::vector<pfTrack> genTracks;
            std::vector<pfTrack> genMuTracks;
            for(int ig = 0; ig < ev.ngpf; ig++) {
              if(abs(ev.gpf_id[ig])!=13) continue;
              TLorentzVector gen;
              gen.SetPtEtaPhiM(ev.gpf_pt[ig], ev.gpf_eta[ig], ev.gpf_phi[ig], gMassMu);
              genTracks.push_back(pfTrack(gen,0,0,0,0,ev.gpf_id[ig],3,true));
              genTracks.back().setMother(ev.gpf_mother[ig]); 
              //cout << genTracks.back().getGenT() << endl; //daug -> mother id -> mother ttbar
              if(abs(ev.gpf_id[ig])==13) genMuTracks.push_back(pfTrack(gen,0,0,0,0, ev.gpf_id[ig],3,true));
            }
            //sort GEN traks
            sort(genMuTracks.begin(), genMuTracks.end(),
                 [] (pfTrack a, pfTrack b) { return a.Pt() > b.Pt(); } );

            for(auto & it : muTracks) { //FIXME reference might not work
              double dR = 0.1; //initial dR
              int best_idx = -1;
              for(auto & itg : genMuTracks) {
                if(it.getPdgId() != itg.getPdgId()) continue; //insure ID and charge
                if(it.getVec().DeltaR(itg.getVec())>dR) continue; //find dR
                //if(((it.Pt()-itg.Pt())/it.Pt())>0.10) continue; //gen and reco less than 10% difference
                dR = it.getVec().DeltaR(itg.getVec());
                best_idx = &itg - &genMuTracks[0]; //get index on current closest gen particle
              }
              if(best_idx<0) { //no gen track matched
                pfmuReject.push_back(it);
              }
              else {
                pfmuMatched.push_back(it);
                pfmuMatched.back().setMother(genMuTracks[best_idx].getMotherId());
                genMuTracks.erase(genMuTracks.begin() + best_idx); //remove gen track so it cannot be matched again!
              }
            }
          }

          /*
          float *pt_chargedB_jpsi = new float[ev.nj]();
          int *sumChBidx = new int[ev.npf]();
          applyTrackingEfficiencySF(ev, pf_eff, piSFB, passBit(runSysts,PI_BIT));
          if(!isData && piSFB[ipf]>=0) { //FIXME
            sumChBidx[ipf]++;
            keep[ev.pf_j[ipf]]++;
            pt_chargedB_jpsi[ev.pf_j[ipf]] += ev.pf_pt[ipf];
            //if(piSFB[ipf]>1 || (piSFB[ipf]==0 && piSFG[ipf]<1)) pt_chargedB[ev.pf_j[ipf]] += ev.pf_pt[ipf];
          }
          int *shouldDrop = new int[ev.npf]();
          applyTrackingEfficiencySF(ev, pf_eff_jpsi, shouldDrop, passBit(runSysts,PI_BIT));
          */
          for(size_t i = 0; i < muTracks.size(); i++) { //0,1 is only ~35% of all J/Psi
            //FIXME for(size_t j = 0; j < muTracks.size(); j++) {
            for(size_t j = i+1; j < muTracks.size(); j++) {
            //int i(0),j(1);
              if(i==j) continue;
              float mass12 = (muTracks[i].getVec()+muTracks[j].getVec()).M();
              if(mass12>2.8 && mass12<3.4) {
                //if(debug) cout << pfmuCands[0].Pt() << " " << pfmuCands[0].Eta() << " " << pfmuCands[0].Phi() << " " << gMassMu << endl;
                //if(debug) cout << pfmuCands[1].Pt() << " " << pfmuCands[1].Eta() << " " << pfmuCands[1].Phi() << " " << gMassMu << endl;
                if(debug) cout << mass12 << endl << endl;
                if(debug) cout << "J/Psi found" << endl;
                if(debug && (mass12>3.0 && mass12<3.2)) cout << "and it's good!" << endl;
                allPlots["massJPsi"+chTag]->Fill(mass12,wgt);
	        allPlots["massJPsi_all"]->Fill(mass12,wgt);

                TLorentzVector jpsi = muTracks[i].getVec() + muTracks[j].getVec();
                if(rbFit && !isData) {
                  float rbWgtJPsi(1.);
                  int idx = (runPeriod.Contains("BCDEF") ? 0 : 1);
                  if(rbFit==1) rbWgtJPsi=(rbWgt->Eval(jpsi.Pt() / jet.getChargedPt(idx)));
                  else if(rbFit==2) rbWgtJPsi=(rbWgts[JPsi]->Eval(jpsi.Pt() / jet.getChargedPt(idx)));
                  if(jet.getChargedPt(idx) == 0) rbWgtJPsi=1.;
                  runBCDEF.SetRbWgt(rbWgtJPsi, JPsi);
                  runGH.SetRbWgt(rbWgtJPsi, JPsi);
                }
                std::vector<pfTrack> tmp_cands = { muTracks[i],muTracks[j] };
                float sf(1.),tWgt(1.);
                /*
                if(!isData) {
                  sf *= sqrt(customSF(tmp_cands[0], "BCDEF"));
                  sf *= sqrt(customSF(tmp_cands[1], "BCDEF"));
                }
                tWgt=runBCDEF.getSF();
                runBCDEF.SetSFs(tWgt*sf);
                runBCDEF.SetPiTrk(sf);
                treeBCDEF.SetPiTrk(sf);
                treeBCDEF.SetSFs(tWgt*sf);
                sf=pisfG;
                tWgt=1.;
                if(!isData) {
                  sf *= sqrt(customSF(tmp_cands[0], "GH"));
                  sf *= sqrt(customSF(tmp_cands[1], "GH"));
                }
                tWgt=runGH.getSF();
                runGH.SetSFs(tWgt*sf);
                runGH.SetPiTrk(sf);
                treeGH.SetSFs(tWgt*sf);
                treeGH.SetPiTrk(sf);
                */
                bool keepRunB(true);
                if(!isData) {
                  /*
                  if(shouldDrop[muTracks[i].getIdx()]<0) keepRunB=false;
                  if(shouldDrop[muTracks[j].getIdx()]<0) keepRunB=false;
                  */
                  if(piSFB[muTracks[i].getIdx()]<0) keepRunB=false;
                  if(piSFB[muTracks[j].getIdx()]<0) keepRunB=false;
                }
                if(keepRunB) {
                runBCDEF.Fill(tmp_cands, leptons, jet, chTag, "jpsi");
                treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "jpsi", ev.event);//, frag);
                treeBCDEF.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);
                }
                runGH.Fill(tmp_cands, leptons, jet, chTag, "jpsi");
                treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "jpsi", ev.event);//, frag);
                treeGH.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);


                if(!isData && pfmuMatched.size() > 1) { //save gen-matched J/Psi
                  runBCDEF.Fill(pfmuMatched, leptons, jet, chTag, "gjpsi");
                  runGH.Fill(pfmuMatched, leptons, jet, chTag, "gjpsi");
                }
                if(!isData && pfmuReject.size() > 1) { //save gen-unmatched J/Psi
                  runBCDEF.Fill(pfmuReject, leptons, jet, chTag, "rgjpsi");
                  runGH.Fill(pfmuReject, leptons, jet, chTag, "rgjpsi");
                }

                runBCDEF.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag, "jpsi");
                runGH.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag, "jpsi");
                //Only run once (i.e. in first jet from collection)
                if(&jet == &kJetsVec.front()) {
                  runBCDEF.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag, "jpsi");
                  runGH.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag, "jpsi");
                }
                runBCDEF.SetSFs(tWgt);
                runGH.SetSFs(tWgt);
                treeBCDEF.SetSFs(tWgt);
                treeGH.SetSFs(tWgt);
                runBCDEF.SetPiTrk(1./sf);
                runGH.SetPiTrk(1./sf);
                treeBCDEF.SetPiTrk(1./sf);
                treeGH.SetPiTrk(1./sf);
              }
              //muTracks.clear();
            } //end j
          } //end i
          //delete[] shouldDrop;
          //evch.nj++;
        }
        if(evch.nmeson>0) cht->Fill();
        evch = {}; //reset just in case (to avoid duplicates)
        if(debug) cout << "J/Psi DONE" << endl;
      }
      //end better J/Psi

      if(htsum<180) continue;
      //Better D meson
      //const float gMassMu(0.1057),gMassK(0.4937),gMassPi(0.1396);
      evch.nj=0;
      if(kalman.isMesonEvent()) { //FIXME can event have BOTH D and J/Psi?
        //copy, not ref, as sum ch is changed per D^0
        for(auto &jet : kJetsVec) {
          //jet.setJchCorrection(j_pt_corrB[jet.getJetIndex()]);
          //if(jet.getPt()>150) continue;
          //if(jet.getChargedPt()>75) continue;
          vector<pfTrack> piTracks,muTracks,piSoftTracks;
          /*
          size_t tmax = 4;
          tmax = jet.getTracks().size() >= tmax ? tmax : jet.getTracks().size();
          */
          for(auto &track : jet.getTracks()) {
            //Only save up to first 4 hardest tracks (sorted by pT already)
            //Only save up to first 4 hardest tracks (sorted by pT already)
            if(abs(track.getMotherId())!=421 && abs(track.getMotherId())!=42113 && abs(track.getMotherId())!=413) continue;
            if(abs(track.getMotherId())==413 && abs(track.getPdgId())==211) piSoftTracks.push_back(track); //save soft pions for D* separately
            //if(abs(track.getPdgId())==211) { piTracks.push_back(track); } //pi and K for D^0 and D*
            else if(abs(track.getMotherId())==42113 && abs(track.getPdgId())==13) { track.setMass(gMassMu); muTracks.push_back(track); } //mu for D^0 + mu (flavor tagging)
            else if(abs(track.getMotherId())==421 && abs(track.getPdgId())==211) { piTracks.push_back(track); } //pi and K for D^0
            //if(abs(track.getPdgId())==13) { cout << endl << ev.event << ": " << track.Pt() << " " << track.Eta() << " " << track.Phi() <<  " " << ev.k_mass[0] << endl; }
            //if(abs(track.getPdgId())==13) { cout << endl << ev.event << ": " << ev.k_pf_pt[0] << " " << ev.k_pf_eta[0] << " " << ev.k_pf_phi[0] <<  " " << ev.k_mass[0] << endl; }
            //if(abs(track.getPdgId())==13) { cout << endl << ev.event << ": " << ev.k_pf_pt[1] << " " << ev.k_pf_eta[1] << " " << ev.k_pf_phi[1] <<  " " << ev.k_mass[1] << endl; }
          }
          if(piTracks.size()<2) continue;
          size_t tmax = 4;
          tmax = piTracks.size() >= tmax ? tmax : piTracks.size();

          std::vector<pfTrack> pfMatched, pfReject, pfMuMatched;
          //Gen-matching
          if(!isData) {
            std::vector<pfTrack> genTracks;
            std::vector<pfTrack> genMuTracks;
            for(int ig = 0; ig < ev.ngpf; ig++) {
              //if(!isJPsiEvent) continue;
              //if(abs(ev.gmeson_id[ig])!=443) continue; //JPsi only
              if(abs(ev.gpf_id[ig])!=13 && abs(ev.gpf_id[ig])!=211) continue;
              TLorentzVector gen;
              gen.SetPtEtaPhiM(ev.gpf_pt[ig], ev.gpf_eta[ig], ev.gpf_phi[ig], gMassMu);
              genTracks.push_back(pfTrack(gen,0,0,0,0,ev.gpf_id[ig],3,true));
              genTracks.back().setMother(ev.gpf_mother[ig]); //daug -> mother id -> mother ttbar
              //cout << genTracks.back().getGenT() << endl;
              if(abs(ev.gpf_id[ig])==13) genMuTracks.push_back(pfTrack(gen,0,0,0,0, ev.gpf_id[ig],3,true));
            }
            //sort GEN tracks
            sort(genTracks.begin(), genTracks.end(),
                 [] (pfTrack a, pfTrack b) { return a.Pt() > b.Pt(); } );

            for(auto & it : piTracks) { //FIXME reference might not work
              double dR = 0.3; //initial dR
              int best_idx = -1;
              for(auto & itg : genTracks) {
                if(it.getPdgId() != itg.getPdgId()) continue; //insure ID and charge
                if(it.getVec().DeltaR(itg.getVec())>dR) continue; //find dR
                if(((it.Pt()-itg.Pt())/it.Pt())>0.10) continue; //gen and reco less than 10% difference
                dR = it.getVec().DeltaR(itg.getVec());
                best_idx = &itg - &genTracks[0]; //get index on current closest gen particle
              }
              if(best_idx<0) { //no gen track matched
                pfReject.push_back(it);
              }
              else {
                //it.setGenT(genMuTracks[best_idx].getGenT());
                pfMatched.push_back(it);
                pfMatched.back().setMother(genTracks[best_idx].getMotherId());
                it.setGenIdx(pfMatched.size()-1);
                genTracks.erase(genTracks.begin() + best_idx); //remove gen track so it cannot be matched again!
              }
            }
            for(auto & it : muTracks) { //FIXME reference might not work
              double dR = 0.3; //initial dR
              int best_idx = -1;
              for(auto & itg : genMuTracks) {
                if(it.getPdgId() != itg.getPdgId()) continue; //insure ID and charge
                if(it.getVec().DeltaR(itg.getVec())>dR) continue; //find dR
                //if(((it.Pt()-itg.Pt())/it.Pt())>0.10) continue; //gen and reco less than 10% difference
                dR = it.getVec().DeltaR(itg.getVec());
                best_idx = &itg - &genMuTracks[0]; //get index on current closest gen particle
              }
              if(best_idx<0) { //no gen track matched
                pfReject.push_back(it);
              }
              else {
                //it.setGenT(genMuTracks[best_idx].getGenT());
                pfMuMatched.push_back(it);
                pfMuMatched.back().setMother(genMuTracks[best_idx].getMotherId());
                it.setGenIdx(pfMuMatched.size()-1);
                genMuTracks.erase(genMuTracks.begin() + best_idx); //remove gen track so it cannot be matched again!
              }
            }
          }
          Jet genJet;
          int genIdx(-1);
          for(auto & gjet : genJetsVec) {
            if(gjet.getVec().DeltaR(jet.getVec())>0.1) continue; //find good dR match
            if(gjet.getVec().DeltaR(jet.getVec()) > genJet.getVec().DeltaR(jet.getVec())) continue; //find tighter dR match
            genJet = gjet;
            genIdx = genJet.getJetIndex();
            break;
          }
          std::vector<float> frag;
          if(genIdx>-1)
            frag = {ev.up[genJet.getJetIndex()],ev.central[genJet.getJetIndex()],ev.down[genJet.getJetIndex()]};
            //frag = {ev.xb[genJet.getJetIndex()]};//,ev.peterson[genJet.getJetIndex()],ev.up[genJet.getJetIndex()],ev.central[genJet.getJetIndex()],ev.down[genJet.getJetIndex()]};

          sort(piTracks.begin(), piTracks.end(),
               [] (pfTrack a, pfTrack b) { return a.M() < b.M(); } );
          //only loop over i<j since mass is assigned in Kalman filter
          //memoize pi/K tracks
          int *trkIdx = new int[ev.npf]();
          //store all PF tracks per jet once
	  for(int ipf = 0; ipf < ev.npf; ipf++) {
            if(ev.pf_j[ipf] != jet.getJetIndex()) continue;
            if(ev.pf_fromPV[ipf]<2) continue;
            evch.pf_pt[ipf]  = ev.pf_pt[ipf];
            evch.pf_eta[ipf] = ev.pf_eta[ipf];
            evch.pf_phi[ipf] = ev.pf_phi[ipf];
            evch.pf_m[ipf]   = ev.pf_m[ipf];
            evch.pf_id[ipf]  = ev.pf_id[ipf];
            evch.npf++;
          }
          //recompute probability for all D^0 pi/K once from pf_eff_pi
          /*
          delete[] piSFB;
          piSFB = new int[ev.npf]();
          applyTrackingEfficiencySF(ev, pf_eff_pi, piSFB);
          */
          int *shouldDrop = new int[piTracks.size()]();
          int *shouldDropK = new int[piTracks.size()]();
          applyTrackingEfficiencySF(piTracks, pf_eff_pi, shouldDrop, passBit(runSysts,PI_BIT));
          size_t *isUsedB = new size_t[piTracks.size()]();
          size_t *isUsedG = new size_t[piTracks.size()]();
          //applyTrackingEfficiencySF(piTracks, pf_eff_k, shouldDropK, passBit(runSysts,PI_BIT));
          for(size_t i = 0; i < piTracks.size(); i++) {
            if(i > tmax) break;
            for(size_t j = i+1; j < piTracks.size(); j++) { //i<j b/c Kalman filter already loops over all i,j
              if(j > tmax) break;
              if(i==j) continue;
              if(abs(piTracks[i].getMotherId())!=421) continue;
              if(abs(piTracks[j].getMotherId())!=421) continue;
              //ensure same Kalman mass (sorting tracks by pT might mess up order in which Kalman masses were saved)
              if(piTracks[i].getKalmanMass() != piTracks[j].getKalmanMass()) continue;
              //if(((piTracks[i].getKalmanMass() - piTracks[j].getKalmanMass()) / piTracks[i].getKalmanMass()) > 0.01) continue;
              //Set mass assmumption
              //piTracks[i].setMass(gMassPi); piTracks[j].setMass(gMassK);
              //Mass already set by Kalman filter
              /*
              sort(pfMatched.begin(), pfMatched.end(),
                   [] (pfTrack a, pfTrack b) { return a.M() < b.M(); } );
              sort(pfMuMatched.begin(), pfMuMatched.end(),
                   [] (pfTrack a, pfTrack b) { return a.M() < b.M(); } );
              */
              //std::cout << piTracks[i].getKalmanMass() << " " << piTracks[j].getKalmanMass() << " " << (piTracks[i].getKalmanMass() - piTracks[j].getKalmanMass()) / piTracks[i].getKalmanMass() << std::endl;
              //Check masses from Kalman class
              if(piTracks[i].M()!=gMassPi) continue;
              if(piTracks[j].M()!=gMassK) continue;
              if(piTracks[i].Pt() < 5) continue;
              //if(piTracks[i].Pt() < 12) continue; //norm GEN vs norm unmatched pi cross at 12 GeV
              if(piTracks[j].Pt() < 1) continue;
              bool cuts(true),keepRunB(true),runG(true);
              cuts &= abs(piTracks[i].Eta()) < 1.5;
              cuts &= abs(piTracks[j].Eta()) < 1.5;
              //cuts &= abs(piTracks[i].getDxy()/piTracks[i].getDxyE()) > 0.5;
              //cuts &= abs(piTracks[j].getDxy()/piTracks[j].getDxyE()) > 0.5;
              /*
              if(abs(piTracks[i].Eta()) > 1.) continue;
              if(abs(piTracks[j].Eta()) > 1.) continue;
              if(abs(piTracks[i].getDxy()/piTracks[i].getDxyE()) < 0.5) continue;
              if(abs(piTracks[j].getDxy()/piTracks[j].getDxyE()) < 0.5) continue;
              */
              TLorentzVector d0 = piTracks[i].getVec() + piTracks[j].getVec();
              /*
              float cosDjet = (d0.Vect()).Dot(jet.getVec().Vect())/(jet.getVec().P() * d0.P());
              if(cosDjet < 0.99) continue;
              */
              float mass12 = (piTracks[i].getVec()+piTracks[j].getVec()).M();
              //istd::cout << piTracks[i].getKalmanMass() << " " << piTracks[j].getKalmanMass() << " " << mass12 << std::endl;
              cuts &= (mass12>1.7 && mass12<2.0);
              //if (mass12>1.65 && mass12<2.0) {
              //if(D0p4.Pt()<20) cuts=false; //FIXME
              //custom SFs
              //float sfB(1.),sfG(1.),tWgt(1.);
              //float sfB(pisfB),sfG(pisfG),tWgt(1.);
              if(trkIdx[piTracks[i].getIdx()]==0) { //only add D^0 pi/K once
                pi_pf_eff->Fill(piTracks[i].Eta(), piTracks[i].Pt());
                trkIdx[piTracks[i].getIdx()]++;
              }
              if(trkIdx[piTracks[j].getIdx()]==0) { //only add D^0 pi/K once
                pi_pf_eff->Fill(piTracks[j].Eta(), piTracks[j].Pt());
                trkIdx[piTracks[j].getIdx()]++;
              }
              if(!isData) {
                //skip if dropped by pi tracking eff
                if(shouldDrop[i]<0) keepRunB=false;
                if(shouldDrop[j]<0) keepRunB=false;
                /*
                if(piSFB[piTracks[i].getIdx()]<0) keepRunB=false;
                if(piSFB[piTracks[j].getIdx()]<0) keepRunB=false;
                */
                //add skipped PF track back into sumch
                if(sumChBidx[piTracks[i].getIdx()] == 0 && sumChBidx[piTracks[j].getIdx()] != 0) {
                  float tmpSumCh(jet.getChargedPt(0));
                  //if(tmpSumCh==0. && jet.getTracks().size()>2) keepRunB=false; //all others were dropped
                  tmpSumCh += piTracks[i].Pt();
                  jet.updateChargedPt(tmpSumCh, pt_chargedG[jet.getJetIndex()]);
                  sumChBidx[piTracks[i].getIdx()]++;
                }
                if(sumChBidx[piTracks[j].getIdx()] == 0 && sumChBidx[piTracks[i].getIdx()] != 0) {
                  float tmpSumCh(jet.getChargedPt(0) + piTracks[j].Pt());
                  //if(tmpSumCh==0.) keepRunB=false; //all others were dropped
                  tmpSumCh += piTracks[j].Pt();
                  jet.updateChargedPt(tmpSumCh, pt_chargedG[jet.getJetIndex()]);
                  sumChBidx[piTracks[j].getIdx()]++;
                }
                if(rbFit && !isData) {
                  float rbWgtD0(1.);
                  int idx = (runPeriod.Contains("BCDEF") ? 0 : 1);
                  if(rbFit==1) rbWgtD0=(rbWgt->Eval(d0.Pt() / jet.getChargedPt(idx)));
                  else if(rbFit==2) rbWgtD0=(rbWgts[D0]->Eval(d0.Pt() / jet.getChargedPt(idx)));
                  if(jet.getChargedPt(idx) == 0.) rbWgtD0=1.;
                  runBCDEF.SetRbWgt(rbWgtD0, D0);
                  runGH.SetRbWgt(rbWgtD0, D0);
                }
                float ptsf = 1., etasf = 1.;
                //customSF(piTracks[i], "BCDEF", ptsf, etasf, pf_eff_pi); //event weight from pi track only
                //customSF(piTracks[j], "BCDEF", ptsf, etasf, pf_eff_pi); //trk weight for K track only
                //piTracks[i].setEtaCorrection(etasf);
                //piTracks[i].setPtCorrection(ptsf);
                if(ptsf > 0) { //pT shape weight only
                allPlots["piwgt_BCDEF"]->Fill(0.,1.0);
                allPlots["piwgt_BCDEF"]->Fill(1.,ptsf);
                }
                //piTracks[j].setPtCorrection(ptsf);
                //sfB /= ptsf;
                if(ptsf > 0) { //pT shape weight only
                allPlots["piwgt_BCDEF"]->Fill(0.,1.0);
                allPlots["piwgt_BCDEF"]->Fill(1.,ptsf);
                }
                /*
                //skip if only added by othe epoch
                if(piSFG[piTracks[i].getIdx()]>0 && piSFB[piTracks[i].getIdx()]==1) sfB = 0.;
                if(piSFG[piTracks[j].getIdx()]>0 && piSFB[piTracks[j].getIdx()]==1) sfB = 0.;
                */
                //tmp_candsB[0].getVec() *= sf;
                //sf *= customSF(tmp_cands[1].Pt(), "BCDEF");
                //allPlots["piwgt_BCDEF"]->Fill(0.,1.0);
                //allPlots["piwgt_BCDEF"]->Fill(1.,sfB);
                //tmp_candsB[1].getVec() *= sf;
              }
              /*
              tWgt=runBCDEF.getSF();
              runBCDEF.SetSFs(tWgt*sfB);
              runBCDEF.SetPiTrk(sfB);
              treeBCDEF.SetSFs(tWgt*sfB);
              treeBCDEF.SetPiTrk(sfB);
              sfG=1.;
              tWgt=1.;
              */
              if(!isData && 0) {
                customSF(piTracks[i], "GH");
                //skip if dropped by pi tracking eff
                /*
                if(piSFG[piTracks[i].getIdx()]<0) sfG = 0.;
                if(piSFG[piTracks[j].getIdx()]<0) sfG = 0.;
                //skip if only added by othe epoch
                if(piSFB[piTracks[i].getIdx()]>0 && piSFG[piTracks[i].getIdx()]==1) sfG = 0.;
                if(piSFB[piTracks[j].getIdx()]>0 && piSFG[piTracks[j].getIdx()]==1) sfG = 0.;
                */
                //tmp_candsG[0].getVec() *= sfG;
                //sfG *= customSF(tmp_cands[1].Pt(), "GH");
                /*
                allPlots["piwgt_GH"]->Fill(0.,1.0);
                allPlots["piwgt_GH"]->Fill(1.,sfG);
                */
                //tmp_candsG[1].getVec() *= sfG;
              }
              std::vector<pfTrack> tmp_cands = { piTracks[i],piTracks[j] };
              TLorentzVector D0p4 = tmp_cands[0].getVec() + tmp_cands[1].getVec();
              /*
              tWgt=runGH.getSF();
              runGH.SetSFs(tWgt*sfG);
              runGH.SetPiTrk(sfG);
              treeGH.SetSFs(tWgt*sfG);
              treeGH.SetPiTrk(sfG);
              */
              if (cuts) {
                //if(debug) cout << pfmuCands[0].Pt() << " " << pfmuCands[0].Eta() << " " << pfmuCands[0].Phi() << " " << gMassMu << endl;
                //if(debug) cout << pfmuCands[1].Pt() << " " << pfmuCands[1].Eta() << " " << pfmuCands[1].Phi() << " " << gMassMu << endl;
                if(debug) cout << mass12 << endl << endl;
                if(debug) cout << "D0 found" << endl;
                if(debug && (mass12>1.8 && mass12<1.9)) cout << "and it's good!" << endl;
                allPlots["massD0"+chTag]->Fill(mass12,wgt);
	        allPlots["massD0_all"]->Fill(mass12,wgt);
                allPlots["mt"+chTag+"_meson"]->Fill(mt,wgt);
	        allPlots["mt_all_meson"]->Fill(mt,wgt);
                if(lightJetsVec.size()>1) {
                  allPlots["massD0_1j"+chTag]->Fill(mass12,wgt);
                  allPlots["massD0_1j_all"]->Fill(mass12,wgt);
                }
                if(lightJetsVec.size()>2) {
                  allPlots["massD0_2j"+chTag]->Fill(mass12,wgt);
                  allPlots["massD0_2j_all"]->Fill(mass12,wgt);
                }
                if(lightJetsVec.size()>3) {
                  allPlots["massD0_3j"+chTag]->Fill(mass12,wgt);
                  allPlots["massD0_3j_all"]->Fill(mass12,wgt);
                }
                if(lightJetsVec.size()>4) {
                  allPlots["massD0_4j"+chTag]->Fill(mass12,wgt);
                  allPlots["massD0_4j_all"]->Fill(mass12,wgt);
                }

                bool keepRunG(true);
                //Check if both pi and K already used in loop (remove ij -> ji duplicates)
                if(isUsedB[i] == j) keepRunB=false;
                if(isUsedG[i] == j) keepRunG=false;
                if(isUsedB[i] == j) std::cout << "DUPLICATE" << std::endl;
                if(isUsedG[i] == j) std::cout << "DUPLICATE" << std::endl;

                if(keepRunB) runBCDEF.Fill(tmp_cands, leptons, jet, chTag, "meson");
                if(keepRunB) treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event);//, tmp_match, frag, genJet);
                if(keepRunG) runGH.Fill(tmp_cands, leptons, jet, chTag, "meson");
                if(keepRunG) treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event);//, tmp_match, frag, genJet);
                if(keepRunB) isUsedB[i] = j;
                isUsedG[i] = j;
                if(isTTbar) {
                  if(keepRunB) runBCDEF.FillGen(tops, chTag, "meson");
                  if(keepRunG) runGH.FillGen(tops, chTag, "meson");
                }
                std::vector<pfTrack> tmp_match;
                if(!isData && pfMatched.size() >1 && 0) {
                  for(size_t i = 0; i < pfMatched.size(); i++) {
                    for(size_t j = 0; j < pfMatched.size(); j++) {
                      tmp_match = { pfMatched[piTracks[i].getGenIdx()],pfMatched[piTracks[j].getGenIdx()] };
                      //void CharmTree::Fill(CharmEvent_t &ev_, std::vector<pfTrack>& pfCands, Leptons lep, Jet jet, TString, TString name, int event, std::vector<pfTrack> genMatch, std::vector<float> frag) 

                      if(keepRunB) runBCDEF.Fill(tmp_match, leptons, jet, chTag, "gmeson");
                      if(keepRunG) runGH.Fill(tmp_match, leptons, jet, chTag, "gmeson");
                    }
                  }
                }
                if(!isData && pfReject.size() > 1 && 0) { //save gen-unmatched J/Psi
                  for(size_t i = 0; i < pfReject.size(); i++) {
                    for(size_t j = 0; j < pfReject.size(); j++) {
                      tmp_match = { pfReject[piTracks[i].getGenIdx()],pfReject[piTracks[j].getGenIdx()] };
                      if(keepRunB) runBCDEF.Fill(pfReject, leptons, jet, chTag, "rgmeson");
                      if(keepRunG) runGH.Fill(pfReject, leptons, jet, chTag, "rgmeson");
                    }
                  }
                }
                /*
                if(genIdx>-1) {  
                  treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event, tmp_match, frag, genJet);
                  if(keepRunG) treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event, tmp_match, frag, genJet);
                }
                else {
                  treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event);//, tmp_match, frag, genJet);
                  if(keepRunG) treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event);//, tmp_match, frag, genJet);
                }
                */

                if(keepRunB) runBCDEF.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag, "meson");
                if(keepRunB) treeBCDEF.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);
                if(keepRunG) runGH.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag, "meson");
                if(keepRunG) treeGH.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);
	        allPlots["HT"+chTag+"_meson"]->Fill(htsum,wgt);
	        allPlots["HT_all_meson"]->Fill(htsum,wgt);
	        allPlots["H"+chTag+"_meson"]->Fill(hsum,wgt);
	        allPlots["H_all_meson"]->Fill(hsum,wgt);
	        if(htbsum) allPlots["HTb"+chTag+"_meson"]->Fill(htbsum,wgt);
	        if(htbsum) allPlots["HTb_all_meson"]->Fill(htbsum,wgt);
	        allPlots["nbj"+chTag+"_meson"]->Fill(nbj,wgt);
	        allPlots["nbj_all_meson"]->Fill(nbj,wgt);
                //Only run once (i.e. in first jet from collection)
                if(&jet == &kJetsVec.front()) {
                  if(keepRunB) runBCDEF.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag, "meson");
                  if(keepRunG) runGH.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag, "meson");
                }
              //D^0 + mu for flavor tagging
              } //end extra cuts within D^0 window
              if(muTracks.size()<1) continue;
              for(auto &track : muTracks) {
                if(abs(track.getMotherId())!=42113) continue;
                if(piTracks[j].getKalmanMass() != track.getKalmanMass()) continue;
                //if(((piTracks[i].getKalmanMass() - track.getKalmanMass()) / piTracks[i].getKalmanMass()) > 0.01) continue;
                //if(!track.highPurity()) continue; //only use high purity tracks FIXME
                //if(!track.trackerMuon() && !track.globalMuon()) continue;
                //if(track.Pt() < 3.0) continue;
                if(debug) cout << "third lepton found" << endl;
                if(piTracks[j].charge()*track.charge()<0) continue; //PDGID 13 is NEGATIVE mu but charge function in pfTrack acconts for this
                //Kaon and lepton have same charge (e.g. b^-1/3 -> c^+2/3 W^- -> c^+2/3 l^- nubar)
                //correct mass assumption
                if(debug) cout << "correct mass assumption" << endl;
                TLorentzVector d0mu = piTracks[i].getVec() + piTracks[j].getVec() + track.getVec();
                if(rbFit && !isData) {
                  float rbWgtD0mu(1.);
                  int idx = (runPeriod.Contains("BCDEF") ? 0 : 1);
                  if(rbFit==1) rbWgtD0mu=(rbWgt->Eval(d0mu.Pt() / jet.getChargedPt(idx)));
                  else if(rbFit==2) rbWgtD0mu=(rbWgts[D0]->Eval(d0mu.Pt() / jet.getChargedPt(idx)));
                  if(jet.getChargedPt(idx) == 0) rbWgtD0mu=1.;
                  runBCDEF.SetRbWgt(rbWgtD0mu, D0mu);
                  runGH.SetRbWgt(rbWgtD0mu, D0mu);
                }
                std::vector<pfTrack> tmp_cands = { piTracks[i],piTracks[j],track };
                //tmp_cands.push_back( track );
                //keepRunB=true;
                if(!isData && piSFB[track.getIdx()]<0) keepRunB=false;
                if(keepRunB) runBCDEF.Fill(tmp_cands, leptons, jet, chTag, "meson");
                runGH.Fill(tmp_cands, leptons, jet, chTag, "meson");
                std::vector<pfTrack> tmp_match;
                /*
                if(!isData && pfMatched.size() > 1 && pfMuMatched.size() > 0) { //save gen-matched J/Psi
                  //std::vector<pfTrack> tmp_matched_cands = { pfMatched[0],pfMatched[1],pfMuMatched[0] };
                  tmp_match = { pfMatched[piTracks[i].getGenIdx()],pfMatched[piTracks[j].getGenIdx()],pfMuMatched[track.getGenIdx()] };
                  runBCDEF.Fill(tmp_match, leptons, jet, chTag, "gmeson");
                  runGH.Fill(tmp_match, leptons, jet, chTag, "gmeson");
                }
                */
                /*
                if(!isData) {
                  treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event, tmp_match);
                  treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event, tmp_match);
                }
                else {
                treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "meson");
                treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "meson");
                }
                */
                /*
                if(genIdx>-1) {
                treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event, tmp_match, frag, genJet);
                treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event, tmp_match, frag, genJet);
                }
                else {
                treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event, tmp_match);
                treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event, tmp_match);
                }
                treeBCDEF.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);
                treeGH.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);
                */
                treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event, tmp_match);
                treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "meson", ev.event, tmp_match);
              }
              if(piSoftTracks.size()<1) continue;
              for(auto &track : piSoftTracks) {
                if(abs(track.getMotherId())!=413) continue;
                //if(piTracks[j].getKalmanMass() != track.getKalmanMass()) continue;
                //if(fabs(mass12-1.864) > 0.05) continue; // tighter mass window cut
                if(fabs(mass12-1.864) > 0.1) continue; // tighter mass window cut
                if( piTracks[j].charge() == track.charge() ) continue;
                  // Kaon and pion have opposite charges
                  // I.e. correct mass assumption
                //std::vector<pfTrack> tmp_cands = { piTracks[i],piTracks[j],track };
                tmp_cands.push_back( track );
                runBCDEF.Fill(tmp_cands, leptons, jet, chTag, "meson");
                runGH.Fill(tmp_cands, leptons, jet, chTag, "meson");
                treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "meson");
                treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "meson");
                treeBCDEF.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);
                treeGH.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);
              } //end D*
              /*
              runBCDEF.SetPiTrk(1./sfB);
              runGH.SetPiTrk(1./sfG);
              treeBCDEF.SetSFs(tWgt);
              runGH.SetSFs(tWgt);
              treeBCDEF.SetSFs(tWgt);
              treeGH.SetSFs(tWgt);
              treeBCDEF.SetPiTrk(1./sfB);
              treeGH.SetPiTrk(1./sfG);
              */
            } //end D^0 j
          } //end D^0 i
          //evch.nj++;
          delete[] trkIdx;
          delete[] shouldDrop;
          delete[] shouldDropK;
        } //end jets loop
        if(evch.nmeson>0) cht->Fill();
        evch = {}; //reset just in case (to avoid duplicates)
        if(debug) cout << "D meosn DONE" << endl;
      }
      //end better D meson

      //test for D^*
      /*
      evch.nj=0;
      if(kalman.isMesonEvent()) { //FIXME can event have BOTH D and J/Psi?
        for(auto &jet : kJetsVec) {
          vector<pfTrack> piTracks,muTracks,piSoftTracks;
          for(auto &track : jet.getTracks()) {
            //Only save up to first 4 hardest tracks (sorted by pT already)
            //Only save up to first 4 hardest tracks (sorted by pT already)
            if(abs(track.getMotherId())!=421 && abs(track.getMotherId())!=42113 && abs(track.getMotherId())!=413) continue;
            if(abs(track.getMotherId())==413 && abs(track.getPdgId())==211) piSoftTracks.push_back(track); //save soft pions for D* separately
            if(abs(track.getPdgId())==211) { piTracks.push_back(track); } //pi and K for D^0 and D*
            else if(abs(track.getMotherId())==42113 && abs(track.getPdgId())==13) { track.setMass(gMassMu); muTracks.push_back(track); } //mu for D^0 + mu (flavor tagging)
            //else if(abs(track.getMotherId())==421 && abs(track.getPdgId())==211) { piTracks.push_back(track); } //pi and K for D^0
            //if(abs(track.getPdgId())==13) { cout << endl << ev.event << ": " << track.Pt() << " " << track.Eta() << " " << track.Phi() <<  " " << ev.k_mass[0] << endl; }
            //if(abs(track.getPdgId())==13) { cout << endl << ev.event << ": " << ev.k_pf_pt[0] << " " << ev.k_pf_eta[0] << " " << ev.k_pf_phi[0] <<  " " << ev.k_mass[0] << endl; }
            //if(abs(track.getPdgId())==13) { cout << endl << ev.event << ": " << ev.k_pf_pt[1] << " " << ev.k_pf_eta[1] << " " << ev.k_pf_phi[1] <<  " " << ev.k_mass[1] << endl; }
          }
          if(piTracks.size()<2) continue;
          size_t tmax = 4;
          tmax = piTracks.size() >= tmax ? tmax : piTracks.size();

          std::vector<pfTrack> pfMatched, pfReject, pfMuMatched;
          //Gen-matching

          sort(piTracks.begin(), piTracks.end(),
               [] (pfTrack a, pfTrack b) { return a.M() < b.M(); } );
          //only loop over i<j since mass is assigned in Kalman filter
          for(size_t i = 0; i < piTracks.size(); i++) {
            if(i > tmax) break;
            for(size_t j = i+1; j < piTracks.size(); j++) {
              if(j > tmax) break;
              if(i==j) continue;
              if(abs(piTracks[i].getMotherId())!=421) continue;
              if(abs(piTracks[j].getMotherId())!=421) continue;
              //ensure same Kalman mass (sorting tracks by pT might mess up order in which Kalman masses were saved)
              if(piTracks[i].getKalmanMass() != piTracks[j].getKalmanMass()) continue;
              if(((piTracks[i].getKalmanMass() - piTracks[j].getKalmanMass()) / piTracks[i].getKalmanMass()) > 0.01) continue;
              //Check masses from Kalman class
              if(piTracks[i].M()!=gMassPi) continue;
              if(piTracks[j].M()!=gMassK) continue;
              if(piTracks[i].Pt() < 5) continue;
              if(piTracks[j].Pt() < 1) continue;
              float mass12 = (piTracks[i].getVec()+piTracks[j].getVec()).M();
              //istd::cout << piTracks[i].getKalmanMass() << " " << piTracks[j].getKalmanMass() << " " << mass12 << std::endl;
              if (mass12>1.65 && mass12<2.0) {
                if(debug) cout << mass12 << endl << endl;
                if(debug) cout << "D0 found" << endl;
                if(debug && (mass12>1.8 && mass12<1.9)) cout << "and it's good!" << endl;
              }
              if(piTracks.size()<1) continue;
              for(auto &track : piSoftTracks) {
                if(abs(track.getMotherId())!=413) continue;
                if(piTracks[j].getKalmanMass() != track.getKalmanMass()) continue;
                if(((piTracks[i].getKalmanMass() - track.getKalmanMass()) / piTracks[i].getKalmanMass()) > 0.01) continue;
                //if(fabs(mass12-1.864) > 0.05) continue; // tighter mass window cut
                if(fabs(mass12-1.864) > 0.1) continue; // tighter mass window cut
                if( piTracks[j].charge() == track.charge() ) continue;
                  // Kaon and pion have opposite charges
                  // I.e. correct mass assumption
                std::vector<pfTrack> tmp_cands = { piTracks[i],piTracks[j],track };
                runBCDEF.Fill(tmp_cands, leptons, jet, chTag, "meson");
                runGH.Fill(tmp_cands, leptons, jet, chTag, "meson");
                treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "meson");
                treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "meson");
                //FIXME
                treeBCDEF.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);
                treeGH.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);
                //FIXME
              } //end D*
            } //end D^0 j
          } //end D^0 i
          //evch.nj++;
        } //end jets loop
        if(evch.nmeson>0) cht->Fill();
        evch = {}; //reset just in case (to avoid duplicates)
        if(debug) cout << "D meosn DONE" << endl;
      }
      */
      /*
      evch.nj=0;
      if(kalman.isMesonEvent()) { //FIXME can event have BOTH D and J/Psi?
        for(auto &jet : kJetsVec) {
          vector<pfTrack> piTracks,muTracks,piSoftTracks;
          for(auto &track : jet.getTracks()) {
            if(abs(track.getMotherId())!=421 && abs(track.getMotherId())!=42113 && abs(track.getMotherId())!=413) continue;
            if(abs(track.getMotherId())==413 && abs(track.getPdgId())==211) piSoftTracks.push_back(track); //save soft pions for D* separately
            if(abs(track.getPdgId())==211) { piTracks.push_back(track); } //pi and K for D^0 and D*
          }
          if(piTracks.size()<2) continue;
          size_t tmax = 4;
          tmax = piTracks.size() >= tmax ? tmax : piTracks.size();

          std::vector<pfTrack> pfMatched, pfReject, pfMuMatched;
          //Gen-matching

          sort(piTracks.begin(), piTracks.end(),
               [] (pfTrack a, pfTrack b) { return a.M() < b.M(); } );
          //only loop over i<j since mass is assigned in Kalman filter
          for(size_t i = 0; i < piTracks.size(); i++) {
            if(i > tmax) break;
            for(size_t j = i+1; j < piTracks.size(); j++) {
              if(j > tmax) break;
              if(i==j) continue;
              if(abs(piTracks[i].getMotherId())!=421) continue;
              if(abs(piTracks[j].getMotherId())!=421) continue;
              //ensure same Kalman mass (sorting tracks by pT might mess up order in which Kalman masses were saved)
              if(piTracks[i].getKalmanMass() != piTracks[j].getKalmanMass()) continue;
              if(((piTracks[i].getKalmanMass() - piTracks[j].getKalmanMass()) / piTracks[i].getKalmanMass()) > 0.01) continue;
              //Check masses from Kalman class
              if(piTracks[i].M()!=gMassPi) continue;
              if(piTracks[j].M()!=gMassK) continue;
              if(piTracks[i].Pt() < 5) continue;
              if(piTracks[j].Pt() < 1) continue;
              float mass12 = (piTracks[i].getVec()+piTracks[j].getVec()).M();
              if (mass12<1.65 || mass12>2.0) continue;
              if(piSoftTracks.size()<3) continue;
              for(auto &track : piSoftTracks) {
                if(abs(track.getMotherId())!=413) continue;
                if(piTracks[j].getKalmanMass() != track.getKalmanMass()) continue;
                if(((piTracks[i].getKalmanMass() - track.getKalmanMass()) / piTracks[i].getKalmanMass()) > 0.01) continue;
                //if(fabs(mass12-1.864) > 0.05) continue; // tighter mass window cut
                if(fabs(mass12-1.864) > 0.1) continue; // tighter mass window cut
                if( piTracks[j].charge() == track.charge() ) continue;
                  // Kaon and pion have opposite charges
                  // I.e. correct mass assumption
                std::vector<pfTrack> tmp_cands = { piTracks[i],piTracks[j],track };
                runBCDEF.Fill(tmp_cands, leptons, jet, chTag, "meson");
                runGH.Fill(tmp_cands, leptons, jet, chTag, "meson");
                treeBCDEF.Fill(evch, tmp_cands, leptons, jet, chTag, "meson");
                treeGH.Fill(evch, tmp_cands, leptons, jet, chTag, "meson");
                //FIXME
                treeBCDEF.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);
                treeGH.Fill(evch, ev.nvtx, htsum, stsum, ev.met_pt[0], lightJetsVec);
                //FIXME
              } //end D*
            } //end D^0 j
          } //end D^0 i
          //evch.nj++;
        } //end jets loop
        if(evch.nmeson>0) cht->Fill();
        evch = {}; //reset just in case (to avoid duplicates)
        if(debug) cout << "D meosn DONE" << endl;
      }
      */
              

      for(size_t ij = 0; ij < kJetsVec.size(); ij++) {

        //if(ij > 1) continue;
        if(ij == 0) {
          runBCDEF.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag, "csv");
          runGH.Fill(1, ev.nvtx, htsum, stsum, ev.met_pt[0], chTag, "csv");

          runBCDEF.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag, "csv");
          runGH.Fill(leptons,lightJetsVec,kJetsVec,allJetsVec, chTag, "csv");
          
          runBCDEF.Fill(leptons, chTag, "csv");
          runGH.Fill(leptons, chTag, "csv");

          //allPlots["nkj"+chTag+"_csv"]->Fill(1,wgt);
          allPlots["nevt"+chTag+"_csv"]->Fill(1,norm);
          allPlots["weight"+chTag+"_csv"]->Fill(wgt,norm);
          allPlots["norm"+chTag+"_csv"]->Fill(norm,norm);
          allPlots["nvtx"+chTag+"_csv"]->Fill(ev.nvtx,wgt);
          allPlots["j_pt"+chTag+"_csv"]->Fill(kJetsVec[0].getVec().Pt(),wgt);
          //allPlots["lj_pt"+chTag+"_csv"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
          allPlots["kj_pt"+chTag+"_csv"]->Fill(kJetsVec[0].getVec().Pt(),wgt);
        }

      }

      kJetsVec.clear();
      delete[] piSFB;
      delete[] sumChBidx;
      delete[] piSFG;
    } // end event

  //close input file
  f->Close();

  //save histos to file  
  TString selPrefix("");  
  if(flavourSplitting!=NOFLAVOURSPLITTING) selPrefix=Form("%d_",flavourSplitting);
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+selPrefix+baseName,"RECREATE");
  fOut->cd();
  if(debug) cout << "writing histograms" << endl;
  treeBCDEF.PdfWeights(evch, fOut);
  treeGH.PdfWeights(evch, fOut);

  for (auto& it : allPlots)  { 
    //if(debug) cout << it.second->GetName() << endl;
    //if(debug) cout << it.second->GetEntries() << endl;

    it.second->SetDirectory(fOut); it.second->Write(); 
    //fOut->cd();
  }
  for (auto& it : allPlots2D)  {
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : allPlots3D)  {
    it.second->SetDirectory(fOut); it.second->Write(); 
  }

  //restore run period for writing
  runBCDEF.CheckRunPeriod("BCDEF");
  runGH.CheckRunPeriod("GH");
  
  runBCDEF.Write();
  runGH.Write();

  pi_pf_eff->SetDirectory(fOut);
  pi_pf_eff->Write();
  if(!isData) {
    epoch = "GH";
    pi_pf_eff->SetName("pi_eta_pt_" + epoch);
    pi_pf_eff->SetTitle("pi_eta_pt_" + epoch);
    pi_pf_eff->Write();
  }
  delete pi_pf_eff;
  if(debug) cout << "writing histograms DONE" << endl;

  if(debug) cout << "writing tree" << endl;
  cht->Write();
  if(debug) cout << "writing tree DONE" << endl;

  if(debug) cout << "closing ROOT file" << endl;
  fOut->Close();

  delete jecParam;
  delete jecUnc;
  delete jres;
  delete jres_sf;
  delete cht;
}

