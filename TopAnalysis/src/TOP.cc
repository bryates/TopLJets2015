#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"
#include "TopLJets2015/TopAnalysis/interface/TOPWidth.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "TopLJets2015/TopAnalysis/interface/OtherFunctions.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "TopLJets2015/TopAnalysis/interface/rochcor2016.h"

#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

using namespace std;

//
/*
Float_t computeMT(TLorentzVector &a, TLorentzVector &b)
{
  return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
}


*/
//
void RunTop(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts,
                 TString era,
                 Bool_t debug=false)
{
  if(debug) cout << "in RunTop" << endl;

  bool isTTbar( filename.Contains("_TTJets") );
  
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
  if(ev.isData) runSysts=false;
  bool requireEletriggerOnly(false);
  if(ev.isData && filename.Contains("SingleElectron")) requireEletriggerOnly=true;
  bool requireMutriggerOnly(false);
  if(ev.isData && filename.Contains("SingleMuon")) requireMutriggerOnly=true;
  bool requireEETriggers(false);
  if(ev.isData && filename.Contains("DoubleEG"))       requireEETriggers=true;
  bool requireMMTriggers(false);
  if(ev.isData && filename.Contains("DoubleMuon"))     requireMMTriggers=true;
  bool requireEMTriggers(false);
  if(ev.isData && filename.Contains("MuonEG"))         requireEMTriggers=true;

  cout << "...producing " << outname << " from " << nentries << " events" << (runSysts ? " syst variations will be considered" : "") << endl;
  
  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  if(!ev.isData)
    {
      if(debug) cout << "loading pileup weight" << endl;
      TString puWgtUrl(era+"/pileupWgts.root");
      gSystem->ExpandPathName(puWgtUrl);
      TFile *fIn=TFile::Open(puWgtUrl);
      for(size_t i=0; i<3; i++)
	{
	  TString grName("pu_nom");
          if(i==1) grName="pu_down";
          if(i==2) grName="pu_up";
	  TGraph *puData=(TGraph *)fIn->Get(grName);
	  Float_t totalData=puData->Integral();
	  TH1 *tmp=(TH1 *)puTrue->Clone("tmp");
	  for(Int_t xbin=1; xbin<=tmp->GetXaxis()->GetNbins(); xbin++)
	    {
	      Float_t yexp=puTrue->GetBinContent(xbin);
	      Double_t xobs,yobs;
	      puData->GetPoint(xbin-1,xobs,yobs);
	      tmp->SetBinContent(xbin, yexp>0 ? yobs/(totalData*yexp) : 0. );
	    }
	  TGraph *gr=new TGraph(tmp);
	  grName.ReplaceAll("pu","puwgts");
	  gr->SetName(grName);
	  puWgtGr.push_back( gr );
	  tmp->Delete();
	}
      TH1 *tmp=(TH1 *)puTrue->Clone("tmp");
      TGraph *gr=new TGraph(tmp);
      gr->SetName("puwgts_nom");
      puWgtGr.push_back( gr );
      tmp->Delete();
    }
    if(debug) cout << "loading pileup weight DONE" << endl;

  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);


  //B-TAG CALIBRATION
  TString btagUncUrl(era+"/btagSFactors.csv");
  gSystem->ExpandPathName(btagUncUrl);
  std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
  TString btagEffExpUrl(era+"/expTageff.root");
  gSystem->ExpandPathName(btagEffExpUrl);
  std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
  BTagSFUtil myBTagSFUtil;
  if(!ev.isData)
    {
      BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "central") );
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "down") ); 
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "up") );

      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "central") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "down") ); 
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "up") );
      
      TFile *beffIn=TFile::Open(btagEffExpUrl);
      expBtagEffPy8["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEffPy8["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEffPy8["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();

      TString btagExpPostFix("");
      if(isTTbar)
	{
	  if(filename.Contains("_herwig")) btagExpPostFix="_herwig";
	  if(filename.Contains("_scaleup")) btagExpPostFix="_scaleup";
	  if(filename.Contains("_scaledown")) btagExpPostFix="_scaledown";
	}
      btagEffExpUrl.ReplaceAll(".root",btagExpPostFix+".root");
      beffIn=TFile::Open(btagEffExpUrl);
      expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();

    }

  //jet energy uncertainties
  TString jecUncUrl(era+"/jecUncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  
  //LIST OF SYSTEMATICS
  
  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",4,0,4);
  std::vector<TString> lfsVec = { "_all", "_e", "_ee", "_em", "_mm", "_m" }; 
  std::vector<TString> cutVec = { "", "_lep", "_lepjets", "_jpsi", "_csv", "_meson" };
  std::vector<TString> wgtVec = { "", "_no_weight" };

  for(int i = 0; i < (int)lfsVec.size(); i++) {
  for(int j = 0; j < (int)cutVec.size(); j++) {
  for(int k = 0; k < (int)wgtVec.size(); k++) {
    TString tag(lfsVec[i]);
    TString cut(cutVec[j]);
    TString weight(wgtVec[k]);
    allPlots["lp_pt_iso"+tag+cut+weight] = new TH1F("lp_pt_iso"+tag+cut+weight,";Lepton P_{T} [GeV] after cleaning;Events / 10 GeV", 20, 0,200);
    allPlots["lp_pt_veto"+tag+cut+weight] = new TH1F("lp_pt_veto"+tag+cut+weight,";Lepton P_{T} [GeV] after veto;Events / 10 GeV", 20, 0,200);
    allPlots["lp_pt_low"+tag+cut+weight] = new TH1F("lp_pt_low"+tag+cut+weight,";Leading Lepton P_{T} [GeV];Events / 1 GeV", 20, 20,40);
    allPlots["lp_pt"+tag+cut+weight] = new TH1F("lp_pt"+tag+cut+weight,";Leading Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["l2p_pt"+tag+cut+weight] = new TH1F("l2p_pt"+tag+cut+weight,";Sub-leading Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_pt"+tag+cut+weight] = new TH1F("dilp_pt"+tag+cut+weight,";Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_m"+tag+cut+weight] = new TH1F("dilp_m"+tag+cut+weight,";M_{ll} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["j_pt"+tag+cut+weight] = new TH1F("j_pt"+tag+cut+weight,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["lj_pt"+tag+cut+weight] = new TH1F("lj_pt"+tag+cut+weight,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["bj_pt"+tag+cut+weight] = new TH1F("bj_pt"+tag+cut+weight,";Leading b Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["j_pt_low"+tag+cut+weight] = new TH1F("j_pt_low"+tag+cut+weight,";Leading light Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["lj_pt_low"+tag+cut+weight] = new TH1F("lj_pt_low"+tag+cut+weight,";Leading light Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["bj_pt_low"+tag+cut+weight] = new TH1F("bj_pt_low"+tag+cut+weight,";Leading b Jet P_{T} [GeV];Events / 1 GeV", 20, 30,50);
    allPlots["nlp"+tag+cut+weight]     = new TH1F("nlp"+tag+cut+weight,";N_{l};Events" ,3,0.,3.);
    allPlots["ndilp"+tag+cut+weight]     = new TH1F("ndilp"+tag+cut+weight,";N_{ll};Events" ,3,0.,3.);
    allPlots["nj"+tag+cut+weight]     = new TH1F("nj"+tag+cut+weight,";N_{jets} (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["nlj"+tag+cut+weight]     = new TH1F("nlj"+tag+cut+weight,";N_{jets} (P_{T} > 30 GeV);Events" ,10,0,10.);
    allPlots["nbj"+tag+cut+weight]     = new TH1F("nbj"+tag+cut+weight,";N_{b-jets} (CSV > 0.8);Events" ,4,1.,5.);
    allPlots["npf"+tag+cut+weight]     = new TH1F("npf"+tag+cut+weight,";N_{pf};Events / 10" ,5,0.,5.);
    allPlots["nstart"+tag+cut+weight]     = new TH1F("jetindex"+tag+cut+weight,";N_{jetindex};Events" ,5,0.,5.);
    allPlots["pfid"+tag+cut+weight]     = new TH1F("pfid"+tag+cut+weight,";PFID;Events" ,440,-220.,220.);
    allPlots["massJPsi"+tag+cut+weight]     = new TH1F("massJPsi"+tag+cut+weight,";M_{ll};Events / 18 MeV" ,50,2.5,3.4);
    allPlots["massJPsiK"+tag+cut+weight]     = new TH1F("massJPsiK"+tag+cut+weight,";M_{llk};Events / 15 MeV" ,100,4.5,6);
    allPlots["massD0"+tag+cut+weight]     = new TH1F("massD0"+tag+cut+weight,";M_{D^{0}};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_lep"+tag+cut+weight]     = new TH1F("massD0_lep"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_mu"+tag+cut+weight]     = new TH1F("massD0_mu"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_e"+tag+cut+weight]     = new TH1F("massD0_ele"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massDsmD0loose"+tag+cut+weight]     = new TH1F("massDsmD0loose"+tag+cut+weight,";M_{K#pi#pi} - M_{K#pi};Events / 0.6 MeV" ,25,0.14,0.17);
    allPlots["massDsmD0"+tag+cut+weight]     = new TH1F("massDsmD0"+tag+cut+weight,";M_{K#pi#pi} - M_{K#pi};Events / 0.6 MeV" ,25,0.14,0.17);
    allPlots["massDs"+tag+cut+weight]     = new TH1F("massDs"+tag+cut+weight,";M_{D^{*}};Events / 10 MeV" ,200,0.,2.0);
    allPlots["pi_pt"+tag+cut+weight] = new TH1F("pi_pt"+tag+cut+weight,";#pi^{#pm} P_{T} [GeV];Events / 5 GeV", 10, 0,50);
    allPlots["MET"+tag+cut+weight] = new TH1F("MET"+tag+cut+weight,";MET [GeV];Events / 20 GeV", 10,0,200);
    allPlots["charge"+tag+cut+weight] = new TH1F("charge"+tag+cut+weight,";Charge(l_{1}*l_{2});Events", 5,-2,2);
    allPlots["csv"+tag+cut+weight] = new TH1F("CSV"+tag+cut+weight,";Jet CSV;Events / 0.1", 10,0,1);
    allPlots["dR"+tag+cut+weight] = new TH1F("dR"+tag+cut+weight,";dR;Events / 0.05", 20,0.0,1.);
    allPlots["pflp_pt"+tag+cut+weight] = new TH1F("pflp_pt"+tag+cut+weight,";PF lepton P_{T} [GeV];Events / 0.2 GeV", 15, 0,3);
    allPlots["massZ"+tag+cut+weight]     = new TH1F("massZ_control"+tag+cut+weight,";M_{ll};Events / 1.0 GeV" ,30,81,111);
    allPlots["chargeZ"+tag+cut+weight]     = new TH1F("chargeZ_control"+tag+cut+weight,";M_{ll};Events / 1.0 GeV" ,5,-2,2);
    allPlots["nevt"+tag+cut+weight]     = new TH1F("nevt"+tag+cut+weight,";N_{events};Events" ,1,1.,2.);
    allPlots["weight"+tag+cut+weight]     = new TH1F("weight"+tag+cut+weight,";N_{events};Events/ 1.0" ,20,0.,2.);
    allPlots["norm"+tag+cut+weight]     = new TH1F("norm"+tag+cut+weight,";N_{events};Events / 1.0" ,2,0.,2.);
    allPlots["nvtx"+tag+cut+weight]     = new TH1F("nvtx"+tag+cut+weight,";N_{events};Events / 1.0" ,50,0.,50.);
    allPlots["chi2"+tag+cut+weight] = new TH1F("normchi2"+tag+cut+weight,";#chi^2/n.d.o.f.;Events", 10,0.,10.);
    allPlots["lp_dxy"+tag+cut+weight] = new TH1F("lp_dxy"+tag+cut+weight,";d_{xy} [cm];Events / 0.01 #mum", 20, 0, 0.2);
    allPlots["lp_dz"+tag+cut+weight] = new TH1F("lp_dz"+tag+cut+weight,";d_{z} [cm];Events / 0.01 #mum", 50, 0, 0.5);
    allPlots["pf_dxy"+tag+cut+weight] = new TH1F("pf_dxy"+tag+cut+weight,";d_{xy} [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dz"+tag+cut+weight] = new TH1F("pf_dz"+tag+cut+weight,";d_{z} [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dxyE"+tag+cut+weight] = new TH1F("pf_dxyE"+tag+cut+weight,";#sigma(d_{xy}) [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dzE"+tag+cut+weight] = new TH1F("pf_dzE"+tag+cut+weight,";#sigma(d_{z}) [cm];Events / 0.02 #mum", 20, 0, 0.1);
    allPlots["pf_dxy_sig"+tag+cut+weight] = new TH1F("pf_dxy_significance"+tag+cut+weight,";d_{xy};Events / 1", 30, 0, 30);
    allPlots["pf_dz_sig"+tag+cut+weight] = new TH1F("pf_dz_significance"+tag+cut+weight,";d_{z};Events / 1", 30, 0, 30);

  }
  }
  }
    allPlots["relIso_m"] = new TH1F("relIso_m",";relIso;Events / 0.05", 20,0,1.);
    allPlots["relIso_e"] = new TH1F("relIso_e",";relIso;Events / 0.05", 20,0,1.);
    allPlots["nevt_iso"] = new TH1F("nevt_iso",";After Isolation;Events", 1,1.,2.);
    allPlots["nevt_veto"] = new TH1F("nevt_veto",";After Veto;Events", 1,1.,2.);
    allPlots["norm_iso"] = new TH1F("norm_iso",";After Isolation;Events / 1.0", 20,0,2.);
    allPlots["norm_veto"] = new TH1F("norm_veto",";After Veto;Events / 1.0", 20,0.,2.);
    allPlots["nvtx_iso"]     = new TH1F("nvtx_iso",";N_{events};Events / 1.0" ,50,0.,50.);
    allPlots["nvtx_veto"]     = new TH1F("nvtx_veto",";N_{events};Events / 1.0" ,50,0.,50.);
    allPlots["chi2_iso"] = new TH1F("normchi2_iso",";#chi^2/n.d.o.f.;Events", 10,0.,10.);
    allPlots["chi2_veto"] = new TH1F("normchi2_veto",";#chi^2/n.d.o.f.;Events", 10,0.,10.);
    allPlots["lp_dxy_iso"] = new TH1F("lp_dxy_iso",";d_{xy} [cm];Events / 0.01 #mum", 20, 0, 0.2);
    allPlots["lp_dxy_veto"] = new TH1F("lp_dxy_veto",";d_{xy} [cm];Events / 0.01 #mum", 20, 0, 0.2);
    allPlots["lp_dz_iso"] = new TH1F("lp_dz_iso",";d_{z} [cm];Events / 0.01 #mum", 50, 0, 0.5);
    allPlots["lp_dz_veto"] = new TH1F("lp_dz_veto",";d_{z} [cm];Events / 0.01 #mum", 50, 0, 0.5);


  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  //for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%5000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));
      //Normalize to XSec and lumi
      float norm(1.0);
      if(!ev.isData) {
        norm =  normH ? normH->GetBinContent(1) : 1.0;
	//update nominal event weight
	if(ev.ttbar_nw>0) norm*=ev.ttbar_w[0];
      }
      allPlots["nevt_all"]->Fill(1,norm);
      allPlots["norm_all"]->Fill(norm,norm);
      allPlots["nvtx_all"]->Fill(ev.nvtx,norm);

      //Basic lepton kinematics
      std::vector<int> tightLeptons,vetoLeptons;
      for(int il=0; il<ev.nl; il++)
	{
          //cout << "in lepton selection" << endl;
          float relIso = ev.l_relIso[il];
	  bool passTightKin((ev.l_pt[il]>20 && fabs(ev.l_eta[il])<2.4 && abs(ev.l_id[il])==13 && relIso<0.25) //muons
               || (ev.l_pt[il]>30 && ((fabs(ev.l_eta[il])<1.479 && relIso<0.0893)  //electrons small eta
                  || (fabs(ev.l_eta[il])>1.479 && fabs(ev.l_eta[il])<2.5 && relIso<0.121)) //electrons medium eta
                  && abs(ev.l_id[il])==11 // TOP mu cut for dilep
                  && (fabs(ev.l_eta[il]<1.4442) && fabs(ev.l_eta[il])>1.5660))); //remove ECAL 1.4442<|eta|<1.5660
          //passTightKin &= (abs(ev.l_id[il])==13 && ev.l_dxy[il]<0.2 && ev.l_dz[il]<0.5); //muon impact parameters

	  bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1); //tight muon ID or tight ID except Iso electron
	  
          //Check veto here, but ONLY for lep+jets (later on in code)
	  bool passVetoIso(  ev.l_id[il]==13 ? relIso<0.24 : true); //FIXME from 7_6_x
          bool passVetoKin(  ev.l_pt[il]>15 && fabs(ev.l_eta[il])<2.4); // TOP veto 10->15

	  //bool passSIP3d(ev.l_ip3dsig[il]<4);
	  //if(channelSelection==21) passSIP3d=true;

	  if(passTightKin && passTightId)// && passSIP3d)
	    {
	      tightLeptons.push_back(il);
	      //else if(passNonIso) tightLeptonsNonIso.push_back(il); //FIXME from 7_6_x
	    }
	  else if(passVetoKin && passVetoIso) vetoLeptons.push_back(il); //FIXME from 7_6_x
	}
      if(debug) cout << "lepton selection DONE" << endl;

      //check if triggers have fired
      //bool hasEETrigger(((ev.elTrigger>>3)&0x1)!=0 || ((ev.elTrigger>>2)&0x1)!=0);//HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v || HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
      //Dielectron
      bool hasEETrigger(((ev.elTrigger>>2)&0x1)!=0);                              //HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
      //Dimuon
      bool hasMMTrigger(((ev.muTrigger>>4)&0x1)!=0);                              //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
      hasMMTrigger |= ((ev.muTrigger>>5)&0x1)!=0;                                 //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v
      //Electron Muon (ME as well)
      bool hasEMTrigger(((ev.elTrigger>>4)&0x1)!=0);                              //HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v 
                                                                                  //(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v  recommended by TOP PAG)
                                                                                  //(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v recommended by TOP PAG)
      //Muon Electron part
      hasEMTrigger |= ((ev.elTrigger>>5)&0x1)!=0;                                 //HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v
                                                                                  //(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
                                                                                  //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v 
                                                                                  //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v recommended by TOP PAG)
      //Single muon
      bool hasMuTrigger(((ev.muTrigger>>2)&0x1)!=0);                              //HLT_IsoMu22_v (24 recommended by TOP PAG)
      hasMuTrigger |= (((ev.muTrigger>>3)&0x1)!=0);                               //HLT_IsoTkMu22_v
      //Single electorn
      bool hasEleTrigger((ev.elTrigger & 0x1)!=0);                                //HLT_Ele27_WPTight_Gsf_v

      //No trigger requirement for MC
      if(!ev.isData)
	{	 
	  hasMuTrigger=true;
	  hasEleTrigger=true;
	  hasEETrigger=true;
	  hasMMTrigger=true;
	  hasEMTrigger=true;
	}
      else
	{
	  if(requireMutriggerOnly && !hasMuTrigger) continue;
	  if(requireEletriggerOnly && !hasEleTrigger) continue;
	  if(requireEETriggers && !hasEETrigger) continue;
	  if(requireMMTriggers && !hasMMTrigger) continue;
	  if(requireEMTriggers && !hasEMTrigger) continue;
	}

      //decide the channel
      if(debug) cout << "decide channel" << endl;
      TString chTag("");
      std::vector<int> selLeptons;
      bool passTightKin(false),passIso(false);
      if(tightLeptons.size()==1 )
	{
          //** Tighter cuts for lepton + jets **
          if(ev.l_id[tightLeptons[0]]==13) { // muon + jets
	    passTightKin = (ev.l_pt[tightLeptons[0]] > 26 && fabs(ev.l_eta[tightLeptons[0]])<2.1); // TOP mu cut for lep+jets
            passIso = (ev.l_relIso[tightLeptons[0]] < 0.15); //TOP mu cut for lep+jets
	    allPlots["lp_pt_iso_m"]->Fill(ev.l_pt[tightLeptons[0]],norm);
            allPlots["relIso_m"]->Fill(ev.l_relIso[tightLeptons[0]],norm);
            //allPlots["chi2_iso"]->Fill(ev.l_chi2norm[tightLeptons[0]],norm);
            //allPlots["dxy_iso"]->Fill(ev.l_dxy[tightLeptons[0]],norm);
            //allPlots["dz_iso"]->Fill(ev.l_dz[tightLeptons[0]],norm);
            if(vetoLeptons.size()==0) {
              allPlots["lp_pt_veto_m"]->Fill(ev.l_pt[tightLeptons[0]],norm);
              //allPlots["chi2_veto"]->Fill(ev.l_chi2norm[tightLeptons[0]],norm);
              //allPlots["dxy_veto"]->Fill(ev.l_dxy[tightLeptons[0]],norm);
              //allPlots["dz_veto"]->Fill(ev.l_dz[tightLeptons[0]],norm);
            }
          }
          else if(ev.l_id[tightLeptons[0]]==11) { // electron + jets
            passTightKin = (ev.l_pt[tightLeptons[0]] > 30); //from TOP-15-005
            passIso = (ev.l_relIso[tightLeptons[0]] < 0.15); //TOP mu cut for lep+jets
            allPlots["lp_pt_iso_e"]->Fill(ev.l_pt[tightLeptons[0]],norm);
            allPlots["relIso_e"]->Fill(ev.l_relIso[tightLeptons[0]],norm);
            if(vetoLeptons.size()==0)
              allPlots["lp_pt_veto_e"]->Fill(ev.l_pt[tightLeptons[0]],norm);
          }

          if(passTightKin && passIso) {
	    selLeptons.push_back( tightLeptons[0] );
            if(debug) cout << "found 1 tight lepton" << endl;
          }
          //no extra isolated leptons
          //if(vetoLeptons.size()>0) continue; //veto only on lep+jets
	}
      if(tightLeptons.size()>=2)
	{
	  selLeptons.push_back(tightLeptons[0]);
	  selLeptons.push_back(tightLeptons[1]);
          if(debug) cout << "found 2 tight leptons" << endl;
	}
      if(debug) if(selLeptons.size()==0) cout << "NO LEPTONS!!" << endl;
      if(selLeptons.size()==0) continue;
      if(selLeptons.size()==1)
	{
	  if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) {
            if(requireEletriggerOnly || !ev.isData) chTag="e"; //Ensure SingleElectorn if data
          }
	  if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger) {
            if(requireMutriggerOnly || !ev.isData) chTag="m"; //Ensure SingleMuon if data
          }
	}
      if(selLeptons.size()==2)
	{
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*11 && hasEETrigger) {
            if(requireEETriggers || !ev.isData) chTag="ee"; //Ensure DoubleEG if data
          }
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==13*13 && hasMMTrigger) {
            if(requireMMTriggers || !ev.isData) chTag="mm"; //Ensure DoubleMuon if data
          }
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*13 && hasEMTrigger) {
            if(requireEMTriggers || !ev.isData) chTag="em"; //Ensure MuonEG if data
          }
	  if(hasMuTrigger && requireEletriggerOnly) chTag="";
	}
      if(chTag=="") continue;
      chTag = "_"+chTag;
      if(debug) cout << "decide channel DONE" << endl;

      //one good lepton either isolated or in the non-isolated sideband or a Z candidate
      Int_t lepIdx=-1;
      Bool_t isZ(false);//,isZPassingSIP3d(false);
      TLorentzVector l1p4,l2p4,dilp4;
      if(selLeptons.size()==1)                                       lepIdx=selLeptons[0];
      //else if (selLeptons.size()==0 && selLeptonsNonIso.size()==1) lepIdx=selLeptonsNonIso[0];
      else if(selLeptons.size()==2)
	{	  
          if(debug) cout << "di-lepton" << endl;
	  l1p4.SetPtEtaPhiM(ev.l_pt[selLeptons[0]],ev.l_eta[selLeptons[0]],ev.l_phi[selLeptons[0]],ev.l_mass[selLeptons[0]]);
	  l2p4.SetPtEtaPhiM(ev.l_pt[selLeptons[1]],ev.l_eta[selLeptons[1]],ev.l_phi[selLeptons[1]],ev.l_mass[selLeptons[1]]);
	  dilp4=l1p4+l2p4;
	  if(ev.l_id[selLeptons[0]]==ev.l_id[selLeptons[1]]          && 
	     ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]]<0 && 
	     fabs(dilp4.M()-91)<15)
	    { 
	      isZ=true; 
	      //isZPassingSIP3d=(ev.l_ip3dsig[0]<4 && ev.l_ip3dsig[1]<4);
	    }
	  lepIdx=selLeptons[0];
          if(debug) cout << "di-lepton DONE" << endl;
	}

      if(lepIdx<0) continue;
      allPlots["nevt_iso"]->Fill(1,norm);
      allPlots["norm_iso"]->Fill(norm,norm);
      allPlots["nvtx_iso"]->Fill(ev.nvtx,norm);
      
      //USE VETO HERE
      if(vetoLeptons.size()>0) continue; //veto only on lep+jets
      allPlots["nevt_veto"]->Fill(1,norm);
      allPlots["norm_veto"]->Fill(norm,norm);
      allPlots["nvtx_veto"]->Fill(ev.nvtx,norm);
      
      //save lepton kinematics
      std::vector<TLorentzVector> leptons;
      for(size_t il=0; il<selLeptons.size(); il++)
	{
	  int lepIdx=selLeptons[il];
	  TLorentzVector lp4;
	  lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);
	  leptons.push_back(lp4);
	}

      //select jets
      Float_t htsum(0);
      TLorentzVector jetDiff(0,0,0,0);
      int nbjets(0),ncjets(0),nljets(0);//,leadingJetIdx(-wgt);
      std::vector<int> resolvedJetIdx;
      std::vector<TLorentzVector> resolvedJetP4;
      std::vector<Jet> bJetsVec, lightJetsVec, allJetsVec;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
	  //jp4=updateJES(jp4,ev.j_rawsf[k],ev.j_area[k],ev.rho,ev.nvtx,jetCorr);

	  //cross clean with respect to leptons deltaR<0.4
          bool overlapsWithLepton(false);
          for(size_t il=0; il<leptons.size(); il++) {
            if(jp4.DeltaR(leptons[il])>0.4) continue;  //Jet is fine
	    overlapsWithLepton=true;                   //Jet ovelaps with an "isolated" lepton, event is bad
          }
          if(overlapsWithLepton) continue;
          if(debug) cout << "Overlap with lepton DONE" << endl;

	  //smear jet energy resolution for MC
	  //jetDiff -= jp4;
	  float genJet_pt(0);
	  if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }
	  //jetDiff += jp4;
	  resolvedJetIdx.push_back(k);
	  resolvedJetP4.push_back(jp4);

	  // re-inforce kinematics cuts
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  //if(leadingJetIdx<0) leadingJetIdx=k;
	  htsum += jp4.Pt();

	  //b-tag
	  if(debug) cout << "Starting b-tagging" << endl;
	  float csv = ev.j_csv[k];	  
	  bool isBTagged(csv>0.800);
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0);
	      if(abs(ev.j_hadflav[k])==4) 
		{ 
		  ncjets++;
		  expEff    = expBtagEff["c"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  nbjets++;
		  expEff    = expBtagEff["b"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
		}
	      else
		{
		  nljets++;
		  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
		}
	      
	      //updated b-tagging decision with the data/MC scale factor
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
	    }
	  if(debug) cout << "b-tagging DONE" << endl;

	  //save jet
          Jet tmpj(jp4, csv, k);
	  for(int ipf = 0; ipf < ev.npf; ipf++) {
	    if(ev.pf_j[ipf] != k) continue; //skip if PF track doesn't belong to current jet
	    if(ev.pf_c[ipf]==0) continue;   //skip if PF track is neutral
	    TLorentzVector tkP4(0,0,0,0);
	    tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);
            pfTrack pftk(tkP4, ev.pf_dxy[ipf], ev.pf_dxyE[ipf], ev.pf_dz[ipf], ev.pf_dzE[ipf], ev.pf_id[ipf]);
	    tmpj.addTrack(pftk,ev.pf_id[ipf]);
	  }
          tmpj.sortTracksByPt();

          if(isBTagged) bJetsVec.push_back(tmpj);
          else lightJetsVec.push_back(tmpj);
          allJetsVec.push_back(tmpj);
	}

      
      //event weight
      float wgt(1.0);

      std::vector<float> puWgts(3,1.0),topPtWgts(2,1.0);
      EffCorrection_t lepSelCorrWgt(1.0,0.0), triggerCorrWgt(1.0,0.0);
      if(debug) cout << "Lepton scale factors" << endl;
      if(!ev.isData)
	{
	  //account for pu weights and effect on normalization
	  allPlots["puwgtctr"]->Fill(0.,1.0);
	  if(debug) cout << "getting puWgts" << endl;
	    for(size_t iwgt=0; iwgt<3; iwgt++)
	      {
	        puWgts[iwgt]=puWgtGr[iwgt]->Eval(ev.putrue);  
	        allPlots["puwgtctr"]->Fill(iwgt+1,puWgts[iwgt]);
	      }
	  if(debug) cout << "getting puWgts DONE!" << endl;
	  //trigger/id+iso efficiency corrections
          if(debug) cout << "calling trigger function" << endl;
          std::vector<int> pdgIds; //vector of IDs for trigger correction function
          for(size_t ilp = 0; ilp < selLeptons.size(); ilp++)
            pdgIds.push_back(ev.l_id[selLeptons[ilp]]);
	  //triggerCorrWgt=lepEffH.getTriggerCorrection(selLeptons,leptons);
	  triggerCorrWgt=lepEffH.getTriggerCorrection(pdgIds,leptons); //FIXME
          if(debug) cout << "calling trigger function DONE!" << endl;
          // ** selLeptons contains only ev_l position, leptons contains p4 **
	  for(size_t il=0; il<selLeptons.size(); il++) {
	    EffCorrection_t selSF=lepEffH.getOfflineCorrection(pdgIds[il],leptons[il].Pt(),leptons[il].Eta());
	    lepSelCorrWgt.second = sqrt( pow(lepSelCorrWgt.first*selSF.second,2)+pow(lepSelCorrWgt.second*selSF.first,2));
            if(debug) cout << "lepSelCorrWgt=" << lepSelCorrWgt.first << endl;
            if(debug) cout << "selSF=" << selSF.first << endl;
	    lepSelCorrWgt.first *= selSF.first;
	  }
	  wgt=triggerCorrWgt.first*lepSelCorrWgt.first*puWgts[0]*norm;
          if(debug) cout << "weight=" << wgt << endl;
          if(debug) cout << "Trigger=" << triggerCorrWgt.first << endl << "Lepton=" << lepSelCorrWgt.first << endl << "PU=" << puWgts[0] << endl << "norm=" << norm  << endl;
          //wgt=1.0;
	}
      if(debug) cout << "Lepton scale factors DONE!" << endl;

      //sort by Pt
      sort(lightJetsVec.begin(),    lightJetsVec.end(),   sortJetsByPt);
      sort(bJetsVec.begin(),    bJetsVec.end(),   sortJetsByPt);
      sort(allJetsVec.begin(),  allJetsVec.end(), sortJetsByPt);

      for(size_t ij = 0; ij < allJetsVec.size(); ij++) {
        float csv = allJetsVec.at(ij).getCSV();
        allPlots["csv"+chTag]->Fill(csv,wgt);
        allPlots["csv_all"]->Fill(csv,wgt);
      }

      for(size_t il=0; il<leptons.size(); il++) {
        for(size_t ij=0; ij<allJetsVec.size(); ij++) {
	  TLorentzVector jp4=allJetsVec[ij].getVec();
          allPlots["dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
          allPlots["dR"+chTag+"_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
        }
      /*
        for(size_t ij=0; ij<bJetsVec.size(); ij++) {
	  TLorentzVector jp4=bJetsVec[ij].getVec();
	  jp4.SetPtEtaPhiM(ev.j_pt[ij],ev.j_eta[ij],ev.j_phi[ij],ev.j_mass[ij]);
          allPlots["dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
          allPlots["dR"+chTag+"_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
        }
        for(size_t ij=0; ij<lightJetsVec.size(); ij++) {
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[ij],ev.j_eta[ij],ev.j_phi[ij],ev.j_mass[ij]);
          allPlots["dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
          allPlots["dR_all"]->Fill(jp4.DeltaR(leptons[il]),wgt);
          allPlots["dR"+chTag+"_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
        }
      */
      }


      //MET and transverse mass
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
      met+=jetDiff;
      met.SetPz(0.); met.SetE(met.Pt());
      //float mt( computeMT(isZ ? dilp4: lp4,met) );

      //simple fill
      bool singleLep(false);
      bool doubleLep(false);
      bool minJets(false);
      if(debug) cout << "sorting jets" << endl;
      if(debug) cout << "starting simple plots" << endl;
      allPlots["nevt"+chTag]->Fill(1,norm);
      allPlots["weight"+chTag]->Fill(wgt,norm);
      allPlots["norm"+chTag]->Fill(norm,norm);
      allPlots["nvtx"+chTag]->Fill(ev.nvtx,norm);
      allPlots["nj"+chTag]->Fill(allJetsVec.size(),wgt);
      allPlots["nlj"+chTag]->Fill(lightJetsVec.size(),wgt);
      allPlots["nbj"+chTag]->Fill(bJetsVec.size(),wgt);
      allPlots["nlj"+chTag+"_no_weight"]->Fill(lightJetsVec.size(),norm);
      allPlots["nbj"+chTag+"_no_weight"]->Fill(bJetsVec.size(),norm);
      allPlots["nlj_all"]->Fill(lightJetsVec.size(),wgt);
      allPlots["nbj_all"]->Fill(bJetsVec.size(),wgt);
      allPlots["nlp"+chTag]->Fill(selLeptons.size(),wgt);
      allPlots["nlp"+chTag+"_no_weight"]->Fill(selLeptons.size(),norm);

      if(leptons.size() > 0) {
        allPlots["lp_pt_low"+chTag]->Fill(leptons[0].Pt(),wgt);
        allPlots["lp_pt"+chTag]->Fill(leptons[0].Pt(),wgt);
        allPlots["lp_pt"+chTag+"_no_weight"]->Fill(leptons[0].Pt(),norm);
        allPlots["lp_pt_all"]->Fill(leptons[0].Pt(),wgt);
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
      if(bJetsVec.size() > 0) {
        allPlots["bj_pt_low"+chTag]->Fill(bJetsVec[0].getVec().Pt(),wgt);
        allPlots["bj_pt"+chTag]->Fill(bJetsVec[0].getVec().Pt(),wgt);
        allPlots["bj_pt"+chTag+"_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
        allPlots["bj_pt_all"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
      }

      //if(selLeptons.size() == 1 && bJetsVec.size() >= 1 && lightJetsVec.size() >= 2) {
      if(debug) cout << "starting lep/jets plots" << endl;
      //Require exactly 1 lepton
      if(selLeptons.size() == 1) {
        if(debug) cout << "single lepton" << endl;
        singleLep = true;
        allPlots["lp_pt"+chTag+"_lep"]->Fill(leptons[0].Pt(),wgt);
        allPlots["lp_pt"+chTag+"_lep_no_weight"]->Fill(leptons[0].Pt(),norm);

        allPlots["lp_pt_all_lep"]->Fill(leptons[0].Pt(),wgt);
        allPlots["MET"+chTag+"_lep"]->Fill(ev.met_pt[0],wgt);
        allPlots["MET"+chTag+"_lep_no_weight"]->Fill(ev.met_pt[0],norm);
        //Require at least 1 b-tagged and at least 2 light jets
        if(bJetsVec.size() >= 1 && lightJetsVec.size() >= 2) {
          if(debug) cout << "jet reqirements" << endl;
          minJets = true;
          std::sort(leptons.begin(), leptons.end(), VecSort);
          allPlots["lp_pt"+chTag+"_lepjets"]->Fill(leptons[0].Pt(),wgt);
          //allPlots["lp_pt"+chTag+"_lepjets_no_weight"]->Fill(leptons[0].Pt(),norm);

          allPlots["lp_pt_all_lepjets"]->Fill(leptons[0].Pt(),wgt);
          allPlots["MET"+chTag+"_lepjets"]->Fill(ev.met_pt[0],wgt);
          //allPlots["MET"+chTag+"_lepjets__no_weight"]->Fill(ev.met_pt[0],norm);
        }
      }
      //Require exactly 2 leptons
      else if(selLeptons.size() == 2) {
        if(debug) cout << "dilepton" << endl;
        doubleLep = true;
        //Z control plot
        if(isZ) {
          allPlots["massZ"+chTag+"_lep"]->Fill(dilp4.M(),wgt);
          allPlots["massZ"+chTag+"_lep_no_weight"]->Fill(dilp4.M(),norm);
        }
        if(abs(dilp4.M()-91)<15)
          allPlots["chargeZ"+chTag]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt); 
        //Exclude Z and low mass and require same falvor dilepton MET > 40 GeV
        if(!isZ && dilp4.M() > 20 &&
          ((ev.l_id[selLeptons[0]]!=ev.l_id[selLeptons[1]]) || (ev.l_id[selLeptons[0]]==ev.l_id[selLeptons[1]] && met.Pt() > 40))) {
          allPlots["ndilp"+chTag+"_lep"]->Fill(selLeptons.size(),wgt);
          allPlots["dilp_pt"+chTag+"_lep"]->Fill(dilp4.Pt(),wgt);
          allPlots["dilp_m"+chTag+"_lep"]->Fill(dilp4.M(),wgt);
          allPlots["ndilp"+chTag+"_lep"+"_no_weight"]->Fill(selLeptons.size(),norm);
          allPlots["dilp_pt"+chTag+"_lep"+"_no_weight"]->Fill(dilp4.Pt(),norm);
          allPlots["dilp_m"+chTag+"_lep"+"_no_weight"]->Fill(dilp4.M(),norm);
          std::sort(leptons.begin(), leptons.end(), VecSort);
          allPlots["lp_pt"+chTag+"_lep"]->Fill(leptons[0].Pt(),wgt);
          allPlots["l2p_pt"+chTag+"_lep"]->Fill(leptons[1].Pt(),wgt);
          allPlots["lp_pt"+chTag+"_lep"+"_no_weight"]->Fill(leptons[0].Pt(),norm);
          allPlots["l2p_pt"+chTag+"_lep"+"_no_weight"]->Fill(leptons[1].Pt(),norm);
          allPlots["MET"+chTag+"_lep"]->Fill(ev.met_pt[0],wgt);
          allPlots["charge"+chTag+"_lep"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
          allPlots["MET"+chTag+"_lep"+"_no_weight"]->Fill(ev.met_pt[0],norm);
          allPlots["charge"+chTag+"_lep"+"_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],norm);
        }
        //Require at least 1 b-tagged and at least 2 light jets
        if(bJetsVec.size() >= 1 && lightJetsVec.size() > 1) {
          if(debug) cout << "jet reqirements" << endl;
          //Z control plot
          if(isZ) {
            allPlots["massZ"+chTag+"_lepjets"]->Fill(dilp4.M(),wgt);
            allPlots["massZ"+chTag+"_lepjets_no_weight"]->Fill(dilp4.M(),norm);
          }
          //Exclude Z mass
          if(isZ) continue;
          //Exclude low mass (M < 20 GeV)
          if(dilp4.M() < 20) continue;
          //Require same falvor dilepton MET > 40 GeV
          if(ev.l_id[selLeptons[0]]==ev.l_id[selLeptons[1]] && met.Pt() < 40) continue; //FIXME
          minJets = true;
          allPlots["ndilp"+chTag+"_lepjets"]->Fill(selLeptons.size(),wgt);
          allPlots["dilp_pt"+chTag+"_lepjets"]->Fill(dilp4.Pt(),wgt);
          allPlots["dilp_m"+chTag+"_lepjets"]->Fill(dilp4.M(),wgt);
          allPlots["ndilp"+chTag+"_lepjets"+"_no_weight"]->Fill(selLeptons.size(),norm);
          allPlots["dilp_pt"+chTag+"_lepjets"+"_no_weight"]->Fill(dilp4.Pt(),norm);
          allPlots["dilp_m"+chTag+"_lepjets"+"_no_weight"]->Fill(dilp4.M(),norm);
          std::sort(leptons.begin(), leptons.end(), VecSort);
          allPlots["lp_pt"+chTag+"_lepjets"]->Fill(leptons[0].Pt(),wgt);
          allPlots["l2p_pt"+chTag+"_lepjets"]->Fill(leptons[1].Pt(),wgt);
          allPlots["lp_pt"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons[0].Pt(),norm);
          allPlots["l2p_pt"+chTag+"_lepjets"+"_no_weight"]->Fill(leptons[1].Pt(),norm);
          allPlots["MET"+chTag+"_lepjets"]->Fill(ev.met_pt[0],wgt);
          allPlots["charge"+chTag+"_lepjets"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
          allPlots["MET"+chTag+"_lepjets"+"_no_weight"]->Fill(ev.met_pt[0],norm);
          allPlots["charge"+chTag+"_lepjets"+"_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],norm);
        }
      }
      if(debug) cout << "simple plots DONE" << endl;


      //Require lep+jets or dilepton
      if(!singleLep && !doubleLep) continue;
      if(debug) cout << "passed lep requirements" << endl;

      allPlots["npf"+chTag+"_lep"]->Fill(ev.npf,wgt);
      //allPlots["npf"+chTag+"_lep"+"_no_weight"]->Fill(ev.npf,norm);
      allPlots["nevt"+chTag+"_lep"]->Fill(1,norm);
      allPlots["weight"+chTag+"_lep"]->Fill(wgt,norm);
      allPlots["norm"+chTag+"_lep"]->Fill(norm,norm);
      allPlots["nvtx"+chTag+"_lep"]->Fill(ev.nvtx,norm);
      allPlots["nevt_all_lep"]->Fill(1,norm);

      allPlots["nj"+chTag+"_lep"]->Fill(allJetsVec.size(),wgt);
      allPlots["nlj"+chTag+"_lep"]->Fill(lightJetsVec.size(),wgt);
      allPlots["nbj"+chTag+"_lep"]->Fill(bJetsVec.size(),wgt);
      //allPlots["nj"+chTag+"_lep"+"_no_weight"]->Fill(lightJetsVec.size(),norm);
      //allPlots["nbj"+chTag+"_lep"+"_no_weight"]->Fill(bJetsVec.size(),norm);
      allPlots["nlp"+chTag+"_lep"]->Fill(selLeptons.size(),wgt);
      //allPlots["nlp"+chTag+"_lep"+"_no_weight"]->Fill(selLeptons.size(),norm);

      if(lightJetsVec.size() > 0 and bJetsVec.size() > 0) {
      allPlots["j_pt"+chTag+"_lep"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
      allPlots["lj_pt"+chTag+"_lep"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
      allPlots["bj_pt"+chTag+"_lep"]->Fill(bJetsVec[0].getVec().Pt(),wgt);

      for(size_t ij = 0; ij < allJetsVec.size(); ij++) {
        float csv = allJetsVec.at(ij).getCSV();
        allPlots["csv"+chTag+"_lep"]->Fill(csv,wgt);
        allPlots["csv_all_lep"]->Fill(csv,wgt);
      }
      }


      //Require b-tagged and light jets
      if(!minJets) continue;
      if(debug) cout << "passed jet requirements" << endl;

      allPlots["npf"+chTag+"_lepjets"]->Fill(ev.npf,wgt);
      allPlots["npf"+chTag+"_lepjets"+"_no_weight"]->Fill(ev.npf,norm);
      allPlots["nevt"+chTag+"_lepjets"]->Fill(1,norm);
      allPlots["weight"+chTag+"_lepjets"]->Fill(wgt,norm);
      allPlots["norm"+chTag+"_lepjets"]->Fill(norm,norm);
      allPlots["nvtx"+chTag+"_lepjets"]->Fill(ev.nvtx,norm);
      allPlots["nevt_all_lepjets"]->Fill(1,norm);

      allPlots["nj"+chTag+"_lepjets"]->Fill(allJetsVec.size(),wgt);
      allPlots["nlj"+chTag+"_lepjets"]->Fill(lightJetsVec.size(),wgt);
      allPlots["nbj"+chTag+"_lepjets"]->Fill(bJetsVec.size(),wgt);
      allPlots["nlp"+chTag+"_lepjets"]->Fill(selLeptons.size(),wgt);

      allPlots["j_pt"+chTag+"_lepjets"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
      allPlots["lj_pt"+chTag+"_lepjets"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
      allPlots["bj_pt"+chTag+"_lepjets"]->Fill(bJetsVec[0].getVec().Pt(),wgt);

      for(size_t ij = 0; ij < allJetsVec.size(); ij++) {
        float csv = allJetsVec.at(ij).getCSV();
        allPlots["csv"+chTag+"_lepjets"]->Fill(csv,wgt);
        allPlots["csv_all_lepjets"]->Fill(csv,wgt);
      }

      //charmed resonance analysis : use only jets with CSV>CSVL, up to two per event
      for(size_t ij = 0; ij < bJetsVec.size(); ij++) {

        if(ij > 1) continue;
        if(ij == 0) {
          //allPlots["nbj"+chTag+"_csv"]->Fill(1,wgt);
          allPlots["nevt"+chTag+"_csv"]->Fill(1,norm);
          allPlots["weight"+chTag+"_csv"]->Fill(wgt,norm);
          allPlots["norm"+chTag+"_csv"]->Fill(norm,norm);
          allPlots["nvtx"+chTag+"_csv"]->Fill(ev.nvtx,norm);
          allPlots["j_pt"+chTag+"_csv"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
          allPlots["lj_pt"+chTag+"_csv"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
          allPlots["bj_pt"+chTag+"_csv"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
        }

        allPlots["csv"+chTag+"_csv"]->Fill(bJetsVec[ij].getCSV(),wgt);
        std::vector<IdTrack> &tracks = bJetsVec[ij].getTracks();

        //J/Psi
        if(debug) cout << "starting J/Psi" << endl;
        const float gMassMu(0.1057),gMassK(0.4937),gMassPi(0.1396);
        std::vector<pfTrack> pfmuCands,kaonCands;
        for(size_t itk = 0; itk < tracks.size(); itk++) {
          if(abs(tracks[itk].second) == 13) {
            TLorentzVector muP4;
            muP4.SetPtEtaPhiM(tracks[itk].first.Pt(), tracks[itk].first.Eta(), tracks[itk].first.Phi(), gMassMu);
            pfTrack pfmu(muP4, tracks[itk].first.getDxy(), tracks[itk].first.getDxyE(), tracks[itk].first.getDz(), tracks[itk].first.getDzE(), tracks[itk].second);
            pfmuCands.push_back(pfmu);
          }
          if(abs(tracks[itk].second) == 211) {
            TLorentzVector kP4;
            pfTrack pfk(kP4, tracks[itk].first.getDxy(), tracks[itk].first.getDxyE(), tracks[itk].first.getDz(), tracks[itk].first.getDzE(), tracks[itk].second);
            kaonCands.push_back(pfk);
          }
        }
    
        if(pfmuCands.size()>1) {
          if(pfmuCands[0].getPfid() != -pfmuCands[1].getPfid()) continue;
          float mass12((pfmuCands[0].getVec() + pfmuCands[1].getVec()).M());
          float mass123( kaonCands.size()>0 ? (pfmuCands[0].getVec()+pfmuCands[1].getVec()+kaonCands[0].getVec()).M() : -1);

          allPlots["massJPsi"+chTag]->Fill(mass12,wgt);
	  allPlots["massJPsi_all"]->Fill(mass12,wgt);
          allPlots["nbj"+chTag+"_jpsi"]->Fill(1,wgt);
          allPlots["bj_pt"+chTag+"_jpsi"]->Fill(bJetsVec[ij].getVec().Pt(),wgt);
          allPlots["lp_pt"+chTag+"_jpsi"]->Fill(leptons[0].Pt(),wgt);
          allPlots["lp_pt_low"+chTag+"_jpsi"]->Fill(leptons[0].Pt(),wgt);
          allPlots["csv"+chTag+"_jpsi"]->Fill(bJetsVec[ij].getCSV(),wgt);
          allPlots["nevt"+chTag+"_jpsi"]->Fill(1,norm);
          allPlots["weight"+chTag+"_jpsi"]->Fill(wgt,norm);
          allPlots["norm"+chTag+"_jpsi"]->Fill(norm,norm);
          allPlots["nvtx"+chTag+"_jpsi"]->Fill(ev.nvtx,norm);

          if(mass12<3.0 || mass12>3.2) continue;

          for(int itk = 0; itk < 2; itk++) {

            for(int i = 0; i < 2; i++) {
              allPlots["pf_dxy"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDxy()),wgt);
              allPlots["pf_dz"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDz()),wgt);
              allPlots["pf_dxyE"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDxyE()),wgt);
              allPlots["pf_dzE"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDzE()),wgt);
              allPlots["pf_dz_sig"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
              allPlots["pf_dxy_sig"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDxy())/abs(pfmuCands[i].getDxyE()),wgt);
              allPlots["pf_dz_sig"+chTag+"_jpsi"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
              allPlots["pf_dxy_all"]->Fill(abs(pfmuCands[i].getDxy()),wgt);
              allPlots["pf_dz_all"]->Fill(abs(pfmuCands[i].getDz()),wgt);
              allPlots["pf_dxyE_all"]->Fill(abs(pfmuCands[i].getDxyE()),wgt);
              allPlots["pf_dzE_all"]->Fill(abs(pfmuCands[i].getDzE()),wgt);
              allPlots["pf_dz_sig_all"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
              allPlots["pf_dxy_sig_all"]->Fill(abs(pfmuCands[i].getDxy())/abs(pfmuCands[i].getDxyE()),wgt);
              allPlots["pf_dz_sig_all"]->Fill(abs(pfmuCands[i].getDz())/abs(pfmuCands[i].getDzE()),wgt);
            }
          }

          if(mass123 > 0) {
            allPlots["massJPsiK"+chTag]->Fill(mass123,wgt);
            allPlots["massJPsiK_all"]->Fill(mass123,wgt);
          }
        }
        if(debug) cout << "J/Psi DONE" << endl;
        //continue; //FIXME

        //D0 and D* 
        if(debug) cout << "Starting D0 and D*" << endl;
        int jetindex = allJetsVec[ij].getJetIndex();
        if(tracks.size() < 3) continue;
        size_t tmax = 4;
        tmax = tracks.size() >= tmax ? tmax : tracks.size();
        for(size_t i = 0; i < tmax; i++)
          for(size_t j = 0; j < tmax; j++) {
            if(i == j) continue;

            //opposite sign
            if(tracks[i].second*tracks[j].second != -211*211) continue;

            TLorentzVector p_track1, p_track2;
            p_track1.SetPtEtaPhiM(tracks[i].first.Pt(), tracks[i].first.Eta(), tracks[i].first.Phi(), gMassPi);
            p_track2.SetPtEtaPhiM(tracks[j].first.Pt(), tracks[j].first.Eta(), tracks[j].first.Phi(), gMassK);
            if(debug) cout << i << ": " << tracks[i].first.Pt() << " " << tracks[i].first.Eta() << " " << tracks[i].first.Phi() << " " << gMassPi << endl;
            if(debug) cout << j << ": " << tracks[j].first.Pt() << " " << tracks[j].first.Eta() << " " << tracks[j].first.Phi() << " " << gMassK << endl << endl;
            float mass12 = (p_track1+p_track2).M();
            if(debug) cout << mass12 << endl;
            allPlots["dR"+chTag+"_meson"]->Fill(p_track1.DeltaR(p_track2), wgt);
            //allPlots["dR"+chTag+"_meson_no_weight"]->Fill(p_track1.DeltaR(p_track2),norm);

            if (mass12>1.65 && mass12<2.0) {
              allPlots["massD0"+chTag]->Fill(mass12,wgt);
              allPlots["massD0_all"]->Fill(mass12,wgt);
              //allPlots["massD0"+chTag+"_no_weight"]->Fill(mass12,norm);
              allPlots["nbj"+chTag+"_meson"]->Fill(1,wgt);
            }

            //looking for lepton
            if(debug) cout << "third lepton" << endl;
            //for(int tk3 = 0; tk3 < ev.npf; tk3++)
            for(size_t k = 0; k < tracks.size(); k++) {
              if(k == i) continue;
              if(k == j) continue;
              if(debug) cout << "third lepton possible" << endl;
            
              if(abs(tracks[k].second) != 13 && abs(tracks[k].second) != 11) continue;
              if(debug) cout << "third lepton found" << endl;

              if(tracks[j].second/abs(tracks[j].second) == -tracks[k].second/abs(tracks[k].second)) {
                //Kaon and lepton have same charge
                //correct mass assumption
                if(debug) cout << "correct mass assumption" << endl;
                allPlots["massD0_lep"+chTag]->Fill(mass12,wgt);
                allPlots["massD0_lep"+chTag+"_no_weight"]->Fill(mass12,norm);

                if(abs(tracks[k].second) == 13)
                  allPlots["massD0_mu"+chTag]->Fill(mass12,wgt);
                if(abs(tracks[k].second) == 11)
                  allPlots["massD0_e"+chTag]->Fill(mass12,wgt);
                if(abs(tracks[k].second) == 13)
                  allPlots["massD0_mu"+chTag+"_no_weight"]->Fill(mass12,norm);
                if(abs(tracks[k].second) == 11)
                  allPlots["massD0_e"+chTag+"_no_weight"]->Fill(mass12,norm);
              }
            }
            //looking for pion
            if(debug) cout << "D*->pi+D0" << endl;
            for(size_t k = 0; k < tracks.size(); k++) {
              if(k == i) continue;
              if(k == j) continue;

              if(abs(tracks[k].second) != 211) continue;
              if(debug) cout << "Pion found" << endl;

              TLorentzVector p_track3, p_cand;
              p_track3.SetPtEtaPhiM(tracks[k].first.Pt(), tracks[k].first.Eta(), tracks[k].first.Phi(), gMassPi);
              if(debug) cout << k << ": " << tracks[k].first.Pt() << " " << tracks[k].first.Eta() << " " << tracks[k].first.Phi() << " " << gMassPi << endl;
              allPlots["pi_pt"+chTag]->Fill(p_track3.Pt(),wgt);
              allPlots["pi_pt"+chTag+"_no_weight"]->Fill(p_track3.Pt(),norm);
              if( tracks[j].second/abs(tracks[j].second) == -tracks[k].second/abs(tracks[k].second) ) {
                // Kaon and pion have opposite charges
                // I.e. correct mass assumption
                if(debug) cout << "correct mass assumption" << endl;

                p_cand = p_track1+p_track2+p_track3;
                allPlots["massDs"+chTag]->Fill(p_cand.M(), wgt);
                allPlots["massDs"+chTag+"_no_weight"]->Fill(p_cand.M(),norm);
                allPlots["massDs_all"]->Fill(p_cand.M(), wgt);

                if(abs(mass12-1.864) < 0.10) { // mass window cut
                  TLorentzVector p_jet;
                  p_jet.SetPtEtaPhiM(ev.j_pt[jetindex], ev.j_eta[jetindex], ev.j_phi[jetindex], 0.);

                  float deltam = p_cand.M() - mass12;

                  allPlots["massDsmD0loose"+chTag]->Fill(deltam, wgt);
                  allPlots["massDsmD0loose"+chTag+"_no_weight"]->Fill(deltam,norm);
                  allPlots["massDsmD0loose_all"]->Fill(deltam, wgt);
                  if(abs(mass12-1.864) < 0.05) { // tighter mass window cut

                    allPlots["massDsmD0"+chTag]->Fill(deltam, wgt);
                    allPlots["massDsmD0"+chTag+"_no_weight"]->Fill(deltam, norm);
                    allPlots["massDsmD0_all"]->Fill(deltam, wgt);
                    allPlots["lp_pt"+chTag+"_meson"]->Fill(leptons[0].Pt(),wgt);
                    allPlots["lp_pt_low"+chTag+"_meson"]->Fill(leptons[0].Pt(),wgt);
                    if(deltam<0.14 || deltam>0.15) continue;
                    allPlots["nevt"+chTag+"_meson"]->Fill(1,norm);
                    allPlots["weight"+chTag+"_meson"]->Fill(wgt,norm);
                    allPlots["norm"+chTag+"_meson"]->Fill(norm,norm);
                    allPlots["nvtx"+chTag+"_meson"]->Fill(ev.nvtx,norm);
                  }
                }
              }
            }
          }
        if(debug) cout << "D0 and D* DONE" << endl;
      }

    }

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

  for (auto& it : allPlots)  { 
    if(debug) cout << it.second->GetName() << endl;
    if(debug) cout << it.second->GetEntries() << endl;

    //fOut->cd( dir );
    it.second->SetDirectory(fOut); it.second->Write(); 
    fOut->cd();
  }
  if(debug) cout << "writing histograms DONE" << endl;
  if(debug) cout << "closing ROOT file" << endl;
  fOut->Close();
}

