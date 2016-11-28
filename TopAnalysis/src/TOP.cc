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
  //bool debug(false);
  //bool isData( filename.Contains("Data") ? true : false);
  
  //READ TREE FROM FILE
  MiniEvent_t ev;
  //TopWidthEvent_t ev;
  TFile *f = TFile::Open(filename);
  //TTree *t = (TTree*)f->Get("twev");
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  //createTopWidthEventTree(t,ev);
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
      //TString puWgtUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data"+era+"/pileupWgts.root");
      TString puWgtUrl(era+"/pileupWgts.root");
      gSystem->ExpandPathName(puWgtUrl);
      TFile *fIn=TFile::Open(puWgtUrl);
      //TGraph *puData=(TGraph *)fIn->Get(grName);
      //Float_t totalData=puData->Integral();
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
  //TString lepEffUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/leptonEfficiencies.root");
  //TString lepEffUrl(era+"/leptonEfficiencies.root");
  //FIXME
  //TString lepEffUrl(era+"/muonEfficiencies.root");
  //gSystem->ExpandPathName(lepEffUrl);
  //std::map<TString,TH2 *> lepEffH;
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);
  /*
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH["m_sel"]=(TH2 *)fIn->Get("m_sel");
      lepEffH["m_trig"]=(TH2 *)fIn->Get("m_trig");      
      for(auto& it : lepEffH) it.second->SetDirectory(0);
      fIn->Close();
    }
  */

  //lepEffUrl="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CutBasedID_TightWP_76X_18Feb.txt_SF2D.root";
  //lepEffUrl=era+"/CutBasedID_TightWP_76X_18Feb.txt_SF2D.root";
  /*
  lepEffUrl=era+"/electronEfficiencies.root";
  gSystem->ExpandPathName(lepEffUrl);
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH["e_sel"]=(TH2 *)fIn->Get("EGamma_SF2D");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
      fIn->Close();
    }
  */

  //B-TAG CALIBRATION
  //TString btagUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CSVv2.csv");
  TString btagUncUrl(era+"/btagSFactors.csv");
  gSystem->ExpandPathName(btagUncUrl);
  std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
  //TString btagEffExpUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/expTageff.root");
  TString btagEffExpUrl(era+"/expTageff.root");
  gSystem->ExpandPathName(btagEffExpUrl);
  std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
  BTagSFUtil myBTagSFUtil;
  //float wgt(1.0);
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
      //expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      //expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      //expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
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

      //wgt=1.0;
      //float norm( normH ? normH->GetBinContent(1) : 1.0);
      //wgt=norm;//lepTriggerSF[0]*lepSelSF[0]*puWeight[0]*norm;
    }

  //JET ENERGY SCALE: https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Summer15_uncertainties
  //jet energy uncertainties
  TString jecUncUrl(era+"/jecUncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  //JetCorrectorParameters *jecParam = new JetCorrectorParameters(jecUncUrl.Data(),"Total");
  //JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty( *jecParam );
  
  //LIST OF SYSTEMATICS
  
  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",4,0,4);
  //addGenScanCounters(allPlots,f); FIXME
  std::vector<TString> lfsVec = { "_all", "_e", "_ee", "_em", "_mm", "_m" }; 
  std::vector<TString> cutVec = { "", "_lep", "_jpsi", "_csv", "_meson" };
  std::vector<TString> wgtVec = { "", "_no_weight" };

  for(int i = 0; i < (int)lfsVec.size(); i++) {
  for(int j = 0; j < (int)cutVec.size(); j++) {
  for(int k = 0; k < (int)wgtVec.size(); k++) {
    TString tag(lfsVec[i]);
    TString cut(cutVec[j]);
    TString weight(wgtVec[k]);
    allPlots["lp_pt"+tag+cut+weight] = new TH1F("lp_pt"+tag+cut+weight,";Leading Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["l2p_pt"+tag+cut+weight] = new TH1F("l2p_pt"+tag+cut+weight,";Sub-leading Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_pt"+tag+cut+weight] = new TH1F("dilp_pt"+tag+cut+weight,";Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_m"+tag+cut+weight] = new TH1F("dilp_m"+tag+cut+weight,";M_{ll} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["j_pt"+tag+cut+weight] = new TH1F("j_pt"+tag+cut+weight,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["bj_pt"+tag+cut+weight] = new TH1F("bj_pt"+tag+cut+weight,";Leading b Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["nlp"+tag+cut+weight]     = new TH1F("nlp"+tag+cut+weight,";N_{l};Events" ,3,0.,3.);
    allPlots["ndilp"+tag+cut+weight]     = new TH1F("ndilp"+tag+cut+weight,";N_{ll};Events" ,3,0.,3.);
    allPlots["nj"+tag+cut+weight]     = new TH1F("nj"+tag+cut+weight,";N_{jets} (P_{T} > 30 GeV);Events" ,8,2.,10.);
    allPlots["nbj"+tag+cut+weight]     = new TH1F("nbj"+tag+cut+weight,";N_{b-jets} (CSV > 0.8);Events" ,4,1.,5.);
    allPlots["npf"+tag+cut+weight]     = new TH1F("npf"+tag+cut+weight,";N_{pf};Events / 10" ,5,0.,5.);
    allPlots["nstart"+tag+cut+weight]     = new TH1F("jetindex"+tag+cut+weight,";N_{jetindex};Events" ,5,0.,5.);
    allPlots["pfid"+tag+cut+weight]     = new TH1F("pfid"+tag+cut+weight,";PFID;Events" ,440,-220.,220.);
/*
    allPlots["massJPsi"+tag+cut+weight]     = new TH1F("massJPsi"+tag+cut+weight,";M_{J/#Psi};Events" ,20,2.,4.);
    allPlots["massD0"+tag+cut+weight]     = new TH1F("massD0"+tag+cut+weight,";M_{jj};Events" ,20,1.,3.);
    allPlots["massDs"+tag+cut+weight]     = new TH1F("massDs"+tag+cut+weight,";M_{D^{*}};Events" ,20,0.,20.);
*/
    //allPlots["massJPsi"+tag+cut+weight]     = new TH1F("massJPsi"+tag+cut+weight,";M_{J/#Psi};Events / 0.01 GeV" ,18,2.5,3.4);
    //allPlots["massJPsi"+tag+cut+weight]     = new TH1F("massJPsi"+tag+cut+weight,";M_{J/#Psi};Events / 0.5 GeV" ,20,0,10);
    allPlots["massJPsi"+tag+cut+weight]     = new TH1F("massJPsi"+tag+cut+weight,";M_{ll};Events / 36 MeV" ,25,2.5,3.4);
    allPlots["massJPsiK"+tag+cut+weight]     = new TH1F("massJPsiK"+tag+cut+weight,";M_{llk};Events / 15 MeV" ,100,4.5,6);
    allPlots["massD0"+tag+cut+weight]     = new TH1F("massD0"+tag+cut+weight,";M_{D^{0}};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_lep"+tag+cut+weight]     = new TH1F("massD0_lep"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_mu"+tag+cut+weight]     = new TH1F("massD0_mu"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_e"+tag+cut+weight]     = new TH1F("massD0_ele"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massDsmD0loose"+tag+cut+weight]     = new TH1F("massDsmD0loose"+tag+cut+weight,";M_{K^{-}#pi^{+}#pi^{+}} - M_{K^{-}#pi^{+}};Events / 0.6 MeV" ,25,0.14,0.17);
    allPlots["massDsmD0"+tag+cut+weight]     = new TH1F("massDsmD0"+tag+cut+weight,";M_{K^{-}#pi^{+}#pi^{+}} - M_{K^{-}#pi^{+}};Events / 0.6 MeV" ,25,0.14,0.17);
    allPlots["massDs"+tag+cut+weight]     = new TH1F("massDs"+tag+cut+weight,";M_{D^{*}};Events / 10 MeV" ,200,0.,2.0);
    allPlots["pi_pt"+tag+cut+weight] = new TH1F("pi_pt"+tag+cut+weight,";#pi^{#pm} P_{T} [GeV];Events / 5 GeV", 10, 0,50);
    allPlots["MET"+tag+cut+weight] = new TH1F("MET"+tag+cut+weight,";MET [GeV];Events / 20 GeV", 10,0,200);
    allPlots["charge"+tag+cut+weight] = new TH1F("charge"+tag+cut+weight,";Charge(l_{1}*l_{2});Events", 5,-2,2);
    allPlots["csv"+tag+cut+weight] = new TH1F("CSV"+tag+cut+weight,";Jet CSV;Events / 0.01", 80,0.8,1);
    allPlots["dR"+tag+cut+weight] = new TH1F("dR"+tag+cut+weight,";dR;Events / 0.05", 20,0.0,1.);
    allPlots["pflp_pt"+tag+cut+weight] = new TH1F("pflp_pt"+tag+cut+weight,";PF lepton P_{T} [GeV];Events / 0.2 GeV", 15, 0,3);
    allPlots["massZ"+tag+cut+weight]     = new TH1F("massZ_control"+tag+cut+weight,";M_{ll};Events / 1.0 GeV" ,30,81,111);
    allPlots["nevt"+tag+cut+weight]     = new TH1F("nevt"+tag+cut+weight,";N_{events};Events" ,1,1.,2.);
    allPlots["pf_dxy"+tag+cut+weight] = new TH1F("pf_dxy"+tag+cut+weight,";d_{xy} [cm];Events / 20 #mum", 100, 0, 0.1);
    allPlots["pf_dz"+tag+cut+weight] = new TH1F("pf_dz"+tag+cut+weight,";d_{z} [cm];Events / 20 #mum", 100, 0, 0.1);
    allPlots["pf_dxyE"+tag+cut+weight] = new TH1F("pf_dxyE"+tag+cut+weight,";#sigma(d_{xy}) [cm];Events / 20 #mum", 100, 0, 0.1);
    allPlots["pf_dzE"+tag+cut+weight] = new TH1F("pf_dzE"+tag+cut+weight,";#sigma(d_{z}) [cm];Events / 20 #mum", 100, 0, 0.1);
    allPlots["pf_dxy_sig"+tag+cut+weight] = new TH1F("pf_dxy_significance"+tag+cut+weight,";d_{xy};Events / 0.3", 100, 0, 30);
    allPlots["pf_dz_sig"+tag+cut+weight] = new TH1F("pf_dz_significance"+tag+cut+weight,";d_{z};Events / 0.3", 100, 0, 30);

  }
  }
  }
    allPlots["relIso_m"] = new TH1F("relIso_m",";relIso;Events / 0.05", 20,0,1.);
    allPlots["relIso_e"] = new TH1F("relIso_e",";relIso;Events / 0.05", 20,0,1.);
    allPlots["nevt_iso"] = new TH1F("nevt_iso",";After Isolation;Events", 1,0,1.);
    allPlots["nevt_veto"] = new TH1F("nevt_veto",";After Veto;Events", 1,0,1.);


  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  //for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }
  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%5000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));
      allPlots["nevt_all"]->Fill(1,1);

      //account for pu weights and effect on normalization
      //float puWeight(1.0);
      if(!ev.isData) 
	{
	  //puWeight=puWgtGr[0]->Eval(ev.putrue);  
          /*
	  allPlots["puwgtctr"]->Fill(0.,1.0);
	  allPlots["puwgtctr"]->Fill(1.,puWeight);
          */
	}

      //select 1 good lepton
      //cout << "entering lepton selection" << endl;
      //std::vector<int> tightLeptonsNonIso;
      std::vector<int> tightLeptons,vetoLeptons;
      for(int il=0; il<ev.nl; il++)
	{
          //cout << "in lepton selection" << endl;
	  bool passTightKin(ev.l_pt[il]>20 && fabs(ev.l_eta[il])<2.4); // TOP mu cut for dilep

	  float relIso(ev.l_relIso[il]);
	  bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);
	  bool passIso( ev.l_id[il]==13 ? relIso<0.25 : (ev.l_pid[il]>>1)&0x1 ); // TOP mu cut for dilep FIXME
	  
	  //bool passNonIso(relIso>0.4); //FIXME from 7_6_x
	  //if( ev.l_id[il]==11 && (passIso || relIso<0.4) ) passNonIso=false; //FIXME from 7_6_x

	  bool passVetoIso(  ev.l_id[il]==13 ? relIso<0.25 : true); //FIXME from 7_6_x
          bool passVetoKin(  ev.l_pt[il]>10 && fabs(ev.l_eta[il])<2.5); // TOP veto

	  //bool passSIP3d(ev.l_ip3dsig[il]<4);
	  //if(channelSelection==21) passSIP3d=true;
	  if( ev.l_id[il] == 13 )
            allPlots["relIso_m"]->Fill(relIso,1);
          else if( ev.l_id[il] == 11)
            allPlots["relIso_e"]->Fill(relIso,1);

	  if(passTightKin && passTightId)// && passSIP3d)
	    {
	      if(passIso)         tightLeptons.push_back(il);
	      //else if(passNonIso) tightLeptonsNonIso.push_back(il); //FIXME from 7_6_x
	    }
	  else if(passVetoKin && passVetoIso) vetoLeptons.push_back(il); //FIXME from 7_6_x
	}
      if(debug) cout << "lepton selection DONE" << endl;

      //check if triggers have fired
      bool hasEETrigger(((ev.elTrigger>>1)&0x1)!=0 || ((ev.elTrigger>>4)&0x1)!=0);
      bool hasMMTrigger(((ev.muTrigger>>2)&0x3)!=0);
      bool hasEMTrigger(((ev.elTrigger>>2)&0x3)!=0);
      bool hasMuTrigger((ev.muTrigger & 0x3)!=0);
      bool hasEleTrigger((ev.elTrigger & 0x1)!=0);
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
	    passTightKin = (ev.l_pt[tightLeptons[0]] > 26 && fabs(ev.l_eta[tightLeptons[0]])<2.1); // TOP mu cut for dilep
            passIso = (ev.l_relIso[tightLeptons[0]] < 0.15); //TOP mu cut for lep+jets
          }
          else if(ev.l_id[tightLeptons[0]]==11) { // electron + jets
            passTightKin = (ev.l_pt[tightLeptons[0]] > 30); //from TOP-15-005
            passIso = (ev.l_relIso[tightLeptons[0]] < 0.15); //TOP mu cut for lep+jets
          }
          //************************************
          if(passTightKin && passIso) {
	    selLeptons.push_back( tightLeptons[0] );
            if(debug) cout << "found 1 tight lepton" << endl;
          }
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
	  if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) chTag="e";
	  if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger)  chTag="m";
	}
      if(selLeptons.size()==2)
	{
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*11 && hasEETrigger) chTag="ee";
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==13*13 && hasMMTrigger) chTag="mm";
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*13)
	    {
	      //if(tightLeptons.size()>=2 && (hasEleTrigger || hasMuTrigger)) chTag="em";
	      if(selLeptons.size()>=2 && (hasEMTrigger)) chTag="em";
	      if(selLeptons.size()==1)
		{
		  if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) chTag="em";
		  if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger) chTag="em";
		}
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
	     fabs(dilp4.M()-91)<10) //&& 
	     //dilp4.Pt()>30)
	    { 
	      isZ=true; 
	      //isZPassingSIP3d=(ev.l_ip3dsig[0]<4 && ev.l_ip3dsig[1]<4);
	    }
	  lepIdx=selLeptons[0];
          if(debug) cout << "di-lepton DONE" << endl;
	}

      if(lepIdx<0) continue;
      allPlots["nevt_iso"]->Fill(1);
      
      //no extra isolated leptons
      if(vetoLeptons.size()>0) continue;
      allPlots["nevt_veto"]->Fill(1);
      
      //apply trigger requirement
      /*
      if(ev.l_id[lepIdx]==13)
	{
	  if(ev.isData  && (ev.muTrigger & 0x3)==0) continue;
	  if(!ev.isData && (ev.muTrigger & 0x3)==0) continue;
	}
      if(ev.l_id[lepIdx]==11)
	{ 
	  if( ((ev.elTrigger>>0)&0x1)==0 ) continue;
	}

      //select according to the lepton id/charge
      Int_t lid=ev.l_id[lepIdx];
      if(isZ) lid=2100+ev.l_id[lepIdx];
      else if(tightLeptonsNonIso.size()==1) lid=100*ev.l_id[lepIdx];
      if(channelSelection!=0)
	{
	  if(channelSelection==21) { if(!isZ) continue; }
	  else                     { if(lid!=channelSelection) continue; }
	}
      if(chargeSelection!=0  && ev.l_charge[lepIdx]!=chargeSelection) continue;
      */

      //lepton kinematics
      if(debug) cout << "checking lepton kinematics" << endl;
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);

      if(debug) cout << "checking lepton kinematics DONE" << endl;

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
      std::vector<TLorentzVector> bJets,lightJets;
      TLorentzVector visSystem(isZ ? dilp4 : lp4);
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

	  //cross clean with respect to leptons 
          bool overlapsWithLepton(false);
          for(size_t il=0; il<leptons.size(); il++) {
            if(jp4.DeltaR(leptons[il])>0.4) continue;
	    overlapsWithLepton=true;
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

	  //require back-to-back configuration with Z
	  //if(isZ && jp4.DeltaPhi(dilp4)<2.7) continue; //FIXME

	  // re-inforce kinematics cuts
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  //if(leadingJetIdx<0) leadingJetIdx=k;
	  htsum += jp4.Pt();
	  if(bJets.size()+lightJets.size()<4) visSystem += jp4;

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
	  if(isBTagged) bJets.push_back(jp4);
	  else          lightJets.push_back(jp4);

          Jet tmpj(jp4, csv, k);
	  for(int ipf = 0; ipf < ev.npf; ipf++) {
	    if(ev.pf_j[ipf] != k) continue;
	    if(ev.pf_c[ipf]==0) continue;
	    TLorentzVector tkP4(0,0,0,0);
	    tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);
            pfTrack pftk(tkP4, ev.pf_dxy[ipf], ev.pf_dxyE[ipf], ev.pf_dz[ipf], ev.pf_dzE[ipf], ev.pf_id[ipf]);
	    //tmpj.addTrack(tkP4,ev.pf_id[ipf]);
	    tmpj.addTrack(pftk,ev.pf_id[ipf]);
            /*
	    tmpj.addTrack(ipf);
	    tmpj.addDxy(ev.pf_dxy[ipf], ev.pf_dxyE[ipf]);
	    tmpj.addDz(ev.pf_dz[ipf], ev.pf_dzE[ipf]);
            */
	  }
          tmpj.sortTracksByPt();

          if(isBTagged) bJetsVec.push_back(tmpj);
          else lightJetsVec.push_back(tmpj);
          allJetsVec.push_back(tmpj);
          /*
          if(isBTagged) allPlots["csv"+chTag]->Fill(csv,wgt);
          if(isBTagged) allPlots["csv_all"]->Fill(csv,wgt);
          */
	}

      
      //event weight
      float wgt(1.0);
      float norm(1.0);
      std::vector<float> puWgts(3,1.0),topPtWgts(2,1.0);
      EffCorrection_t lepSelCorrWgt(1.0,0.0), triggerCorrWgt(1.0,0.0);
      if(debug) cout << "Lepton scale factors" << endl;
      if(!ev.isData)
	{
	  //update lepton selection scale factors, if found
	  //float lepTriggerSF(1.0),lepSelSF(1.0);
          //FIXME

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
          std::vector<int> pdgIds;
          for(size_t ilp = 0; ilp < selLeptons.size(); ilp++)
            pdgIds.push_back(ev.l_id[selLeptons[ilp]]);
	  //triggerCorrWgt=lepEffH.getTriggerCorrection(selLeptons,leptons);
	  triggerCorrWgt=lepEffH.getTriggerCorrection(pdgIds,leptons);
          if(debug) cout << "calling trigger function DONE!" << endl;
          // ** selLeptons contains only ev_l position, leptons contains p4 **
	  for(size_t il=0; il<selLeptons.size(); il++) {
	    EffCorrection_t selSF=lepEffH.getOfflineCorrection(ev.l_id[selLeptons[il]],leptons[il].Pt(),leptons[il].Eta());
	    lepSelCorrWgt.second = sqrt( pow(lepSelCorrWgt.first*selSF.second,2)+pow(lepSelCorrWgt.second*selSF.first,2));
            if(debug) cout << "lepSelCorrWgt=" << lepSelCorrWgt.first << endl;
            if(debug) cout << "selSF=" << selSF.first << endl;
	    lepSelCorrWgt.first *= selSF.first;
            //if(lepSelCorrWgt.first < 0.7) lepSelCorrWgt.first = selSF.first;
	   }
          /*
	  for(UInt_t il=0; il<leptons.size(); il++)
	    {
	      TString prefix(abs(ev.l_id[il])==11 ? "e" : "m");
	      float minEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmin() ), maxEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmax()-0.01 );
	      float etaForEff=TMath::Max(TMath::Min(float(fabs(leptons[il].Eta())),maxEtaForEff),minEtaForEff);
	      Int_t etaBinForEff=lepEffH[prefix+"_sel"]->GetXaxis()->FindBin(etaForEff);
	      
	      float minPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmin() ), maxPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmax()-0.01 );
	      float ptForEff=TMath::Max(TMath::Min(float(leptons[il].Pt()),maxPtForEff),minPtForEff);
	      Int_t ptBinForEff=lepEffH[prefix+"_sel"]->GetYaxis()->FindBin(ptForEff);
		  		  
	      lepSelSF=(lepEffH[prefix+"_sel"]->GetBinContent(etaBinForEff,ptBinForEff));
	      
	      if(il!=0) continue;
	      if(prefix=="m")
		{
		  lepTriggerSF=(lepEffH[prefix+"_trig"]->GetBinContent(etaBinForEff,ptBinForEff));
		
		}
	    }
          */
	  
	  //https://indico.cern.ch/event/539804/contributions/2196937/attachments/1291290/1923328/TopTriggers2016v2.pdf
	  /*
	  if(era.Contains("2015"))
	    {
	      if(chTag=="EE") lepTriggerSF=0.930;
	      if(chTag=="MM") lepTriggerSF=0.894;
	      if(chTag=="EM") lepTriggerSF=0.931;
	    }
	  else
	    {
	      if(chTag=="EE") lepTriggerSF=0.93;
              if(chTag=="MM") lepTriggerSF=0.87;
              if(chTag=="EM") lepTriggerSF=0.88;
	    }
          */

	  //update nominal event weight
	  //float norm( normH ? normH->GetBinContent(1) : 1.0);
	  norm =  normH ? normH->GetBinContent(1) : 1.0;
	  //wgt=lepTriggerSF*lepSelSF*puWeight*norm;
	  wgt=triggerCorrWgt.first*lepSelCorrWgt.first*puWgts[0]*norm;
	  if(ev.ttbar_nw>0) wgt*=ev.ttbar_w[0];
          if(debug) cout << "weight=" << wgt << endl;
          if(debug) cout << "Trigger=" << triggerCorrWgt.first << endl << "Lepton=" << lepSelCorrWgt.first << endl << "PU=" << puWgts[0] << endl << "norm=" << norm  << endl;
          if(filename.Contains("_WJets")) cout << "Trigger=" << triggerCorrWgt.first << endl << "Lepton=" << lepSelCorrWgt.first << endl << "PU=" << puWgts[0] << endl << "norm=" << norm  << endl << "ttbar_w[0]=" << ev.ttbar_w[0] << endl << "wgt=" << wgt << endl;
          for(size_t il = 0; il < leptons.size(); il++) {
            if(!filename.Contains("_WJets")) continue;
            cout << "pT: " << leptons[il].Pt() << endl;
          }
          //wgt=1.0;
	}
      if(debug) cout << "Lepton scale factors DONE!" << endl;

      //sort by Pt
      sort(bJetsVec.begin(),    bJetsVec.end(),   sortJetsByPt);
      sort(allJetsVec.begin(),  allJetsVec.end(), sortJetsByPt);

      for(size_t ij = 0; ij < bJetsVec.size(); ij++) {
        float csv = bJetsVec.at(ij).getCSV();
        allPlots["csv"+chTag]->Fill(csv,wgt);
        allPlots["csv_all"]->Fill(csv,wgt);
      }

      if(bJetsVec.size()==0) continue;
      for(size_t il=0; il<leptons.size(); il++) {
        for(size_t ij=0; ij<bJets.size(); ij++) {
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[ij],ev.j_eta[ij],ev.j_phi[ij],ev.j_mass[ij]);
          allPlots["dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
          allPlots["dR"+chTag+"_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
        }
        for(size_t ij=0; ij<lightJets.size(); ij++) {
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[ij],ev.j_eta[ij],ev.j_phi[ij],ev.j_mass[ij]);
          allPlots["dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
          allPlots["dR_all"]->Fill(jp4.DeltaR(leptons[il]),wgt);
          allPlots["dR"+chTag+"_no_weight"]->Fill(jp4.DeltaR(leptons[il]),norm);
        }
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
      if(debug) cout << "sorting jets" << endl;
      std::sort(lightJets.begin(), lightJets.end(), VecSort);
      std::sort(bJets.begin(), bJets.end(), VecSort);
      if(debug) cout << "starting simple plots" << endl;
      if(selLeptons.size() == 1 && bJetsVec.size() >= 1 && lightJetsVec.size() >= 4) {
        singleLep = true;
        allPlots["nj"+chTag]->Fill(lightJetsVec.size(),wgt);
        allPlots["nbj"+chTag]->Fill(bJetsVec.size(),wgt);
        allPlots["nj"+chTag+"_no_weight"]->Fill(lightJetsVec.size(),norm);
        allPlots["nbj"+chTag+"_no_weight"]->Fill(bJetsVec.size(),norm);
        allPlots["nj_all"]->Fill(lightJetsVec.size(),wgt);
        allPlots["nbj_all"]->Fill(bJetsVec.size(),wgt);
        allPlots["nlp"+chTag]->Fill(selLeptons.size(),wgt);
        allPlots["nlp"+chTag+"_no_weight"]->Fill(selLeptons.size(),norm);
        allPlots["nevt"+chTag]->Fill(1,wgt);
        std::sort(leptons.begin(), leptons.end(), VecSort);
        allPlots["lp_pt"+chTag]->Fill(leptons[0].Pt(),wgt);
        allPlots["lp_pt"+chTag+"_no_weight"]->Fill(leptons[0].Pt(),norm);
        allPlots["lp_pt_all"]->Fill(leptons[0].Pt(),wgt);
        allPlots["j_pt"+chTag]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
        allPlots["bj_pt"+chTag]->Fill(bJetsVec[0].getVec().Pt(),wgt);
        allPlots["j_pt_all"]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
        allPlots["bj_pt_all"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
        allPlots["MET"+chTag]->Fill(ev.met_pt[0],wgt);
        allPlots["j_pt"+chTag+"_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
        allPlots["bj_pt"+chTag+"_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
        allPlots["MET"+chTag+"_no_weight"]->Fill(ev.met_pt[0],norm);
        for(int i = 0; i < ev.npf; i++) {
          allPlots["pfid"+chTag]->Fill(ev.pf_id[i],wgt);
        }
      }
      else if(selLeptons.size() == 2 && bJets.size() >= 1 && lightJetsVec.size() > 1) {
        if(isZ) {
          allPlots["massZ"+chTag]->Fill(dilp4.M(),wgt);
          allPlots["massZ"+chTag+"_no_weight"]->Fill(dilp4.M(),norm);
        }
        if(isZ) continue;
        if(dilp4.M() < 10) continue; // && ev.l_charge[selLeptons[0]]!=ev.l_charge[selLeptons[1]]) continue;
        if(ev.l_id[selLeptons[0]]==ev.l_id[selLeptons[1]] && met.Pt() < 40) continue; //FIXME
        doubleLep = true;
        allPlots["nj"+chTag]->Fill(lightJetsVec.size(),wgt);
        allPlots["nbj"+chTag]->Fill(bJetsVec.size(),wgt);
        allPlots["ndilp"+chTag]->Fill(selLeptons.size(),wgt);
        allPlots["dilp_pt"+chTag]->Fill(dilp4.Pt(),wgt);
        allPlots["dilp_m"+chTag]->Fill(dilp4.M(),wgt);
        allPlots["nj"+chTag+"_no_weight"]->Fill(lightJetsVec.size(),norm);
        allPlots["nbj"+chTag+"_no_weight"]->Fill(bJetsVec.size(),norm);
        allPlots["ndilp"+chTag+"_no_weight"]->Fill(selLeptons.size(),norm);
        allPlots["dilp_pt"+chTag+"_no_weight"]->Fill(dilp4.Pt(),norm);
        allPlots["dilp_m"+chTag+"_no_weight"]->Fill(dilp4.M(),norm);
        std::sort(leptons.begin(), leptons.end(), VecSort);
        allPlots["lp_pt"+chTag]->Fill(leptons[0].Pt(),wgt);
        allPlots["l2p_pt"+chTag]->Fill(leptons[1].Pt(),wgt);
        allPlots["lp_pt"+chTag+"_no_weight"]->Fill(leptons[0].Pt(),norm);
        allPlots["l2p_pt"+chTag+"_no_weight"]->Fill(leptons[1].Pt(),norm);
        allPlots["j_pt"+chTag]->Fill(lightJetsVec[0].getVec().Pt(),wgt);
        allPlots["bj_pt"+chTag]->Fill(bJetsVec[0].getVec().Pt(),wgt);
        allPlots["j_pt"+chTag+"_no_weight"]->Fill(lightJetsVec[0].getVec().Pt(),norm);
        allPlots["bj_pt"+chTag+"_no_weight"]->Fill(bJetsVec[0].getVec().Pt(),norm);
        allPlots["MET"+chTag]->Fill(ev.met_pt[0],wgt);
        allPlots["charge"+chTag]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
        allPlots["MET"+chTag+"_no_weight"]->Fill(ev.met_pt[0],norm);
        allPlots["charge"+chTag+"_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],norm);
        for(int i = 0; i < ev.npf; i++) {
          allPlots["pfid"+chTag]->Fill(ev.pf_id[i],wgt);
        }
      }
      if(debug) cout << "simple plots DONE" << endl;


      if(!singleLep && !doubleLep) continue;

      allPlots["npf"+chTag+"_lep"]->Fill(ev.npf,wgt);
      allPlots["npf"+chTag+"_lep"+"_no_weight"]->Fill(ev.npf,norm);
      allPlots["nevt"+chTag+"_lep"]->Fill(1,wgt);
      allPlots["nevt_all_lep"]->Fill(1,wgt);

      //charmed resonance analysis : use only jets with CSV>CSVL, up to two per event
      for(size_t ij = 0; ij < bJetsVec.size(); ij++) {

        if(ij > 1) continue;
        if(ij == 0) {
          allPlots["nbj"+chTag+"_csv"]->Fill(1,wgt);
          allPlots["nevt"+chTag+"_csv"]->Fill(1,wgt);
        }
        allPlots["bj_pt"+chTag+"_csv"]->Fill(bJetsVec[ij].getVec().Pt(),wgt);

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
          allPlots["csv"+chTag+"_jpsi"]->Fill(bJetsVec[ij].getCSV(),wgt);
          allPlots["nevt"+chTag+"_jpsi"]->Fill(1,wgt);
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

          if(filename.Contains("_WJets"))
            cout << endl << mass12 << " " << wgt << endl;
          if(mass123 > 0) {
            allPlots["massJPsiK"+chTag]->Fill(mass123,wgt);
            allPlots["massJPsiK_all"]->Fill(mass123,wgt);
          }
          /*
          if(mass123 > 0) kaonCands.erase(kaonCands.begin());
          pfmuCands.erase(pfmuCands.begin()+1);
          pfmuCands.erase(pfmuCands.begin());
          */
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
                    if(filename.Contains("_WJets")) {
                      cout << endl << deltam << " " << wgt << endl;
                      cout << "pi1: " << p_track1.Pt() << endl;
                      cout << "K: " << p_track2.Pt() << endl;
                      cout << "pi2: " << p_track3.Pt() << endl;
                    }

                    allPlots["massDsmD0"+chTag]->Fill(deltam, wgt);
                    allPlots["massDsmD0"+chTag+"_no_weight"]->Fill(deltam, norm);
                    allPlots["massDsmD0_all"]->Fill(deltam, wgt);
                    if(deltam<0.14 || deltam>0.15) continue;
/*
                    allPlots["pf_dxy"+chTag+"_meson"]->Fill(abs(tracks[k].first.getDxy()),wgt);
                    allPlots["pf_dz"+chTag+"_meson"]->Fill(abs(tracks[k].first.getDz()),wgt);
                    allPlots["pf_dxyE"+chTag+"_meson"]->Fill(abs(tracks[k].first.getDxyE()),wgt);
                    allPlots["pf_dzE"+chTag+"_meson"]->Fill(abs(tracks[k].first.getDzE()),wgt);
                    allPlots["pf_dxy_sig"+chTag+"_meson"]->Fill(abs(tracks[k].first.getDxy())/abs(tracks[k].first.getDxyE()),wgt);
                    allPlots["pf_dz_sig"+chTag+"_meson"]->Fill(abs(tracks[k].first.getDz())/abs(tracks[k].first.getDzE()),wgt);
                    allPlots["pf_dxy_all"]->Fill(abs(tracks[k].first.getDxy()),wgt);
                    allPlots["pf_dz_all"]->Fill(abs(tracks[k].first.getDz()),wgt);
                    allPlots["pf_dxyE_all"]->Fill(abs(tracks[k].first.getDxyE()),wgt);
                    allPlots["pf_dzE_all"]->Fill(abs(tracks[k].first.getDzE()),wgt);
                    allPlots["pf_dxy_sig_all"]->Fill(abs(tracks[k].first.getDxy())/abs(tracks[k].first.getDxyE()),wgt);
                    allPlots["pf_dz_sig_all"]->Fill(abs(tracks[k].first.getDz())/abs(tracks[k].first.getDzE()),wgt);
*/
                    allPlots["nevt"+chTag+"_meson"]->Fill(1,wgt);
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
