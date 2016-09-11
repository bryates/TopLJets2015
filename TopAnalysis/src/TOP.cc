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
  float wgt(1.0);
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
    allPlots["massJPsi"+tag+cut+weight]     = new TH1F("massJPsi"+tag+cut+weight,";M_{ll};Events / 9 MeV" ,100,2.5,3.4);
    allPlots["massJPsiK"+tag+cut+weight]     = new TH1F("massJPsiK"+tag+cut+weight,";M_{llk};Events / 15 MeV" ,100,4.5,6);
    allPlots["massD0"+tag+cut+weight]     = new TH1F("massD0"+tag+cut+weight,";M_{D^{0}};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_lep"+tag+cut+weight]     = new TH1F("massD0_lep"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_mu"+tag+cut+weight]     = new TH1F("massD0_mu"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massD0_e"+tag+cut+weight]     = new TH1F("massD0_ele"+tag+cut+weight,";M_{K#pi};Events / 3 MeV" ,100,1.7,2.0);
    allPlots["massDsmD0loose"+tag+cut+weight]     = new TH1F("massDsmD0loose"+tag+cut+weight,";M_{K^{-}#pi^{+}#pi^{+}} - M_{K^{-}#pi^{+}};Events / 0.3 MeV" ,100,0.14,0.17);
    allPlots["massDsmD0"+tag+cut+weight]     = new TH1F("massDsmD0"+tag+cut+weight,";M_{K^{-}#pi^{+}#pi^{+}} - M_{K^{-}#pi^{+}};Events / 0.3 MeV" ,100,0.14,0.17);
    allPlots["massDs"+tag+cut+weight]     = new TH1F("massDs"+tag+cut+weight,";M_{D^{*}};Events / 10 MeV" ,200,0.,2.0);
    allPlots["pi_pt"+tag+cut+weight] = new TH1F("pi_pt"+tag+cut+weight,";#pi^{#pm} P_{T} [GeV];Events / 5 GeV", 10, 0,50);
    allPlots["MET"+tag+cut+weight] = new TH1F("MET"+tag+cut+weight,";MET [GeV];Events / 20 GeV", 10,0,200);
    allPlots["charge"+tag+cut+weight] = new TH1F("charge"+tag+cut+weight,";Charge(l_{1}*l_{2});Events", 5,-2,2);
    allPlots["csv"+tag+cut+weight] = new TH1F("CSV"+tag+cut+weight,";Jet CSV;Events / 0.01", 80,0.8,1);
    allPlots["dR"+tag+cut+weight] = new TH1F("dR"+tag+cut+weight,";dR;Events / 0.05", 20,0.0,1.);
    allPlots["pflp_pt"+tag+cut+weight] = new TH1F("pflp_pt"+tag+cut+weight,";PF lepton P_{T} [GeV];Events / 0.2 GeV", 15, 0,3);
    allPlots["massZ"+tag+cut+weight]     = new TH1F("massZ_control"+tag+cut+weight,";M_{ll};Events / 1.0 GeV" ,30,81,111);
    allPlots["nevt"+tag+cut+weight]     = new TH1F("nevt"+tag+cut+weight,";N_{events};Events" ,1,1.,2.);

  }
  }
  }
    allPlots["relIso_m"] = new TH1F("relIso_m",";relIso;Events / 0.05", 20,0,1.);
    allPlots["relIso_e"] = new TH1F("relIso_e",";relIso;Events / 0.05", 20,0,1.);


  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  //for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }
  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%5000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

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
      //std::vector<int> tightLeptonsNonIso, vetoLeptons;
      std::vector<int> tightLeptons,looseLeptons;
      for(int il=0; il<ev.nl; il++)
	{
          //cout << "in lepton selection" << endl;
	  //bool passTightKin(ev.l_pt[il]>30 && fabs(ev.l_eta[il])<2.1);
	  bool passTightKin(ev.l_pt[il]>20 && fabs(ev.l_eta[il])<2.5);
	  //bool passVetoKin(ev.l_pt[il]>10 && fabs(ev.l_eta[il])<2.5); //FIXME from 7_6_x
	  float relIso(ev.l_relIso[il]);
	  bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);
	  //bool passTightIso( ev.l_id[il]==13 ? relIso<0.15 : (ev.l_pid[il]>>1)&0x1 ); //FIXME from 7_6_x
	  bool passIso( ev.l_id[il]==13 ? relIso<0.15 : (ev.l_pid[il]>>1)&0x1 ); //FIXME from 7_6_x
	  //bool passNonIso(relIso>0.4); //FIXME from 7_6_x
	  //if( ev.l_id[il]==11 && (passIso || relIso<0.4) ) passNonIso=false; //FIXME from 7_6_x

	  //bool passVetoIso(  ev.l_id[il]==13 ? relIso<0.25 : true); //FIXME from 7_6_x

	  //bool passSIP3d(ev.l_ip3dsig[il]<4);
	  //if(channelSelection==21) passSIP3d=true;
	  if( ev.l_id[il] == 13 )
            allPlots["relIso_m"]->Fill(relIso,wgt);
          else if( ev.l_id[il] == 11)
            allPlots["relIso_e"]->Fill(relIso,wgt);

	  if(passTightKin && passTightId)// && passSIP3d)
	    {
	      if(passIso)         tightLeptons.push_back(il);
	      //else if(passNonIso) tightLeptonsNonIso.push_back(il); //FIXME from 7_6_x
	    }
	  //else if(passVetoKin && passVetoIso) vetoLeptons.push_back(il); //FIXME from 7_6_x
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
      if(tightLeptons.size()==1 )//&& looseLeptons.size()==0)
	{
	  selLeptons.push_back( tightLeptons[0] );
          if(debug) cout << "found 1 tight lepton" << endl;
	}
      if(tightLeptons.size()>=2)
	{
	  selLeptons.push_back(tightLeptons[0]);
	  selLeptons.push_back(tightLeptons[1]);
          if(debug) cout << "found 2 tight leptons" << endl;
	}
      if(tightLeptons.size()==1 && looseLeptons.size()>=1)
	{
	  selLeptons.push_back(tightLeptons[0]);
	  selLeptons.push_back(looseLeptons[0]);	  
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
	      if(tightLeptons.size()>=2 && (hasEMTrigger)) chTag="em";
	      if(tightLeptons.size()==1)
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
      if(tightLeptons.size()==1)                                       lepIdx=tightLeptons[0];
      //else if (tightLeptons.size()==0 && tightLeptonsNonIso.size()==1) lepIdx=tightLeptonsNonIso[0];
      else if(tightLeptons.size()==2)
	{	  
          if(debug) cout << "di-lepton" << endl;
	  l1p4.SetPtEtaPhiM(ev.l_pt[tightLeptons[0]],ev.l_eta[tightLeptons[0]],ev.l_phi[tightLeptons[0]],ev.l_mass[tightLeptons[0]]);
	  l2p4.SetPtEtaPhiM(ev.l_pt[tightLeptons[1]],ev.l_eta[tightLeptons[1]],ev.l_phi[tightLeptons[1]],ev.l_mass[tightLeptons[1]]);
	  dilp4=l1p4+l2p4;
	  if(ev.l_id[tightLeptons[0]]==ev.l_id[tightLeptons[1]]          && 
	     ev.l_charge[tightLeptons[0]]*ev.l_charge[tightLeptons[1]]<0 && 
	     fabs(dilp4.M()-91)<10) //&& 
	     //dilp4.Pt()>30)
	    { 
	      isZ=true; 
	      //isZPassingSIP3d=(ev.l_ip3dsig[0]<4 && ev.l_ip3dsig[1]<4);
	    }
	  lepIdx=tightLeptons[0];
          if(debug) cout << "di-lepton DONE" << endl;
	}

      if(lepIdx<0) continue;
      
      //no extra isolated leptons
      //if(vetoLeptons.size()>0) continue;
      
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
      std::vector<Jet> bJetsVec, allJetsVec;
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
          /*
	  float genJet_pt(0);
	  if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
          */
          /*
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }
          */
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
	    tmpj.addTrack(tkP4,ev.pf_id[ipf]);
	  }
	  tmpj.sortTracksByPt();

          if(isBTagged) bJetsVec.push_back(tmpj);
          allJetsVec.push_back(tmpj);
          if(isBTagged) allPlots["csv"+chTag]->Fill(csv,wgt);
          if(isBTagged) allPlots["csv_all"]->Fill(csv,wgt);
	}

      
      //event weight
      //float wgt(1.0);
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
	  triggerCorrWgt=lepEffH.getTriggerCorrection(selLeptons,leptons);
          if(debug) cout << "calling trigger function DONE!" << endl;
	  for(size_t il=0; il<selLeptons.size(); il++) {
	    EffCorrection_t selSF=lepEffH.getOfflineCorrection(selLeptons[il],leptons[il].Pt(),leptons[il].Eta());
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
	  float norm( normH ? normH->GetBinContent(1) : 1.0);
	  //wgt=lepTriggerSF*lepSelSF*puWeight*norm;
	  wgt=triggerCorrWgt.first*lepSelCorrWgt.first*puWgts[0]*norm;
	  if(ev.ttbar_nw>0) wgt*=ev.ttbar_w[0];
          if(debug) cout << "weight=" << wgt << endl;
          if(debug) cout << "Trigger=" << triggerCorrWgt.first << endl << "Lepton=" << lepSelCorrWgt.first << endl << "PU=" << puWgts[0] << endl << "norm=" << norm  << endl;
          //wgt=1.0;
	}
      if(debug) cout << "Lepton scale factors DONE!" << endl;

      //sort by Pt
      sort(bJetsVec.begin(),    bJetsVec.end(),   sortJetsByPt);
      sort(allJetsVec.begin(),  allJetsVec.end(), sortJetsByPt);

      if(bJetsVec.size()==0) continue;
      for(size_t il=0; il<leptons.size(); il++) {
        for(size_t ij=0; ij<bJets.size(); ij++) {
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[ij],ev.j_eta[ij],ev.j_phi[ij],ev.j_mass[ij]);
          allPlots["dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
          allPlots["dR"+chTag+"_no_weight"]->Fill(jp4.DeltaR(leptons[il]),1);
        }
        for(size_t ij=0; ij<lightJets.size(); ij++) {
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[ij],ev.j_eta[ij],ev.j_phi[ij],ev.j_mass[ij]);
          allPlots["dR"+chTag]->Fill(jp4.DeltaR(leptons[il]),wgt);
          allPlots["dR_all"]->Fill(jp4.DeltaR(leptons[il]),wgt);
          allPlots["dR"+chTag+"_no_weight"]->Fill(jp4.DeltaR(leptons[il]),1);
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
      std::sort(lightJets.begin(), lightJets.end(), VecSort);
      std::sort(bJets.begin(), bJets.end(), VecSort);
      if(debug) cout << "starting simple plots" << endl;
      if(selLeptons.size() == 1 && bJets.size() >= 1 && lightJets.size() >= 4) {
        singleLep = true;
        allPlots["nj"+chTag]->Fill(lightJets.size(),wgt);
        allPlots["nbj"+chTag]->Fill(bJets.size(),wgt);
        allPlots["nj"+chTag+"_no_weight"]->Fill(lightJets.size(),1);
        allPlots["nbj"+chTag+"_no_weight"]->Fill(bJets.size(),1);
        allPlots["nj_all"]->Fill(lightJets.size(),wgt);
        allPlots["nbj_all"]->Fill(bJets.size(),wgt);
        allPlots["nlp"+chTag]->Fill(selLeptons.size(),wgt);
        allPlots["nlp"+chTag+"_no_weight"]->Fill(selLeptons.size(),1);
        allPlots["nevt"+chTag]->Fill(1,wgt);
        std::sort(leptons.begin(), leptons.end(), VecSort);
        allPlots["lp_pt"+chTag]->Fill(leptons[0].Pt(),wgt);
        allPlots["lp_pt"+chTag+"_no_weight"]->Fill(leptons[0].Pt(),1);
        allPlots["lp_pt_all"]->Fill(leptons[0].Pt(),wgt);
        if(debug) cout << "sorting jets" << endl;
	std::sort(lightJets.begin(), lightJets.end(), VecSort);
        std::sort(bJets.begin(), bJets.end(), VecSort);
        /*
        for (auto it : lightJets)
          allPlots["j_pt"+chTag]->Fill(it.Pt(),wgt);
        for (auto it : bJets)
          allPlots["bj_pt"+chTag]->Fill(it.Pt(),wgt);
        */
        allPlots["j_pt"+chTag]->Fill(lightJets[0].Pt(),wgt);
        allPlots["bj_pt"+chTag]->Fill(bJets[0].Pt(),wgt);
        allPlots["j_pt_all"]->Fill(lightJets[0].Pt(),wgt);
        allPlots["bj_pt_all"]->Fill(bJets[0].Pt(),wgt);
        allPlots["MET"+chTag]->Fill(ev.met_pt[0],wgt);
        allPlots["j_pt"+chTag+"_no_weight"]->Fill(lightJets[0].Pt(),1);
        allPlots["bj_pt"+chTag+"_no_weight"]->Fill(bJets[0].Pt(),1);
        allPlots["MET"+chTag+"_no_weight"]->Fill(ev.met_pt[0],1);
        for(int i = 0; i < ev.npf; i++) {
          allPlots["pfid"+chTag]->Fill(ev.pf_id[i],wgt);
        }
      }
      else if(selLeptons.size() == 2 && bJets.size() >= 1 && lightJets.size() > 1) {
        if(isZ) {
          allPlots["massZ"+chTag]->Fill(dilp4.M(),wgt);
          allPlots["massZ"+chTag+"_no_weight"]->Fill(dilp4.M(),1);
        }
        if(isZ) continue;
        if(dilp4.M() < 10) continue; // && ev.l_charge[selLeptons[0]]!=ev.l_charge[selLeptons[1]]) continue;
        if(ev.l_id[selLeptons[0]]==ev.l_id[selLeptons[1]] && met.Pt() < 40) continue; //FIXME
        doubleLep = true;
        allPlots["nj"+chTag]->Fill(lightJets.size(),wgt);
        allPlots["nbj"+chTag]->Fill(bJets.size(),wgt);
        allPlots["ndilp"+chTag]->Fill(selLeptons.size(),wgt);
        allPlots["dilp_pt"+chTag]->Fill(dilp4.Pt(),wgt);
        allPlots["dilp_m"+chTag]->Fill(dilp4.M(),wgt);
        allPlots["nj"+chTag+"_no_weight"]->Fill(lightJets.size(),1);
        allPlots["nbj"+chTag+"_no_weight"]->Fill(bJets.size(),1);
        allPlots["ndilp"+chTag+"_no_weight"]->Fill(selLeptons.size(),1);
        allPlots["dilp_pt"+chTag+"_no_weight"]->Fill(dilp4.Pt(),1);
        allPlots["dilp_m"+chTag+"_no_weight"]->Fill(dilp4.M(),1);
        std::sort(leptons.begin(), leptons.end(), VecSort);
        allPlots["lp_pt"+chTag]->Fill(leptons[0].Pt(),wgt);
        allPlots["l2p_pt"+chTag]->Fill(leptons[1].Pt(),wgt);
        allPlots["lp_pt"+chTag+"_no_weight"]->Fill(leptons[0].Pt(),1);
        allPlots["l2p_pt"+chTag+"_no_weight"]->Fill(leptons[1].Pt(),1);
        if(debug) cout << "sorting jets" << endl;
	std::sort(lightJets.begin(), lightJets.end(), VecSort);
        std::sort(bJets.begin(), bJets.end(), VecSort);
        /*
        for (auto it : lightJets)
          allPlots["j_pt"+chTag]->Fill(it.Pt(),wgt);
        for (auto it : bJets)
          allPlots["bj_pt"+chTag]->Fill(it.Pt(),wgt);
        */
        allPlots["j_pt"+chTag]->Fill(lightJets[0].Pt(),wgt);
        allPlots["bj_pt"+chTag]->Fill(bJets[0].Pt(),wgt);
        allPlots["j_pt"+chTag+"_no_weight"]->Fill(lightJets[0].Pt(),1);
        allPlots["bj_pt"+chTag+"_no_weight"]->Fill(bJets[0].Pt(),1);
        TLorentzVector leadlp41,leadlp42;
        /*
        leadlp41.SetPtEtaPhiM(ev.l_pt[selLeptons[0]],ev.l_eta[selLeptons[0]],ev.l_phi[selLeptons[0]],ev.l_mass[selLeptons[0]]);
        leadlp42.SetPtEtaPhiM(ev.l_pt[selLeptons[1]],ev.l_eta[selLeptons[1]],ev.l_phi[selLeptons[1]],ev.l_mass[selLeptons[1]]);
        allPlots["leadL_pt"+chTag]->Fill(ev.l_pt[selLeptons[0]],wgt);
        allPlots["leadL_pt"+chTag]->Fill(ev.l_pt[selLeptons[1]],wgt);
        */
        allPlots["MET"+chTag]->Fill(ev.met_pt[0],wgt);
        allPlots["charge"+chTag]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],wgt);
        allPlots["MET"+chTag+"_no_weight"]->Fill(ev.met_pt[0],1);
        allPlots["charge"+chTag+"_no_weight"]->Fill(ev.l_charge[selLeptons[0]]*ev.l_charge[selLeptons[1]],1);
        for(int i = 0; i < ev.npf; i++) {
          allPlots["pfid"+chTag]->Fill(ev.pf_id[i],wgt);
        }
      }
      if(debug) cout << "simple plots DONE" << endl;


      if(!singleLep && !doubleLep) continue;

      allPlots["npf"+chTag+"_lep"]->Fill(ev.npf,wgt);
      allPlots["npf"+chTag+"_lep"+"_no_weight"]->Fill(ev.npf,1);
      allPlots["nevt"+chTag+"_lep"]->Fill(1,wgt);
      allPlots["nevt_all_lep"]->Fill(1,wgt);
      //float lpt(0), bpt(0);
      //int pfid = 0;
      /*
      //Track candidates
      for(int it = 0; it < ev.npf; it++) {
        if(ev.pf_pt[it] > pfpt) pfpt = ev.pf_pt[it];
        if(abs(ev.pf_id[it]) == 13 || abs(ev.pf_id[it]) == 11) {
          allPlots["nlp"+chTag+"_lep"]->Fill(1,wgt);
          allPlots["lp_pt"+chTag+"_lep"]->Fill(ev.pf_pt[it],wgt);
        }
        else if(abs(ev.pf_id[it]) == 211) {
          allPlots["nbj"+chTag+"_lep"]->Fill(1,wgt);
          allPlots["bj_pt"+chTag+"_lep"]->Fill(ev.pf_pt[it],wgt);
        }
      }
      */

      //charmed resonance analysis : use only jets with CSV>CSVL, up to two per event
      //sort(allJetsVec.begin(),    allJetsVec.end(),     sortJetsByCSV);
      for(size_t ij = 0; ij < bJetsVec.size(); ij++) {

        if(ij > 1) continue;
        if(bJetsVec[ij].getCSV()<0.460) continue;

        std::vector<IdTrack> &tracks = bJetsVec[ij].getTracks();

        //J/Psi
        if(debug) cout << "starting J/Psi" << endl;
        float gMassMu(0.1057),gMassK(0.4937);//,gMassPi(0.1396);
        std::vector<TLorentzVector> pfmuCands,kaonCands;
        for(size_t itk = 0; itk < tracks.size(); itk++) {
          if(abs(tracks[itk].second) == 13) {
            TLorentzVector muP4;
            muP4.SetPtEtaPhiM(tracks[itk].first.Pt(), tracks[itk].first.Eta(), tracks[itk].first.Phi(), gMassMu);
            pfmuCands.push_back(muP4);
          }
          if(abs(tracks[itk].second) == 211) {
            TLorentzVector kP4;
            kP4.SetPtEtaPhiM( tracks[itk].first.Pt(), tracks[itk].first.Eta(), tracks[itk].first.Phi(), gMassK);
            kaonCands.push_back(kP4);
          }
    
          if(pfmuCands.size()>1) {
            float mass12((pfmuCands[0] + pfmuCands[1]).M());
            float mass123( kaonCands.size()>0 ? (pfmuCands[0]+pfmuCands[1]+kaonCands[0]).M() : -1);
            allPlots["massJPsi"+chTag]->Fill(mass12,wgt);
	    allPlots["massJPsi_all"]->Fill(mass12,wgt);
            if(mass123 > 0) {
              allPlots["massJPsiK"+chTag]->Fill(mass123,wgt);
              allPlots["massJPsiK_all"]->Fill(mass123,wgt);
            }
          }
        }
        if(debug) cout << "J/Psi DONE" << endl;
        //continue; //FIXME

        //D0 and D* 
        if(debug) cout << "Starting D0 and D*" << endl;
        //nstart = firstTrackIndex(jetindex,&tracks);
        //if((tracks.size() - nstart) < 3) continue;
        int jetindex = allJetsVec[ij].getJetIndex();
        if(tracks.size() < 3) continue;
        for(size_t i = 0; i < 3; i++)
          //for(size_t j = i+1; j < i+2; j++)
          for(size_t j = 0; j < 3; j++) {
            if(i == j) continue;
            /*
            int tk1 = get<2>(tracks.at(i));
            int tk2 = get<2>(tracks.at(j));
            */

            /*
            allPlots["npf"+chTag+"_meson"]->Fill(ev.npf,wgt);
            allPlots["npf"+chTag+"_meson_no_weight"]->Fill(ev.npf,1);
            for(unsigned int it = 0; it < jetPt.size(); it++)
              bpt = get<2>(jetPt.at(it)) > bpt ? get<2>(jetPt.at(it)) : bpt;
            allPlots["bj_pt"+chTag+"_meson"]->Fill(bpt,wgt);
            allPlots["bj_pt"+chTag+"_meson_no_weight"]->Fill(bpt,1);
            allPlots["pfid"+chTag+"_meson"]->Fill(ev.pf_id[tk1],wgt);
            allPlots["nbj"+chTag+"_meson"]->Fill(1,wgt);
            allPlots["nbj"+chTag+"_meson_no_weight"]->Fill(1,1);
            */

            //opposite sign
            //if(ev.pf_id[tk1]*ev.pf_id[tk2] != -211*211) continue;
            if(tracks[i].second*tracks[j].second != -211*211) continue;

            const float gMassK  = 0.4937;
            const float gMassPi = 0.1396;
          
            //p_track1.SetPtEtaPhiM(ev.pf_pt[tk1], ev.pf_eta[tk1], ev.pf_phi[tk1], gMassPi);
            //p_track2.SetPtEtaPhiM(ev.pf_pt[tk2], ev.pf_eta[tk2], ev.pf_phi[tk2], gMassK);
            TLorentzVector p_track1, p_track2;
            p_track1.SetPtEtaPhiM(tracks[i].first.Pt(), tracks[i].first.Eta(), tracks[i].first.Phi(), gMassPi);
            p_track2.SetPtEtaPhiM(tracks[j].first.Pt(), tracks[j].first.Eta(), tracks[j].first.Phi(), gMassK);
            if(debug) cout << tracks[i].first.Pt() << " " << tracks[i].first.Eta() << " " << tracks[i].first.Phi() << " " << gMassPi << endl;
            if(debug) cout << tracks[j].first.Pt() << " " << tracks[j].first.Eta() << " " << tracks[j].first.Phi() << " " << gMassK << endl << endl;
            float mass12 = (p_track1+p_track2).M();
            if(debug) cout << mass12 << endl;
            allPlots["dR"+chTag+"_meson"]->Fill(p_track1.DeltaR(p_track2), wgt);
            allPlots["dR"+chTag+"_meson_no_weight"]->Fill(p_track1.DeltaR(p_track2), 1);

            if (mass12>1.65 && mass12<2.0) {
              allPlots["massD0"+chTag]->Fill(mass12,wgt);
              allPlots["massD0_all"]->Fill(mass12,wgt);
              allPlots["massD0"+chTag+"_no_weight"]->Fill(mass12,1);
            }

            //looking for lepton
            if(debug) cout << "third lepton" << endl;
            //for(int tk3 = 0; tk3 < ev.npf; tk3++)
            for(size_t k = 0; k < tracks.size(); k++) {
              //int tk3 = get<2>(tracks.at(k));
              //if(ev.pf_j[tk3] != jetindex) continue;
              //if(tk3 == tk1) continue;
              //if(tk3 == tk2) continue;
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
                allPlots["massD0_lep"+chTag+"_no_weight"]->Fill(mass12,1);

                if(abs(tracks[k].second) == 13)
                  allPlots["massD0_mu"+chTag]->Fill(mass12,wgt);
                if(abs(tracks[k].second) == 11)
                  allPlots["massD0_e"+chTag]->Fill(mass12,wgt);
                if(abs(tracks[k].second) == 13)
                  allPlots["massD0_mu"+chTag+"_no_weight"]->Fill(mass12,1);
                if(abs(tracks[k].second) == 11)
                  allPlots["massD0_e"+chTag+"_no_weight"]->Fill(mass12,1);

              }
            }
            //looking for pion
            if(debug) cout << "D*->pi+D0" << endl;
            //for(int tk3 = 0; tk3 < ev.npf; tk3++)
            for(size_t k = 0; k < tracks.size(); k++) {
              //int tk3 = get<2>(tracks.at(k));
              //if(ev.pf_j[tk3] != jetindex) continue;
              //if(tk3 == tk1) continue;
              //if(tk3 == tk2) continue;
              if(k == i) continue;
              if(k == j) continue;

              if(abs(tracks[k].second) != 211) continue;
              if(debug) cout << "Pion found" << endl;

              TLorentzVector p_track3, p_cand;
              p_track3.SetPtEtaPhiM(tracks[k].first.Pt(), tracks[k].first.Eta(), tracks[k].first.Phi(), gMassPi);
              allPlots["pi_pt"+chTag]->Fill(p_track3.Pt(),wgt);
              allPlots["pi_pt"+chTag+"_no_weight"]->Fill(p_track3.Pt(),1);
              if( tracks[j].second/abs(tracks[j].second) == -tracks[k].second/abs(tracks[k].second) ) {
                // Kaon and pion have opposite charges
                // I.e. correct mass assumption
                if(debug) cout << "correct mass assumption" << endl;
                
                p_cand = p_track1+p_track2+p_track3;
                allPlots["massDs"+chTag]->Fill(p_cand.M(), wgt);
                allPlots["massDs"+chTag+"_no_weight"]->Fill(p_cand.M(), 1);
                allPlots["massDs_all"]->Fill(p_cand.M(), wgt);

                if(abs(mass12-1.864) < 0.10) { // mass window cut
                  TLorentzVector p_jet;
                  p_jet.SetPtEtaPhiM(ev.j_pt[jetindex], ev.j_eta[jetindex], ev.j_phi[jetindex], 0.);

                  //float hardpt = std::max(ev.pf_pt[tk3], std::max(ev.pf_pt[tk1], ev.pf_pt[tk2]));
                  //float softpt = std::min(ev.pf_pt[tk3], std::min(ev.pf_pt[tk1], ev.pf_pt[tk2]));
                  float deltam = p_cand.M() - mass12;

                  if(filename.Contains("_WJets"))
                    cout << deltam << " " << wgt << endl;
                  allPlots["massDsmD0loose"+chTag]->Fill(deltam, wgt);
                  allPlots["massDsmD0loose"+chTag+"_no_weight"]->Fill(deltam, 1);
                  allPlots["massDsmD0loose_all"]->Fill(deltam, wgt);
                  if(abs(mass12-1.864) < 0.05) { // tighter mass window cut
                      //FillCharmTree(413,  jetindex, tk1, gMassPi, tk2, gMassK, tk3, gMassPi);
                      //FillCharmTree(-413, jetindex, deltam, p_cand, p_jet, hardpt, softpt);
                      allPlots["massDsmD0"+chTag]->Fill(deltam, wgt);
                      allPlots["massDsmD0"+chTag+"_no_weight"]->Fill(deltam, 1);
                      allPlots["massDsmD0_all"]->Fill(deltam, wgt);
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
/*
  for(auto it : lfsVec) {
    it.Remove(0,1);
    fOut->mkdir(it);
  }
*/
  for (auto& it : allPlots)  { 
    if(debug) cout << it.second->GetName() << endl;
    if(debug) cout << it.second->GetEntries() << endl;
/*
    TString dir = it.first;
    dir.Remove(0,dir.Last('_')+1);
*/
    //fOut->cd( dir );
    it.second->SetDirectory(fOut); it.second->Write(); 
    fOut->cd();
  }
  if(debug) cout << "writing histograms DONE" << endl;
  if(debug) cout << "closing ROOT file" << endl;
  fOut->Close();
}

