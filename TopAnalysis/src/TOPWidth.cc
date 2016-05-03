#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"


using namespace std;

//
void RunTOPWidth(TString filename,
              TString outname,
              Int_t channelSelection, 
              Int_t chargeSelection, 
              FlavourSplitting flavourSplitting,
              TH1F *normH,
              Bool_t runSysts = false,
              Bool_t interpolate = false,
              int nIntrps = 0,
              float minIntrpWid = 1,
              float maxIntrpWid = 2) 
{

  bool isTTbar( filename.Contains("_TTJets") ); // TODO: use?

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);
  if(ev.isData) runSysts=false;
  cout << "...producing " << outname << " from " << nentries << " events" << (runSysts ? " syst variations will be considered" : "") << endl;

  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;

  //for data only get the lumi per run map
  std::map<Int_t,Float_t> lumiMap;
  if(ev.isData) lumiMap=lumiPerRun();

  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  if(!ev.isData)
    {
      TString puWgtUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/pileupWgts.root");
      gSystem->ExpandPathName(puWgtUrl);
      TFile *fIn=TFile::Open(puWgtUrl);
      if(fIn)
	{
	  puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_nom") );
	  puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_down") );
	  puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_up") );
	  fIn->Close();
	}
    }

  //LEPTON EFFICIENCIES
  TString lepEffUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/leptonEfficiencies.root");
  gSystem->ExpandPathName(lepEffUrl);
  std::map<TString,TH2 *> lepEffH;
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH["m_sel"]=(TH2 *)fIn->Get("m_sel");
      lepEffH["m_trig"]=(TH2 *)fIn->Get("m_trig");      
      for(auto& it : lepEffH) it.second->SetDirectory(0);
      fIn->Close();
    }

  lepEffUrl="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CutBasedID_TightWP_fromTemplates_withSyst_Final.txt_SF2D.root";
  gSystem->ExpandPathName(lepEffUrl);
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH["e_sel"]=(TH2 *)fIn->Get("EGamma_SF2D");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
      fIn->Close();
    }

  //B-TAG CALIBRATION
  TString btagUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CSVv2.csv");
  gSystem->ExpandPathName(btagUncUrl);
  std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
  TString btagEffExpUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/expTageff.root");
  gSystem->ExpandPathName(btagEffExpUrl);
  std::map<TString, TGraphAsymmErrors *> expBtagEff;
  BTagSFUtil myBTagSFUtil;
  if(!ev.isData)
    {
      BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "central") );
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "down") ); 
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "up") );

      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "comb", "central") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "comb", "down") ); 
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "comb", "up") );
      
      TFile *beffIn=TFile::Open(btagEffExpUrl);
      expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();
    }

  //JET ENERGY SCALE: https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Summer15_uncertainties
  TString jecUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/Summer15_25nsV7_DATA_UncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  //FactorizedJetCorrector *jetCorr=getFactorizedJetEnergyCorrector("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/jectxt",!ev.isData);
  std::vector<TString> jecUncSrcs;
  std::vector<JetCorrectionUncertainty*> jecUncs;
  if(runSysts)
    {
      jecUncSrcs.push_back("AbsoluteStat");        
      jecUncSrcs.push_back("AbsoluteScale");        
      jecUncSrcs.push_back("AbsoluteFlavMap");        
      jecUncSrcs.push_back("AbsoluteMPFBias");        
      jecUncSrcs.push_back("Fragmentation");        
      jecUncSrcs.push_back("SinglePionECAL");  
      jecUncSrcs.push_back("SinglePionHCAL"); 
      jecUncSrcs.push_back("TimeEta");
      jecUncSrcs.push_back("TimePt");
      jecUncSrcs.push_back("RelativeJEREC1");  
      jecUncSrcs.push_back("RelativeJEREC2"); 
      jecUncSrcs.push_back("RelativeJERHF");
      jecUncSrcs.push_back("RelativePtBB");    
      jecUncSrcs.push_back("RelativePtEC1");  
      jecUncSrcs.push_back("RelativePtEC2");   
      jecUncSrcs.push_back("RelativePtHF");  
      jecUncSrcs.push_back("RelativeFSR");
      jecUncSrcs.push_back("RelativeStatEC"); 
      jecUncSrcs.push_back("RelativeStatHF");
      jecUncSrcs.push_back("PileUpDataMC");    
      jecUncSrcs.push_back("PileUpPtRef");    
      jecUncSrcs.push_back("PileUpPtBB");     
      jecUncSrcs.push_back("PileUpPtEC1");      
      jecUncSrcs.push_back("PileUpPtEC2");      
      jecUncSrcs.push_back("PileUpPtHF");    
      jecUncSrcs.push_back("FlavorPureGluon"); 
      jecUncSrcs.push_back("FlavorPureQuark");
      jecUncSrcs.push_back("FlavorPureCharm"); 
      jecUncSrcs.push_back("FlavorPureBottom");
      for(size_t i=0; i<jecUncSrcs.size(); i++)
	{
	  JetCorrectorParameters *p = new JetCorrectorParameters(jecUncUrl.Data(), jecUncSrcs[i].Data());
	  jecUncs.push_back( new JetCorrectionUncertainty(*p) );
	}
    }

  //LIST OF SYSTEMATICS
  Int_t nGenSysts(0);
  std::vector<TString> expSysts;
  if(runSysts)
    {
      if(normH)
	for(Int_t xbin=1; xbin<=normH->GetNbinsX(); xbin++) 
	  nGenSysts += (normH->GetBinContent(xbin)!=0);	
      
      expSysts=jecUncSrcs;
      expSysts.push_back("JER");
      expSysts.push_back("Pileup");
      expSysts.push_back("MuTrigger");
      expSysts.push_back("MuEfficiency");
      expSysts.push_back("MuScale");
      expSysts.push_back("EleTrigger");
      expSysts.push_back("EleEfficiency");
      expSysts.push_back("EleScale");
      expSysts.push_back("BtagEff");
      expSysts.push_back("CtagEff");
      expSysts.push_back("LtagEff");
      expSysts.push_back("UncMET");
      expSysts.push_back("topPt");
      
      cout << "\t..." << nGenSysts << "/" << expSysts.size() << " generator level/experimental systematics will be considered" << endl;
    }


  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  std::map<TString, TH2 *> all2dPlots;
  std::vector<TString> lfsVec = { "E", "EE", "EM", "MM", "M" }; // eac hardcoded
  //allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",2,0,2);
//  for(Int_t ij=1; ij<=4; ij++)        {
//  for(Int_t itag=-1; itag<=2; itag++) {
  for(size_t ilfs=0; ilfs<lfsVec.size(); ilfs++)   { 
//	  if(itag>ij) continue;
	  TString tag(/*itag<0 ? Form("%s",lfsVec.at(ilfs).Data()) :*/ Form("%s",lfsVec.at(ilfs).Data()));
          
	  if(lumiMap.size()) allPlots["mlbwa_"+tag+"_ratevsrun"] = new TH1F("mlbwa_"+tag+"_ratevsrun",";Run number; Events/pb",lumiMap.size(),0,lumiMap.size());
	  //Int_t runCtr(0);
	  //for(std::map<Int_t,Float_t>::iterator it=lumiMap.begin(); it!=lumiMap.end(); it++,runCtr++)
	    //allPlots["mlbwa_"+tag+"_ratevsrun"]->GetXaxis()->SetBinLabel(runCtr+1,Form("%d",it->first));
	  //allPlots["mlbwa_"+tag+"_L_Pt"]     = new TH1F("mlbwa_"+tag+"_L_Pt",";Transverse momentum [GeV];Events" ,20,0.,300.);
	  //allPlots["mlbwa_"+tag+"_L_Eta"]    = new TH1F("mlbwa_"+tag+"_L_Eta",";Pseudo-rapidity;Events" ,12,0.,3.);
	  //allPlots["mlbwa_"+tag+"_J_Pt"]     = new TH1F("mlbwa_"+tag+"_J_Pt",";Transverse momentum [GeV];Events" ,20,0.,300.);
	  //allPlots["mlbwa_"+tag+"_J_Eta"]    = new TH1F("mlbwa_"+tag+"_J_Eta",";Pseudo-rapidity;Events" ,12,0.,3.);
	  //allPlots["mlbwa_"+tag+"_mttbar"]   = new TH1F("mlbwa_"+tag+"_mttbar",";#sqrt{#hat{s}} [GeV];Events" ,50,0.,1000.);
	  //allPlots["mlbwa_"+tag+"_TMass"]    = new TH1F("mlbwa_"+tag+"_TMass",";Transverse Mass [GeV];Events" ,20,0.,200.);
	  //allPlots["mlbwa_"+tag+"_mTop"]     = new TH1F("mlbwa_"+tag+"_mTop",";Top quark mass [GeV];Events" ,50,0.,500.);
	  allPlots["mlbwa_"+tag+"_Count"]    = new TH1F("mlbwa_"+tag+"_Count","Event present;Events" ,25,0.,250.);
	  allPlots["mlbwa_"+tag+"_MET"]      = new TH1F("mlbwa_"+tag+"_MET",";Missing transverse energy [GeV];Events" ,10,0.,200.);
	  if(/*itag==-1*/ true)
	    {
	      allPlots["mlbwa_"+tag+"_B_Num"] = new TH1F("mlbwa_"+tag+"_B_Num",";Category;Events" ,3,0.,3.);
	      allPlots["mlbwa_"+tag+"_B_Num"]->GetXaxis()->SetBinLabel(1, "=0b");
	      allPlots["mlbwa_"+tag+"_B_Num"]->GetXaxis()->SetBinLabel(2, "=1b");
	      allPlots["mlbwa_"+tag+"_B_Num"]->GetXaxis()->SetBinLabel(3, "#geq2b");
	    }
	  
      // mlb plots
      allPlots["mlbwa_"+tag+"_Mlb_min"]  = new TH1F("mlbwa_"+tag+"_Mlb_min",";min Mass(lepton,b) [GeV];Events",70,0,350);
	  allPlots["mlbwa_"+tag+"_Mlb_inc"]  = new TH1F("mlbwa_"+tag+"_Mlb_inc",";inc Mass(lepton,b) [GeV];Events",70,0,350);
	  allPlots["mlbwa_"+tag+"_Mlb_mdr"]  = new TH1F("mlbwa_"+tag+"_Mlb_mdr",";mdr Mass(lepton,b) [GeV];Events",70,0,350);

      allPlots["mlbwa_"+tag+"_lbmatch_min"]  = new TH1F("mlbwa_"+tag+"_lbmatch_min",";Category;Events",3,0,3);
      allPlots["mlbwa_"+tag+"_lbmatch_min"]->GetXaxis()->SetBinLabel(1, "correct");
      allPlots["mlbwa_"+tag+"_lbmatch_min"]->GetXaxis()->SetBinLabel(2, "incorrect");
      allPlots["mlbwa_"+tag+"_lbmatch_min"]->GetXaxis()->SetBinLabel(3, "mismatch");

      allPlots["mlbwa_"+tag+"_lbmatch_mdr"]  = new TH1F("mlbwa_"+tag+"_lbmatch_mdr",";Category;Events",3,0,3);
      allPlots["mlbwa_"+tag+"_lbmatch_mdr"]->GetXaxis()->SetBinLabel(1, "correct");
      allPlots["mlbwa_"+tag+"_lbmatch_mdr"]->GetXaxis()->SetBinLabel(2, "incorrect");
      allPlots["mlbwa_"+tag+"_lbmatch_mdr"]->GetXaxis()->SetBinLabel(3, "mismatch");

      allPlots["mlbwa_"+tag+"_lbmatch_inc"]  = new TH1F("mlbwa_"+tag+"_lbmatch_inc",";Category;Events",3,0,3);
      allPlots["mlbwa_"+tag+"_lbmatch_inc"]->GetXaxis()->SetBinLabel(1, "correct");
      allPlots["mlbwa_"+tag+"_lbmatch_inc"]->GetXaxis()->SetBinLabel(2, "incorrect");
      allPlots["mlbwa_"+tag+"_lbmatch_inc"]->GetXaxis()->SetBinLabel(3, "mismatch");

      all2dPlots["mlbwa_"+tag+"_Tcmp_min"] = new TH2F("mlbwa_"+tag+"_Tcmp_min", 
                                                      "top Mass(lepton,b) [GeV];min mass [GeV]",
                                                      70, 0, 350, 75, 0, 375);
      all2dPlots["mlbwa_"+tag+"_Tcmp_mdr"] = new TH2F("mlbwa_"+tag+"_Tcmp_mdr", 
                                                      "top Mass(lepton,b) [GeV];mdr mass [GeV]",
                                                      70, 0, 350, 75, 0, 375);
      all2dPlots["mlbwa_"+tag+"_Tcmp_inc"] = new TH2F("mlbwa_"+tag+"_Tcmp_inc", 
                                                      "top Mass(lepton,b) [GeV];inc mass [GeV]",
                                                      70, 0, 350, 75, 0, 375);

	  //shape analysis
	  if(runSysts)
	    {
	      //gen level systematics
	      if(normH)
		{
		  all2dPlots["metptshapes_"+tag+"_gen"]                   
		    = new TH2F("metptshapes_"+tag+"_gen", ";Missing transverse energy [GeV];Events" ,    10,0.,200., nGenSysts,0,nGenSysts);
		  all2dPlots["mlbshapes_"+tag+"_gen"]               
		    = new TH2F("mlbshapes_"+tag+"_gen", ";min Mass(lepton,b) [GeV];Events" , 25,0.,250., nGenSysts,0,nGenSysts);
		  all2dPlots["RMPFshapes_"+tag+"_gen"]               
		    = new TH2F("RMPFshapes_"+tag+"_gen", ";R_{MPF};Events" , 20,0.,2., nGenSysts,0,nGenSysts);
		  if(/*itag==-1*/ true) 
		    all2dPlots["nbtagsshapes_"+tag+"_gen"] 
		      = new TH2F("nbtagsshapes_"+tag+"_gen", ";Category;Events" , 3, 0.,3.,   nGenSysts,0,nGenSysts);		  
		  for(Int_t igen=0; igen<nGenSysts; igen++)
		    {
		      TString label(normH->GetXaxis()->GetBinLabel(igen+1));
		      all2dPlots["metptshapes_"+tag+"_gen"]    ->GetYaxis()->SetBinLabel(igen+1,label);
		      all2dPlots["mlbshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,label);
		      all2dPlots["RMPFshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,label);
		      if(/*itag!=-1*/ false) continue;
		      all2dPlots["nbtagsshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,label);
		    }
		}
	      
	      //experimental systematics
	      Int_t nExpSysts=expSysts.size();
	      if(nExpSysts>0)
		{
		  all2dPlots["metptshapes_"+tag+"_exp"]                  
		    = new TH2F("metptshapes_"+tag+"_exp",  ";Missing transverse energy [GeV];Events" ,   10,0.,200., 2*nExpSysts,0,2*nExpSysts);
		  all2dPlots["mlbshapes_"+tag+"_exp"]              
		    = new TH2F("mlbshapes_"+tag+"_exp", ";min Mass(lepton,b) [GeV];Events" ,25,0.,250., 2*nExpSysts,0,2*nExpSysts);
		  all2dPlots["RMPFshapes_"+tag+"_exp"]               
		    = new TH2F("RMPFshapes_"+tag+"_exp", ";R_{MPF};Events" , 20,0.,2., 2*nExpSysts,0,2*nExpSysts);
		  if(/*itag==-1*/ true) 
		    all2dPlots["nbtagsshapes_"+tag+"_exp"] 
		      = new TH2F("nbtagsshapes_"+tag+"_exp", ";Category;Events" ,  3,0.,3.,    2*nExpSysts,0,2*nExpSysts);
		  for(Int_t isyst=0; isyst<nExpSysts; isyst++)
		    {
		      for(Int_t ivar=0; ivar<2; ivar++)
			{
			  TString label(expSysts[isyst] + (ivar==0 ? "Down" : "Up"));
			  all2dPlots["metptshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
			  all2dPlots["mlbshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
			  all2dPlots["RMPFshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
			  if(/*itag!=-1*/ false) continue;
			  all2dPlots["nbtagsshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
			}
		    }
		}
	    }
    }
    //}
    //}

  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%5000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      //select 1 good lepton
      std::vector<int> tightLeptonsIso, tightLeptonsNonIso, vetoLeptons;
      for(int il=0; il<ev.nl; il++)
	{
	  bool passTightKin(ev.l_pt[il]>30 && fabs(ev.l_eta[il])<2.1);
	  bool passVetoKin(ev.l_pt[il]>10 && fabs(ev.l_eta[il])<2.5);
	  bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);
	  float relIso(ev.l_relIso[il]);
	  bool passIso( ev.l_id[il]==13 ? relIso<0.15 : (ev.l_pid[il]>>1)&0x1 );
	  bool passNonIso(relIso>0.4);
	  if( ev.l_id[il]==11 && (passIso || relIso<0.4) ) passNonIso=false;
	  bool passVetoIso(  ev.l_id[il]==13 ? relIso<0.25 : true);
	  bool passSIP3d(ev.l_ip3dsig[il]<4);
	  if(channelSelection==21) passSIP3d=true;

	  if(passTightKin && passTightId && passSIP3d)
	    {
	      if(passIso)         tightLeptonsIso.push_back(il);
	      else if(passNonIso) tightLeptonsNonIso.push_back(il);
	    }
	  else if(passVetoKin && passVetoIso) vetoLeptons.push_back(il);
	}


      //one good lepton either isolated or in the non-isolated sideband or a Z candidate
      Int_t lepIdx=-1;
      //Bool_t isZ(false); //isZPassingSIP3d(false);
      TLorentzVector l1p4,l2p4,dilp4;
      TString lfs;
      if(tightLeptonsIso.size()==1) {
          lepIdx=tightLeptonsIso[0];

          // E/M selection
          if(ev.l_id[lepIdx] == 11)       lfs = lfsVec.at(0); // eac hardcoded
          else if (ev.l_id[lepIdx] == 13) lfs = lfsVec.at(4); // eac hardcoded
          else continue;
      } else if (tightLeptonsIso.size()==0 && tightLeptonsNonIso.size()==1) {
          lepIdx=tightLeptonsNonIso[0];

          // E/M selection
          if(ev.l_id[lepIdx] == 11)       lfs = lfsVec.at(0); // eac hardcoded
          else if (ev.l_id[lepIdx] == 13) lfs = lfsVec.at(4); // eac hardcoded
          else continue;
      } else if(tightLeptonsIso.size()==2) {
	    l1p4.SetPtEtaPhiM(ev.l_pt[tightLeptonsIso[0]],
                          ev.l_eta[tightLeptonsIso[0]],
                          ev.l_phi[tightLeptonsIso[0]],
                          ev.l_mass[tightLeptonsIso[0]]);
	    l2p4.SetPtEtaPhiM(ev.l_pt[tightLeptonsIso[1]],
                          ev.l_eta[tightLeptonsIso[1]],
                          ev.l_phi[tightLeptonsIso[1]],
                          ev.l_mass[tightLeptonsIso[1]]);
	    dilp4=l1p4+l2p4;
    //	if(ev.l_id[tightLeptonsIso[0]]==ev.l_id[tightLeptonsIso[1]]          && 
    //	   ev.l_charge[tightLeptonsIso[0]]*ev.l_charge[tightLeptonsIso[1]]<0 && 
    //	   fabs(dilp4.M()-91)<10 && 
    //	   dilp4.Pt()>30)
    //	  {
	//        isZ=true; 
	//        //isZPassingSIP3d=(ev.l_ip3dsig[0]<4 && ev.l_ip3dsig[1]<4);
    //	  }
	    lepIdx=tightLeptonsIso[0];
        
        // EE, EM, MM selection
        if(ev.l_id[lepIdx] == 11 && ev.l_id[tightLeptonsIso[1]] == 11)
            lfs = lfsVec.at(1); // eac hardcoded
        else if(ev.l_id[lepIdx] == 11 && ev.l_id[tightLeptonsIso[1]] == 13)
            lfs = lfsVec.at(2); // eac hardcoded
        else if(ev.l_id[lepIdx] == 13 && ev.l_id[tightLeptonsIso[1]] == 11)
            lfs = lfsVec.at(2); // eac hardcoded
        else if(ev.l_id[lepIdx] == 13 && ev.l_id[tightLeptonsIso[1]] == 13)
            lfs = lfsVec.at(3); // eac hardcoded
        else continue;
   	  }

      if(lepIdx<0) continue;
      
      //no extra isolated leptons
      if(vetoLeptons.size()>0) continue;
            
      //apply trigger requirement
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
     // Int_t lid=ev.l_id[lepIdx];
     // if(isZ) lid=2100+ev.l_id[lepIdx];
     // else if(tightLeptonsNonIso.size()==1) lid=100*ev.l_id[lepIdx];
      if(channelSelection!=0)
	{
//	  if(channelSelection==21) { if(!isZ) continue; }
//	  else                     { if(lid!=channelSelection) continue; }
      if(channelSelection==13 && lfs != "M") continue;
      if(channelSelection==11 && lfs != "E") continue;
      if(channelSelection==13*11 && lfs != "EM") continue;
      if(channelSelection==13*13 && lfs != "MM") continue;
      if(channelSelection==11*11 && lfs != "EE") continue;
	}
//      if(chargeSelection!=0  && ev.l_charge[lepIdx]!=chargeSelection) continue;

      //lepton kinematics
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);

      //select jets
      Float_t htsum(0);
      TLorentzVector jetDiff(0,0,0,0);
      std::vector<TLorentzVector> bJets,lightJets;
      std::vector<int> bJetIDs;
      std::vector<float> gen_tmass;
      std::vector<int> gen_tids;
      //TLorentzVector visSystem(isZ ? dilp4 : lp4);
      int nbjets(0),ncjets(0),nljets(0),leadingJetIdx(-1);
      std::vector<int> resolvedJetIdx;
      std::vector<TLorentzVector> resolvedJetP4;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
	  //jp4=updateJES(jp4,ev.j_rawsf[k],ev.j_area[k],ev.rho,ev.nvtx,jetCorr);

	  //cross clean with respect to leptons 
	  if(jp4.DeltaR(lp4)<0.4) continue;
	  //if(isZ && jp4.DeltaR(l2p4)<0.4)continue;

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
//	  if(isZ && jp4.DeltaPhi(dilp4)<2.7) continue;

	  // re-inforce kinematics cuts
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  if(leadingJetIdx<0) leadingJetIdx=k;
	  htsum += jp4.Pt();
	  //if(bJets.size()+lightJets.size()<4) visSystem += jp4;

	  //b-tag
	  float csv = ev.j_csv[k];	  
	  bool isBTagged(csv>0.890);
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0);
	      if(abs(ev.j_hadflav[k])==4) 
		{ 
		  ncjets++;
		  expEff    = expBtagEff["c"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  nbjets++;
		  expEff    = expBtagEff["b"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		}
	      else
		{
		  nljets++;
		  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
          jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		}
	      
	      //updated b-tagging decision with the data/MC scale factor
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
	    }

	  //save jet
	  if(isBTagged) {
        bJets.push_back(jp4);
        // TODO: add bJetID array
        bJetIDs.push_back(ev.j_hadflav[k]);
      } else lightJets.push_back(jp4);
	}

      //check if flavour splitting was required
//      if(!ev.isData)
//	{
//	  if(flavourSplitting!=FlavourSplitting::NOFLAVOURSPLITTING)
//	    {
//	      if(flavourSplitting==FlavourSplitting::BSPLITTING)         { if(nbjets==0)    continue; }
//	      else if(flavourSplitting==FlavourSplitting::CSPLITTING)    { if(ncjets==0 || nbjets!=0)    continue; }
//	      else if(flavourSplitting==FlavourSplitting::UDSGSPLITTING) { if(nljets==0 || ncjets!=0 || nbjets!=0) continue; }
//	    }
//	}

      //MET and transverse mass
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
      met+=jetDiff;
      met.SetPz(0.); met.SetE(met.Pt());
//      float mt( computeMT(isZ ? dilp4: lp4,met) );
//
//      //compute neutrino kinematics
//      neutrinoPzComputer.SetMET(met);
//      neutrinoPzComputer.SetLepton(isZ ? dilp4 : lp4);
//      
//      float nupz=neutrinoPzComputer.Calculate();
//      TLorentzVector neutrinoHypP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));
//      visSystem+=neutrinoHypP4;
//
//      //balancing variable and ISR control variable
//      float RMPF(0.0),alpha(0.0);
//      TVector2 metT(met.Px(),met.Py());
//      TVector2 visT(isZ ? dilp4.Px() : lp4.Px(), isZ ? dilp4.Py() : lp4.Py());
//      RMPF=1.0+(metT*visT)/visT.Mod2();
//      for(size_t ij=0; ij<resolvedJetIdx.size(); ij++)
//	{
//	  Int_t k=resolvedJetIdx[ij];
//	  if(k==leadingJetIdx) continue;
//	  if(ev.j_pt[k]<15 || fabs(ev.j_eta[k])>3.0) continue;
//	  alpha=ev.j_pt[k]/(isZ ? dilp4.Pt() : lp4.Pt());
//	}

      //event weight
      float wgt(1.0);
      std::vector<float> puWeight(3,1.0),lepTriggerSF(3,1.0),lepSelSF(3,1.0), topPtWgt(3,1.0);
      if(!ev.isData)
	{
	  //update lepton selection scale factors, if found
	  TString prefix("m");
	  if(lfs == "E" || lfs== "EE") prefix="e"; // QUESTION: is this right?
	  if(lepEffH.find(prefix+"_sel")!=lepEffH.end())// && !isZ)
	    {
	      for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)tightLeptonsIso.size()); il++)
		{
		  Int_t ilIdx=tightLeptonsIso[il];
		  float minEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmin() ), maxEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmax()-0.01 );
		  float etaForEff=TMath::Max(TMath::Min(float(fabs(ev.l_eta[ilIdx])),maxEtaForEff),minEtaForEff);
		  Int_t etaBinForEff=lepEffH[prefix+"_sel"]->GetXaxis()->FindBin(etaForEff);
		  		  
		  float minPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmin() ), maxPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmax()-0.01 );
		  float ptForEff=TMath::Max(TMath::Min(float(ev.l_pt[ilIdx]),maxPtForEff),minPtForEff);
		  Int_t ptBinForEff=lepEffH[prefix+"_sel"]->GetYaxis()->FindBin(ptForEff);
		  		  
		  float selSF(lepEffH[prefix+"_sel"]->GetBinContent(etaBinForEff,ptBinForEff));
		  float selSFUnc(lepEffH[prefix+"_sel"]->GetBinError(etaBinForEff,ptBinForEff));
		  lepSelSF[0]*=selSF;      lepSelSF[1]*=(selSF-selSFUnc);       lepSelSF[2]*=(selSF+selSFUnc);
		  
		  float trigSF(1.0), trigSFUnc(0.03);
		  if(prefix=="m")
		    {
		      trigSF=(lepEffH[prefix+"_trig"]->GetBinContent(etaBinForEff,ptBinForEff));
		      trigSFUnc=(lepEffH[prefix+"_trig"]->GetBinError(etaBinForEff,ptBinForEff));
		    }

		  lepTriggerSF[0]*=trigSF; lepTriggerSF[1]*=(trigSF-trigSFUnc); lepTriggerSF[2]*=(trigSF+trigSFUnc);
		}
	    }

	  Int_t ntops(0);
	  float ptsf(1.0);
      // TODO: create array/ vector pair of top mass variables
	  for(Int_t igen=0; igen<ev.ngtop; igen++)
	    {
	      ntops++;
	      gen_tmass.push_back(ev.gtop_m[igen]);
	      gen_tids.push_back(ev.gtop_id[igen]);	      
	      ptsf *= TMath::Exp(0.156-0.00137*ev.gtop_pt[igen]);
	    }
	  if(ptsf>0 && ntops==2)
	    {
	      ptsf=TMath::Sqrt(ptsf);
	      topPtWgt[1]=1./ptsf;
	      topPtWgt[2]=ptsf;
	    }
	  
	  //update pileup weights, if found
	  if(puWgtGr.size())
	    {
	      puWeight[0]=puWgtGr[0]->Eval(ev.putrue);  
	      puWeight[1]=puWgtGr[1]->Eval(ev.putrue); 
	      puWeight[2]=puWgtGr[2]->Eval(ev.putrue);
	    }
	  
	  //update nominal event weight
	  float norm( normH ? normH->GetBinContent(1) : 1.0);
	  wgt=lepTriggerSF[0]*lepSelSF[0]*puWeight[0]*norm;
	  if(ev.ttbar_nw>0) wgt*=ev.ttbar_w[0];
	} else {
	for(Int_t igen=0; igen<ev.ngtop; igen++)
	  {
	    if(abs(ev.gtop_id[igen])!=6) continue;

	    gen_tmass.push_back(ev.gtop_m[igen]);
	    gen_tids.push_back(ev.gtop_id[igen]);
	  }
      }

      //nominal selection control histograms
      //int nJetsCat=TMath::Min((int)(bJets.size()+lightJets.size()),(int)4);
      //int nBtagsCat=TMath::Min((int)(bJets.size()),(int)2);
      std::vector<TString> catsToFill(1,Form("%s", lfs.Data()));
      //catsToFill[1]+= Form("%dT",nBtagsCat);
      //allPlots["puwgtctr"]->Fill(0.,puWeight[0]!=0 ? wgt/puWeight[0] : 0.);
      //allPlots["puwgtctr"]->Fill(1.,wgt);
      // if dilepton, at least 1 bjet TODO
      // if single-lep, at least 1 bjet + 3 jets total
      if((bJets.size() > 1 && (lfs == "EE" || lfs == "EM" || lfs == "MM")) ||
         ( bJets.size()+lightJets.size()>=3 && (lfs == "E" || lfs == "M")))
	{
	  //std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);
      // QUESTION: why run twice here?
	  for(size_t icat=0; icat<1; icat++)
	    {
	      TString tag=catsToFill[icat];

          //if(rIt!=lumiMap.end()) {
              //Int_t runCtr=std::distance(lumiMap.begin(),rIt);
              //allPlots["mlbwa_"+tag+"_ratevsrun"]->Fill(runCtr,1.e+6/rIt->second);
          //}

	      //allPlots["mlbwa_"+tag+"_L_Pt"]->Fill(isZ ? dilp4.Pt() : lp4.Pt(),wgt);
	      //allPlots["mlbwa_"+tag+"_L_Eta"]->Fill(isZ ? dilp4.Eta() : lp4.Eta(),wgt);
	      //allPlots["mlbwa_"+tag+"_J_Pt"]->Fill(ev.j_pt[leadingJetIdx],wgt);
	      //allPlots["mlbwa_"+tag+"_J_Eta"]->Fill(fabs(ev.j_eta[leadingJetIdx]),wgt);
	      //allPlots["mlbwa_"+tag+"_metphi"]->Fill(ev.met_phi[0],wgt);
	      //allPlots["mlbwa_"+tag+"_mttbar"]->Fill(visSystem.M(),wgt);
	      //allPlots["mlbwa_"+tag+"_TMass"]->Fill(mt,wgt);
	      allPlots["mlbwa_"+tag+"_Count"]->Fill(1, wgt);
	      allPlots["mlbwa_"+tag+"_MET"]->Fill(ev.met_pt[0],wgt); // metpt
	      allPlots["mlbwa_"+tag+"_B_Num"]->Fill(bJets.size(),wgt);
 
          // mlb storage
          //if(bJets.size()) {
          int ibjMax = (bJets.size() > 2 ? 2 : bJets.size());
          int ilpMax = (tightLeptonsIso.size() > 2 ? 2 : tightLeptonsIso.size());
          int itpMax = (gen_tmass.size() > 2 ? 2 : gen_tmass.size());

          float minmlb    = std::numeric_limits<float>::infinity();
          float minmlb_dr = std::numeric_limits<float>::infinity();

          int tbJetID_min(-1), tbJetID_mdr(-1),
              tLepChg_min(-1), tLepChg_mdr(-1);

          //std::cout << "\nFor file " << filename << "\n\t- ntops = " << gen_tmass.size()
          //          << "\n\t- ibjMax = " << ibjMax
          //          << "\n\t- ilpMax = " << ilpMax
          //          << "\n\t- isTTbar = " << isTTbar
          //          << "\n\t- is data = " << ev.isData
          //          << std::endl;

          for(int ibjet = 0; ibjet < ibjMax; ibjet++) {
              TLorentzVector tbJet = bJets[ibjet];
              int tbJetID = bJetIDs[ibjet];

              for(int ilp = 0; ilp < ilpMax; ilp++) {
                  TLorentzVector tLep(ev.l_pt[tightLeptonsIso[ilp]],
                          ev.l_eta[tightLeptonsIso[ilp]],
                          ev.l_phi[tightLeptonsIso[ilp]],
                          ev.l_mass[tightLeptonsIso[ilp]]);
                  float tmlb=(tbJet+tLep).M(); 
                  int tLepChg = ev.l_charge[ilp];

                  // min mlb calculations: minimum mass over all combinations
                  if(tmlb < minmlb) {
                      minmlb = tmlb;
                      tbJetID_min = tbJetID;
                      tLepChg_min = tLepChg;
                  }

                  // mdr mlb calculations: minimum mass when dR matched within 2.0
                  if(tbJet.DeltaR(tLep) <= 2.0 && tmlb < minmlb_dr) {
                      minmlb_dr = tmlb;
                      tbJetID_mdr = tbJetID;
                      tLepChg_mdr = tLepChg;
                  }

                  // loop over the top quarks

                  for(int itop = 0; itop < itpMax; itop++) {
                      float tmass   = gen_tmass.at(itop);
                      int   topID   = gen_tids.at(itop);
                      
                      // store the correctness of the inclusive entry
                      int lepTop_inc(-1), bJetTop_inc(-1);
                      if(abs(tbJetID) == 5 && topID * tbJetID > 0) bJetTop_inc = itop;
                      if(topID * ev.l_charge[tightLeptonsIso[ilp]] > 0) lepTop_inc = itop;
                      float pairTag_inc( bJetTop_inc<0 || lepTop_inc<0 ? 2.5 : // unmatched | cor | wrng 
                                                                   (bJetTop_inc==lepTop_inc ? 0.5 : 1.5 ));
                      allPlots["mlbwa_"+tag+"_lbmatch_inc"]->Fill(pairTag_inc); 
                      
                      // inc mlb calculations: add everything
                      allPlots["mlbwa_"+tag+"_Mlb_inc"]->Fill(tmlb,wgt); 
                      all2dPlots["mlbwa_"+tag+"_Tcmp_inc"]->Fill(tmlb,tmass,wgt);

                      // cheap trick: don't store min/mdr unless we're at the end of lep/bjet loops
                      // this avoids having to loop a second time through the top arrays
                      if(ilp == (ilpMax - 1) && ibjet == (ibjMax - 1)) {
                          
                          // check if our min pairing was correct
                          int lepTop_min(-1), bJetTop_min(-1);
                          if(abs(tbJetID_min) == 5 && topID * tbJetID_min > 0) bJetTop_min = itop;
                          if(topID * ev.l_charge[tightLeptonsIso[tLepChg_min]] > 0) lepTop_min = itop;
                          float pairTag_min( bJetTop_min<0 || lepTop_min<0 ? 2.5 :   // mismatch
                                              (bJetTop_min==lepTop_min ? 0.5 : 1.5 )); // cor | wrng
                          allPlots["mlbwa_"+tag+"_lbmatch_min"]->Fill(pairTag_min); 
                      
                          // add minmlb to our plots
                          if(!std::isinf(minmlb) && !std::isnan(minmlb)) {
                            allPlots["mlbwa_"+tag+"_Mlb_min"]->Fill(minmlb,wgt); 
                            all2dPlots["mlbwa_"+tag+"_Tcmp_min"]->Fill(minmlb,tmass,wgt);
                          }

                          // check if our mdr pairing was correct
                          int lepTop_mdr(-1), bJetTop_mdr(-1);
                          if(abs(tbJetID_mdr) == 5 && topID * tbJetID_mdr > 0) bJetTop_min = itop;
                          if(topID * ev.l_charge[tightLeptonsIso[tLepChg_mdr]] > 0) lepTop_min = itop;

                          float pairTag_mdr( bJetTop_mdr<0 || lepTop_mdr<0 ? 2.5 :   // mismatch
                                              (bJetTop_mdr==lepTop_mdr ? 0.5 : 1.5 )); // cor | wrng
                          allPlots["mlbwa_"+tag+"_lbmatch_mdr"]->Fill(pairTag_mdr); 

                          // add minmlb_dr to our plots (if there is such a value)
                          if(!std::isinf(minmlb_dr) && !std::isnan(minmlb_dr)) {
                              allPlots["mlbwa_"+tag+"_Mlb_mdr"]->Fill(minmlb_dr,wgt); 
                              all2dPlots["mlbwa_"+tag+"_Tcmp_mdr"]->Fill(minmlb_dr,tmass,wgt);
                          }
                      }
                  }
              }
          }



          //}
        }
	}

      //ANALYSIS WITH SYSTEMATICS
      if(!runSysts) continue;
      
      //gen weighting systematics
//      if(bJets.size()+lightJets.size()>=1 && normH)
//	{
//	  float mlb(bJets.size() ? (bJets[0]+lp4).M() : 0.);
//	  if(bJets.size()>1) mlb=TMath::Min( (float) mlb, (float)(bJets[1]+lp4).M() );
//	  
//	  for(Int_t igen=0; igen<nGenSysts; igen++)
//	    {
//	      for(size_t icat=0; icat<2; icat++)
//		{
//		  float newWgt = wgt*ev.ttbar_w[igen]/ev.ttbar_w[0];
//
//		  //for signal we only consider shapes and acceptance effects as it will be fit
//		  if(isTTbar) 
//		    newWgt *= normH->GetBinContent(igen+1)/normH->GetBinContent(1);
//
//		  TString tag=catsToFill[icat];	 
//		  all2dPlots["metptshapes_"+tag+"_gen"]->Fill(ev.met_pt[0],igen,newWgt);
//		  all2dPlots["mlbshapes_"+tag+"_gen"]->Fill(mlb,igen,newWgt);
//		  if(alpha<0.3) all2dPlots["RMPFshapes_"+tag+"_gen"]->Fill(RMPF,igen,newWgt);
//		  if(icat==0) all2dPlots["nbtagsshapes_"+tag+"_gen"]->Fill(bJets.size(),igen,newWgt);
//		}
//	    }
//	}
//
//      //experimental systematics
//      for(size_t ivar=0; ivar<expSysts.size(); ivar++)
//	{
//	  TString varName=expSysts[ivar];
//	  bool updateJEn(ivar<=jecUncSrcs.size());
//	  bool updateBtag(varName.Contains("tagEff"));
//
//	  //update jet selection, if needed
//	  for(int isign=0; isign<2; isign++)
//	    {
//	      //weight systematics
//	      float newWgt(wgt);
//	      if(varName=="Pileup" && puWeight[0]!=0) newWgt*=(isign==0 ? puWeight[1]/puWeight[0] : puWeight[2]/puWeight[0]);
//	      if(
//		 (varName=="MuTrigger" && (lid==13 || lid==1300 || lid==2113)) ||
//		 (varName=="EleTrigger" && (lid==11 || lid==1100 || lid==2111))
//		 )
//		newWgt *= (isign==0 ? lepTriggerSF[1]/lepTriggerSF[0] : lepTriggerSF[2]/lepTriggerSF[0]);
//	      if(
//		 (varName=="MuEfficiency" && (lid==13 ||  lid==1300 || lid==2113)) ||
//		 (varName=="EleEfficiency" && (lid==11 || lid==1100 || lid==2111))
//		 )
//		newWgt *= (isign==0 ? lepSelSF[1]/lepSelSF[0] : lepSelSF[2]/lepSelSF[0]);
//	      if(
//		 varName=="topPt"
//		 )
//		newWgt *= topPtWgt[1+isign];
//
//	      //lepton scale systematics
//	      TLorentzVector varlp4(lp4),varl2p4(0,0,0,0),vardilp4(0,0,0,0);
//	      if(isZ)
//		varl2p4.SetPtEtaPhiM(ev.l_pt[tightLeptonsIso[1]],ev.l_eta[tightLeptonsIso[1]],ev.l_phi[tightLeptonsIso[1]],ev.l_mass[tightLeptonsIso[1]]);		
//	      if(
//		 (varName=="MuScale" && ( lid==13 || lid==1300 || lid==2113)) ||
//		 (varName=="EleScale" && (lid==11 || lid==1100 || lid==2111))
//		 )
//		{
//		  varlp4 = (1.0+(isign==0?-1.:1.)*getLeptonEnergyScaleUncertainty(lid,lp4.Pt(),lp4.Eta()) ) *lp4;
//		  if(isZ) varl2p4 = (1.0+(isign==0?-1.:1.)*getLeptonEnergyScaleUncertainty(lid,varl2p4.Pt(),varl2p4.Eta()) ) *varl2p4;
//		}
//	      if(varlp4.Pt()<30 || fabs(varlp4.Eta())>2.1) continue;
//	      if(isZ) 
//		{
//		  if(varl2p4.Pt()<30 || fabs(varl2p4.Eta())>2.1) continue;
//		  vardilp4=varlp4+varl2p4;
//		  if(fabs(vardilp4.M()-91)>10) continue;
//		  if(vardilp4.Pt()<30) continue;
//		}
//
//	      //jets
//	      std::vector<TLorentzVector> varBJets,varLightJets;
//	      TLorentzVector jetDiff(0,0,0,0),jetSum(0,0,0,0);
//	      if(!(updateJEn || updateBtag)) { varBJets=bJets; varLightJets=lightJets; }
//	      else
//		{
//		  for (size_t ij=0; ij<resolvedJetIdx.size();ij++)
//		    {
//		      int k(resolvedJetIdx[ij]);
//		      int jflav( abs(ev.j_hadflav[k]) ),jflavForJES( abs(ev.j_flav[k]) );
//		      
//		      //check kinematics
//		      TLorentzVector jp4;
//		      jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
//		      //jp4=updateJES(jp4,ev.j_rawsf[k],ev.j_area[k],ev.rho,ev.nvtx,jetCorr);
//		      jetDiff -= jp4;
//		      
//		      //smear jet energy resolution for MC
//		      float genJet_pt(ev.genj_pt[k]); 
//		      if(!ev.isData && genJet_pt>0) 
//			{
//			  std::vector<float> jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt);
//			  if(varName=="JER") jp4 *= jerSmear[isign+1];
//			  else               jp4 *= jerSmear[0];
//			}
//		      
//		      //update jet energy correction
//		      if(updateJEn && varName!="JER")
//			{
//			  jecUncs[ivar]->setJetEta(fabs(jp4.Eta()));
//			  jecUncs[ivar]->setJetPt(jp4.Pt());
//			  if(
//			     (jflavForJES==4 && varName=="FlavorPureCharm") || 
//			     (jflavForJES==5 && varName=="FlavorPureBottom") || 
//			     (jflavForJES==21 && varName=="FlavorPureGluon") || 
//			     (jflavForJES<4 && varName=="FlavorPureQuark") || 
//			     !varName.Contains("FlavorPure") 
//			     )
//			    {
//			      float jecuncVal=jecUncs[ivar]->getUncertainty(true);
//			      jp4 *= (1.0+(isign==0?-1.:1.)*jecuncVal);
//			    }
//			}
//		      
//		      jetDiff += jp4;
//		      jetSum += jp4;
//
//		      //re-apply kinematics
//		      if(jp4.Pt()<30) continue;
//		      if(fabs(jp4.Eta()) > 2.4) continue;
//		      
//		      //b-tag
//		      float csv = ev.j_csv[k];	  
//		      bool isBTagged(csv>0.890);
//		      if(!ev.isData)
//			{
//			  float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
//			  float expEff(1.0), jetBtagSF(1.0);
//			  if(jflav==4) 
//			    { 
//			      expEff        = expBtagEff["c"]->Eval(jptForBtag); 
//			      int idx(0);
//			      if(varName=="CtagEff") idx=(isign==0 ? 1 : 2);
//			      jetBtagSF     = sfbReaders[idx]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
//			    }
//			  else if(jflav==5)
//			    { 
//			      expEff=expBtagEff["b"]->Eval(jptForBtag); 
//			      int idx(0);
//			      if(varName=="BtagEff") idx=(isign==0 ? 1 : 2);
//			      jetBtagSF=sfbReaders[idx]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
//			    }
//			  else
//			    {
//			      expEff=expBtagEff["udsg"]->Eval(jptForBtag);
//			      int idx(0);
//			      if(varName=="LtagEff") idx=(isign==0 ? 1 : 2);
//			      jetBtagSF=sflReaders[idx]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
//			    }
//	      
//			  //updated b-tagging decision with the data/MC scale factor
//			  myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
//			}
//
//		      //save jet
//		      if(isBTagged) varBJets.push_back(jp4);
//		      else          varLightJets.push_back(jp4);
//		    }
//		}
//	      
//	      //MET variations
//	      TLorentzVector varMet(met);
//	      varMet += jetDiff;
//	      varMet += (varlp4-lp4);
//	      if(varName=="UncMET") 
//		{
//		  TLorentzVector uncMet(met+varlp4+jetSum);
//		  uncMet *= (1.0 + (isign==0?-1.:1.)*0.1);
//		  varMet = uncMet-varlp4-jetSum;
//		}
//	      varMet.SetPz(0.); varMet.SetE(varMet.Pt());
//	      // float varmt(computeMT(varlp4,varMet) );
//
//	      if(varBJets.size()+varLightJets.size()<1) continue;
//
//	      //ready to fill histograms
//	      float mlb(varBJets.size() ? (varBJets[0]+varlp4).M() : 0.);
//	      if(varBJets.size()>1) mlb=TMath::Min( (float) mlb, (float)(varBJets[1]+varlp4).M() );
//
//	      //balancing variable
//	      float varRMPF(0);
//	      TVector2 metT(varMet.Px(),varMet.Py());
//	      TVector2 visT(isZ ? vardilp4.Px() : varlp4.Px(), isZ ? vardilp4.Py() : varlp4.Py());
//	      varRMPF=1.0+(metT*visT)/visT.Mod2();
//	
//	      //update categories
//	      //int nvarJetsCat=TMath::Min((int)(varBJets.size()+varLightJets.size()),(int)4);
//	      //int nvarBtagsCat=TMath::Min((int)(varBJets.size()),(int)2);
//	      std::vector<TString> varcatsToFill(2,Form("%s", lfs.Data()));
//	      //varcatsToFill[1]+= Form("%dT",nvarBtagsCat);
//
//	      for(size_t icat=0; icat<2; icat++)
//                {
//		  TString tag=varcatsToFill[icat];
//		  all2dPlots["metptshapes_"+tag+"_exp"]->Fill(varMet.Pt(),2*ivar+isign,newWgt);
//		  all2dPlots["mlbshapes_"+tag+"_exp"]->Fill(mlb,2*ivar+isign,newWgt);
//		  all2dPlots["RMPFshapes_"+tag+"_exp"]->Fill(varRMPF,2*ivar+isign,newWgt);
//		  if(icat!=0) continue;
//		  all2dPlots["nbtagsshapes_"+tag+"_exp"]->Fill(varBJets.size(),2*ivar+isign,newWgt);
//		}
//	    }
//	}
    }


  // APPLY INTERPOLATION 
  if(interpolate && !ev.isData && isTTbar) {
      std::vector<TString> mlbVec = {"min", "mdr", "inc"};

      // loop over final states and kinds of mlb
      for(size_t imlb = 0; imlb < mlbVec.size(); imlb++) {
      for(size_t ilfs = 0; ilfs < lfsVec.size(); ilfs++) {
        TString mlb = mlbVec.at(imlb);
        TString lfs = lfsVec.at(ilfs);

        // project 2dplot onto tmass axis

        TH2F *th2 = (TH2F*) all2dPlots["mlbwa_"+lfs+"_Tcmp_"+mlb];
        TH1F *tmass_hist = (TH1F*) all2dPlots["mlbwa_"+lfs+"_Tcmp_"+mlb]->ProjectionY("mlbwa_"+lfs+"_"+mlb+"_wts");
        
        // fit projected plots with rbw
        
        TF1 *bw = new TF1("bw", (TString("[0]*([1]*[2]*sqrt([1]*[1]*([1]*[1]+[2]*[2]))")
                                 +TString("/sqrt([1]*[1]+sqrt([1]*[1]*([1]*[1]+[2]*[2]))))")
                                 +TString("/(TMath::Power(x*x-[1]*[1],2)+TMath::Power([1]*[2],2))")).Data(),
                           160, 200);
        bw->SetParName(0,"N");
        bw->SetParName(1,"m_{t}");
        bw->SetParName(2,"#Gamma_{t}");

        bw->SetParLimits(1,160,200);
        bw->SetParLimits(2,0.01,10);

        bw->SetParameter(0,1);
        bw->SetParameter(1,172);
        bw->SetParameter(2,1.3);

        tmass_hist->Fit(bw,"", "", 160, 200);

        //loop over interpolations
        for(int iInt = 0; iInt <= nIntrps + 1; iInt++) {

            // setup histos for storage
            float intrpWid = iInt*(maxIntrpWid - minIntrpWid)/(nIntrps + 1) + minIntrpWid;

            char intrpWidStr[256];
            sprintf(intrpWidStr, "mlbwa_%s_%.2f_Mlb_%s", lfs.Data(), intrpWid, mlb.Data());
            TString intrpWidTStr = TString(intrpWidStr).ReplaceAll(".", "p");

            TH1F* intrp_mlb_histo = new TH1F(intrpWidTStr, intrpWidTStr, 70, 0, 350);
            
            // calculate BW for interpolated width
            TF1 *bwInterp = (TF1*) bw->Clone();
            bwInterp->SetParameter(2, intrpWid);

            // loop over mlb bins
            for(int mlb_bin = 0; mlb_bin < 70; mlb_bin++) {
                float weighted_mlb(0);

                // loop over tmass bins
                for(int tmass_bin = 0; tmass_bin < 75; tmass_bin++) {
                   float wt_tmass = th2->GetYaxis()->GetBinCenter(tmass_bin);
                   float unweighted_mlb = th2->GetBinContent(mlb_bin,tmass_bin);
                    
                   if(bw->Eval(wt_tmass) != 0) {
                       weighted_mlb += unweighted_mlb * bwInterp->Eval(wt_tmass) / bw->Eval(wt_tmass);
                   }
                }

                intrp_mlb_histo->AddBinContent(mlb_bin, weighted_mlb);
            }

            // store weighted histograms
            allPlots[intrpWidTStr] = intrp_mlb_histo;
        } 
        
      }}


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
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : all2dPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}
