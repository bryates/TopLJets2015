#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/ReadTree.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

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

Float_t computeMT(TLorentzVector &a, TLorentzVector &b)
{
  return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
}


void ReadTree(TString filename,
	      TString outname,
	      Int_t channelSelection, 
	      Int_t chargeSelection, 
	      FlavourSplitting flavourSplitting,
	      TH1F *normH, 
	      Bool_t runSysts)
{

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
  TString jecUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/Summer15_25nsV6M3_DATA_UncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
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
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",2,0,2);
  for(Int_t ij=1; ij<=4; ij++)
    {
      for(Int_t itag=-1; itag<=2; itag++)
	{
	  if(itag>ij) continue;
	  TString tag(itag<0 ? Form("%dj",ij) : Form("%dj%dt",ij,itag));
	  if(lumiMap.size()) allPlots["ratevsrun_"+tag] = new TH1F("ratevsrun_"+tag,";Run number; Events/pb",lumiMap.size(),0,lumiMap.size());
	  Int_t runCtr(0);
	  for(std::map<Int_t,Float_t>::iterator it=lumiMap.begin(); it!=lumiMap.end(); it++,runCtr++)
	    allPlots["ratevsrun_"+tag]->GetXaxis()->SetBinLabel(runCtr+1,Form("%d",it->first));
	  allPlots["lpt_"+tag]        = new TH1F("lpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
	  allPlots["lsip3d_"+tag]     = new TH1F("lsip3d_"+tag,";3d impact parameter significance;Events" ,40,0.,4.);
	  allPlots["lreliso_"+tag]     = new TH1F("lreliso_"+tag,";Relative isolation;Events" ,25,0.,0.5);
	  allPlots["leta_"+tag]       = new TH1F("leta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["jpt_"+tag]        = new TH1F("jpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
	  allPlots["jeta_"+tag]       = new TH1F("jeta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["ht_"+tag]         = new TH1F("ht_"+tag,";H_{T} [GeV];Events",40,0,800);
	  allPlots["csv_"+tag]        = new TH1F("csv_"+tag,";CSV discriminator;Events",100,0,1.0);
	  allPlots["nvtx_"+tag]       = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events" ,50,0.,50.);
	  allPlots["metpt_"+tag]      = new TH1F("metpt_"+tag,";Missing transverse energy [GeV];Events" ,10,0.,200.);
	  allPlots["metphi_"+tag]     = new TH1F("metphi_" + tag,";MET #phi [rad];Events" ,50,-3.2,3.2);
	  allPlots["mttbar_"+tag]     = new TH1F("mttbar_"+tag,";#sqrt{#hat{s}} [GeV];Events" ,50,0.,1000.);
	  allPlots["mt_"+tag]         = new TH1F("mt_"+tag,";Transverse Mass [GeV];Events" ,20,0.,200.);
	  allPlots["minmlb_"+tag]     = new TH1F("minmlb_"+tag,";min Mass(lepton,b) [GeV];Events" ,25,0.,250.);
	  allPlots["passSIP3d_"+tag]     = new TH1F("passSIP3d_"+tag,";Pass SIP3d requirement;Events" ,2,0.,2.);
	  allPlots["RMPF_"+tag]       = new TH1F("RMPF_"+tag,";R_{MPF};Events" ,20,0.,2.);
	  allPlots["alpha_"+tag]      = new TH1F("alpha_"+tag,";#alpha;Events" ,20,0.,2.);
	  if(itag==-1)
	    {
	      allPlots["nbtags_"+tag]     = new TH1F("nbtags_"+tag,";Category;Events" ,3,0.,3.);
	      allPlots["nbtags_"+tag]->GetXaxis()->SetBinLabel(1, "=0b");
	      allPlots["nbtags_"+tag]->GetXaxis()->SetBinLabel(2, "=1b");
	      allPlots["nbtags_"+tag]->GetXaxis()->SetBinLabel(3, "#geq2b");
	    }

	  //shape analysis
	  if(runSysts)
	    {
	      //gen level systematics
	      if(normH)
		{
		  all2dPlots["metptshapes_"+tag+"_gen"]                   
		    = new TH2F("metptshapes_"+tag+"_gen", ";Missing transverse energy [GeV];Events" ,    10,0.,200., nGenSysts,0,nGenSysts);
		  all2dPlots["minmlbshapes_"+tag+"_gen"]               
		    = new TH2F("minmlbshapes_"+tag+"_gen", ";min Mass(lepton,b) [GeV];Events" , 25,0.,250., nGenSysts,0,nGenSysts);
		  all2dPlots["RMPFshapes_"+tag+"_gen"]               
		    = new TH2F("RMPFshapes_"+tag+"_gen", ";R_{MPF};Events" , 20,0.,2., nGenSysts,0,nGenSysts);
		  if(itag==-1) 
		    all2dPlots["nbtagsshapes_"+tag+"_gen"] 
		      = new TH2F("nbtagsshapes_"+tag+"_gen", ";Category;Events" , 3, 0.,3.,   nGenSysts,0,nGenSysts);		  
		  for(Int_t igen=0; igen<nGenSysts; igen++)
		    {
		      TString label(normH->GetXaxis()->GetBinLabel(igen+1));
		      all2dPlots["metptshapes_"+tag+"_gen"]    ->GetYaxis()->SetBinLabel(igen+1,label);
		      all2dPlots["minmlbshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,label);
		      all2dPlots["RMPFshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,label);
		      if(itag!=-1) continue;
		      all2dPlots["nbtagsshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,label);
		    }
		}
	      
	      //experimental systematics
	      Int_t nExpSysts=expSysts.size();
	      if(nExpSysts>0)
		{
		  all2dPlots["metptshapes_"+tag+"_exp"]                  
		    = new TH2F("metptshapes_"+tag+"_exp",  ";Missing transverse energy [GeV];Events" ,   10,0.,200., 2*nExpSysts,0,2*nExpSysts);
		  all2dPlots["minmlbshapes_"+tag+"_exp"]              
		    = new TH2F("minmlbshapes_"+tag+"_exp", ";min Mass(lepton,b) [GeV];Events" ,25,0.,250., 2*nExpSysts,0,2*nExpSysts);
		  all2dPlots["RMPFshapes_"+tag+"_exp"]               
		    = new TH2F("RMPFshapes_"+tag+"_exp", ";R_{MPF};Events" , 20,0.,2., 2*nExpSysts,0,2*nExpSysts);
		  if(itag==-1) 
		    all2dPlots["nbtagsshapes_"+tag+"_exp"] 
		      = new TH2F("nbtagsshapes_"+tag+"_exp", ";Category;Events" ,  3,0.,3.,    2*nExpSysts,0,2*nExpSysts);
		  for(Int_t isyst=0; isyst<nExpSysts; isyst++)
		    {
		      for(Int_t ivar=0; ivar<2; ivar++)
			{
			  TString label(expSysts[isyst] + (ivar==0 ? "Down" : "Up"));
			  all2dPlots["metptshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
			  all2dPlots["minmlbshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
			  all2dPlots["RMPFshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
			  if(itag!=-1) continue;
			  all2dPlots["nbtagsshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
			}
		    }
		}
	    }
	}
    }

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
      Bool_t isZ(false),isZPassingSIP3d(false);
      TLorentzVector l1p4,l2p4,dilp4;
      if(tightLeptonsIso.size()==1)                                       lepIdx=tightLeptonsIso[0];
      else if (tightLeptonsIso.size()==0 && tightLeptonsNonIso.size()==1) lepIdx=tightLeptonsNonIso[0];
      else if(tightLeptonsIso.size()==2)
	{	  
	  l1p4.SetPtEtaPhiM(ev.l_pt[tightLeptonsIso[0]],ev.l_eta[tightLeptonsIso[0]],ev.l_phi[tightLeptonsIso[0]],ev.l_mass[tightLeptonsIso[0]]);
	  l2p4.SetPtEtaPhiM(ev.l_pt[tightLeptonsIso[1]],ev.l_eta[tightLeptonsIso[1]],ev.l_phi[tightLeptonsIso[1]],ev.l_mass[tightLeptonsIso[1]]);
	  dilp4=l1p4+l2p4;
	  if(ev.l_id[tightLeptonsIso[0]]==ev.l_id[tightLeptonsIso[1]]          && 
	     ev.l_charge[tightLeptonsIso[0]]*ev.l_charge[tightLeptonsIso[1]]<0 && 
	     fabs(dilp4.M()-91)<10 && 
	     dilp4.Pt()>30)
	    { 
	      isZ=true; 
	      isZPassingSIP3d=(ev.l_ip3dsig[0]<4 && ev.l_ip3dsig[1]<4);
	    }
	  lepIdx=tightLeptonsIso[0];
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
      Int_t lid=ev.l_id[lepIdx];
      if(isZ) lid=2100+ev.l_id[lepIdx];
      else if(tightLeptonsNonIso.size()==1) lid=100*ev.l_id[lepIdx];
      if(channelSelection!=0)
	{
	  if(channelSelection==21) { if(!isZ) continue; }
	  else                     { if(lid!=channelSelection) continue; }
	}
      if(chargeSelection!=0  && ev.l_charge[lepIdx]!=chargeSelection) continue;

      //lepton kinematics
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);

      //select jets
      Float_t htsum(0);
      TLorentzVector jetDiff(0,0,0,0);
      std::vector<TLorentzVector> bJets,lightJets;
      TLorentzVector visSystem(isZ ? dilp4 : lp4);
      int nbjets(0),ncjets(0),nljets(0),leadingJetIdx(-1);
      std::vector<int> resolvedJetIdx;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

	  //cross clean with respect to leptons 
	  if(jp4.DeltaR(lp4)<0.4) continue;
	  if(isZ && jp4.DeltaR(l2p4)<0.4)continue;

	  //smear jet energy resolution for MC
	  jetDiff -= jp4;
	  float genJet_pt(ev.genj_pt[k]); 
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }
	  jetDiff += jp4;
	  resolvedJetIdx.push_back(k);

	  //require back-to-back configuration with Z
	  if(isZ && jp4.DeltaPhi(dilp4)<2.7) continue;

	  // re-inforce kinematics cuts
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  if(leadingJetIdx<0) leadingJetIdx=k;
	  htsum += jp4.Pt();
	  if(bJets.size()+lightJets.size()<4) visSystem += jp4;

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
	  if(isBTagged) bJets.push_back(jp4);
	  else          lightJets.push_back(jp4);
	}

      //check if flavour splitting was required
      if(!ev.isData)
	{
	  if(flavourSplitting!=FlavourSplitting::NOFLAVOURSPLITTING)
	    {
	      if(flavourSplitting==FlavourSplitting::BSPLITTING)         { if(nbjets==0)    continue; }
	      else if(flavourSplitting==FlavourSplitting::CSPLITTING)    { if(ncjets==0 || nbjets!=0)    continue; }
	      else if(flavourSplitting==FlavourSplitting::UDSGSPLITTING) { if(nljets==0 || ncjets!=0 || nbjets!=0) continue; }
	    }
	}

      //MET and transverse mass
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
      met+=jetDiff;
      met.SetPz(0.); met.SetE(met.Pt());
      float mt( computeMT(isZ ? dilp4: lp4,met) );

      //compute neutrino kinematics
      neutrinoPzComputer.SetMET(met);
      neutrinoPzComputer.SetLepton(isZ ? dilp4 : lp4);
      
      float nupz=neutrinoPzComputer.Calculate();
      TLorentzVector neutrinoHypP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));
      visSystem+=neutrinoHypP4;

      //balancing variable and ISR control variable
      float RMPF(0.0),alpha(0.0);
      TVector2 metT(met.Px(),met.Py());
      TVector2 visT(isZ ? dilp4.Px() : lp4.Px(), isZ ? dilp4.Py() : lp4.Py());
      RMPF=1.0+(metT*visT)/visT.Mod2();
      for(size_t ij=0; ij<resolvedJetIdx.size(); ij++)
	{
	  Int_t k=resolvedJetIdx[ij];
	  if(k==leadingJetIdx) continue;
	  if(ev.j_pt[k]<15 || fabs(ev.j_eta[k])>3.0) continue;
	  alpha=ev.j_pt[k]/(isZ ? dilp4.Pt() : lp4.Pt());
	}
      
      //event weight
      float wgt(1.0);
      std::vector<float> puWeight(3,1.0),lepTriggerSF(3,1.0),lepSelSF(3,1.0), topPtWgt(3,1.0);
      if(!ev.isData)
	{
	  //update lepton selection scale factors, if found
	  TString prefix("m");
	  if(lid==11 || lid==1100) prefix="e";
	  if(lepEffH.find(prefix+"_sel")!=lepEffH.end() && !isZ)
	    {
	      for(size_t il=0; il<tightLeptonsIso.size(); il++)
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
		  
		  float trigSF(lepEffH[prefix+"_trig"]->GetBinContent(etaBinForEff,ptBinForEff));
		  float trigSFUnc(lepEffH[prefix+"_trig"]->GetBinError(etaBinForEff,ptBinForEff));

		  lepTriggerSF[0]*=trigSF; lepTriggerSF[1]*=(trigSF-trigSFUnc); lepTriggerSF[2]*=(trigSF+trigSFUnc);
		  lepSelSF[0]*=selSF;      lepSelSF[1]*=(selSF-selSFUnc);       lepSelSF[2]*=(selSF+selSFUnc);
		}
	    }
	  
	  Int_t ntops(0);
	  float ptsf(1.0);
	  for(Int_t igen=0; igen<ev.ngenHardProc; igen++)
	    {
	      if(abs(ev.ghp_id[igen])!=6) continue;
	      ntops++;
	      ptsf *= TMath::Exp(0.156-0.00137*ev.ghp_pt[igen]);
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
	}

      //nominal selection control histograms
      int nJetsCat=TMath::Min((int)(bJets.size()+lightJets.size()),(int)4);
      int nBtagsCat=TMath::Min((int)(bJets.size()),(int)2);
      std::vector<TString> catsToFill(2,Form("%dj",nJetsCat));
      catsToFill[1]+= Form("%dt",nBtagsCat);
      allPlots["puwgtctr"]->Fill(0.,puWeight[0]!=0 ? wgt/puWeight[0] : 0.);
      allPlots["puwgtctr"]->Fill(1.,wgt);
      if(bJets.size()+lightJets.size()>=1)
	{
	  std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);	  
	  for(size_t icat=0; icat<2; icat++)
	    {
	      TString tag=catsToFill[icat];
	      if(rIt!=lumiMap.end()) 
		{
		  Int_t runCtr=std::distance(lumiMap.begin(),rIt);
		  allPlots["ratevsrun_"+tag]->Fill(runCtr,1.e+6/rIt->second);
		}

	      allPlots["lpt_"+tag]->Fill(isZ ? dilp4.Pt() : lp4.Pt(),wgt);
	      allPlots["leta_"+tag]->Fill(isZ ? dilp4.Eta() : lp4.Eta(),wgt);
	      allPlots["lsip3d_"+tag]->Fill(ev.l_ip3dsig[lepIdx],wgt);
	      allPlots["passSIP3d_"+tag]->Fill(isZPassingSIP3d,wgt);
	      allPlots["lreliso_"+tag]->Fill(ev.l_relIso[lepIdx],wgt);
	      allPlots["jpt_"+tag]->Fill(ev.j_pt[ leadingJetIdx ],wgt);
	      allPlots["jeta_"+tag]->Fill(fabs(ev.j_eta[ leadingJetIdx ]),wgt);
	      allPlots["csv_"+tag]->Fill(ev.j_csv[ leadingJetIdx ],wgt);
	      allPlots["nvtx_"+tag]->Fill(ev.nvtx,wgt);
	      allPlots["metpt_"+tag]->Fill(ev.met_pt[0],wgt);
	      allPlots["metphi_"+tag]->Fill(ev.met_phi[0],wgt);
	      allPlots["mt_"+tag]->Fill(mt,wgt);
	      allPlots["mttbar_"+tag]->Fill(visSystem.M(),wgt);
	      allPlots["ht_"+tag]->Fill(htsum,wgt);
	      allPlots["alpha_"+tag]->Fill(alpha,wgt);
	      if(alpha<0.3) allPlots["RMPF_"+tag]->Fill(RMPF,wgt);
	      if(icat==0) allPlots["nbtags_"+tag]->Fill(bJets.size(),wgt);
	     
	      if(bJets.size())
		{
		  float mlb=(bJets[0]+(isZ ? dilp4 : lp4)).M();		  
		  if(bJets.size()>1) mlb=TMath::Min( (float) mlb, (float)(bJets[1]+(isZ ? dilp4 : lp4)).M() );
		  allPlots["minmlb_"+tag]->Fill(mlb,wgt);		 		 
		}	  
	    }
	}

      //ANALYSIS WITH SYSTEMATICS
      if(!runSysts) continue;
      
      //gen weighting systematics
      if(bJets.size()+lightJets.size()>=1 && normH)
	{
	  float mlb(bJets.size() ? (bJets[0]+lp4).M() : 0.);
	  if(bJets.size()>1) mlb=TMath::Min( (float) mlb, (float)(bJets[1]+lp4).M() );
	  
	  for(Int_t igen=0; igen<nGenSysts; igen++)
	    {
	      for(size_t icat=0; icat<2; icat++)
		{
		  float newWgt = wgt*ev.ttbar_w[igen]/ev.ttbar_w[0];
		  TString tag=catsToFill[icat];	 
		  all2dPlots["metptshapes_"+tag+"_gen"]->Fill(ev.met_pt[0],igen,newWgt);
		  all2dPlots["minmlbshapes_"+tag+"_gen"]->Fill(mlb,igen,newWgt);
		  if(alpha<0.3) all2dPlots["RMPFshapes_"+tag+"_gen"]->Fill(RMPF,igen,newWgt);
		  if(icat==0) all2dPlots["nbtagsshapes_"+tag+"_gen"]->Fill(bJets.size(),igen,newWgt);
		}
	    }
	}

      //experimental systematics
      for(size_t ivar=0; ivar<expSysts.size(); ivar++)
	{
	  TString varName=expSysts[ivar];
	  bool updateJEn(ivar<=jecUncSrcs.size());
	  bool updateBtag(varName.Contains("tagEff"));

	  //update jet selection, if needed
	  for(int isign=0; isign<2; isign++)
	    {
	      //weight systematics
	      float newWgt(wgt);
	      if(varName=="Pileup" && puWeight[0]!=0) newWgt*=(isign==0 ? puWeight[1]/puWeight[0] : puWeight[2]/puWeight[0]);
	      if(
		 (varName=="MuTrigger" && (lid==13 || lid==1300 || lid==2113)) ||
		 (varName=="EleTrigger" && (lid==11 || lid==1100 || lid==2111))
		 )
		newWgt *= (isign==0 ? lepTriggerSF[1]/lepTriggerSF[0] : lepTriggerSF[2]/lepTriggerSF[0]);
	      if(
		 (varName=="MuEfficiency" && (lid==13 ||  lid==1300 || lid==2113)) ||
		 (varName=="EleEfficiency" && (lid==11 || lid==1100 || lid==2111))
		 )
		newWgt *= (isign==0 ? lepSelSF[1]/lepSelSF[0] : lepSelSF[2]/lepSelSF[0]);
	      if(
		 varName=="topPt"
		 )
		newWgt *= topPtWgt[1+isign];

	      //lepton scale systematics
	      TLorentzVector varlp4(lp4),varl2p4(0,0,0,0),vardilp4(0,0,0,0);
	      if(isZ)
		varl2p4.SetPtEtaPhiM(ev.l_pt[tightLeptonsIso[1]],ev.l_eta[tightLeptonsIso[1]],ev.l_phi[tightLeptonsIso[1]],ev.l_mass[tightLeptonsIso[1]]);		
	      if(
		 (varName=="MuScale" && ( lid==13 || lid==1300 || lid==2113)) ||
		 (varName=="EleScale" && (lid==11 || lid==1100 || lid==2111))
		 )
		{
		  varlp4 = (1.0+(isign==0?-1.:1.)*getLeptonEnergyScaleUncertainty(lid,lp4.Pt(),lp4.Eta()) ) *lp4;
		  if(isZ) varl2p4 = (1.0+(isign==0?-1.:1.)*getLeptonEnergyScaleUncertainty(lid,varl2p4.Pt(),varl2p4.Eta()) ) *varl2p4;
		}
	      if(varlp4.Pt()<30 || fabs(varlp4.Eta())>2.1) continue;
	      if(isZ) 
		{
		  if(varl2p4.Pt()<30 || fabs(varl2p4.Eta())>2.1) continue;
		  vardilp4=varlp4+varl2p4;
		  if(fabs(vardilp4.M()-91)>10) continue;
		  if(vardilp4.Pt()<30) continue;
		}

	      //jets
	      std::vector<TLorentzVector> varBJets,varLightJets;
	      TLorentzVector jetDiff(0,0,0,0),jetSum(0,0,0,0);
	      if(!(updateJEn || updateBtag)) { varBJets=bJets; varLightJets=lightJets; }
	      else
		{
		  for (size_t ij=0; ij<resolvedJetIdx.size();ij++)
		    {
		      int k(resolvedJetIdx[ij]);
		      int jflav( abs(ev.j_hadflav[k]) ),jflavForJES( abs(ev.j_flav[k]) );
		      
		      //check kinematics
		      TLorentzVector jp4;
		      jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
		      jetDiff -= jp4;
		      
		      //smear jet energy resolution for MC
		      float genJet_pt(ev.genj_pt[k]); 
		      if(!ev.isData && genJet_pt>0) 
			{
			  std::vector<float> jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt);
			  if(varName=="JER") jp4 *= jerSmear[isign+1];
			  else               jp4 *= jerSmear[0];
			}
		      
		      //update jet energy correction
		      if(updateJEn && varName!="JER")
			{
			  jecUncs[ivar]->setJetEta(fabs(jp4.Eta()));
			  jecUncs[ivar]->setJetPt(jp4.Pt());
			  if(
			     (jflavForJES==4 && varName=="FlavorPureCharm") || 
			     (jflavForJES==5 && varName=="FlavorPureBottom") || 
			     (jflavForJES==21 && varName=="FlavorPureGluon") || 
			     (jflavForJES<4 && varName=="FlavorPureQuark") || 
			     !varName.Contains("FlavorPure") 
			     )
			    {
			      float jecuncVal=jecUncs[ivar]->getUncertainty(true);
			      jp4 *= (1.0+(isign==0?-1.:1.)*jecuncVal);
			    }
			}
		      
		      jetDiff += jp4;
		      jetSum += jp4;

		      //re-apply kinematics
		      if(jp4.Pt()<30) continue;
		      if(fabs(jp4.Eta()) > 2.4) continue;
		      
		      //b-tag
		      float csv = ev.j_csv[k];	  
		      bool isBTagged(csv>0.890);
		      if(!ev.isData)
			{
			  float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
			  float expEff(1.0), jetBtagSF(1.0);
			  if(jflav==4) 
			    { 
			      expEff        = expBtagEff["c"]->Eval(jptForBtag); 
			      int idx(0);
			      if(varName=="CtagEff") idx=(isign==0 ? 1 : 2);
			      jetBtagSF     = sfbReaders[idx]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
			    }
			  else if(jflav==5)
			    { 
			      expEff=expBtagEff["b"]->Eval(jptForBtag); 
			      int idx(0);
			      if(varName=="BtagEff") idx=(isign==0 ? 1 : 2);
			      jetBtagSF=sfbReaders[idx]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
			    }
			  else
			    {
			      expEff=expBtagEff["udsg"]->Eval(jptForBtag);
			      int idx(0);
			      if(varName=="LtagEff") idx=(isign==0 ? 1 : 2);
			      jetBtagSF=sflReaders[idx]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
			    }
	      
			  //updated b-tagging decision with the data/MC scale factor
			  myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
			}

		      //save jet
		      if(isBTagged) varBJets.push_back(jp4);
		      else          varLightJets.push_back(jp4);
		    }
		}
	      
	      //MET variations
	      TLorentzVector varMet(met);
	      varMet += jetDiff;
	      varMet += (varlp4-lp4);
	      if(varName=="UncMET") 
		{
		  TLorentzVector uncMet(met+varlp4+jetSum);
		  uncMet *= (1.0 + (isign==0?-1.:1.)*0.1);
		  varMet = uncMet-varlp4-jetSum;
		}
	      varMet.SetPz(0.); varMet.SetE(varMet.Pt());
	      // float varmt(computeMT(varlp4,varMet) );

	      if(varBJets.size()+varLightJets.size()<1) continue;

	      //ready to fill histograms
	      float mlb(varBJets.size() ? (varBJets[0]+varlp4).M() : 0.);
	      if(varBJets.size()>1) mlb=TMath::Min( (float) mlb, (float)(varBJets[1]+varlp4).M() );

	      //balancing variable
	      float varRMPF(0);
	      TVector2 metT(varMet.Px(),varMet.Py());
	      TVector2 visT(isZ ? vardilp4.Px() : varlp4.Px(), isZ ? vardilp4.Py() : varlp4.Py());
	      varRMPF=1.0+(metT*visT)/visT.Mod2();
	
	      //update categories
	      int nvarJetsCat=TMath::Min((int)(varBJets.size()+varLightJets.size()),(int)4);
	      int nvarBtagsCat=TMath::Min((int)(varBJets.size()),(int)2);
	      std::vector<TString> varcatsToFill(2,Form("%dj",nvarJetsCat));
	      varcatsToFill[1]+= Form("%dt",nvarBtagsCat);

	      for(size_t icat=0; icat<2; icat++)
                {
		  TString tag=varcatsToFill[icat];
		  all2dPlots["metptshapes_"+tag+"_exp"]->Fill(varMet.Pt(),2*ivar+isign,newWgt);
		  all2dPlots["minmlbshapes_"+tag+"_exp"]->Fill(mlb,2*ivar+isign,newWgt);
		  all2dPlots["RMPFshapes_"+tag+"_exp"]->Fill(varRMPF,2*ivar+isign,newWgt);
		  if(icat!=0) continue;
		  all2dPlots["nbtagsshapes_"+tag+"_exp"]->Fill(varBJets.size(),2*ivar+isign,newWgt);
		}
	    }
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
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : all2dPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}

//
std::map<Int_t,Float_t> lumiPerRun()
{
  std::map<Int_t,Float_t> lumiMap;
  lumiMap[ 256630 ]=   972453.859  ;
  lumiMap[ 256673 ]=    5777.795   ;
  lumiMap[ 256674 ]=   96620.611   ;
  lumiMap[ 256675 ]=  7596277.091  ;
  lumiMap[ 256676 ]=  9548596.716  ;
  lumiMap[ 256677 ]=  16055929.005 ;
  lumiMap[ 256801 ]=  9165130.512  ;
  lumiMap[ 256842 ]=   17179.759   ;
  lumiMap[ 256843 ]=  38623648.354 ;
  lumiMap[ 256866 ]=   60646.099   ;
  lumiMap[ 256867 ]=  4728164.411  ;
  lumiMap[ 256868 ]=  23371474.596 ;
  lumiMap[ 256869 ]=  1592810.667  ;
  lumiMap[ 256926 ]=  1571058.735  ;
  lumiMap[ 256941 ]=  9007928.067  ;
  lumiMap[ 257461 ]=  3194841.970  ;
  lumiMap[ 257531 ]=  8830131.818  ;
  lumiMap[ 257599 ]=  5193986.242  ;
  lumiMap[ 257613 ]=  78942504.398 ;
  lumiMap[ 257614 ]=   882267.423  ;
  lumiMap[ 257645 ]=  65301655.020 ;
  lumiMap[ 257682 ]=  13710324.756 ;
  lumiMap[ 257722 ]=   853623.834  ;
  lumiMap[ 257723 ]=  6248262.776  ;
  lumiMap[ 257735 ]=   542314.185  ;
  lumiMap[ 257751 ]=  28193282.293 ;
  lumiMap[ 257804 ]=   221266.820  ;
  lumiMap[ 257805 ]=  17835928.295 ;
  lumiMap[ 257816 ]=  25362711.952 ;
  lumiMap[ 257819 ]=  15747569.138 ;
  lumiMap[ 257968 ]=  17516821.210 ;
  lumiMap[ 257969 ]=  40746097.178 ;
  lumiMap[ 258129 ]=  6048266.862  ;
  lumiMap[ 258136 ]=  3759607.514  ;
  lumiMap[ 258157 ]=  4049174.674  ;
  lumiMap[ 258158 ]= 109407017.749 ;
  lumiMap[ 258159 ]=  26405897.556 ;
  lumiMap[ 258177 ]= 110399397.530 ;
  lumiMap[ 258211 ]=  6705307.183  ;
  lumiMap[ 258213 ]=  12097115.770 ;
  lumiMap[ 258214 ]=  15902014.828 ;
  lumiMap[ 258215 ]=   430560.536  ;
  lumiMap[ 258287 ]=  13905008.423 ;
  lumiMap[ 258403 ]=  16154078.120 ;
  lumiMap[ 258425 ]=  10663266.828 ;
  lumiMap[ 258426 ]=   786886.125  ;
  lumiMap[ 258427 ]=  8265152.264  ;
  lumiMap[ 258428 ]=  12092990.844 ;
  lumiMap[ 258432 ]=   291643.938  ;
  lumiMap[ 258434 ]=  31742537.170 ;
  lumiMap[ 258440 ]=  46151016.104 ;
  lumiMap[ 258444 ]=  2140077.077  ;
  lumiMap[ 258445 ]=  16854431.566 ;
  lumiMap[ 258446 ]=  7668549.265  ;
  lumiMap[ 258448 ]=  36534661.251 ;
  lumiMap[ 258655 ]=   400033.787  ;
  lumiMap[ 258656 ]=  26619578.792 ;
  lumiMap[ 258694 ]=  15960985.116 ;
  lumiMap[ 258702 ]=  30438445.207 ;
  lumiMap[ 258703 ]=  32463087.020 ;
  lumiMap[ 258705 ]=  7903450.408  ;
  lumiMap[ 258706 ]=  53978046.823 ;
  lumiMap[ 258712 ]=  35567148.443 ;
  lumiMap[ 258713 ]=  10456680.471 ;
  lumiMap[ 258714 ]=  4285875.372  ;
  lumiMap[ 258741 ]=  4899388.331  ;
  lumiMap[ 258742 ]=  65057022.311 ;
  lumiMap[ 258745 ]=  14597123.777 ;
  lumiMap[ 258749 ]=  47585277.273 ;
  lumiMap[ 258750 ]=  15156434.748 ;
  lumiMap[ 259626 ]=  11169319.844 ;
  lumiMap[ 259637 ]=  15450881.949 ;
  lumiMap[ 259681 ]=  1961953.523  ;
  lumiMap[ 259683 ]=  7494008.728  ;
  lumiMap[ 259685 ]=  54416576.426 ;
  lumiMap[ 259686 ]=  26080434.224 ;
  lumiMap[ 259721 ]=  12072683.502 ;
  lumiMap[ 259809 ]=  14229830.368 ;
  lumiMap[ 259810 ]=  9772545.665  ;
  lumiMap[ 259811 ]=  7373602.434  ;
  lumiMap[ 259813 ]=   733457.718  ;
  lumiMap[ 259817 ]=   353071.431  ;
  lumiMap[ 259818 ]=  13027383.967 ;
  lumiMap[ 259820 ]=   937057.379  ;
  lumiMap[ 259821 ]=  4734961.135  ;
  lumiMap[ 259822 ]=  32487481.427 ;
  lumiMap[ 259861 ]=  6369581.222  ;
  lumiMap[ 259862 ]=  44470716.726 ;
  lumiMap[ 259884 ]=  6510093.864  ;
  lumiMap[ 259890 ]=  9352026.765  ;
  lumiMap[ 259891 ]=  9243934.111  ;
  lumiMap[ 260373 ]=  10572653.190 ;
  lumiMap[ 260424 ]=  65639766.663 ;
  lumiMap[ 260425 ]=  23157369.009 ;
  lumiMap[ 260426 ]=  43093150.340 ;
  lumiMap[ 260427 ]=  15638701.311 ;
  lumiMap[ 260431 ]=  34141975.488 ;
  lumiMap[ 260532 ]=  67694389.854 ;
  lumiMap[ 260533 ]=  1084296.728  ;
  lumiMap[ 260534 ]=  1965785.974  ;
  lumiMap[ 260536 ]=  14172830.429 ;
  lumiMap[ 260538 ]=  21910218.595 ;
  lumiMap[ 260541 ]=  1783295.494  ;
  lumiMap[ 260575 ]=  1755955.214  ;
  lumiMap[ 260576 ]=  16233038.140 ;
  lumiMap[ 260577 ]=  8223890.985  ;
  lumiMap[ 260593 ]=  34830922.812 ;
  lumiMap[ 260627 ]= 171231412.138 ;
  return lumiMap;
};

// 
float getLeptonEnergyScaleUncertainty(int l_id,float l_pt,float l_eta)
{
  float unc(0.02);
  
  // electron uncertainties for 8 TeV cf. AN-14-145   
  if(abs(l_id)==11 || abs(l_id)==1100 || abs(l_id)==2111)
    {
      float par0(-2.27e-02), par1(-7.01e-02), par2(-3.71e-04);
      if (fabs(l_eta) > 0.8 && fabs(l_eta)<1.5)
        {
          par0 = -2.92e-02;
          par1 = -6.59e-02;
          par2 = -7.22e-04;
        }
      else if(fabs(l_eta)>1.5)
        {
          par0 = -2.27e-02;
          par1 = -7.01e-02;
          par2 = -3.71e-04;
        }
      unc=fabs(par0 * TMath::Exp(par1 * l_pt) + par2);
    }

  return unc;
}

//Sources :  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt)
{
  std::vector<float> res(3,1.0);

  float ptSF(1.0), ptSF_err(0.0);
  if(TMath::Abs(eta)<0.8)       { ptSF=1.061; ptSF_err = 0.023; }
  else if(TMath::Abs(eta)<1.3)  { ptSF=1.088; ptSF_err = 0.029; }
  else if(TMath::Abs(eta)<1.9)  { ptSF=1.106; ptSF_err = 0.030; }
  else if(TMath::Abs(eta)<2.5)  { ptSF=1.126; ptSF_err = 0.094; }
  else if(TMath::Abs(eta)<3.0)  { ptSF=1.343; ptSF_err = 0.123; }
  else if(TMath::Abs(eta)<3.2)  { ptSF=1.303; ptSF_err = 0.111; }
  else if(TMath::Abs(eta)<5.0)  { ptSF=1.320; ptSF_err = 0.286; }

  res[0] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF)*(pt-genjpt)))/pt;
  res[1] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF-ptSF_err)*(pt-genjpt)))/pt;
  res[2] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF+ptSF_err)*(pt-genjpt)))/pt;
  
  return res;
}
