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
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib,     BTagEntry::OP_MEDIUM, "mujets", "central") );
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib,   BTagEntry::OP_MEDIUM, "mujets", "up") );
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "down") ); 
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib,     BTagEntry::OP_MEDIUM, "comb", "central") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib,   BTagEntry::OP_MEDIUM, "comb", "up") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "comb", "down") ); 
      
      TFile *beffIn=TFile::Open(btagEffExpUrl);
      expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();
    }

  //JET ENERGY SCALE: https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Summer15_uncertainties
  TString jecUncUrl("${CMSSW_BASE}/src/UserCode/BJetEnergyPeak/data/Summer15_25nsV6M3_DATA_UncertaintySources_AK4PFchs.txt");
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

  //experimental systematics
  std::vector<TString> expSysts;
  if(runSysts)
    {
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
    }


  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  std::map<TString, TH2 *> all2dPlots;
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
	  allPlots["lsip3d_"+tag]     = new TH1F("lsip3d_"+tag,";3d impact parameter significance;Events" ,40,0.,20.);
	  allPlots["lchreliso_"+tag]  = new TH1F("lchreliso_"+tag,";Charged hadron relative isolation;Events" ,25,0.,0.2);
	  allPlots["leta_"+tag]       = new TH1F("leta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["jpt_"+tag]        = new TH1F("jpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
	  allPlots["jeta_"+tag]       = new TH1F("jeta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["ht_"+tag]         = new TH1F("ht_"+tag,";H_{T} [GeV];Events",40,0,800);
	  allPlots["csv_"+tag]        = new TH1F("csv_"+tag,";CSV discriminator;Events",100,0,1.0);
	  allPlots["nvtx_"+tag]       = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events" ,50,0.,50.);
	  allPlots["met_"+tag]        = new TH1F("metpt_"+tag,";Missing transverse energy [GeV];Events" ,20,0.,300.);
	  allPlots["metphi_"+tag]     = new TH1F("metphi_" + tag,";MET #phi [rad];Events" ,50,-3.2,3.2);
	  allPlots["mttbar_"+tag]     = new TH1F("mttbar_"+tag,";#sqrt{#hat{s}} [GeV];Events" ,50,0.,1000.);
	  allPlots["mt_"+tag]         = new TH1F("mt_"+tag,";Transverse Mass [GeV];Events" ,20,0.,200.);
	  allPlots["minmlb_"+tag]     = new TH1F("minmlb_"+tag,";min Mass(lepton,b) [GeV];Events" ,25,0.,250.);
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
	      //gen weight systematics
	      Int_t nGenSysts=normH->GetNbinsX();
	      all2dPlots["mtshapes"+tag+"_gen"]                   = new TH2F("mtshapes"+tag+"_gen",      ";Transverse Mass [GeV];Events" ,    20,0.,200., nGenSysts,0,nGenSysts);
	      all2dPlots["minmlbshapes"+tag+"_gen"]               = new TH2F("minmlbshapes"+tag+"_gen",  ";min Mass(lepton,b) [GeV];Events" , 25,0.,250., nGenSysts,0,nGenSysts);
	      if(itag==-1) all2dPlots["nbtagsshapes_"+tag+"_gen"] = new TH2F("nbtagsshapes_"+tag+"_gen", ";Category;Events" ,                 3, 0.,3.,   nGenSysts,0,nGenSysts);
	      for(Int_t igen=0; igen<nGenSysts; igen++)
		{
		  all2dPlots["mtshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,normH->GetXaxis()->GetBinLabel(igen+1));
		  all2dPlots["minmlbshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,normH->GetXaxis()->GetBinLabel(igen+1));
		  all2dPlots["nbtagsshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,normH->GetXaxis()->GetBinLabel(igen+1));
		}
	      
	      //experimental systematics
	      Int_t nExpSysts=expSysts.size();
	      all2dPlots["mtshapes_"+tag+"_exp"]                  = new TH2F("mtshapes_"+tag+"_exp",     ";Transverse Mass [GeV];Events" ,   20,0.,200., 2*nExpSysts,0,2*nExpSysts);
	      all2dPlots["minmlbshapes_"+tag+"_exp"]              = new TH2F("minmlbshapes_"+tag+"_exp", ";min Mass(lepton,b) [GeV];Events" ,25,0.,250., 2*nExpSysts,0,2*nExpSysts);
	      if(itag==-1) all2dPlots["nbtagsshapes_"+tag+"_exp"] = new TH2F("nbtagsshapes_"+tag+"_exp", ";Category;Events" ,                3,0.,3.,    2*nExpSysts,0,2*nExpSysts);
	      for(Int_t isyst=0; isyst<nExpSysts; isyst++)
		{
		  for(Int_t ivar=0; ivar<2; ivar++)
		    {
		      TString label(expSysts[isyst] + (ivar==0 ? "Up" : "Down"));
		      all2dPlots["mtshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
		      all2dPlots["minmlbshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
		      if(itag==-1) all2dPlots["nbtagsshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
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
      printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));
      
      //base kinematics
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt,ev.l_eta,ev.l_phi,ev.l_mass);
      if(lp4.Pt()<30 || fabs(lp4.Eta())>2.1) continue;

      //select according to the lepton id/charge
      if(channelSelection!=0)
	{
	  if(abs(ev.l_id)!=abs(channelSelection)) continue;
	  if(channelSelection==1300)
	    {
	      float relchIso = ev.l_chargedHadronIso/ev.l_pt;
	      if(relchIso<0.4) continue;
	    }
	}

      if(chargeSelection!=0 &&  ev.l_charge!=chargeSelection) continue;

      //apply trigger requirement
      if((abs(ev.l_id) == 13 || abs(ev.l_id) == 1300))
	{
	  if(ev.isData  && ((ev.muTrigger>>0)&0x1)==0) continue;
	  if(!ev.isData && ((ev.muTrigger>>2)&0x1)==0) continue;
	}
      if((abs(ev.l_id) == 11 || abs(ev.l_id) == 1100) && ((ev.elTrigger>>0)&0x1)==0) continue;

      //MET and transverse mass
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);
      float mt( computeMT(lp4,met) );
      
      //compute neutrino kinematics
      neutrinoPzComputer.SetMET(met);
      neutrinoPzComputer.SetLepton(lp4);
      float nupz=neutrinoPzComputer.Calculate();
      TLorentzVector neutrinoHypP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));
      
      //select jets
      Float_t htsum(0);
      std::vector<TLorentzVector> bJets,lightJets;
      TLorentzVector visSystem(lp4+neutrinoHypP4);
      int nbjets(0),ncjets(0),nljets(0),leadingJetIdx(-1);
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

	  //smear jet energy resolution for MC
	  float genJet_pt(ev.genj_pt[k]); 
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }

	  //cross clean with respect to leptons and re-inforce kinematics cuts
	  if(jp4.DeltaR(lp4)<0.4) continue;
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  if(leadingJetIdx<0) leadingJetIdx=k;
	  htsum     += jp4.Pt();
	  if(bJets.size()+lightJets.size()<=4) visSystem += jp4;

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
		  expEff        = expBtagEff["c"]->Eval(jptForBtag); 
		  jetBtagSF     = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  nbjets++;
		  expEff=expBtagEff["b"]->Eval(jptForBtag); 
		  jetBtagSF=sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		}
	      else
		{
		  nljets++;
		  expEff=expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF=sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
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
      

      //event weight
      float wgt(1.0);
      std::vector<float> puWeight(3,1.0),lepTriggerSF(3,1.0),lepSelSF(3,1.0);
      if(!ev.isData)
	{
	  lepTriggerSF = getLeptonTriggerScaleFactor  (ev.l_id,ev.l_pt,ev.l_eta,ev.isData);
	  lepSelSF     = getLeptonSelectionScaleFactor(ev.l_id,ev.l_pt,ev.l_eta,ev.isData);
	  if(puWgtGr.size())
	    {
	      puWeight[0]=puWgtGr[0]->Eval(ev.putrue);  
	      puWeight[1]=puWgtGr[1]->Eval(ev.putrue); 
	      puWeight[2]=puWgtGr[2]->Eval(ev.putrue);
	    }
	  float norm( normH ? normH->GetBinContent(1) : 1.0);
	  wgt=lepTriggerSF[0]*lepSelSF[0]*puWeight[0]*norm;
	  if(ev.ttbar_nw>0) wgt*=ev.ttbar_w[0];
	}

      //nominal selection control histograms
      if(bJets.size()+lightJets.size()>=1)
	{
	  int nJetsCat=TMath::Min((int)(bJets.size()+lightJets.size()),(int)4);
	  int nBtagsCat=TMath::Min((int)(bJets.size()),(int)2);
	  std::vector<TString> catsToFill(2,Form("%dj",nJetsCat));
	  catsToFill[1]+= Form("%dt",nBtagsCat);

	  std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);	  

	  for(size_t icat=0; icat<2; icat++)
	    {
	      TString tag=catsToFill[icat];
	      if(rIt!=lumiMap.end()) 
		{
		  Int_t runCtr=std::distance(lumiMap.begin(),rIt);
		  allPlots["ratevsrun_"+tag]->Fill(runCtr,1.e+6/rIt->second);
		}

	      allPlots["lpt_"+tag]->Fill(ev.l_pt,wgt);
	      allPlots["lsip3d_"+tag]->Fill(ev.l_ip3dsig,wgt);
	      allPlots["lchreliso_"+tag]->Fill(ev.l_chargedHadronIso/ev.l_pt,wgt);
	      allPlots["leta_"+tag]->Fill(ev.l_eta,wgt);
	      allPlots["jpt_"+tag]->Fill(ev.j_pt[ leadingJetIdx ],wgt);
	      allPlots["jeta_"+tag]->Fill(fabs(ev.j_eta[ leadingJetIdx ]),wgt);
	      allPlots["csv_"+tag]->Fill(ev.j_csv[ leadingJetIdx ],wgt);
	      allPlots["nvtx_"+tag]->Fill(ev.nvtx,wgt);
	      allPlots["met_"+tag]->Fill(ev.met_pt,wgt);
	      allPlots["metphi_"+tag]->Fill(ev.met_phi,wgt);
	      allPlots["mt_"+tag]->Fill(mt,wgt);
	      allPlots["mttbar_"+tag]->Fill(visSystem.M(),wgt);
	      allPlots["ht_"+tag]->Fill(htsum,wgt);
	      if(icat==0) allPlots["nbtags_"+tag]->Fill(bJets.size(),wgt);
	      
	      if(bJets.size())
		{
		  float mlb=(bJets[0]+lp4).M();
		  if(bJets.size()>1) mlb=TMath::Min( (float) mlb, (float)(bJets[1]+lp4).M() );
		  allPlots["minmlb_"+tag]->Fill(mlb,wgt);
		}	  
	    }
	}


      if(!runSysts) continue;

      /*
	// work in progress
		  allPlots["minmlb_qcdScaleDown_"+tag]->Fill(minMlb,wgtQCDScaleLo);
		  allPlots["minmlb_qcdScaleUp_"+tag]->Fill(minMlb,wgtQCDScaleHi);
		  allPlots["minmlb_puUp_"+tag]->Fill(minMlb,wgtPuUp);
		  allPlots["minmlb_puDown_"+tag]->Fill(minMlb,wgtPuDown);
		  allPlots["minmlb_muEffUp_"+tag]->Fill(minMlb,wgtMuEffUp);
		  allPlots["minmlb_muEffDown_"+tag]->Fill(minMlb,wgtMuEffDown);
		  allPlots["minmlb_eEffUp_"+tag]->Fill(minMlb,wgtElEffUp);
		  allPlots["minmlb_eEffDown_"+tag]->Fill(minMlb,wgtElEffDown);
		  allPlots["minmlb_umetUp_"+tag]->Fill(minMlb,wgt);
		  allPlots["minmlb_umetDown_"+tag]->Fill(minMlb,wgt);


	      allPlots["mt_qcdScaleDown_"+tag]->Fill(mt,wgtQCDScaleLo);
	      allPlots["mt_qcdScaleUp_"+tag]->Fill(mt,wgtQCDScaleHi);
	      allPlots["mt_puUp_"+tag]->Fill(mt,wgtPuUp);
	      allPlots["mt_puDown_"+tag]->Fill(mt,wgtPuDown);
	      allPlots["mt_muEffUp_"+tag]->Fill(mt,wgtMuEffUp);
	      allPlots["mt_muEffDown_"+tag]->Fill(mt,wgtMuEffDown);
	      allPlots["mt_eEffUp_"+tag]->Fill(mt,wgtElEffUp);
	      allPlots["mt_eEffDown_"+tag]->Fill(mt,wgtElEffDown);	      
	      allPlots["mt_umetUp_"+tag]->Fill(mtUMetup,wgtElEffUp);
	      allPlots["mt_umetDown_"+tag]->Fill(mtUMetdown,wgtElEffDown);	      


	  allPlots["njetsnbtags_nom"]->Fill(binToFill,wgt);
	  allPlots["njetsnbtags_qcdScaleDown"]->Fill(binToFill,wgtQCDScaleLo);
	  allPlots["njetsnbtags_qcdScaleUp"]->Fill(binToFill,wgtQCDScaleHi);
	  allPlots["njetsnbtags_puUp"]->Fill(binToFill,wgtPuUp);
	  allPlots["njetsnbtags_puDown"]->Fill(binToFill,wgtPuDown);
	  allPlots["njetsnbtags_muEffUp"]->Fill(binToFill,wgtMuEffUp);
	  allPlots["njetsnbtags_muEffDown"]->Fill(binToFill,wgtMuEffDown);
	  allPlots["njetsnbtags_umetUp"]->Fill(binToFill,wgt);
	  allPlots["njetsnbtags_umetDown"]->Fill(binToFill,wgt);


      float wgtPuUp      (norm*lepSF[0]                                *puWeight[1]);
      float wgtPuDown    (norm*lepSF[0]                                *puWeight[2]);
      float wgtMuEffUp   (norm*(abs(ev.l_id)==13 ? lepSF[1] : lepSF[0])*puWeight[0]);
      float wgtMuEffDown (norm*(abs(ev.l_id)==13 ? lepSF[2] : lepSF[0])*puWeight[0]);
      float wgtElEffUp   (norm*(abs(ev.l_id)==11 ? lepSF[1] : lepSF[0])*puWeight[0]);
      float wgtElEffDown (norm*(abs(ev.l_id)==11 ? lepSF[2] : lepSF[0])*puWeight[0]);
      float wgtQCDScaleLo(wgt),wgtQCDScaleHi(wgt);
      if(genWgtMode!=NOGENWGT && !ev.isData) 
	{
	  wgt           *= ev.ttbar_w[0];
	  wgtPuUp       *= ev.ttbar_w[0];
	  wgtPuDown     *= ev.ttbar_w[0];
	  wgtMuEffUp    *= ev.ttbar_w[0];
	  wgtMuEffDown  *= ev.ttbar_w[0];
	  wgtElEffUp    *= ev.ttbar_w[0];
	  wgtElEffDown  *= ev.ttbar_w[0];
	  wgtQCDScaleLo *= ev.ttbar_w[0];
	  wgtQCDScaleHi *= ev.ttbar_w[0];
	}
      if(isTTbar)
	{
	  wgtQCDScaleLo   = wgt*(normH->GetBinContent(10)/norm)*(ev.ttbar_w[9]/ev.ttbar_w[0]);
	  wgtQCDScaleHi   = wgt*(normH->GetBinContent(6)/norm)*(ev.ttbar_w[5]/ev.ttbar_w[0]);	 
	}


      TLorentzVector metJESup( met+(jetSum-jetSumJESup)),metJESdown(met+(jetSum-jetSumJESdown));
      TLorentzVector metJERup( met+(jetSum-jetSumJERup)),metJERdown(met+(jetSum-jetSumJERdown));
      TLorentzVector metUMetdown(0.9*met-0.1*(jetSum+lp4)),metUMetup(1.1*met+0.1*(jetSum+lp4));
      

      float mtJESup( computeMT(lp4,metJESup) ), mtJESdown( computeMT(lp4,metJESdown) );
      float mtJERup( computeMT(lp4,metJESup) ), mtJERdown( computeMT(lp4,metJERdown) );
      float mtUMetup( computeMT(lp4,metUMetup) ), mtUMetdown( computeMT(lp4,metUMetdown) );


	  //jet energy resolution
	  std::vector<float> jerScale(3,1.0);
	  float genJet_pt(ev.genj_pt[k]);
	  if(!ev.isData && genJet_pt>0) jerScale=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt);

	  //FIXME
	  //jet energy scale variations
	  //jecUnc->setJetEta(fabs(jp4.Eta()));
	  //jecUnc->setJetPt(jp4.Pt());
	  double unc = 0; //jecUnc->getUncertainty(true);    

	  
	  //readout the b-tagging scale factors for this jet
	  float csv = ev.j_csv[k];	  
	  bool isBTagged(csv>0.890),isBTaggedUp(isBTagged),isBTaggedDown(isBTagged);
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0), jetBtagSFUp(1.0), jetBtagSFDown(1.0);
	      if(abs(ev.j_hadflav[k])==4) 
		{ 
		  expEff        = expEff_c->Eval(jptForBtag); 
		  jetBtagSF     = btagSFbReader.eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSFUp   = btagSFbupReader.eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSFDown = btagSFbdownReader.eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  expEff=expEff_b->Eval(jptForBtag); 
		  jetBtagSF=btagSFbReader.eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSFUp=btagSFbupReader.eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSFDown=btagSFbdownReader.eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		}
	      else
		{
		  expEff=expEff_udsg->Eval(jptForBtag);
                  jetBtagSF=btagSFlReader.eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
                  jetBtagSFUp=btagSFlupReader.eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
                  jetBtagSFDown=btagSFldownReader.eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		}
	      
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
	      myBTagSFUtil.modifyBTagsWithSF(isBTaggedUp,  jetBtagSFUp,   expEff);
	      myBTagSFUtil.modifyBTagsWithSF(isBTaggedDown,jetBtagSFDown, expEff);
	    }
	  
	  //select the jet
	  if(jp4.Pt()>30)          
	    { 
	      selJetsIdx.push_back(k);
	      nJets++;      
	      jetSum += jp4;       
	      htsum += jp4.Pt();
	      if(abs(ev.j_hadflav[k])==4)      ncJets++;
	      else if(abs(ev.j_hadflav[k])==5) nbJets++;
	      else nudsgJets++;
	      nBtags += isBTagged;
	      if(!ev.isData)
		{
		  if(abs(ev.j_hadflav[k])==4)
		    {
		      nBtagsBeffLo   += isBTagged;        nBtagsBeffHi   += isBTagged;
		      nBtagsCeffLo   += isBTaggedDown;    nBtagsCeffHi   += isBTaggedUp;
		      nBtagsMistagLo += isBTagged;        nBtagsMistagHi += isBTagged;
		    }
		  else if(abs(ev.j_hadflav[k])==5)
		    {
		      nBtagsBeffLo += isBTaggedDown;	  nBtagsBeffHi += isBTaggedUp;
		      nBtagsCeffLo += isBTagged;	  nBtagsCeffHi += isBTagged;
		      nBtagsMistagLo += isBTagged;	  nBtagsMistagHi += isBTagged;
		    }
		  else
		    {
		      nBtagsBeffLo += isBTagged;           nBtagsBeffHi += isBTagged;
		      nBtagsCeffLo += isBTagged;           nBtagsCeffHi += isBTagged;
		      nBtagsMistagLo += isBTaggedDown;     nBtagsMistagHi += isBTaggedUp;
		    }
		}

	      if(isBTagged) minMlb = TMath::Min(minMlb,(Float_t)(jp4+lp4).M()); 
	      if(abs(ev.j_hadflav[k])==4 || abs(ev.j_hadflav[k])==5)
		{
		  if(isBTaggedDown) minMlbBeffLo=TMath::Min(minMlbBeffLo,(Float_t)(jp4+lp4).M());
		  if(isBTaggedUp)   minMlbBeffHi=TMath::Min(minMlbBeffHi,(Float_t)(jp4+lp4).M());
		} 
	      else
		{
		  if(isBTaggedDown) minMlbLeffLo=TMath::Min(minMlbLeffLo,(Float_t)(jp4+lp4).M());
		  if(isBTaggedUp)   minMlbLeffHi=TMath::Min(minMlbLeffHi,(Float_t)(jp4+lp4).M());
		}
	    }
	  if((jp4.Pt())*(1+unc)>30) 
	    { 
	      nJetsJESHi++;
	      jetSumJESup += (1+unc)*jp4;  
	      if(isBTagged) minMlbJESup   = TMath::Min(minMlbJESup,(Float_t)((1+unc)*jp4+lp4).M());
	    }
	  if((jp4.Pt())*(1-unc)>30) 
	    { 
	      nJetsJESLo++; 
	      jetSumJESdown += (1-unc)*jp4;
	      if(isBTagged) minMlbJESdown = TMath::Min(minMlbJESdown,(Float_t)((1-unc)*jp4+lp4).M());
	    } 
	  if(jerScale[1]*jp4.Pt()>30) 
	    { 
	      nJetsJERLo++; 
	      jetSumJERdown += jerScale[1]*jp4;
	      if(isBTagged)  minMlbJERdown = TMath::Min(minMlbJERdown,(Float_t)( jerScale[1]*jp4+lp4).M()); 
	    }
	  if(jerScale[2]*jp4.Pt()>30) 
	    { 
	      nJetsJERHi++; 
	      jetSumJERup += jerScale[2]*jp4;
	      if(isBTagged)  minMlbJERup   = TMath::Min(minMlbJERup,(Float_t)( jerScale[2]*jp4+lp4).M());
	    }
	}

      if(nJetsJESHi>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJESHi,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=0;//getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_jesUp"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots["mt_jesUp_"+tag]->Fill(mtJESup,wgt);
	  allPlots["minmlb_jesUp_"+tag]->Fill(minMlbJESup,wgt);
	}
      if(nJetsJESLo>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJESLo,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=0;//getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_jesDown"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots["mt_jesDown_"+tag]->Fill(mtJESdown,wgt);
	  allPlots["minmlb_jesDown_"+tag]->Fill(minMlbJESdown,wgt);
	}
      if(nJetsJERHi>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJERHi,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=0;//getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_jerUp"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots["mt_jerUp_"+tag]->Fill(mtJERup,wgt);
	  allPlots["minmlb_jerUp_"+tag]->Fill(minMlbJERup,wgt);
	}
      if(nJetsJERLo>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJERLo,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=0;//getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_jerDown"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots["mt_jerDown_"+tag]->Fill(mtJERdown,wgt);
	  allPlots["minmlb_jerDown_"+tag]->Fill(minMlbJERdown,wgt);
	}
      if(nJets>=1)
	{
	  int nJetsCat=TMath::Min((int)nJets,(int)4);
          
	  int nBtagsCat=TMath::Min((int)nBtagsBeffLo,(int)2);
          int binToFill=0;//getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_beffDown"]->Fill(binToFill,wgt); 
	  
	  nBtagsCat=TMath::Min((int)nBtagsBeffHi,(int)2);
	  binToFill=0;//getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots["njetsnbtags_beffUp"]->Fill(binToFill,wgt); 

	  nBtagsCat=TMath::Min((int)nBtagsMistagLo,(int)2);
	  binToFill=0;//getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots["njetsnbtags_mistagDown"]->Fill(binToFill,wgt); 

	  nBtagsCat=TMath::Min((int)nBtagsMistagHi,(int)2);
	  binToFill=0;//getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots["njetsnbtags_mistagUp"]->Fill(binToFill,wgt); 

	  cout << nBtagsCat << " " << nJetsCat << endl;
	}
      */

    }
  
  //close input file
  f->Close();

  //save histos to file  
  TString selPrefix("");  
  if(flavourSplitting!=NOFLAVOURSPLITTING) selPrefix=Form("%d_",flavourSplitting);
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+selPrefix+baseName,"RECREATE");
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}

//
std::map<Int_t,Float_t> lumiPerRun()
{
  std::map<Int_t,Float_t> lumiMap;
  lumiMap[256630]=948417.609;
  lumiMap[256673]=5534.408;
  lumiMap[256674]=92567.485;
  lumiMap[256675]=7282554.291;
  lumiMap[256676]=9172872.886;
  lumiMap[256677]=15581756.928;
  lumiMap[256801]=8830347.675;
  lumiMap[256842]=16510.582;
  lumiMap[256843]=37367788.087;
  lumiMap[256866]=58250.406;
  lumiMap[256867]=4546508.033;
  lumiMap[256868]=22542014.201;
  lumiMap[256869]=1539580.832;
  lumiMap[256926]=1499855.808;
  lumiMap[256941]=8741029.381;
  lumiMap[257461]=3057928.782;
  lumiMap[257531]=8418598.194;
  lumiMap[257599]=4940876.751;
  lumiMap[257613]=75639153.224;
  lumiMap[257614]=850778.922;
  lumiMap[257645]=62520503.888;
  lumiMap[257682]=13053256.987;
  lumiMap[257722]=810552.138;
  lumiMap[257723]=5941442.106;
  lumiMap[257735]=521278.124;
  lumiMap[257751]=27029514.967;
  lumiMap[257804]=210956.374;
  lumiMap[257805]=17038078.687;
  lumiMap[257816]=24328019.178;
  lumiMap[257819]=15147148.510;
  lumiMap[257968]=16769109.914;
  lumiMap[257969]=39179793.996;
  lumiMap[258129]=5813530.480;
  lumiMap[258136]=3617731.160;
  lumiMap[258157]=3866329.715;
  lumiMap[258158]=105692715.515;
  lumiMap[258159]=25632955.799;
  lumiMap[258177]=101938042.657;
  lumiMap[258211]=6371404.543;
  lumiMap[258213]=11525115.399;
  lumiMap[258214]=15172855.551;
  lumiMap[258215]=411364.505;
  lumiMap[258287]=13271542.328;
  lumiMap[258403]=15612888.766;
  lumiMap[258425]=10170405.992;
  lumiMap[258426]=751812.067;
  lumiMap[258427]=7901302.746;
  lumiMap[258428]=11578208.046;
  lumiMap[258432]=279944.749;
  lumiMap[258434]=30536738.787;
  lumiMap[258440]=44624917.855;
  lumiMap[258444]=2079367.112;
  lumiMap[258445]=16403375.902;
  lumiMap[258446]=7474593.995;
  lumiMap[258448]=35705866.510;
  lumiMap[258655]=383658.834;
  lumiMap[258656]=25581933.798;
  lumiMap[258694]=15219679.888;
  lumiMap[258702]=29093248.654;
  lumiMap[258703]=31138680.065;
  lumiMap[258705]=7604760.367;
  lumiMap[258706]=52122692.407;
  lumiMap[258712]=34495799.123;
  lumiMap[258713]=10164347.291;
  lumiMap[258714]=4168945.356;
  lumiMap[258741]=4446752.908;
  lumiMap[258742]=59299810.293;
  lumiMap[258745]=20966777.757;
  lumiMap[258749]=44752319.500;
  lumiMap[258750]=14330984.460;
  lumiMap[259626]=10597934.251;
  lumiMap[259637]=14838753.649;
  lumiMap[259681]=1866936.539;
  lumiMap[259683]=7150122.168;
  lumiMap[259685]=52046691.339;
  lumiMap[259686]=25072787.031;
  lumiMap[259721]=11233249.877;
  lumiMap[259809]=13415540.397;
  lumiMap[259810]=9213536.606;
  lumiMap[259811]=6958703.081;
  lumiMap[259813]=696989.195;
  lumiMap[259817]=338781.444;
  lumiMap[259818]=12345301.694;
  lumiMap[259820]=11764360.572;
  lumiMap[259821]=13847071.055;
  lumiMap[259822]=31229541.823;
  lumiMap[259861]=6085783.718;
  lumiMap[259862]=42639175.757;
  lumiMap[259884]=6228050.692;
  lumiMap[259890]=8940297.429;
  lumiMap[259891]=8847861.645;
  lumiMap[260373]=10082072.543;
  lumiMap[260424]=63075300.274;
  lumiMap[260425]=22346537.163;
  lumiMap[260426]=41726987.678;
  lumiMap[260427]=15191168.895;
  lumiMap[260431]=33250906.885;
  lumiMap[260532]=61132427.407;
  lumiMap[260533]=1059991.553;
  lumiMap[260534]=28416797.584;
  lumiMap[260536]=13137987.544;
  lumiMap[260538]=20326040.651;
  lumiMap[260541]=1663588.353;
  lumiMap[260575]=1681737.625;
  lumiMap[260576]=15468612.251;
  lumiMap[260577]=7858938.175;
  lumiMap[260593]=33318468.559;
  lumiMap[260627]=164681593.193;
  return lumiMap;
};

//
std::vector<float> getLeptonTriggerScaleFactor(int id,float pt,float eta,bool isData)
{
  std::vector<float> lepTrigSF(3,1.0);
  return lepTrigSF;
}

//
std::vector<float> getLeptonSelectionScaleFactor(int id,float pt,float eta,bool isData)
{
  std::vector<float> lepSelSF(3,1.0);
  return lepSelSF;
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
