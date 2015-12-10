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
	  allPlots["lreliso_"+tag]     = new TH1F("lreliso_"+tag,";Relative isolation;Events" ,25,0.,0.5);
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
	      //gen level systematics
	      if(normH)
		{
		  all2dPlots["mtshapes_"+tag+"_gen"]                   
		    = new TH2F("mtshapes_"+tag+"_gen", ";Transverse Mass [GeV];Events" ,    20,0.,200., nGenSysts,0,nGenSysts);
		  all2dPlots["minmlbshapes_"+tag+"_gen"]               
		    = new TH2F("minmlbshapes_"+tag+"_gen", ";min Mass(lepton,b) [GeV];Events" , 25,0.,250., nGenSysts,0,nGenSysts);
		  if(itag==-1) 
		    all2dPlots["nbtagsshapes_"+tag+"_gen"] 
		      = new TH2F("nbtagsshapes_"+tag+"_gen", ";Category;Events" , 3, 0.,3.,   nGenSysts,0,nGenSysts);		  
		  for(Int_t igen=0; igen<nGenSysts; igen++)
		    {
		      TString label(normH->GetXaxis()->GetBinLabel(igen+1));
		      all2dPlots["mtshapes_"+tag+"_gen"]    ->GetYaxis()->SetBinLabel(igen+1,label);
		      all2dPlots["minmlbshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,label);
		      if(itag!=-1) continue;
		      all2dPlots["nbtagsshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,label);
		    }
		}
	      
	      //experimental systematics
	      Int_t nExpSysts=expSysts.size();
	      if(nExpSysts>0)
		{
		  all2dPlots["mtshapes_"+tag+"_exp"]                  
		    = new TH2F("mtshapes_"+tag+"_exp",  ";Transverse Mass [GeV];Events" ,   20,0.,200., 2*nExpSysts,0,2*nExpSysts);
		  all2dPlots["minmlbshapes_"+tag+"_exp"]              
		    = new TH2F("minmlbshapes_"+tag+"_exp", ";min Mass(lepton,b) [GeV];Events" ,25,0.,250., 2*nExpSysts,0,2*nExpSysts);
		  if(itag==-1) 
		    all2dPlots["nbtagsshapes_"+tag+"_exp"] 
		      = new TH2F("nbtagsshapes_"+tag+"_exp", ";Category;Events" ,  3,0.,3.,    2*nExpSysts,0,2*nExpSysts);
		  for(Int_t isyst=0; isyst<nExpSysts; isyst++)
		    {
		      for(Int_t ivar=0; ivar<2; ivar++)
			{
			  TString label(expSysts[isyst] + (ivar==0 ? "Down" : "Up"));
			  all2dPlots["mtshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
			  all2dPlots["minmlbshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
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
      
      /*
      //base kinematics
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt,ev.l_eta,ev.l_phi,ev.l_mass);
      if(lp4.Pt()<30 || fabs(lp4.Eta())>2.1) continue;
      float relIsoDeltaBeta((ev.l_chargedHadronIso
			     +max(0.,ev.l_neutralHadronIso+ev.l_photonIso-0.5*ev.l_puChargedHadronIso))/ev.l_pt);

      //select according to the lepton id/charge
      if(channelSelection!=0)
	{
	  if(abs(ev.l_id)!=abs(channelSelection)) continue;
	  if(channelSelection==1300)
	    {
	      if(relIsoDeltaBeta<0.2) continue;
	      //float relchIso = ev.l_chargedHadronIso/ev.l_pt;
	      //if(relchIso<0.4) continue;
	    }
	}

      if(chargeSelection!=0 &&  ev.l_charge!=chargeSelection) continue;

      //apply trigger requirement
      if((abs(ev.l_id) == 13 || abs(ev.l_id) == 1300))
	{
	  if(ev.isData  && ((ev.muTrigger>>2)&0x1)==0) continue;
	  if(!ev.isData && ((ev.muTrigger>>0)&0x1)==0) continue;
	}
      if((abs(ev.l_id) == 11 || abs(ev.l_id) == 1100) && ((ev.elTrigger>>0)&0x1)==0) continue;
      
      //select jets
      Float_t htsum(0);
      TLorentzVector jetDiff(0,0,0,0);
      std::vector<TLorentzVector> bJets,lightJets;
      TLorentzVector visSystem(lp4);
      int nbjets(0),ncjets(0),nljets(0),leadingJetIdx(-1);
      std::vector<int> resolvedJetIdx;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
	  if(jp4.DeltaR(lp4)<0.4) continue;

	  resolvedJetIdx.push_back(k);
	  jetDiff -= jp4;

	  //smear jet energy resolution for MC
	  float genJet_pt(ev.genj_pt[k]); 
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }
	  jetDiff += jp4;

	  //cross clean with respect to leptons and re-inforce kinematics cuts

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

      //MET and transverse mass
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);
      met+=jetDiff;
      met.SetPz(0.); met.SetE(met.Pt());
      
      float mt( computeMT(lp4,met) );
      
      //compute neutrino kinematics
      neutrinoPzComputer.SetMET(met);
      neutrinoPzComputer.SetLepton(lp4);
      float nupz=neutrinoPzComputer.Calculate();
      TLorentzVector neutrinoHypP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));
      visSystem+=neutrinoHypP4;

      //event weight
      float wgt(1.0);
      std::vector<float> puWeight(3,1.0),lepTriggerSF(3,1.0),lepSelSF(3,1.0);
      if(!ev.isData)
	{
	  //update lepton selection scale factors, if found
	  TString prefix("m");
	  if(abs(ev.l_id)==11 || abs(ev.l_id)==1100) prefix="e";
	  if(lepEffH.find(prefix+"_sel")!=lepEffH.end())
	    {
	      float minEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmin() ), maxEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmax()-0.01 );
	      float etaForEff=TMath::Max(TMath::Min(fabs(ev.l_eta),maxEtaForEff),minEtaForEff);
	      Int_t etaBinForEff=lepEffH[prefix+"_sel"]->GetXaxis()->FindBin(etaForEff);
	      
	      float minPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmin() ), maxPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmax()-0.01 );
	      float ptForEff=TMath::Max(TMath::Min(fabs(ev.l_pt),maxPtForEff),minPtForEff);
	      Int_t ptBinForEff=lepEffH[prefix+"_sel"]->GetYaxis()->FindBin(ptForEff);

	      float selSF(lepEffH[prefix+"_sel"]->GetBinContent(etaBinForEff,ptBinForEff));
	      float selSFUnc(lepEffH[prefix+"_sel"]->GetBinError(etaBinForEff,ptBinForEff));

	      float trigSF(lepEffH[prefix+"_trig"]->GetBinContent(etaBinForEff,ptBinForEff));
	      float trigSFUnc(lepEffH[prefix+"_trig"]->GetBinError(etaBinForEff,ptBinForEff));

	      lepTriggerSF[0]=trigSF; lepTriggerSF[1]=trigSF-trigSFUnc; lepTriggerSF[2]=trigSF+trigSFUnc;
	      lepSelSF[0]=selSF;      lepSelSF[1]=selSF-selSFUnc;       lepSelSF[2]=selSF+selSFUnc;

	      if(trigSF<0.5 || selSF<0.5)
		{
		  cout << ev.l_id << " " << ev.l_pt << " " << ptForEff << " " << ev.l_eta << " " << etaForEff << endl;
		  cout << trigSF << " " << trigSF-trigSFUnc << " " << trigSF+trigSFUnc << endl;
		  cout << selSF << " " << selSF-selSFUnc << " " << selSF+selSFUnc << endl;
		}
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

	      allPlots["lpt_"+tag]->Fill(ev.l_pt,wgt);
	      allPlots["lsip3d_"+tag]->Fill(ev.l_ip3dsig,wgt);
	      allPlots["lreliso_"+tag]->Fill(relIsoDeltaBeta,wgt);
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
      
      //gen weighting systematics
      if(bJets.size()+lightJets.size()>=1 && normH)
	{
	  float mlb(bJets.size() ? (bJets[0]+lp4).M() : 0.);
	  if(bJets.size()>1) mlb=TMath::Min( (float) mlb, (float)(bJets[1]+lp4).M() );
	  
	  for(Int_t igen=0; igen<nGenSysts; igen++)
	    {
	      for(size_t icat=0; icat<2; icat++)
		{
		  float newWgt = wgt*(normH->GetBinContent(igen+1)*ev.ttbar_w[igen])/(normH->GetBinContent(1)*ev.ttbar_w[0]);		  
		  TString tag=catsToFill[icat];	 
		  all2dPlots["mtshapes_"+tag+"_gen"]->Fill(mt,igen,newWgt);
		  all2dPlots["minmlbshapes_"+tag+"_gen"]->Fill(mlb,igen,newWgt);
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
		 (varName=="MuTrigger" && (abs(ev.l_id)==13 || abs(ev.l_id)==1300)) ||
		 (varName=="EleTrigger" && (abs(ev.l_id)==11 || abs(ev.l_id)==1300))
		 )
		newWgt *= (isign==0 ? lepTriggerSF[1]/lepTriggerSF[0] : lepTriggerSF[2]/lepTriggerSF[0]);
	      if(
		 (varName=="MuEfficiency" && (abs(ev.l_id)==13 || abs(ev.l_id)==1300)) ||
		 (varName=="EleEfficiency" && (abs(ev.l_id)==11 || abs(ev.l_id)==1300))
		 )
		newWgt *= (isign==0 ? lepSelSF[1]/lepSelSF[0] : lepSelSF[2]/lepSelSF[0]);

	      //lepton scale systematics
	      TLorentzVector varlp4(lp4);
	      if(
		 (varName=="MuScale" && (abs(ev.l_id)==13 || abs(ev.l_id)==1300)) ||
		 (varName=="EleScale" && (abs(ev.l_id)==11 || abs(ev.l_id)==1300))
		 )
		varlp4 = (1.0+(isign==0?-1.:1.)*getLeptonEnergyScaleUncertainty(abs(ev.l_id),lp4.Pt(),lp4.Eta()) ) *lp4;
	      if(varlp4.Pt()<30 || fabs(varlp4.Eta())>2.1) continue;

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
	      float varmt(computeMT(varlp4,varMet) );

	      if(varBJets.size()+varLightJets.size()<1) continue;

	      //ready to fill histograms
	      float mlb(varBJets.size() ? (varBJets[0]+varlp4).M() : 0.);
	      if(varBJets.size()>1) mlb=TMath::Min( (float) mlb, (float)(varBJets[1]+varlp4).M() );

	      //update categories
	      int nvarJetsCat=TMath::Min((int)(varBJets.size()+varLightJets.size()),(int)4);
	      int nvarBtagsCat=TMath::Min((int)(varBJets.size()),(int)2);
	      std::vector<TString> varcatsToFill(2,Form("%dj",nvarJetsCat));
	      varcatsToFill[1]+= Form("%dt",nvarBtagsCat);

	      for(size_t icat=0; icat<2; icat++)
                {
		  TString tag=varcatsToFill[icat];
		  all2dPlots["mtshapes_"+tag+"_exp"]->Fill(varmt,2*ivar+isign,newWgt);
		  all2dPlots["minmlbshapes_"+tag+"_exp"]->Fill(mlb,2*ivar+isign,newWgt);
		  if(icat!=0) continue;
		  all2dPlots["nbtagsshapes_"+tag+"_exp"]->Fill(varBJets.size(),2*ivar+isign,newWgt);
		}
	    }
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
float getLeptonEnergyScaleUncertainty(int l_id,float l_pt,float l_eta)
{
  float unc(0.02);
  
  // electron uncertainties for 8 TeV cf. AN-14-145   
  if(abs(l_id)==11)
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
