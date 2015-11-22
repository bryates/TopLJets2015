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
	      TH1F *normH, 
	      Bool_t isTTbar,
	      FlavourSplitting flavourSplitting,
	      GenWeightMode genWgtMode,
	      Bool_t runSysts)
{

  std::vector<TString> systs(1,"nom");

  bool isMC(filename.Contains("MC13TeV"));
  if(!isMC) runSysts=false;

  //setup pileup weighting
  TGraph *puWgtGr=0, *puUpWgtGr=0, *puDownWgtGr=0;
  if(isMC)
    {
      TString puWgtUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/pileupWgts.root");
      gSystem->ExpandPathName(puWgtUrl);
      TFile *fIn=TFile::Open(puWgtUrl);
      if(fIn){
	puWgtGr     = (TGraph *)fIn->Get("puwgts_nom");
	puDownWgtGr = (TGraph *)fIn->Get("puwgts_down");
	puUpWgtGr   = (TGraph *)fIn->Get("puwgts_up");
	fIn->Close();
      }
    }

  // setup b-tag calibration readers
  TString btagUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CSVv2.csv");
  gSystem->ExpandPathName(btagUncUrl);
  BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
  BTagCalibrationReader btagSFbReader(&btvcalib,     BTagEntry::OP_MEDIUM, "mujets", "central");
  BTagCalibrationReader btagSFbupReader(&btvcalib,   BTagEntry::OP_MEDIUM, "mujets", "up");
  BTagCalibrationReader btagSFbdownReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "down"); 
  BTagCalibrationReader btagSFlReader(&btvcalib,     BTagEntry::OP_MEDIUM, "comb", "central");
  BTagCalibrationReader btagSFlupReader(&btvcalib,   BTagEntry::OP_MEDIUM, "comb", "up");
  BTagCalibrationReader btagSFldownReader(&btvcalib, BTagEntry::OP_MEDIUM, "comb", "down"); 

  TString btagEffExpUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/expTageff.root");
  gSystem->ExpandPathName(btagEffExpUrl);
  TFile *beffIn=TFile::Open(btagEffExpUrl);
  TGraphAsymmErrors *expEff_b=(TGraphAsymmErrors *)beffIn->Get("b");
  TGraphAsymmErrors *expEff_c=(TGraphAsymmErrors *)beffIn->Get("c");
  TGraphAsymmErrors *expEff_udsg=(TGraphAsymmErrors *)beffIn->Get("udsg");
  BTagSFUtil myBTagSFUtil;


  //jES uncertainties
  TString jesSrcNames[] = {"Absolute",        "HighPtExtra",    "SinglePionECAL",  "SinglePionHCAL", "Time",
			   "RelativeJEREC1",  "RelativeJEREC2", "RelativeJERHF",
			   "RelativePtBB",    "RelativePtEC1",  "RelativePtEC2",   "RelativePtHF",  "RelativeFSR",
			   "RelativeStatEC2", "RelativeStatHF",
			   "PileUpDataMC",    "PileUpPtBB",     "PileUpPtEC",      "PileUpPtHF", "PileUpBias",
			   "FlavorPureGluon", "FlavorPureQuark","FlavorPureCharm", "FlavorPureBottom"};
  size_t nJecSrc=sizeof(jesSrcNames)/sizeof(TString);
  std::map<Int_t,JetCorrectionUncertainty*> jecUncs;
  TString jecUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/Summer15_25nsV6M3_DATA_UncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  if(runSysts)
    {
      systs.push_back("muFUp");        systs.push_back("muFDown");
      systs.push_back("muRUp");        systs.push_back("muRDown");
      systs.push_back("muRmuFUp");     systs.push_back("muRmuFDown");
      
      for(size_t i=0; i<nJecSrc; i++)
	{
	  systs.push_back(jesSrcNames[i]+"Up");     
	  systs.push_back(jesSrcNames[i]+"Down");    

	  JetCorrectorParameters *p = new JetCorrectorParameters(jecUncUrl.Data(), jesSrcNames[i].Data());
	  jecUncs[JES_Absolute+i]=new JetCorrectionUncertainty(*p);
	}
    }

  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;

  //book histograms
  std::map<TString, TH1 *> allPlots;
  std::map<TString, TH2 *> all2dPlots;

  std::map<Int_t,Float_t> lumiMap=lumiPerRun();
  for(Int_t ij=1; ij<=4; ij++)
    {
      for(Int_t itag=-1; itag<=2; itag++)
	{
	  if(itag>ij) continue;
	  TString tag(itag<0 ? Form("%dj",ij) : Form("%dj%dt",ij,itag));
	  allPlots["ratevsrun_"+tag] = new TH1F("ratevsrun_"+tag,";Run number; Events/pb",lumiMap.size(),0,lumiMap.size());
	  Int_t runCtr(0);
	  for(std::map<Int_t,Float_t>::iterator it=lumiMap.begin(); it!=lumiMap.end(); it++,runCtr++)
	    allPlots["ratevsrun_"+tag]->GetXaxis()->SetBinLabel(runCtr+1,Form("%d",it->first));
	  allPlots["lpt_"+tag]        = new TH1F("lpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
	  allPlots["lsip3d_"+tag]     = new TH1F("lsip3d_"+tag,";3d impact parameter significance;Events" ,40,0.,20.);
	  allPlots["lchiso_"+tag]     = new TH1F("lchiso_"+tag,";Charged hadron isolation [GeV];Events" ,25,0.,50.);
	  allPlots["lchreliso_"+tag]  = new TH1F("lchreliso_"+tag,";Charged hadron relative isolation;Events" ,25,0.,0.2);
	  allPlots["leta_"+tag]       = new TH1F("leta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["jpt_"+tag]        = new TH1F("jpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
	  allPlots["jeta_"+tag]       = new TH1F("jeta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["ht_"+tag]         = new TH1F("ht_"+tag,";H_{T} [GeV];Events",40,0,800);
	  allPlots["csv_"+tag]        = new TH1F("csv_"+tag,";CSV discriminator;Events",100,0,1.0);
	  allPlots["nvtx_"+tag]       = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events" ,50,0.,50.);
	  allPlots["met_"+tag]        = new TH1F("metpt_"+tag,";Missing transverse energy [GeV];Events" ,20,0.,300.);
	  allPlots["metphi_"+tag]     = new TH1F("metphi_" + tag,";MET #phi [rad];Events" ,50,-3.2,3.2);
	  allPlots["mttbar_"+tag]     = new TH1F("mttbar_"+tag,";#sqrt{#hat{s}} [GeV];Events" ,100,300.,700.);
	  allPlots["mt_"+tag]         = new TH1F("mt_"+tag,";Transverse Mass [GeV];Events" ,20,0.,200.);
	  allPlots["minmlb_"+tag]     = new TH1F("minmlb_"+tag,";min Mass(lepton,b) [GeV];Events" ,25,0.,250.);
	  allPlots["minmlblowmttbar_"+tag]   = new TH1F("minmlblowmttbar_"+tag,";min Mass(lepton,b) [GeV];Events" ,25,0.,250.);
	  allPlots["minmlbhighmttbar_"+tag]  = new TH1F("minmlbhighmttbar_"+tag,";min Mass(lepton,b) [GeV];Events" ,25,0.,250.);

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
	      all2dPlots["mtshapes_"+tag]     = new TH2F("mtshapes_"+tag,";Transverse Mass [GeV];Events" ,20,0.,200., 2*LASTSYST,0,2*LASTSYST);
	      all2dPlots["minmlbshapes_"+tag] = new TH2F("minmlbshapes_"+tag,";min Mass(lepton,b) [GeV];Events" ,25,0.,250.,2*LASTSYST,0,2*LASTSYST);
	      if(itag==-1) all2dPlots["nbtagsshapes_"+tag] = new TH2F("nbtagsshapes_"+tag,";Category;Events" ,3,0.,3.,2*LASTSYST,0,2*LASTSYST);
	      for(size_t isyst=0; isyst<LASTSYST; isyst++)
		{
		  TString baseLabel(getSystematicsLabel(isyst));
		  for(size_t ivar=0; ivar<2; ivar++)
		    {
		      TString label(baseLabel + (ivar==0 ? "Up" : "Down"));
		      all2dPlots["mtshapes_"+tag]->GetYaxis()->SetBinLabel(isyst+1, label);
		      all2dPlots["minmlbshapes_"+tag]->GetYaxis()->SetBinLabel(isyst+1,label);
		      if(itag==-1) all2dPlots["nbtagsshapes_"+tag]->GetYaxis()->SetBinLabel(isyst+1,label);
		    }
		}
	    }
	}
    }

  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  
  //read tree from file
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev);

  //loop over events
  Int_t nentries(t->GetEntriesFast());
  cout << "...producing " << outname << " from " << nentries << " events" << (runSysts ? " syst variations will be considered" : "") << endl;
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
      if((abs(ev.l_id) == 13 || abs(ev.l_id) == 1300) && ((ev.muTrigger>>0)&0x1)==0) continue;
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
	  if(jp4.DeltaR(lp4)<0.4) continue;
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  
	  if(leadingJetIdx<0) leadingJetIdx=k;
	  htsum     += jp4.Pt();
	  visSystem += jp4;

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
		  expEff        = expEff_c->Eval(jptForBtag); 
		  jetBtagSF     = btagSFbReader.eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  nbjets++;
		  expEff=expEff_b->Eval(jptForBtag); 
		  jetBtagSF=btagSFbReader.eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		}
	      else
		{
		  nljets++;
		  expEff=expEff_udsg->Eval(jptForBtag);
                  jetBtagSF=btagSFlReader.eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
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
      std::vector<float> puWeight(3,1.0),lepSF(3,1.0);
      if(!ev.isData)
	{
	  lepSF=getLeptonSelectionScaleFactor(ev.l_id,ev.l_pt,ev.l_eta,ev.isData);
	  if(puWgtGr)
	    {
	      puWeight[0]=puWgtGr->Eval(ev.putrue);  
	      puWeight[1]=puUpWgtGr->Eval(ev.putrue); 
	      puWeight[2]=puDownWgtGr->Eval(ev.putrue);
	    }
	  float norm( normH ? normH->GetBinContent(1) : 1.0);
	  wgt=norm*lepSF[0]*puWeight[0];
	  if(genWgtMode!=NOGENWGT) wgt *= ev.ttbar_w[0];
	}


      //nominal selection control histograms
      if(bJets.size()+lightJets.size()>=1)
	{
	  std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);

	  int nJetsCat=TMath::Min((int)(bJets.size()+lightJets.size()),(int)4);
	  int nBtagsCat=TMath::Min((int)(bJets.size()),(int)2);

	  std::vector<TString> catsToFill(2,Form("%dj",nJetsCat));
	  catsToFill[1]+= Form("%dt",nBtagsCat);
	  
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
	      allPlots["lchiso_"+tag]->Fill(ev.l_chargedHadronIso,wgt);
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
	      allPlots["nbtags_"+tag]->Fill(bJets.size(),wgt);

	      if(bJets.size())
		{
		  float mlb=(bJets[0]+lp4).M();
		  if(bJets.size()>1) mlb=TMath::Min( (float) mlb, (float)(bJets[1]+lp4).M() );
		  allPlots["minmlb_nom_"+tag]->Fill(mlb,wgt);
		  if(visSystem.M()<375) allPlots["minmlblowmttbar_"+tag]->Fill(mlb,wgt);
		  else                  allPlots["minmlbhighmttbar_"+tag]->Fill(mlb,wgt);
		}	  
	    }
	}


      if(!runSysts) continue;

      /*

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
TString getSystematicsLabel(int varIdx)
{
  switch(varIdx)
    {
    case NOMINAL: return "nominal"; break;
    case PU: return "pileup"; break;                 
    case MUTRIGGER: return "#mu trigger"; break;          
    case MUEFF: return "#mu sel."; break;              
    case MUSCALE: return "#mu scale"; break;            
    case ETRIGGER: return "e trigger"; break;           
    case EEFF: return "e sel."; break;               
    case ESCALE: return "e scale"; break;             
    case BEFF: return "SF_{b}"; break;               
    case CEFF: return "SF_{c}"; break;               
    case LEFF: return "SF_{udsg}"; break;               
    case UMET: return "unc. MET"; break;               
    case JER: return "JER"; break;                
    case JES_Absolute: return "Absolute"; break;
    case JES_HighPtExtra: return "HighPtExtra"; break;
    case JES_SinglePionECAL: return "SinglePionECAL"; break;
    case JES_SinglePionHCAL: return "SinglePionHCAL"; break;
    case JES_Time: return "Time"; break;
    case JES_RelativeJEREC1: return "RelativeJEREC1"; break;
    case JES_RelativeJEREC2: return "RelativeJEREC2"; break;
    case JES_RelativeJERHF: return "RelativeJERHF"; break;
    case JES_RelativePtBB: return "RelativePtBB"; break;
    case JES_RelativePtEC1: return "RelativePtEC1"; break;
    case JES_RelativePtEC2: return "RelativePtEC2"; break;
    case JES_RelativePtHF: return "RelativePtHF"; break;
    case JES_RelativeFSR: return "RelativeFSR"; break;
    case JES_RelativeStatEC2: return "RelativeStatEC2"; break;
    case JES_RelativeStatHF: return "RelativeStatHF"; break;
    case JES_PileUpDataMC: return "PileUpDataMC"; break;
    case JES_PileUpPtBB: return "PileUpPtBB"; break;
    case JES_PileUpPtEC: return "PileUpPtEC"; break;
    case JES_PileUpPtHF: return "PileUpPtHF"; break;
    case JES_PileUpBias: return "PileUpBias"; break;
    case JES_FlavorPureGluon: return "FlavorPureGluon"; break;
    case JES_FlavorPureQuark: return "FlavorPureQuark"; break;
    case JES_FlavorPureCharm: return "FlavorPureCharm"; break;
    case JES_FlavorPureBottom: return "FlavorPureBottom"; break;
    case QCD_MUF: return "#mu_{F}"; break;
    case QCD_MUR: return "#mu_{R}"; break;
    case QCD_MURMUF: return "#mu_{R}#mu_{F}"; break;
    default: return ""; break;
    }

  return "";
}

//
std::vector<float> getLeptonSelectionScaleFactor(int id,float pt,float eta,bool isData)
{
  std::vector<float> lepSelSF(3,1.0);
  if(isData) return lepSelSF;

  std::pair<float,float>res(1.0,0.0);

  //electrons
  if(abs(id)==11)
    {
      if (fabs(eta)<0.8)
	{
	  if (pt<30)      { res.first=0.927; res.second=0.073; }
	  else if (pt<40) { res.first=0.975; res.second=0.018; }
	  else if (pt<50) { res.first=0.962; res.second=0.036; }
	  else            { res.first=0.955; res.second=0.022; }
	}
      else if (fabs(eta)<1.5)
	{
	  if (pt<30)      { res.first=0.891; res.second=0.074; }
	  else if (pt<40) { res.first=0.965; res.second=0.020; }
	  else if (pt<50) { res.first=0.968; res.second=0.018; }
	  else            { res.first=0.955; res.second=0.018; }
	}
      else
	{
	  if (pt<30)      { res.first=0.956; res.second=0.059; }
	  else if (pt<40) { res.first=0.995; res.second=0.018; }
	  else if (pt<50) { res.first=0.993; res.second=0.019; }
	  else            { res.first=0.985; res.second=0.023; }
	}
    }

  //muons
  if (abs(id)==13)
    {
      if (fabs(eta)<0.9)
	{
	  if (pt<30)      { res.first=1.003; res.second=0.019; }
	  else if (pt<40) { res.first=1.014; res.second=0.015; }
	  else if (pt<50) { res.first=1.001; res.second=0.014; }
	  else            { res.first=0.983; res.second=0.014; }
	}
      else if(fabs(eta)<1.2)
	{
	  if (pt<30)      { res.first=0.993; res.second=0.019; }
	  else if (pt<40) { res.first=0.994; res.second=0.015; }
	  else if (pt<50) { res.first=0.980; res.second=0.014; }
	  else            { res.first=0.987; res.second=0.015; }
	}
      else
	{
	  if (pt<30)      { res.first=1.023; res.second=0.028; }
	  else if (pt<40) { res.first=0.994; res.second=0.014; }
	  else if (pt<50) { res.first=0.996; res.second=0.014; }
	  else            { res.first=0.979; res.second=0.014; }
	}
    }

  lepSelSF[0]=res.first;
  lepSelSF[1]=res.first+res.second;
  lepSelSF[2]=res.first-res.second;
  return lepSelSF;
}

//Sources
//  Assuming nominal JER but uncertainties from Run I
//  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt)
{
  std::vector<float> res(3,1.0);

  float ptSF(1.0), ptSF_err(0.06);
  if(TMath::Abs(eta)<0.5) 
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.012,2)+pow(0.5*(0.062+0.061),2));
    }
  else if(TMath::Abs(eta)<1.1)
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.012,2)+pow(0.5*(0.056+0.055),2));
    }
  else if(TMath::Abs(eta)<1.7)
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.017,2)+pow(0.5*(0.063+0.062),2));
    }
  else if(TMath::Abs(eta)<2.3)
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    }
  else
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.127,2)+pow(0.5*(0.155+0.153),2));
    }

  res[0] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF)*(pt-genjpt)))/pt;
  res[1] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF-ptSF_err)*(pt-genjpt)))/pt;
  res[2] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF+ptSF_err)*(pt-genjpt)))/pt;
  
  return res;
}
