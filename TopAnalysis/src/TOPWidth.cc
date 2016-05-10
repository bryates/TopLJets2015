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
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"


#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"

using namespace std;


//
void RunTopWidth(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts)
{

  bool isTTbar( filename.Contains("_TTJets") );

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);
  cout << "...producing " << outname << " from " << nentries << " events" << endl;


  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  if(!ev.isData)
    {
      TString puWgtUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/pileupWgts.root");
      gSystem->ExpandPathName(puWgtUrl);
      TFile *fIn=TFile::Open(puWgtUrl);
      TGraph *puData=(TGraph *)fIn->Get("pu_nom");
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
      gr->SetName("puwgts_nom");
      puWgtGr.push_back( gr );
      tmp->Delete();
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

  lepEffUrl="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CutBasedID_TightWP_76X_18Feb.txt_SF2D.root";
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
  std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
  BTagSFUtil myBTagSFUtil;
  if(!ev.isData)
    {
      BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "central") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "central") );

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

  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  std::map<TString, TH2 *> all2dPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",2,0,2);

  std::vector<TString> lfsVec = { "E", "EE", "EM", "MM", "M" }; 
  for(size_t ilfs=0; ilfs<lfsVec.size(); ilfs++)   
    { 
      TString tag(lfsVec[ilfs]);
      allPlots["nvtx_"+tag]  = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events",30,0,30);
      allPlots["njets_"+tag]  = new TH1F("njets_"+tag,";Jet multiplicity;Events",5,0,5);

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
    }
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%10000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      //select 1 good lepton
      std::vector<int> tightLeptons,looseLeptons;
      for(int il=0; il<ev.nl; il++)
	{
	  bool passTightKin(ev.l_pt[il]>30 && fabs(ev.l_eta[il])<2.1);
	  bool passLooseKin(ev.l_pt[il]>10 && fabs(ev.l_eta[il])<2.5);
	  bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);
	  float relIso(ev.l_relIso[il]);
	  bool passTightIso( ev.l_id[il]==13 ? relIso<0.15 : (ev.l_pid[il]>>1)&0x1 );
	  bool passLooseIso(  ev.l_id[il]==13 ? relIso<0.25 : true);
	  if(passTightKin && passTightId && passTightIso) tightLeptons.push_back(il);
	  else if(passLooseKin && passLooseIso)           looseLeptons.push_back(il);
	}
      
      //check if triggers have fired
      bool hasMuTrigger((ev.muTrigger & 0x3)!=0);
      bool hasEleTrigger((ev.elTrigger & 0x1)!=0);

      //
      std::vector<int> selLeptons;
      if(tightLeptons.size()==1 && looseLeptons.size()==0)
	{
	  selLeptons.push_back( tightLeptons[0] );
	}
      if(tightLeptons.size()>=2)
	{
	  selLeptons.push_back(tightLeptons[0]);
	  selLeptons.push_back(tightLeptons[1]);
	}
      if(tightLeptons.size()==1 && looseLeptons.size()>=1)
	{
	  selLeptons.push_back(tightLeptons[0]);
	  selLeptons.push_back(looseLeptons[0]);
	}
      if(selLeptons.size()==0) continue;
      TString chTag("");
      if(selLeptons.size()==1)
	{
	  if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) chTag="E";
	  if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger)  chTag="M";
	}
      if(selLeptons.size()==2)
	{
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*11 && hasEleTrigger) chTag="EE";
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==13*13 && hasMuTrigger) chTag="MM";
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*13 && (hasEleTrigger || hasMuTrigger)) chTag="EM";
	}
      if(chTag=="") continue;

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
      std::vector<TLorentzVector> bJets,lightJets;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

	  //cross clean with leptons
	  bool overlapsWithLepton(false);
	  for(size_t il=0; il<leptons.size(); il++)
	    {
	      if(jp4.DeltaR(leptons[il])>0.4) continue;
	      overlapsWithLepton=true;
	    }
	  if(overlapsWithLepton) continue;

	  //smear jet energy resolution for MC
	  float genJet_pt(0);
	  if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }

	  // re-inforce kinematics cuts
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  //b-tag
	  float csv = ev.j_csv[k];	  
	  bool isBTagged(csv>0.800);
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0);
	      if(abs(ev.j_hadflav[k])==4) 
		{ 	
		  expEff    = expBtagEff["c"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  expEff    = expBtagEff["b"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
		}
	      else
		{
		  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
		}
	      
	      //updated b-tagging decision with the data/MC scale factor
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
	    }

	  //save jet
	  if(isBTagged) bJets.push_back(jp4);
	  else          lightJets.push_back(jp4);
	}

      //at least 2 b-jets
      if(bJets.size()<2) continue;

      //event weight
      float wgt(1.0),puWeight(1.0);
      if(!ev.isData)
	{
	  //update lepton selection scale factors, if found
	  float lepTriggerSF(1.0),lepSelSF(1.0);
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

	  
	  //pileup weight
	  puWeight=puWgtGr[0]->Eval(ev.putrue);  
	  
	  //update nominal event weight
	  float norm( normH ? normH->GetBinContent(1) : 1.0);
	  wgt=lepTriggerSF*lepSelSF*puWeight*norm;
	  if(ev.ttbar_nw>0) wgt*=ev.ttbar_w[0];
	}
      
      //nominal selection control histograms
      allPlots["puwgtctr"]->Fill(0.,puWeight!=0 ? wgt/puWeight : 0.);
      allPlots["puwgtctr"]->Fill(1.,wgt);

      allPlots["nvtx_"+chTag]->Fill(ev.nvtx,wgt);
      allPlots["njets_"+chTag]->Fill(lightJets.size(),wgt);
    }
  
  //close input file
  f->Close();

  //save histos to file  
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : all2dPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}

