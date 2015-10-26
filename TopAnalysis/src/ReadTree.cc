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

#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"

using namespace std;

Int_t getBtagCatBinToFill(Int_t nBtags,Int_t nJets)
{
  Int_t nJetsBin(nJets>4 ? 4 : nJets);

  Int_t binToFill(nBtags>=2?2:nBtags);
  binToFill+=3*(nJetsBin-1);
  return binToFill;
}

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
	      TGraph *puWgtGr,TGraph *puUpWgtGr,TGraph *puDownWgtGr)
{

  TString systs[]={"nom",
		   "puUp","puDown",
		   "muEffUp","muEffDown",
		   "eEffUp","eEffDown",
		   "qcdScaleDown","qcdScaleUp",
		   "umetDown", "umetUp",
		   "jesDown","jesUp",
		   "jerDown","jerUp",
		   "beffDown","beffUp",
		   "mistagDown","mistagUp"};



  //book histograms
  std::map<TString, TH1 *> allPlots;

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
	  allPlots["lpt_"+tag]  = new TH1F("lpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
	  allPlots["lsip3d_"+tag]  = new TH1F("lsip3d_"+tag,";3d impact parameter significance;Events" ,40,0.,20.);
	  allPlots["lchiso_"+tag]  = new TH1F("lchiso_"+tag,";Charged hadron isolation [GeV];Events" ,25,0.,50.);
	  allPlots["lchreliso_"+tag]  = new TH1F("lchreliso_"+tag,";Charged hadron relative isolation;Events" ,25,0.,0.2);
	  allPlots["leta_"+tag] = new TH1F("leta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["jpt_"+tag]  = new TH1F("jpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
	  allPlots["jeta_"+tag] = new TH1F("jeta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["ht_"+tag]   = new TH1F("ht_"+tag,";H_{T} [GeV];Events",40,0,800);
	  allPlots["csv_"+tag]  = new TH1F("csv_"+tag,";CSV discriminator;Events",100,0,1.0);
	  allPlots["nvtx_"+tag] = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events" ,50,0.,50.);
	  allPlots["met_"+tag]  = new TH1F("metpt_"+tag,";Missing transverse energy [GeV];Events" ,20,0.,300.);
	  allPlots["metphi_"+tag] = new TH1F("metphi_" + tag,";MET #phi [rad];Events" ,50,-3.2,3.2);
	  for(size_t isyst=0; isyst<sizeof(systs)/sizeof(TString); isyst++)
	    {
	      allPlots["mt_"+systs[isyst]+"_"+tag]     = new TH1F("mt_"+systs[isyst]+"_"+tag,";Transverse Mass [GeV];Events" ,20,0.,200.);
	      allPlots["minmlb_"+systs[isyst]+"_"+tag] = new TH1F("minmlb_"+systs[isyst]+"_"+tag,";min Mass(lepton,b) [GeV];Events" ,25,0.,250.);
	    }
	}
    }

  //category counting
  for(size_t i=0; i<sizeof(systs)/sizeof(TString); i++)
    {
      allPlots["njetsnbtags_"+systs[i]] = new TH1F("njetsnbtags_"+systs[i],";Category;Events" ,12,0.,12.);
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(1, "1j,=0b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(2, "1j,=1b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(3, "");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(4, "2j,=0b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(5, "2j,=1b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(6, "2j,#geq2b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(7, "3j,=0b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(8, "3j,=1b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(9, "3j,#geq2b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(10,"4j,=0b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(11,"4j,=1b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(12,"4j,#geq2b");
    }

  for (auto& it : allPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }


  //jet uncertainty parameterization
  TString jecUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/Summer15_25nsV5_DATA_Uncertainty_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(jecUncUrl.Data());
  
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

  //read tree from file
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev);

  //loop over events
  Int_t nentries(t->GetEntriesFast());
  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  for (Int_t i=0;i<nentries;i++)
    {
      t->GetEntry(i);
      printf ("\r [%3.0f/100] done",100.*(float)(i)/(float)(nentries));

      //base kinematics
      if(ev.l_pt<30 || fabs(ev.l_eta)>2.1) continue;

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
      
      //lepton kinematics
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt,ev.l_eta,ev.l_phi,ev.l_mass);

      //select jets
      uint32_t nJets(0),  nJetsJESLo(0),   nJetsJESHi(0),   nJetsJERLo(0),     nJetsJERHi(0);
      TLorentzVector jetSum(0,0,0,0), jetSumJESup(jetSum), jetSumJESdown(jetSum), jetSumJERup(jetSum), jetSumJERdown(jetSum);
      Float_t htsum(0);
      Int_t nudsgJets(0),ncJets(0), nbJets(0);
      uint32_t nBtags(0), nBtagsBeffLo(0), nBtagsBeffHi(0), nBtagsMistagLo(0), nBtagsMistagHi(0);
      Float_t minMlb(99999.);
      Float_t minMlbJESup(minMlb),  minMlbJESdown(minMlb);
      Float_t minMlbJERup(minMlb),  minMlbJERdown(minMlb);
      Float_t minMlbBeffHi(minMlb), minMlbBeffLo(minMlb);
      Float_t minMlbLeffHi(minMlb), minMlbLeffLo(minMlb);
      std::vector<int> selJetsIdx;           
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
	  float csv = ev.j_csv[k];	  

	  if(fabs(jp4.Eta()) > 2.4) continue;

	  //jet energy scale variations
	  jecUnc->setJetEta(fabs(jp4.Eta()));
	  jecUnc->setJetPt(jp4.Pt());
	  double unc = jecUnc->getUncertainty(true);    
	  
	  //jet energy resolution
	  std::vector<float> jerScale(3,1.0);
	  float genJet_pt(ev.genj_pt[k]);
	  if(!ev.isData && genJet_pt>0) jerScale=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt);
	  
	  //readout the b-tagging scale factors for this jet
	  bool isBTagged(csv>0.890),isBTaggedUp(isBTagged),isBTaggedDown(isBTagged);
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0), jetBtagSFUp(1.0), jetBtagSFDown(1.0);
	      if(abs(ev.j_hadflav[k])==4) 
		{ 
		  expEff=expEff_c->Eval(jptForBtag); 
		  jetBtagSF=btagSFbReader.eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSFUp=btagSFbupReader.eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSFDown=btagSFbdownReader.eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
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
	      if(ev.isData)
		{
		  nBtagsBeffLo += isBTagged;
		  nBtagsBeffHi += isBTagged;
		  nBtagsMistagLo += isBTagged;
		  nBtagsMistagHi += isBTagged;
		}
	      else
		{
		  if(abs(ev.j_hadflav[k])==4|| abs(ev.j_hadflav[k])==5)
		    {
		      nBtagsBeffLo += isBTaggedDown;
		      nBtagsBeffHi += isBTaggedUp;
		      nBtagsMistagLo += isBTagged;
		      nBtagsMistagHi += isBTagged;
		    }
		  else
		    {
		      nBtagsBeffLo += isBTagged;
		      nBtagsBeffHi += isBTagged;
		      nBtagsMistagLo += isBTaggedDown;
		      nBtagsMistagHi += isBTaggedUp;
		    }
		}
	      if(isBTagged) minMlb        = TMath::Min(minMlb,(Float_t)(jp4+lp4).M()); 
	      if(abs(ev.j_hadflav[k])==4 || abs(ev.j_hadflav[k])==5)
		{
		  if(isBTaggedDown) minMlbBeffLo=TMath::Min(minMlbBeffLo,(Float_t)(jp4+lp4).M());
		  if(isBTaggedUp) minMlbBeffHi=TMath::Min(minMlbBeffHi,(Float_t)(jp4+lp4).M());
		} 
	      else
		{
		  if(isBTaggedDown) minMlbLeffLo=TMath::Min(minMlbLeffLo,(Float_t)(jp4+lp4).M());
		  if(isBTaggedUp) minMlbLeffHi=TMath::Min(minMlbLeffHi,(Float_t)(jp4+lp4).M());
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

      //varied MET
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);
      TLorentzVector metJESup( met+(jetSum-jetSumJESup)),metJESdown(met+(jetSum-jetSumJESdown));
      TLorentzVector metJERup( met+(jetSum-jetSumJERup)),metJERdown(met+(jetSum-jetSumJERdown));
      TLorentzVector metUMetdown(0.9*met-0.1*(jetSum+lp4)),metUMetup(1.1*met+0.1*(jetSum+lp4));
      
      //varied MT
      float mt( computeMT(lp4,met) );
      float mtJESup( computeMT(lp4,metJESup) ), mtJESdown( computeMT(lp4,metJESdown) );
      float mtJERup( computeMT(lp4,metJESup) ), mtJERdown( computeMT(lp4,metJERdown) );
      float mtUMetup( computeMT(lp4,metUMetup) ), mtUMetdown( computeMT(lp4,metUMetdown) );
      
      //check if flavour splitting was required
      if(flavourSplitting!=FlavourSplitting::NOFLAVOURSPLITTING)
	{
	  if(flavourSplitting==FlavourSplitting::BSPLITTING)         { if(nbJets==0)    continue; }
	  else if(flavourSplitting==FlavourSplitting::CSPLITTING)    { if(ncJets==0 || nbJets!=0)    continue; }
	  else if(flavourSplitting==FlavourSplitting::UDSGSPLITTING) { if(nudsgJets==0 || ncJets!=0 || nbJets!=0) continue; }
	}
      
      //generator level weights to apply
      std::vector<float> lepSF=getLeptonSelectionScaleFactor(ev.l_id,ev.l_pt,ev.l_eta,ev.isData);
      std::vector<float> puWeight(3,1.0);
      if(!ev.isData && puWgtGr)
	{
	  puWeight[0]=puWgtGr->Eval(ev.putrue);  puWeight[1]=puUpWgtGr->Eval(ev.putrue); puWeight[2]=puDownWgtGr->Eval(ev.putrue);
	}

      float norm=normH ? normH->GetBinContent(1) : 1.0;
      float wgt          (norm*lepSF[0]                                *puWeight[0]);
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

      //nominal selection
      if(nJets>=1)
	{
	  int nJetsCat=TMath::Min((int)nJets,(int)4);
	  int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  
	  allPlots["njetsnbtags_nom"]->Fill(binToFill,wgt);
	  allPlots["njetsnbtags_qcdScaleDown"]->Fill(binToFill,wgtQCDScaleLo);
	  allPlots["njetsnbtags_qcdScaleUp"]->Fill(binToFill,wgtQCDScaleHi);
	  allPlots["njetsnbtags_puUp"]->Fill(binToFill,wgtPuUp);
	  allPlots["njetsnbtags_puDown"]->Fill(binToFill,wgtPuDown);
	  allPlots["njetsnbtags_muEffUp"]->Fill(binToFill,wgtMuEffUp);
	  allPlots["njetsnbtags_muEffDown"]->Fill(binToFill,wgtMuEffDown);
	  allPlots["njetsnbtags_umetUp"]->Fill(binToFill,wgt);
	  allPlots["njetsnbtags_umetDown"]->Fill(binToFill,wgt);
	  
	  std::vector<TString> catsToFill(2,Form("%dj",nJetsCat));
	  catsToFill[1]+= Form("%dt",nBtagsCat);
	  
	  std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);
	  for(size_t icat=0; icat<2; icat++)
	    {
	      TString tag=catsToFill[icat];
	      if(rIt!=lumiMap.end()) {
		Int_t runCtr=std::distance(lumiMap.begin(),rIt);
		allPlots["ratevsrun_"+tag]->Fill(runCtr,1.e+6/rIt->second);
	      }
	      allPlots["lpt_"+tag]->Fill(ev.l_pt,wgt);
	      allPlots["lsip3d_"+tag]->Fill(ev.l_ip3dsig,wgt);
	      allPlots["lchiso_"+tag]->Fill(ev.l_chargedHadronIso,wgt);
	      allPlots["lchreliso_"+tag]->Fill(ev.l_chargedHadronIso/ev.l_pt,wgt);
	      allPlots["leta_"+tag]->Fill(ev.l_eta,wgt);
	      allPlots["jpt_"+tag]->Fill(ev.j_pt[ selJetsIdx[0] ],wgt);
	      allPlots["jeta_"+tag]->Fill(fabs(ev.j_eta[ selJetsIdx[0] ]),wgt);
	      allPlots["csv_"+tag]->Fill(ev.j_csv[ selJetsIdx[0] ],wgt);
	      allPlots["ht_"+tag]->Fill(htsum,wgt);
	      allPlots["nvtx_"+tag]->Fill(ev.nvtx,wgt);
	      allPlots["met_"+tag]->Fill(ev.met_pt,wgt);
	      allPlots["metphi_"+tag]->Fill(ev.met_phi,wgt);
	      allPlots["mt_nom_"+tag]->Fill(mt,wgt);
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

	      if(nBtagsCat>0)
		{
		  allPlots["minmlb_nom_"+tag]->Fill(minMlb,wgt);
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
		}	  
	    }
	}

      if(nJetsJESHi>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJESHi,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
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
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
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
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
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
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_jerDown"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots["mt_jerDown_"+tag]->Fill(mtJERdown,wgt);
	  allPlots["minmlb_jerDown_"+tag]->Fill(minMlbJERdown,wgt);
	}
      if(nJets>=1)
	{
	  int nJetsCat=TMath::Min((int)nJets,(int)4);
          
	  int nBtagsCat=TMath::Min((int)nBtagsBeffLo,(int)2);
          int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_beffDown"]->Fill(binToFill,wgt); 
	  
	  nBtagsCat=TMath::Min((int)nBtagsBeffHi,(int)2);
	  binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots["njetsnbtags_beffUp"]->Fill(binToFill,wgt); 

	  nBtagsCat=TMath::Min((int)nBtagsMistagLo,(int)2);
	  binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots["njetsnbtags_mistagDown"]->Fill(binToFill,wgt); 

	  nBtagsCat=TMath::Min((int)nBtagsMistagHi,(int)2);
	  binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots["njetsnbtags_mistagUp"]->Fill(binToFill,wgt); 
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
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}

//
std::map<Int_t,Float_t> lumiPerRun()
{
  std::map<Int_t,Float_t> toReturn;
toReturn[ 256630 ]=  948417.609 ;
toReturn[ 256673 ]=   5534.408  ;
toReturn[ 256674 ]=  92567.485  ;
toReturn[ 256675 ]= 7099162.193 ;
toReturn[ 256676 ]= 9172872.886 ;
toReturn[ 256677 ]= 15581756.928;
toReturn[ 256801 ]= 8830347.675 ;
toReturn[ 256842 ]=  16510.582  ;
toReturn[ 256843 ]= 37131085.338;
toReturn[ 256866 ]=  58250.406  ;
toReturn[ 256867 ]= 4546508.033 ;
toReturn[ 256868 ]= 22542014.201;
toReturn[ 256869 ]= 1539580.832 ;
toReturn[ 256926 ]= 1499855.808 ;
toReturn[ 256941 ]= 8741029.381 ;
toReturn[ 257461 ]= 3057928.782 ;
toReturn[ 257531 ]= 8418598.194 ;
toReturn[ 257599 ]= 4940876.751 ;
toReturn[ 257613 ]= 75519819.209;
toReturn[ 257614 ]=  850778.922 ;
toReturn[ 257645 ]= 62388624.946;
toReturn[ 257682 ]= 13053256.987;
toReturn[ 257722 ]=  719350.314 ;
toReturn[ 257723 ]= 5941442.106 ;
toReturn[ 257735 ]=  521278.124 ;
toReturn[ 257751 ]= 27029514.967;
toReturn[ 257804 ]=  210956.374 ;
toReturn[ 257805 ]= 17038078.687;
toReturn[ 257816 ]= 24328019.178;
toReturn[ 257819 ]= 15147148.510;
toReturn[ 257968 ]= 16769109.914;
toReturn[ 257969 ]= 39179793.996;
toReturn[ 258129 ]= 5813530.480 ;
toReturn[ 258136 ]= 3617731.160 ;
toReturn[ 258157 ]= 3866329.715 ;
toReturn[ 258158 ]=105571609.093;
toReturn[ 258159 ]= 25531210.007;
toReturn[ 258177 ]=101938042.657;
toReturn[ 258211 ]= 6371404.543 ;
toReturn[ 258213 ]= 11238447.671;
toReturn[ 258214 ]= 15172855.551;
toReturn[ 258215 ]=  411364.505 ;
toReturn[ 258287 ]= 12966493.641;
toReturn[ 258403 ]= 15612888.766;
toReturn[ 258425 ]= 10170405.992;
toReturn[ 258426 ]=  751812.067 ;
toReturn[ 258427 ]= 7901302.746 ;
toReturn[ 258428 ]= 11578208.046;
toReturn[ 258432 ]=  279944.749 ;
toReturn[ 258434 ]= 30536738.787;
toReturn[ 258440 ]= 44624917.855;
toReturn[ 258444 ]= 2079367.112 ;
toReturn[ 258445 ]= 16403375.902;
toReturn[ 258446 ]= 7474593.995 ;
toReturn[ 258448 ]= 35705866.510;
toReturn[ 258655 ]=  383658.834 ;
toReturn[ 258656 ]= 25581933.798;
toReturn[ 258694 ]= 15219679.888;
toReturn[ 258702 ]= 29093248.654;
toReturn[ 258703 ]= 31138680.065;
toReturn[ 258705 ]= 7604760.367 ;
toReturn[ 258706 ]= 52122692.407;
toReturn[ 258712 ]= 34495799.123;
toReturn[ 258713 ]= 10164347.291;
toReturn[ 258714 ]= 4168945.356 ;
toReturn[ 258741 ]= 4446752.908 ;
toReturn[ 258742 ]= 59299810.293;
toReturn[ 258745 ]= 20966777.757;
toReturn[ 258749 ]= 44752319.500;
toReturn[ 258750 ]= 14330984.460;
  return toReturn;
};

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
