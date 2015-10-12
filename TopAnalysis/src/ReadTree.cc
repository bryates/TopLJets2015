#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/ReadTree.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"

using namespace std;


Int_t getSecVtxBinToFill(Float_t firstVtxMass,Float_t secondVtxMass,Int_t nJets,Int_t nsvtxMassBins,Float_t minSvtxMass,Float_t maxSvtxMass)
{
  Int_t nJetsBin(nJets>4 ? 4 : nJets);
  Int_t nVtx(0); nVtx+=(firstVtxMass>0); nVtx += (secondVtxMass>0); 
  Int_t secvtxBinToFill(0);
  if(nVtx==1) secvtxBinToFill=(Int_t)nsvtxMassBins*(TMath::Max(TMath::Min(firstVtxMass,maxSvtxMass),minSvtxMass)-minSvtxMass)/maxSvtxMass+1;
  if(nVtx==2) secvtxBinToFill=(Int_t)nsvtxMassBins*(TMath::Max(TMath::Min(secondVtxMass,maxSvtxMass),minSvtxMass)-minSvtxMass)/maxSvtxMass+nsvtxMassBins+1;
  secvtxBinToFill += (nJetsBin-1)*(2*nsvtxMassBins+1);
  return secvtxBinToFill;
}

Int_t getBtagCatBinToFill(Int_t nBtags,Int_t nJets)
{
  Int_t nJetsBin(nJets>4 ? 4 : nJets);

  Int_t binToFill(nBtags>=2?2:nBtags);
  binToFill+=3*(nJetsBin-1);
  return binToFill;
}

void ReadTree(TString filename,
	      TString outname,
	      Int_t channelSelection, 
	      Int_t chargeSelection, 
	      Float_t norm, 
	      Bool_t isTTbar,
	      FlavourSplitting flavourSplitting,
	      GenWeightMode genWgtMode,
	      TGraph *puWgtGr,TGraph *puUpWgtGr,TGraph *puDownWgtGr)
{
  //book histograms
  std::map<TString, TH1 *> allPlots;

  std::map<Int_t,Float_t> lumiMap=lumiPerRun();
  for(Int_t ij=1; ij<=4; ij++)
    {
      TString tag(Form("%dj",ij));
      allPlots["ratevsrun_"+tag] = new TH1F("ratevsrun_"+tag,";Run number; Events/pb",lumiMap.size(),0,lumiMap.size());
      Int_t runCtr(0);
      for(std::map<Int_t,Float_t>::iterator it=lumiMap.begin(); it!=lumiMap.end(); it++,runCtr++)
	allPlots["ratevsrun_"+tag]->GetXaxis()->SetBinLabel(runCtr+1,Form("%d",it->first));
      allPlots["lpt_"+tag]  = new TH1F("lpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
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
      allPlots["mt_"+tag] = new TH1F("mt_"+tag,";Transverse Mass [GeV];Events" ,100,0.,200.);
      allPlots["secvtxmass_"+tag] = new TH1F("secvtxmass_"+tag,";SecVtx Mass [GeV];Events" ,10,0.,5.);
      allPlots["secvtx3dsig_"+tag] = new TH1F("secvtx3dsig_"+tag,";SecVtx 3D sig;Events" ,10,0.,100.);
    }

  Int_t nsvtxMassBins=allPlots["secvtxmass_4j"]->GetXaxis()->GetNbins(); 
  Float_t maxSvtxMass=allPlots["secvtxmass_4j"]->GetXaxis()->GetXmax();
  Float_t minSvtxMass=allPlots["secvtxmass_4j"]->GetXaxis()->GetXmin();
  allPlots["catcountSecVtx"] = new TH1F("catcountSecVtx",";SecVtx Mass [GeV];Events",4*(2*nsvtxMassBins+1),0.,4*(2*nsvtxMassBins+1));
  for(Int_t njets=1; njets<=4; njets++)
    {
      Int_t startBin=(njets-1)*(2*nsvtxMassBins+1)+1;
      allPlots["catcountSecVtx"]->GetXaxis()->SetBinLabel(startBin,Form("%dj,0v",njets));
      allPlots["catcountSecVtx"]->GetXaxis()->SetBinLabel(startBin+1,"0");
      allPlots["catcountSecVtx"]->GetXaxis()->SetBinLabel(startBin+6,Form("%dj,1v",njets));
      allPlots["catcountSecVtx"]->GetXaxis()->SetBinLabel(startBin+nsvtxMassBins,"5");
      allPlots["catcountSecVtx"]->GetXaxis()->SetBinLabel(startBin+1+nsvtxMassBins,"0");
      allPlots["catcountSecVtx"]->GetXaxis()->SetBinLabel(startBin+6+nsvtxMassBins,Form("%dj,2v",njets));
      allPlots["catcountSecVtx"]->GetXaxis()->SetBinLabel(startBin+2*nsvtxMassBins,"5");
    }
  allPlots["catcount"] = new TH1F("catcount",";Category;Events" ,12,0.,12.);
  allPlots["catcount"]->GetXaxis()->SetBinLabel(1, "1j,=0b");
  allPlots["catcount"]->GetXaxis()->SetBinLabel(2, "1j,=1b");
  allPlots["catcount"]->GetXaxis()->SetBinLabel(3, "");
  allPlots["catcount"]->GetXaxis()->SetBinLabel(4, "2j,=0b");
  allPlots["catcount"]->GetXaxis()->SetBinLabel(5, "2j,=1b");
  allPlots["catcount"]->GetXaxis()->SetBinLabel(6, "2j,#geq2b");
  allPlots["catcount"]->GetXaxis()->SetBinLabel(7, "3j,=0b");
  allPlots["catcount"]->GetXaxis()->SetBinLabel(8, "3j,=1b");
  allPlots["catcount"]->GetXaxis()->SetBinLabel(9, "3j,#geq2b");
  allPlots["catcount"]->GetXaxis()->SetBinLabel(10,"4j,=0b");
  allPlots["catcount"]->GetXaxis()->SetBinLabel(11,"4j,=1b");
  allPlots["catcount"]->GetXaxis()->SetBinLabel(12,"4j,#geq2b");
  TString systs[]={"puUp","puDown","muEffUp","muEffDown","eEffUp","eEffDown","qcdScaleDown","qcdScaleUp","hdampScaleDown","hdampScaleUp","jesDown","jesUp","jerDown","jerUp","beffDown","beffUp","mistagDown","mistagUp"};
  for(size_t i=0; i<sizeof(systs)/sizeof(TString); i++)
    {
      allPlots["catcount_"+systs[i]]       = (TH1 *)allPlots["catcount"]->Clone("catcount_"+systs[i]);
      allPlots["catcountSecVtx_"+systs[i]] = (TH1 *)allPlots["catcountSecVtx"]->Clone("catcountSecVtx_"+systs[i]);
    }

  for (auto& it : allPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }


  //jet uncertainty parameterization
  TString jecUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/Summer15_50nsV5_DATA_Uncertainty_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(jecUncUrl.Data());
  
  // setup calibration readers
  TString btagUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CSVv2.csv");
  gSystem->ExpandPathName(btagUncUrl);

  /*
  BTagCalibration calib("csvv2", btagUncUrl.Data());
  BTagCalibrationReader btvreader(&calib,               // calibration instance
				  BTagEntry::OP_LOOSE,  // operating point
				  "comb",               // measurement type
				  "central");           // systematics type
  BTagCalibrationReader btvreader_up(&calib, BTagEntry::OP_LOOSE, "comb", "up");  // sys up
  BTagCalibrationReader btvreader_do(&calib, BTagEntry::OP_LOOSE, "comb", "down");  // sys down
  */

  //read tree from file
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev);

  //loop over events
  Int_t nentries(t->GetEntriesFast());
  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  for (Int_t i=0;i<t->GetEntriesFast();i++)
    {
      t->GetEntry(i);
      
      //select according to the lepton id/charge
      if(channelSelection!=0)
	{
	  if(abs(ev.l_id)!=abs(channelSelection)) continue;
	  if(channelSelection==1300 && (ev.l_chargedHadronIso/ev.l_pt<0.06 || ev.l_chargedHadronIso/ev.l_pt>0.2)) continue;
	}
      if(chargeSelection!=0 &&  ev.l_charge!=chargeSelection) continue;

      //apply trigger requirement
      if((abs(ev.l_id) == 13 || abs(ev.l_id) == 1300) && ((ev.muTrigger>>0)&0x1)==0) continue;
      if((abs(ev.l_id) == 11 || abs(ev.l_id) == 1100) && ((ev.elTrigger>>0)&0x1)==0) continue;
      
      //select jets
      Int_t nudsgJets(0),ncJets(0), nbJets(0);
      uint32_t nJets(0), nJetsJESLo(0),nJetsJESHi(0), nJetsJERLo(0), nJetsJERHi(0);
      uint32_t nBtags(0), nBtagsBeffLo(0), nBtagsBeffHi(0), nBtagsMistagLo(0),nBtagsMistagHi(0);
      Float_t htsum(0),firstVtxMass(0.),firstVtxLxySig(-9999.),secondVtxMass(0.),secondVtxLxySig(-9999.);
      std::vector<int> selJetsIdx;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  float jet_pt  = ev.j_pt[k], jet_eta = ev.j_eta[k], csv = ev.j_csv[k];    
	  if(fabs(jet_eta) > 2.5) continue;
	  if(jet_pt > 30)
	    {
	      nJets ++;
	      selJetsIdx.push_back(k);
	      htsum += jet_pt;
	      bool isBTagged(csv>0.890);
	      if(ev.j_vtx3DSig[k]>firstVtxLxySig)
		{
		  secondVtxLxySig=firstVtxLxySig;
		  secondVtxMass=firstVtxMass;
		  firstVtxLxySig=ev.j_vtx3DSig[k];
		  firstVtxMass=ev.j_vtxmass[k];
		}
	      else if(ev.j_vtx3DSig[k]>secondVtxLxySig)
		{
		  secondVtxLxySig=ev.j_vtx3DSig[k];
		  secondVtxMass=ev.j_vtxmass[k];
		}

	      //BTagEntry::JetFlavor btagFlav( BTagEntry::FLAV_UDSG  );
	      if(abs(ev.j_hadflav[k])==4)      { /*btagFlav=BTagEntry::FLAV_C;*/ ncJets++; }
	      else if(abs(ev.j_hadflav[k])==5) { /*btagFlav=BTagEntry::FLAV_B;*/ nbJets++; }
	      else nudsgJets++;


	      //readout the b-tagging scale factors for this jet
	      /*
		BTagEntry::JetFlavor btagFlav( BTagEntry::FLAV_UDSG  );
		if(abs(ev.j_hadflav[k])==4) btagFlav=BTagEntry::FLAV_C;
		if(abs(ev.j_hadflav[k])==5) btagFlav=BTagEntry::FLAV_B;	
		float jetBtagSF(1.0), jetBtagSFUp(1.0), jetBtagSFDown(1.0);
		if (jet_pt < 1000.) {
		jetBtagSF = btvreader.eval(btagFlav, jet_eta, jet_pt);
		jetBtagSFUp = btvreader_up.eval(btagFlav, jet_eta, jet_pt);
		jetBtagSFDown = btvreader_do.eval(btagFlav, jet_eta, jet_pt);
		}
	      */
	
	      nBtags += isBTagged;
	      nBtagsBeffLo += isBTagged;
	      nBtagsBeffHi += isBTagged;
	      nBtagsMistagLo += isBTagged;
	      nBtagsMistagHi += isBTagged;
	    }

	  //jet energy scale variations
	  jecUnc->setJetEta(fabs(jet_eta));
	  jecUnc->setJetPt(jet_pt);
	  double unc = jecUnc->getUncertainty(true);    
	  if((jet_pt)*(1+unc)>30) nJetsJESHi++;
	  if((jet_pt)*(1-unc)>30) nJetsJESLo++;
	  
	  //jet energy resolution
	  std::vector<float> jerScale(3,1.0);
	  float genJet_pt(ev.genj_pt[k]);
	  if(!ev.isData && genJet_pt>0)
	      jerScale=getJetResolutionScales(jet_pt,jet_eta,genJet_pt);
	  if(jerScale[1]*jet_pt>30) nJetsJERLo++;
	  if(jerScale[2]*jet_pt>30) nJetsJERHi++;
	}

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
      float wgt          (norm*lepSF[0]                                *puWeight[0]);
      float wgtPuUp      (norm*lepSF[0]                                *puWeight[1]);
      float wgtPuDown    (norm*lepSF[0]                                *puWeight[2]);
      float wgtMuEffUp   (norm*(abs(ev.l_id)==13 ? lepSF[1] : lepSF[0])*puWeight[0]);
      float wgtMuEffDown (norm*(abs(ev.l_id)==13 ? lepSF[2] : lepSF[0])*puWeight[0]);
      float wgtElEffUp   (norm*(abs(ev.l_id)==11 ? lepSF[1] : lepSF[0])*puWeight[0]);
      float wgtElEffDown (norm*(abs(ev.l_id)==11 ? lepSF[2] : lepSF[0])*puWeight[0]);

      float wgtQCDScaleLo(wgt),wgtQCDScaleHi(wgt),wgthdampScaleLo(wgt),wgthdampScaleHi(wgt);
      if(genWgtMode!=NOGENWGT && !ev.isData) 
	{
	  wgt          *= ev.ttbar_w[0];
	  wgtPuUp      *= ev.ttbar_w[0];
	  wgtPuDown    *= ev.ttbar_w[0];
	  wgtMuEffUp   *= ev.ttbar_w[0];
	  wgtMuEffDown *= ev.ttbar_w[0];
	  wgtElEffUp   *= ev.ttbar_w[0];
	  wgtElEffDown *= ev.ttbar_w[0];
	}
      if(isTTbar)
	{
	  wgtQCDScaleLo   = wgt*ev.ttbar_w[9]/ev.ttbar_w[0];
	  wgtQCDScaleHi   = wgt*ev.ttbar_w[5]/ev.ttbar_w[0];
	  wgthdampScaleLo = wgt*ev.ttbar_w[ev.ttbar_nw-17]/ev.ttbar_w[0];
	  wgthdampScaleHi = wgt*ev.ttbar_w[ev.ttbar_nw-9]/ev.ttbar_w[0];
	}
      
      //main histogram for xsec extraction
      int binToFill=getBtagCatBinToFill(nBtags,nJets);
      int secvtxBinToFill=getSecVtxBinToFill(firstVtxMass,secondVtxMass,nJets,nsvtxMassBins, minSvtxMass,maxSvtxMass);
      if(nJets>=1)
	{
	  allPlots["catcountSecVtx"]->Fill(secvtxBinToFill,wgt);
	  allPlots["catcountSecVtx_qcdScaleDown"]->Fill(secvtxBinToFill,wgtQCDScaleLo);
	  allPlots["catcountSecVtx_qcdScaleUp"]->Fill(secvtxBinToFill,wgtQCDScaleHi);
	  allPlots["catcountSecVtx_hdampScaleDown"]->Fill(secvtxBinToFill,wgthdampScaleLo);
	  allPlots["catcountSecVtx_hdampScaleUp"]->Fill(secvtxBinToFill,wgthdampScaleHi);
	  allPlots["catcountSecVtx_puUp"]->Fill(secvtxBinToFill,wgtPuUp);
	  allPlots["catcountSecVtx_puDown"]->Fill(secvtxBinToFill,wgtPuDown);
	  allPlots["catcountSecVtx_muEffUp"]->Fill(secvtxBinToFill,wgtMuEffUp);
	  allPlots["catcountSecVtx_muEffDown"]->Fill(secvtxBinToFill,wgtMuEffDown);
	  allPlots["catcountSecVtx_eEffUp"]->Fill(secvtxBinToFill,wgtElEffUp);
	  allPlots["catcountSecVtx_eEffDown"]->Fill(secvtxBinToFill,wgtElEffDown);

	  allPlots["catcount"]->Fill(binToFill,wgt);
	  allPlots["catcount_qcdScaleDown"]->Fill(binToFill,wgtQCDScaleLo);
	  allPlots["catcount_qcdScaleUp"]->Fill(binToFill,wgtQCDScaleHi);
	  allPlots["catcount_hdampScaleDown"]->Fill(binToFill,wgthdampScaleLo);
	  allPlots["catcount_hdampScaleUp"]->Fill(binToFill,wgthdampScaleHi);
	  allPlots["catcount_puUp"]->Fill(binToFill,wgtPuUp);
	  allPlots["catcount_puDown"]->Fill(binToFill,wgtPuDown);
	  allPlots["catcount_muEffUp"]->Fill(binToFill,wgtMuEffUp);
	  allPlots["catcount_muEffDown"]->Fill(binToFill,wgtMuEffDown);
	  allPlots["catcount_eEffUp"]->Fill(binToFill,wgtElEffUp);
	  allPlots["catcount_eEffDown"]->Fill(binToFill,wgtElEffDown);
	}

      binToFill=getBtagCatBinToFill(nBtags,nJetsJESHi);
      secvtxBinToFill=getSecVtxBinToFill(firstVtxMass,secondVtxMass,nJetsJESHi,nsvtxMassBins, minSvtxMass,maxSvtxMass);
      if(nJetsJESHi>=1)
	{
	  allPlots["catcount_jesUp"]->Fill(binToFill,wgt);
	  allPlots["catcountSecVtx_jesUp"]->Fill(secvtxBinToFill,wgt);
	}

      binToFill=getBtagCatBinToFill(nBtags,nJetsJESLo);
      secvtxBinToFill=getSecVtxBinToFill(firstVtxMass,secondVtxMass,nJetsJESLo,nsvtxMassBins, minSvtxMass,maxSvtxMass);
      if(nJetsJESLo>=1) 
	{
	  allPlots["catcount_jesDown"]->Fill(binToFill,wgt);
	  allPlots["catcountSecVtx_jesDown"]->Fill(secvtxBinToFill,wgt);
	}

      binToFill=getBtagCatBinToFill(nBtags,nJetsJERHi);
      secvtxBinToFill=getSecVtxBinToFill(firstVtxMass,secondVtxMass,nJetsJERHi,nsvtxMassBins, minSvtxMass,maxSvtxMass);
      if(nJetsJERHi>=1) 
	{
	  allPlots["catcount_jerUp"]->Fill(binToFill,wgt);
	  allPlots["catcountSecVtx_jerUp"]->Fill(secvtxBinToFill,wgt);
	}

      binToFill=getBtagCatBinToFill(nBtags,nJetsJERLo);
      secvtxBinToFill=getSecVtxBinToFill(firstVtxMass,secondVtxMass,nJetsJERLo,nsvtxMassBins, minSvtxMass,maxSvtxMass);
      if(nJetsJERHi>=1) 
	{
	  allPlots["catcount_jerDown"]->Fill(binToFill,wgt);
	  allPlots["catcountSecVtx_jerDown"]->Fill(secvtxBinToFill,wgt);
	}
  
      binToFill=getBtagCatBinToFill(nBtagsBeffLo,nJets);
      if(nJets>=1) allPlots["catcount_beffDown"]->Fill(binToFill,wgt); 

      binToFill=getBtagCatBinToFill(nBtagsBeffHi,nJets);
      if(nJets>=1) allPlots["catcount_beffUp"]->Fill(binToFill,wgt); 
      
      binToFill=getBtagCatBinToFill(nBtagsMistagLo,nJets);
      if(nJets>=1) allPlots["catcount_mistagDown"]->Fill(binToFill,wgt); 

      binToFill=getBtagCatBinToFill(nBtagsMistagHi,nJets);
      if(nJets>=1) allPlots["catcount_mistagUp"]->Fill(binToFill,wgt); 

      //control histograms for the nominal selection only
      if(nJets<1) continue;
      
      TString tag("");
      if(nJets>4) tag="4";
      else tag += nJets;
      tag+="j";
      std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);
      if(rIt!=lumiMap.end()) {
	Int_t runCtr=std::distance(lumiMap.begin(),rIt);
	allPlots["ratevsrun_"+tag]->Fill(runCtr,1.e+6/rIt->second);
      }

      allPlots["lpt_"+tag]->Fill(ev.l_pt,wgt);
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
      allPlots["mt_"+tag]->Fill(ev.mt,wgt);
      if(firstVtxMass>0) 
	{
	  allPlots["secvtxmass_"+tag]->Fill(firstVtxMass,wgt);
	  allPlots["secvtx3dsig_"+tag]->Fill(firstVtxLxySig,wgt);
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
  toReturn[256630]=948417.609   ;
  toReturn[256673]=5534.408     ;
  toReturn[256674]=92567.485    ;
  toReturn[256675]=7282554.291  ;
  toReturn[256676]=9172872.886  ;
  toReturn[256677]=15581756.928 ;
  toReturn[256729]=66555084.031 ;
  toReturn[256734]=7199959.798  ;
  toReturn[256801]=8830347.675  ;
  toReturn[256842]=16510.582    ;
  toReturn[256843]=37367788.087 ;
  toReturn[256866]=58250.406    ;
  toReturn[256867]=4546508.033  ;
  toReturn[256868]=22542014.201 ;
  toReturn[256869]=1539580.832  ;
  toReturn[256926]=1499855.808  ;
  toReturn[256941]=8741029.381  ;
  toReturn[257394]=1630030.328  ;
  toReturn[257395]=701401.393   ;
  toReturn[257461]=3057928.782  ;
  toReturn[257531]=8418598.194  ;
  toReturn[257599]=4940876.751  ;
  toReturn[257613]=75639153.224 ;
  toReturn[257614]=850778.922   ;
  toReturn[257645]=62520503.888 ;
  toReturn[257682]=13053256.987 ;
  toReturn[257722]=810552.138   ;
  toReturn[257723]=5941442.106  ;
  toReturn[257735]=521278.124   ;
  toReturn[257751]=27029514.967 ;
  toReturn[257804]=210956.374   ;
  toReturn[257805]=17038078.687 ;
  toReturn[257816]=24328019.178 ;
  toReturn[257819]=15147148.51  ;
  toReturn[257968]=16769109.914 ;
  toReturn[257969]=39179793.996 ;
  toReturn[258129]=5813530.48   ;
  toReturn[258136]=3617731.16   ;
  toReturn[258157]=3866329.715  ;
  toReturn[258158]=29505332.962 ;
  toReturn[258159]=25687049.652 ;
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
