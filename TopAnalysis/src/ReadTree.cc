#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/ReadTree.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <vector>
#include <iostream>
#include <algorithm>


using namespace std;

void ReadTree(TString filename,
	      TString outname,
	      Int_t channelSelection, 
	      Int_t chargeSelection, 
	      Float_t norm, 
	      Bool_t isTTbar,
	      FlavourSplitting flavourSplitting,
	      GenWeightMode genWgtMode)
{
  //book histograms
  std::map<TString, TH1 *> allPlots;
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
  allPlots["catcount_qcdScaleDown"]   = (TH1 *)allPlots["catcount"]->Clone("catcount_qcdScaleDown");
  allPlots["catcount_qcdScaleUp"]     = (TH1 *)allPlots["catcount"]->Clone("catcount_qcdScaleUp");
  allPlots["catcount_hdampScaleDown"] = (TH1 *)allPlots["catcount"]->Clone("catcount_hdampScaleDown");
  allPlots["catcount_hdampScaleUp"]   = (TH1 *)allPlots["catcount"]->Clone("catcount_hdampScaleUp");
  allPlots["catcount_jesDown"]        = (TH1 *)allPlots["catcount"]->Clone("catcount_jesDown");
  allPlots["catcount_jesUp"]          = (TH1 *)allPlots["catcount"]->Clone("catcount_jesUp");
  allPlots["catcount_jerDown"]        = (TH1 *)allPlots["catcount"]->Clone("catcount_jerDown");
  allPlots["catcount_jerUp"]          = (TH1 *)allPlots["catcount"]->Clone("catcount_jerUp");
  allPlots["catcount_beffDown"]       = (TH1 *)allPlots["catcount"]->Clone("catcount_beffDown");
  allPlots["catcount_beffUp"]         = (TH1 *)allPlots["catcount"]->Clone("catcount_beffUp");
  allPlots["catcount_mistagDown"]     = (TH1 *)allPlots["catcount"]->Clone("catcount_mistagDown");
  allPlots["catcount_mistagUp"]       = (TH1 *)allPlots["catcount"]->Clone("catcount_mistagUp");
  for(Int_t ij=1; ij<=4; ij++)
    {
      TString tag(Form("%dj",ij));
      allPlots["lpt_"+tag]  = new TH1F("lpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
      allPlots["leta_"+tag] = new TH1F("leta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
      allPlots["jpt_"+tag]  = new TH1F("jpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
      allPlots["jeta_"+tag] = new TH1F("jeta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
      allPlots["ht_"+tag]   = new TH1F("ht_"+tag,";H_{T} [GeV];Events",20,0,400);
      allPlots["csv_"+tag]  = new TH1F("csv_"+tag,";CSV discriminator;Events",100,-1.2,1.2);
      allPlots["nvtx_"+tag] = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events" ,50,0.,50.);
      allPlots["met_"+tag]  = new TH1F("metpt_"+tag,";Missing transverse energy [GeV];Events" ,20,0.,300.);
      allPlots["metphi_"+tag] = new TH1F("metphi_" + tag,";MET #phi [rad];Events" ,50,-3.2,3.2);
      allPlots["mt_"+tag] = new TH1F("mt_"+tag,";Transverse Mass [GeV];Events" ,100,0.,200.);
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
      if(channelSelection!=0 && abs(ev.l_id)!=abs(channelSelection)) continue;
      if(chargeSelection!=0 &&  ev.l_charge!=chargeSelection) continue;

      //apply trigger requirement
      if(abs(ev.l_id) == 13 && ((ev.muTrigger>>0)&0x1)==0) continue;
      if(abs(ev.l_id) == 11 && ((ev.elTrigger>>0)&0x1)==0) continue;
      
      //select jets
      Int_t nudsgJets(0),ncJets(0), nbJets(0);
      uint32_t nJets(0), nJetsJESLo(0),nJetsJESHi(0), nJetsJERLo(0), nJetsJERHi(0);
      uint32_t nBtags(0), nBtagsBeffLo(0), nBtagsBeffHi(0), nBtagsMistagLo(0),nBtagsMistagHi(0);
      Float_t htsum(0);
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
	  float JERCor_UP(1.0),JERCor_DOWN(1.0);
	  
	  //NEEDS GEN JETS IN THE TREE!!!!!!!!!!
	  //float genJet_pt(ev.genj_pt[k]);
	  //if(genJet_pt>0)
	  //{
	  //	JERCor = getJERfactor(jet_pt, jet_eta, genJet_pt);
	  //	JERCor_UP = getJERfactor_up(jet_pt, jet_eta, genJet_pt);
	  //  JERCor_DOWN = getJERfactor_down(jet_pt, jet_eta, genJet_pt);
	  // }
	  if(JERCor_UP*jet_pt>30) nJetsJERHi++;
	  if(JERCor_DOWN*jet_pt>30) nJetsJERLo++;
	}

      //check if flavour splitting was required
      if(flavourSplitting!=FlavourSplitting::NOFLAVOURSPLITTING)
      {
      	if(flavourSplitting==FlavourSplitting::BSPLITTING)         { if(nbJets==0)    continue; }
      	else if(flavourSplitting==FlavourSplitting::CSPLITTING)    { if(ncJets==0 || nbJets!=0)    continue; }
      	else if(flavourSplitting==FlavourSplitting::UDSGSPLITTING) { if(nudsgJets==0 || ncJets!=0 || nbJets!=0) continue; }
      }
	
      //generator level weights to apply
      float wgt(norm),wgtQCDScaleLo(norm),wgtQCDScaleHi(norm),wgthdampScaleLo(norm),wgthdampScaleHi(norm);
      if(genWgtMode==FULLWEIGHT) wgt *= ev.ttbar_w[0];
      if(genWgtMode==ONLYSIGN)   wgt *= (ev.ttbar_w[0]>0 ? +1.0 : -1.0)*norm;
      if(isTTbar)
	{
	  wgtQCDScaleLo   = wgt*ev.ttbar_w[9]/ev.ttbar_w[0];
	  wgtQCDScaleHi   = wgt*ev.ttbar_w[5]/ev.ttbar_w[0];
	  wgthdampScaleLo = wgt*ev.ttbar_w[ev.ttbar_nw-17]/ev.ttbar_w[0];
	  wgthdampScaleHi = wgt*ev.ttbar_w[ev.ttbar_nw-9]/ev.ttbar_w[0];
	}
      
      //main histogram for xsec extraction
      int binToFill(nBtags>=2?2:nBtags);
      binToFill+=3*(nJets-1);
      if(nJets>=1)
	{
	  allPlots["catcount"]->Fill(binToFill,wgt);
	  allPlots["catcount_qcdScaleDown"]->Fill(binToFill,wgtQCDScaleLo);
	  allPlots["catcount_qcdScaleUp"]->Fill(binToFill,wgtQCDScaleHi);
	  allPlots["catcount_hdampScaleDown"]->Fill(binToFill,wgthdampScaleLo);
	  allPlots["catcount_hdampScaleUp"]->Fill(binToFill,wgthdampScaleHi);
	}
      binToFill=(nBtags>=2?2:nBtags);
      binToFill+=3*(nJetsJESHi-1);
      if(nJetsJESHi>=1) allPlots["catcount_jesUp"]->Fill(binToFill,wgt);

      binToFill=(nBtags>=2?2:nBtags);
      binToFill+=3*(nJetsJESLo-1);
      if(nJetsJESLo>=1) allPlots["catcount_jesDown"]->Fill(binToFill,wgt);

      binToFill=(nBtags>=2?2:nBtags);
      binToFill+=3*(nJetsJERHi-1);
      if(nJetsJERHi>=1) allPlots["catcount_jerUp"]->Fill(binToFill,wgt);

      binToFill=(nBtags>=2?2:nBtags);
      binToFill+=3*(nJetsJERLo-1);
      if(nJetsJERLo>=1) allPlots["catcount_jerDown"]->Fill(binToFill,wgt);
  
      binToFill=(nBtagsBeffLo>=2?2:nBtagsBeffLo);
      binToFill+=3*(nJets-1);
      if(nJets>=1) allPlots["catcount_beffDown"]->Fill(binToFill,wgt); 

      binToFill=(nBtagsBeffHi>=2?2:nBtagsBeffHi);
      binToFill+=3*(nJets-1);
      if(nJets>=1) allPlots["catcount_beffUp"]->Fill(binToFill,wgt); 
      
      binToFill=(nBtagsMistagLo>=2?2:nBtagsMistagLo);
      binToFill+=3*(nJets-1);
      if(nJets>=1) allPlots["catcount_mistagDown"]->Fill(binToFill,wgt); 

      binToFill=(nBtagsMistagHi>=2?2:nBtagsMistagHi);
      binToFill+=3*(nJets-1);
      if(nJets>=1) allPlots["catcount_mistagUp"]->Fill(binToFill,wgt); 

      //control histograms for the nominal selection only
      if(nJets<1) continue;
      
      TString tag("");
      if(nJets>4) tag="4";
      else tag += nJets;
      tag+="j";

      allPlots["lpt_"+tag]->Fill(ev.l_pt,wgt);
      allPlots["leta_"+tag]->Fill(ev.l_eta,wgt);
      allPlots["jpt_"+tag]->Fill(ev.j_pt[ selJetsIdx[0] ],wgt);
      allPlots["jeta_"+tag]->Fill(fabs(ev.j_eta[ selJetsIdx[0] ]),wgt);
      allPlots["csv_"+tag]->Fill(ev.j_csv[ selJetsIdx[0] ],wgt);
      allPlots["ht_"+tag]->Fill(htsum,wgt);
      allPlots["nvtx_"+tag]->Fill(ev.nvtx,wgt);
      allPlots["met_"+tag]->Fill(ev.met_pt,wgt);
      allPlots["metphi_"+tag]->Fill(ev.met_phi,wgt);
      allPlots["mt_"+tag]->Fill(ev.mt,wgt);
    }

  //close input file
  f->Close();

  //save histos to file  
  TString selPrefix("");  
  if(flavourSplitting!=NOFLAVOURSPLITTING) selPrefix=Form("%d_",flavourSplitting);
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+selPrefix+baseName,"RECREATE");
  for (auto& it : allPlots)  { it.second->SetDirectory(fOut); it.second->Write(); }
  fOut->Close();
}
