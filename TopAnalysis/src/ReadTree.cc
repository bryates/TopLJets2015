#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include "UserCode/TopAnalysis/interface/MiniEvent.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include <vector>
#include <iostream>
#include <algorithm>

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"


void ReadTree(TString filename,TString output,int chToSelect,bool isTTbar)
{
  gROOT->Reset();

  TH1F *cutflow = new TH1F("cutflow",";Cut;Events" ,6,0.,6.);
  cutflow->GetXaxis()->SetBinLabel(1,"preselected");
  cutflow->GetXaxis()->SetBinLabel(2,"#geq 2j");
  cutflow->GetXaxis()->SetBinLabel(3,"#geq 3j");
  cutflow->GetXaxis()->SetBinLabel(4,"#geq 4j");
  cutflow->GetXaxis()->SetBinLabel(5,"#geq 1b-tag");
  cutflow->GetXaxis()->SetBinLabel(6,"#geq 2b-tags");

  TH1F *bjetcutflow = new TH1F("bjetcutflow",";Category;Events" ,9,0.,9.);
  bjetcutflow->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

  std::map<TString, TH1 *> systVars;
  systVars["qcdScaleDown"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_qcdScaleDown");
  systVars["qcdScaleUp"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_qcdScaleUp");
  systVars["hdampScaleDown"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_hdampScaleDown");
  systVars["hdampScaleUp"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_hdampScaleUp");
  systVars["jesDown"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_jesDown");
  systVars["jesUp"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_jesUp");
  systVars["jerDown"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_jerDown");
  systVars["jerUp"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_jerUp");
  systVars["beffDown"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_beffDown");
  systVars["beffUp"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_beffUp");
  systVars["mistagDown"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_mistagDown");
  systVars["mistagUp"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_mistagUp");

// Lepton pt 
  TH1F *leppt_2j = new TH1F("leppt_2j",";Transverse momentum [GeV];Events" ,20,0.,300.);
  TH1F *leppt_3j = (TH1F *)leppt_2j->Clone("leppt_3j");
  TH1F *leppt_4j = (TH1F *)leppt_2j->Clone("leppt_4j");

// Lepton eta
  TH1F *lepeta_2j = new TH1F("lepeta_2j",";Pseudo-rapidity;Events" ,12,0.,3.);
  TH1F *lepeta_3j = (TH1F *)lepeta_2j->Clone("lepeta_3j");
  TH1F *lepeta_4j = (TH1F *)lepeta_2j->Clone("lepeta_4j");

// Lepton phi
  TH1F *lepphi_2j = new TH1F("lepphi_2j",";#phi [rad];Events" ,50,-3.2,3.2);
  TH1F *lepphi_3j = (TH1F *)lepphi_2j->Clone("lepphi_3j");
  TH1F *lepphi_4j = (TH1F *)lepphi_2j->Clone("lepphi_4j");

// Lepton transverse mass
  TH1F *leptmass_2j = new TH1F("leptmass_2j",";Transverse Mass [GeV];Events" ,100,0.,200.);
  TH1F *leptmass_3j = (TH1F *)leptmass_2j->Clone("leptmass_3j");
  TH1F *leptmass_4j = (TH1F *)leptmass_2j->Clone("leptmass_4j");

// Jet pt
  TH1F *jetpt_2j = new TH1F("jetpt_2j",";Transverse momentum [GeV];Events" ,20,0.,300.);
  TH1F *jetpt_3j = (TH1F *)jetpt_2j->Clone("jetpt_3j");
  TH1F *jetpt_4j = (TH1F *)jetpt_2j->Clone("jetpt_4j");

// Jet eta  
  TH1F *jeteta_2j = new TH1F("jeteta_2j",";Pseudo-rapidity;Events" ,12,0.,3.);
  TH1F *jeteta_3j = (TH1F *)jeteta_2j->Clone("jeteta_3j");
  TH1F *jeteta_4j = (TH1F *)jeteta_2j->Clone("jeteta_4j");

// CSV
  TH1F *jetcsv_2j = new TH1F("jetcsv_2j",";CSV discriminator;Events" ,100,-1.2,1.2);
  TH1F *jetcsv_3j     = (TH1F *)jetcsv_2j->Clone("jetcsv_3j");
  TH1F *jetcsv_4j     = (TH1F *)jetcsv_2j->Clone("jetcsv_4j");

// numvertices
  TH1F *numvertices_2j = new TH1F("numvertices_2j",";Vertex multiplicity;Events" ,50,0.,50.);
  TH1F *numvertices_3j     = (TH1F *)numvertices_2j->Clone("numvertices_3j");
  TH1F *numvertices_4j     = (TH1F *)numvertices_2j->Clone("numvertices_4j");

// MET pt
  TH1F *metpt_2j = new TH1F("metpt_2j",";Missing transverse energy [GeV];Events" ,20,0.,300.);
  TH1F *metpt_3j     = (TH1F *)metpt_2j->Clone("metpt_3j");
  TH1F *metpt_4j     = (TH1F *)metpt_2j->Clone("metpt_4j");

// MET phi
  TH1F *metphi_2j = new TH1F("metphi_2j",";MET #phi [rad];Events" ,50,-3.2,3.2);
  TH1F *metphi_3j     = (TH1F *)metphi_2j->Clone("metphi_3j");
  TH1F *metphi_4j     = (TH1F *)metphi_2j->Clone("metphi_4j");

// MET transverse mass
  TH1F *mettmass_2j = new TH1F("mettmass_2j",";Transverse Mass [GeV];Events" ,100,0.,200.);
  TH1F *mettmass_3j     = (TH1F *)mettmass_2j->Clone("mettmass_3j");
  TH1F *mettmass_4j     = (TH1F *)mettmass_2j->Clone("mettmass_4j");

  //jet uncertainty parameterization
  TString jecUncUrl("${CMSSW_BASE}/src/UserCode/TopAnalysis/jec/Summer15_50nsV5_DATA_Uncertainty_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(jecUncUrl.Data());

  // setup calibration readers
  TString btagUncUrl("${CMSSW_BASE}/src/UserCode/TopAnalysis/jec/CSVv2.csv");
  BTagCalibration calib("csvv2", btagUncUrl.Data());
  BTagCalibrationReader btvreader(&calib,               // calibration instance
			       BTagEntry::OP_LOOSE,  // operating point
			       "comb",               // measurement type
			       "central");           // systematics type
  BTagCalibrationReader btvreader_up(&calib, BTagEntry::OP_LOOSE, "comb", "up");  // sys up
  BTagCalibrationReader btvreader_do(&calib, BTagEntry::OP_LOOSE, "comb", "down");  // sys down
  
  //read tree from file
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);

  //get the original number of events in the dataset
  TH1F *origCutFlow=(TH1F *)f->Get("demo/cutflow");
  if(origCutFlow) {
    Float_t origEvents=origCutFlow->GetBinContent(1);
    cutflow->SetBinContent(1,origEvents);
    }

  //get the tree
  TTree *t = (TTree*)f->Get("demo/AnaTree");
  attachToMiniEventTree(t,ev);

  //fill histograms, loop over all entries
  Int_t nentries = (Int_t)t->GetEntriesFast();
  for (Int_t i=0;i<nentries;i++){
  t->GetEntry(i);

  //select according to the lepton id
  int lepton_id(ev.l_id);	
  if(chToSelect>0 && lepton_id!=chToSelect) continue;
  //FIXME 
  //  if(abs(lepton_id) == 13 && !ev.muTrigger) continue;
  //  if(abs(lepton_id) == 11 && !ev.elTrigger) continue;

  //select jets
  uint32_t nJets(0), nJetsJESLo(0),nJetsJESHi(0), nJetsJERLo(0), nJetsJERHi(0);
  uint32_t nBtags(0), nBtagsBeffLo(0), nBtagsBeffHi(0), nBtagsMistagLo(0),nBtagsMistagHi(0);
  for (int k=0; k<ev.nj;k++){
    //check pt and eta of this jet
    float jet_pt  = ev.j_pt[k], jet_eta = ev.j_eta[k], csv = ev.j_csv[k];    
    if(fabs(jet_eta) > 2.5) continue;
    
    if(jet_pt > 30)
      {
	nJets ++;
	bool isBTagged(csv>0.890);
	
	//readout the b-tagging scale factors for this jet
	/*
	BTagEntry::JetFlavor btagFlav( BTagEntry::FLAV_UDSG  );
	if(abs(ev.j_flav[k])==4) btagFlav=BTagEntry::FLAV_C;
	if(abs(ev.j_flav[k])==5) btagFlav=BTagEntry::FLAV_B;	
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

	//
      }
      
    ////////////////////////////////////////////////////////////////////////////////////////

    //jet energy scale variations
    jecUnc->setJetEta(fabs(jet_eta));
    jecUnc->setJetPt(jet_pt);
    double unc = jecUnc->getUncertainty(true);    
    if((jet_pt)*(1+unc)>30) nJetsJESHi++;
    if((jet_pt)*(1-unc)>30) nJetsJESLo++;
    
    //jet energy resolution
    //    float JERCor(1.0),JERCor_UP(1.0),JERCor_DOWN(1.0);
    
    //NEEDS GEN JETS IN THE TREE!!!!!!!!!!
    //float genJet_pt(ev.genj_pt[k]);
    //if(genJet_pt>0)
    //{
    //	JERCor = getJERfactor(jet_pt, jet_eta, genJet_pt);
    //	JERCor_UP = getJERfactor_up(jet_pt, jet_eta, genJet_pt);
    //  JERCor_DOWN = getJERfactor_down(jet_pt, jet_eta, genJet_pt);
    // }
    //    if(JERCor_UP*jet_pt>30) nJetsJERHi++;
    //    if(JERCor_DOWN*jet_pt>30) nJetsJERLo++;
  }
  

  float wgt(1.0),wgtQCDScaleLo(1.0),wgtQCDScaleHi(1.0),wgthdampScaleLo(1.0),wgthdampScaleHi(1.0);
  if(isTTbar) 
    {
      wgt             = ev.ttbar_w[0];
      wgtQCDScaleLo   = wgt*ev.ttbar_w[9]/ev.ttbar_w[0];
      wgtQCDScaleHi   = wgt*ev.ttbar_w[5]/ev.ttbar_w[0];
      wgthdampScaleLo = wgt*ev.ttbar_w[ev.ttbar_nw-17]/ev.ttbar_w[0];
      wgthdampScaleHi = wgt*ev.ttbar_w[ev.ttbar_nw-9]/ev.ttbar_w[0];
    }
  
  //fill cutflow histos
  if(nJets >= 2 )               cutflow->Fill(1,wgt);
  if(nJets >= 3 )               cutflow->Fill(2,wgt);
  if(nJets >= 4 )               cutflow->Fill(3,wgt);
  if(nJets >= 4 && nBtags >= 1) cutflow->Fill(4,wgt);
  if(nJets >= 4 && nBtags >= 2) cutflow->Fill(5,wgt);
  
  //main histogram for xsec extraction
  int binToFill(nBtags>=2?2:nBtags);
  binToFill+=3*(nJets-2);
  if(nJets>=2)
    {
      bjetcutflow->Fill(binToFill,wgt);
      systVars["qcdScaleDown"]->Fill(binToFill,wgtQCDScaleLo);
      systVars["qcdScaleUp"]->Fill(binToFill,wgtQCDScaleHi);
      systVars["hdampScaleDown"]->Fill(binToFill,wgthdampScaleLo);
      systVars["hdampScaleUp"]->Fill(binToFill,wgthdampScaleHi);
    }

  binToFill=(nBtags>=2?2:nBtags);
  binToFill+=3*(nJetsJESHi-2);
  if(nJetsJESHi>=2) systVars["jesUp"]->Fill(binToFill,wgt);

  binToFill=(nBtags>=2?2:nBtags);
  binToFill+=3*(nJetsJESLo-2);
  if(nJetsJESLo>=2) systVars["jesDown"]->Fill(binToFill,wgt);

  binToFill=(nBtags>=2?2:nBtags);
  binToFill+=3*(nJetsJERHi-2);
  if(nJetsJERHi>=2) systVars["jerUp"]->Fill(binToFill,wgt);

  binToFill=(nBtags>=2?2:nBtags);
  binToFill+=3*(nJetsJERLo-2);
  if(nJetsJERLo>=2) systVars["jerDown"]->Fill(binToFill,wgt);
  
  binToFill=(nBtagsBeffLo>=2?2:nBtagsBeffLo);
  binToFill+=3*(nJets-2);
  if(nJets>=2) systVars["beffDown"]->Fill(binToFill,wgt); 

  binToFill=(nBtagsBeffHi>=2?2:nBtagsBeffHi);
  binToFill+=3*(nJets-2);
  if(nJets>=2) systVars["beffUp"]->Fill(binToFill,wgt); 

  binToFill=(nBtagsMistagLo>=2?2:nBtagsMistagLo);
  binToFill+=3*(nJets-2);
  if(nJets>=2) systVars["mistagDown"]->Fill(binToFill,wgt); 

  binToFill=(nBtagsMistagHi>=2?2:nBtagsMistagHi);
  binToFill+=3*(nJets-2);
  if(nJets>=2) systVars["mistagUp"]->Fill(binToFill,wgt); 

  //control histograms for the nominal selection only
  if(nJets<2) continue;

  TH1F *lepptH=leppt_2j, *lepetaH=lepeta_2j,*lepphiH=lepphi_2j,*mtH=leptmass_2j, *jetptH=jetpt_2j, *jetetaH=jeteta_2j, *jetcsvH=jetcsv_2j, *numvtxH=numvertices_2j, *metptH=metpt_2j, *metphiH=metphi_2j, *mettmassH=mettmass_2j;
  
  if(nJets==3){
    lepptH=leppt_3j; lepetaH=lepeta_3j; lepphiH=lepphi_3j; mtH=leptmass_3j; jetptH=jetpt_3j; jetetaH=jeteta_3j; jetcsvH=jetcsv_3j; numvtxH=numvertices_3j; metptH=metpt_3j; metphiH=metphi_3j; mettmassH=mettmass_3j;
  }
  if(nJets>=4){
    lepptH=leppt_4j; lepetaH=lepeta_4j; lepphiH=lepphi_4j; mtH=leptmass_4j; jetptH=jetpt_4j; jetetaH=jeteta_4j; jetcsvH=jetcsv_4j; numvtxH=numvertices_4j; metptH=metpt_4j; metphiH=metphi_4j; mettmassH=mettmass_4j;
  }

  lepptH->Fill(ev.l_pt,wgt);
  lepetaH->Fill(ev.l_eta,wgt);
  lepphiH->Fill(ev.l_phi,wgt);
  mtH->Fill(ev.l_tmass,wgt);

  jetptH->Fill(ev.j_pt[0],wgt);
  jetetaH->Fill(fabs(ev.j_eta[0]),wgt);
  jetcsvH->Fill(fabs(ev.j_csv[0]),wgt);
  numvtxH->Fill(ev.nvtx,wgt);

  metptH->Fill(ev.met_pt,wgt);
  metphiH->Fill(ev.met_phi,wgt);
  mettmassH->Fill(ev.mt,wgt);

  }
  
  
  //close file
  
  f->Close();
  
  //open output file
  
  TFile *fOut=TFile::Open(output+"/"+gSystem->BaseName(filename),"RECREATE");

  cutflow->SetDirectory(fOut);
  cutflow->Write();

  bjetcutflow->SetDirectory(fOut);
  bjetcutflow->Write();
  for(std::map<TString, TH1 *>::iterator it=systVars.begin(); it!=systVars.end(); it++) it->second->Write();

  // 2-Jets
  leppt_2j->SetDirectory(fOut);
  leppt_2j->Write();
  lepeta_2j->SetDirectory(fOut);
  lepeta_2j->Write();
  lepphi_2j->SetDirectory(fOut);
  lepphi_2j->Write();
  leptmass_2j->SetDirectory(fOut);
  leptmass_2j->Write();
  
  jetpt_2j->SetDirectory(fOut);
  jetpt_2j->Write();
  jeteta_2j->SetDirectory(fOut);
  jeteta_2j->Write();

  jetcsv_2j->SetDirectory(fOut);
  jetcsv_2j->Write();

  numvertices_2j->SetDirectory(fOut);
  numvertices_2j->Write();

  metpt_2j->SetDirectory(fOut);
  metpt_2j->Write();

  metphi_2j->SetDirectory(fOut);
  metphi_2j->Write();

  mettmass_2j->SetDirectory(fOut);
  mettmass_2j->Write();

// 3-Jets
  leppt_3j->Write();
  lepeta_3j->Write();
  lepphi_3j->Write();
  leptmass_3j->Write();
  
  jetpt_3j->Write();
  jeteta_3j->Write();
  jetcsv_3j->Write();
  numvertices_3j->Write();

  metpt_3j->Write();
  metphi_3j->Write();
  mettmass_3j->Write();

  // 4-Jets
  leppt_4j->Write();
  lepeta_4j->Write();
  lepphi_4j->Write();
  leptmass_4j->Write();
  
  jetpt_4j->Write();
  jeteta_4j->Write();
  jetcsv_4j->Write();
  numvertices_4j->Write();

  metpt_4j->Write();
  metphi_4j->Write();
  mettmass_4j->Write();

  fOut->Close();
}


void RunOverSamples(TString inDir,TString output, int chToSelect){
  TString files[]={
    "DYJetsToLL_M50_SEP22.root",
    "ST_s_channel_SEP22.root",
    "ST_t_channel_SEP22.root",
    "ST_tW_antitop_SEP22.root",
    "ST_tW_top_SEP22.root",
    "TT_TuneCUETP8M1_SEP22.root",
    "WJetsToLNu_TuneCUETP8M1_SEP22.root",
    "QCD_Pt_80to120_EMEnriched_SEP22.root",
    "QCD_Pt_120to170_EMEnriched_SEP22.root",
    "QCD_Pt_170to300_EMEnriched_SEP22.root",
    "QCD_Pt_300toInf_EMEnriched_SEP22.root",
    "singleEle_PromptReco2015C_SEP25.root"
    //"singleEle_2015D_SEP25.root",
    "singlemu_2015C_SEP25.root"
    //"singlemu_2015D_SEP25.root"
  };
  
  for(size_t i=0; i<sizeof(files)/sizeof(TString); i++){
    bool isTTbar(false);
    if(files[i].Contains("TT_")) isTTbar=true;

    if(files[i].Contains("singlemu") && chToSelect==11) continue;
    if(files[i].Contains("singleEle") && chToSelect==13) continue;

    TString url(files[i]);
    if(inDir!="") url=inDir+"/"+files[i];
    std::cout << "Processing " << url << " " << output << " " << chToSelect << " " << isTTbar << std::endl;
    ReadTree(url,output,chToSelect,isTTbar);
  }
}
