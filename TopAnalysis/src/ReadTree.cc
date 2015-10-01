#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include "UserCode/TopAnalysis/interface/MiniEvent.h"
 
#include <vector>
#include <iostream>
#include <algorithm>

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "FWCore/Framework/interface/ESHandle.h"

using namespace edm;

void ReadTree(TString filename,TString output,int chToSelect,bool isTTbar){

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

  std::map<TString, TH1F *> systVars;
  systVars["qcdScaleLo"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_qcdScaleLo");
  systVars["qcdScaleHi"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_qcdScaleHi");
  systVars["hdampScaleLo"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_hdampScaleLo");
  systVars["hdampScaleHi"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_hdampScaleHi");
  systVars["jesLo"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_jesLo");
  systVars["jesHi"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_jesHi");
  systVars["jerLo"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_jerLo");
  systVars["jerHi"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_jerHi");
  systVars["beffLo"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_beffLo");
  systVars["beffHi"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_beffHi");
  systVars["mistagLo"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_mistagLo");
  systVars["mistagHi"] = (TH1 *)bjetcutflow->Clone("bjetcutflow_mistagHi");

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
 JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("Jec11_V12_Uncertainty_AK5PF.txt");

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
  if(abs(lepton_id) == 13 && !ev.muTrigger) continue;
  if(abs(lepton_id) == 11 && !ev.elTrigger) continue;

  //select jets
  uint32_t nJets(0), nJetsJESLo(0),nJetsJESHi(0), nJetsJERLo(0), nJetsJERHi(0);
  uint32_t nBtags(0), nBtagsBeffLo(0), nBtagsBeffHi(0), nBtagsMistagLo(0),nBtagsMistagHi(0);
  for (int k=0; k<ev.nj;k++){
    //check pt and eta of this jet
    float jet_pt  = ev.j_pt[k], jet_eta = ev.j_eta[k], csv = ev.j_csv[k];
    if(fabs(jet_eta) < 2.5) continue;
    
    if(jet_pt > 30)
      {
	nJets ++;
	bool isBTagged(csv>0.890);
	
	//update csv decision based on the flavour of the jet
	bool isHeavyFlavor(abs(ev.j_flav[k])==5 || abs(ev.j_flav[k])==4);
	
	//TODO USING B-TAGGING CODE!!!!!
	
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
    float JERCor(1.0),JERCor_UP(1.0),JERCor_DOWN(1.0);
    
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

  float wgt(1.0),wgtQCDScaleLo(1.0),wgtQCDScaleHi(1.0),wgthdampScaleLo(1.0),wgthdampScaleHi(1.0);
  if(isTTbar) 
  {
  	wgt=ev.ttbar_w[0];
  	wgtQCDScaleLo=wgt*ttbar_w[9]/ttbar_w[0];
  	wgtQCDScaleHi=wgt*ttbar_w[5]/ttbar_w[0];
  	wgthdampScaleLo=wgt*ttbar_w[ttbar_nw-17]/ttbar_w[0];
  	wgthdampScaleHi=wgt*wgt*ttbar_w[ttbar_nw-9]/ttbar_w[0];
  }
  
  //fill cutflow histos
  if(nJets >= 2 )               cutflow->Fill(1,wgt);
  if(nJets >= 3 )               cutflow->Fill(2,wgt);
  if(nJets >= 4 )               cutflow->Fill(3,wgt);
  if(nJets >= 4 && nBtags >= 1) cutflow->Fill(4,wgt);
  if(nJets >= 4 && nBtags >= 2) cutflow->Fill(5,wgt);
  
  //main histogram for xsec extraction
  int binToFill(nBtags>=2?2:nBtags);
  binToFill+=3*nJets;
  if(nJets>=2){
	bjetcutflow->Fill(0,wgt);
  	systVars["qcdScaleLo"]->Fill(binToFill,wgtQCDScaleLo);
	systVars["qcdScaleHi"]->Fill(binToFill,wgtQCDScaleHi);
  	systVars["hdampScaleLo"]->Fill(binToFill,wgthdampScaleLo);
  	systVars["hdampScaleHi"]->Fill(binToFill,wgthdampScaleHi);
  }

  binToFill=(nBtags>=2?2:nBtags);
  binToFill+=3*nJetsJESHi;
  if(nJetsJESHi>=2) systVars["jesHi"]->Fill(binToFill,wgt);

  binToFill=(nBtags>=2?2:nBtags);
  binToFill+=3*nJetsJESLo;
  if(nJetsJESLo>=2) systVars["jesLo"]->Fill(binToFill,wgt);
  
  binToFill=(nBtags>=2?2:nBtags);
  binToFill+=3*nJetsJERHi;
  if(nJetsJERHi>=2) systVars["jerHi"]->Fill(binToFill,wgt);

  binToFill=(nBtags>=2?2:nBtags);
  binToFill+=3*nJetsJERLo;
  if(nJetsJERLo>=2) systVars["jerLo"]->Fill(binToFill,wgt);
  
  binToFill=(nBtagsBeffLo>=2?2:nBtagsBeffLo);
  binToFill+=3*nJets;
  if(nJets>=2) systVars["beffLo"]->Fill(binToFill,wgt); 

  binToFill=(nBtagsBeffHi>=2?2:nBtagsBeffHi);
  binToFill+=3*nJets;
  if(nJets>=2) systVars["beffHi"]->Fill(binToFill,wgt); 

  binToFill=(nBtagsMistagLo>=2?2:nBtagsMistagLo);
  binToFill+=3*nJets;
  if(nJets>=2) systVars["mistagLo"]->Fill(binToFill,wgt); 

  binToFill=(nBtagsMistagHi>=2?2:nBtagsMistagHi);
  binToFill+=3*nJets;
  if(nJets>=2) systVars["mistagHi"]->Fill(binToFill,wgt); 

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
  TFile *fOut=TFile::Open(output+"/"+filename,"RECREATE");
  cutflow->Write();
  bjetcutflow->Write();
  for(std::map<TString, TH1F *>::iterator it=systVars.begin(); it!=systVars.end(); it++) it->second->Write();

  runlumi_2j->Write();
  runlumi_3j->Write();
  runlumi_4j->Write();

  // 2-Jets
  leppt_2j->Write();
  lepeta_2j->Write();
  lepphi_2j->Write();
  leptmass_2j->Write();
  
  jetpt_2j->Write();
  jeteta_2j->Write();
  jetcsv_2j->Write();
  numvertices_2j->Write();

  metpt_2j->Write();
  metphi_2j->Write();
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


void RunOverSamples(TString output, int chToSelect){
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
//  "singleEle_2015D_SEP25.root",
//  "singlemu_2015C_SEP25.root"
//  "singlemu_2015D_SEP25.root"
  };
  
  for(size_t i=0; i<sizeof(files)/sizeof(TString); i++){
    bool isTTbar(false);
    if(files[i].Contains("TT_")) isTTbar=true;
    ReadTree(files[i],output,chToSelect,isTTbar);
    }
}
