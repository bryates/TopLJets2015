#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include "UserCode/TopAnalysis/interface/MiniEvent.h"
 
#include <vector>
#include <iostream>
#include <algorithm>


void ReadTree(TString filename,TString output,int chToSelect){

  gROOT->Reset();

  TH1F *cutflow = new TH1F("cutflow",";Cut;Events" ,6,0.,6.);
  cutflow->GetXaxis()->SetBinLabel(1,"preselected");
  cutflow->GetXaxis()->SetBinLabel(2,"#geq 2j");
  cutflow->GetXaxis()->SetBinLabel(3,"#geq 3j");
  cutflow->GetXaxis()->SetBinLabel(4,"#geq 4j");
  cutflow->GetXaxis()->SetBinLabel(5,"#geq 1b-tag");
  cutflow->GetXaxis()->SetBinLabel(6,"#geq 2b-tags");

  TH1F *bjetcutflow = new TH1F("bjetcutflow","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

// Lepton pt 
  TH1F *leppt_2j_leading = new TH1F("leppt_2j_leading",";ch;Events" ,20,0.,300.);
  TH1F *leppt_3j_leading = (TH1F *)leppt_2j_leading->Clone("leppt_3j_leading");
  TH1F *leppt_4j_leading = (TH1F *)leppt_2j_leading->Clone("leppt_4j_leading");
  leppt_3j_leading   ->SetName("leppt_3j_leading");
  leppt_4j_leading   ->SetName("leppt_4j_leading");

// Lepton eta
  TH1F *lepeta_2j_leading = new TH1F("lepeta_2j_leading",";ch;Events" ,12,0.,3.);
  TH1F *lepeta_3j_leading = (TH1F *)lepeta_2j_leading->Clone("lepeta_3j_leading");
  TH1F *lepeta_4j_leading = (TH1F *)lepeta_2j_leading->Clone("lepeta_4j_leading");
  lepeta_3j_leading   ->SetName("lepeta_3j_leading");
  lepeta_4j_leading   ->SetName("lepeta_4j_leading");

// Lepton phi
  TH1F *lepphi_2j_leading = new TH1F("lepphi_2j_leading",";ch;Events" ,50,-3.2,3.2);
  TH1F *lepphi_3j_leading = (TH1F *)lepphi_2j_leading->Clone("lepphi_3j_leading");
  TH1F *lepphi_4j_leading = (TH1F *)lepphi_2j_leading->Clone("lepphi_4j_leading");
  lepphi_3j_leading   ->SetName("lepphi_3j_leading");
  lepphi_4j_leading   ->SetName("lepphi_4j_leading");

// Lepton transverse mass
  TH1F *leptmass_2j_leading = new TH1F("leptmass_2j_leading",";Transverse Mass;Events" ,100,0.,200.);
  TH1F *leptmass_3j_leading = (TH1F *)leptmass_2j_leading->Clone("leptmass_3j_leading");
  TH1F *leptmass_4j_leading = (TH1F *)leptmass_2j_leading->Clone("leptmass_4j_leading");
  leptmass_3j_leading   ->SetName("leptmass_3j_leading");
  leptmass_4j_leading   ->SetName("leptmass_4j_leading");

// Jet pt
  TH1F *jetpt_2j_leading = new TH1F("jetpt_2j_leading",";pt;Events" ,20,0.,300.);
  TH1F *jetpt_3j_leading = (TH1F *)jetpt_2j_leading->Clone("jetpt_3j_leading");
  TH1F *jetpt_4j_leading = (TH1F *)jetpt_2j_leading->Clone("jetpt_4j_leading");
  jetpt_3j_leading   ->SetName("jetpt_3j_leading");
  jetpt_4j_leading   ->SetName("jetpt_4j_leading");

// Jet eta  
  TH1F *jeteta_2j_leading = new TH1F("jeteta_2j_leading",";eta;Events" ,12,0.,3.);
  TH1F *jeteta_3j_leading = (TH1F *)jeteta_2j_leading->Clone("jeteta_3j_leading");
  TH1F *jeteta_4j_leading = (TH1F *)jeteta_2j_leading->Clone("jeteta_4j_leading");
  jeteta_3j_leading   ->SetName("jeteta_3j_leading");
  jeteta_4j_leading   ->SetName("jeteta_4j_leading");

// CSV
  TH1F *jetcsv_2j_leading = new TH1F("jetcsv_2j_leading",";csv;Events" ,100,-1.2,1.2);
  TH1F *jetcsv_3j_leading     = (TH1F *)jetcsv_2j_leading->Clone("jetcsv_3j_leading");
  TH1F *jetcsv_4j_leading     = (TH1F *)jetcsv_2j_leading->Clone("jetcsv_4j_leading");
  jetcsv_3j_leading     -> SetName("jetcsv_3j_leading");
  jetcsv_4j_leading     -> SetName("jetcsv_4j_leading");

// numvertices
  TH1F *numvertices_2j_leading = new TH1F("numvertices_2j_leading",";vertices;Events" ,25,0.,50.);
  TH1F *numvertices_3j_leading     = (TH1F *)numvertices_2j_leading->Clone("numvertices_3j_leading");
  TH1F *numvertices_4j_leading     = (TH1F *)numvertices_2j_leading->Clone("numvertices_4j_leading");
  numvertices_3j_leading     -> SetName("numvertices_3j_leading");
  numvertices_4j_leading     -> SetName("numvertices_4j_leading");

// MET pt
  TH1F *metpt_2j_leading = new TH1F("metpt_2j_leading",";ch;Events" ,20,0.,300.);
  TH1F *metpt_3j_leading     = (TH1F *)metpt_2j_leading->Clone("metpt_3j_leading");
  TH1F *metpt_4j_leading     = (TH1F *)metpt_2j_leading->Clone("metpt_4j_leading");
  metpt_3j_leading     -> SetName("metpt_3j_leading");
  metpt_4j_leading     -> SetName("metpt_4j_leading");

// MET phi
  TH1F *metphi_2j_leading = new TH1F("metphi_2j_leading",";ch;Events" ,50,-3.2,3.2);
  TH1F *metphi_3j_leading     = (TH1F *)metphi_2j_leading->Clone("metphi_3j_leading");
  TH1F *metphi_4j_leading     = (TH1F *)metphi_2j_leading->Clone("metphi_4j_leading");
  metphi_3j_leading     -> SetName("metphi_3j_leading");
  metphi_4j_leading     -> SetName("metphi_4j_leading");

// MET transverse mass
  TH1F *mettmass_2j_leading = new TH1F("mettmass_2j_leading",";Transverse Mass;Events" ,100,0.,200.);
  TH1F *mettmass_3j_leading     = (TH1F *)mettmass_2j_leading->Clone("mettmass_3j_leading");
  TH1F *mettmass_4j_leading     = (TH1F *)mettmass_2j_leading->Clone("mettmass_4j_leading");
  mettmass_3j_leading     -> SetName("mettmass_3j_leading");
  mettmass_4j_leading     -> SetName("mettmass_4j_leading");
  

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

  //std::cout <<"Electron Triiger :" <<ev.elTrigger<<" : Muon Trigger "<<ev.muTrigger <<std::endl;
  //select jets
  uint32_t nJets(0), nBtags(0); //,nJets30(0) ;
  int lepton_id(ev.l_id); //Nvertices(ev.nvtx), lepton_id(ev.l_id); 	

  for (int k=0; k<ev.nj;k++){
    //check pt and eta of this jet
    float jet_pt  = ev.j_pt[k], jet_eta = ev.j_eta[k], csv = ev.j_csv[k];
    if(jet_pt > 30 && fabs(jet_eta) < 2.5) 
      {
	nJets ++;
	if(csv>0.890) nBtags ++;
      }
  }
  
  //select according to the lepton id
  if(chToSelect>0 && lepton_id!=chToSelect) continue;
  if(nJets<2) continue;

  float wgt(1.0);

  //fill cutflow histos
  if(nJets >= 2 )               cutflow->Fill(1,wgt);
  if(nJets >= 3 )               cutflow->Fill(2,wgt);
  if(nJets >= 4 )               cutflow->Fill(3,wgt);
  
  if(nJets >= 4 && nBtags >= 1) cutflow->Fill(4,wgt);
  if(nJets >= 4 && nBtags >= 2) cutflow->Fill(5,wgt);
  
  if(nJets == 2 && nBtags == 0 ) bjetcutflow->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow->Fill(1,wgt);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow->Fill(2,wgt);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow->Fill(3,wgt);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow->Fill(4,wgt);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow->Fill(5,wgt);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow->Fill(6,wgt);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow->Fill(7,wgt);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow->Fill(8,wgt);

  
  TH1F *lepptH=leppt_2j_leading, *lepetaH=lepeta_2j_leading,*lepphiH=lepphi_2j_leading,*mtH=leptmass_2j_leading, *jetptH=jetpt_2j_leading, *jetetaH=jeteta_2j_leading, *jetcsvH=jetcsv_2j_leading, *numvtxH=numvertices_2j_leading, *metptH=metpt_2j_leading, *metphiH=metphi_2j_leading, *mettmassH=mettmass_2j_leading;
  if(nJets==3){
    lepptH=leppt_3j_leading; lepetaH=lepeta_3j_leading; lepphiH=lepphi_3j_leading; mtH=leptmass_3j_leading; jetptH=jetpt_3j_leading; jetetaH=jeteta_3j_leading; jetcsvH=jetcsv_3j_leading; numvtxH=numvertices_3j_leading; metptH=metpt_3j_leading, metphiH=metphi_3j_leading, mettmassH=mettmass_3j_leading;
  }
  if(nJets>=4){
    lepptH=leppt_4j_leading; lepetaH=lepeta_4j_leading; lepphiH=lepphi_4j_leading; mtH=leptmass_4j_leading; jetptH=jetpt_4j_leading; jetetaH=jeteta_4j_leading; jetcsvH=jetcsv_4j_leading; numvtxH=numvertices_4j_leading; metptH=metpt_4j_leading, metphiH=metphi_4j_leading, mettmassH=mettmass_4j_leading;

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

  // 2-Jets
  leppt_2j_leading->Write();
  lepeta_2j_leading->Write();
  lepphi_2j_leading->Write();
  leptmass_2j_leading->Write();
  
  jetpt_2j_leading->Write();
  jeteta_2j_leading->Write();
  jetcsv_2j_leading->Write();
  numvertices_2j_leading->Write();

  metpt_2j_leading->Write();
  metphi_2j_leading->Write();
  mettmass_2j_leading->Write();

// 3-Jets
  leppt_3j_leading->Write();
  lepeta_3j_leading->Write();
  lepphi_3j_leading->Write();
  leptmass_3j_leading->Write();
  
  jetpt_3j_leading->Write();
  jeteta_3j_leading->Write();
  jetcsv_3j_leading->Write();
  numvertices_3j_leading->Write();

  metpt_3j_leading->Write();
  metphi_3j_leading->Write();
  mettmass_3j_leading->Write();

// 4-Jets
  leppt_4j_leading->Write();
  lepeta_4j_leading->Write();
  lepphi_4j_leading->Write();
  leptmass_4j_leading->Write();
  
  jetpt_4j_leading->Write();
  jeteta_4j_leading->Write();
  jetcsv_4j_leading->Write();
  numvertices_4j_leading->Write();

  metpt_4j_leading->Write();
  metphi_4j_leading->Write();
  mettmass_4j_leading->Write();

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
  "singleEle_28Aug2015_SEP25.root",
  "singleEle_PromptReco2015C_SEP25.root",
  "singleEle_2015D_SEP25.root",
  "singlemu_2015C_SEP25.root",
  "singlemu_2015D_SEP25.root"
  };
  
  for(size_t i=0; i<sizeof(files)/sizeof(TString); i++){
    ReadTree(files[i],output,chToSelect);
    }
}
