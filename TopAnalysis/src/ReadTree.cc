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
 
// JEC Uncertainty
  TH1F *upjec_2j_leading  = new TH1F("upjec_2j_leading",";Transverse momentum [GeV];# jets", 100, 0., 250.);
  TH1F *upjec_3j_leading = (TH1F *)TH1F *upjec_2j_leading -> Clone("upjec_3j_leading");
  TH1F *upjec_4j_leading = (TH1F *)TH1F *upjec_2j_leading -> Clone("upjec_4j_leading");
  TH1F *upjec_3j_leading ->SetName("upjec_3j_leading");
  TH1F *upjec_4j_leading ->SetName("upjec_4j_leading");

  TH1F *downjec_2j_leading  = new TH1F("downjec_2j_leading",";Transverse momentum [GeV];# jets", 100, 0., 250.);
  TH1F *downjec_3j_leading = (TH1F *)TH1F *downjec_2j_leading -> Clone("downjec_3j_leading");
  TH1F *downjec_4j_leading = (TH1F *)TH1F *downjec_2j_leading -> Clone("downjec_4j_leading");
  TH1F *downjec_3j_leading ->SetName("downjec_3j_leading");
  TH1F *downjec_4j_leading ->SetName("downjec_4j_leading");

// Run & Luminosity
  TH2F *lumi_2j_leading = new TH2F("lumi_2j_leading",";Run Number;Number of jets" ,10,0.,10.,10,0.,10.);
  TH2F *lumi_3j_leading = (TH2F *)lumi_2j_leading->Clone("lumi_3j_leading");
  TH2F *lumi_4j_leading = (TH2F *)lumi_2j_leading->Clone("lumi_4j_leading");
  lumi_3j_leading -> SetName("lumi_3j_leading");
  lumi_4j_leading -> SetName("lumi_4j_leading");

  TH1F *runlumi_2j = new TH1F("runlumi_2j",";Run;Events" ,18,0.,19.);
  runlumi_2j->GetXaxis()->SetBinLabel(1,"231");
  runlumi_2j->GetXaxis()->SetBinLabel(2,"232");
  runlumi_2j->GetXaxis()->SetBinLabel(3,"790");
  runlumi_2j->GetXaxis()->SetBinLabel(4,"852");
  runlumi_2j->GetXaxis()->SetBinLabel(5,"879");
  runlumi_2j->GetXaxis()->SetBinLabel(6,"906");
  runlumi_2j->GetXaxis()->SetBinLabel(7,"907");
  runlumi_2j->GetXaxis()->SetBinLabel(8,"914");
  runlumi_2j->GetXaxis()->SetBinLabel(9,"630");
  runlumi_2j->GetXaxis()->SetBinLabel(10,"673");
  runlumi_2j->GetXaxis()->SetBinLabel(11,"674");
  runlumi_2j->GetXaxis()->SetBinLabel(12,"675");
  runlumi_2j->GetXaxis()->SetBinLabel(13,"677");
  runlumi_2j->GetXaxis()->SetBinLabel(14,"729");
  runlumi_2j->GetXaxis()->SetBinLabel(15,"734");
  runlumi_2j->GetXaxis()->SetBinLabel(16,"801");
  runlumi_2j->GetXaxis()->SetBinLabel(17,"842");
  runlumi_2j->GetXaxis()->SetBinLabel(18,"843");
/*  runlumi_2j->GetXaxis()->SetBinLabel(19,"231");
  runlumi_2j->GetXaxis()->SetBinLabel(20,"232");
  runlumi_2j->GetXaxis()->SetBinLabel(21,"790");
  runlumi_2j->GetXaxis()->SetBinLabel(22,"852");
  runlumi_2j->GetXaxis()->SetBinLabel(23,"879");
  runlumi_2j->GetXaxis()->SetBinLabel(24,"906");
  runlumi_2j->GetXaxis()->SetBinLabel(25,"907");
  runlumi_2j->GetXaxis()->SetBinLabel(26,"914");
  runlumi_2j->GetXaxis()->SetBinLabel(27,"630");
  runlumi_2j->GetXaxis()->SetBinLabel(28,"673");
  runlumi_2j->GetXaxis()->SetBinLabel(29,"674");
  runlumi_2j->GetXaxis()->SetBinLabel(30,"675");
  runlumi_2j->GetXaxis()->SetBinLabel(31,"677");
  runlumi_2j->GetXaxis()->SetBinLabel(32,"729");
  runlumi_2j->GetXaxis()->SetBinLabel(33,"734");
  runlumi_2j->GetXaxis()->SetBinLabel(34,"801");
  runlumi_2j->GetXaxis()->SetBinLabel(35,"842");
  runlumi_2j->GetXaxis()->SetBinLabel(36,"843");
*/
  TH1F *runlumi_3j = new TH1F("runlumi_3j",";Run;Events" ,18,0.,19.);
  runlumi_3j->GetXaxis()->SetBinLabel(1,"231");
  runlumi_3j->GetXaxis()->SetBinLabel(2,"232");
  runlumi_3j->GetXaxis()->SetBinLabel(3,"790");
  runlumi_3j->GetXaxis()->SetBinLabel(4,"852");
  runlumi_3j->GetXaxis()->SetBinLabel(5,"879");
  runlumi_3j->GetXaxis()->SetBinLabel(6,"906");
  runlumi_3j->GetXaxis()->SetBinLabel(7,"907");
  runlumi_3j->GetXaxis()->SetBinLabel(8,"914");
  runlumi_3j->GetXaxis()->SetBinLabel(9,"630");
  runlumi_3j->GetXaxis()->SetBinLabel(10,"673");
  runlumi_3j->GetXaxis()->SetBinLabel(11,"674");
  runlumi_3j->GetXaxis()->SetBinLabel(12,"675");
  runlumi_3j->GetXaxis()->SetBinLabel(13,"677");
  runlumi_3j->GetXaxis()->SetBinLabel(14,"729");
  runlumi_3j->GetXaxis()->SetBinLabel(15,"734");
  runlumi_3j->GetXaxis()->SetBinLabel(16,"801");
  runlumi_3j->GetXaxis()->SetBinLabel(17,"842");
  runlumi_3j->GetXaxis()->SetBinLabel(18,"843");
/*  runlumi_3j->GetXaxis()->SetBinLabel(19,"231");
  runlumi_3j->GetXaxis()->SetBinLabel(20,"232");
  runlumi_3j->GetXaxis()->SetBinLabel(21,"790");
  runlumi_3j->GetXaxis()->SetBinLabel(22,"852");
  runlumi_3j->GetXaxis()->SetBinLabel(23,"879");
  runlumi_3j->GetXaxis()->SetBinLabel(24,"906");
  runlumi_3j->GetXaxis()->SetBinLabel(25,"907");
  runlumi_3j->GetXaxis()->SetBinLabel(26,"914");
  runlumi_3j->GetXaxis()->SetBinLabel(27,"630");
  runlumi_3j->GetXaxis()->SetBinLabel(28,"673");
  runlumi_3j->GetXaxis()->SetBinLabel(29,"674");
  runlumi_3j->GetXaxis()->SetBinLabel(30,"675");
  runlumi_3j->GetXaxis()->SetBinLabel(31,"677");
  runlumi_3j->GetXaxis()->SetBinLabel(32,"729");
  runlumi_3j->GetXaxis()->SetBinLabel(33,"734");
  runlumi_3j->GetXaxis()->SetBinLabel(34,"801");
  runlumi_3j->GetXaxis()->SetBinLabel(35,"842");
  runlumi_3j->GetXaxis()->SetBinLabel(36,"843");
*/
 TH1F *runlumi_4j = new TH1F("runlumi_4j",";Run;Events" ,18,0.,19.);
  runlumi_4j->GetXaxis()->SetBinLabel(1,"231");
  runlumi_4j->GetXaxis()->SetBinLabel(2,"232");
  runlumi_4j->GetXaxis()->SetBinLabel(3,"790");
  runlumi_4j->GetXaxis()->SetBinLabel(4,"852");
  runlumi_4j->GetXaxis()->SetBinLabel(5,"879");
  runlumi_4j->GetXaxis()->SetBinLabel(6,"906");
  runlumi_4j->GetXaxis()->SetBinLabel(7,"907");
  runlumi_4j->GetXaxis()->SetBinLabel(8,"914");
  runlumi_4j->GetXaxis()->SetBinLabel(9,"630");
  runlumi_4j->GetXaxis()->SetBinLabel(10,"673");
  runlumi_4j->GetXaxis()->SetBinLabel(11,"674");
  runlumi_4j->GetXaxis()->SetBinLabel(12,"675");
  runlumi_4j->GetXaxis()->SetBinLabel(13,"677");
  runlumi_4j->GetXaxis()->SetBinLabel(14,"729");
  runlumi_4j->GetXaxis()->SetBinLabel(15,"734");
  runlumi_4j->GetXaxis()->SetBinLabel(16,"801");
  runlumi_4j->GetXaxis()->SetBinLabel(17,"842");
  runlumi_4j->GetXaxis()->SetBinLabel(18,"843");
/*  runlumi_4j->GetXaxis()->SetBinLabel(19,"231");
  runlumi_4j->GetXaxis()->SetBinLabel(20,"232");
  runlumi_4j->GetXaxis()->SetBinLabel(21,"790");
  runlumi_4j->GetXaxis()->SetBinLabel(22,"852");
  runlumi_4j->GetXaxis()->SetBinLabel(23,"879");
  runlumi_4j->GetXaxis()->SetBinLabel(24,"906");
  runlumi_4j->GetXaxis()->SetBinLabel(25,"907");
  runlumi_4j->GetXaxis()->SetBinLabel(26,"914");
  runlumi_4j->GetXaxis()->SetBinLabel(27,"630");
  runlumi_4j->GetXaxis()->SetBinLabel(28,"673");
  runlumi_4j->GetXaxis()->SetBinLabel(29,"674");
  runlumi_4j->GetXaxis()->SetBinLabel(30,"675");
  runlumi_4j->GetXaxis()->SetBinLabel(31,"677");
  runlumi_4j->GetXaxis()->SetBinLabel(32,"729");
  runlumi_4j->GetXaxis()->SetBinLabel(33,"734");
  runlumi_4j->GetXaxis()->SetBinLabel(34,"801");
  runlumi_4j->GetXaxis()->SetBinLabel(35,"842");
  runlumi_4j->GetXaxis()->SetBinLabel(36,"843");
*/ 
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
//  if(abs(lepton_id) == 13 && !ev.muTrigger) continue;
//  if(abs(lepton_id) == 11 && !ev.elTrigger) continue;

  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
//  iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);
 
   jecUnc->setJetEta(fabs(ev.j_eta[0]));
   jecUnc->setJetPt(ev.j_pt[0]);
   double unc = jecUnc->getUncertainty(true);
   double up_jecpt = (ev.j_pt[0])*(1+unc);
   double down_jecpt = (ev.j_pt[0])*(1-unc);
 
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

  if(nJets == 2 && ev.run == 254231)  runlumi_2j->Fill(1,wgt/0.029);
  if(nJets == 2 && ev.run == 254232)  runlumi_2j->Fill(2,wgt/0.097);
  if(nJets == 2 && ev.run == 254790)  runlumi_2j->Fill(3,wgt/10.36);
  if(nJets == 2 && ev.run == 254852)  runlumi_2j->Fill(4,wgt/0.82);
  if(nJets == 2 && ev.run == 254879)  runlumi_2j->Fill(5,wgt/1.66);
  if(nJets == 2 && ev.run == 254906)  runlumi_2j->Fill(6,wgt/1.48);
  if(nJets == 2 && ev.run == 254907)  runlumi_2j->Fill(7,wgt/1.02);
  if(nJets == 2 && ev.run == 254914)  runlumi_2j->Fill(8,wgt/0.87);
  if(nJets == 2 && ev.run == 256630)  runlumi_2j->Fill(9,wgt/0.95);
  if(nJets == 2 && ev.run == 256673)  runlumi_2j->Fill(10,wgt/0.005);
  if(nJets == 2 && ev.run == 256674)  runlumi_2j->Fill(11,wgt/0.09);
  if(nJets == 2 && ev.run == 256675)  runlumi_2j->Fill(12,wgt/7.28);
  if(nJets == 2 && ev.run == 256677)  runlumi_2j->Fill(13,wgt/15.58);
  if(nJets == 2 && ev.run == 256729)  runlumi_2j->Fill(14,wgt/66.55);
  if(nJets == 2 && ev.run == 256734)  runlumi_2j->Fill(15,wgt/7.2);
  if(nJets == 2 && ev.run == 256801)  runlumi_2j->Fill(16,wgt/8.83);
  if(nJets == 2 && ev.run == 256842)  runlumi_2j->Fill(17,wgt/0.02);
  if(nJets == 2 && ev.run == 256843)  runlumi_2j->Fill(18,wgt/17.48);
/*  if(nJets == 2 && ev.run == 254231)  runlumi_2j->Fill(19,wgt/0.029);
  if(nJets == 2 && ev.run == 254232)  runlumi_2j->Fill(20,wgt/0.097);
  if(nJets == 2 && ev.run == 254790)  runlumi_2j->Fill(21,wgt/10.36);
  if(nJets == 2 && ev.run == 254852)  runlumi_2j->Fill(22,wgt/0.82);
  if(nJets == 2 && ev.run == 254879)  runlumi_2j->Fill(23,wgt/1.66);
  if(nJets == 2 && ev.run == 254906)  runlumi_2j->Fill(24,wgt/1.48);
  if(nJets == 2 && ev.run == 254907)  runlumi_2j->Fill(25,wgt/1.02);
  if(nJets == 2 && ev.run == 254914)  runlumi_2j->Fill(26,wgt/0.87);
  if(nJets == 2 && ev.run == 256630)  runlumi_2j->Fill(27,wgt/0.95);
  if(nJets == 2 && ev.run == 256673)  runlumi_2j->Fill(28,wgt/0.005);
  if(nJets == 2 && ev.run == 256674)  runlumi_2j->Fill(29,wgt/0.09);
  if(nJets == 2 && ev.run == 256675)  runlumi_2j->Fill(30,wgt/7.28);
  if(nJets == 2 && ev.run == 256677)  runlumi_2j->Fill(31,wgt/15.58);
  if(nJets == 2 && ev.run == 256729)  runlumi_2j->Fill(32,wgt/66.55);
  if(nJets == 2 && ev.run == 256734)  runlumi_2j->Fill(33,wgt/7.2);
  if(nJets == 2 && ev.run == 256801)  runlumi_2j->Fill(34,wgt/8.83);
  if(nJets == 2 && ev.run == 256842)  runlumi_2j->Fill(35,wgt/0.02);
  if(nJets == 2 && ev.run == 256843)  runlumi_2j->Fill(36,wgt/17.48);
*/
  if(nJets == 3 && ev.run == 254231)  runlumi_3j->Fill(1,wgt/0.029);
  if(nJets == 3 && ev.run == 254232)  runlumi_3j->Fill(2,wgt/0.097);
  if(nJets == 3 && ev.run == 254790)  runlumi_3j->Fill(3,wgt/10.36);
  if(nJets == 3 && ev.run == 254852)  runlumi_3j->Fill(4,wgt/0.82);
  if(nJets == 3 && ev.run == 254879)  runlumi_3j->Fill(5,wgt/1.66);
  if(nJets == 3 && ev.run == 254906)  runlumi_3j->Fill(6,wgt/1.48);
  if(nJets == 3 && ev.run == 254907)  runlumi_3j->Fill(7,wgt/1.02);
  if(nJets == 3 && ev.run == 254914)  runlumi_3j->Fill(8,wgt/0.87);
  if(nJets == 3 && ev.run == 256630)  runlumi_3j->Fill(9,wgt/0.95);
  if(nJets == 3 && ev.run == 256673)  runlumi_3j->Fill(10,wgt/0.005);
  if(nJets == 3 && ev.run == 256674)  runlumi_3j->Fill(11,wgt/0.09);
  if(nJets == 3 && ev.run == 256675)  runlumi_3j->Fill(12,wgt/7.28);
  if(nJets == 3 && ev.run == 256677)  runlumi_3j->Fill(13,wgt/15.58);
  if(nJets == 3 && ev.run == 256729)  runlumi_3j->Fill(14,wgt/66.55);
  if(nJets == 3 && ev.run == 256734)  runlumi_3j->Fill(15,wgt/7.2);
  if(nJets == 3 && ev.run == 256801)  runlumi_3j->Fill(16,wgt/8.83);
  if(nJets == 3 && ev.run == 256842)  runlumi_3j->Fill(17,wgt/0.02);
  if(nJets == 3 && ev.run == 256843)  runlumi_3j->Fill(18,wgt/17.48);
/*  if(nJets == 3 && ev.run == 254231)  runlumi_3j->Fill(19,wgt/0.029);
  if(nJets == 3 && ev.run == 254232)  runlumi_3j->Fill(20,wgt/0.097);
  if(nJets == 3 && ev.run == 254790)  runlumi_3j->Fill(21,wgt/10.36);
  if(nJets == 3 && ev.run == 254852)  runlumi_3j->Fill(22,wgt/0.82);
  if(nJets == 3 && ev.run == 254879)  runlumi_3j->Fill(23,wgt/1.66);
  if(nJets == 3 && ev.run == 254906)  runlumi_3j->Fill(24,wgt/1.48);
  if(nJets == 3 && ev.run == 254907)  runlumi_3j->Fill(25,wgt/1.02);
  if(nJets == 3 && ev.run == 254914)  runlumi_3j->Fill(26,wgt/0.87);
  if(nJets == 3 && ev.run == 256630)  runlumi_3j->Fill(27,wgt/0.95);
  if(nJets == 3 && ev.run == 256673)  runlumi_3j->Fill(28,wgt/0.005);
  if(nJets == 3 && ev.run == 256674)  runlumi_3j->Fill(29,wgt/0.09);
  if(nJets == 3 && ev.run == 256675)  runlumi_3j->Fill(30,wgt/7.28);
  if(nJets == 3 && ev.run == 256677)  runlumi_3j->Fill(31,wgt/15.58);
  if(nJets == 3 && ev.run == 256729)  runlumi_3j->Fill(32,wgt/66.55);
  if(nJets == 3 && ev.run == 256734)  runlumi_3j->Fill(33,wgt/7.2);
  if(nJets == 3 && ev.run == 256801)  runlumi_3j->Fill(34,wgt/8.83);
  if(nJets == 3 && ev.run == 256842)  runlumi_3j->Fill(35,wgt/0.02);
  if(nJets == 3 && ev.run == 256843)  runlumi_3j->Fill(36,wgt/17.48);
*/
  if(nJets >= 4 && ev.run == 254231)  runlumi_4j->Fill(1,wgt/0.029);
  if(nJets >= 4 && ev.run == 254232)  runlumi_4j->Fill(2,wgt/0.097);
  if(nJets >= 4 && ev.run == 254790)  runlumi_4j->Fill(3,wgt/10.36);
  if(nJets >= 4 && ev.run == 254852)  runlumi_4j->Fill(4,wgt/0.82);
  if(nJets >= 4 && ev.run == 254879)  runlumi_4j->Fill(5,wgt/1.66);
  if(nJets >= 4 && ev.run == 254906)  runlumi_4j->Fill(6,wgt/1.48);
  if(nJets >= 4 && ev.run == 254907)  runlumi_4j->Fill(7,wgt/1.02);
  if(nJets >= 4 && ev.run == 254914)  runlumi_4j->Fill(8,wgt/0.87);
  if(nJets >= 4 && ev.run == 256630)  runlumi_4j->Fill(9,wgt/0.95);
  if(nJets >= 4 && ev.run == 256673)  runlumi_4j->Fill(10,wgt/0.005);
  if(nJets >= 4 && ev.run == 256674)  runlumi_4j->Fill(11,wgt/0.09);
  if(nJets >= 4 && ev.run == 256675)  runlumi_4j->Fill(12,wgt/7.28);
  if(nJets >= 4 && ev.run == 256677)  runlumi_4j->Fill(13,wgt/15.58);
  if(nJets >= 4 && ev.run == 256729)  runlumi_4j->Fill(14,wgt/66.55);
  if(nJets >= 4 && ev.run == 256734)  runlumi_4j->Fill(15,wgt/7.2);
  if(nJets >= 4 && ev.run == 256801)  runlumi_4j->Fill(16,wgt/8.83);
  if(nJets >= 4 && ev.run == 256842)  runlumi_4j->Fill(17,wgt/0.02);
  if(nJets >= 4 && ev.run == 256843)  runlumi_4j->Fill(18,wgt/17.48);
/*  if(nJets >= 4 && ev.run == 254231)  runlumi_4j->Fill(19,wgt/0.029);
  if(nJets >= 4 && ev.run == 254232)  runlumi_4j->Fill(20,wgt/0.097);
  if(nJets >= 4 && ev.run == 254790)  runlumi_4j->Fill(21,wgt/10.36);
  if(nJets >= 4 && ev.run == 254852)  runlumi_4j->Fill(22,wgt/0.82);
  if(nJets >= 4 && ev.run == 254879)  runlumi_4j->Fill(23,wgt/1.66);
  if(nJets >= 4 && ev.run == 254906)  runlumi_4j->Fill(24,wgt/1.48);
  if(nJets >= 4 && ev.run == 254907)  runlumi_4j->Fill(25,wgt/1.02);
  if(nJets >= 4 && ev.run == 254914)  runlumi_4j->Fill(26,wgt/0.87);
  if(nJets >= 4 && ev.run == 256630)  runlumi_4j->Fill(27,wgt/0.95);
  if(nJets >= 4 && ev.run == 256673)  runlumi_4j->Fill(28,wgt/0.005);
  if(nJets >= 4 && ev.run == 256674)  runlumi_4j->Fill(29,wgt/0.09);
  if(nJets >= 4 && ev.run == 256675)  runlumi_4j->Fill(30,wgt/7.28);
  if(nJets >= 4 && ev.run == 256677)  runlumi_4j->Fill(31,wgt/15.58);
  if(nJets >= 4 && ev.run == 256729)  runlumi_4j->Fill(32,wgt/66.55);
  if(nJets >= 4 && ev.run == 256734)  runlumi_4j->Fill(33,wgt/7.2);
  if(nJets >= 4 && ev.run == 256801)  runlumi_4j->Fill(34,wgt/8.83);
  if(nJets >= 4 && ev.run == 256842)  runlumi_4j->Fill(35,wgt/0.02);
  if(nJets >= 4 && ev.run == 256843)  runlumi_4j->Fill(36,wgt/17.48);
*/

  TH1F *lepptH=leppt_2j_leading, *lepetaH=lepeta_2j_leading,*lepphiH=lepphi_2j_leading,*mtH=leptmass_2j_leading, *jetptH=jetpt_2j_leading, *jetetaH=jeteta_2j_leading, *jetcsvH=jetcsv_2j_leading, *numvtxH=numvertices_2j_leading, *metptH=metpt_2j_leading, *metphiH=metphi_2j_leading, *mettmassH=mettmass_2j_leading, *upjecH==upjec_2j_leading, *downjecH=downjec_2j_leading;
  TH2F *lumiH=lumi_2j_leading;
  if(nJets==3){
    lepptH=leppt_3j_leading; lepetaH=lepeta_3j_leading; lepphiH=lepphi_3j_leading; mtH=leptmass_3j_leading; jetptH=jetpt_3j_leading; jetetaH=jeteta_3j_leading; jetcsvH=jetcsv_3j_leading; numvtxH=numvertices_3j_leading; metptH=metpt_3j_leading; metphiH=metphi_3j_leading; mettmassH=mettmass_3j_leading; upjecH==upjec_3j_leading; downjecH=downjec_3j_leading;
lumiH=lumi_3j_leading;
  }
  if(nJets>=4){
    lepptH=leppt_4j_leading; lepetaH=lepeta_4j_leading; lepphiH=lepphi_4j_leading; mtH=leptmass_4j_leading; jetptH=jetpt_4j_leading; jetetaH=jeteta_4j_leading; jetcsvH=jetcsv_4j_leading; numvtxH=numvertices_4j_leading; metptH=metpt_4j_leading; metphiH=metphi_4j_leading; mettmassH=mettmass_4j_leading; upjecH==upjec_4j_leading; downjecH=downjec_4j_leading; lumiH=lumi_4j_leading;

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

  upjecH->Fill(up_jecpt);
  downjecH->Fill(down_jecpt);
  lumiH->Fill(ev.run,ev.nj);
  }
  
  
  //close file
  
  f->Close();
  
  //open output file
  TFile *fOut=TFile::Open(output+"/"+filename,"RECREATE");
  cutflow->Write();
  bjetcutflow->Write();
  runlumi_2j->Write();
  runlumi_3j->Write();
  runlumi_4j->Write();

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
  lumi_2j_leading->Write();

  upjec_2j_leading->Write();
  downjec_2j_leading->Write();

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
  lumi_3j_leading->Write();

  upjec_3j_leading->Write();
  downjec_3j_leading->Write();

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
  lumi_4j_leading->Write();

  upjec_4j_leading->Write();
  downjec_4j_leading->Write();

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
    ReadTree(files[i],output,chToSelect);
    }
}
