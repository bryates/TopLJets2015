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

// QCD Scale
  TH1F *bjetcutflow_qcdScaleLo = new TH1F("bjetcutflow_qcdScaleLo","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_qcdScaleLo->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_qcdScaleLo->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_qcdScaleLo->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_qcdScaleLo->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_qcdScaleLo->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_qcdScaleLo->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_qcdScaleLo->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_qcdScaleLo->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_qcdScaleLo->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

  TH1F *bjetcutflow_qcdScaleHi = new TH1F("bjetcutflow_qcdScaleHi","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_qcdScaleHi->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_qcdScaleHi->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_qcdScaleHi->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_qcdScaleHi->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_qcdScaleHi->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_qcdScaleHi->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_qcdScaleHi->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_qcdScaleHi->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_qcdScaleHi->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

// HDAMP Scale
  TH1F *bjetcutflow_hdampLo = new TH1F("bjetcutflow_hdampLo","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_hdampLo->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_hdampLo->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_hdampLo->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_hdampLo->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_hdampLo->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_hdampLo->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_hdampLo->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_hdampLo->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_hdampLo->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

  TH1F *bjetcutflow_hdampHi = new TH1F("bjetcutflow_hdampHi","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_hdampHi->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_hdampHi->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_hdampHi->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_hdampHi->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_hdampHi->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_hdampHi->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_hdampHi->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_hdampHi->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_hdampHi->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

// JES
  TH1F *bjetcutflow_jesLo = new TH1F("bjetcutflow_jesLo","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_jesLo->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_jesLo->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_jesLo->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_jesLo->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_jesLo->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_jesLo->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_jesLo->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_jesLo->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_jesLo->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

  TH1F *bjetcutflow_jesHi = new TH1F("bjetcutflow_jesHi","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_jesHi->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_jesHi->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_jesHi->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_jesHi->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_jesHi->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_jesHi->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_jesHi->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_jesHi->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_jesHi->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

// JER
  TH1F *bjetcutflow_jerLo = new TH1F("bjetcutflow_jerLo","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_jerLo->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_jerLo->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_jerLo->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_jerLo->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_jerLo->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_jerLo->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_jerLo->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_jerLo->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_jerLo->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

  TH1F *bjetcutflow_jerHi = new TH1F("bjetcutflow_jerHi","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_jerHi->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_jerHi->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_jerHi->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_jerHi->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_jerHi->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_jerHi->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_jerHi->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_jerHi->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_jerHi->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

// bTag
  TH1F *bjetcutflow_btagLo = new TH1F("bjetcutflow_btagLo","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_btagLo->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_btagLo->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_btagLo->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_btagLo->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_btagLo->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_btagLo->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_btagLo->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_btagLo->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_btagLo->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

  TH1F *bjetcutflow_btagHi = new TH1F("bjetcutflow_btagHi","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_btagHi->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_btagHi->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_btagHi->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_btagHi->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_btagHi->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_btagHi->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_btagHi->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_btagHi->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_btagHi->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

// misTag
  TH1F *bjetcutflow_mistagLo = new TH1F("bjetcutflow_mistagLo","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_mistagLo->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_mistagLo->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_mistagLo->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_mistagLo->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_mistagLo->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_mistagLo->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_mistagLo->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_mistagLo->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_mistagLo->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

  TH1F *bjetcutflow_mistagHi = new TH1F("bjetcutflow_mistagHi","bjet cutflow ;Cut;Events" ,9,0.,9.);
  bjetcutflow_mistagHi->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow_mistagHi->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow_mistagHi->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow_mistagHi->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow_mistagHi->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow_mistagHi->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow_mistagHi->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow_mistagHi->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow_mistagHi->GetXaxis()->SetBinLabel(9,"4j,#geq2b");

// Lepton pt 
  TH1F *leppt_2j = new TH1F("leppt_2j",";ch;Events" ,20,0.,300.);
  TH1F *leppt_3j = (TH1F *)leppt_2j->Clone("leppt_3j");
  TH1F *leppt_4j = (TH1F *)leppt_2j->Clone("leppt_4j");
  leppt_3j   ->SetName("leppt_3j");
  leppt_4j   ->SetName("leppt_4j");

// Lepton eta
  TH1F *lepeta_2j = new TH1F("lepeta_2j",";ch;Events" ,12,0.,3.);
  TH1F *lepeta_3j = (TH1F *)lepeta_2j->Clone("lepeta_3j");
  TH1F *lepeta_4j = (TH1F *)lepeta_2j->Clone("lepeta_4j");
  lepeta_3j   ->SetName("lepeta_3j");
  lepeta_4j   ->SetName("lepeta_4j");

// Lepton phi
  TH1F *lepphi_2j = new TH1F("lepphi_2j",";ch;Events" ,50,-3.2,3.2);
  TH1F *lepphi_3j = (TH1F *)lepphi_2j->Clone("lepphi_3j");
  TH1F *lepphi_4j = (TH1F *)lepphi_2j->Clone("lepphi_4j");
  lepphi_3j   ->SetName("lepphi_3j");
  lepphi_4j   ->SetName("lepphi_4j");

// Lepton transverse mass
  TH1F *leptmass_2j = new TH1F("leptmass_2j",";Transverse Mass;Events" ,100,0.,200.);
  TH1F *leptmass_3j = (TH1F *)leptmass_2j->Clone("leptmass_3j");
  TH1F *leptmass_4j = (TH1F *)leptmass_2j->Clone("leptmass_4j");
  leptmass_3j   ->SetName("leptmass_3j");
  leptmass_4j   ->SetName("leptmass_4j");

// Jet pt
  TH1F *jetpt_2j = new TH1F("jetpt_2j",";pt;Events" ,20,0.,300.);
  TH1F *jetpt_3j = (TH1F *)jetpt_2j->Clone("jetpt_3j");
  TH1F *jetpt_4j = (TH1F *)jetpt_2j->Clone("jetpt_4j");
  jetpt_3j   ->SetName("jetpt_3j");
  jetpt_4j   ->SetName("jetpt_4j");

// Jet eta  
  TH1F *jeteta_2j = new TH1F("jeteta_2j",";eta;Events" ,12,0.,3.);
  TH1F *jeteta_3j = (TH1F *)jeteta_2j->Clone("jeteta_3j");
  TH1F *jeteta_4j = (TH1F *)jeteta_2j->Clone("jeteta_4j");
  jeteta_3j   ->SetName("jeteta_3j");
  jeteta_4j   ->SetName("jeteta_4j");

// CSV
  TH1F *jetcsv_2j = new TH1F("jetcsv_2j",";csv;Events" ,100,-1.2,1.2);
  TH1F *jetcsv_3j     = (TH1F *)jetcsv_2j->Clone("jetcsv_3j");
  TH1F *jetcsv_4j     = (TH1F *)jetcsv_2j->Clone("jetcsv_4j");
  jetcsv_3j     -> SetName("jetcsv_3j");
  jetcsv_4j     -> SetName("jetcsv_4j");

// numvertices
  TH1F *numvertices_2j = new TH1F("numvertices_2j",";vertices;Events" ,25,0.,50.);
  TH1F *numvertices_3j     = (TH1F *)numvertices_2j->Clone("numvertices_3j");
  TH1F *numvertices_4j     = (TH1F *)numvertices_2j->Clone("numvertices_4j");
  numvertices_3j     -> SetName("numvertices_3j");
  numvertices_4j     -> SetName("numvertices_4j");

// MET pt
  TH1F *metpt_2j = new TH1F("metpt_2j",";ch;Events" ,20,0.,300.);
  TH1F *metpt_3j     = (TH1F *)metpt_2j->Clone("metpt_3j");
  TH1F *metpt_4j     = (TH1F *)metpt_2j->Clone("metpt_4j");
  metpt_3j     -> SetName("metpt_3j");
  metpt_4j     -> SetName("metpt_4j");

// MET phi
  TH1F *metphi_2j = new TH1F("metphi_2j",";ch;Events" ,50,-3.2,3.2);
  TH1F *metphi_3j     = (TH1F *)metphi_2j->Clone("metphi_3j");
  TH1F *metphi_4j     = (TH1F *)metphi_2j->Clone("metphi_4j");
  metphi_3j     -> SetName("metphi_3j");
  metphi_4j     -> SetName("metphi_4j");

// MET transverse mass
  TH1F *mettmass_2j = new TH1F("mettmass_2j",";Transverse Mass;Events" ,100,0.,200.);
  TH1F *mettmass_3j     = (TH1F *)mettmass_2j->Clone("mettmass_3j");
  TH1F *mettmass_4j     = (TH1F *)mettmass_2j->Clone("mettmass_4j");
  mettmass_3j     -> SetName("mettmass_3j");
  mettmass_4j     -> SetName("mettmass_4j");

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
  if(abs(lepton_id) == 13 && !ev.muTrigger) continue;
  if(abs(lepton_id) == 11 && !ev.elTrigger) continue;

// JES
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
//  iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
//  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

   JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("Jec11_V12_Uncertainty_AK5PF.txt");

   jecUnc->setJetEta(fabs(ev.j_eta[0]));
   jecUnc->setJetPt(ev.j_pt[0]);
   double unc = jecUnc->getUncertainty(true);
   double jesHi = (ev.j_pt[0])*(1+unc);
   double jesLo = (ev.j_pt[0])*(1-unc);

// b-tagging 
   const reco::Candidate *genjet = j.genParton();
        bool bjet(false),cjet(false),lightjets(false);
        if(genjet){
  bjet=(abs(genjet->pdgId())==5);
        cjet=(abs(genjet->pdgId())==4);
  //for u, d, s, g 
        lightjets=(abs(genjet->pdgId())== 1 || abs(genjet->pdgId())==2 || abs(genjet->pdgId())==3 || abs(genjet->pdgId())==21);
        }

  math::XYZTLorentzVector vvv = j.p4();
  if (doJER && !Data) {
    (j.genJet() && vvv.pt()>10)? JERCor = getJERfactor(j.pt(), j.eta(), j.genJet()->pt()): JERCor = 1.0;
vvv *= JERCor;
jer_pt = vvv.pt(), jer_eta = vvv.eta(), jer_phi = vvv.phi();
}
if (doJER && !Data) {
(j.genJet() && vvv.pt()>10)? JERCor_UP = getJERfactor_up(j.pt(), j.eta(), j.genJet()->pt()): JERCor_UP = 1.0;
vvv *= JERCor_UP;
jerup_pt = vvv.pt(), jerup_eta = vvv.eta(), jerup_phi = vvv.phi();
}
if (doJER && !Data) {
(j.genJet() && vvv.pt()>10)? JERCor_DOWN = getJERfactor_down(j.pt(), j.eta(), j.genJet()->pt()): JERCor_DOWN = 1.0;
vvv *= JERCor_DOWN;
jerdown_pt = vvv.pt(), jerdown_eta = vvv.eta(), jerdown_phi = vvv.phi();
}

for( pat::JetCollection::const_iterator jetIter = jets->begin(); jetIter != jets->end(); ++jetIter ){

//        float csvm = jetIter->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
math::XYZTLorentzVector vvv = jetIter->p4();
if (doJER && !Data) {
// get the JER correction factor for this jet
(jetIter->genJet() && vvv.pt()>10)? JERCor = getJERfactor(jetIter->pt(), jetIter->eta(), jetIter->genJet()->pt()): JERCor = 1.0;
vvv *= JERCor;
jer_pt = vvv.pt(), jer_eta = vvv.eta(), jer_phi = vvv.phi();
}
if (doJER && !Data) {
// get the JER correction factor for this jet
(jetIter->genJet() && vvv.pt()>10)? JERCor_UP = getJERfactor_up(jetIter->pt(), jetIter->eta(), jetIter->genJet()->pt()): JERCor_UP = 1.0;
vvv *= JERCor_UP;
jerup_pt = vvv.pt(), jerup_eta = vvv.eta(), jerup_phi = vvv.phi();
}
if (doJER && !Data) {
// get the JER correction factor for this jet
(jetIter->genJet() && vvv.pt()>10)? JERCor_DOWN = getJERfactor_down(jetIter->pt(), jetIter->eta(), jetIter->genJet()->pt()): JERCor_DOWN = 1.0;
vvv *= JERCor_DOWN;
jerdown_pt = vvv.pt(), jerdown_eta = vvv.eta(), jerdown_phi = vvv.phi();
}

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

// QCD Scale
  if(nJets == 2 && nBtags == 0 ) bjetcutflow_qcdScaleLo->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_qcdScaleLo->Fill(1,wgt*ttbar_w[9]/ttbar_w[0]);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_qcdScaleLo->Fill(2,wgt*ttbar_w[9]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_qcdScaleLo->Fill(3,wgt*ttbar_w[9]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_qcdScaleLo->Fill(4,wgt*ttbar_w[9]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_qcdScaleLo->Fill(5,wgt*ttbar_w[9]/ttbar_w[0]);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_qcdScaleLo->Fill(6,wgt*ttbar_w[9]/ttbar_w[0]);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_qcdScaleLo->Fill(7,wgt*ttbar_w[9]/ttbar_w[0]);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_qcdScaleLo->Fill(8,wgt*ttbar_w[9]/ttbar_w[0]);

  if(nJets == 2 && nBtags == 0 ) bjetcutflow_qcdScaleHi->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_qcdScaleHi->Fill(1,wgt*ttbar_w[5]/ttbar_w[0]);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_qcdScaleHi->Fill(2,wgt*ttbar_w[5]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_qcdScaleHi->Fill(3,wgt*ttbar_w[5]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_qcdScaleHi->Fill(4,wgt*ttbar_w[5]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_qcdScaleHi->Fill(5,wgt*ttbar_w[5]/ttbar_w[0]);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_qcdScaleHi->Fill(6,wgt*ttbar_w[5]/ttbar_w[0]);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_qcdScaleHi->Fill(7,wgt*ttbar_w[5]/ttbar_w[0]);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_qcdScaleHi->Fill(8,wgt*ttbar_w[5]/ttbar_w[0]);

// HDAMP
  if(nJets == 2 && nBtags == 0 ) bjetcutflow_hdampLo->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_hdampLo->Fill(1,wgt*ttbar_w[ttbar_nw-17]/ttbar_w[0]);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_hdampLo->Fill(2,wgt*ttbar_w[ttbar_nw-17]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_hdampLo->Fill(3,wgt*ttbar_w[ttbar_nw-17]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_hdampLo->Fill(4,wgt*ttbar_w[ttbar_nw-17]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_hdampLo->Fill(5,wgt*ttbar_w[ttbar_nw-17]/ttbar_w[0]);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_hdampLo->Fill(6,wgt*ttbar_w[ttbar_nw-17]/ttbar_w[0]);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_hdampLo->Fill(7,wgt*ttbar_w[ttbar_nw-17]/ttbar_w[0]);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_hdampLo->Fill(8,wgt*ttbar_w[ttbar_nw-17]/ttbar_w[0]);

  if(nJets == 2 && nBtags == 0 ) bjetcutflow_hdampHi->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_hdampHi->Fill(1,wgt*ttbar_w[ttbar_nw-9]/ttbar_w[0]);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_hdampHi->Fill(2,wgt*ttbar_w[ttbar_nw-9]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_hdampHi->Fill(3,wgt*ttbar_w[ttbar_nw-9]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_hdampHi->Fill(4,wgt*ttbar_w[ttbar_nw-9]/ttbar_w[0]);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_hdampHi->Fill(5,wgt*ttbar_w[ttbar_nw-9]/ttbar_w[0]);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_hdampHi->Fill(6,wgt*ttbar_w[ttbar_nw-9]/ttbar_w[0]);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_hdampHi->Fill(7,wgt*ttbar_w[ttbar_nw-9]/ttbar_w[0]);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_hdampHi->Fill(8,wgt*ttbar_w[ttbar_nw-9]/ttbar_w[0]);

// JES
  if(nJets == 2 && nBtags == 0 ) bjetcutflow_jesLo->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_jesLo->Fill(1,wgt);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_jesLo->Fill(2,wgt);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_jesLo->Fill(3,wgt);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_jesLo->Fill(4,wgt);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_jesLo->Fill(5,wgt);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_jesLo->Fill(6,wgt);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_jesLo->Fill(7,wgt);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_jesLo->Fill(8,wgt);

  if(nJets == 2 && nBtags == 0 ) bjetcutflow_jesHi->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_jesHi->Fill(1,wgt);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_jesHi->Fill(2,wgt);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_jesHi->Fill(3,wgt);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_jesHi->Fill(4,wgt);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_jesHi->Fill(5,wgt);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_jesHi->Fill(6,wgt);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_jesHi->Fill(7,wgt);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_jesHi->Fill(8,wgt);

// JER
  if(nJets == 2 && nBtags == 0 ) bjetcutflow_jerLo->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_jerLo->Fill(1,wgt);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_jerLo->Fill(2,wgt);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_jerLo->Fill(3,wgt);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_jerLo->Fill(4,wgt);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_jerLo->Fill(5,wgt);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_jerLo->Fill(6,wgt);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_jerLo->Fill(7,wgt);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_jerLo->Fill(8,wgt);

  if(nJets == 2 && nBtags == 0 ) bjetcutflow_jerHi->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_jerHi->Fill(1,wgt);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_jerHi->Fill(2,wgt);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_jerHi->Fill(3,wgt);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_jerHi->Fill(4,wgt);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_jerHi->Fill(5,wgt);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_jerHi->Fill(6,wgt);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_jerHi->Fill(7,wgt);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_jerHi->Fill(8,wgt);

// bTag
  if(nJets == 2 && nBtags == 0 ) bjetcutflow_btagLo->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_btagLo->Fill(1,wgt);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_btagLo->Fill(2,wgt);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_btagLo->Fill(3,wgt);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_btagLo->Fill(4,wgt);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_btagLo->Fill(5,wgt);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_btagLo->Fill(6,wgt);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_btagLo->Fill(7,wgt);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_btagLo->Fill(8,wgt);

  if(nJets == 2 && nBtags == 0 ) bjetcutflow_btagHi->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_btagHi->Fill(1,wgt);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_btagHi->Fill(2,wgt);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_btagHi->Fill(3,wgt);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_btagHi->Fill(4,wgt);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_btagHi->Fill(5,wgt);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_btagHi->Fill(6,wgt);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_btagHi->Fill(7,wgt);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_btagHi->Fill(8,wgt);

// misTag
  if(nJets == 2 && nBtags == 0 ) bjetcutflow_mistagLo->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_mistagLo->Fill(1,wgt);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_mistagLo->Fill(2,wgt);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_mistagLo->Fill(3,wgt);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_mistagLo->Fill(4,wgt);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_mistagLo->Fill(5,wgt);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_mistagLo->Fill(6,wgt);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_mistagLo->Fill(7,wgt);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_mistagLo->Fill(8,wgt);

  if(nJets == 2 && nBtags == 0 ) bjetcutflow_mistagHi->Fill(0);	
  if(nJets == 2 && nBtags == 1 ) bjetcutflow_mistagHi->Fill(1,wgt);	
  if(nJets == 2 && nBtags == 2 ) bjetcutflow_mistagHi->Fill(2,wgt);	
  if(nJets >= 3 && nBtags == 0 ) bjetcutflow_mistagHi->Fill(3,wgt);	
  if(nJets >= 3 && nBtags == 1 ) bjetcutflow_mistagHi->Fill(4,wgt);	
  if(nJets >= 3 && nBtags >= 2 ) bjetcutflow_mistagHi->Fill(5,wgt);	
  if(nJets >= 4 && nBtags == 0 ) bjetcutflow_mistagHi->Fill(6,wgt);
  if(nJets >= 4 && nBtags == 1 ) bjetcutflow_mistagHi->Fill(7,wgt);
  if(nJets >= 4 && nBtags >= 2 ) bjetcutflow_mistagHi->Fill(8,wgt);


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
  bjetcutflow_qcdScaleLo->Write();
  bjetcutflow_qcdScaleHi->Write();
  bjetcutflow_hdampLo->Write();
  bjetcutflow_hdampHi->Write();
  bjetcutflow_jesLo->Write();
  bjetcutflow_jesHi->Write();
  bjetcutflow_jerLo->Write();
  bjetcutflow_jerHi->Write();
  bjetcutflow_btagLo->Write();
  bjetcutflow_btagHi->Write();
  bjetcutflow_mistagLo->Write();
  bjetcutflow_mistagHi->Write();

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
    ReadTree(files[i],output,chToSelect);
    }
}
