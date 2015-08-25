#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include "UserCode/TopAnalysis/interface/MiniEvent.h"
 
#include <vector>
#include <iostream>
#include <algorithm>

bool sortBySignificance(std::pair<int,std::pair<float,float> > a, std::pair<int,std::pair<float,float> > b){
  if( a.second.first>0 || b.second.first>0 ) return (a.second.first>b.second.first);
  return (a.second.second>b.second.second);
}

void ReadTree(TString filename,TString output,bool useChPt,int minVtx, int maxVtx){

  gROOT->Reset();

  TH1F *cutflow = new TH1F("cutflow",";Cut;Events" ,6,0.,6.);
  cutflow->GetXaxis()->SetBinLabel(1,"preselected");
  cutflow->GetXaxis()->SetBinLabel(2,"#geq 2j");
  cutflow->GetXaxis()->SetBinLabel(3,"#geq 3j");
  cutflow->GetXaxis()->SetBinLabel(4,"#geq 4j");
  cutflow->GetXaxis()->SetBinLabel(5,"#geq 1b-tag");
  cutflow->GetXaxis()->SetBinLabel(6,"#geq 2b-tags");

  TH1F *bjetcutflow = new TH1F("bjetcutflow",";Cut;Events" ,9,0.,9.);
  bjetcutflow->GetXaxis()->SetBinLabel(1,"2j,=0b");
  bjetcutflow->GetXaxis()->SetBinLabel(2,"2j,=1b");
  bjetcutflow->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  bjetcutflow->GetXaxis()->SetBinLabel(4,"3j,=0b");
  bjetcutflow->GetXaxis()->SetBinLabel(5,"3j,=1b");
  bjetcutflow->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  bjetcutflow->GetXaxis()->SetBinLabel(7,"4j,=0b");
  bjetcutflow->GetXaxis()->SetBinLabel(8,"4j,=1b");
  bjetcutflow->GetXaxis()->SetBinLabel(9,"4j,#geq2b");
  
  TH1F *disvtx_2j_leading     = new TH1F("disvtx_2j_leading",";lxy;Events" ,25,0.,10.);
  TH1F *disvtx_2j_nextleading = (TH1F *)disvtx_2j_leading->Clone("disvtx_2j_nextleading");
  TH1F *disvtx_3j_leading     = (TH1F *)disvtx_2j_leading->Clone("disvtx_3j_leading");
  TH1F *disvtx_3j_nextleading = (TH1F *)disvtx_2j_leading->Clone("disvtx_3j_nextleading");
  TH1F *disvtx_4j_leading     = (TH1F *)disvtx_2j_leading->Clone("disvtx_4j_leading");
  TH1F *disvtx_4j_nextleading = (TH1F *)disvtx_2j_leading->Clone("disvtx_4j_nextleading");
  disvtx_2j_nextleading->SetName("disvtx_2j_nextleading");
  disvtx_3j_leading->SetName("disvtx_3j_leading");    
  disvtx_3j_nextleading->SetName("disvtx_3j_nextleading");
  disvtx_4j_leading->SetName("disvtx_4j_leading");
  disvtx_4j_nextleading->SetName("disvtx_4j_nextleading");
  
  TH1F *lxyz_sig_2j_leading     = new TH1F("lxyz_sig_2j_leading",";lxyz_sig;Events" ,25,0.,20.);
  TH1F *lxyz_sig_2j_nextleading = (TH1F *)lxyz_sig_2j_leading->Clone("lxyz_sig_2j_nextleading");
  TH1F *lxyz_sig_3j_leading     = (TH1F *)lxyz_sig_2j_leading->Clone("lxyz_sig_3j_leading");
  TH1F *lxyz_sig_3j_nextleading = (TH1F *)lxyz_sig_2j_leading->Clone("lxyz_sig_3j_nextleading");
  TH1F *lxyz_sig_4j_leading     = (TH1F *)lxyz_sig_2j_leading->Clone("lxyz_sig_4j_leading");
  TH1F *lxyz_sig_4j_nextleading = (TH1F *)lxyz_sig_2j_leading->Clone("lxyz_sig_4j_nextleading");
  lxyz_sig_2j_nextleading -> SetName("lxyz_sig_2j_nextleading");
  lxyz_sig_3j_leading     -> SetName("lxyz_sig_3j_leading");
  lxyz_sig_3j_nextleading -> SetName("lxyz_sig_3j_nextleading");
  lxyz_sig_4j_leading     -> SetName("lxyz_sig_4j_leading");
  lxyz_sig_4j_nextleading -> SetName("lxyz_sig_4j_nextleading");

  TH1F *vertexmass_2j_leading     = new TH1F("vertexmass_2j_leading",";vertexmass;Events" ,25,0.,6.);
  TH1F *vertexmass_2j_nextleading = (TH1F *)vertexmass_2j_leading->Clone("vertexmass_2j_nextleading");
  TH1F *vertexmass_3j_leading     = (TH1F *)vertexmass_2j_leading->Clone("vertexmass_3j_leading");
  TH1F *vertexmass_3j_nextleading = (TH1F *)vertexmass_2j_leading->Clone("vertexmass_3j_nextleading");
  TH1F *vertexmass_4j_leading     = (TH1F *)vertexmass_2j_leading->Clone("vertexmass_4j_leading");
  TH1F *vertexmass_4j_nextleading = (TH1F *)vertexmass_2j_leading->Clone("vertexmass_4j_nextleading");
  vertexmass_2j_nextleading -> SetName("vertexmass_2j_nextleading");
  vertexmass_3j_leading     -> SetName("vertexmass_3j_leading");
  vertexmass_3j_nextleading -> SetName("vertexmass_3j_nextleading");
  vertexmass_4j_leading     -> SetName("vertexmass_4j_leading");
  vertexmass_4j_nextleading -> SetName("vertexmass_4j_nextleading");

  TH1F *numvertices_2j_leading     = new TH1F("numvertices_2j_leading",";vertices;Events" ,100,0.,100.);
  TH1F *numvertices_2j_nextleading = (TH1F *)numvertices_2j_leading->Clone("numvertices_2j_nextleading");
  TH1F *numvertices_3j_leading     = (TH1F *)numvertices_2j_leading->Clone("numvertices_3j_leading");
  TH1F *numvertices_3j_nextleading = (TH1F *)numvertices_2j_leading->Clone("numvertices_3j_nextleading");
  TH1F *numvertices_4j_leading     = (TH1F *)numvertices_2j_leading->Clone("numvertices_4j_leading");
  TH1F *numvertices_4j_nextleading = (TH1F *)numvertices_2j_leading->Clone("numvertices_4j_nextleading");
  numvertices_2j_nextleading -> SetName("numvertices_2j_nextleading");
  numvertices_3j_leading     -> SetName("numvertices_3j_leading");   
  numvertices_3j_nextleading -> SetName("numvertices_3j_nextleading");
  numvertices_4j_leading     -> SetName("numvertices_4j_leading");
  numvertices_4j_nextleading -> SetName("numvertices_4j_nextleading");

  TH1F *jetpt_2j_leading     = new TH1F("jetpt_2j_leading",";pt;Events" ,50,0.,300.);
  TH1F *jetpt_2j_nextleading = (TH1F *)jetpt_2j_leading->Clone("jetpt_2j_nextleading");
  TH1F *jetpt_3j_leading     = (TH1F *)jetpt_2j_leading->Clone("jetpt_3j_leading");
  TH1F *jetpt_3j_nextleading = (TH1F *)jetpt_2j_leading->Clone("jetpt_3j_nextleading");
  TH1F *jetpt_4j_leading     = (TH1F *)jetpt_2j_leading->Clone("jetpt_4j_leading");
  TH1F *jetpt_4j_nextleading = (TH1F *)jetpt_2j_leading->Clone("jetpt_4j_nextleading");
  jetpt_2j_nextleading -> SetName("jetpt_2j_nextleading");
  jetpt_3j_leading     -> SetName("jetpt_3j_leading");   
  jetpt_3j_nextleading -> SetName("jetpt_3j_nextleading");
  jetpt_4j_leading     -> SetName("jetpt_4j_leading");
  jetpt_4j_nextleading -> SetName("jetpt_4j_nextleading");

  TH1F *jeteta_2j_leading     = new TH1F("jeteta_2j_leading",";eta;Events" ,25,0.,3.);
  TH1F *jeteta_2j_nextleading = (TH1F *)jeteta_2j_leading->Clone("jeteta_2j_nextleading");
  TH1F *jeteta_3j_leading     = (TH1F *)jeteta_2j_leading->Clone("jeteta_3j_leading");
  TH1F *jeteta_3j_nextleading = (TH1F *)jeteta_2j_leading->Clone("jeteta_3j_nextleading");
  TH1F *jeteta_4j_leading     = (TH1F *)jeteta_2j_leading->Clone("jeteta_4j_leading");
  TH1F *jeteta_4j_nextleading = (TH1F *)jeteta_2j_leading->Clone("jeteta_4j_nextleading");
  jeteta_2j_nextleading -> SetName("jeteta_2j_nextleading");
  jeteta_3j_leading     -> SetName("jeteta_3j_leading");   
  jeteta_3j_nextleading -> SetName("jeteta_3j_nextleading");
  jeteta_4j_leading     -> SetName("jeteta_4j_leading");
  jeteta_4j_nextleading -> SetName("jeteta_4j_nextleading");
  
  TH1F *lepCh_2j_leading     = new TH1F("lepCh_2j_leading",";ch;Events" ,6,-2.,2.);
  TH1F *lepCh_2j_nextleading = (TH1F *)lepCh_2j_leading->Clone("lepCh_2j_nextleading");
  TH1F *lepCh_3j_leading     = (TH1F *)lepCh_2j_leading->Clone("lepCh_3j_leading");
  TH1F *lepCh_3j_nextleading = (TH1F *)lepCh_2j_leading->Clone("lepCh_3j_nextleading");
  TH1F *lepCh_4j_leading     = (TH1F *)lepCh_2j_leading->Clone("lepCh_4j_leading");
  TH1F *lepCh_4j_nextleading = (TH1F *)lepCh_2j_leading->Clone("lepCh_4j_nextleading");
  lepCh_2j_nextleading -> SetName("lepCh_2j_nextleading");
  lepCh_3j_leading     -> SetName("lepCh_3j_leading");   
  lepCh_3j_nextleading -> SetName("lepCh_3j_nextleading");
  lepCh_4j_leading     -> SetName("lepCh_4j_leading");
  lepCh_4j_nextleading -> SetName("lepCh_4j_nextleading");
  
  TH1F *wmass_2j_leading     = new TH1F("wmass_2j_leading",";Wmass;Events" ,8,0.,300.);
  TH1F *wmass_2j_nextleading = (TH1F *)wmass_2j_leading->Clone("wmass_2j_nextleading");
  TH1F *wmass_3j_leading     = (TH1F *)wmass_2j_leading->Clone("wmass_3j_leading");
  TH1F *wmass_3j_nextleading = (TH1F *)wmass_2j_leading->Clone("wmass_3j_nextleading");
  TH1F *wmass_4j_leading     = (TH1F *)wmass_2j_leading->Clone("wmass_4j_leading");
  TH1F *wmass_4j_nextleading = (TH1F *)wmass_2j_leading->Clone("wmass_4j_nextleading");
  wmass_2j_nextleading -> SetName("wmass_2j_nextleading");
  wmass_3j_leading     -> SetName("wmass_3j_leading");   
  wmass_3j_nextleading -> SetName("wmass_3j_nextleading");
  wmass_4j_leading     -> SetName("wmass_4j_leading");
  wmass_4j_nextleading -> SetName("wmass_4j_nextleading");

  TH1F *wchmass_2j_leading     = new TH1F("wchmass_2j_leading",";chWmass;Events" ,8,0.,300.);
  TH1F *wchmass_2j_nextleading = (TH1F *)wchmass_2j_leading->Clone("wchmass_2j_nextleading");
  TH1F *wchmass_3j_leading     = (TH1F *)wchmass_2j_leading->Clone("wchmass_3j_leading");
  TH1F *wchmass_3j_nextleading = (TH1F *)wchmass_2j_leading->Clone("wchmass_3j_nextleading");
  TH1F *wchmass_4j_leading     = (TH1F *)wchmass_2j_leading->Clone("wchmass_4j_leading");
  TH1F *wchmass_4j_nextleading = (TH1F *)wchmass_2j_leading->Clone("wchmass_4j_nextleading");
  wchmass_2j_nextleading -> SetName("wchmass_2j_nextleading");
  wchmass_3j_leading     -> SetName("wchmass_3j_leading");   
  wchmass_3j_nextleading -> SetName("wchmass_3j_nextleading");
  wchmass_4j_leading     -> SetName("wchmass_4j_leading");
  wchmass_4j_nextleading -> SetName("wchmass_4j_nextleading");
  
  TH1F *neutiso_2j     = new TH1F("neutiso_2j",";neutiso;Events" ,8,0.,20.);
//  TH1F *neutiso_2j_nextleading = (TH1F *)neutiso_2j_leading->Clone("neutiso_2j_nextleading");
  TH1F *neutiso_3j     = (TH1F *)neutiso_2j->Clone("neutiso_3j");
//  TH1F *neutiso_3j_nextleading = (TH1F *)neutiso_2j_leading->Clone("neutiso_3j_nextleading");
  TH1F *neutiso_4j     = (TH1F *)neutiso_2j->Clone("neutiso_4j");
//  TH1F *neutiso_4j_nextleading = (TH1F *)neutiso_2j_leading->Clone("neutiso_4j_nextleading");

  TH1F *photoniso_2j     = new TH1F("photoniso_2j",";photoniso;Events" ,15,0.,30.);
//  TH1F *photoniso_2j_nextleading = (TH1F *)photoniso_2j_leading->Clone("photoniso_2j_nextleading");
  TH1F *photoniso_3j     = (TH1F *)photoniso_2j->Clone("photoniso_3j");
//  TH1F *photoniso_3j_nextleading = (TH1F *)photoniso_2j_leading->Clone("photoniso_3j_nextleading");
  TH1F *photoniso_4j     = (TH1F *)photoniso_2j->Clone("photoniso_4j");
//  TH1F *photoniso_4j_nextleading = (TH1F *)photoniso_2j_leading->Clone("photoniso_4j_nextleading");

  TH1F *sumiso_2j     = new TH1F("sumiso_2j",";sumiso;Events" ,8,0.,20.);
//  TH1F *sumiso_2j_nextleading = (TH1F *)sumiso_2j_leading->Clone("sumiso_2j_nextleading");
  TH1F *sumiso_3j     = (TH1F *)sumiso_2j->Clone("sumiso_3j");
//  TH1F *sumiso_3j_nextleading = (TH1F *)sumiso_2j_leading->Clone("sumiso_3j_nextleading");
  TH1F *sumiso_4j     = (TH1F *)sumiso_2j->Clone("sumiso_4j");
//  TH1F *sumiso_4j_nextleading = (TH1F *)sumiso_2j_leading->Clone("sumiso_4j_nextleading");

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

      //select jets
      uint32_t nJets(0), nBtags(0),nJets30(0) ;
      int Nvertices(ev.nvtx); 	
      std::vector< std::pair<int,std::pair<float,float> > > vlxyz_sig;
      for (int k=0; k<ev.nj;k++){
        //check pt and eta of this jet
        float pt  = ev.j_pt[k];
        float chPt  = ev.j_chpt[k];
        float eta = ev.j_eta[k];
        //int Nvertices = ev.nvtx;
        float csv = ev.j_csv[k];
        float lxyz=ev.j_vtx3DVal[k];
        float lxyz_sig= ev.j_vtx3DSig[k];
        //float mass_Wbos= ev.mt;
        //float chmass_Wbos= ev.chmt;

	  if (useChPt) {
	    if(chPt > 15 && fabs(eta) < 2.5 ){
		    nJets++;
		  if (chPt > 30) nJets30 ++;
		  if(lxyz>0) {
        nBtags++;
        std::pair<float,float> jetKinematics(lxyz_sig,chPt);
        vlxyz_sig.push_back( std::pair<int,std::pair<float,float> >(k,jetKinematics) );
	  	    }   
	      }
	    }
	  else {
	      if (pt > 30 && fabs(eta) < 2.5){
        nJets++; 
        nJets30++;	
		  //b-discriminator cut https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
      if (csv>0.890) {
		      nBtags++;
		      std::pair<float,float> jetKinematics(csv,pt);
		      vlxyz_sig.push_back( std::pair<int,std::pair<float,float> >(k,jetKinematics) );
		      }
		    }
	    }
	  }

  std::sort(vlxyz_sig.begin(),vlxyz_sig.end(),sortBySignificance);
	  
  //filter on the number of vertices
  if(Nvertices < minVtx || Nvertices >= maxVtx) continue;
  //if(lep_charge < lepCh_neg || lep_charge >= lepCh_pos) continue;
  
  //fill cutflow histos
  if(nJets>=2 && nJets30>0)               cutflow->Fill(1);
  if(nJets>=3 && nJets30>0)               cutflow->Fill(2);
  if(nJets>=4 && nJets30>0)               cutflow->Fill(3);
  if(nJets>=4 && nJets30>0 && nBtags >=1) cutflow->Fill(4);
  if(nJets>=4 && nJets30>0 && nBtags >=2) cutflow->Fill(5);
  if(nJets==2 && nJets30>0 && nBtags ==0) bjetcutflow->Fill(0);	
  if(nJets==2 && nJets30>0 && nBtags ==1) bjetcutflow->Fill(1);	
  if(nJets==2 && nJets30>0 && nBtags >=2) bjetcutflow->Fill(2);	
  if(nJets==3 && nJets30>0 && nBtags ==0) bjetcutflow->Fill(3);	
  if(nJets==3 && nJets30>0 && nBtags ==1) bjetcutflow->Fill(4);	
  if(nJets==3 && nJets30>0 && nBtags >=2) bjetcutflow->Fill(5);	
  if(nJets>=4 && nJets30>0 && nBtags ==0) bjetcutflow->Fill(6);
  if(nJets>=4 && nJets30>0 && nBtags ==1) bjetcutflow->Fill(7);
  if(nJets>=4 && nJets30>0 && nBtags >=2) bjetcutflow->Fill(8);
  
  //show b-tagging discriminator distributions
  for(size_t v=0; v<vlxyz_sig.size(); v++) {
	  if(v>1) break;
	  int k=vlxyz_sig[v].first;
	  float pt  = useChPt ?  ev.j_chpt[k] : ev.j_pt[k];
	  float eta = ev.j_eta[k];
	  float lxyz=ev.j_vtx3DVal[k];
	  float lxyz_sig= ev.j_vtx3DSig[k];
	  float vtxmass = ev.j_vtxmass[k]; 
	  int numvertex = ev.nvtx;
	  int lep_charge = ev.l_charge;
    float mass_Wbos= ev.mt;
    float chmass_Wbos= ev.chmt;
	  float neut_iso = ev.l_neutralHadronIso;
	  float photon_iso = ev.l_photonIso;
	  float sum_iso = ev.l_neutralHadronIso + ev.l_photonIso;
	  int lepton_id = ev.l_id;
	  //most significantly displaced vertex
	  if(v==0 && lepton_id == 13) {
	    if(nJets ==2 && nJets30>0){
	    if(lxyz>0){
          disvtx_2j_leading->Fill(lxyz);
          lxyz_sig_2j_leading->Fill(lxyz_sig);
          vertexmass_2j_leading->Fill(vtxmass);
          numvertices_2j_leading->Fill(numvertex);
          }
	      jetpt_2j_leading->Fill(pt);
	      jeteta_2j_leading->Fill(fabs(eta));
	      lepCh_2j_leading->Fill(lep_charge);
	      wmass_2j_leading->Fill(mass_Wbos);
	      wchmass_2j_leading->Fill(chmass_Wbos);
	    //neutiso_2j_leading->Fill(neut_iso);
	    //photoniso_2j_leading->Fill(photon_iso);
	    //sumiso_2j_leading->Fill(sum_iso);
	    }
	    
  if(nJets ==3 && nJets30>0){
    if(lxyz>0){
       disvtx_3j_leading->Fill(lxyz);
       lxyz_sig_3j_leading->Fill(lxyz_sig);
       vertexmass_3j_leading->Fill(vtxmass);
       numvertices_3j_leading->Fill(numvertex);
      }
    jetpt_3j_leading->Fill(pt);
    jeteta_3j_leading->Fill(fabs(eta));
    lepCh_3j_leading->Fill(lep_charge);
    wmass_3j_leading->Fill(mass_Wbos);
    wchmass_3j_leading->Fill(chmass_Wbos);
 // neutiso_3j_leading->Fill(neut_iso);
      //  photoniso_3j_leading->Fill(photon_iso);
      //  sumiso_3j_leading->Fill(sum_iso);
	    }

	    if(nJets >=4 && nJets30>0){
	      if(lxyz>0){
          disvtx_4j_leading->Fill(lxyz);
          lxyz_sig_4j_leading->Fill(lxyz_sig);
          vertexmass_4j_leading->Fill(vtxmass);
          numvertices_4j_leading->Fill(numvertex);
          }
        
        jetpt_4j_leading->Fill(pt);
        jeteta_4j_leading->Fill(fabs(eta));
	      lepCh_4j_leading->Fill(lep_charge);
	      wmass_4j_leading->Fill(mass_Wbos);
        wchmass_4j_leading->Fill(chmass_Wbos);
//	    neutiso_4j_leading->Fill(neut_iso);
//      photoniso_4j_leading->Fill(photon_iso);
//      sumiso_4j_leading->Fill(sum_iso);
	    }
	  }
	  
	  //second most significantly displaced vertex
	  if(v==1 && lepton_id == 13) {
	    if(nJets ==2 && nJets30>0){
	      if(lxyz>0){
		      disvtx_2j_nextleading->Fill(lxyz);
		      lxyz_sig_2j_nextleading->Fill(lxyz_sig);
		      vertexmass_2j_nextleading->Fill(vtxmass);
		      numvertices_2j_nextleading->Fill(numvertex);
	        }
	      jetpt_2j_nextleading->Fill(pt);
	      jeteta_2j_nextleading->Fill(fabs(eta));
        lepCh_2j_nextleading->Fill(lep_charge);
	      wmass_2j_nextleading->Fill(mass_Wbos);
        wchmass_2j_nextleading->Fill(chmass_Wbos);
//	    neutiso_2j_nextleading->Fill(neut_iso);
//      photoniso_2j_nextleading->Fill(photon_iso);
//      sumiso_2j_nextleading->Fill(sum_iso);
	    }

   if(nJets ==3 && nJets30>0){
   if(lxyz>0){
      disvtx_3j_nextleading->Fill(lxyz);
      lxyz_sig_3j_nextleading->Fill(lxyz_sig);
      vertexmass_3j_nextleading->Fill(vtxmass);
		  numvertices_3j_nextleading->Fill(numvertex);
      }
     jetpt_3j_nextleading->Fill(pt);
     jeteta_3j_nextleading->Fill(fabs(eta));
	   lepCh_3j_nextleading->Fill(lep_charge);
	   wmass_3j_nextleading->Fill(mass_Wbos);
     wchmass_3j_nextleading->Fill(chmass_Wbos);
//   neutiso_3j_nextleading->Fill(neut_iso);
//   photoniso_3j_nextleading->Fill(photon_iso);
//   sumiso_3j_nextleading->Fill(sum_iso);
	   }
	    
	   if(nJets >=4 && nJets30>0){
	      if(lxyz>0){
        disvtx_4j_nextleading->Fill(lxyz);
        lxyz_sig_4j_nextleading->Fill(lxyz_sig);
        vertexmass_4j_nextleading->Fill(vtxmass);
        numvertices_4j_nextleading->Fill(numvertex);
	      }
	      jetpt_4j_nextleading->Fill(pt);
	      jeteta_4j_nextleading->Fill(fabs(eta));
	      lepCh_4j_nextleading->Fill(lep_charge);
	      wmass_4j_nextleading->Fill(mass_Wbos);
        wchmass_4j_nextleading->Fill(chmass_Wbos);
//      neutiso_4j_nextleading->Fill(neut_iso);
//      photoniso_4j_nextleading->Fill(photon_iso);
//      sumiso_4j_nextleading->Fill(sum_iso);
	      }
	    }
      if(nJets == 2 && nJets30 > 0 && lepton_id == 13) neutiso_2j->Fill(neut_iso);
      if(nJets == 2 && nJets30 > 0 && lepton_id == 13) photoniso_2j->Fill(photon_iso);
      if(nJets == 2 && nJets30 > 0 && lepton_id == 13) sumiso_2j->Fill(sum_iso);

      if(nJets == 3 && nJets30 > 0 && lepton_id == 13) neutiso_3j->Fill(neut_iso);
      if(nJets == 3 && nJets30 > 0 && lepton_id == 13) photoniso_3j->Fill(photon_iso);
      if(nJets == 3 && nJets30 > 0 && lepton_id == 13) sumiso_3j->Fill(sum_iso);

      if(nJets >= 4 && nJets30 > 0 && lepton_id == 13) neutiso_4j->Fill(neut_iso);
      if(nJets >= 4 && nJets30 > 0 && lepton_id == 13) photoniso_4j->Fill(photon_iso);
      if(nJets >= 4 && nJets30 > 0 && lepton_id == 13) sumiso_4j->Fill(sum_iso);
        }
      }
      
  
  //close file
  
  f->Close();
  
  //open output file
  TFile *fOut=TFile::Open(output+"/"+filename,"RECREATE");
  cutflow->Write();
  bjetcutflow->Write();

  disvtx_2j_leading->Write();
  lxyz_sig_2j_leading->Write();
  vertexmass_2j_leading->Write();
  numvertices_2j_leading->Write();
  jetpt_2j_leading->Write();
  jeteta_2j_leading->Write();
  lepCh_2j_leading->Write();
  wmass_2j_leading->Write();
  wchmass_2j_leading->Write();
  neutiso_2j->Write();
  photoniso_2j->Write();
  sumiso_2j->Write();

  disvtx_3j_leading->Write();
  lxyz_sig_3j_leading->Write();
  vertexmass_3j_leading->Write();
  numvertices_3j_leading->Write();
  jetpt_3j_leading->Write();
  jeteta_3j_leading->Write();
  lepCh_3j_leading->Write();
  wmass_3j_leading->Write();
  wchmass_3j_leading->Write();
  neutiso_3j->Write();
  photoniso_3j->Write();
  sumiso_3j->Write();

  disvtx_4j_leading->Write();
  lxyz_sig_4j_leading->Write();
  vertexmass_4j_leading->Write();
  numvertices_4j_leading->Write();
  jetpt_4j_leading->Write();
  jeteta_4j_leading->Write();
  lepCh_4j_leading->Write();
  wmass_4j_leading->Write();
  wchmass_4j_leading->Write();
  neutiso_4j->Write();
  photoniso_4j->Write();
  sumiso_4j->Write();

  disvtx_2j_nextleading->Write();
  lxyz_sig_2j_leading->Write();
  vertexmass_2j_nextleading->Write();
  numvertices_2j_nextleading->Write();
  jetpt_2j_nextleading->Write();
  jeteta_2j_nextleading->Write();
  lepCh_2j_nextleading->Write();
  wmass_2j_nextleading->Write();
  wchmass_2j_nextleading->Write();
//  neutiso_2j_nextleading->Write();
//  photoniso_2j_nextleading->Write();
//  sumiso_2j_nextleading->Write();

  disvtx_3j_nextleading->Write();
  lxyz_sig_3j_leading->Write();
  vertexmass_3j_nextleading->Write();
  numvertices_3j_nextleading->Write();
  jetpt_3j_nextleading->Write();
  jeteta_3j_nextleading->Write();
  lepCh_3j_nextleading->Write();
  wmass_3j_nextleading->Write();
  wchmass_3j_nextleading->Write();
//  neutiso_3j_nextleading->Write();
//  photoniso_3j_nextleading->Write();
//  sumiso_3j_nextleading->Write();

  disvtx_4j_nextleading->Write();
  lxyz_sig_4j_nextleading->Write();
  vertexmass_4j_nextleading->Write();
  numvertices_4j_nextleading->Write();
  jetpt_4j_nextleading->Write();
  jeteta_4j_nextleading->Write();
  lepCh_4j_nextleading->Write();
  wmass_4j_nextleading->Write();
  wchmass_4j_nextleading->Write();
//  neutiso_4j_nextleading->Write();
//  photoniso_4j_nextleading->Write();
//  sumiso_4j_nextleading->Write();

 
  fOut->Close();
}


void RunOverSamples(TString output, bool useChPt,int minVtx, int maxVtx){
  TString files[]={
  "DYJetsToLL_M50_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8.root",
  "ST_tW_antitop_5f_DS_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M1.root",
  "ST_tW_top_5f_DS_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M1.root",
  "ST_t_channel_4f_leptonDecays_13TeV_amcatnlo_pythia8_TuneCUETP8M1.root",
  "WJetsToLNu_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8.root"
  };
  
  
  for(size_t i=0; i<sizeof(files)/sizeof(TString); i++){
    ReadTree(files[i],output,useChPt,minVtx,maxVtx);
    }
}
