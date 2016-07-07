#ifndef _topwidth_h_
#define _topwidth_h_

#include "TH1.h"
#include "TString.h"
#include "TFile.h"

#include <map>

struct TopWidthEvent_t
{
  Int_t cat,nw,nl,nj,nt;
  Float_t weight[100];
  Float_t l_pt[2],l_eta[2],l_phi[2],l_m[2],l_les[2];
  Int_t l_id[2];
  Float_t gl_pt[2],gl_eta[2],gl_phi[2],gl_m[2];
  Int_t gl_id[2];
  Float_t j_pt[50],j_eta[50],j_phi[50],j_m[50],j_jer[50],j_jes[50];
  Int_t j_btag[50];
  Float_t gj_pt[50],gj_eta[50],gj_phi[50],gj_m[50];
  Int_t gj_flav[50],gj_hadflav[50];
  Float_t t_pt[10],t_eta[10],t_phi[10],t_m[10];
  Int_t t_id[10];
  Float_t met_pt,met_phi;
};

void createTopWidthEventTree(TTree *t,TopWidthEvent_t &twev);
void resetTopWidthEvent(TopWidthEvent_t &twev);
void addGenScanCounters(std::map<TString, TH1 *> &plotColl,TFile *fIn);
void RunTopWidth(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts,
		 TString era);
#endif
