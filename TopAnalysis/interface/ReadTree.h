#ifndef _readtree_h_
#define _readtree_h_

#include <TString.h>
#include <TGraph.h>
#include <TH1F.h>
#include "TLorentzVector.h"

#include <map>
#include <vector>

enum FlavourSplitting {NOFLAVOURSPLITTING=0, UDSGSPLITTING=1, CSPLITTING=4, BSPLITTING=5 };
enum GenWeightMode { NOGENWGT=0, GENWEIGHT=1 };
enum SystematicVariations { NOMINAL=0,
			    PU,                 
			    MUTRIGGER,          
			    MUEFF,              
			    MUSCALE,            
			    ETRIGGER,           
			    EEFF,               
			    ESCALE,             
			    BEFF,               
			    CEFF,               
			    LEFF,               
			    UMET,               
			    JER,                
			    JES_Absolute,
			    JES_HighPtExtra,
			    JES_SinglePionECAL,
			    JES_SinglePionHCAL,
			    JES_Time,
			    JES_RelativeJEREC1,
			    JES_RelativeJEREC2,
			    JES_RelativeJERHF,
			    JES_RelativePtBB,
			    JES_RelativePtEC1,
			    JES_RelativePtEC2,
			    JES_RelativePtHF,
			    JES_RelativeFSR,
			    JES_RelativeStatEC2,
			    JES_RelativeStatHF,
			    JES_PileUpDataMC,
			    JES_PileUpPtBB,
			    JES_PileUpPtEC,
			    JES_PileUpPtHF,
			    JES_PileUpBias,
			    JES_FlavorPureGluon,
			    JES_FlavorPureQuark,
			    JES_FlavorPureCharm,
			    JES_FlavorPureBottom,
			    QCD_MUF,
			    QCD_MUR,
			    QCD_MURMUF,
			    LASTSYST};

TString getSystematicsLabel(int varIdx);
Float_t computeMT(TLorentzVector &a, TLorentzVector &b);
void ReadTree(TString filename,
	      TString outname,
	      Int_t channelSelection=13,
	      Int_t chargeSelection=0,
	      TH1F *normH=0,
	      Bool_t isTTbar=false,
	      FlavourSplitting flavourSplitting=NOFLAVOURSPLITTING,
	      GenWeightMode genWgtMode=NOGENWGT,
	      Bool_t runSysts=false);
std::map<Int_t,Float_t> lumiPerRun();
std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt);
std::vector<float> getLeptonSelectionScaleFactor(int l_id,float l_pt,float l_eta,bool isData);

#endif
