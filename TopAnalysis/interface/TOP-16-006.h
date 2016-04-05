#ifndef _top16006_h_
#define _top16006_h_

#include <TString.h>
#include <TGraph.h>
#include <TH1F.h>
#include "TLorentzVector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include <map>
#include <vector>

enum FlavourSplitting {NOFLAVOURSPLITTING=0, UDSGSPLITTING=1, CSPLITTING=4, BSPLITTING=5 };

FactorizedJetCorrector *getFactorizedJetEnergyCorrector(TString,bool);
float getLeptonEnergyScaleUncertainty(int l_id,float l_pt,float l_eta);
Float_t computeMT(TLorentzVector &a, TLorentzVector &b);
void RunTop16006(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts);
std::map<Int_t,Float_t> lumiPerRun();
std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt);
float getLeptonEnergyScaleUncertainty(int l_id,float l_pt,float l_eta);
#endif
