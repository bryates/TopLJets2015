#ifndef _correction_tools_h_
#define _correction_tools_h_

#include <vector>
#include <string>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "TGraphAsymmErrors.h"

//pileup weighting
std::vector<TGraph *> getPileupWeights(TString era,TH1 *genPU);

//apply jet energy resolutions
MiniEvent_t smearJetEnergies(MiniEvent_t ev, std::string option = "central");

//see working points in https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco
MiniEvent_t addBTagDecisions(MiniEvent_t ev,float wp=0.8484);

//details in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
MiniEvent_t updateBTagDecisions(MiniEvent_t ev, 
				std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> &btvsfReaders,
				std::map<BTagEntry::JetFlavor, TGraphAsymmErrors*> &expBtagEff, 
				std::map<BTagEntry::JetFlavor, TGraphAsymmErrors*> &expBtagEffPy8, 
				BTagSFUtil *myBTagSFUtil, 
				std::string option = "central");

//details in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> getBTVcalibrationReaders(TString era,BTagEntry::OperatingPoint btagOP=BTagEntry::OP_MEDIUM);

//the expections are created with the script scripts/saveExpectedBtagEff.py (cf README)
std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> readExpectedBtagEff(TString era,TString btagExpPostFix="");

#endif
