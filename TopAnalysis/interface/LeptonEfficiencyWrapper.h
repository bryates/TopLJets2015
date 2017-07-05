#ifndef _LeptonEfficiencyWrapper_h_
#define _LeptonEfficiencyWrapper_h_

#include "TString.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <map>

#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"

typedef std::pair<float,float> EffCorrection_t;

class LeptonEfficiencyWrapper 
{
 public:
  LeptonEfficiencyWrapper(bool isData,TString era,TString runPeriod, bool debug=false);
  EffCorrection_t getTriggerCorrection(Leptons leptons);
  //EffCorrection_t getTriggerCorrection(std::vector<int> pdgId,std::vector<TLorentzVector> leptons);
  EffCorrection_t getOfflineCorrection(Particle lep, int nvtx);
  //EffCorrection_t getOfflineCorrection(int pdgId,float pt,float eta, TString runPeriod);
  ~LeptonEfficiencyWrapper();  
 private:
  void init(TString era,TString runPeriod);
  bool debug_;
  int era_;
  TString runPeriod_;
  std::map<TString,TH2 *> lepEffH_;
  std::map<TString,TGraphAsymmErrors *> lepEffGr_;
  std::vector<std::vector<float>> eeSFHigh;
  std::vector<std::vector<float>> eeSFLow;
};

#endif
