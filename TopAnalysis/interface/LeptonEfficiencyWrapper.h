#ifndef _LeptonEfficiencyWrapper_h_
#define _LeptonEfficiencyWrapper_h_

#include "TString.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <map>

typedef std::pair<float,float> EffCorrection_t;

class LeptonEfficiencyWrapper 
{
 public:
  LeptonEfficiencyWrapper(bool isData,TString era,TString runPeriod);
  EffCorrection_t getTriggerCorrection(std::vector<int> pdgId,std::vector<TLorentzVector> leptons);
  EffCorrection_t getOfflineCorrection(int pdgId,float pt,float eta, TString runPeriod);
  ~LeptonEfficiencyWrapper();  
 private:
  void init(TString era,TString runPeriod);
  int era_;
  std::map<TString,TH2 *> lepEffH_;
  std::map<TString,TGraphAsymmErrors *> lepEffGr_;
};

#endif
