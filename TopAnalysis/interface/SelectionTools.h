#ifndef _Selection_Tools_h
#define _Selection_Tools_h

#include <vector>
#include <map>

#include "TString.h"

enum type { SingleMuon, SingleElectron,
              DoubleMuon, DoubleElectron,
              ElectronMuon, MC, None };

class Trigger{
 public:
  Trigger(int, int, bool debug=false);
  ~Trigger();
  void setMuonTrigger(int);
  void setElectronTrigger(int);
  void setDataType(TString);
  void addRequiredMuonTrigger(TString);
  void addRequiredMuonTrigger(std::vector<TString>);
  void addRequiredElectronTrigger(TString);
  void addRequiredElectronTrigger(std::vector<TString>);
  void deleteRequiredMuonTrigger(TString);
  void deleteRequiredElectronTrigger(TString);
  bool triggerFired(TString);
  bool muonFired();
  bool electronFired();
  bool isSingleMuonEvent();
  bool isSingleElectronEvent();
  bool isMuonFile();
  bool isElectronFile();

 private:
  bool debug_;
  bool isData_;
  type dataType_;
  std::map<TString, bool> triggers_;
  std::vector<TString> requiredMuonTriggers_;
  std::vector<TString> requiredElectronTriggers_;

};

#endif
