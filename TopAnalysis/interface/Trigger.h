#ifndef _Trigger_h
#define _Trigger_h

#include <vector>
#include <map>

#include "TString.h"

enum triggerType { SingleMuon, SingleElectron,
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
  void addRequiredDoubleMuonTrigger(TString);
  void addRequiredDoubleMuonTrigger(std::vector<TString>);
  void addRequiredDoubleElectronTrigger(TString);
  void addRequiredDoubleElectronTrigger(std::vector<TString>);
  void addRequiredEMTrigger(TString);
  void addRequiredEMTrigger(std::vector<TString>);
  void deleteRequiredMuonTrigger(TString);
  void deleteRequiredElectronTrigger(TString);
  void deleteRequiredDoubleMuonTrigger(TString);
  bool triggerFired(TString);
  bool muonFired();
  bool electronFired();
  bool doubleMuonFired();
  bool doubleElectronFired();
  bool EMFired();
  bool isSingleMuonEvent();
  bool isSingleElectronEvent();
  bool isDoubleMuonEvent();
  bool isDoubleElectronEvent();
  bool isEMEvent();
  bool isMuonFile();
  bool isElectronFile();
  bool isDoubleMuonFile();
  bool isDoubleElectronFile();
  bool isEMFile();

 private:
  bool debug_;
  bool isData_;
  triggerType dataType_;
  std::map<TString, bool> triggers_;
  std::vector<TString> requiredMuonTriggers_;
  std::vector<TString> requiredElectronTriggers_;
  std::vector<TString> requiredDoubleMuonTriggers_;
  std::vector<TString> requiredDoubleElectronTriggers_;
  std::vector<TString> requiredEMTriggers_;

};

#endif
