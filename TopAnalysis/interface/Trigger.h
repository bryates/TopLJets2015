#ifndef _Trigger_h
#define _Trigger_h

#include <vector>
#include <map>

#include "TString.h"

#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"

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
  bool triggerHasDZ(TString); //Check if a DZ version exists
  bool muonFired();
  bool electronFired();
  bool doubleMuonFired();
  bool doubleElectronFired();
  bool EMFired();
  bool isSingleMuonEvent();
  bool isSingleMuonEvent(Leptons);
  bool isSingleElectronEvent();
  bool isSingleElectronEvent(Leptons);
  bool isDoubleMuonEvent();
  bool isDoubleMuonEvent(Leptons);
  bool isDoubleElectronEvent();
  bool isDoubleElectronEvent(Leptons);
  bool isEMEvent();
  bool isEMEvent(Leptons);
  bool isMuonFile();
  bool isElectronFile();
  bool isDoubleMuonFile();
  bool isDoubleElectronFile();
  bool isEMFile();
  bool isMCFile();

 private:
  bool debug_;
  bool isData_;
  bool isH_;
  triggerType dataType_;
  std::map<TString, bool> triggers_;
  std::vector<TString> requiredMuonTriggers_;
  std::vector<TString> requiredElectronTriggers_;
  std::vector<TString> requiredDoubleMuonTriggers_;
  std::vector<TString> requiredDoubleElectronTriggers_;
  std::vector<TString> requiredEMTriggers_;

};

#endif
