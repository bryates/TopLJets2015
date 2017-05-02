#include <iostream>
#include "TopLJets2015/TopAnalysis/interface/Trigger.h"

std::vector<TString> muTriggers = {"HLT_IsoMu24_v", "HLT_IsoTkMu24_v","HLT_IsoMu22_v", "HLT_IsoTkMu22_v", "HLT_IsoMu22_eta2p1_v", "HLT_IsoTkMu22_eta2p1_v", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"};
std::vector<TString> eleTriggers = {"HLT_Ele27_WPTight_Gsf_v", "HLT_Ele32_WPTight_Gsf_v", "HLT_Ele32_eta2p1_WPTight_Gsf_v", "HLT_Ele25_eta2p1_WPTight_Gsf_v", "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"}; 

int count(TString st, TString sub) {
  int c = 0;
  int index = 0;
  while((index = st.Index(sub))!=-1) {
    st.Remove(index,sub.Length());
    c++;
  }
  return c;
}

Trigger::Trigger(int muTrig, int eleTrig, bool debug) {
  debug_ = debug;
  dataType_ = None;
  setMuonTrigger(muTrig);
  setElectronTrigger(eleTrig);
}

Trigger::~Trigger() {};

void Trigger::setMuonTrigger(int muTrig) {
  if(debug_) std::cout << "====== Parsing triggers ======" << std::endl;
  for(size_t imu = 0; imu < muTriggers.size(); imu++) {
    triggers_[muTriggers[imu]] = ((muTrig>>imu)&0x1);
    if(debug_) std::cout << ((muTrig>>imu)&0x1) << " " << muTriggers[imu] << std::endl;
  }
  if(debug_) std::cout << "=== Parsing triggers done! ===" << std::endl;
  /*
  addRequiredMuonTrigger("HLT_IsoMu24_v");
  addRequiredMuonTrigger("HLT_IsoTkMu24_v");
  */
}

void Trigger::setElectronTrigger(int eleTrig) {
  for(size_t iel = 0; iel < eleTriggers.size(); iel++) {
    triggers_[eleTriggers[iel]] = ((eleTrig>>iel)&0x1);
    if(debug_) std::cout << ((eleTrig>>iel)&0x1) << " " << eleTriggers[iel] << std::endl;
  }
}

void Trigger::setDataType(TString fileName) {
  isData_ = fileName.Contains("Data");
  if(!isData_) {
    dataType_ = MC;
    return;
  }
  if(fileName.Contains("SingleMuon")) {
    dataType_ = SingleMuon;
    return;
  }
  if(fileName.Contains("SingleElectron")) {
    dataType_ = SingleElectron;
    return;
  }
  if(fileName.Contains("DoubleMuon")) {
    dataType_ = DoubleMuon;
    return;
  }
  if(fileName.Contains("DoubleEG")) {
    dataType_ = DoubleElectron;
    return;
  }
  if(fileName.Contains("MuonEG")) {
    dataType_ = ElectronMuon;
    return;
  }
}

void Trigger::addRequiredMuonTrigger(TString trigger) {
  if(std::find(requiredMuonTriggers_.begin(), requiredMuonTriggers_.end(), trigger) != requiredMuonTriggers_.end()) return;
  if(std::find(muTriggers.begin(), muTriggers.end(), trigger) == muTriggers.end()) {
    if(debug_) std::cout << trigger << " does not exist!" << std::endl;
    return;
  }
  requiredMuonTriggers_.push_back(trigger);
}

void Trigger::addRequiredElectronTrigger(TString trigger) {
  if(std::find(requiredElectronTriggers_.begin(), requiredElectronTriggers_.end(), trigger) != requiredElectronTriggers_.end()) return;
  if(std::find(eleTriggers.begin(), eleTriggers.end(), trigger) == eleTriggers.end()) {
    if(debug_) std::cout << trigger << " does not exist!" << std::endl;
    return;
  }
  requiredElectronTriggers_.push_back(trigger);
}

void Trigger::addRequiredDoubleMuonTrigger(TString trigger) {
  if(std::find(requiredDoubleMuonTriggers_.begin(), requiredDoubleMuonTriggers_.end(), trigger) != requiredDoubleMuonTriggers_.end()) return;
  if(std::find(muTriggers.begin(), muTriggers.end(), trigger) == muTriggers.end()) {
    if(debug_) std::cout << trigger << " does not exist!" << std::endl;
    return;
  }
  requiredDoubleMuonTriggers_.push_back(trigger);
}

void Trigger::addRequiredDoubleElectronTrigger(TString trigger) {
  if(std::find(requiredDoubleElectronTriggers_.begin(), requiredDoubleElectronTriggers_.end(), trigger) != requiredDoubleElectronTriggers_.end()) return;
  if(std::find(eleTriggers.begin(), eleTriggers.end(), trigger) == eleTriggers.end()) {
    if(debug_) std::cout << trigger << " does not exist!" << std::endl;
    return;
  }
  requiredDoubleElectronTriggers_.push_back(trigger);
}

void Trigger::addRequiredEMTrigger(TString trigger) {
  if(std::find(requiredEMTriggers_.begin(), requiredEMTriggers_.end(), trigger) != requiredEMTriggers_.end()) return;
  if(std::find(eleTriggers.begin(), eleTriggers.end(), trigger) == eleTriggers.end()) {
    if(debug_) std::cout << trigger << " does not exist!" << std::endl;
    return;
  }
  requiredEMTriggers_.push_back(trigger);
}

void Trigger::addRequiredMuonTrigger(std::vector<TString> triggers) {
  for(auto itMu = triggers.begin(); itMu != triggers.end(); itMu++)
    addRequiredMuonTrigger(*itMu);
}

void Trigger::addRequiredDoubleMuonTrigger(std::vector<TString> triggers) {
  for(auto itMu = triggers.begin(); itMu != triggers.end(); itMu++)
    addRequiredDoubleMuonTrigger(*itMu);

}

void Trigger::deleteRequiredMuonTrigger(TString trigger) {
  requiredMuonTriggers_.erase(std::remove(requiredMuonTriggers_.begin(), requiredMuonTriggers_.end(), trigger), requiredMuonTriggers_.end());
}

void Trigger::addRequiredElectronTrigger(std::vector<TString> triggers) {
  for(auto itEl = triggers.begin(); itEl != triggers.end(); itEl++)
    addRequiredElectronTrigger(*itEl);
}

void Trigger::addRequiredDoubleElectronTrigger(std::vector<TString> triggers) {
  for(auto itEl = triggers.begin(); itEl != triggers.end(); itEl++)
    addRequiredDoubleElectronTrigger(*itEl);
}

void Trigger::addRequiredEMTrigger(std::vector<TString> triggers) {
  for(auto itEl = triggers.begin(); itEl != triggers.end(); itEl++)
    addRequiredEMTrigger(*itEl);
}

void Trigger::deleteRequiredElectronTrigger(TString trigger) {
  requiredElectronTriggers_.erase(std::remove(requiredElectronTriggers_.begin(), requiredElectronTriggers_.end(), trigger), requiredElectronTriggers_.end());
}

bool Trigger::triggerFired(TString triggerName) {
  if(debug_) {
    std::cout << triggers_[triggerName] << " " << triggerName << std::endl;
  }
  return triggers_[triggerName];
}

bool Trigger::muonFired() {
  if(debug_) std::cout << "Checking if required trigger(s) fired" << std::endl;
  for(size_t itrig = 0; itrig < requiredMuonTriggers_.size(); itrig++) {
    if(triggerFired(requiredMuonTriggers_[itrig])) return true;
  }
  if(debug_) std::cout << "No required triggers fired" << std::endl;
  return false;
}

bool Trigger::electronFired() {
  if(debug_) std::cout << "Checking if required trigger(s) fired" << std::endl;
  for(size_t itrig = 0; itrig < requiredElectronTriggers_.size(); itrig++)
    if(triggerFired(requiredElectronTriggers_[itrig])) return true;
  return false;
}

bool Trigger::doubleMuonFired() {
  if(debug_) std::cout << "Checking if required trigger(s) fired" << std::endl;
  for(size_t itrig = 0; itrig < requiredDoubleMuonTriggers_.size(); itrig++) {
    if(triggerFired(requiredDoubleMuonTriggers_[itrig])) return true;
  }
  if(debug_) std::cout << "No required triggers fired" << std::endl;
  return false;
}

bool Trigger::doubleElectronFired() {
  if(debug_) std::cout << "Checking if required trigger(s) fired" << std::endl;
  for(size_t itrig = 0; itrig < requiredDoubleElectronTriggers_.size(); itrig++) {
    if(triggerFired(requiredDoubleElectronTriggers_[itrig])) return true;
  }
  if(debug_) std::cout << "No required triggers fired" << std::endl;
  return false;
}

bool Trigger::EMFired() {
  if(debug_) std::cout << "Checking if required trigger(s) fired" << std::endl;
  for(size_t itrig = 0; itrig < requiredEMTriggers_.size(); itrig++) {
    if(triggerFired(requiredEMTriggers_[itrig])) return true;
  }
  if(debug_) std::cout << "No required triggers fired" << std::endl;
  return false;
}

bool Trigger::isMuonFile() {
  return dataType_ == SingleMuon;
}

bool Trigger::isElectronFile() {
  return dataType_ == SingleElectron;
}

bool Trigger::isDoubleMuonFile() {
  return dataType_ == DoubleMuon;
}

bool Trigger::isDoubleElectronFile() {
  return dataType_ == DoubleElectron;
}

bool Trigger::isEMFile() {
  return dataType_ ==  ElectronMuon;
}

bool Trigger::isMCFile() {
  return dataType_ ==  MC;
}

bool Trigger::isSingleMuonEvent() {
  if(dataType_ == None) {
    if(debug_) std::cout << "No type is set!" << std::endl;
    return false;
  }
  if(isData_ && isMuonFile()) {
    if(debug_) std::cout << "Event is" << (isMuonFile() ? " " : " not ")
                         << "from a single muon file" << std::endl;
  }
  else if(isData_) return false;
  return muonFired();
}

//Check for Single Muon trigger based on Leptons class
bool Trigger::isSingleMuonEvent(Leptons leps) {
  if(leps.size() != 1) return false;
  if(abs(leps[0].getPdgId())!=13) return false;
  if(isMCFile()) return true;
  return isSingleMuonEvent();
}

bool Trigger::isSingleElectronEvent() {
  if(dataType_ == None) {
    if(debug_) std::cout << "No type is set!" << std::endl;
    return false;
  }
  if(isData_ && isElectronFile()) {
    if(debug_) std::cout << "Event is" << (isElectronFile() ? " " : " not ")
                         << "from a single electron file" << std::endl;
  }
  else if(isData_) return false;
  return electronFired();
}

//Check for Single Electron trigger based on Leptons class
bool Trigger::isSingleElectronEvent(Leptons leps) {
  if(leps.size() != 1) return false;
  if(abs(leps[0].getPdgId())!=11) return false;
  if(isMCFile()) return true;
  return isSingleElectronEvent();
}

bool Trigger::isDoubleMuonEvent() {
  if(dataType_ == None) {
    if(debug_) std::cout << "No type is set!" << std::endl;
    return false;
  }
  if(isData_ && isDoubleMuonFile()) {
    if(debug_) std::cout << "Event is" << (isDoubleMuonFile() ? " " : " not ")
                         << "from a double muon file" << std::endl;
  }
  else if(isData_) return false;
  return doubleMuonFired();
}

//Check for Double Muon trigger based on Leptons class
bool Trigger::isDoubleMuonEvent(Leptons leps) {
  if(leps.size() < 2) return false;
  if(abs(leps[0].getPdgId()*leps[1].getPdgId())!=13*13) return false;
  if(isMCFile()) return true;
  return isDoubleMuonEvent();
}

bool Trigger::isDoubleElectronEvent() {
  if(dataType_ == None) {
    if(debug_) std::cout << "No type is set!" << std::endl;
    return false;
  }
  if(isData_ && isDoubleElectronFile()) {
    if(debug_) std::cout << "Event is" << (isDoubleElectronFile() ? " " : " not ")
                         << "from a double electron file" << std::endl;
  }
  else if(isData_) return false;
  return doubleElectronFired();
}

//Check for Double Electron trigger based on Leptons class
bool Trigger::isDoubleElectronEvent(Leptons leps) {
  if(leps.size() < 2) return false;
  if(abs(leps[0].getPdgId()*leps[1].getPdgId())!=11*11) return false;
  if(isMCFile()) return true;
  return isDoubleElectronEvent();
}

bool Trigger::isEMEvent() {
  if(dataType_ == None) {
    if(debug_) std::cout << "No type is set!" << std::endl;
    return false;
  }
  if(isData_ && isEMFile()) {
    if(debug_) std::cout << "Event is" << (isEMFile() ? " " : " not ")
                         << "from an EG file" << std::endl;
  }
  else if(isData_) return false;
  return EMFired();
}

//Check for EM trigger based on Leptons class
bool Trigger::isEMEvent(Leptons leps) {
  if(leps.size() < 2) return false;
  if(abs(leps[0].getPdgId()*leps[1].getPdgId())!=11*13) return false;
  if(isMCFile()) return true;
  return isEMEvent();
}

