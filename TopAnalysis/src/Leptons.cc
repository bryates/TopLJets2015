#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"

Leptons::Leptons(enum particleType reqType, bool debug) {
  reqType_ = reqType;
  debug_  = debug;
} 

//Leptons::Leptons(bool veto, bool debug) : veto_(veto), debug_(debug) {} 
Leptons::~Leptons() {};

void Leptons::setMinPt(double minPt) { minPt_ = minPt; }
void Leptons::setMaxEta(double maxEta) { maxEta_ = maxEta; }
void Leptons::setMaxRelIso(double maxRelIso) { maxRelIso_ = maxRelIso; }

void Leptons::changeMinPt(double minPt) { 
  minPt_ = minPt; 
  for(auto itP = leps_.begin(); itP != leps_.end(); itP++) {
    if(itP->Pt() < minPt_ && debug_) std::cout << "Removing particle below new pT. (" << itP->Pt() << " < " << minPt_ << ")" << std::endl;
    if(itP->Pt() < minPt_)
      leps_.erase(itP--);
  }
}

void Leptons::changeMaxEta(double maxEta) { 
  maxEta_ = maxEta; 
  for(auto itP = leps_.begin(); itP != leps_.end(); itP++) {
    if(itP->Eta() > maxEta_ && debug_) std::cout << "Removing particle above new eta. (" << itP->Eta() << " > " << maxEta_ << ")" << std::endl;
    if(itP->Eta() > maxEta_)
      leps_.erase(itP--);
  }
}
void Leptons::changeMaxRelIso(double maxRelIso) { 
  maxRelIso_ = maxRelIso; 
  for(auto itP = leps_.begin(); itP != leps_.end(); itP++) {
    if(itP->getRelIso() > maxRelIso_ && debug_) std::cout << "Removing particle above new relIso. (" << itP->getRelIso() << " > " << maxRelIso_ << ")" << std::endl;
    if(itP->getRelIso() > maxRelIso_)
      leps_.erase(itP--);
  }
}

void Leptons::changeParticleType(enum particleType reqType) { 
  reqType_ = reqType;
  for(auto itP = leps_.begin(); itP != leps_.end(); itP++) {
    if(itP->getType() < reqType_ && debug_) std::cout << "Removing particle with old type. (" << itP->getType() << ")" << std::endl;
    if(itP->getType() < reqType_)
      leps_.erase(itP--);
  }
}

void Leptons::addParticle(Particle p) {
  //check if particle passes required cuts
  if(debug_) std::cout << "Checking if particle pT greter than required. (" << p.Pt() << " > " << minPt_ << ")" << std::endl;
  if(p.Pt() < minPt_) return;
  if(debug_) std::cout << "pT passed!!" << std::endl;
  if(debug_) std::cout << "Checking if particle eta less than required. (" << p.Eta() << " < " << maxEta_ << ")" << std::endl;
  if(p.Eta() > maxEta_) return;
  if(debug_) std::cout << "eta passed!!" << std::endl;
  if(debug_) std::cout << "Checking if particle relIso less than required. (" << p.getRelIso() << " < " << maxRelIso_ << ")" << std::endl;
  if(p.getRelIso() > maxRelIso_) return;
  if(debug_) std::cout << "relIso passed!!" << std::endl;
  if(debug_) std::cout << "Checking if particle type is correct." << std::endl;
  if(p.getType() < reqType_) return;
  if(debug_) std::cout << "Type is good!!" << std::endl;

  //add good particles to collection
  leps_.push_back(p);
}

void Leptons::combineLeptons(Leptons lep) {
  leps_.insert(leps_.end(), lep.leps_.begin(), lep.leps_.end());
}

int Leptons::getSize() {
  int s = leps_.size();
  if(debug_) std::cout << "There " << (s != 1 ? "are " : "is ") << s << " lepton" << (s != 1 ? "s " : " ") << "in this collection" << std::endl;
  return s;
}

void Leptons::sortLeptonsByPt() {
  sort(leps_.begin(),leps_.end(), [](Particle i, Particle j) { return i.Pt() > j.Pt() ; } );
  if(leps_.size() > 2) leps_.resize(2); //keep only 2 hardest leptons if more
}

Particle& Leptons::getElement(int pos) {
  return leps_[pos];
}

Particle& Leptons::operator[](std::size_t idx) {
  return getElement(idx);
}

