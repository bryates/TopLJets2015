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
  if(p.getType() != reqType_) return;
  if(debug_) std::cout << "Type is good!!" << std::endl;

  //add good particles to collection
  leps_.push_back(p);
}

int Leptons::getSize() {
  int s = leps_.size();
  if(debug_) std::cout << "There " << (s > 1 ? "are " : "is ") << s << " lepton" << (s > 1 ? "s " : " ") << "in this collection" << std::endl;
  return s;
}

Particle& Leptons::getElement(int pos) {
  return leps_[pos];
}
