#ifndef _Leptons_h
#define _Leptons_h
#include <vector>
#include <iostream>
#include "TopLJets2015/TopAnalysis/interface/Particle.h"

class Leptons {
 public:
  Leptons(enum particleType reqType,bool debug_=false);
  Leptons(enum particleType reqType,enum particleType maxType,bool debug_=false);
  ~Leptons();
  void setMinPt(double);
  void setMaxEta(double);
  void setMaxRelIso(double);
  void setMaxType(enum particleType);
  void changeMinPt(double);
  void changeMaxEta(double);
  void changeMaxRelIso(double);
  void changeParticleType(enum particleType);
  void addParticle(Particle);
  void combineLeptons(Leptons);
  void sortLeptonsByPt();
  int getSize();
  inline size_t size() { return getSize(); }
  Particle& getElement(int);
  Particle& operator[](std::size_t);

 private:
  bool debug_;
  double minPt_;
  double maxEta_;
  double maxRelIso_;
  particleType reqType_;
  particleType maxType_;
  std::vector<Particle> leps_;

};

#endif
