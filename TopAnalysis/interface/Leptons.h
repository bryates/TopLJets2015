#ifndef _Leptons_h
#define _Leptons_h
#include <vector>
#include <iostream>
#include "TopLJets2015/TopAnalysis/interface/Particle.h"

class Leptons {
 public:
  Leptons(enum particleType reqType,bool debug_=false);
  ~Leptons();
  void setMinPt(double);
  void setMaxEta(double);
  void setMaxRelIso(double);
  void addParticle(Particle);
  int getSize();
  Particle& getElement(int);

 private:
  bool debug_;
  double minPt_;
  double maxEta_;
  double maxRelIso_;
  particleType reqType_;
  std::vector<Particle> leps_;

};

#endif
