#ifndef _Particle_h
#define _Particle_h
#include "TLorentzVector.h"

enum particleType { Veto, Loose, Medium, TightNoIso, Tight };

//Defind Particle class
//
class Particle {
 public:
  Particle(float pt, float eta, float phi, float mass, int pdgId, float relIso, int pid);
  Particle();
  ~Particle();
  //isMuon (if abs(pdgId)==13)
  bool isMuon();
  //isElectron (if abs(pdgId)==11)
  bool isElectron();
  TLorentzVector& getVec();
  float Pt();
  float Eta();
  float Phi();
  float Mass();
  int charge();
  particleType& getType();
  float& getRelIso();
  int& getPdgId();

 private:
  TLorentzVector p4_;
  int pdgId_;
  particleType t_; //Store type (tight/medium/loose)
  float relIso_;
  int pid_;

};

#endif
