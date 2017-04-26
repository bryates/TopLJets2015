#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TLorentzVector.h"

Particle::Particle(float pt, float eta, float phi, float mass, int pdgId, float relIso, int pid) {
  p4_.SetPtEtaPhiM(pt, eta, phi, mass);
  pdgId_ = pdgId;
  relIso_ = relIso;
  pid_ = pid;
  //tight muon ID or tight ID except Iso electron
  if(isMuon() && (pid_&0x1)) t_ = Tight;
  else if(isElectron() && (pid_>>1)&0x1) t_ = Tight;
  else if(isElectron() && (pid_>>2)&0x1) t_ = TightNoIso;
  //Veto if muon isLoose or electron VetoId
  else if(isMuon() && (pid_>>2)&0x1) t_ = Veto;
  else if(isElectron() && (pid_)&0x1) t_ = Veto;
}

Particle::Particle() {} ;

Particle::~Particle() {} ;

bool Particle::isMuon() {
  return (abs(pdgId_) == 13);
}
bool Particle::isElectron() {
  return (abs(pdgId_) == 11);
}
TLorentzVector& Particle::getVec() {
  return p4_;
}
float Particle::Pt() {
  return p4_.Pt();
}
float Particle::Eta() {
  return p4_.Eta();
}
float Particle::Phi() {
  return p4_.Phi();
}
float Particle::Mass() {
  return p4_.M();
}
int Particle::charge() {
  return abs(pdgId_)/pdgId_;
}
particleType& Particle::getType() {
  return t_;
}
float& Particle::getRelIso() {
  return relIso_;
}
int& Particle::getPdgId() {
  return pdgId_;
}

