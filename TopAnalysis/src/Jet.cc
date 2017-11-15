#include <TLorentzVector.h>
#include <vector>
#include "TopLJets2015/TopAnalysis/interface/Jet.h"

pfTrack::pfTrack(TLorentzVector p4, float dxy, float dxyE, float dz, float dzE, int pfid, int quality, bool highPurity) : vec_(p4), dxy_(dxy), dxyE_(dxyE), dz_(dz), dzE_(dzE), pfid_(pfid), quality_(quality), highPurity_(highPurity)
{ }
pfTrack::pfTrack(TLorentzVector p4, int pfid) : vec_(p4), dxy_(0), dxyE_(0), dz_(0), dzE_(0), pfid_(pfid)
{ }
pfTrack::pfTrack(TLorentzVector p4,float k_mass, float l3d, float sigmal3d, float chi2, float vtxProb, int pfid, int motherId) : vec_(p4),k_mass_(k_mass), l3d_(l3d), sigmal3d_(sigmal3d), chi2_(chi2), vtxProb_(vtxProb), pfid_(pfid), motherId_(motherId)
{ }
pfTrack::~pfTrack() {};
int pfTrack::getPfid() { return pfid_ ; }
int pfTrack::getMotherId() { return motherId_; }
int pfTrack::charge() { return pfid_ / abs(pfid_); }
int pfTrack::getQuality() { return quality_; }
int pfTrack::getGenT() { return genT_; }
bool pfTrack::highPurity() { return highPurity_; }
bool pfTrack::globalMuon() { return globalMuon_; }
bool pfTrack::trackerMuon() { return trackerMuon_; }
float pfTrack::Pt() { return vec_.Pt() ; }
float pfTrack::Eta() { return vec_.Eta() ; }
float pfTrack::Phi() { return vec_.Phi() ; }
float pfTrack::P() { return vec_.P() ; }
float pfTrack::Pz() { return vec_.Pz() ; }
float pfTrack::M() { return vec_.M() ; }
float pfTrack::DeltaR(pfTrack &rhs) { return vec_.DeltaR(rhs.getVec()) ; }
float pfTrack::getDxy() { return dxy_ ; }
float pfTrack::getDxyE() { return dxyE_ ; }
float pfTrack::getDz() { return dz_ ; }
float pfTrack::getDzE() { return dzE_ ; }
float pfTrack::getL3D() { return l3d_; }
float pfTrack::getSigmaL3D() { return sigmal3d_; }
TLorentzVector &pfTrack::getVec() { return vec_ ; }
void pfTrack::setPfid(int pfid) { pfid_ = pfid ; }
void pfTrack::setMass(float mass) { vec_.SetPtEtaPhiM(Pt(), Eta(), Phi(), mass); }
void pfTrack::setGlobalMuon(bool globalMuon) { globalMuon_ = globalMuon; }
void pfTrack::setTrackerMuon(bool trackerMuon) { trackerMuon_ = trackerMuon; }
void pfTrack::setGenT(int genT) { genT_ = genT; }
void pfTrack::print() { std:: cout << "pdgId=" << getPdgId() << " pT=" << Pt() << " eta=" << Eta() << " phi=" << Phi() << " mass=" << M() << std::endl; }

Jet::Jet(TLorentzVector p4, float csv, int idx) : p4_(p4), csv_(csv), idx_(idx) { }
Jet::Jet(TLorentzVector p4, float csv, int idx, float chargedPt, float PFPt) : p4_(p4), csv_(csv), idx_(idx), chargedPt_(chargedPt), PFPt_(PFPt) { }
Jet::Jet(TLorentzVector p4, float csv, int idx, float chargedPt, float PFPt, int genJet) : p4_(p4), csv_(csv), idx_(idx), chargedPt_(chargedPt), PFPt_(PFPt), genJet_(genJet) { }
Jet::Jet(TLorentzVector p4, float csv, int idx, float chargedPt,  float chargedPz, float chargedP, float PFPt, float PFPz, float PFP, int genJet) : p4_(p4), csv_(csv), idx_(idx), chargedPt_(chargedPt), chargedPz_(chargedPz), chargedP_(chargedP), PFPt_(PFPt), PFPz_(PFPz), PFP_(PFP), genJet_(genJet) { }

/*
Jet::Jet(TLorentzVector p4, float csv, int idx) {
  p4_ = p4;
  csv_ = csv;
  idx_ = idx;
}
*/

Jet::~Jet() {};

//void Jet::addTrack(TLorentzVector p4, int pfid) { trks_.push_back( IdTrack(p4,pfid) ); }
//void Jet::addTrack(pfTrack pf, int pfid) { trks_.push_back( IdTrack(pf,pfid) ); }
void Jet::addTrack(pfTrack pf) { trks_.push_back( pf ); }
void Jet::addTrack(int idx) { index_.push_back( idx ) ; }
void Jet::addDxy(float dxy, float dxyE) { dxy_.push_back( dxy ) ; dxyE_.push_back( dxyE ) ; }
void Jet::addDz(float dz, float dzE) { dz_.push_back( dz ) ; dzE_.push_back( dzE) ; }
TLorentzVector &Jet::getVec() { return p4_; }
float &Jet::getCSV() { return csv_; }
float Jet::getPt() { return p4_.Pt(); }
float &Jet::getChargedPt() { return chargedPt_; }
float &Jet::getPFPt() { return PFPt_; }
float Jet::getP() { return p4_.P(); }
float &Jet::getChargedP() { return chargedP_; }
float &Jet::getPFP() { return PFP_; }
float Jet::getPz() { return p4_.Pz(); }
float &Jet::getChargedPz() { return chargedPz_; }
float &Jet::getPFPz() { return PFPz_; }
int &Jet::getGenJet() { return genJet_; }
int &Jet::getJetIndex() { return idx_; }
int &Jet::getIndex(int idx) { return index_[idx]; }
float &Jet::getDxy(int idx) { return dxy_[idx]; }
float &Jet::getDz(int idx) { return dz_[idx]; }
float &Jet::getDxyE(int idx) { return dxyE_[idx]; }
float &Jet::getDzE(int idx) { return dzE_[idx]; }
std::vector<pfTrack> &Jet::getTracks() { return trks_; }

//void Jet::sortTracksByPt() { sort(trks_.begin(),trks_.end(), sortIdTracksByPt); }
void Jet::sortTracksByPt() { sort(trks_.begin(),trks_.end(), [](pfTrack i, pfTrack j) { return i.Pt() > j.Pt() ; } ); }

