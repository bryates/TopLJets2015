#ifndef _other_functions_h_
#define _other_functions_h_

#include <TLorentzVector.h>
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include <vector>
#include "TMath.h"

Int_t npf, pf_j[5000];
Double_t pi = TMath::Pi();

int firstTrackIndex(int jetindex) {
    // Find index of the first track in this jet
    int result = 0;
    for (result = 0; result < npf; ++result) {
        if (pf_j[result]==jetindex) {
            break; // at this point, result is correct
        }
    }
    return result;
}

int firstTrackIndex(int jetindex, std::vector<std::tuple<int,int,float>> *jets) {
    // Find index of the first track in this jet
    int result = 0;
    for (result = 0; result < (int)jets->size(); ++result) {
        if (std::get<1>(jets->at(result))==jetindex) {
            break; // at this point, result is correct
        }
    }
    return result;
}

int firstTrackIndex(int jetindex, std::vector<std::tuple<int,float,float>> *jets) {
    // Find index of the first track in this jet
    int result = 0;
    for (result = 0; result < (int)jets->size(); ++result) {
        if (std::get<0>(jets->at(result))==jetindex) {
            break; // at this point, result is correct
        }
    }
    return result;
}

double DR(TLorentzVector &v1, TLorentzVector &v2) {
    double deta = v1.Eta() - v2.Eta();
    double dphi = v1.Phi() - v2.Phi();
    dphi = dphi < pi ? dphi : pi - dphi;
    return TMath::Sqrt( deta*deta+dphi*dphi );
}

bool VecSort(TLorentzVector j1, TLorentzVector j2) { return j1.Pt() > j2.Pt(); }
//bool sortJetTuple(std::tuple<int,float> i, std::tuple<int,float> j) { return std::get<0>(i) > std::get<0>(j); }
bool sortJetTuple(std::tuple<int,int,float> i, std::tuple<int,int,float> j) { return std::get<1>(i) < std::get<1>(j) || ( std::get<1>(i) == std::get<1>(j) && std::get<2>(i) > std::get<2>(j) ); }
bool sortJetCSVTuple(std::tuple<int,float,float> i, std::tuple<int,float,float> j) { return std::get<1>(i) < std::get<1>(j) || ( std::get<1>(i) == std::get<1>(j) && std::get<2>(i) > std::get<2>(j) ); }

class pfTrack {

 public:
  pfTrack(TLorentzVector p4, float dxy, float dxyE, float dz, float dzE, int pfid);
  ~pfTrack();
  float Pt();
  float Eta();
  float Phi();
  int getPfid();
  float getDxy();
  float getDxyE();
  float getDz();
  float getDzE();
  TLorentzVector &getVec();
  void setPfid(int pfid);

 private:
  TLorentzVector vec_;
  float dxy_;
  float dxyE_;
  float dz_;
  float dzE_;
  int pfid_;

};

pfTrack::pfTrack(TLorentzVector p4, float dxy, float dxyE, float dz, float dzE, int pfid) : vec_(p4), dxy_(dxy), dxyE_(dxyE), dz_(dz), dzE_(dzE), pfid_(pfid)
{ }
pfTrack::~pfTrack() {};
int pfTrack::getPfid() { return pfid_ ; }
float pfTrack::Pt() { return vec_.Pt() ; }
float pfTrack::Eta() { return vec_.Eta() ; }
float pfTrack::Phi() { return vec_.Phi() ; }
float pfTrack::getDxy() { return dxy_ ; }
float pfTrack::getDxyE() { return dxyE_ ; }
float pfTrack::getDz() { return dz_ ; }
float pfTrack::getDzE() { return dzE_ ; }
TLorentzVector &pfTrack::getVec() { return vec_ ; }
void pfTrack::setPfid(int pfid) { pfid_ = pfid ; }

//typedef std::pair<TLorentzVector,int> IdTrack;
typedef std::pair<pfTrack,int> IdTrack;

class Jet {

 public:
  Jet(TLorentzVector p4, float csv, int idx);
  ~Jet();
  //void addTrack(TLorentzVector p4, int pfid);
  void addTrack(pfTrack pf, int pfid);
  void addTrack(int idx);
  void addDxy(float dxy, float dxyE);
  void addDz(float dz, float dzE);
  void addDz(int idx);
  TLorentzVector &getVec();
  float &getCSV();
  int &getJetIndex();
  int &getIndex(int idx);
  float &getDxy(int idx);
  float &getDz(int idx);
  float &getDxyE(int idx);
  float &getDzE(int idx);
  std::vector<IdTrack> &getTracks();
  IdTrack getTrack(int idx,float mass);
  void sortTracksByPt();
 
 private:
  TLorentzVector p4_;
  float csv_;
  int jetindex_;
  int idx_;
  std::vector<IdTrack> trks_;
  std::vector<int> index_;
  std::vector<float> dxy_;
  std::vector<float> dxyE_;
  std::vector<float> dz_;
  std::vector<float> dzE_;

};

Jet::Jet(TLorentzVector p4, float csv, int idx) : p4_(p4), csv_(csv), idx_(idx) { }

/*
Jet::Jet(TLorentzVector p4, float csv, int idx) {
  p4_ = p4;
  csv_ = csv;
  idx_ = idx;
}
*/

Jet::~Jet() {};

//void Jet::addTrack(TLorentzVector p4, int pfid) { trks_.push_back( IdTrack(p4,pfid) ); }
void Jet::addTrack(pfTrack pf, int pfid) { trks_.push_back( IdTrack(pf,pfid) ); }
void Jet::addTrack(int idx) { index_.push_back( idx ) ; }
void Jet::addDxy(float dxy, float dxyE) { dxy_.push_back( dxy ) ; dxyE_.push_back( dxyE ) ; }
void Jet::addDz(float dz, float dzE) { dz_.push_back( dz ) ; dzE_.push_back( dzE) ; }
TLorentzVector &Jet::getVec() { return p4_; }
float &Jet::getCSV() { return csv_; }
int &Jet::getJetIndex() { return idx_; }
int &Jet::getIndex(int idx) { return index_[idx]; }
float &Jet::getDxy(int idx) { return dxy_[idx]; }
float &Jet::getDz(int idx) { return dz_[idx]; }
float &Jet::getDxyE(int idx) { return dxyE_[idx]; }
float &Jet::getDzE(int idx) { return dzE_[idx]; }
std::vector<IdTrack> &Jet::getTracks() { return trks_; }

bool sortJetsByPt(Jet i, Jet j)  { return i.getVec().Pt() > j.getVec().Pt(); }
bool sortJetsByCSV(Jet i, Jet j) { return i.getCSV() > j.getCSV(); }
bool sortIdTracksByPt(IdTrack i, IdTrack j)  { return i.first.Pt() > j.first.Pt(); }

void Jet::sortTracksByPt() { sort(trks_.begin(),trks_.end(), sortIdTracksByPt); }

#endif
