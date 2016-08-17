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

typedef std::pair<TLorentzVector,int> IdTrack;

class Jet {

 public:
  Jet(TLorentzVector p4, float csv, int idx);
  ~Jet();
  void addTrack(TLorentzVector p4, int pfid);
  TLorentzVector &getVec();
  float &getCSV();
  int &getJetIndex();
  std::vector<IdTrack> &getTracks();
  IdTrack getTrack(int idx,float mass);
  void sortTracksByPt();
 
 private:
  TLorentzVector p4_;
  float csv_;
  int jetindex_;
  int idx_;
  std::vector<IdTrack> trks_;

};

Jet::Jet(TLorentzVector p4, float csv, int idx) {
  p4 = p4_;
  csv = csv_;
  idx = idx_;
}

Jet::~Jet() {};

void Jet::addTrack(TLorentzVector p4, int pfid) { trks_.push_back( IdTrack(p4,pfid) ); }
TLorentzVector &Jet::getVec() { return p4_; }
float &Jet::getCSV() { return csv_; }
std::vector<IdTrack> &Jet::getTracks() { return trks_; }

bool sortJetsByPt(Jet i, Jet j)  { return i.getVec().Pt() > j.getVec().Pt(); }
bool sortJetsByCSV(Jet i, Jet j) { return i.getCSV() > j.getCSV(); }
bool sortIdTracksByPt(IdTrack i, IdTrack j)  { return i.first.Pt() > j.first.Pt(); }

void Jet::sortTracksByPt() { sort(trks_.begin(),trks_.end(), sortIdTracksByPt); }

#endif
