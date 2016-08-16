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

class Jet {

 public:
  Jet(TLorentzVector vec_, float csv_);
  ~Jet();
  TLorentzVector getVec();
  float getCSV();
  int getJetIndex();
  void setVec(TLorentzVector vec_);
  void setCSV(float csv_);
  int setJetIndex(int jetindex_);
 
 private:
  TLorentzVector vec;
  float csv;
  int jetindex;

};

Jet::Jet(TLorentzVector vec_, float csv_) {
  vec = vec_;
  csv = csv_;
}

Jet::~Jet() {};

#endif
