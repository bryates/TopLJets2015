#ifndef StdPlots_h
#define StdPlots_h

#include <string>
#include <iostream>
#include <map>
#include "TString.h"
#include "TH1F.h"
#include "TTree.h"
#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"
#include "TopLJets2015/TopAnalysis/interface/Jet.h"

enum CharmMeson { JPsi=0, D0, D0mu };

class StdPlots {

 public:
  TH1F *h_particles, *h_pt, *h_eta, *h_phi;
 
  StdPlots(TString, TString, bool debug=false);
  ~StdPlots();
  //inline friend StdPlots operator+(const StdPlots&, const StdPlots&);
  //StdPlots Combine(const StdPlots&);
  friend StdPlots Combine(const StdPlots&, const StdPlots&);
  inline friend StdPlots operator+(const StdPlots &lhs, const StdPlots &rhs) { return Combine(lhs,rhs); };
  inline float getWgt() { return norm_ * sfs_.first * puWgt_ * top_pt_wgt_ * tracker_wgt_ * pi_wgt_.first; };
  //inline float getUnc() { return sqrt( pow(sfs_.second,2) + pow(pi_wgt_.second,2) ); }
  inline std::pair<float,float> getUnc() { return std::pair<float, float>(sqrt( pow(sfs_.second.first,2) + pow(pi_wgt_.second.first,2) ),
                                                            sqrt( pow(sfs_.second.second,2) + pow(pi_wgt_.second.second,2) )); }
  void Fill(Leptons, TString, TString name="");
  void FillGen(std::vector<Particle>, TString, TString name="");
  //void Fill(std::vector<Jet> lightJetsVec, std::vector<Jet> kJetsVec, std::vector<Jet> allJetsVec, TString, TString name="");
  void Fill(std::vector<Jet>&, std::vector<Jet>&, std::vector<Jet>&, TString, TString name="");
  void Fill(Leptons&, std::vector<Jet>&, std::vector<Jet>&, std::vector<Jet>&, TString, TString name="");
  void Fill(double nevt, double nvtx, double HT, double ST, double MET, TString chTag, TString name="");
  void Fill(Double_t weight, Int_t N, Double_t pt, Double_t eta, Double_t phi);
  void Fill(std::vector<pfTrack>&, TString, TString name=""); //e.g. J/Psi plots
  void Fill(std::vector<pfTrack>&, Leptons, TString, TString name=""); //e.g. J/Psi plots+mu
  void Fill(std::vector<pfTrack>&, Jet, TString, TString name=""); //e.g. J/Psi plots+mu
  void Fill(std::vector<pfTrack>&, Leptons, Jet, TString, TString name="");
  void SetNorm(float);
  void SetSFs(float);
  void SetSFs(float,float);
  void SetSFs(float,float,float);
  void SetPuWgt(float);
  void SetTopPtWgt(float);
  void SetTrackerWgt(float);
  void SetPiWgt(float,float);
  void SetPiWgt(float,float,float);
  void SetRbWgt(float,CharmMeson);
  void CheckRunPeriod(TString);
  void Write();

 private:
  std::map<TString, TH1 *> allPlots;
  TString runPeriod_;
  TString name_;
  bool debug_;
  bool isGood_;
  float norm_;
  //std::pair <float,float> sfs_;
  std::pair <float,std::pair<float,float>> sfs_;
  float puWgt_;
  float top_pt_wgt_;
  float tracker_wgt_;
  //std::pair <float,float> pi_wgt_;
  std::pair <float,std::pair<float,float>> pi_wgt_;
  //std::vector<float> top_pt_wgt_vec;
  float rbWgt_[3];

};


#endif
