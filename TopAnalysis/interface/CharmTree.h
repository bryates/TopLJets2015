#ifndef CharmTree_h
#define CharmTree_h

#include <string>
#include <iostream>
#include "TString.h"
#include "TTree.h"
#include "TopLJets2015/TopAnalysis/interface/CharmEvent.h"
#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"
#include "TopLJets2015/TopAnalysis/interface/Jet.h"

class CharmTree {

 public:
  CharmTree(TTree *t, TString, TString, bool debug=false);
  ~CharmTree();
  inline float getWgt() { return norm_ * sfs_.first * puWgt_ * top_pt_wgt_; };
  void Fill(Leptons, TString, TString name="");
  void FillGen(std::vector<Particle>, TString, TString name="");
  void Fill(std::vector<Jet> lightJetsVec, std::vector<Jet> bJetsVec, std::vector<Jet> allJetsVec, TString, TString name="");
  void Fill(Double_t weight, Int_t N, Double_t pt, Double_t eta, Double_t phi);
  void Fill(std::vector<pfTrack>&, TString, TString name=""); //e.g. J/Psi plots
  void Fill(std::vector<pfTrack>&, Leptons, TString, TString name=""); //e.g. J/Psi plots+mu
  void Fill(std::vector<pfTrack>&, Jet, TString, TString name=""); //e.g. J/Psi plots+mu
  void Fill(CharmEvent_t &ev_, double nvtx, double HT, double ST, double MET, std::vector<Jet> lightJets);
  void Fill(CharmEvent_t &ev_, std::vector<pfTrack>&, Leptons, Jet, TString, TString name="", int event=0, std::vector<pfTrack> genMatch={}, std::vector<float> frag={1.,1.,1.,1.});
  void SetNorm(float);
  void SetSFs(float, float unc=0.);
  void SetPuWgt(float);
  void SetTopPtWgt(float);
  void SetTrackerWgt(float tracker_wgt);
  void SetPiWgt(float pi_wgt, float unc);
  void SetLumi(float);
  void CheckRunPeriod(TString);
  void Write();

 private:
  std::map<TString, TH1 *> allPlots;
  TString runPeriod_;
  TString name_;
  bool debug_;
  bool isGood_;
  float norm_;
  std::pair <float,float> sfs_;
  float puWgt_;
  float top_pt_wgt_;
  float tracker_wgt_;
  float lumi_;
  std::pair <float,float> pi_wgt_;
  //CharmEvent_t ev_;

};


#endif
