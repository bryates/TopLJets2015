#ifndef StdPlots_h
#define StdPlots_h

#include <string>
#include <iostream>
#include <map>
#include "TString.h"
#include "TH1F.h"
#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"
#include "TopLJets2015/TopAnalysis/interface/Jet.h"

class StdPlots {

 public:
  TH1F *h_particles, *h_pt, *h_eta, *h_phi;
 
  StdPlots(TString, TString, bool debug=false);
  ~StdPlots();
  void Fill(Leptons, TString, TString name="");
  void Fill(std::vector<Jet> lightJetsVec, std::vector<Jet> bJetsVec, std::vector<Jet> allJetsVec, TString, TString name="");
  void Fill(double nevt, double nvtx, double HT, double ST, double MET, TString chTag, TString name="");
  void Fill(Double_t weight, Int_t N, Double_t pt, Double_t eta, Double_t phi);
  void Fill(std::vector<pfTrack>, TString, TString name=""); //e.g. J/Psi plots
  void Fill(std::vector<pfTrack>, Leptons, TString, TString name=""); //e.g. J/Psi plots+mu
  void Fill(std::vector<pfTrack>, Jet, TString, TString name=""); //e.g. J/Psi plots+mu
  void Fill(std::vector<pfTrack>, Leptons, Jet, TString, TString name="");
  void SetNorm(float);
  void SetSFs(float);
  void SetPuWgt(float);
  void Write();

 private:
  std::map<TString, TH1 *> allPlots;
  TString runPeriod_;
  bool debug_;
  bool isGood_;
  float norm_;
  float sfs_;
  float puWgt_;


};


#endif
