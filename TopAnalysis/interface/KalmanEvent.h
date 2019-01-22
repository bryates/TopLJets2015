#ifndef _KalmanEvent_h
#define _KalmanEvent_h
#include <TTree.h>
#include <TGraphAsymmErrors.h>
#include <TRandom3.h>
#include <vector>
#include "TopLJets2015/TopAnalysis/interface/Jet.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

class KalmanEvent
{
 public:
  KalmanEvent(bool debug);
  ~KalmanEvent();
  void loadEvent(const MiniEvent_t &ev);
  void buildJets();
  inline int getNJPsi() { return njpsi_; }
  inline int getNMeson() { return nmeson_; }
  inline int getEvent() { return event_; }
  inline bool isJPsiEvent() { return njpsi_ > 0; }
  inline bool isMesonEvent() { return nmeson_ > 0; }
  inline bool isGoodEvent() { return nmeson_ > 0; }
  bool isGoodJet(int idx);
  std::vector<Jet> getJets();

 private:
  MiniEvent_t ev_;
  bool debug_;
  int event_;
  int row_;
  int njpsi_;
  int nmeson_;
  float vtxProb_;
  float chi2_;
  float l3dsig_;
  float opang_;
  float csv_;
  std::vector<Jet> jets_ {};

};

#endif
