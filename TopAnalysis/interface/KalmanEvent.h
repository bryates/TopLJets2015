#ifndef _KalmanEvent_h
#define _KalmanEvent_h
#include <TTree.h>
#include <vector>
#include "TopLJets2015/TopAnalysis/interface/Jet.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

class KalmanEvent
{
 public:
  KalmanEvent(TTree *tree);
  KalmanEvent(bool debug);
  ~KalmanEvent();
  void loadTree(TTree *tree);
  void loadEvent(int event);
  void loadEvent(const MiniEvent_t &ev);
  void buildJets();
  inline int getNJPsi() { return njpsi_; }
  inline int getNMeson() { return nmeson_; }
  inline int getEvent() { return event_; }
  inline bool isJPsiEvent() { return njpsi_ > 0; }
  inline bool isGoodEvent() { return nmeson_ > 0; }
  std::vector<Jet> getJets();

 private:
  TTree *tree_;
  MiniEvent_t ev_;
  bool debug_;
  int event_;
  int row_;
  int njpsi_;
  int nmeson_;
  float vtxProb_;
  float chi2_;
  std::vector<Jet> jets_;

};

#endif
