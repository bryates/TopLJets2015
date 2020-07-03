#ifndef _BFragmentationAnalyzerUtils_h_
#define _BFragmentationAnalyzerUtils_h_

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#define IS_BHADRON_PDGID(id) ( ((abs(id)/100)%10 == 5) || (abs(id) >= 5000 && abs(id) <= 5999) )
#define IS_CHADRON_PDGID(id) ( ((abs(id)/100)%10 == 4) || (abs(id) >= 4000 && abs(id) <= 4999) )
#define IS_NEUTRINO_PDGID(id) ( (abs(id) == 12) || (abs(id) == 14) || (abs(id) == 16) )

struct JetFragInfo_t
{
  float xb,xb_charged,pt,charmId,motherId,eta,phi;
  int leadTagId;
  bool hasSemiLepDecay,hasTauSemiLepDecay,hasCharm;
  int nbtags,nctags,ntautags;
  std::vector< std::vector<double> > meson;
  std::vector<int> mesonId;
};

JetFragInfo_t analyzeJet(const reco::GenJet &genJet,float tagScale=1.0E+20);

#endif
