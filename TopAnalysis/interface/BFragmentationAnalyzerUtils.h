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
  float xb;
  int leadTagId;
  bool hasSemiLepDecay,hasTauSemiLepDecay;
  int nbtags,nctags,ntautags;
};

JetFragInfo_t analyzeJet(const reco::GenJet &genJet,float tagScale=1.0E+20);

#endif
