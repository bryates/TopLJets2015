#ifndef GENERATORTOOLS_h
#define GENERATORTOOLS_h

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

#include <vector>
#include <string>
#include <map>

#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TF1.h"

//available theory systs
typedef std::pair<TString, float> WeightSysts_t;
std::vector< WeightSysts_t > getWeightSysts(TFile *, TString sample="TTJets2016");
std::vector< WeightSysts_t > getPartonShowerWeightSysts(TFile *f);


#endif
