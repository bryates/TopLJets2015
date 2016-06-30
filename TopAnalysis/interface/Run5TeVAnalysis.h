#ifndef _run5tevanalysis_h_
#define _run5tevanalysis_h_

#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"

void Run5TeVAnalysis(TString filename,
		     TString outname,
		     Int_t channelSelection, 
		     Int_t chargeSelection, 
		     FlavourSplitting flavourSplitting,
		     TH1F *normH, 
		     Bool_t runSysts,
		     TString era);

#endif
