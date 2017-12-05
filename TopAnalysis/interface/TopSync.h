#ifndef _topsync_h_
#define _topsync_h_

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"

#include "TH1.h"
#include "TString.h"
#include "TFile.h"

void RunTopSync(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts,
		 TString era,
                 TString runPeriod,
                 Bool_t debug);

#endif
