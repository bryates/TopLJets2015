#ifndef _top16006_h_
#define _top16006_h_

#include <TString.h>
#include <TGraph.h>
#include <TH1F.h>
#include "TLorentzVector.h"
#include <map>
#include <vector>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"

void RunTop16006(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts,
		 TString era);
#endif
