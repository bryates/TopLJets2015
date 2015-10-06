#ifndef _readtree_h_
#define _readtree_h_

#include <TString.h>

void ReadTree(TString filename,
	      TString outDir,
	      Int_t channelSelection=13,
	      Int_t chargeSelection=0,
	      Float_t norm=1.0,
	      Bool_t isTTbar=false);


#endif
