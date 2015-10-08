#ifndef _readtree_h_
#define _readtree_h_

#include <TString.h>

enum FlavourSplitting {NOFLAVOURSPLITTING=0, UDSGSPLITTING=1, CSPLITTING=4, BSPLITTING=5 };
enum GenWeightMode { NOGENWGT=0, ONLYSIGN=1, FULLWEIGHT=2 };

void ReadTree(TString filename,
	      TString outname,
	      Int_t channelSelection=13,
	      Int_t chargeSelection=0,
	      Float_t norm=1.0,
	      Bool_t isTTbar=false,
	      FlavourSplitting flavourSplitting=NOFLAVOURSPLITTING,
	      GenWeightMode genWgtMode=NOGENWGT);


#endif
