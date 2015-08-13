#ifndef _readtree_h_
#define _readtree_h_
#include <iostream>
#include <TString.h>

void ReadTree(TString filename="chargediso_QCD_1000_MuEnriched_PU20bx25.root",TString output="plots",bool useChPt=true,int minVtx=-1, int maxVtx=9999999);
void RunOverSamples(TString output="plots/",bool useChPt=true,int minVtx=-1, int maxVtx=9999999);

#endif
