#ifndef _readtree_h_
#define _readtree_h_
#include <iostream>
#include <TString.h>

void ReadTree(TString filename="chargediso_QCD_1000_MuEnriched_PU20bx25.root",TString output="plots",int chToSelect=13);
void RunOverSamples(TString output="plots/",int chToSelect=13);

#endif
