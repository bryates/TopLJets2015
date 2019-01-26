#ifndef _fragevent_h_
#define _fragevent_h_

#include "TTree.h"

struct FragEvent_t
{
  FragEvent_t()
  {
    nB=0;
  }

  Bool_t isData;
  Int_t run,event,lumi;
  Int_t nB;

  //Fragmentation
  Float_t xb[500],xbc[500],pt[500];
  Int_t id[500];
};

void createFragEventTree(TTree *t,FragEvent_t &ev);
void attachToFragEventTree(TTree *t, FragEvent_t &ev,bool full=false);

#endif
