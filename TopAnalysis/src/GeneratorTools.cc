#include "TopLJets2015/TopAnalysis/interface/GeneratorTools.h"
#include "TH1F.h"

std::vector< WeightSysts_t > getWeightSysts(TFile *f,TString sample) {
  std::vector< WeightSysts_t > systsOfInterest;

  TH1 *h=(TH1 *) f->Get("analysis/generator_initrwgt");
  if(h==0) return systsOfInterest;

  for(Int_t xbin=1; xbin<=h->GetNbinsX(); xbin++) {
    TString label=h->GetXaxis()->GetBinLabel(xbin);
    if(label.Length()==0) continue;

    if(sample=="TTJets2016") {
      //NNPDF30
      /*
      if(label.Contains("PDF set = 260")) {
        Int_t start=label.Index("PDF set = 260")+13;
        TString id=label(start,3);
        systsOfInterest.push_back( WeightSysts_t("PDF"+id, xbin-1) );
      }
      */
      if(label.Contains("PDF set = 260001")) systsOfInterest.push_back( WeightSysts_t("PDF001", xbin-1) );
      if(label.Contains("PDF set = 265000")) systsOfInterest.push_back( WeightSysts_t("PDF101", xbin-1) );
      if(label.Contains("PDF set = 266000")) systsOfInterest.push_back( WeightSysts_t("PDF102", xbin-1) );
    }
  }

  return systsOfInterest;
}

std::vector< WeightSysts_t > getPartonShowerWeightSysts(TFile *f) {
  std::vector< WeightSysts_t > systsOfInterest;

  TH1 *h=(TH1 *) f->Get("analysis/genlumiwgts");
  if(h==0) return systsOfInterest;

  for(Int_t xbin=1; xbin<=h->GetNbinsX(); xbin++) {
    TString label=h->GetXaxis()->GetBinLabel(xbin);
    if(label.Length()==0) continue;
    if(label=="nominal" || label=="Baseline") continue;
    systsOfInterest.push_back( WeightSysts_t(label,xbin-1) );
  }

  return systsOfInterest;
}
