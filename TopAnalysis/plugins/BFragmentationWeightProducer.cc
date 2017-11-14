#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TopLJets2015/TopAnalysis/interface/BFragmentationAnalyzerUtils.h"

#include <memory>
#include <string>
#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

//#define IS_BHADRON_PDGID(id) ( ((abs(id)/100)%10 == 5) || (abs(id) >= 5000 && abs(id) <= 5999) )
//#define IS_NEUTRINO_PDGID(id) ( (abs(id) == 12) || (abs(id) == 14) || (abs(id) == 16) )
#define IS_JPSI_PDGID(id) ( (abs(id) == 443) )
#define IS_BQUARK_PDGID(id) ( abs(id) == 5 )



class FragmentationWeightProducer : public edm::EDAnalyzer {

 public:
  FragmentationWeightProducer(const edm::ParameterSet &);
  virtual void analyze(const edm::Event &, const edm::EventSetup &);
  virtual void endJob() override;

 private:
  std::string fragModel;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;
  std::map<std::string, TGraph *> wgtGr_;
  TTree *data_;
  Int_t nB_;
  Float_t xb_[100], model_[100];

};
 
FragmentationWeightProducer::FragmentationWeightProducer(const edm::ParameterSet &cfg) :
  fragModel(cfg.getParameter< std::string >("fragModel")),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("pseudoTop:jets"))),
  prunedGenParticlesToken_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles")))
{
  std::string weights[]={"upFrag","centralFrag","downFrag","PetersonFrag","semilepbrUp","semilepbrDown"};

  //readout weights from file and declare them for the producer
  edm::FileInPath fp = cfg.getParameter<edm::FileInPath>("cfg");
  TFile *fIn=TFile::Open(fp.fullPath().c_str());
  for(size_t i=0; i<sizeof(weights)/sizeof(std::string); i++)
    {
      //produces<edm::ValueMap<float> >(weights[i]);
      TGraph *gr=(TGraph *)fIn->Get(weights[i].c_str());  
      if(gr==0) continue;
      wgtGr_[weights[i]]=gr;
    }
  fIn->Close();

  //Load TFile Service
  edm::Service<TFileService> fs;
  if(!fs)
    throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  //sumary tree for xb and fragmentation re-weighting
  data_ = fs->make<TTree>("FragTree", "FragTree");
  data_->Branch("nB",    &nB_,    "nB/I");
  data_->Branch("fragModel",   model_,  "fragModel[nB]/F");
  data_->Branch("xb",    xb_,    "xb[nB]F");

}

void FragmentationWeightProducer::endJob() {
  data_->Draw("fragModel:xb","xb>0 && xb<2","goff");
  TGraph *g = new TGraph(data_->GetSelectedRows(),data_->GetV2(),data_->GetV1());
  g->SetName(TString(fragModel));
  g->Draw("AP");
  g->Write();
}

void FragmentationWeightProducer::analyze(const edm::Event &evt, const edm::EventSetup &setup) {

  //gen jets
  using namespace edm;
  edm::Handle<std::vector<reco::GenJet>> genJets;
  evt.getByToken(genJetsToken_,genJets);
  //Fragmentation re-weighting
  nB_ = 0;
  for(auto genJet : *genJets) {
    //map the gen particles which are clustered in this jet
    JetFragInfo_t jinfo=analyzeJet(genJet);
    
    //evaluate the weight to an alternative fragmentation model (if a tag id is available)
    if(jinfo.leadTagId != 0) {
      model_[nB_] = wgtGr_[fragModel]->Eval(jinfo.xb);
      xb_[nB_] = jinfo.xb;
      nB_++;
    }
    /*
    else {
      model_[nB_] = 1;
    }
    */
  }

  //Fill ntuple
  if(nB_) data_->Fill();

}


DEFINE_FWK_MODULE(FragmentationWeightProducer);
