#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "TopLJets2015/TopAnalysis/interface/BFragmentationAnalyzerUtils.h"

#include "TTree.h"
#include "TGraph.h"

using namespace std;


class FragmentationAnalyzer : public edm::EDAnalyzer {
public:
  explicit FragmentationAnalyzer(const edm::ParameterSet&);
  ~FragmentationAnalyzer();  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  virtual void endRun(const edm::Run&,const edm::EventSetup&);  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  std::vector<int> hadronList_;
  std::map<std::string, TH1F *> histos_;
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken_;
  TTree *data_;
  Int_t nB_;
  Float_t xb_[100], model_[100];
};


//
FragmentationAnalyzer::FragmentationAnalyzer(const edm::ParameterSet& iConfig) :
  hadronList_(iConfig.getParameter<std::vector<int> >("hadronList")),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("pseudoTop:jets")))
{
  //prepare monitoring histograms
  size_t nhadrons=hadronList_.size()+1;
  histos_["semilepbr"]    = fs->make<TH1F>("semilepbr", ";B-hadron;BR",nhadrons,0,nhadrons);
  histos_["semilepbrinc"] = fs->make<TH1F>("semilepbrinc", ";B-hadron;BR",nhadrons,0,nhadrons);
  for(size_t i=0; i<nhadrons; i++)
    {
      std::string name("inc");
      if(i) { char buf[20]; sprintf(buf,"%d", hadronList_[i-1]); name=buf; }
      histos_["semilepbr"]->GetXaxis()->SetBinLabel(i+1,name.c_str());
      histos_["semilepbrinc"]->GetXaxis()->SetBinLabel(i+1,name.c_str());
      histos_["xb_"+name] = fs->make<TH1F>(("xb_"+name).c_str(), (name+";x_{b}=p_{T}(B)/p_{T}(jet); Jets").c_str(), 100, 0, 2);
    }
  for(auto it : histos_) it.second->Sumw2();

  //Load TFile Service
  if(!fs)
    throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  //sumary tree for xb and fragmentation re-weighting
  data_ = fs->make<TTree>("FragTree", "FragTree");
  data_->Branch("nB",    &nB_,    "nB/I");
  data_->Branch("fragModel",   model_,  "fragModel[nB]/F");
  data_->Branch("xb",    xb_,    "xb[nB]F");
}


//
FragmentationAnalyzer::~FragmentationAnalyzer()
{
}


//
void FragmentationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //
  edm::Handle<std::vector<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsToken_,genJets);  
  nB_ = 0;
  for(auto genJet : *genJets)
    {
      //map the gen particles which are clustered in this jet
      JetFragInfo_t jinfo=analyzeJet(genJet);

      //skip if not B
      int absid(abs(jinfo.leadTagId));
      if(!IS_BHADRON_PDGID(absid)) continue;

      //inclusive histos
      if(jinfo.hasSemiLepDecay) 
	{
	  if(!jinfo.hasTauSemiLepDecay)
	    histos_["semilepbr"]->Fill(0);
	  histos_["semilepbrinc"]->Fill(0);
	}
      histos_["xb_inc"]->Fill(jinfo.xb);
      
      //exclusive histograms
      std::vector<int>::iterator hit=std::find(hadronList_.begin(), hadronList_.end(), absid);
      if(hit!=hadronList_.end())
	{
	  char buf[20]; 
	  sprintf(buf,"xb_%d", *hit);
	  if(jinfo.hasSemiLepDecay)
	    {
	      if(!jinfo.hasTauSemiLepDecay) 
		histos_["semilepbr"]->Fill(1+hit-hadronList_.begin());
	      histos_["semilepbrinc"]->Fill(1+hit-hadronList_.begin());
	    }
	  histos_[buf]->Fill(jinfo.xb);
	}
      xb_[nB_] = jinfo.xb_charged;
      nB_++;
    }
  //Fill ntuple
  if(nB_) data_->Fill();
}

//
void FragmentationAnalyzer::beginJob(){ }

//
void  FragmentationAnalyzer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup){ }

//
void FragmentationAnalyzer::endJob(){
  data_->Draw("xb","xb>0 && xb<2","goff");
  TGraph *g = new TGraph(data_->GetSelectedRows(),data_->GetV2(),data_->GetV1());
  g->SetName("xb_{charged}");
  g->Draw("AP");
  g->Write();
}

//
void FragmentationAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FragmentationAnalyzer);
