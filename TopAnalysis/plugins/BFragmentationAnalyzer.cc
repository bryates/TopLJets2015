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
#include "TopLJets2015/TopAnalysis/interface/FragEvent.h"

#include <memory>
#include <string>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//#include "Utilities/General/interface/FileInPath.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "TTree.h"
#include "TGraph.h"

using namespace std;
#define IS_BQUARK_PDGID(id) ( abs(id) == 5 )


class FragmentationAnalyzer : public edm::EDAnalyzer {
public:
  explicit FragmentationAnalyzer(const edm::ParameterSet&);
  ~FragmentationAnalyzer();  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  virtual void endRun(const edm::Run&,const edm::EventSetup&);  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void genAnalysis(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;

  std::vector<int> hadronList_;
  std::map<std::string, TH1F *> histos_;
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken_;
  std::vector<int> numEntries_;
  /*
  TTree *data_;
  FragEvent_t ev_;
  */
};


//
FragmentationAnalyzer::FragmentationAnalyzer(const edm::ParameterSet& iConfig) :
  hadronList_(iConfig.getParameter<std::vector<int> >("hadronList")),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("pseudoTop:jets"))),
  numEntries_(iConfig.getParameter<std::vector<int> >("numEntries"))
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
      histos_["xb_semilep"+name] = fs->make<TH1F>(("xb_semilep"+name).c_str(), (name+";x_{b}=p_{T}(B)/p_{T}(jet); Jets").c_str(), 100, 0, 2);
    }
  for(auto it : histos_) it.second->Sumw2();
  //produces<edm::ValueMap<float> >("xb");

  //Load TFile Service
  if(!fs)
    throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  //sumary tree for xb and fragmentation re-weighting
  //data_ = fs->make<TTree>("FragTree", "FragTree");
  /*
  data_ = fs->make<TTree>("fragTree","fragTree");
  createFragEventTree(data_,ev_);
  data_->Branch("nB",    &nB_,    "nB/I");
  data_->Branch("xb",    xb_,    "xb[nB]/F");
  data_->Branch("xbc",   xbc_,   "xbc[nB]/F");
  data_->Branch("id",    id_,    "id[nB]/F");
  data_->Branch("pt",    pt_,    "pt[nB]/F");
  */
}


//
FragmentationAnalyzer::~FragmentationAnalyzer()
{
}


//
void FragmentationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  genAnalysis(iEvent, iSetup);
  if(histos_["xb_semilepinc"]->GetEntries() > numEntries_[0]) endJob();
  if(histos_["xb_semilepinc"]->GetEntries() > numEntries_[0]) return;
}
void FragmentationAnalyzer::genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //
  /*
  ev_.run     = iEvent.id().run();
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event(); 
  */

  edm::Handle<std::vector<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsToken_,genJets);  
  /*
  ev_.nB = 0;
  std::map< std::string, std::vector<float> > jetWeights;
  jetWeights["xb"]=std::vector<float>();
  */
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
          //xb for semi-leptonic ttbar decay
          histos_["xb_semilepinc"]->Fill(jinfo.xb);
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
      /*
      if(jinfo.hasCharm) {
        ev_.xb[ev_.nB] = jinfo.xb;
        //jetWeights["xb"].push_back(jinfo.xb);
        ev_.xbc[ev_.nB] = jinfo.xb_charged;
        ev_.id[ev_.nB] = jinfo.charmId;
        ev_.pt[ev_.nB] = jinfo.pt;
        ev_.nB++;
      }
      */
    }
  //Fill ntuple
  //if(ev_.nB) data_->Fill();

  //put in event
  /*
  for(auto it : jetWeights) {
    auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    typename edmValueMap<float>::Filler filler(*valMap);
    filler.insert(genJets, it.second.begin(), it.second.end());
    filler.fill();
    iEvent.put(valMap, it.first);
  }
  */
}

//
void FragmentationAnalyzer::beginJob(){ }

//
void  FragmentationAnalyzer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup){ }

//
void FragmentationAnalyzer::endJob() {
  /*
  data_->Draw("xb:pt","","goff");
  TGraph *g = new TGraph(data_->GetSelectedRows(),data_->GetV2(),data_->GetV1());
  g->SetName("xb");
  g->SetTitle("x_{b}");
  g->Draw("AP");
  g->Write();
  */
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
