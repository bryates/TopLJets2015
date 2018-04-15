#include <memory>
#include <string>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Utilities/General/interface/FileInPath.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "TopLJets2015/TopAnalysis/interface/BFragmentationAnalyzerUtils.h"

#include "TFile.h"
#include "TGraph.h"

using namespace std;


class BFragmentationProducer : public edm::stream::EDProducer<> {
   public:
      explicit BFragmentationProducer(const edm::ParameterSet&);
      ~BFragmentationProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

  edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken_;
  std::map<std::string, TGraph *> wgtGr_;
};

//
BFragmentationProducer::BFragmentationProducer(const edm::ParameterSet& iConfig):
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("pseudoTop:jets")))
{
  std::string weights[]={"upFrag","downFrag","semilepbrUp","semilepbrDown"};
  //std::string weights[]={"upFrag","centralFrag","downFrag","PetersonFrag","semilepbrUp","semilepbrDown"};

  //readout weights from file and declare them for the producer
  edm::FileInPath fp = iConfig.getParameter<edm::FileInPath>("cfg");
  TFile *fIn=TFile::Open(fp.fullPath().c_str());
  produces<edm::ValueMap<float> >("xb");
  for(size_t i=0; i<sizeof(weights)/sizeof(std::string); i++)
    {
      produces<edm::ValueMap<float> >(weights[i]);
      TGraph *gr=(TGraph *)fIn->Get(weights[i].c_str());  
      if(gr==0) continue;
      wgtGr_[weights[i]]=gr;
    }
  fIn->Close();
}

//
BFragmentationProducer::~BFragmentationProducer()
{
}

//
void BFragmentationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<std::vector<reco::GenJet> > genJets;
   iEvent.getByToken(genJetsToken_,genJets);
   std::map< std::string, std::vector<float> > jetWeights;
   jetWeights["upFrag"]=std::vector<float>();
   //jetWeights["centralFrag"]=std::vector<float>();
   jetWeights["downFrag"]=std::vector<float>();
   //jetWeights["PetersonFrag"]=std::vector<float>();
   jetWeights["semilepbrUp"]=std::vector<float>();
   jetWeights["semilepbrDown"]=std::vector<float>();
   for(auto genJet : *genJets)
     {
       //map the gen particles which are clustered in this jet
       JetFragInfo_t jinfo=analyzeJet(genJet);
       
       //evaluate the weight to an alternative fragmentation model (if a tag id is available)
       if(jinfo.leadTagId != 0)
       {
        jetWeights["upFrag"].push_back(wgtGr_["upFrag"]->Eval(jinfo.xb));
        //jetWeights["centralFrag"].push_back(wgtGr_["centralFrag"]->Eval(jinfo.xb));
        jetWeights["downFrag"].push_back(wgtGr_["downFrag"]->Eval(jinfo.xb));
        //jetWeights["PetersonFrag"].push_back(wgtGr_["PetersonFrag"]->Eval(jinfo.xb));
       }
       else
       {
        jetWeights["upFrag"].push_back(1.);
        //jetWeights["centralFrag"].push_back(1.);
        jetWeights["downFrag"].push_back(1.);
        //jetWeights["PetersonFrag"].push_back(1.);
       }

       float semilepbrUp(1.0),semilepbrDown(1.0);
       int absBid(abs(jinfo.leadTagId));
       if(absBid==511 || absBid==521 || absBid==531 || absBid==5122)
	   {
	    int bid( jinfo.hasSemiLepDecay ? absBid : -absBid);
	    semilepbrUp=wgtGr_["semilepbrUp"]->Eval(bid);
	    semilepbrDown=wgtGr_["semilepbrDown"]->Eval(bid);
	   }
       jetWeights["semilepbrUp"].push_back(semilepbrUp);
       jetWeights["semilepbrDown"].push_back(semilepbrDown);

       /*
       if(IS_BHADRON_PDGID(absBid))
	 cout << genJet.pt()       << " GeV jet matched with " << absBid 
	      << " upFrag: "       << wgtGr_["upFrag"]->Eval(jinfo.xb)
	      << " centralFrag: "  << wgtGr_["centralFrag"]->Eval(jinfo.xb)
	      << " downFrag: "     << wgtGr_["downFrag"]->Eval(jinfo.xb)
	      << " PetersonFrag: " << wgtGr_["PetersonFrag"]->Eval(jinfo.xb)
	      << " semilep+: "     << semilepbrUp 
	      << " semilep-: "     << semilepbrDown << endl;
       */
     }

   //put in event
   for(auto it : jetWeights)
     {
       auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
       typename edm::ValueMap<float>::Filler filler(*valMap);
       filler.insert(genJets, it.second.begin(), it.second.end());
       filler.fill();
       iEvent.put(valMap, it.first);
     }
}

//
void BFragmentationProducer::beginStream(edm::StreamID)
{
}

//
void BFragmentationProducer::endStream() {
}

//
void BFragmentationProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BFragmentationProducer);
