/ -*- C++ -*-
//
// Package:    UserCode/MiniAnalyzer
// Class:      MiniAnalyzer
// 
/**\class MiniAnalyzer MiniAnalyzer.cc Test/MiniAnalyzer/plugins/MiniAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Qamar Ul Hassan
//         Created:  Sun, 13 Jul 2014 06:22:18 GMT
//
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "UserCode/TopAnalysis/interface/MiniEvent.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//#include "RecoBTag/PerformanceMeasurements/interface/JetInfoBranches.h"
//#include "UserCode/TopAnalysis/interface/EventInfoBranches.h"
//#include "RecoBTag/PerformanceMeasurements/interface/BookHistograms.h"
//#include "RecoBTag/PerformanceMeasurements/interface/MVAEvaluator.h"

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"

#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>
#include <iostream>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;
using namespace pat; 

//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAnalyzer(const edm::ParameterSet&);
  ~MiniAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  bool doFiducialAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
      
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  //edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<edm::View<pat::Electron>  >  electronToken_;
  //edm::EDGetToken electronToken_; 
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;

  //Electron Decisions
  edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;

  // MVA_ID decisions objects and MVA values and categories (optional)
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumMVAIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightMVAIdMapToken_;  
  
  edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;
                  
  std::unordered_map<std::string,TH1F*> histContainer_;
  std::unordered_map<std::string,TH2F*> histContainer2d_; 

  PFJetIDSelectionFunctor pfjetIDLoose_;

  bool DEBUG_;
  //bool storeEventInfo_;
  
  TTree *tree_;
  MiniEvent_t ev_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//


//
// constructors and destructor
//
MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig) :
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  //electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  eleMediumMVAIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumMVAIdMap"))),
  eleTightMVAIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightMVAIdMap"))),
  mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
  mvaCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap"))),
  pfjetIDLoose_( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE ),
  DEBUG_(iConfig.getUntrackedParameter<bool>("DEBUG",false))
{
  //now do what ever initialization is needed
  electronToken_ = mayConsume<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  //  storeEventInfo_ = iConfig.getParameter<bool>("storeEventInfo");
}


MiniAnalyzer::~MiniAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
bool MiniAnalyzer::doFiducialAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel(genParticleCollectionName_, genParticles);
    
    //require only one lepton (can be from tau, if tau not from hadron)
    int nLeptons(0);
    float lphi(0), leta(0)
    for (size_t i = 0; i < genParticles->size(); ++i) {
      const GenParticle & genIt = (*genParticles)[i];
      if(!genIt->isPromptFinalState() && !genIt->isDirectPromptTauDecayProductFinalState()) continue;
      int ID = abs(genIt.pdgId());
      if(ID!=11 && ID!=13) continue;
      if(genIt->pt()<20 || fabs(genIt->eta())>2.5) continue;
      nLeptons++;
      lphi=genIt->phi();
      leta=genIt->eta();
    }
    if(nLeptons!=1) return false;
    
    //require 2 jets not overlapping with lepton
    int njets(0);
    edm::Handle< std::vector<reco::GenJet> > genJets;
    iEvent.getByLabel("ak4GenJets", genJets);
    for(std::vector<reco::GenJet>::const_iterator genJet=genJets->begin(); genJet!=genJets->end(); genJet++)
     {
       if(genJet->pt()<20 || fabs(genJet->eta())>2.5) continue;
       float dR=deltaR(genJet->eta(),genJet->phi(),leta,lphi);
       if(dR<0.4) continue;
       njets++;
     }
    if(njets<2) return false;
    return true;
}

// ------------ method called for each event  ------------
void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  /* 
     id1_pdfInfo_.clear();
     id2_pdfInfo_.clear();
     x1_pdfInfo_.clear();
     x2_pdfInfo_.clear();
     scalePDF_pdfInfo_.clear();
     ptHat_=0;
     mcWeight_=0;
  */
  bool is_Data = false;
  bool is_MC = false;

  ev_.isFiducial = true; 
  if(iEvent.isRealData()) is_Data = true;
  else {
  	is_MC = true;
        ev_.isFiducial = doFiducialAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup); 
  }
 
  ev_.run     = iEvent.id().run();
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event();
  ev_.isData  = is_Data;
  ev_.isMC    = is_MC;

  TriggerResults tr; 
  edm::Handle<edm::TriggerResults> h_trigRes;
  iEvent.getByToken(triggerBits_, h_trigRes);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  // Get the electron ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream. If you need more info, check with the EGM group.
  
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);

  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);

  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);

  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(eleTightIdMapToken_ ,tight_id_decisions);

  // Get MVA values and categories (optional)
  edm::Handle<edm::ValueMap<bool> > medium_mva_id_decisions;
  iEvent.getByToken(eleMediumMVAIdMapToken_,medium_mva_id_decisions); 
  
  edm::Handle<edm::ValueMap<bool> > tight_mva_id_decisions; 
  iEvent.getByToken(eleTightMVAIdMapToken_,tight_mva_id_decisions);

  edm::Handle<edm::ValueMap<float> > mvaValues;
  iEvent.getByToken(mvaValuesMapToken_,mvaValues); 
  
  edm::Handle<edm::ValueMap<int> > mvaCategories;
  iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories); 

  tr = *h_trigRes;
  std::vector<string> triggerList;
  Service<service::TriggerNamesService> tns;
  histContainer_["cutflow"]->Fill(0);
  histContainer_["ecutflow"]->Fill(0);
  histContainer_["mucutflow"]->Fill(0);
  
  //TRIGGER
  bool mu_trigger = false;
  bool el_trigger = false;
  string electronTrigger;
  string muonTrigger;
  if (is_Data){
    electronTrigger =  "Ele27_eta2p1_WPLoose_Gsf";
    muonTrigger = "IsoMu24_eta2p1";
  }  
  else if (is_MC){
    electronTrigger = "Ele27_eta2p1_WP75_Gsf";
    muonTrigger = "IsoMu24_eta2p1";
  }
 
  bool foundNames = tns->getTrigPaths(tr,triggerList);
  if (!foundNames) std::cout << "Could not get trigger names!\n";
  if (tr.size()!=triggerList.size()) std::cout << "ERROR: length of names and paths not the same: "
					       << triggerList.size() << "," << tr.size() << endl;
  // dump trigger list at first event
  for (unsigned int i=0; i< tr.size(); i++) 
    {
      if( !tr[i].accept() ) continue;
      mu_trigger |= (triggerList[i].find(muonTrigger)!=string::npos);
      el_trigger |= (triggerList[i].find(electronTrigger)!=string::npos);
    }
  if(mu_trigger||el_trigger) histContainer_["cutflow"]->Fill(1);
  else return;
  
  if(mu_trigger) histContainer_["mucutflow"]->Fill(1);
  if(el_trigger) histContainer_["ecutflow"]->Fill(1);
  ev_.muTrigger = mu_trigger;
  ev_.elTrigger = el_trigger;

  //GENERATOR LEVEL INFO
  ev_.ttbar_nw=0;
  if(!is_Data)
    {
      edm::Handle<GenEventInfoProduct> evt;
      iEvent.getByLabel("generator","", evt);
      if(evt.isValid())
	{
	  ev_.ttbar_allmepartons   = evt->nMEPartons();
	  ev_.ttbar_matchmepartons = evt->nMEPartonsFiltered();
	  ev_.ttbar_w[0]           = evt->weight();
	  ev_.ttbar_nw++;
	}
      edm::Handle<LHEEventProduct> evet;
      iEvent.getByLabel("externalLHEProducer","", evet);

      if(evet.isValid())
	{
	  double asdd=evet->originalXWGTUP();
	  for(unsigned int i=0  ; i<evet->weights().size();i++){
	    double asdde=evet->weights()[i].wgt;
	    ev_.ttbar_w[ev_.ttbar_nw]=ev_.ttbar_w[0]*asdde/asdd;
	    ev_.ttbar_nw++;
	  }
	}
    }   
 

  //VERTICES
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &primVtx = vertices->front();
  ev_.nvtx=vertices->size();
  
  //RHO
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  float rho=*rhoH;
  ev_.rho=rho;
	    
  // MUONS
  // cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  histContainer_["nrecomuons"]->Fill(muons->size());
  std::vector<const pat::Muon *> selectedMuons,selectedNonIsoMuons,vetoMuons,nonIsoVetoMuons;        
  for (const pat::Muon &mu : *muons) 
    { 
      //kinematics
      bool passPt( mu.pt() > 30 );
      bool passVetoPt( mu.pt() > 15 );
      bool passEta(fabs(mu.eta()) < 2.1 );
      bool passVetoEta(fabs(mu.eta()) < 2.4 );
    
      //ID
      bool isMedium(muon::isMediumMuon(mu));
      bool isTight(muon::isTightMuon(mu,primVtx));
      float dz(fabs( mu.vertex().z()- primVtx.z())); 

      //isolation
      float relchIso = (mu.chargedHadronIso() + std::max(0., mu.neutralHadronIso() + mu.photonIso() - 0.5*mu.puChargedHadronIso()))/mu.pt(); 
      bool passIso( relchIso  < 0.05 );
      bool passVetoIso( relchIso  < 0.1 );

      //N-1 plots
      if(mu_trigger && passPt && passEta && isTight)
	{
	  histContainer_["muonchreliso"]->Fill(relchIso);
	  histContainer_["muonchiso"]->Fill(mu.chargedHadronIso());
	  histContainer_["muonneuthadiso"]->Fill(mu.neutralHadronIso());
	  histContainer_["muonphotoniso"]->Fill(mu.photonIso());
	  histContainer_["muonpuchiso"]->Fill(mu.puChargedHadronIso());
	  if(passIso)
	    {
	      histContainer_["muondb"]->Fill(mu.dB());
	      histContainer_["muondz"]->Fill(dz);
	    }
	}
      if(mu_trigger &&isTight && passIso)
	{
	  if(passEta) histContainer_["muonpt"]->Fill(mu.pt()); //N-1 plot
	  if(passPt) histContainer_["muoneta"]->Fill(fabs(mu.eta()));	    
	}

      //MUON SELECTION
      if(passPt && passEta) selectedNonIsoMuons.push_back( &mu );
      else if(passVetoPt && passVetoEta && isMedium && passVetoIso) nonIsoVetoMuons.push_back( &mu )
      if(passPt &&  passEta && isTight && passIso)
	selectedMuons.push_back( &mu );	
      else if(passVetoPt && passVetoEta && isMedium && passVetoIso)
	vetoMuons.push_back( &mu );
    }		
  histContainer_["nselmuons"]->Fill(selectedMuons.size());
 
  // ELECTRONS
  // cf. https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
  edm::Handle<edm::View<pat::Electron> >    electrons;
  iEvent.getByToken(electronToken_, electrons);
  histContainer_["nrecoelectrons"]->Fill(electrons->size());
  std::vector<const pat::Electron *> selectedElectrons, selectedNonIsoElectrons, vetoElectrons, nonIsoVetoElectrons;
  Int_t nele(0);
  for (const pat::Electron &el : *electrons) 
    {        
      const auto e = electrons->ptrAt(nele); 
      nele++;

      //kinematics cuts
      bool passPt(e->pt() > 30.0);
      bool passVetoPt(e->pt() > 15.0);
      bool passEta(fabs(e->eta()) < 2.5 && (fabs(e->superCluster()->eta()) < 1.4442 || fabs(e->superCluster()->eta()) > 1.5660));

      //look up id decisions
      bool passVetoId = (*veto_id_decisions)[e];
      //      bool passLooseId = (*loose_id_decisions)[e];
      //bool passMediumId = (*medium_id_decisions)[e];
      bool passTightId  = (*tight_id_decisions)[e];
      //bool passMediumMVAId = (*medium_mva_id_decisions)[e];
      //bool passTightMVAId  = (*tight_mva_id_decisions)[e];

      //N-1 plots
      if(el_trigger && passTightId ){
	if(passEta) histContainer_["electronpt"]->Fill(e->pt());  
	if(passPt)  histContainer_["electroneta"]->Fill(fabs(e->eta()));  
      }
 
      if(el_trigger && passPt && passEta && passTightId){
	histContainer_["electronchiso"]->Fill(e->chargedHadronIso());
	histContainer_["electronneuthadiso"]->Fill(e->neutralHadronIso());
	histContainer_["electronphotoniso"]->Fill(e->photonIso());
	histContainer_["electronpuchiso"]->Fill(e->puChargedHadronIso());
      }

      bool passTightIdNoIso(passTightId);
      if(passPt && passEta && passTightIdNoIso )
	selectedNonIsoElectrons.push_back(&el);
      else if(passVetoPt && passEta && passVetoId)
	nonIsoVetoElectrons.push_back(&el);

      if(passPt && passEta && passTightId )
	selectedElectrons.push_back(&el);
      else if(passVetoPt && passEta && passVetoId)
	vetoElectrons.push_back(&el);
    }
  histContainer_["nselelectrons"]->Fill(selectedElectrons.size());

  //require only 1 tight lepton in the event
  bool passLepton((mu_trigger && selectedMuons.size()==1) || (el_trigger && selectedElectrons.size()==1));
  bool passNonIsoLepton((mu_trigger && selectedNonIsoMuons.size()==1) || (el_trigger && selectedNonIsoElectrons.size()==1));
  if(passLepton) histContainer_["cutflow"]->Fill(2);  
  if(mu_trigger && (selectedMuons.size()==1 || selectedNonIsoMuons.size()==1))
    {
      if(selectedMuons.size()==1) histContainer_["mucutflow"]->Fill(2);
      const pat::Muon *mu=(selectedMuons.size()==1 ? selectedMuons[0] : selectedNonIsoMuons[0]);
      ev_.l_id=(selectedMuons.size()==1 ? 13 :  1300);
      ev_.l_charge=mu->charge();
      ev_.l_pt=mu->pt();
      ev_.l_eta=mu->eta();
      ev_.l_phi=mu->phi();
      ev_.l_tmass=mu->mt();
      ev_.l_chargedHadronIso=mu->chargedHadronIso();
      ev_.l_neutralHadronIso=mu->neutralHadronIso();
      ev_.l_photonIso=mu->photonIso();
      ev_.l_puChargedHadronIso=mu->puChargedHadronIso();
    }
  if(el_trigger && (selectedElectrons.size()==1 || selectedNonIsoElectrons.size()))
    {
      if(selectedElectrons.size()==1)  histContainer_["ecutflow"]->Fill(2);
      const pat::Electron *el=(selectedElectrons.size()==1 ? selectedElectrons[0] : selectedNonIsoElectrons[0]);
      ev_.l_id=(selectedElectrons.size()==1 ? 11 :  1100);
      ev_.l_charge=el->charge();
      ev_.l_pt=el->pt();
      ev_.l_eta=el->eta();
      ev_.l_phi=el->phi();
      ev_.l_tmass=el->mt();
      ev_.l_chargedHadronIso=el->chargedHadronIso();
      ev_.l_neutralHadronIso=el->neutralHadronIso();
      ev_.l_photonIso=el->photonIso();
      ev_.l_puChargedHadronIso=el->puChargedHadronIso();
    }

  //require no other leptons in the event
  bool passVetoLepton(passLepton && (vetoElectrons.size()+vetoMuons.size()));
  bool passNonIsoVetoLepton(passNonIsoLepton && (vetoNonIsoElectrons.size()+vetoNonIsoMuons.size()));
  if(passVeto) histContainer_["cutflow"]->Fill(3);
  if(mu_trigger && passVetoLepton)        histContainer_["mucutflow"]->Fill(3);  
  if(el_trigger && passNonIsoVetoLepton)  histContainer_["ecutflow"]->Fill(3);
  
  // JETS
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  std::vector<const pat::Jet *> selectedJets;
  histContainer_["nrecojets"]->Fill(jets->size());
  ev_.nj=0;
  int nSecVtx=0;
  float pt=0, eta =0.;
  for (const pat::Jet &j : *jets) 
    {  
      //kinematics
      pt = j.pt();
      eta = j.eta();
      if(pt<30 || fabs(eta)>2.5) continue;

      //cross clean to selected iso leptons
      float dR2lepton=9999.
      for(size_t il=0; il<selectedMuons.size(); il++) dR2lepton = TMath::Min(deltaR2lepton,deltaR(j,*(selectedMuons[il])));
      for(size_t il=0; il<selectedElectrons.size(); il++) dR2lepton = TMath::Min(deltaR2lepton,deltaR(j,*(selectedElectrons[il])));
      if(dR2lepton<0.4) continue;
      
      // PF jet ID
      pat::strbitset retpf = pfjetIDLoose_.getBitTemplate();
      retpf.set(false);
      bool passLoose=pfjetIDLoose_( j, retpf );
      if(!passLoose) continue;

      //parton matched to the jet
      bool isLightFromW(false),isB(false);
      const reco::Candidate *genParton = j.genParton();
      if(!is_Data)
	{
	  if(genParton)
	    {
	      isB=(abs(genParton->pdgId())==5);
	      if(genParton->mother())
		isLightFromW=(abs(genParton->mother()->pdgId())==24);
	    }
	}

      float csv=j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      float svtxmass=j.userFloat("vtxMass");
      int vtxNtracks=j.userFloat("vtxNtracks");
      float vtx3DVal=j.userFloat("vtx3DVal");
      float vtx3DSig=j.userFloat("vtx3DSig");

      //control plots
      histContainer_["jetpt"]->Fill(pt);
      if(isB)
	histContainer_["bjetpt"]->Fill(pt);
      else if(isLightFromW)
	histContainer_["lightjetpt"]->Fill(pt);
      else
	histContainer_["otherjetpt"]->Fill(pt);
      histContainer_["jeteta"]->Fill(fabs(j.eta())); 
      histContainer_["jetcsv"]->Fill(csv);
      if(svtxmass>0)
	{
	  nSecVtx++;
	  histContainer_["jetsecvtxmass"]->Fill(svtxmass);
	  histContainer_["jetvtxNtracks"]->Fill(vtxNtracks);
	  histContainer_["jetvtx3DVal"]->Fill(vtx3DVal);
	  histContainer_["jetvtx3DSig"]->Fill(vtx3DSig);
	}
      histContainer_["jetpileupid"]->Fill(j.userFloat("pileupJetId:fullDiscriminant"));

      selectedJets.push_back( &j );
      const reco::GenJet *genJet=j.genJet(); 
      ev_.j_pt[ev_.nj]=j.pt();
      ev_.j_mass[ev_.nj]=j.mass();
      ev_.j_eta[ev_.nj]=j.eta();
      ev_.j_phi[ev_.nj]=j.phi();
      ev_.genj_pt[ev_.nj]=genJet ? genJet->pt() : 0;
      ev_.genj_mass[ev_.nj]=genJet ? genJet->mass() : 0;
      ev_.genj_eta[ev_.nj]=genJet ? genJet->eta() : 0;
      ev_.genj_phi[ev_.nj]=genJet ?  genJet->phi() : 0;
      ev_.j_csv[ev_.nj]=csv;
      ev_.j_vtxmass[ev_.nj]=svtxmass;
      ev_.j_vtxNtracks[ev_.nj]=vtxNtracks;
      ev_.j_vtx3DVal[ev_.nj]=vtx3DVal;
      ev_.j_vtx3DSig[ev_.nj]=vtx3DSig;
      ev_.j_puid[ev_.nj]=j.userFloat("pileupJetId:fullDiscriminant");
      ev_.j_flav[ev_.nj]=j.partonFlavour();
      ev_.j_pid[ev_.nj]=genParton ? genParton->pdgId() : 0;
      ev_.nj++;
    }
  if(passVetoLepton) histContainer_["nseljets"]->Fill(selectedJets.size());

  //
  // MET
  //
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  float metpt = mets->at(0).pt();
  float metphi = mets->at(0).phi();
  float dphi_met_lepton = deltaPhi(ev_.l_phi, metphi); // use the function to restrict to the 0,pi range
  float mt=sqrt(2*ev_.l_pt*metpt*(1-cos(dphi_met_lepton)));
  //charged met
  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);
  TLorentzVector chMet(0,0,0,0);
  float chHT(0);
  for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
    //require not to be associated to other PVs and to be charged
    const pat::PackedCandidate &pf = (*pfs)[i];
    if (pf.fromPV() == 0 || pf.charge()== 0) continue;
    chMet -= TLorentzVector(pf.px(),pf.py(),0,pf.pt());
    chHT += pf.pt();
  }
  float dphi_chmet_lepton = deltaPhi(ev_.l_phi, chMet.Phi()); // use the function to restrict to the 0,pi range
  float chmt=sqrt(2*ev_.l_pt*chMet.Pt()*(1-cos(dphi_chmet_lepton)));

  //save to ttree
  ev_.met_pt=metpt;
  ev_.met_phi=metphi;
  ev_.mt=mt;
  ev_.chmet_pt=chMet.Pt();
  ev_.chmet_phi=chMet.Phi();
  ev_.chmt=chmt;

  //save tree if event is interesting
  if((passVetoLeptons || passNonIsoVetoLeptons) && selectedJets.size()>=2) tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzer::beginJob(){
  edm::Service<TFileService> fs;
  histContainer_["cutflow"]   = fs->make<TH1F>("cutflow",    ";Selection cut;Events", 5, 0., 5.); 
  histContainer_["ecutflow"]  = fs->make<TH1F>("ecutflow",   ";Selection cut;Events", 5, 0., 5.); 
  histContainer_["mucutflow"] = fs->make<TH1F>("mucutflow",  ";Selection cut;Events", 5, 0., 5.); 
  TString steps[]={"N(Reco.)","N(Trig.)","N(=1 good lep)","N(=0 loose lep)","N(#geq2 jets)"};
  for(size_t i=0; i<sizeof(steps)/sizeof(TString); i++){
    histContainer_["cutflow"]   -> GetXaxis()->SetBinLabel(i+1,steps[i]);
    histContainer_["ecutflow"]  -> GetXaxis()->SetBinLabel(i+1,steps[i]);
    histContainer_["mucutflow"] -> GetXaxis()->SetBinLabel(i+1,steps[i]);
  }


  histContainer_["nrecomuons"]     = fs->make<TH1F>("nrecomuons",      ";# reconstructed muons; Events",             10, 0., 10.);
  histContainer_["muonpt"]         = fs->make<TH1F>("muonpt",          ";Transverse momentum [GeV];# muons",         100, 0., 300.);
  histContainer_["muoneta"]        = fs->make<TH1F>("muoneta",         ";Pseudo-rapidity;#muons ",                   100, 0, 3.);
  histContainer_["muondb"]         = fs->make<TH1F>("muondb",          ";d_{0} [cm];# muons",         100, 0.,0.3);
  histContainer_["muondz"]         = fs->make<TH1F>("muondz",          ";|d_{z}| [cm];# muons",       100, 0.,0.6);
  histContainer_["muonchreliso"]   = fs->make<TH1F>("muonchreliso",       ";Relative charged hadron isolation;# muons ",  100, 0, 0.1);
  histContainer_["muonchiso"]      = fs->make<TH1F>("muonchiso",       ";Charged hadron isolation [GeV];# muons ",  100, 0, 10);
  histContainer_["muonpuchiso"]    = fs->make<TH1F>("muonpuchiso",     ";Pileup charged hadron isolation [GeV];# muons ",  100, 0, 10);
  histContainer_["muonneuthadiso"] = fs->make<TH1F>("muonneuthadiso",  ";Neutral hadron isolation [GeV];# muons ",  100, 0, 10);
  histContainer_["muonphotoniso"]  = fs->make<TH1F>("muonphotoniso",   ";Photon isolation [GeV];# muons ",  100, 0, 10);
  histContainer_["nselmuons"]      = fs->make<TH1F>("nselmuons",       ";# reconstructed muons;Events",              5, 0.,5.);

  histContainer_["nrecoelectrons"]    = fs->make<TH1F>("nrecoelectrons",  ";# reconstructed electrons;Events",         5, 0., 5.);
  histContainer_["electronpt"]        = fs->make<TH1F>("electronpt",      ";Transverse momentum [GeV];# electrons",    100, 0., 300.);
  histContainer_["electroneta"]       = fs->make<TH1F>("electroneta",     ";Pseudo-rapidity;#electrons",               100, 0., 3.);
  histContainer_["electronchiso"]      = fs->make<TH1F>("electronchiso",       ";Charged hadron isolation [GeV];# electrons ",  100, 0, 10);
  histContainer_["electronpuchiso"]    = fs->make<TH1F>("electronpuchiso",     ";Pileup charged hadron isolation [GeV];# electrons ",  100, 0, 10);
  histContainer_["electronneuthadiso"] = fs->make<TH1F>("electronneuthadiso",  ";Neutral hadron isolation [GeV];# electrons ",  100, 0, 10);
  histContainer_["electronphotoniso"]  = fs->make<TH1F>("electronphotoniso",   ";Photon isolation [GeV];# electrons ",  100, 0, 10);
  histContainer_["nselelectrons"]     = fs->make<TH1F>("nselelectrons",   ";# selected electrons;Events", 5, 0., 5.);

  histContainer_["nrecojets"]      = fs->make<TH1F>("nrecojets",   ";#reconstructed jets;Events", 50, 0., 50.);
  histContainer_["jetpt"]     = fs->make<TH1F>("jetpt",       ";Transverse momentum [GeV];# jets", 100, 0., 250.);
  histContainer_["bjetpt"]     = fs->make<TH1F>("bjetpt",       ";Transverse momentum [GeV];# jets", 100, 0., 250.);
  histContainer_["lightjetpt"]     = fs->make<TH1F>("lightjetpt",       ";Transverse momentum [GeV];# jets", 100, 0., 250.);
  histContainer_["otherjetpt"]     = fs->make<TH1F>("otherjetpt",       ";Transverse momentum [GeV];# jets", 100, 0., 250.);
  histContainer_["jeteta"]         = fs->make<TH1F>("jeteta",      ";Pseudo-rapidity;# jets", 100, 0., 3.);
  histContainer_["jetcsv"]         = fs->make<TH1F>("jetcsv",      ";Combined secondary vertes;# jets", 100, -1.2, 1.2);
  histContainer_["jetpileupid"]    = fs->make<TH1F>("jetpileupid", ";Pileup jet id;#jets", 100, -1.2, 1.2);
  histContainer_["jetsecvtxmass"]  = fs->make<TH1F>("jetvtxMass", ";Secondary vertex mass [GeV];#jets", 100, 0., 6.);
  histContainer_["jetvtxNtracks"]  = fs->make<TH1F>("jetvtxNtracks", ";Vertex Tracks;#jets", 6, 0., 6.);  
  histContainer_["jetvtx3DVal"]    = fs->make<TH1F>("jetvtx3DVal", ";vtx3DVal [cm];#jets", 100, -5., 5.);
  histContainer_["jetvtx3DSig"]    = fs->make<TH1F>("jetvtx3DSig", ";vtx3DSig;#jets", 100, 0., 5.);
  histContainer_["nseljets"]       = fs->make<TH1F>("nseljets",    ";#selected jets;Events", 10, 0., 10.);
  }


  //instruct ROOT to compute the uncertainty from the square root of weights
  //http://root.cern.ch/root/html/TH1.html#TH1:Sumw2
  for(std::unordered_map<std::string,TH1F*>::iterator it=histContainer_.begin();   it!=histContainer_.end();   it++) it->second->Sumw2();
  
  //create a tree for the selected events
  tree_ = fs->make<TTree>("AnaTree", "AnaTree");
  createMiniEventTree(tree_,ev_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  MiniAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  MiniAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  MiniAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  MiniAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
