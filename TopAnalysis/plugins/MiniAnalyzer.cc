//
// -*- C++ -*-
//
// Package:    TopLJets2015/TopAnalysis
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
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
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
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/MyIPTools.h"
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

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

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
  virtual void endRun(const edm::Run & iRun, edm::EventSetup const & iSetup);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  Int_t doFiducialAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;


  // member data 
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron>  >  electronToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > qgToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;

  //Electron Decisions
  edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleTightIdFullInfoMapToken_;
                  
  std::unordered_map<std::string,TH1F*> histContainer_;

  PFJetIDSelectionFunctor pfjetIDLoose_;

  std::vector<std::string> muTriggersToUse_, elTriggersToUse_;

  //ttbar event classifier
  edm::EDGetTokenT<int> genTtbarIdToken_;
  
  TTree *tree_;
  MiniEvent_t ev_;

  edm::Service<TFileService> fs;

  bool saveTree_;
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
  jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
  qgToken_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  eleTightIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleTightIdFullInfoMap"))),
  pfjetIDLoose_( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE ),
  genTtbarIdToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("genTtbarId"))),
  saveTree_( iConfig.getParameter<bool>("saveTree") )
{
  //now do what ever initialization is needed
  electronToken_ = mayConsume<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  elTriggersToUse_ = iConfig.getParameter<std::vector<std::string> >("elTriggersToUse");
  muTriggersToUse_ = iConfig.getParameter<std::vector<std::string> >("muTriggersToUse");

  for(Int_t igenjet=0; igenjet<5; igenjet++)
    {
      TString tag("fidcounter"); tag+=igenjet;
      histContainer_[tag.Data()] = fs->make<TH1F>(tag,    ";Variation;Events", 1000, 0., 1000.); 
    }
  histContainer_["counter"]   = fs->make<TH1F>("counter",    ";Counter;Events",2,0,2);
  for(std::unordered_map<std::string,TH1F*>::iterator it=histContainer_.begin();   it!=histContainer_.end();   it++) it->second->Sumw2();

  //create a tree for the selected events
  if(saveTree_)
    {
      tree_ = fs->make<TTree>("data","data");
      createMiniEventTree(tree_,ev_);
    }
}


MiniAnalyzer::~MiniAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
Int_t MiniAnalyzer::doFiducialAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("prunedGenParticles", genParticles);
    
    //require only one lepton (can be from tau, if tau not from hadron)
    int nLeptons(0);
    float lphi(0), leta(0);
    for (size_t i = 0; i < genParticles->size(); ++i) {
      const GenParticle & genIt = (*genParticles)[i];

      if(genIt.isHardProcess())
	{
	  ev_.ghp_id[ ev_.ngenHardProc ] = genIt.pdgId();
	  ev_.ghp_pt[ ev_.ngenHardProc ] = genIt.pt();
	  ev_.ghp_eta[ ev_.ngenHardProc ] = genIt.eta();
	  ev_.ghp_phi[ ev_.ngenHardProc ] = genIt.phi();
	  ev_.ghp_m[ ev_.ngenHardProc ] = genIt.mass();
	  ev_.ngenHardProc++;
	}

      if(!genIt.isPromptFinalState() && !genIt.isDirectPromptTauDecayProductFinalState()) continue;
      int ID = abs(genIt.pdgId());
      if(ID!=11 && ID!=13) continue;
      if(genIt.pt()<20 || fabs(genIt.eta())>2.5) continue;
      nLeptons++;
      lphi=genIt.phi();
      leta=genIt.eta();
    }
    if(nLeptons!=1) return 0;
    
    //require 1 jets not overlapping with lepton
    edm::Handle< std::vector<reco::GenJet> > genJets;
    iEvent.getByLabel("ak4GenJetsCustom", genJets);
    for(std::vector<reco::GenJet>::const_iterator genJet=genJets->begin(); genJet!=genJets->end(); genJet++)
     {
       if(genJet->pt()<20 || fabs(genJet->eta())>2.5) continue;
       float dR=deltaR(genJet->eta(),genJet->phi(),leta,lphi);
       if(dR<0.4) continue;
       ev_.genj_pt[ev_.ngenj]=genJet->pt();
       ev_.genj_eta[ev_.ngenj]=genJet->eta();
       ev_.genj_phi[ev_.ngenj]=genJet->phi();
       ev_.genj_mass[ev_.ngenj]=genJet->mass();
       ev_.ngenj++;
     }
    return ev_.ngenj;
}

// ------------ method called for each event  ------------
void MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //count events analyzed
  histContainer_["counter"]->Fill(0);

  //EVENT HEADER
  ev_.run     = iEvent.id().run();
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event(); 
  ev_.isData=iEvent.isRealData();

  //VERTICES
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &primVtx = vertices->front();
  reco::VertexRef primVtxRef(vertices,0);

  ev_.nvtx=vertices->size();
  if(ev_.nvtx==0) return;
  
  //RHO
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  float rho=*rhoH;
  ev_.rho=rho;

  //GENERATOR LEVEL INFO
  ev_.isFiducial = true;  
  ev_.ttbar_nw=0;
  ev_.me_np=0;
  ev_.ngenj=0;
  ev_.ttbar_genId=0;
  ev_.ngenHardProc=0;
  if(!ev_.isData)
    {
      edm::Handle<int> genTtbarIdHandle;
      iEvent.getByToken(genTtbarIdToken_, genTtbarIdHandle);
      if(genTtbarIdHandle.isValid()) ev_.ttbar_genId=*genTtbarIdHandle;

      Int_t ngenJets = doFiducialAnalysis(iEvent,iSetup);
      if(ngenJets<1) ev_.isFiducial=false;

      edm::Handle<GenEventInfoProduct> evt;
      iEvent.getByLabel("generator","", evt);
      if(evt.isValid())
	{
	  ev_.ttbar_allmepartons   = evt->nMEPartons();
	  ev_.ttbar_matchmepartons = evt->nMEPartonsFiltered();
	  ev_.ttbar_w[0]           = evt->weight();
	  ev_.ttbar_nw++;
	}
      histContainer_["counter"]->Fill(1,ev_.ttbar_w[0]);

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
	  
	  const lhef::HEPEUP &hepeup=evet->hepeup();
	  ev_.me_id=hepeup.IDPRUP;
	  for(int ip=0;ip<hepeup.NUP; ip++)
	    {
	      ev_.me_pid[ev_.me_np]=hepeup.IDUP[ip];
	      ev_.me_px[ev_.me_np]=hepeup.PUP[ip][0];
	      ev_.me_py[ev_.me_np]=hepeup.PUP[ip][1];
	      ev_.me_pz[ev_.me_np]=hepeup.PUP[ip][2];
	      ev_.me_mass[ev_.me_np]=hepeup.PUP[ip][4];
	      ev_.me_np++;
	    }

	}
      
      for(Int_t igenjet=0; igenjet<5; igenjet++)
	{
	  TString tag("fidcounter"); tag+=igenjet;
	  histContainer_[tag.Data()]->Fill(0.,ev_.ttbar_w[0]);
	  if(igenjet<=ngenJets)
	    {
	      for(Int_t iw=1; iw<ev_.ttbar_nw; iw++)
		histContainer_[tag.Data()]->Fill((float)iw,ev_.ttbar_w[iw]);
	    }
	}

      edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
      iEvent.getByLabel("slimmedAddPileupInfo", PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator ipu;
      for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) {
	if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
	ev_.pu=ipu->getPU_NumInteractions();
	ev_.putrue=ipu->getTrueNumInteractions();
      }
    }
 
  //TRIGGER INFORMATION
  edm::Handle<edm::TriggerResults> h_trigRes;
  iEvent.getByToken(triggerBits_, h_trigRes);
  std::vector<string> triggerList;
  Service<service::TriggerNamesService> tns;
  tns->getTrigPaths(*h_trigRes,triggerList);
  ev_.muTrigger=0;
  ev_.elTrigger=0;
  for (unsigned int i=0; i< h_trigRes->size(); i++) 
    {
      //std::cout << triggerList[i] << std::endl;
      if( !(*h_trigRes)[i].accept() ) continue;
      for(size_t imu=0; imu<muTriggersToUse_.size(); imu++)
	{
	  if (triggerList[i].find(muTriggersToUse_[imu])==string::npos) continue;
	  ev_.muTrigger |= (1 << imu);
	}
      for(size_t iel=0; iel<elTriggersToUse_.size(); iel++) 
	{
	  if (triggerList[i].find(elTriggersToUse_[iel])==string::npos)continue;
	  ev_.elTrigger |= (1 << iel);
	}
    }
  bool passMuTrigger(ev_.isData ? ev_.muTrigger!=0 : true);
  bool passElTrigger(ev_.isData ? ev_.elTrigger!=0 : true);  
  if(!passMuTrigger && !passElTrigger) return;
 

  //MUON SELECTION: cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  std::vector<const pat::Muon *> selectedMuons,selectedNonIsoMuons,vetoMuons,vetoNonIsoMuons;        
  for (const pat::Muon &mu : *muons) 
    { 
      //kinematics
      bool passPt( mu.pt() > 30 );
      bool passVetoPt( mu.pt() > 10 );
      bool passEta(fabs(mu.eta()) < 2.1 );
      bool passVetoEta(fabs(mu.eta()) < 2.4 );
    
      //ID
      bool isMedium(muon::isMediumMuon(mu));
      bool isTight(muon::isTightMuon(mu,primVtx));

      //isolation
      //float relchIso = mu.chargedHadronIso()/mu.pt(); 
      float relIsoDeltaBeta( (mu.pfIsolationR04().sumChargedHadronPt + max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt))/mu.pt() ); 
      bool passIso( relIsoDeltaBeta  < 0.15 );
      bool passVetoIso( relIsoDeltaBeta  < 0.25 );

      //MUON SELECTION
      if(passPt && passEta)                                           selectedNonIsoMuons.push_back( &mu );
      else if(passVetoPt && passVetoEta && isMedium && passVetoIso)   vetoNonIsoMuons.push_back( &mu );
      if(passPt &&  passEta && isTight && passIso)                    selectedMuons.push_back( &mu );	
      else if(passVetoPt && passVetoEta && isMedium && passVetoIso)   vetoMuons.push_back( &mu );
    }
  
  // ELECTRON SELECTION: cf. https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(eleTightIdMapToken_ ,tight_id_decisions);
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > tight_id_cutflow_data;
  iEvent.getByToken(eleTightIdFullInfoMapToken_,tight_id_cutflow_data);
  edm::Handle<edm::View<pat::Electron> >    electrons;
  iEvent.getByToken(electronToken_, electrons);  
  std::vector<const pat::Electron *> selectedElectrons, selectedNonIsoElectrons, vetoElectrons, vetoNonIsoElectrons;
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
      bool passTightId  = (*tight_id_decisions)[e];
      vid::CutFlowResult fullCutFlowData = (*tight_id_cutflow_data)[e];
      bool passTightIdNoIso(true);
      int ncuts = fullCutFlowData.cutFlowSize();
      for(int icut = 0; icut<ncuts; icut++)
	{
	  if(icut==9 && fullCutFlowData.getCutResultByIndex(icut) ) passTightIdNoIso=false;
	  if(icut!=9 && !fullCutFlowData.getCutResultByIndex(icut)) passTightIdNoIso=false;
	}
      
      if(passPt && passEta && passTightIdNoIso )   selectedNonIsoElectrons.push_back(&el);
      else if(passVetoPt && passEta && passVetoId) vetoNonIsoElectrons.push_back(&el);

      if(passPt && passEta && passTightId )        selectedElectrons.push_back(&el);
      else if(passVetoPt && passEta && passVetoId) vetoElectrons.push_back(&el);
    }

  //require only 1 tight lepton in the event
  bool passLepton((passMuTrigger && selectedMuons.size()==1) || (passElTrigger && selectedElectrons.size()==1));
  bool passNonIsoLepton((passMuTrigger && selectedNonIsoMuons.size()==1) || (passElTrigger && selectedNonIsoElectrons.size()==1));
  if(!passLepton && !passNonIsoLepton) return;
  if(passMuTrigger && (selectedMuons.size()==1 || selectedNonIsoMuons.size()==1))
    {
      const pat::Muon *mu=(selectedMuons.size()==1 ? selectedMuons[0] : selectedNonIsoMuons[0]);
      const reco::GenParticle * gen=mu->genLepton(); 
      ev_.isPromptFinalState = gen ? gen->isPromptFinalState() : false;
      ev_.isDirectPromptTauDecayProductFinalState = gen ? gen->isDirectPromptTauDecayProductFinalState() : false;
      ev_.l_id=(selectedMuons.size()==1 ? 13 :  1300);
      ev_.l_charge=mu->charge();
      ev_.l_pt=mu->pt();
      ev_.l_eta=mu->eta();
      ev_.l_phi=mu->phi();
      ev_.l_mass=mu->mass();
      ev_.l_chargedHadronIso=mu->pfIsolationR04().sumChargedHadronPt;
      ev_.l_neutralHadronIso=mu->pfIsolationR04().sumNeutralHadronEt;
      ev_.l_photonIso=mu->pfIsolationR04().sumPhotonEt;
      ev_.l_puChargedHadronIso=mu->pfIsolationR04().sumPUPt;
      ev_.l_ip3d = -9999.;
      ev_.l_ip3dsig = -9999;
      if(mu->innerTrack().get())
	{
	  std::pair<bool,Measurement1D> ip3dRes = getImpactParameter<reco::TrackRef>(mu->innerTrack(), primVtxRef, iSetup, true);
	  ev_.l_ip3d    = ip3dRes.second.value();
	  ev_.l_ip3dsig = ip3dRes.second.significance();
	}
    }
  if(passElTrigger && (selectedElectrons.size()==1 || selectedNonIsoElectrons.size()))
    {     
      const pat::Electron *el=(selectedElectrons.size()==1 ? selectedElectrons[0] : selectedNonIsoElectrons[0]);
      const reco::GenParticle * gen=el->genLepton(); 
      ev_.isPromptFinalState = gen ? gen->isPromptFinalState() : false;
      ev_.isDirectPromptTauDecayProductFinalState = gen ? gen->isDirectPromptTauDecayProductFinalState() : false;
      ev_.l_id=(selectedElectrons.size()==1 ? 11 :  1100);
      ev_.l_charge=el->charge();
      ev_.l_pt=el->pt();
      ev_.l_eta=el->eta();
      ev_.l_phi=el->phi();
      ev_.l_mass=el->mass();
      ev_.l_chargedHadronIso=el->chargedHadronIso();
      ev_.l_neutralHadronIso=el->neutralHadronIso();
      ev_.l_photonIso=el->photonIso();
      ev_.l_puChargedHadronIso=el->puChargedHadronIso();
      ev_.l_ip3d = -9999.;
      ev_.l_ip3dsig = -9999;
      if(el->gsfTrack().get())
	{
	  std::pair<bool,Measurement1D> ip3dRes = getImpactParameter<reco::GsfTrackRef>(el->gsfTrack(), primVtxRef, iSetup, true);
	  ev_.l_ip3d    = ip3dRes.second.value();
	  ev_.l_ip3dsig = ip3dRes.second.significance();
	}
    }

  //require no other leptons in the event
  bool passVetoLepton(vetoElectrons.size()+vetoMuons.size()==0);
  bool passNonIsoVetoLepton(vetoNonIsoElectrons.size()+vetoNonIsoMuons.size()==0);
  if( !passVetoLepton && !passNonIsoVetoLepton) return;
  
  // JETS
  ev_.nj=0;
  ev_.ngen=0;
  ev_.npf=0;
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetToken_,jets);
  edm::Handle<edm::ValueMap<float> > qgHandle; 
  iEvent.getByToken(qgToken_, qgHandle);  
  int ncentralJets(0);
  for(auto j = jets->begin();  j != jets->end(); ++j)
    {
      //kinematics
      if(j->pt()<20 || fabs(j->eta())>4.7) continue;

      //cross clean to selected tight leptons
      float dR2lepton=9999.;
      for(size_t il=0; il<selectedMuons.size(); il++)     dR2lepton = TMath::Min(dR2lepton,(float)deltaR(*j,*(selectedMuons[il])));
      for(size_t il=0; il<selectedElectrons.size(); il++) dR2lepton = TMath::Min(dR2lepton,(float)deltaR(*j,*(selectedElectrons[il])));
      if(dR2lepton<0.4) continue;
      
      // PF jet ID
      pat::strbitset retpf = pfjetIDLoose_.getBitTemplate();
      retpf.set(false);
      bool passLoose=pfjetIDLoose_( *j, retpf );
      if(!passLoose) continue;

      bool isCentral(false);
      if(j->pt()>25 && fabs(j->eta())<2.5) isCentral=true;
      ncentralJets+=isCentral;      

      //save jet
      const reco::Candidate *genParton = j->genParton();
      const reco::GenJet *genJet=j->genJet(); 
      ev_.j_area[ev_.nj]=j->jetArea();
      ev_.j_pt[ev_.nj]=j->pt();
      ev_.j_mass[ev_.nj]=j->mass();
      ev_.j_eta[ev_.nj]=j->eta();
      ev_.j_phi[ev_.nj]=j->phi();
      ev_.genj_pt[ev_.nj]=genJet ? genJet->pt() : 0;
      ev_.genj_mass[ev_.nj]=genJet ? genJet->mass() : 0;
      ev_.genj_eta[ev_.nj]=genJet ? genJet->eta() : 0;
      ev_.genj_phi[ev_.nj]=genJet ?  genJet->phi() : 0;
      ev_.j_csv[ev_.nj]=j->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      ev_.j_vtxpx[ev_.nj]=j->userFloat("vtxPx");
      ev_.j_vtxpy[ev_.nj]=j->userFloat("vtxPy");
      ev_.j_vtxpz[ev_.nj]=j->userFloat("vtxPz");
      ev_.j_vtxmass[ev_.nj]=j->userFloat("vtxMass");
      ev_.j_vtxNtracks[ev_.nj]=j->userFloat("vtxNtracks");
      ev_.j_vtx3DVal[ev_.nj]=j->userFloat("vtx3DVal");
      ev_.j_vtx3DSig[ev_.nj]=j->userFloat("vtx3DSig");
      ev_.j_puid[ev_.nj]=j->userFloat("pileupJetId:fullDiscriminant");
      edm::RefToBase<pat::Jet> jetRef(edm::Ref<edm::View<pat::Jet> >(jets, j - jets->begin()) );
      ev_.j_qg[ev_.nj] = (*qgHandle)[jetRef];
      ev_.j_flav[ev_.nj]=j->partonFlavour();
      ev_.j_hadflav[ev_.nj]=j->hadronFlavour();
      ev_.j_pid[ev_.nj]=genParton ? genParton->pdgId() : 0;
      
      //save gen/PF candidates for this jet
      if(isCentral)
	{
	  if(genJet)
	    {
	      for(size_t igen=0; igen<genJet->numberOfDaughters(); igen++)
		{
		  const reco::Candidate*gp=genJet->daughter(igen);
		  ev_.g_j[ev_.ngen]      = ev_.nj;
		  ev_.g_id[ev_.ngen]     = gp->pdgId();
		  ev_.g_charge[ev_.ngen] = gp->charge();
		  ev_.g_px[ev_.ngen]     = gp->px();
		  ev_.g_py[ev_.ngen]     = gp->py();
		  ev_.g_pz[ev_.ngen]     = gp->pz();
		  ev_.ngen++;
		}
	    }

	  for(size_t ipf=0; ipf<j->numberOfDaughters(); ipf++)
	    {
	      const reco::Candidate *pf=j->daughter(ipf);
	      ev_.pf_j[ev_.npf]      = ev_.nj;
	      ev_.pf_id[ev_.npf]     = pf->pdgId();
	      ev_.pf_charge[ev_.npf] = pf->charge();
	      ev_.pf_px[ev_.npf]     = pf->px();
	      ev_.pf_py[ev_.npf]     = pf->py();
	      ev_.pf_pz[ev_.npf]     = pf->pz();
	      ev_.npf++;
	    }
	}
      
      //store jet
      ev_.nj++;
    }
  if(ncentralJets==0) return;


  // MET
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  float metpt = mets->at(0).pt();
  float metphi = mets->at(0).phi();
  float dphi_met_lepton = deltaPhi(ev_.l_phi, metphi); 
  float mt=sqrt(2*ev_.l_pt*metpt*(1-cos(dphi_met_lepton)));
  ev_.met_pt=metpt;
  ev_.met_phi=metphi;
  ev_.mt=mt;

  //save tree if event is interesting
  if(saveTree_) tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzer::beginJob(){
}

//
void
MiniAnalyzer::endRun(const edm::Run & iRun, edm::EventSetup const & iSetup)
{
  try{
    
    edm::Handle<LHERunInfoProduct> lheruninfo;
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    iRun.getByLabel( "externalLHEProducer", lheruninfo );
    
    LHERunInfoProduct myLHERunInfoProduct = *(lheruninfo.product());
    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); 
	 iter!=myLHERunInfoProduct.headers_end(); 
	 iter++)
      {
	std::string tag("generator");
	if(iter->tag()!="") tag+="_"+iter->tag();
		
	std::vector<std::string> lines = iter->lines();
	std::vector<std::string> prunedLines;
	for (unsigned int iLine = 0; iLine<lines.size(); iLine++) 
	  {
	    if(lines.at(iLine)=="") continue;
	    if(lines.at(iLine).find("weightgroup")!=std::string::npos) continue;
	    prunedLines.push_back( lines.at(iLine) );
	  }
	
	if(histContainer_.find(tag)==histContainer_.end()) 
	  histContainer_[tag]=fs->make<TH1F>(tag.c_str(),tag.c_str(),prunedLines.size(),0,prunedLines.size());
	for (unsigned int iLine = 0; iLine<prunedLines.size(); iLine++) 
	  histContainer_[tag]->GetXaxis()->SetBinLabel(iLine+1,prunedLines.at(iLine).c_str());  
      }
  }
  catch(...){
    std::cout << "Failed to retrieve LHERunInfoProduct" << std::endl;
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAnalyzer::endJob() 
{
}

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
