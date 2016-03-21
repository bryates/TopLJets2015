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
  float getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
			 const reco::Candidate* ptcl,  
                         float r_iso_min, float r_iso_max, float kt_scale,
                         bool charged_only);

  // member data 
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron>  >  electronToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_, puppiMetToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorevtToken_;
  edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>  > genJetsToken_;
  edm::EDGetTokenT< edm::ValueMap<double> > puppiWeightToken_,puppiWeightNoLepToken_;

  //Electron Decisions
  edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleTightIdFullInfoMapToken_;

  std::unordered_map<std::string,TH1F*> histContainer_;

  PFJetIDSelectionFunctor pfjetIDLoose_;

  std::vector<std::string> muTriggersToUse_, elTriggersToUse_;

  
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
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  puppiMetToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("puppimets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  generatorevtToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator",""))),
  generatorlheToken_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
  puToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),
  genParticleToken_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"))),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJets"))),
  puppiWeightToken_(consumes< edm::ValueMap<double> >(edm::InputTag("PuppiWithLepton"))),
  puppiWeightNoLepToken_(consumes< edm::ValueMap<double> >(edm::InputTag("PuppiWithoutLepton"))),
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  eleTightIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleTightIdFullInfoMap"))),
  pfjetIDLoose_( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE ),
  
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
  histContainer_["counter"] = fs->make<TH1F>("counter", ";Counter;Events",2,0,2);
  histContainer_["pu"]      = fs->make<TH1F>("pu",      ";Pileup observed;Events",100,0,100);
  histContainer_["putrue"]  = fs->make<TH1F>("putrue",  ";Pileup true;Events",100,0,100);
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
    iEvent.getByToken(genParticleToken_,genParticles);
    
    //require only one lepton (can be from tau, if tau not from hadron)
    std::vector<std::pair<float,float> > letaphi;
    for (size_t i = 0; i < genParticles->size(); ++i) 
      {
	const GenParticle & genIt = (*genParticles)[i];
	if(!genIt.isPromptFinalState() && !genIt.isDirectPromptTauDecayProductFinalState()) continue;
	int ID = abs(genIt.pdgId());
	if(ID!=11 && ID!=13) continue;
	if(genIt.pt()<25 || fabs(genIt.eta())>2.5) continue;
	letaphi.push_back(std::pair<float,float>(genIt.eta(),genIt.phi()));
      }
    
    //require 1 jets not overlapping with leptons
    edm::Handle< std::vector<reco::GenJet> > genJets;
    iEvent.getByToken(genJetsToken_,genJets);
    std::map<const reco::Candidate *,int> jetConstsMap;
    for(std::vector<reco::GenJet>::const_iterator genJet=genJets->begin(); genJet!=genJets->end(); genJet++)
     {
       if(genJet->pt()<20) continue;

       //cross clean with leptons
       float minDR(99999.);
       for(size_t il=0; il<letaphi.size(); il++)
	 {
	   minDR=min((float)minDR,
		     (float)deltaR(genJet->eta(),genJet->phi(),letaphi[il].first,letaphi[il].second));
	 }
       if(minDR<0.4) continue;

       //map the gen particles which are clustered in this jet
       std::vector< const reco::Candidate * > jconst=genJet->getJetConstituentsQuick ();
       for(size_t ijc=0; ijc <jconst.size(); ijc++)
	 jetConstsMap[ jconst[ijc] ] = ev_.ngenj;
       
       ev_.genj_pt[ev_.ngenj]   = genJet->pt();
       ev_.genj_eta[ev_.ngenj]  = genJet->eta();
       ev_.genj_phi[ev_.ngenj]  = genJet->phi();
       ev_.genj_mass[ev_.ngenj] = genJet->mass();       
       ev_.ngenj++;
     }
    
    //store final state particles
    for (size_t i = 0; i < genParticles->size(); ++i)
      {
	const GenParticle & genIt = (*genParticles)[i];
	if(genIt.isHardProcess())
	  {
	    ev_.ghp_id[ ev_.ngenHardProc ]  = genIt.pdgId();
	    ev_.ghp_pt[ ev_.ngenHardProc ]  = genIt.pt();
	    ev_.ghp_eta[ ev_.ngenHardProc ] = genIt.eta();
	    ev_.ghp_phi[ ev_.ngenHardProc ] = genIt.phi();
	    ev_.ghp_m[ ev_.ngenHardProc ]   = genIt.mass();
	    ev_.ngenHardProc++;
	  }

	//final state particles
	if(genIt.status()==1)
	  {
	    ev_.g_id[ev_.ngen]     = genIt.pdgId();

	    //check if clustered in jet
	    ev_.g_genj[ev_.ngen]=-1;
	    for(std::map<const reco::Candidate *,int>::iterator it=jetConstsMap.begin();
		it!=jetConstsMap.end();
		it++)
	      {
		if(it->first->pdgId()==genIt.pdgId()) continue;
		if(deltaR( *(it->first), genIt)>0.01) continue; 
		ev_.g_genj[ev_.ngen]=it->second;
		break;
	      }	    
	    ev_.g_charge[ev_.ngen] = genIt.charge();
	    ev_.g_pt[ev_.ngen]     = genIt.pt();
	    ev_.g_eta[ev_.ngen]    = genIt.eta();
	    ev_.g_phi[ev_.ngen]    = genIt.phi();
	    ev_.ngen++;
	  }
      }

    return letaphi.size();
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
  ev_.isData  = iEvent.isRealData();
  
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
  ev_.isFiducial = 0;
  ev_.ttbar_nw=0;
  ev_.ngenj=0; 
  ev_.ngen=0;
  ev_.ngenHardProc=0;
  
  if(!ev_.isData)
    {
       Int_t nleptons = doFiducialAnalysis(iEvent,iSetup);
      ev_.isFiducial= (nleptons>0 ? ev_.ngenj : -1);
       
      edm::Handle<GenEventInfoProduct> evt;
      iEvent.getByToken( generatorToken_,evt);
      if(evt.isValid())
	{
	  ev_.ttbar_allmepartons   = evt->nMEPartons();
	  ev_.ttbar_matchmepartons = evt->nMEPartonsFiltered();
	  ev_.ttbar_w[0]           = evt->weight();
	  ev_.ttbar_nw++;
	}
      histContainer_["counter"]->Fill(1,ev_.ttbar_w[0]);

      edm::Handle<LHEEventProduct> evet;
      iEvent.getByToken(generatorlheToken_, evet);
      if(evet.isValid())
	{
	  double asdd=evet->originalXWGTUP();
	  for(unsigned int i=0  ; i<evet->weights().size();i++){
	    double asdde=evet->weights()[i].wgt;
	    ev_.ttbar_w[ev_.ttbar_nw]=ev_.ttbar_w[0]*asdde/asdd;
	    ev_.ttbar_nw++;
	  }
	}

      for(Int_t igenjet=0; igenjet<5; igenjet++)
	{
	  TString tag("fidcounter"); tag+=igenjet;
	  histContainer_[tag.Data()]->Fill(0.,ev_.ttbar_w[0]);
	  if(igenjet<=ev_.ngenj && nleptons>0)
	    {
	      for(Int_t iw=1; iw<ev_.ttbar_nw; iw++)
		histContainer_[tag.Data()]->Fill((float)iw,ev_.ttbar_w[iw]);
	    }
	}

      edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
      iEvent.getByToken(puToken_,PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator ipu;
      for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) 
	{
	  if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
	  ev_.pu=ipu->getPU_NumInteractions();
	  ev_.putrue=ipu->getTrueNumInteractions();
	  histContainer_["pu"]->Fill(ev_.pu);
	  histContainer_["putrue"]->Fill(ev_.putrue);
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
 
  //PF candidates
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfToken_,pfcands);
  edm::Handle< edm::ValueMap<double> > puppiWeight,puppiWeightNoLep;
  iEvent.getByToken( puppiWeightToken_, puppiWeight);
  iEvent.getByToken( puppiWeightNoLepToken_, puppiWeightNoLep);

  //
  //LEPTON SELECTION
  //
  int nLeadLeptons(0);
  ev_.nl=0;

  //MUON SELECTION: cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  for (const pat::Muon &mu : *muons) 
    { 
      //kinematics
      bool passPt( mu.pt() > 10 );
      bool passEta(fabs(mu.eta()) < 2.4 );
      if(!passPt || !passEta) continue;

      //ID
      bool isMedium(muon::isMediumMuon(mu));
      bool isTight(muon::isTightMuon(mu,primVtx));
      if(!isMedium) continue;

      //count
      nLeadLeptons+=( isTight && mu.pt()>25);

      //save it
      const reco::GenParticle * gen=mu.genLepton(); 
      ev_.isPromptFinalState[ev_.nl] = gen ? gen->isPromptFinalState() : false;
      ev_.isDirectPromptTauDecayProductFinalState[ev_.nl] = gen ? gen->isDirectPromptTauDecayProductFinalState() : false;
      ev_.l_id[ev_.nl]=13;
      ev_.l_charge[ev_.nl]=mu.charge();
      ev_.l_pt[ev_.nl]=mu.pt();
      ev_.l_eta[ev_.nl]=mu.eta();
      ev_.l_phi[ev_.nl]=mu.phi();
      ev_.l_mass[ev_.nl]=mu.mass();
      ev_.l_pid[ev_.nl]=(isMedium | (isTight<<1));
      ev_.l_chargedHadronIso[ev_.nl]=mu.pfIsolationR04().sumChargedHadronPt;
      ev_.l_miniIso[ev_.nl]=getMiniIsolation(pfcands,&mu,0.05,0.2, 10., false);
      ev_.l_relIso[ev_.nl]=(mu.pfIsolationR04().sumChargedHadronPt + max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt))/mu.pt();
      ev_.l_ip3d[ev_.nl] = -9999.;
      ev_.l_ip3dsig[ev_.nl] = -9999;
      if(mu.innerTrack().get())
	{
	  std::pair<bool,Measurement1D> ip3dRes = getImpactParameter<reco::TrackRef>(mu.innerTrack(), primVtxRef, iSetup, true);
	  ev_.l_ip3d[ev_.nl]    = ip3dRes.second.value();
	  ev_.l_ip3dsig[ev_.nl] = ip3dRes.second.significance();
	}
      ev_.nl++;     
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
  Int_t nele(0);
  for (const pat::Electron &el : *electrons) 
    {        
      const auto e = electrons->ptrAt(nele); 

      //kinematics cuts
      bool passPt(e->pt() > 15.0);
      bool passEta(fabs(e->eta()) < 2.5 && (fabs(e->superCluster()->eta()) < 1.4442 || fabs(e->superCluster()->eta()) > 1.5660));
      if(!passPt || !passEta) continue;
     
      //look up id decisions
      bool passVetoId = (*veto_id_decisions)[e];
      bool passTightId  = (*tight_id_decisions)[e];
      vid::CutFlowResult fullCutFlowData = (*tight_id_cutflow_data)[e];
      bool passTightIdExceptIso(true);
      for(size_t icut = 0; icut<fullCutFlowData.cutFlowSize(); icut++)
 	{
	  if(icut!=9 && !fullCutFlowData.getCutResultByIndex(icut)) passTightIdExceptIso=false;
	}

      //count
      nLeadLeptons += (passTightIdExceptIso && e->pt()>25);

      //save the electron
      const reco::GenParticle * gen=el.genLepton(); 
      ev_.isPromptFinalState[ev_.nl] = gen ? gen->isPromptFinalState() : false;
      ev_.isDirectPromptTauDecayProductFinalState[ev_.nl] = gen ? gen->isDirectPromptTauDecayProductFinalState() : false;
      ev_.l_id[ev_.nl]=11;
      ev_.l_pid[ev_.nl]=(passVetoId | (passTightId<<1) | (passTightIdExceptIso<<2));
      ev_.l_charge[ev_.nl]=e->charge();
      ev_.l_pt[ev_.nl]=e->pt();
      ev_.l_eta[ev_.nl]=e->eta();
      ev_.l_phi[ev_.nl]=e->phi();
      ev_.l_mass[ev_.nl]=e->mass();
      ev_.l_miniIso[ev_.nl]=getMiniIsolation(pfcands,&el,0.05, 0.2, 10., false);
      //ev_.l_relIso[ev_.nl]=fullCutFlowData.getValueCutUpon(9);
      ev_.l_relIso[ev_.nl]=(el.chargedHadronIso()+ max(0., el.neutralHadronIso() + el.photonIso()  - 0.5*el.puChargedHadronIso()))/el.pt();     
      ev_.l_chargedHadronIso[ev_.nl]=el.chargedHadronIso();
      ev_.l_ip3d[ev_.nl] = -9999.;
      ev_.l_ip3dsig[ev_.nl] = -9999;
      if(el.gsfTrack().get())
	{
	  std::pair<bool,Measurement1D> ip3dRes = getImpactParameter<reco::GsfTrackRef>(el.gsfTrack(), primVtxRef, iSetup, true);
	  ev_.l_ip3d[ev_.nl]    = ip3dRes.second.value();
	  ev_.l_ip3dsig[ev_.nl] = ip3dRes.second.significance();
	}
      ev_.nl++;
      if(passTightIdExceptIso && e->pt()>30) nele++;
    }

  //require at least 1 tight lepton in the event
  if(nLeadLeptons==0) return;

  // JETS
  ev_.nj=0;
  ev_.npf=0;
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetToken_,jets);
  std::vector<const reco::Candidate *>clustCands;
  for(auto j = jets->begin();  j != jets->end(); ++j)
    {
      //kinematics
      if(j->pt()<20 || fabs(j->eta())>4.7) continue;
      
      // PF jet ID
      pat::strbitset retpf = pfjetIDLoose_.getBitTemplate();
      retpf.set(false);
      bool passLoose=pfjetIDLoose_( *j, retpf );
      if(!passLoose) continue;
     
      //save jet
      const reco::Candidate *genParton = j->genParton();
      const reco::GenJet *genJet=j->genJet(); 
      ev_.j_area[ev_.nj]=j->jetArea();
      ev_.j_rawsf[ev_.nj]=j->correctedJet("Uncorrected").pt()/j->pt();
      ev_.j_pt[ev_.nj]=j->pt();
      ev_.j_mass[ev_.nj]=j->mass();
      ev_.j_eta[ev_.nj]=j->eta();
      ev_.j_phi[ev_.nj]=j->phi();

      Int_t matchGenJ(-1);
      if(genJet)
	{
	  float minDR(99999.);
	  for(Int_t igenj=0; igenj<ev_.ngenj; igenj++)   
	    {
	      float dR=deltaR(genJet->eta(),genJet->phi(),ev_.genj_eta[igenj],ev_.genj_phi[igenj]);
	      if(minDR<dR) continue;
	      minDR=dR;
	      matchGenJ=igenj;
	    }
	}
      ev_.j_ngj[ev_.nj]=matchGenJ;	
      ev_.j_csv[ev_.nj]=j->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      ev_.j_cvsl[ev_.nj]=j->bDiscriminator("pfCombinedCvsLJetTags");
      ev_.j_cvsb[ev_.nj]=j->bDiscriminator("pfCombinedCvsBJetTags");
      ev_.j_vtxpx[ev_.nj]=j->userFloat("vtxPx");
      ev_.j_vtxpy[ev_.nj]=j->userFloat("vtxPy");
      ev_.j_vtxpz[ev_.nj]=j->userFloat("vtxPz");
      ev_.j_vtxmass[ev_.nj]=j->userFloat("vtxMass");
      ev_.j_vtxNtracks[ev_.nj]=j->userFloat("vtxNtracks");
      ev_.j_vtx3DVal[ev_.nj]=j->userFloat("vtx3DVal");
      ev_.j_vtx3DSig[ev_.nj]=j->userFloat("vtx3DSig");
      ev_.j_puid[ev_.nj]=j->userFloat("pileupJetId:fullDiscriminant");      
      ev_.j_flav[ev_.nj]=j->partonFlavour();
      ev_.j_hadflav[ev_.nj]=j->hadronFlavour();
      ev_.j_pid[ev_.nj]=genParton ? genParton->pdgId() : 0;
      ev_.nj++;

      //save all PF candidates central jet
      if(fabs(j->eta())>2.5) continue;
      for(size_t ipf=0; ipf<j->numberOfDaughters(); ipf++)
	{
	  const reco::Candidate *pf=j->daughter(ipf);
	  clustCands.push_back(pf);
	  if(ev_.npf>=5000) continue;
	  ev_.pf_j[ev_.npf]           = ev_.nj-1;
	  ev_.pf_id[ev_.npf]         = pf->pdgId();
	  ev_.pf_charge[ev_.npf] = pf->charge();
	  ev_.pf_pt[ev_.npf]         = pf->pt();
	  ev_.pf_eta[ev_.npf]       = pf->eta();
	  ev_.pf_phi[ev_.npf]       = pf->phi();
	  //ev_.pf_puppiWgt[ev_.npf]=(*puppiWeight)[ *pf ];
	  // ev_.pf_puppiWgtNoLep[ev_.npf]=(*puppiWeightNoLep)[pf];	  
	  ev_.npf++;
	}
    }
      
  // MET
  ev_.nmet=2;
  for(int i=0; i<2; i++)
    {
        edm::Handle<pat::METCollection> mets;
	if(i==0) iEvent.getByToken(metToken_, mets);
	if(i==1) iEvent.getByToken(puppiMetToken_, mets);
	ev_.met_pt[i]=mets->at(0).pt();
	ev_.met_phi[i]=mets->at(0).phi();
    }

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

//-------------
//cf. https://twiki.cern.ch/twiki/bin/view/CMS/MiniIsolationSUSY
float MiniAnalyzer::getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
				     const reco::Candidate* ptcl,  
				     float r_iso_min, float r_iso_max, float kt_scale,
				     bool charged_only) 
{

    if (ptcl->pt()<5.) return 99999.;

    float deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
      if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  
    } else {
      //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    float iso_nh(0.), iso_ch(0.), iso_ph(0.), iso_pu(0.);
    float ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    float r_iso = (float)TMath::Max((float)r_iso_min,
				    (float)TMath::Min((float)r_iso_max, (float)(kt_scale/ptcl->pt())));
    for (const pat::PackedCandidate &pfc : *pfcands) {
      if (abs(pfc.pdgId())<7) continue;
      
      float dr = deltaR(pfc, *ptcl);
      if (dr > r_iso) continue;
      
      //////////////////  NEUTRALS  /////////////////////////
      if (pfc.charge()==0){
        if (pfc.pt()>ptThresh) {
          /////////// PHOTONS ////////////
          if (abs(pfc.pdgId())==22) {
            if(dr < deadcone_ph) continue;
            iso_ph += pfc.pt();
	    /////////// NEUTRAL HADRONS ////////////
          } else if (abs(pfc.pdgId())==130) {
            if(dr < deadcone_nh) continue;
            iso_nh += pfc.pt();
          }
        }
        //////////////////  CHARGED from PV  /////////////////////////
      } else if (pfc.fromPV()>1){
        if (abs(pfc.pdgId())==211) {
          if(dr < deadcone_ch) continue;
          iso_ch += pfc.pt();
        }
        //////////////////  CHARGED from PU  /////////////////////////
      } else {
        if (pfc.pt()>ptThresh){
          if(dr < deadcone_pu) continue;
          iso_pu += pfc.pt();
        }
      }
    }
    float iso(0.);
    if (charged_only){
      iso = iso_ch;
    } else {
      iso = iso_ph + iso_nh;
      iso -= 0.5*iso_pu;
      if (iso>0) iso += iso_ch;
      else iso = iso_ch;
    }
    iso = iso/ptcl->pt();

    return iso;
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
