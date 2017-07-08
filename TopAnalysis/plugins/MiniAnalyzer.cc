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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

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

#include "TopLJets2015/TopAnalysis/interface/rochcor2016.h"

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

// For PF skimming
bool pfSort(std::pair<int,double> pf1, std::pair<int,double> pf2) { return pf1.second > pf2.second; }

class MiniAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAnalyzer(const edm::ParameterSet&);
  ~MiniAnalyzer();  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  virtual void endRun(const edm::Run&,const edm::EventSetup&);  
private:
  virtual void beginJob() override;
  void genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void recAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  float getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
			 const reco::Candidate* ptcl,  
                         float r_iso_min, float r_iso_max, float kt_scale,
                         bool charged_only);

  // member data 
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorevtToken_;
  edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;
  edm::EDGetTokenT<LHERunInfoProduct> generatorRunInfoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>  > genLeptonsToken_,   genJetsToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> pseudoTopToken_;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron>  >  electronToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_, puppiMetToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  
  //Electron Decisions
  edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleTightIdFullInfoMapToken_;

  std::unordered_map<std::string,TH1F*> histContainer_;

  //muon rochester corrections
  rochcor2016 *rochcor_;

  PFJetIDSelectionFunctor pfjetIDLoose_;

  std::vector<std::string> muTriggersToUse_, elTriggersToUse_;

  bool saveTree_,savePF_;
  TTree *tree_;
  MiniEvent_t ev_;

  //a counter for generator level scans (e.g. sms scans)
  std::map<TString, float>  genScanCounter_;

  edm::Service<TFileService> fs;
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
  generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  generatorevtToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator",""))),
  generatorlheToken_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
  generatorRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>({"externalLHEProducer"})),
  puToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),  
  genLeptonsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("pseudoTop:leptons"))),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("pseudoTop:jets"))),
  genParticlesToken_(consumes<pat::PackedGenParticleCollection>(edm::InputTag("packedGenParticles"))),
  prunedGenParticlesToken_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"))),
  pseudoTopToken_(consumes<reco::GenParticleCollection>(edm::InputTag("pseudoTop"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),  
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  puppiMetToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("puppimets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  eleTightIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleTightIdFullInfoMap"))),
  pfjetIDLoose_( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE ),  
  saveTree_( iConfig.getParameter<bool>("saveTree") ),
  savePF_( iConfig.getParameter<bool>("savePF") )
{
  //now do what ever initialization is needed
  electronToken_ = mayConsume<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  elTriggersToUse_ = iConfig.getParameter<std::vector<std::string> >("elTriggersToUse");
  muTriggersToUse_ = iConfig.getParameter<std::vector<std::string> >("muTriggersToUse");

  //start the rochester correction tool
  rochcor_=new rochcor2016(2016);
 
  //  usesResource("TFileService");

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
void MiniAnalyzer::genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //
  // PILEUP
  //
  edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
  iEvent.getByToken(puToken_,PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator ipu;
  for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) 
    {
      if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
      ev_.pu=ipu->getPU_NumInteractions();
      ev_.putrue=ipu->getTrueNumInteractions();
    }
  histContainer_["pu"]->Fill(ev_.pu);
  histContainer_["putrue"]->Fill(ev_.putrue);
  
  //
  // GENERATOR WEIGHTS
  //
  ev_.ttbar_nw=0;
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
  
  //alternative weights for systematics
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
      
  //
  // GENERATOR LEVEL EVENT
  //
  ev_.ng=0; ev_.ngjets=0; ev_.ngbjets=0;
  edm::Handle<std::vector<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsToken_,genJets);  
  std::map<const reco::Candidate *,int> jetConstsMap;
  for(auto genJet = genJets->begin();  genJet != genJets->end(); ++genJet)
    {
      //map the gen particles which are clustered in this jet
      std::vector< const reco::Candidate * > jconst=genJet->getJetConstituentsQuick ();
      for(size_t ijc=0; ijc <jconst.size(); ijc++)
	jetConstsMap[ jconst[ijc] ] = ev_.ng;
      
      ev_.g_id[ev_.ng]   = genJet->pdgId();
      ev_.g_pt[ev_.ng]   = genJet->pt();
      ev_.g_eta[ev_.ng]  = genJet->eta();
      ev_.g_phi[ev_.ng]  = genJet->phi();
      ev_.g_m[ev_.ng]    = genJet->mass();       
      ev_.ng++;
      
      //gen level selection
      if(genJet->pt()>25 && fabs(genJet->eta())<2.5)
	{
	  ev_.ngjets++;
	  if(abs(genJet->pdgId())) ev_.ngbjets++;
	}
    }

  //leptons
  ev_.ngleptons=-0;
  edm::Handle<std::vector<reco::GenJet> > dressedLeptons;
  iEvent.getByToken(genLeptonsToken_,dressedLeptons);
  for(auto genLep = dressedLeptons->begin();  genLep != dressedLeptons->end(); ++genLep)
    {
      //map the gen particles which are clustered in this lepton
      std::vector< const reco::Candidate * > jconst=genLep->getJetConstituentsQuick ();
      for(size_t ijc=0; ijc <jconst.size(); ijc++)
	jetConstsMap[ jconst[ijc] ] = ev_.ng;
      
      ev_.g_pt[ev_.ng]   = genLep->pt();
      ev_.g_id[ev_.ng]   = genLep->pdgId();
      ev_.g_eta[ev_.ng]  = genLep->eta();
      ev_.g_phi[ev_.ng]  = genLep->phi();
      ev_.g_m[ev_.ng]    = genLep->mass();       
      ev_.ng++;

      //gen level selection
      if(genLep->pt()>20 && fabs(genLep->eta())<2.5) ev_.ngleptons++;
    }
  
  
  //final state particles 
  ev_.ngpf=0;
  edm::Handle<pat::PackedGenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_,genParticles);
  for (size_t i = 0; i < genParticles->size(); ++i)
    {
      const pat::PackedGenParticle & genIt = (*genParticles)[i];

      //this shouldn't be needed according to the workbook
      //if(genIt.status()!=1) continue;
      if(genIt.pt()<0.5 || fabs(genIt.eta())>2.5) continue;
      
      ev_.gpf_id[ev_.ngpf]     = genIt.pdgId();
      ev_.gpf_c[ev_.ngpf]      = genIt.charge();
      ev_.gpf_g[ev_.ngpf]=-1;
      for(std::map<const reco::Candidate *,int>::iterator it=jetConstsMap.begin();
	  it!=jetConstsMap.end();
	  it++)
	{
	  if(it->first->pdgId()!=genIt.pdgId()) continue;
	  if(deltaR( *(it->first), genIt)>0.01) continue; 
	  ev_.gpf_g[ev_.ngpf]=it->second;
	  break;
	}
      ev_.gpf_pt[ev_.ngpf]     = genIt.pt();
      ev_.gpf_eta[ev_.ngpf]    = genIt.eta();
      ev_.gpf_phi[ev_.ngpf]    = genIt.phi();
      ev_.gpf_m[ev_.ngpf]      = genIt.mass();
      ev_.ngpf++;    
    }

  //J/Psi decay
  //Meson decay (J/Psi or D0)
  ev_.ngjpsi=0;
  ev_.ngmeson=0;
  //prurned also used for top/stop
  edm::Handle<reco::GenParticleCollection> prunedGenParticles;
  iEvent.getByToken(prunedGenParticlesToken_,prunedGenParticles);
  for(size_t i = 0; i < prunedGenParticles->size(); i++) {
    const reco::GenParticle &genIt = (*prunedGenParticles)[i];
    int absid=abs(genIt.pdgId());
    if(absid!=443 && absid!=421 && absid!=413) continue;
    if(genIt.numberOfDaughters()<2) continue;
    //if(genIt.daughter(0)->pdgId()*genIt.daughter(1)->pdgId()!=-13*13 &&
    //   genIt.daughter(0)->pdgId()*genIt.daughter(1)->pdgId()!=-211*321 &&
    //   genIt.daughter(0)->pdgId()*genIt.daughter(1)->pdgId()!=-211*321*-211) continue;
    //cout << "New meson found: ngmeson = " << ev_.ngmeson << endl; 
    //cout << "daughter found" << endl;
    bool JPsiDaughter = false;
    bool D0Daughter = false;
    bool DsDaughter = false;
    ev_.ngmeson_daug=0;
    for(size_t ipf=0; ipf<genIt.numberOfDaughters(); ipf++) {
      //cout << "loading daughter" << endl;
      const reco::Candidate *daug=genIt.daughter(ipf);
      if(abs(daug->pdgId())==13 && absid==443) JPsiDaughter = true;
      else if(abs(daug->pdgId())==211 && absid==421) D0Daughter = true;
      else if(abs(daug->pdgId())==321 && absid==421) D0Daughter = true;
      else if(abs(daug->pdgId())==211 && absid==413) DsDaughter = true;
      else if(abs(daug->pdgId())==321 && absid==413) DsDaughter = true;
      else continue;
      cout << "event = " << ev_.event << " ngmeson = " << ev_.ngmeson << " : pdgId() = " << genIt.pdgId() << " : daughter n = " << ipf <<  " : ngmeson_daughter = " << ev_.ngmeson_daug << " : daug->pdgId() = " << daug->pdgId() << endl;
      /*
      if(JPsiDaughter) cout << "J/Psi" << endl;
      else if(D0Daughter) cout << "D0" << endl;
      cout << "daugther is muon: pT, eta, phi" << endl;
      cout << daug->pdgId() << endl;
      cout << daug->pt() << ", ";
      cout << daug->eta() << ", ";
      cout << daug->phi() << endl;
      */
      //2*ev_.ngjpsi to account for 2 duaghters per J/Psi
      ev_.gmeson_daug_id[ev_.ngmeson_daug] = daug->pdgId();
      ev_.gmeson_daug_pt[ev_.ngmeson_daug] = daug->pt();
      ev_.gmeson_daug_eta[ev_.ngmeson_daug] = daug->eta();
      ev_.gmeson_daug_phi[ev_.ngmeson_daug] = daug->phi();
      ev_.gmeson_daug_meson_index[ev_.ngmeson_daug] = ev_.ngmeson;
      //cout << "i: " << i << " " << ev_.event << " " << ev_.ngmeson << " " << ev_.ngmeson_daug << " " << ev_.ngmeson << endl;

      /*
      ev_.gjpsi_mu_dxy[ev_.ngjpsi]  = daug->dxy(primVtx.position());
      ev_.gjpsi_mu_dxyE[ev_.ngjpsi]  = daug->dxyE();
      ev_.gjpsi_mu_dz[ev_.ngjpsi]  = daug->dz(primVtx.position());
      ev_.gjpsi_mu_dzE[ev_.ngjpsi]  = daug->dzE();
      */

      //Find t(tbar) mother
      while(abs(daug->pdgId()) != 6) {
        if(daug->mother() == 0) break;
        //int charge = daug->charge();
        //if(!charge) charge = 1;
        //cout << "PdgId= " << daug->pdgId()*charge << endl;
        daug = daug->mother();
        //charge = daug->charge();
        //if(!charge) charge = 1;
        //cout << "Mother PdgId= " << daug->pdgId() << endl;
        if(abs(daug->pdgId()) == 2212) break;
        if(abs(daug->pdgId()) == 22) break;
      }
      ev_.gmeson_mother_id[ev_.ngmeson_daug] = daug->pdgId(); //*daug->charge();
      ev_.ngmeson_daug++;
    }
    if(!JPsiDaughter && !D0Daughter && !DsDaughter) continue;
    /*
    cout << "J/Psi: pT, eta, phi, M" << endl;
    cout << genIt.pt() << ", ";
    cout << genIt.eta() << ", ";
    cout << genIt.phi() << ", ";
    cout << genIt.mass() << endl;
    */
    ev_.gmeson_id[ev_.ngmeson]     = genIt.pdgId();
    ev_.gmeson_pt[ev_.ngmeson]     = genIt.pt();
    ev_.gmeson_eta[ev_.ngmeson]    = genIt.eta();
    ev_.gmeson_phi[ev_.ngmeson]    = genIt.phi();
    ev_.gmeson_m[ev_.ngmeson]      = genIt.mass();
    ev_.gmeson_daug_dR[ev_.ngmeson]  = reco::deltaR(genIt.daughter(0)->eta(),genIt.daughter(0)->phi(),genIt.daughter(1)->eta(),genIt.daughter(1)->phi());
    ev_.gmeson_index[ev_.ngmeson] = ev_.ngmeson;
    if(JPsiDaughter) ev_.ngjpsi++;
    //if(JPsiDaughter || D0Daughter) ev_.ngmeson++;
    ev_.ngmeson++;
  }

  //top or stop quarks (lastCopy)
  //edm::Handle<reco::GenParticleCollection> prunedGenParticles;
  //iEvent.getByToken(prunedGenParticlesToken_,prunedGenParticles);
  ev_.ngtop=0; 
  float mStop(-1),mNeutralino(-1);
  for (size_t i = 0; i < prunedGenParticles->size(); ++i)
    {
      const reco::GenParticle & genIt = (*prunedGenParticles)[i];
      int absid=abs(genIt.pdgId());
      if(absid!=6 && absid!=1000006 && absid!=1000022) continue;
      if(!genIt.isLastCopy()) continue;      

      if(absid==1000006) mStop=genIt.mass();
      if(absid==1000022) mNeutralino=genIt.mass();

      ev_.gtop_id[ ev_.ngtop ]  = genIt.pdgId();
      ev_.gtop_pt[ ev_.ngtop ]  = genIt.pt();
      ev_.gtop_eta[ ev_.ngtop ] = genIt.eta();
      ev_.gtop_phi[ ev_.ngtop ] = genIt.phi();
      ev_.gtop_m[ ev_.ngtop ]   = genIt.mass();
      ev_.ngtop++;
    }

  //check if this is a SMS scan
  if(mStop>0 && mNeutralino>0)
    {
      TString key(Form("mstop_%d_mchi0_%d",int(10*mStop),int(10*mNeutralino)));
      if(genScanCounter_.find(key)==genScanCounter_.end()) genScanCounter_[key]=0.;
      genScanCounter_[key]+=ev_.ttbar_w[0];
    }
    

  //pseudo-tops 
  edm::Handle<reco::GenParticleCollection> pseudoTop;
  iEvent.getByToken(pseudoTopToken_,pseudoTop);
  for (size_t i = 0; i < pseudoTop->size(); ++i)
    {
      const GenParticle & genIt = (*pseudoTop)[i];
      ev_.gtop_id[ ev_.ngtop ]  = genIt.pdgId()*1000;
      ev_.gtop_pt[ ev_.ngtop ]  = genIt.pt();
      ev_.gtop_eta[ ev_.ngtop ] = genIt.eta();
      ev_.gtop_phi[ ev_.ngtop ] = genIt.phi();
      ev_.gtop_m[ ev_.ngtop ]   = genIt.mass();
      ev_.ngtop++;
    }
  
  //fiducial counters
  for(Int_t igenjet=0; igenjet<5; igenjet++)
    {
      TString tag("fidcounter"); tag+=igenjet;
      histContainer_[tag.Data()]->Fill(0.,ev_.ttbar_w[0]);
      if(igenjet<=ev_.ngjets && ev_.ngleptons>0)
	{
	  for(Int_t iw=1; iw<ev_.ttbar_nw; iw++)
	    histContainer_[tag.Data()]->Fill((float)iw,ev_.ttbar_w[iw]);
	}
    }
}


//
void MiniAnalyzer::recAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
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

  //
  //LEPTON SELECTION
  //
  ev_.nl=0; ev_.nleptons=0;

  //MUON SELECTION: cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  for (const pat::Muon &mu : *muons) 
    { 
      //correct the 4-momentum
      TLorentzVector p4;
      p4.SetPtEtaPhiM(mu.pt(),mu.eta(),mu.phi(),mu.mass());
      float qter(1.0);
      try {
        if(iEvent.isRealData()) {
	  rochcor_->momcor_data(p4, mu.charge(), 0, qter);
        }
        else {
          int ntrk=mu.innerTrack()->hitPattern().trackerLayersWithMeasurement();
	  rochcor_->momcor_data(p4, mu.charge(), ntrk, qter);
        }
      }
      catch(...) {
        //probably no inner track...
      }

      //kinematics
      bool passPt( mu.pt() > 10 );
      bool passEta(fabs(mu.eta()) < 2.4 );
      if(!passPt || !passEta) continue;

      //ID
      bool isMedium(muon::isMediumMuon(mu));
      bool isTight(muon::isTightMuon(mu,primVtx));
      bool isLoose(muon::isLooseMuon(mu));
      bool isSoft(muon::isSoftMuon(mu,primVtx));
      if(!isMedium) continue;

      //save it
      const reco::GenParticle * gen=mu.genLepton(); 
      ev_.isPromptFinalState[ev_.nl] = gen ? gen->isPromptFinalState() : false;
      ev_.isDirectPromptTauDecayProductFinalState[ev_.nl] = gen ? gen->isDirectPromptTauDecayProductFinalState() : false;
      ev_.l_id[ev_.nl]=13;
      ev_.l_g[ev_.nl]=-1;
      for(int ig=0; ig<ev_.ng; ig++)
	{
	  if(abs(ev_.g_id[ig])!=ev_.l_id[ev_.nl]) continue;
	  if(deltaR( mu.eta(),mu.phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
	  ev_.l_g[ev_.nl]=ig;
	  break;
	}	 
      ev_.l_charge[ev_.nl]=mu.charge();
      ev_.l_pt[ev_.nl]=mu.pt();
      ev_.l_eta[ev_.nl]=mu.eta();
      ev_.l_phi[ev_.nl]=mu.phi();
      ev_.l_mass[ev_.nl]=mu.mass();
      ev_.l_pid[ev_.nl]=(isTight | (isMedium<<1) | (isLoose<<2) | (isSoft<<3));
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
      if(mu.outerTrack().isNonnull()) {
        ev_.l_dxy[ev_.nl]=mu.dB();
        ev_.l_dxyE[ev_.nl]=mu.edB();
      }
      if(mu.innerTrack().isNonnull()) {
        ev_.l_dxy[ev_.nl]=mu.dB();
        ev_.l_dxyE[ev_.nl]=mu.edB();
        ev_.l_dz[ev_.nl]=mu.innerTrack()->dz(primVtx.position());
        ev_.l_global[ev_.nl]=mu.isGlobalMuon();
        ev_.l_pf[ev_.nl]=mu.isPFMuon();
        ev_.l_nValTrackerHits[ev_.nl] = mu.innerTrack()->hitPattern().numberOfValidTrackerHits();
        ev_.l_nValPixelHits[ev_.nl] = mu.innerTrack()->hitPattern().numberOfValidPixelHits();
        ev_.l_nMatchedStations[ev_.nl] = mu.numberOfMatchedStations();
        ev_.l_pixelLayerWithMeasurement[ev_.nl]    = mu.innerTrack()->hitPattern().pixelLayersWithMeasurement();
        ev_.l_trackerLayersWithMeasurement[ev_.nl] = mu.innerTrack()->hitPattern().trackerLayersWithMeasurement();

        //Medium
        ev_.l_validFraction[ev_.nl] = mu.innerTrack()->validFraction();
        ev_.l_chi2LocalPosition[ev_.nl] = mu.combinedQuality().chi2LocalPosition;
        ev_.l_trkKink[ev_.nl] = mu.combinedQuality().trkKink;
      }

      if (mu.globalTrack().isNonnull()) {
        ev_.l_chi2norm[ev_.nl]=mu.normChi2();
        ev_.l_global[ev_.nl]=mu.isGlobalMuon();
        ev_.l_pf[ev_.nl]=mu.isPFMuon();
        ev_.l_globalTrackNumberOfValidHits[ev_.nl] = mu.globalTrack()->hitPattern().numberOfValidMuonHits();
        ev_.l_nMatchedStations[ev_.nl] = mu.numberOfMatchedStations();
      }
      ev_.nl++;    
      ev_.nleptons += ( isTight && mu.pt()>25); 
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
      nele++;

      //kinematics cuts
      bool passPt(el.pt() > 15.0);
      bool passEta(fabs(el.eta()) < 2.5 && (fabs(el.superCluster()->eta()) < 1.4442 || fabs(el.superCluster()->eta()) > 1.5660));
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

      //save the electron
      const reco::GenParticle * gen=el.genLepton(); 
      ev_.isPromptFinalState[ev_.nl] = gen ? gen->isPromptFinalState() : false;
      ev_.isDirectPromptTauDecayProductFinalState[ev_.nl] = gen ? gen->isDirectPromptTauDecayProductFinalState() : false;
      ev_.l_id[ev_.nl]=11;
      ev_.l_g[ev_.nl]=-1;
      for(int ig=0; ig<ev_.ng; ig++)
	{
	  if(abs(ev_.g_id[ig])!=ev_.l_id[ev_.nl]) continue;
	  if(deltaR( el.eta(),el.phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
	  ev_.l_g[ev_.nl]=ig;
	  break;
	}	 
      ev_.l_pid[ev_.nl]=(passVetoId | (passTightId<<1) | (passTightIdExceptIso<<2));
      ev_.l_charge[ev_.nl]=el.charge();
      ev_.l_pt[ev_.nl]=el.pt();
      ev_.l_eta[ev_.nl]=el.eta();
      ev_.l_phi[ev_.nl]=el.phi();
      ev_.l_mass[ev_.nl]=el.mass();
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
      ev_.nleptons += (passTightIdExceptIso && el.pt()>25);
    }

  // JETS
  ev_.nj=0; 
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetToken_,jets);
  std::vector< std::pair<const reco::Candidate *,int> > clustCands;
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
      ev_.j_area[ev_.nj]=j->jetArea();
      ev_.j_rawsf[ev_.nj]=j->correctedJet("Uncorrected").pt()/j->pt();
      ev_.j_pt[ev_.nj]=j->pt();
      ev_.j_mass[ev_.nj]=j->mass();
      ev_.j_eta[ev_.nj]=j->eta();
      ev_.j_phi[ev_.nj]=j->phi();
      ev_.j_g[ev_.nj]=-1;
      for(int ig=0; ig<ev_.ng; ig++)
	{
	  if(abs(ev_.g_id[ig])==11 || abs(ev_.g_id[ig])==13) continue;
	  if(deltaR( j->eta(),j->phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
	  ev_.j_g[ev_.nj]=ig;
	  break;
	}	 
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

      ev_.j_pt_charged[ev_.nj]=0; //pT of charged tracks only
      //save all PF candidates central jet
      if(fabs(j->eta())>2.5) continue;
      ev_.j_pt_pf[ev_.nj]=0;
      for(size_t ipf=0; ipf<j->numberOfDaughters(); ipf++)
	{
	  const reco::Candidate *pf=j->daughter(ipf);
	  clustCands.push_back(std::pair<const reco::Candidate *,int>(pf,ev_.nj));
          ev_.j_pt_charged[ev_.nj] += abs(pf->pt()*pf->charge()); //pT of charged tracks only (neutral particles weighted by 0)
          ev_.j_pt_pf[ev_.nj] += pf->pt(); //pT of charged tracks only (neutral particles weighted by 0)
	}

      ev_.nj++;
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

  //PF candidates
  //Skim
  //Save all muons an 6 hardest others
  //This block only counts the number of PF types
  ev_.npf=0;
  std::vector<pair<int,double>> pfCand;
  std::map<int,pair<int,int>> nPFJet; //pf_j, nMu, nKPi;
  std::pair<int,double> hard = std::pair<int,double>(-1,0); //save hadest PF
  for(auto pf = pfcands->begin();  pf != pfcands->end(); ++pf) {
    if(ev_.npf>=5000) continue;

    ev_.pf_j[ev_.npf] = -1;
    for(size_t i=0; i<clustCands.size(); i++) {
      if(pf->pdgId()!=clustCands[i].first->pdgId()) continue;
      if(deltaR(*pf,*(clustCands[i].first))>0.01) continue;
      ev_.pf_j[ev_.npf]=clustCands[i].second;
      break;
    }

    //extra requirements for unclustered PF candidates
    if(ev_.pf_j[ev_.npf]==-1) continue;

    //Only save hardest if not already saved
    //e.g. not a muon or not in the top 6 cands
    if(ev_.pf_pt[ev_.npf] > hard.second)
      hard = std::pair<int,double>(ev_.npf,ev_.pf_pt[ev_.npf]);

    std::pair<int,int> mupik = std::pair<int,int>(0,0);
    if(fabs(pf->pdgId())==13) mupik.first++;
    else if(fabs(pf->pdgId())==211) mupik.second++;
    else continue;

    //If this jet is already in the map, increment with the new paticles
    if(nPFJet.find(ev_.pf_j[ev_.npf]) != nPFJet.end()) {
      nPFJet[ev_.pf_j[ev_.npf]].first += mupik.first;
      nPFJet[ev_.pf_j[ev_.npf]].second += mupik.second;
      /*
      cout << "jet " << ev_.pf_j[ev_.npf] << endl;
      cout << "mu " << nPFJet[ev_.pf_j[ev_.npf]].first << endl;
      cout << "pi/K " << nPFJet[ev_.pf_j[ev_.npf]].second << endl;
      */
    }
    //Otherwise, add this jet to the map
    else nPFJet[ev_.pf_j[ev_.npf]] = mupik;
    pfCand.push_back(std::pair<int,double>(ev_.npf,pf->pt()));
    ev_.pf_jnpf[ev_.pf_j[ev_.npf]]++;
    if(pf->trackHighPurity()) ev_.pf_jnhppf[ev_.pf_j[ev_.npf]]++;
    ev_.npf++;

  }

  sort(pfCand.begin(),pfCand.end(), [](std::pair<int,double> i, std::pair<int,double> j) { return i.second > j.second ; } );
  pfCand.resize(6); //keep only 6 hardest non muon PF candidates

  //Now actually put the desired jets in the ntuple
  ev_.npf=0;
  for(auto pf = pfcands->begin();  pf != pfcands->end(); ++pf)
    {
      if(ev_.npf>=5000) continue;
      //int npf = ev_.npf;
      //if(!(fabs(pf->pdgId())==13 || std::any_of(pfCand.begin(), pfCand.end(),
      //                              [npf](std::pair<int,double>& elem) {return elem.first == npf;} ))) continue;

      ev_.pf_j[ev_.npf] = -1;
      for(size_t i=0; i<clustCands.size(); i++)
	{
	  if(pf->pdgId()!=clustCands[i].first->pdgId()) continue;
	  if(deltaR(*pf,*(clustCands[i].first))>0.01) continue;
	  ev_.pf_j[ev_.npf]=clustCands[i].second;
	  break;
	}
      if(ev_.pf_j[ev_.npf]==-1) continue;

      //only save jets with 2+ mu (+/- 13) or 2+ K/pi (+/- 211) or hardest pT track (if not a mu or K)
      if(ev_.pf_id[ev_.npf])
      if(ev_.pf_j[ev_.npf]>=0 &&
        (nPFJet[ev_.pf_j[ev_.npf]].first < 2 &&
         nPFJet[ev_.pf_j[ev_.npf]].second < 2) &&
         ev_.npf != hard.first) continue;
      
      //extra requirements for unclustered PF candidates
      //ONLY save mu/pi/K
      //if(fabs(ev_.pf_id[ev_.npf])!=13 && fabs(ev_.pf_id[ev_.npf])!=211) continue;
      if(ev_.pf_j[ev_.npf]==-1)
	{
	  if(pf->charge()==0) continue;
	  if(pf->fromPV()<2 && fabs(pf->pdgId())!=13) continue;
	  if(pf->pt()<0.5 || fabs(pf->eta())>2.5) continue;
	  if(pf->puppiWeight()<0.01) continue; // && fabs(pf->pdgId())!=13) continue;
	}
      if(ev_.pf_j[ev_.npf]==-1) continue; // skip unclustered PF candidates
      
      ev_.pf_fromPV[ev_.npf]   = pf->fromPV();
      ev_.pf_id[ev_.npf]       = pf->pdgId();
      ev_.pf_c[ev_.npf]        = pf->charge();
      ev_.pf_pt[ev_.npf]       = pf->pt();
      ev_.pf_eta[ev_.npf]      = pf->eta();
      ev_.pf_phi[ev_.npf]      = pf->phi();
      ev_.pf_m[ev_.npf]        = pf->mass();
      ev_.pf_puppiWgt[ev_.npf] = pf->puppiWeight();      

      ev_.pf_dxy[ev_.npf]      = pf->dxy(primVtx.position());
      ev_.pf_dxyE[ev_.npf]     = pf->dxyError();
      ev_.pf_dz[ev_.npf]       = pf->dz(primVtx.position());
      ev_.pf_dzE[ev_.npf]      = pf->dzError();

      ev_.pf_highPurity[ev_.npf] = pf->trackHighPurity();
      ev_.pf_quality[ev_.npf] = -1;
      ev_.pf_quality[ev_.npf] = pf->pseudoTrack().qualityMask();
      ev_.pf_muon[ev_.npf] = pf->isMuon();
      ev_.pf_standAloneMuon[ev_.npf] = pf->isStandAloneMuon();
      ev_.pf_globalMuon[ev_.npf] = pf->isGlobalMuon();
      ev_.pf_trackerMuon[ev_.npf] = pf->isTrackerMuon();

      ev_.npf++;
    }
}

// ------------ method called for each event  ------------
void MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  histContainer_["counter"]->Fill(0);

  //analyze the event
  if(!iEvent.isRealData()) genAnalysis(iEvent,iSetup);
  recAnalysis(iEvent,iSetup);
  
  //save event if at least one lepton at gen or reco level
  if((ev_.ngleptons==0 && ev_.nleptons==0) || !saveTree_) return;  
  ev_.run     = iEvent.id().run();
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event(); 
  ev_.isData  = iEvent.isRealData();
  if(!savePF_) { ev_.ngpf=0; ev_.npf=0; }
  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzer::beginJob(){
}

//
void 
MiniAnalyzer::endRun(const edm::Run& iRun,
		     const EventSetup& iSetup) 
{
  try{

    cout << "[MiniAnalyzer::endRun]" << endl;

    //save histograms with the counts per point in the gen scan (if available)
    for(auto it=genScanCounter_.begin(); it!=genScanCounter_.end(); it++)
      {
	TString key(it->first);
	float counts(it->second);
	histContainer_[key.Data()]=fs->make<TH1F>(key,key,1,0,1);
	histContainer_[key.Data()]->SetBinContent(1,counts);
      }

    edm::Handle<LHERunInfoProduct> lheruninfo;
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    iRun.getByToken(generatorRunInfoToken_, lheruninfo );
    //iRun.getByLabel( "externalLHEProducer", lheruninfo );

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
	  {
	    std::cout << "Starting histo for " << tag << std::endl;
	    histContainer_[tag]=fs->make<TH1F>(tag.c_str(),tag.c_str(),prunedLines.size(),0,prunedLines.size());
	  }
	for (unsigned int iLine = 0; iLine<prunedLines.size(); iLine++) 
	  histContainer_[tag]->GetXaxis()->SetBinLabel(iLine+1,prunedLines.at(iLine).c_str());  
      }
  }
  catch(std::exception &e){
    std::cout << e.what() << endl
	      << "Failed to retrieve LHERunInfoProduct" << std::endl;
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
  std::cout << "[MiniAnalyzer::endJob]" << endl;
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
