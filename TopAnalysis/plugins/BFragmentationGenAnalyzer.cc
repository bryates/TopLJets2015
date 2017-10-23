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

#include <memory>
#include <string>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TTree.h"

#define IS_BHADRON_PDGID(id) ( ((abs(id)/100%10) == 5) || (abs(id) >= 5000 && abs(id) <= 5999) )
#define IS_NEUTRINO_PDGID(id) ( (abs(id) == 12) || (abs(id) == 14) || (abs(id) == 16) )



class FragmentationAnalyzer : public edm::EDAnalyzer {

 public:
  FragmentationAnalyzer(const edm::ParameterSet &);
  virtual void analyze(const edm::Event &, const edm::EventSetup &);

 private:
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;
  std::map<std::string, TH1F*> hists;
  TTree *data_;
  Int_t nB_, nL_, nJPsi_, Bid_[100], Lid_[100];
  Float_t Bpt_[100], Beta_[100], Bphi_[100], Bm_[100];
  Float_t Bjpt_[100], Bjeta_[100], Bjphi_[100], Bjm_[100];
  Float_t JPsipt_[100], JPsieta_[100], JPsiphi_[100], JPsim_[100];
  Float_t Lpt_[100], Leta_[100], Lphi_[100], Lm_[100];

};
 
FragmentationAnalyzer::FragmentationAnalyzer(const edm::ParameterSet &cfg) :
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("pseudoTop:jets"))),
  prunedGenParticlesToken_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles")))
{
  //Load TFile Service
  edm::Service<TFileService> fs;
  if(!fs)
    throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  hists["genBHadronJPsiDecay"] = fs->make<TH1F>("genBHadronJPsiDecay", "genBHadronJPsiDecay", 2, 0, 2);
  hists["genBHadronNuDecay"] = fs->make<TH1F>("genBHadronDecay", "genBHadronDecay", 2, 0, 2);
  hists["genBHadronPtFraction"] = fs->make<TH1F>("genBHadronPtFraction", "genBHadronPtFraction", 2, 0, 2);
  for(auto & it : hists) it.second->Sumw2();

  //sumary tree for B, J/Psi, and leptons
  data_ = fs->make<TTree>("FragTree", "FragTree");
  data_->Branch("nB",    &nB_,    "nB/I");
  data_->Branch("Bid",    Bid_,   "Bid[nB]/I");
  data_->Branch("Bpt",    Bpt_,   "Bpt[nB]/F");
  data_->Branch("Beta",   Beta_,  "Beta[nB]/F");
  data_->Branch("Bphi",   Bphi_,  "Bphi[nB]/F");
  data_->Branch("Bm",     Bm_,    "Bm[nB]/F");
  data_->Branch("Bjpt",   Bjpt_,  "Bjpt[nB]/F");
  data_->Branch("Bjeta",  Bjeta_, "Bjeta[nB]/F");
  data_->Branch("Bjphi",  Bjphi_, "Bjphi[nB]/F");
  data_->Branch("Bjm",    Bjm_,   "Bjm[nB]/F");
  data_->Branch("nJPsi",  &nJPsi_,   "nJPsi/I");
  data_->Branch("JPsipt",  JPsipt_,  "JPsipt[nJPsi]/F");
  data_->Branch("JPsieta", JPsieta_, "JPsieta[nJPsi]/F");
  data_->Branch("JPsiphi", JPsiphi_, "JPsiphi[nJPsi]/F");
  data_->Branch("JPsim",   JPsim_,   "JPsim[nJPsi]/F");
  data_->Branch("nL",   &nL_,   "nL/I");
  data_->Branch("Lid",   Lid_,  "Lid[nL]/F");
  data_->Branch("Lpt",   Lpt_,  "Lpt[nL]/F");
  data_->Branch("Leta",  Leta_, "Leta[nL]/F");
  data_->Branch("Lphi",  Lphi_, "Lphi[nL]/F");
  data_->Branch("Lm",    Lm_,   "Lm[nL]/F");

}

void FragmentationAnalyzer::analyze(const edm::Event &evt, const edm::EventSetup &setup) {
  std::vector<int> leptons,bHadrons,JPsi;

  //gen jets
  edm::Handle<std::vector<reco::GenJet>> genJets;
  evt.getByToken(genJetsToken_,genJets);

  //gen particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(prunedGenParticlesToken_,genParticles);
  for(size_t i = 0; i  <  genParticles->size(); ++i) {
    const reco::GenParticle &p = (*genParticles)[i];
    const reco::Candidate *mother = p.mother();
    int absid = p.pdgId();
    if(p.pt() == 0) continue;

    //iso leptons
    if(absid==11 || absid==13) {
      if(p.pt()<20 || fabs(p.eta())>2.5 || mother==0 || abs(mother->pdgId())!=24) continue;
      leptons.push_back(i);
    }

    //B hadrons
    if(!IS_BHADRON_PDGID(absid)) continue;

    int n = p.numberOfDaughters();
    if(n<2) continue;

    bool hasBDaughter = false;
    bool hasNuDaughter = false;
    bool hasJPsiDaughter = false;
    //std::vector<int> tmpJPsi;
    for(int j = 0; j < n; ++j) {
      const reco::Candidate *d = p.daughter(j);
      int daugId = d->pdgId();
      //J/Psi
      /*
      if(abs(daugId)==443) {
        hasJPsiDaughter = true;
        tmpJPsi.push_back(j);
      }
      */

      if(IS_BHADRON_PDGID(daugId)) {
        hasBDaughter = true;
        break;
      }
      if(IS_NEUTRINO_PDGID(daugId)) hasNuDaughter = true;
    }
    if(hasBDaughter) continue;

    bHadrons.push_back(i);
    //for(auto & ij : tmpJPsi) JPsi.push_back(ij); //Add only J/Psi found in B hadron (not e.g. B*)
    for(int j = 0; j < n; ++j) {
      const reco::Candidate *d = p.daughter(j);
      int daugId = d->pdgId();
      //J/Psi
      if(abs(daugId)==443) {
        hasJPsiDaughter = true;
        JPsi.push_back(j);
      }
    }

    hists["genBHadronJPsiDecay"]->Fill(hasJPsiDaughter);

    //Weakly decaying B hadron
    hists["genBHadronNuDecay"]->Fill(hasNuDaughter);

    //Fragmentation distribution
    for(auto ijet : *genJets) {
      if(p.pt()==0 || ijet.pt()==0) continue;
      double dr = deltaR(p,ijet);

      //Simple dR match of hadron an GenJet
      if(dr < 0.5) {
        double xb = p.pt()/ijet.pt();
        hists["genBHadronPtFraction"]->Fill(xb);
        break;
      }
    }
  }

  if(bHadrons.size()==0) return;

  //for b-hadron studies
  nB_ = 0;
  nJPsi_=0;
  for(auto ijet : *genJets) {
    //fiducial cuts
    if(ijet.pt()<20 || fabs(ijet.eta()>2.5)) continue;

    //quality cuts
    if(ijet.hadEnergy()/ijet.energy()<0.05) continue;
    if(ijet.emEnergy()/ijet.energy()>0.95) continue;
    //if(ijet.getGenConstituents().size()<2) continue;

    //check if a B hadron can be matched
    bool isBjet(false);
    //for(auto & k : bHadrons) {
    for(size_t k=0; k<bHadrons.size(); ++k) {
      const reco::GenParticle &bhad = (*genParticles)[bHadrons[k]];
      float dR=deltaR(bhad,ijet);
      if(dR>0.5) continue;

      Bid_[nB_]  = bhad.pdgId();
      Bpt_[nB_]  = bhad.pt();
      Beta_[nB_] = bhad.eta();
      Bphi_[nB_] = bhad.phi();
      Bm_[nB_]   = bhad.mass();
      isBjet = true;
      break;
    }
    if(!isBjet) continue;

    //Looking for J/Psi in B jet
    for(size_t k=0; k<JPsi.size(); ++k) {
      const reco::GenParticle &jpsi = (*genParticles)[JPsi[k]];
      float dR=deltaR(jpsi,ijet);
      if(dR>0.5) continue;

      std::cout << jpsi.pt() << std::endl;
      JPsipt_[nJPsi_]  = jpsi.pt();
      JPsieta_[nJPsi_] = jpsi.eta();
      JPsiphi_[nJPsi_] = jpsi.phi();
      JPsim_[nJPsi_]   = jpsi.mass();
      nJPsi_++;
    }

    Bjpt_[nB_]  = ijet.pt();
    Bjeta_[nB_] = ijet.eta();
    Bjphi_[nB_] = ijet.phi();
    Bjm_[nB_]   = ijet.mass();
    nB_++;
  }

  //Leptons from W may be from semileptonic hadron decays (at least in Sherpa)
  nL_ = 0;
  //for(auto & ilep : leptons) {
  for(size_t ilep=0; ilep<leptons.size(); ++ilep) {
    const reco::GenParticle &lep = (*genParticles)[leptons[ilep]];
    float minDR(9999.);
    //for(auto & ib : bHadrons) {
    for(size_t ib=0; ib<bHadrons.size(); ++ib) {
      const reco::GenParticle &bhad = (*genParticles)[bHadrons[ib]];
      float dR=deltaR(lep,bhad);
      if(minDR<dR) minDR=dR;
    }
    if(minDR<0.5) continue;
    Lpt_[nL_]  = lep.pt();
    Leta_[nL_] = lep.eta();
    Lphi_[nL_] = lep.phi();
    Lm_[nL_]   = lep.mass();
    nL_++;
  }

  //Fill ntuple
  if(nL_ && nB_) data_->Fill();

}


DEFINE_FWK_MODULE(FragmentationAnalyzer);
