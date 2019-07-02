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
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"

//#define IS_BHADRON_PDGID(id) ( ((abs(id)/100)%10 == 5) || (abs(id) >= 5000 && abs(id) <= 5999) )
//#define IS_NEUTRINO_PDGID(id) ( (abs(id) == 12) || (abs(id) == 14) || (abs(id) == 16) )
#define IS_JPSI_PDGID(id) ( (abs(id) == 443) )
#define IS_BQUARK_PDGID(id) ( abs(id) == 5 )



class FragmentationGenAnalyzer : public edm::EDAnalyzer {

 public:
  FragmentationGenAnalyzer(const edm::ParameterSet &);
  virtual void analyze(const edm::Event &, const edm::EventSetup &);
  virtual void fragAnalyze(const edm::Event &, const edm::EventSetup &);
  virtual void endJob() override;

 private:
  std::string fragModel;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;
  std::map<std::string, TH1F*> hists;
  std::map<std::string, TGraph *> wgtGr_;
  std::map< std::string, std::vector<float> > jetWeights;
  TTree *data_;
  Int_t nb_, nB_, nM_, nL_, nJPsi_, Bid_[100], Lid_[100], JPsiPrompt_[100];
  Float_t bpt_[100], xb_[100], beta_[100], bphi_[100];
  Float_t fBpt_[100], Bpt_[100], Beta_[100], Bphi_[100], Bm_[100];
  Float_t Bjpt_[100], Bjeta_[100], Bjphi_[100], Bjm_[100];
  Float_t model_[100];
  Float_t JPsipt_[100], JPsieta_[100], JPsiphi_[100], JPsim_[100], JPsidxy_[100], JPsidz_[100];
  Float_t Lpt_[100], Leta_[100], Lphi_[100], Lm_[100];

};
 
FragmentationGenAnalyzer::FragmentationGenAnalyzer(const edm::ParameterSet &cfg) :
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

  hists["genBHadronJPsiDecay"] = fs->make<TH1F>("genBHadronJPsiDecay", "genBHadronJPsiDecay", 2, 0, 2);
  hists["genBHadronNuDecay"] = fs->make<TH1F>("genBHadronDecay", "genBHadronDecay", 2, 0, 2);
  hists["genBHadronPtFraction"] = fs->make<TH1F>("genBHadronPtFraction", "genBHadronPtFraction", 2, 0, 2);
  for(auto & it : hists) it.second->Sumw2();
  jetWeights["upFrag"]=std::vector<float>();
  jetWeights["centralFrag"]=std::vector<float>();
  jetWeights["downFrag"]=std::vector<float>();
  jetWeights["PetersonFrag"]=std::vector<float>();
  jetWeights["semilepbrUp"]=std::vector<float>();
  jetWeights["semilepbrDown"]=std::vector<float>();

  //sumary tree for B, J/Psi, and leptons
  data_ = fs->make<TTree>("FragTree", "FragTree");
  data_->Branch("nb",    &nb_,    "nb/I");
  data_->Branch("bpt",    bpt_,   "bpt[nb]/F");
  data_->Branch("beta",   beta_,  "beta[nb]/F");
  data_->Branch("bphi",   bphi_,  "bphi[nb]/F");
  data_->Branch("nB",    &nB_,    "nB/I");
  data_->Branch("nM",    &nM_,    "nM/I");
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
  data_->Branch("JPsiPrompt",  JPsiPrompt_,   "JPsiPrompt[nJPsi]/I");
  data_->Branch("JPsipt",  JPsipt_,  "JPsipt[nJPsi]/F");
  data_->Branch("JPsieta", JPsieta_, "JPsieta[nJPsi]/F");
  data_->Branch("JPsiphi", JPsiphi_, "JPsiphi[nJPsi]/F");
  data_->Branch("JPsim",   JPsim_,   "JPsim[nJPsi]/F");
  data_->Branch("JPsidxy", JPsidxy_, "JPsidxy[nJPsi]/F");
  data_->Branch("JPsidz",  JPsidz_,  "JPsidz[nJPsi]/F");
  data_->Branch("nL",   &nL_,   "nL/I");
  data_->Branch("Lid",   Lid_,  "Lid[nL]/F");
  data_->Branch("Lpt",   Lpt_,  "Lpt[nL]/F");
  data_->Branch("Leta",  Leta_, "Leta[nL]/F");
  data_->Branch("Lphi",  Lphi_, "Lphi[nL]/F");
  data_->Branch("Lm",    Lm_,   "Lm[nL]/F");
  data_->Branch("fragModel",   model_,  "fragModel[nM]/F");
  data_->Branch("xb",    xb_,    "xb[nM]F");

}

void FragmentationGenAnalyzer::analyze(const edm::Event &evt, const edm::EventSetup &setup) {
  fragAnalyze(evt, setup);

  //Fill ntuple
  if((nL_ && nB_) || nM_) data_->Fill();
}

void FragmentationGenAnalyzer::endJob() {
  data_->Draw("fragModel:xb","","goff");
  TGraph *g = new TGraph(data_->GetSelectedRows(),data_->GetV2(),data_->GetV1());
  g->SetName(TString(fragModel));
  g->Draw("AP");
  g->Write();
}

void FragmentationGenAnalyzer::fragAnalyze(const edm::Event &evt, const edm::EventSetup &setup) {
  std::vector<int> leptons,bHadrons;//,JPsi;

  //gen jets
  using namespace edm;
  edm::Handle<std::vector<reco::GenJet>> genJets;
  evt.getByToken(genJetsToken_,genJets);
  //playing with fragmentation re-weighting
  nM_ = 0;
  for(auto genJet : *genJets) {
    //map the gen particles which are clustered in this jet
    JetFragInfo_t jinfo=analyzeJet(genJet);
    
    //evaluate the weight to an alternative fragmentation model (if a tag id is available)
    if(jinfo.leadTagId != 0) {
      model_[nM_] = wgtGr_[fragModel]->Eval(jinfo.xb);
      xb_[nM_] = jinfo.xb;
      nM_++;
    }
    /*
    else {
      model_[nM_] = 1;
    }
    */
  }
  data_->Fill();

  //gen particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(prunedGenParticlesToken_,genParticles);
  nb_=0;
  for(size_t i = 0; i  <  genParticles->size(); ++i) {
    const reco::GenParticle &p = (*genParticles)[i];
    const reco::Candidate *mother = p.mother();
    int absid = p.pdgId();
    if(p.pt() == 0) continue;

    //b quark
    if(IS_BQUARK_PDGID(absid)) {
      bpt_[nb_]  = p.pt();
      beta_[nb_] = p.eta();
      bphi_[nb_] = p.phi();
      nb_++;
    }

    //Prompt J/Psi
    /*
    if(IS_JPSI_PDGID(absid)) {
      const reco::Candidate *tmpMother = p.mother();
      while(abs(tmpMother->pdgId()) != 5 && abs(tmpMother->pdgId()) != 22 && abs(tmpMother->pdgId()) != 2212) {
        if(tmpMother->mother() == 0) break;
        tmpMother = tmpMother->mother();
      }
      if(abs(tmpMother->pdgId())==2212) {
        JPsipt_[nJPsi_]  = p.pt();
        JPsieta_[nJPsi_] = p.eta();
        JPsiphi_[nJPsi_] = p.phi();
        JPsim_[nJPsi_]   = p.mass();
        JPsidxy_[nJPsi_] = sqrt(pow(p.vx(),2)+pow(p.vy(),2));
        JPsidz_[nJPsi_]  = p.vz();
        JPsiPrompt_[nJPsi_] = 1;
        nJPsi_++;
      }
    }
    */

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
      if(IS_BHADRON_PDGID(daugId)) {
        hasBDaughter = true;
        break;
      }
      if(IS_NEUTRINO_PDGID(daugId)) hasNuDaughter = true;
    }
    if(hasBDaughter) continue;

    bHadrons.push_back(i);

    //J/Psi
    for(int j = 0; j < n; ++j) {
      const reco::Candidate *d = p.daughter(j);
      if(IS_JPSI_PDGID(d->pdgId())) {
        hasJPsiDaughter = true;
        //JPsi.push_back(j);
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
    //if(ijet.getGenConstituents().size()<2) continue; //FIXME causes unexpected crashes (GenJet constituent is not of the type GenParticle)

    //check if a B hadron can be matched
    bool isBjet(false);
    int b=-1;
    for(size_t k=0; k<bHadrons.size(); ++k) {
      const reco::GenParticle &bhad = (*genParticles)[bHadrons[k]];
      float dR=deltaR(bhad,ijet);
      if(dR>0.5) continue;
      if(bhad.pt()/ijet.pt()>2) continue;

      Bid_[nB_]  = bhad.pdgId();
      Bpt_[nB_]  = bhad.pt();
      Beta_[nB_] = bhad.eta();
      Bphi_[nB_] = bhad.phi();
      Bm_[nB_]   = bhad.mass();
      isBjet = true;
      b=k;
      break;
    }
    if(!isBjet) continue;

    //Looking for J/Psi in B jet
    //for(size_t k=0; k<JPsi.size(); ++k) {
    const reco::GenParticle &bhad = (*genParticles)[bHadrons[b]];
    //const reco::GenParticle &jpsi = (*genParticles)[JPsi[b]];
    for(size_t k = 0; k  < bhad.numberOfDaughters(); ++k) {
      const reco::Candidate *jpsi = bhad.daughter(k);
      if(!IS_JPSI_PDGID(jpsi->pdgId())) continue;
      //float dR=deltaR(jpsi,ijet);
      //if(dR>0->5) continue;
      /*
      const reco::Candidate *tmpMother = bhad.mother();
      while(abs(tmpMother->pdgId()) != 5 && abs(tmpMother->pdgId()) != 22 && abs(tmpMother->pdgId()) != 2212) {
        if(tmpMother->mother() == 0) break;
        tmpMother = tmpMother->mother();
      }
      std::cout << tmpMother->pdgId() << std::endl;
      if(abs(tmpMother->pdgId())!=5) continue;
      */

      JPsipt_[nJPsi_]  = jpsi->pt();
      JPsieta_[nJPsi_] = jpsi->eta();
      JPsiphi_[nJPsi_] = jpsi->phi();
      JPsim_[nJPsi_]   = jpsi->mass();
      JPsidxy_[nJPsi_] = sqrt(pow(jpsi->vx(),2)+pow(jpsi->vy(),2));
      JPsidz_[nJPsi_]  = jpsi->vz();
      JPsiPrompt_[nJPsi_] = 0;
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

}


DEFINE_FWK_MODULE(FragmentationGenAnalyzer);
