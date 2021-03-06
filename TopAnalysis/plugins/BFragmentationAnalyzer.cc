#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "TLorentzVector.h"

#include "TopLJets2015/TopAnalysis/interface/BFragmentationAnalyzerUtils.h"
#include "TopLJets2015/TopAnalysis/interface/FragEvent.h"

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
  std::vector<double> hadronUncDzb_;
  std::vector<double> hadronUncDz_;
  std::map<std::string, TH1F *> histos_;
  std::map<std::string, TH2F *> histos2_;
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;
  std::vector<int> numEntries_;
  /*
  TTree *data_;
  FragEvent_t ev_;
  */
};


//
FragmentationAnalyzer::FragmentationAnalyzer(const edm::ParameterSet& iConfig) :
  hadronList_(iConfig.getParameter<std::vector<int> >("hadronList")),
  hadronUncDzb_(iConfig.getParameter<std::vector<double> >("hadronUncDzb")),
  hadronUncDz_(iConfig.getParameter<std::vector<double> >("hadronUncDz")),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("pseudoTop:jets"))),
  genParticlesToken_(consumes<std::vector<reco::GenParticle> >(edm::InputTag("pseudoTop"))),
  prunedGenParticlesToken_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"))),
  numEntries_(iConfig.getParameter<std::vector<int> >("numEntries"))
{
  //prepare monitoring histograms
  size_t nhadrons=hadronList_.size()+1;
  histos_["semilepbr"]    = fs->make<TH1F>("semilepbr", ";B-hadron;BR",nhadrons,0,nhadrons);
  histos2_["semilepbr"]    = fs->make<TH2F>("semilepbr2d", ";B-hadron;BR",nhadrons,0,nhadrons,nhadrons,0,nhadrons);
  histos_["semilepbrinc"] = fs->make<TH1F>("semilepbrinc", ";B-hadron;BR",nhadrons,0,nhadrons);
  histos_["bzunc"] = fs->make<TH1F>("bzunc", ";B^{0}-hadron;Unc",nhadrons,0,nhadrons);
  histos_["bpmunc"] = fs->make<TH1F>("bpmunc", ";B^{#pm}-hadron;Unc",nhadrons,0,nhadrons);
  histos_["bsunc"] = fs->make<TH1F>("bsunc", ";B^{0}_{s}-hadron;Unc",nhadrons,0,nhadrons);
  histos_["lbunc"] = fs->make<TH1F>("lbunc", ";#Lambda^{0}_{s}-hadron;Unc",nhadrons,0,nhadrons);
  for(size_t i=0; i<nhadrons; i++)
    {
      std::string name("inc");
      if(i) { char buf[20]; sprintf(buf,"%d", hadronList_[i-1]); name=buf; }
      histos_["semilepbr"]->GetXaxis()->SetBinLabel(i+1,name.c_str());
      histos_["semilepbrinc"]->GetXaxis()->SetBinLabel(i+1,name.c_str());
      histos_["bzunc"]->GetXaxis()->SetBinLabel(i+1,name.c_str());
      histos_["bpmunc"]->GetXaxis()->SetBinLabel(i+1,name.c_str());
      histos_["bsunc"]->GetXaxis()->SetBinLabel(i+1,name.c_str());
      histos_["lbunc"]->GetXaxis()->SetBinLabel(i+1,name.c_str());
      histos_["xb_"+name] = fs->make<TH1F>(("xb_"+name).c_str(), (name+";x_{b}=p_{T}(B)/p_{T}(jet); Jets").c_str(), 100, 0, 2);
      histos_["massDs"+name] = fs->make<TH1F>(("massDs"+name).c_str(), "massDs", 1000, 1.9, 2.2);
      histos_["massDs_ss"+name] = fs->make<TH1F>(("massDs_ss"+name).c_str(), "massDs_ss", 1000, 1.9, 2.2);
      histos_["massDsmD0"+name] = fs->make<TH1F>(("massDsmD0"+name).c_str(), "massDsmD0", 1000, 0.14, 0.17);
      histos_["massDsmD0_ss"+name] = fs->make<TH1F>(("massDsmD0_ss"+name).c_str(), "massDsmD0_ss", 1000, 0.14, 0.17);
      histos_["xb_semilep"+name] = fs->make<TH1F>(("xb_semilep"+name).c_str(), ("xb_semilep"+name+";#it{x}_{b}=#it{p}_{T}(B)/#it{p}_{T}(jet); Jets").c_str(), 100, 0, 2);
      if(i==0) histos_["xb_semilep421"] = fs->make<TH1F>(("xb_semilep421"), ("xb_semilep421;#it{x}_{b}=#it{p}_{T}(B)/#it{p}_{T}(jet); Jets"), 100, 0, 2);
      histos_["xb_semilepDzb"+name] = fs->make<TH1F>(("xb_semilepDzb"+name).c_str(), (name+";#it{x}_{b}=#it{p}_{T}(B)/#it{p}_{T}(jet); Jets").c_str(), 100, 0, 2);
      histos_["xb_semilepDzbu"+name] = fs->make<TH1F>(("xb_semilepDzbu"+name).c_str(), (name+";#it{x}_{b}=#it{p}_{T}(B)/#it{p}_{T}(jet); Jets").c_str(), 100, 0, 2);
      histos_["xb_semilepDzbd"+name] = fs->make<TH1F>(("xb_semilepDzbd"+name).c_str(), (name+";#it{x}_{b}=#it{p}_{T}(B)/#it{p}_{T}(jet); Jets").c_str(), 100, 0, 2);
      histos_["xb_semilepDzu"+name] = fs->make<TH1F>(("xb_semilepDzu"+name).c_str(), (name+";#it{x}_{b}=#it{p}_{T}(B)/#it{p}_{T}(jet); Jets").c_str(), 100, 0, 2);
      histos_["xb_semilepDzd"+name] = fs->make<TH1F>(("xb_semilepDzd"+name).c_str(), (name+";#it{x}_{b}=#it{p}_{T}(B)/#it{p}_{T}(jet); Jets").c_str(), 100, 0, 2);
      histos_["z_semilep"+name] = fs->make<TH1F>(("z_semilep"+name).c_str(), (name+";z=#it{p}_{T}(B)/#it{p}_{T}(b); Jets").c_str(), 100, 0, 2);
      histos2_["zVxb_semilep"+name] = fs->make<TH2F>(("zVxb_semilep"+name).c_str(), (name+";z=#it{p}_{T}(B)/#it{p}_{T}(b); #it{x}_{b}=#it{p}_{T}(B)/#it{p}_{T}(jet)").c_str(), 100, 0, 2, 100, 0, 2);
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
  edm::Handle<std::vector<reco::GenParticle> > genParticles;
  iEvent.getByToken(genJetsToken_,genJets);  
  iEvent.getByToken(genParticlesToken_,genParticles);  
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
          if(jinfo.charmId==-421)
            histos_["xb_semilepDzbinc"]->Fill(jinfo.xb);

          float wgt = 1.;
          //BR up
          if(jinfo.leadTagId==511 && jinfo.charmId==-421)
            wgt = (1. + hadronUncDzb_[0]/100);
          if(jinfo.leadTagId==521 && jinfo.charmId==-421)
            wgt = (1. + hadronUncDzb_[1]/100);
          if(jinfo.leadTagId==531 && jinfo.charmId==-421)
            wgt = (1. + hadronUncDzb_[2]/100);
          if(jinfo.leadTagId==5122 && jinfo.charmId==-421)
            wgt = (1. + hadronUncDzb_[3]/100);

          if(jinfo.leadTagId==511 && jinfo.charmId==421)
            wgt = (1. + hadronUncDz_[0]/100);
          if(jinfo.leadTagId==521 && jinfo.charmId==421)
            wgt = (1. + hadronUncDz_[1]/100);
          if(jinfo.leadTagId==531 && jinfo.charmId==421)
            wgt = (1. + hadronUncDz_[2]/100);
          if(jinfo.leadTagId==5122 && jinfo.charmId==421)
            wgt = (1. + hadronUncDz_[3]/100);

          histos_["xb_semilepDzbuinc"]->Fill(jinfo.xb, wgt);
          histos_["xb_semilepDzuinc"]->Fill(jinfo.xb, wgt);

          if(jinfo.leadTagId==511 && jinfo.charmId==-421)
            histos_["xb_semilepDzbu511"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==511 && jinfo.charmId==421)
            histos_["xb_semilepDzu511"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==521 && jinfo.charmId==-421)
            histos_["xb_semilepDzbu521"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==521 && jinfo.charmId==421)
            histos_["xb_semilepDzu521"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==531 && jinfo.charmId==-421)
            histos_["xb_semilepDzbu531"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==531 && jinfo.charmId==421)
            histos_["xb_semilepDzu531"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==5122 && jinfo.charmId==-421)
            histos_["xb_semilepDzbu5122"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==5122 && jinfo.charmId==421)
            histos_["xb_semilepDzu5122"]->Fill(jinfo.xb, wgt);
          if(jinfo.motherId==413 && jinfo.charmId==-421)
            histos_["xb_semilepDzbu413"]->Fill(jinfo.xb, 1. + hadronUncDzb_[4]/100);

          if(jinfo.motherId==413 && jinfo.charmId==421)
            histos_["xb_semilepDzu413"]->Fill(jinfo.xb, 1. + hadronUncDz_[4]/100);
          if(jinfo.motherId==423 && jinfo.charmId==-421)
            histos_["xb_semilepDzbu423"]->Fill(jinfo.xb, 1. + hadronUncDzb_[5]/100);
          if(jinfo.motherId==423 && jinfo.charmId==421)
            histos_["xb_semilepDzu423"]->Fill(jinfo.xb, 1. + hadronUncDz_[5]/100);
          if(jinfo.motherId==431 && jinfo.charmId==-421)
            histos_["xb_semilepDzbu431"]->Fill(jinfo.xb, 1. + hadronUncDzb_[6]/100);
          if(jinfo.motherId==431 && jinfo.charmId==421)
            histos_["xb_semilepDzu431"]->Fill(jinfo.xb, 1. + hadronUncDz_[6]/100);

          //BR down
          if(jinfo.leadTagId==511 && jinfo.charmId==-421)
            wgt = (1. - hadronUncDzb_[0]/100);
          if(jinfo.leadTagId==521 && jinfo.charmId==-421)
            wgt = (1. - hadronUncDzb_[1]/100);
          if(jinfo.leadTagId==531 && jinfo.charmId==-421)
            wgt = (1. - hadronUncDzb_[2]/100);
          if(jinfo.leadTagId==5122 && jinfo.charmId==-421)
            wgt = (1. - hadronUncDzb_[3]/100);

          if(jinfo.leadTagId==511 && jinfo.charmId==421)
            wgt = (1. - hadronUncDz_[0]/100);
          if(jinfo.leadTagId==521 && jinfo.charmId==421)
            wgt = (1. - hadronUncDz_[1]/100);
          if(jinfo.leadTagId==531 && jinfo.charmId==421)
            wgt = (1. - hadronUncDz_[2]/100);
          if(jinfo.leadTagId==5122 && jinfo.charmId==421)
            wgt = (1. - hadronUncDz_[3]/100);

          histos_["xb_semilepDzbdinc"]->Fill(jinfo.xb, wgt);
          histos_["xb_semilepDzdinc"]->Fill(jinfo.xb, wgt);

          if(jinfo.leadTagId==511 && jinfo.charmId==-421)
            histos_["xb_semilepDzbd511"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==511 && jinfo.charmId==421)
            histos_["xb_semilepDzd511"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==521 && jinfo.charmId==-421)
            histos_["xb_semilepDzbd521"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==521 && jinfo.charmId==421)
            histos_["xb_semilepDzd521"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==531 && jinfo.charmId==-421)
            histos_["xb_semilepDzbd531"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==531 && jinfo.charmId==421)
            histos_["xb_semilepDzd531"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==5122 && jinfo.charmId==-421)
            histos_["xb_semilepDzbd5122"]->Fill(jinfo.xb, wgt);
          if(jinfo.leadTagId==5122 && jinfo.charmId==421)
            histos_["xb_semilepDzd5122"]->Fill(jinfo.xb, wgt);
          if(jinfo.motherId==423 && jinfo.charmId==-421)
            histos_["xb_semilepDzbu423"]->Fill(jinfo.xb, 1. - hadronUncDzb_[5]/100);
          if(jinfo.motherId==423 && jinfo.charmId==421)
            histos_["xb_semilepDzu423"]->Fill(jinfo.xb, 1. - hadronUncDz_[5]/100);
          if(jinfo.motherId==431 && jinfo.charmId==-421)
            histos_["xb_semilepDzbu431"]->Fill(jinfo.xb, 1. - hadronUncDzb_[6]/100);
          if(jinfo.motherId==431 && jinfo.charmId==421)
            histos_["xb_semilepDzu431"]->Fill(jinfo.xb, 1. - hadronUncDz_[6]/100);

          //if(jinfo.motherId>0 && jinfo.charmId==421) std::cout << jinfo.motherId << std::endl;

          if(jinfo.leadTagId==511)
            histos_["xb_semilep511"]->Fill(jinfo.xb);
          if(jinfo.leadTagId==521)
            histos_["xb_semilep521"]->Fill(jinfo.xb);
          if(jinfo.leadTagId==531)
            histos_["xb_semilep531"]->Fill(jinfo.xb);
          if(jinfo.leadTagId==5122)
            histos_["xb_semilep5122"]->Fill(jinfo.xb);
	}
      histos_["xb_inc"]->Fill(jinfo.xb);

      //find parent quark
      /*
      edm::Handle<reco::GenParticleCollection> genParticles;
      iEvent.getByToken(genParticlesToken_,genParticles);
      for(size_t i = 0; i < genParticles->size(); i++) {
        const reco::GenParticle &genIt = (*genParticles)[i];
        float deta = genIt.eta() - jinfo.eta;
        float dphi = std::abs(genIt.phi() - jinfo.phi);
        if (dphi > M_PI)
          dphi -= (2 * M_PI);
	if((deta * deta + dphi * dphi) > 0.01*0.01) continue;
        for(size_t ipf=0; ipf<genIt.numberOfDaughters(); ipf++) {
          const reco::Candidate *daug=genIt.daughter(ipf);
          while(abs(daug->pdgId()) > 6) {
            if(daug->mother() == 0) break;
            daug = daug->mother();
          }
          //if(abs(daug->pdgId()) <= 6) std::cout << "Parent quark found! PDGID= " << daug->pdgId() << std::endl;
          if(jinfo.hasSemiLepDecay)  {
            histos_["z_semilepinc"]->Fill(jinfo.pt/daug->pt());
            histos_["zVxb_semilepinc"]->Fill(jinfo.pt/daug->pt(),jinfo.xb);
          }
        }
      }
      */
      
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
      const float gMassMu(0.1057),gMassK(0.4937),gMassPi(0.1396);
      for(size_t i = 0; i < jinfo.meson.size(); i++) {
        if(abs(jinfo.charmId)!=421) continue;
        for(size_t j = 0; j < jinfo.meson.size(); j++) {
          if(i==j) continue;
          if(jinfo.mesonId[i]==jinfo.mesonId[j]) continue;
          TLorentzVector p4p;
          TLorentzVector p4k;
          p4p.SetPtEtaPhiM(jinfo.meson[i][0], jinfo.meson[i][1], jinfo.meson[i][2], gMassPi);
          p4k.SetPtEtaPhiM(jinfo.meson[j][0], jinfo.meson[j][1], jinfo.meson[j][2], gMassK);
          float mass12 = (p4p+p4k).M();
          if(mass12<1.7 || mass12>2.0) continue;
          histos_["xb_semilep421"]->Fill(jinfo.xb);
          for(size_t k = 0; k < jinfo.meson.size(); k++) {
            if(i==k) continue;
            //std::cout << jinfo.mesonId[i] << " " << jinfo.mesonId[j] << " " << jinfo.mesonId[k] << std::endl;
            if(jinfo.mesonId[i]==jinfo.mesonId[k]) continue;
            TLorentzVector p4p2;
            p4k.SetPtEtaPhiM(jinfo.meson[k][0], jinfo.meson[k][1], jinfo.meson[k][2], gMassK);
            float mass123 = (p4p+p4k+p4p2).M();// - mass12;
            //if(jinfo.mesonId[i]==jinfo.mesonId[k]) histos_["massDsinc"]->Fill(mass123);
            //else histos_["massDs_ssinc"]->Fill(mass123);
            //std::cout << "mass123=" << mass123 << std::endl;
	    if (mass123>1.9 && mass123<2.2) {
              if(jinfo.mesonId[i]==jinfo.mesonId[k]) {
                histos_["massDsinc"]->Fill(mass123);
                histos_["massDsmD0inc"]->Fill(mass123-mass12);
              }
              else {
                histos_["massDs_ssinc"]->Fill(mass123);
                histos_["massDsmD0_ssinc"]->Fill(mass123-mass12);
              }
            //if(mass123>0.14 && mass123<0.17) {
              //std::cout << mass123 << std::endl;
              histos_["xb_semilep413"]->Fill(jinfo.xb);

              if(jinfo.charmId==-421)
                histos_["xb_semilepDzbu413"]->Fill(jinfo.xb, 1. + hadronUncDzb_[4]/100);
              if(jinfo.charmId==421)
                histos_["xb_semilepDzu413"]->Fill(jinfo.xb, 1. + hadronUncDz_[4]/100);
              if(jinfo.charmId==-421)
                histos_["xb_semilepDzbd413"]->Fill(jinfo.xb, 1. - hadronUncDzb_[4]/100);
              if(jinfo.charmId==421)
                histos_["xb_semilepDzd413"]->Fill(jinfo.xb, 1. - hadronUncDz_[4]/100);
            }
          }
        }
      }
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
