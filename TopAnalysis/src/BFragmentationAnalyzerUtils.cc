#include "TopLJets2015/TopAnalysis/interface/BFragmentationAnalyzerUtils.h"
#include <TRandom3.h>
#define IS_JPSI_PDGID(id) ( (abs(id) == 443) )
#define IS_D0_PDGID(id)   ( (abs(id) == 421) )
#define IS_Ds_PDGID(id)   ( (abs(id) == 413) )

//
JetFragInfo_t analyzeJet(const reco::GenJet &genJet,float tagScale)
{
  //loop over the constituents to analyze the jet leading pT tag and the neutrinos
  std::vector< const reco::Candidate * > jconst=genJet.getJetConstituentsQuick();
  const reco::Candidate *leadTagConst=0;
  reco::Candidate::LorentzVector nup4(0,0,0,0);
  bool hasSemiLepDecay(false),hasTauNeutrino(false),hasCharm(false);
  bool hasDspi0(false),hasDsgamma(false);
  int nbtags(0),nctags(0),ntautags(0);
  float pt_charged(0),pt_rand(0),charmId(-1),motherId(-1),bpt,charmPt(-1);
  std::vector< std::vector<double> > meson;
  std::vector<int> mesonId;
  bool hasD0(false), hasPi(false);
  std::vector<TLorentzVector> D0;
  std::vector<TLorentzVector> pi;
  std::vector<short> pdgId;
  TRandom3 *rand = new TRandom3(0);
  for(size_t ijc=0; ijc <jconst.size(); ijc++) 
    {
      const reco::Candidate *par=jconst[ijc];
      int absid=abs(par->pdgId());

      //account for neutrinos for the total energy estimation and check 
      //which ones are coming from B hadron decays
      if(par->status()==1 && IS_NEUTRINO_PDGID(absid)) 
	{
	  nup4 += par->p4()*tagScale;
	  int motherid( abs(par->mother()->pdgId()) );
	  if(absid==16) hasTauNeutrino=true;
	  if(IS_BHADRON_PDGID(motherid)) hasSemiLepDecay=true;
	}
      if(par->status()==1) { //final state sum charged pT
        pt_charged += abs(par->pt() * par->charge());
        if(rand->Uniform(1)<0.02 && par->charge() != 0)
          pt_rand += par->pt();
      }
      if(absid>400 && absid<500) {
      if(!IS_JPSI_PDGID(absid) && !IS_D0_PDGID(absid) && !IS_Ds_PDGID(absid))
        motherId= absid;
        //std::cout << motherId << std::endl;
      }
      if(absid==211) {
        TLorentzVector p4;
        p4.SetPtEtaPhiM(par->pt(), par->eta(), par->phi(), 0.1396);
        pi.push_back(p4);
        pdgId.push_back(par->pdgId());
        hasPi = true;
        meson.push_back(std::vector<double>({par->pt(),par->eta(),par->phi(),par->mass()}));
        mesonId.push_back(par->pdgId());
      }
      else if(absid==321) {
        TLorentzVector p4;
        p4.SetPtEtaPhiM(par->pt(), par->eta(), par->phi(), 0.4937);
        pi.push_back(p4);
        pdgId.push_back(par->pdgId());
      }
      if(absid==111 && IS_Ds_PDGID(motherId)) hasDspi0 = true;
      if(absid==22 && IS_Ds_PDGID(motherId)) hasDsgamma = true;
      if(par->status()!=2) continue;
      
      //count number of tags
      if(absid==15)               ntautags++;
      if(IS_BHADRON_PDGID(absid)) nbtags++;
      if(IS_CHADRON_PDGID(absid)) nctags++;
      if(IS_JPSI_PDGID(absid) || IS_D0_PDGID(absid) || IS_Ds_PDGID(absid)) {
        hasCharm = true;
        if(par->pt()*tagScale > charmPt) { //hardest charmPt
          charmId = par->pdgId();
          charmPt = par->pt() * tagScale;
        }
        if(IS_D0_PDGID(absid)) {
          TLorentzVector p4;
          p4.SetPtEtaPhiM(par->pt()*tagScale, par->eta(), par->phi(), par->mass()*tagScale);
          D0.push_back(p4);
          hasD0 = true;
        }
      }
      /*
      if(par->mother() != nullptr) {
        //if(absid==411 || absid==10411 || absid==10421 || absid==413 || absid==423 || absid==10413 || absid==10423 || absid==20413 || absid==20423 || absid==415 || absid==425 || absid==431 || absid==10431 || absid==433 || absid==10433 || absid==20433 || absid==435)
            //std::cout << par->pdgId() << std::endl;
        if(IS_D0_PDGID(absid)) std::cout << absid << "->" << int(par->mother()->pdgId()) << std::endl;
      }
      */


      //save leading pT tag
      if(leadTagConst && leadTagConst->pt()>par->pt()) continue;
      leadTagConst=par;
    }
      
  //fill the jet info
  JetFragInfo_t jinfo;
  reco::Candidate::LorentzVector totalP4(genJet.p4()+nup4);
  jinfo.xb              = leadTagConst ? (leadTagConst->pt()*tagScale)/totalP4.pt() : -1;
  jinfo.xb_charged      = (leadTagConst && pt_charged>0) ? (leadTagConst->pt()*tagScale)/pt_charged : -1;
  jinfo.xb_rand      = (leadTagConst && pt_rand>0) ? (leadTagConst->pt()*tagScale)/(totalP4.pt() - pt_rand) : -1;
  jinfo.pt              = leadTagConst ? (leadTagConst->pt()*tagScale) : -1;
  jinfo.eta             = leadTagConst ? leadTagConst->eta() : 999;
  jinfo.phi             = leadTagConst ? leadTagConst->phi() : 999;
  jinfo.leadTagId       = leadTagConst ? leadTagConst->pdgId() : 0;
  jinfo.hasSemiLepDecay = hasSemiLepDecay;
  jinfo.hasTauSemiLepDecay = hasTauNeutrino;
  jinfo.nbtags          = nbtags;
  jinfo.nctags          = nctags;
  jinfo.ntautags        = ntautags;
  jinfo.hasCharm        = hasCharm;
  jinfo.charmId         = charmId;
  jinfo.xb_charm        = hasCharm ? charmPt/totalP4.pt() : -1;
  jinfo.xb_charm_rand = (hasCharm && pt_rand>0) ? charmPt/(totalP4.pt() - pt_rand) : -1;
  jinfo.xb_charged_charm = (hasCharm && pt_charged>0) ? charmPt/pt_charged : -1;
  jinfo.xb_charged_charm_rand = (hasCharm && pt_rand>0) ? charmPt/pt_charged : -1;
  jinfo.motherId = motherId;
  jinfo.meson = meson;
  jinfo.mesonId = mesonId;
  jinfo.hasDspi0 = hasDspi0;
  jinfo.hasDsgamma = hasDsgamma;
  jinfo.hasDs = hasD0 && hasPi;
  jinfo.D0 = D0;
  jinfo.pi = pi;
  jinfo.pdgId = pdgId;

  return jinfo;
}
