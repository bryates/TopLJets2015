#include "TopLJets2015/TopAnalysis/interface/BFragmentationAnalyzerUtils.h"
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
  int nbtags(0),nctags(0),ntautags(0);
  float pt_charged(0),charmId(-1),motherId(-1),bpt;
  std::vector< std::vector<double> > meson;
  std::vector<int> mesonId;
  for(size_t ijc=0; ijc <jconst.size(); ijc++) 
    {
      const reco::Candidate *par=jconst[ijc];
      int absid=abs(par->pdgId());

      //account for neutrinos for the total energy estimation and check 
      //which ones are coming from B hadron decays
      if(abs(par->pdgId())==5)
            std::cout << "b quark pT " << bpt << " status= " << par->status() << std::endl;
      if(par->status()==1 && IS_NEUTRINO_PDGID(absid)) 
	{
	  nup4 += par->p4()*tagScale;
	  int motherid( abs(par->mother()->pdgId()) );
	  if(absid==16) hasTauNeutrino=true;
	  if(IS_BHADRON_PDGID(motherid)) hasSemiLepDecay=true;
	}
      if(par->status()==1) { //final state sum charged pT
        pt_charged += abs(par->pt() * par->charge());
      }
      if(absid>400 && absid<500) {
      if(!IS_JPSI_PDGID(absid) && !IS_D0_PDGID(absid) && !IS_Ds_PDGID(absid))
        motherId= absid;
      }
      if(absid==211) {
        meson.push_back(std::vector<double>({par->pt(),par->eta(),par->phi(),par->mass()}));
        mesonId.push_back(par->pdgId());
      }
      if(par->status()!=2) continue;
      
      //count number of tags
      if(absid==15)               ntautags++;
      if(IS_BHADRON_PDGID(absid)) nbtags++;
      if(IS_CHADRON_PDGID(absid)) nctags++;
      if(IS_JPSI_PDGID(absid) || IS_D0_PDGID(absid) || IS_Ds_PDGID(absid)) {
        hasCharm = true;
        charmId = par->pdgId();
      }


      //save leading pT tag
      if(leadTagConst && leadTagConst->pt()>par->pt()) continue;
      leadTagConst=par;
    }
      
  //fill the jet info
  JetFragInfo_t jinfo;
  reco::Candidate::LorentzVector totalP4(genJet.p4()+nup4);
  jinfo.xb              = leadTagConst ? (leadTagConst->pt()*tagScale)/totalP4.pt() : -1;
  jinfo.xb_charged      = (leadTagConst && pt_charged>0) ? (leadTagConst->pt()*tagScale)/pt_charged : -1;
  jinfo.pt              = leadTagConst ? (leadTagConst->pt()*tagScale) : -1;
  jinfo.eta             = leadTagConst ? leadTagConst->eta() : 999;
  jinfo.phi             = leadTagConst ? leadTagConst->phi() : 999;
  jinfo.leadTagId       = leadTagConst ? leadTagConst->pdgId() : 0;
  jinfo.hasSemiLepDecay = hasSemiLepDecay;
  jinfo.hasTauSemiLepDecay = hasTauNeutrino;
  jinfo.nbtags          = nbtags;
  jinfo.nctags          = nctags;
  jinfo.ntautags        = ntautags;
  jinfo.hasCharm = hasCharm;
  jinfo.charmId = charmId;
  jinfo.motherId = motherId;
  jinfo.meson = meson;
  jinfo.mesonId = mesonId;

  return jinfo;
}
