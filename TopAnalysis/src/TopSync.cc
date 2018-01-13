#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CharmEvent.h"
#include "TopLJets2015/TopAnalysis/interface/TopSync.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"
#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"

//#include "TopLJets2015/TopAnalysis/interface/OtherFunctions.h"
#include "TopLJets2015/TopAnalysis/interface/Trigger.h"
#include "TopLJets2015/TopAnalysis/interface/Particle.h"
#include "TopLJets2015/TopAnalysis/interface/Leptons.h"
#include "TopLJets2015/TopAnalysis/interface/Jet.h"
#include "TopLJets2015/TopAnalysis/interface/StdPlots.h"
#include "TopLJets2015/TopAnalysis/interface/CharmTree.h"
#include "TopLJets2015/TopAnalysis/interface/KalmanEvent.h"

#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

using namespace std;

void RunTopSync(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts,
                 TString era,
                 TString runPeriod,
                 Bool_t debug=false)
{
  if(debug) cout << "in RunTopSync" << endl;

  bool isTTbar( filename.Contains("_TTJets") );


  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << (runSysts ? " syst variations will be considered" : "") << endl;
  
  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  allPlots["nevt_mm"]     = new TH1F("nevt_mm",";N_{events};Events" ,5,0.,5.);
  allPlots["nevt_ee"]     = new TH1F("nevt_ee",";N_{events};Events" ,5,0.,5.);
  allPlots["nevt_em"]     = new TH1F("nevt_em",";N_{events};Events" ,5,0.,5.);
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      //evch = {};
      t->GetEntry(iev);
      if(iev%5000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      allPlots["nevt_mm"]->Fill(0);
      allPlots["nevt_ee"]->Fill(0);
      allPlots["nevt_em"]->Fill(0);

      //Basic lepton kinematics
      std::vector<int> tightLeptons,vetoLeptons;
      Leptons Muons(Tight,debug);
      Leptons Electrons(Tight,debug);

      Muons.setMinPt(20);
      Muons.setMaxEta(2.4);
      Muons.setMaxRelIso(0.15);

      Electrons.setMinPt(20);
      Electrons.setMaxEta(2.4);
      Electrons.setMaxRelIso(0.0678);

      for(int il=0; il<ev.nl; il++)
	{
          //cout << "in lepton selection" << endl;
          Particle p(ev.l_pt[il], ev.l_eta[il], ev.l_phi[il], ev.l_mass[il], ev.l_id[il]*ev.l_charge[il], ev.l_relIso[il], ev.l_pid[il]);
          if(p.isMuon()) Muons.addParticle(p); //only accepts tight
          else if(p.isElectron()) Electrons.addParticle(p);
	}

      if(debug) cout << "lepton selection DONE" << endl;
      Leptons leptons(Tight,debug);
      leptons.combineLeptons(Muons);
      leptons.combineLeptons(Electrons);
      if(debug) cout << "sorting leptons" << endl;
      leptons.sortLeptonsByPt();
      if(leptons.size()<2) continue; // di-lepton events
      int ch = leptons[0].getPdgId()*leptons[1].getPdgId();
      if(ch > 0) continue; // opposite sign only
      TString chTag("");
      if(ch==-13*13) chTag = "_mm";
      else if(ch==-11*11) chTag = "_ee";
      else if(ch==-11*13) chTag = "_em";
      else continue;
      if((leptons[0].getVec()+leptons[1].getVec()).M()<20) continue; // mass > 20 GeV

      allPlots["nevt"+chTag]->Fill(1);
      if(fabs((leptons[0].getVec()+leptons[1].getVec()).M()-91)<15) continue; // mZ veto
      allPlots["nevt"+chTag]->Fill(2);

      std::vector<Jet> kJetsVec, lightJetsVec, allJetsVec;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

	  //cross clean with respect to leptons deltaR<0.4
          bool overlapsWithLepton(false);
          for(size_t il=0; il<leptons.size(); il++) {
            if(jp4.DeltaR(leptons[il].getVec())>0.4) continue;  //Jet is fine
	    overlapsWithLepton=true;                   //Jet ovelaps with an "isolated" lepton, event is bad
          }
          if(overlapsWithLepton) continue;
          if(debug) cout << "Overlap with lepton DONE" << endl;

	  //smear jet energy resolution for MC
	  //jetDiff -= jp4;
	  float genJet_pt(0);
	  if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }
	  //jetDiff += jp4;

	  // re-inforce kinematics cuts
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  //if(leadingJetIdx<0) leadingJetIdx=k;

	  //save jet
          Jet tmpj(jp4, 0, k, ev.j_pt_charged[k], ev.j_pz_charged[k], ev.j_p_charged[k], ev.j_pt_pf[k], ev.j_pz_pf[k], ev.j_p_pf[k], ev.j_g[k]); //Store pt of charged and total PF tracks and gen matched index
          allJetsVec.push_back(tmpj);
	}

      //MET and transverse mass
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
      met.SetPz(0.); met.SetE(met.Pt());

      if((ch==-13*13 || ch==-11*11) && met.Pt()<40) continue;
      allPlots["nevt"+chTag]->Fill(3);

      if(allJetsVec.size()<2) continue; // 2+ jets
      allPlots["nevt"+chTag]->Fill(4);
    }

  //close input file
  f->Close();

  //save histos to file  
  TString selPrefix("");  
  TString selPostfix("");  
  if(flavourSplitting!=NOFLAVOURSPLITTING) selPrefix=Form("%d_",flavourSplitting);
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+selPrefix+baseName,"RECREATE");
  fOut->cd();
  if(debug) cout << "writing histograms" << endl;

  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  if(debug) cout << "writing histograms DONE" << endl;
  fOut->Close();
}
