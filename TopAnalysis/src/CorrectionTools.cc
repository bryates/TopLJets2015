#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"



TH2 *pf_eff_;
//
std::vector<TGraph *> getPileupWeights(TString era,TH1 *genPU)
{
  std::vector<TGraph *>puWgtGr;
  if(genPU==0) return  puWgtGr;

  if(genPU->GetNbinsX()==1000) genPU->Rebin(10);
  genPU->Scale(1./genPU->Integral());

  //readout the pileup weights and take the ratio of data/MC
  TString puWgtUrl(era+"/pileupWgts.root");
  gSystem->ExpandPathName(puWgtUrl);
  TFile *fIn=TFile::Open(puWgtUrl);
  for(size_t i=0; i<3; i++)
    {
      TString grName("pu_nom");
      if(i==1) grName="pu_down";
      if(i==2) grName="pu_up";
      TGraph *puData=(TGraph *)fIn->Get(grName);
      Float_t totalData=puData->Integral();
      TH1 *tmp=(TH1 *)genPU->Clone("tmp");
      for(Int_t xbin=1; xbin<=tmp->GetXaxis()->GetNbins(); xbin++)
	{
	  Float_t yexp=genPU->GetBinContent(xbin);
	  Double_t xobs,yobs;
	  puData->GetPoint(xbin-1,xobs,yobs);
	  tmp->SetBinContent(xbin, yexp>0 ? yobs/(totalData*yexp) : 0. );
	}
      TGraph *gr=new TGraph(tmp);
      grName.ReplaceAll("pu","puwgts");
      gr->SetName(grName);
      puWgtGr.push_back( gr );
      tmp->Delete();
    }
  return puWgtGr;
}


//apply jet energy resolutions
MiniEvent_t smearJetEnergies(MiniEvent_t ev, std::string option) {
  for (int k = 0; k < ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    //smear jet energy resolution for MC
    float genJet_pt(0);
    if(ev.j_g[k]>-1) genJet_pt = ev.g_pt[ ev.j_g[k] ];
    if(!ev.isData && genJet_pt>0) {
      int smearIdx(0);
      if(option=="up") smearIdx=1;
      if(option=="down") smearIdx=2;
      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt)[smearIdx];
      jp4 *= jerSmear;
      ev.j_pt[k]   = jp4.Pt();
      ev.j_eta[k]  = jp4.Eta();
      ev.j_phi[k]  = jp4.Phi();
      ev.j_mass[k] = jp4.M();
    }
  }
  
  return ev;
}

void smearJetEnergyStochastic(TLorentzVector &jp4, double resolution, std::string option)
{
  TRandom *random = new TRandom3(0);
  int smearIdx(0);
  if(option=="up") smearIdx=1;
  if(option=="down") smearIdx=2;
  float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),0.)[smearIdx];
  float jerFactor = 1 + random->Gaus(0, resolution) * sqrt(std::max(pow(jerSmear, 2) - 1., 0.));
  jp4 *= jerFactor;
  delete random;
}

//see working points in https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco
/*
MiniEvent_t addBTagDecisions(MiniEvent_t ev,float wp) {
  for (int k = 0; k < ev.nj; k++) {
    ev.j_btag[k] = (ev.j_csv[k] > wp);
  }
  
  return ev;
}
*/


//details in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
MiniEvent_t updateBTagDecisions(MiniEvent_t ev, 
				std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> &btvsfReaders,
				std::map<BTagEntry::JetFlavor, TGraphAsymmErrors*> &expBtagEff, 
				std::map<BTagEntry::JetFlavor, TGraphAsymmErrors*> &expBtagEffPy8, 
				BTagSFUtil *myBTagSFUtil, 
				std::string option) {
  for (int k = 0; k < ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    //bool isBTagged(ev.j_btag[k]);
    bool isBTagged(ev.j_csv[k]>0.848);
    if(!ev.isData) {
      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
      float expEff(1.0), jetBtagSF(1.0);
      
      BTagEntry::JetFlavor hadFlav=BTagEntry::FLAV_UDSG;
      if(abs(ev.j_hadflav[k])==4) hadFlav=BTagEntry::FLAV_C;
      if(abs(ev.j_hadflav[k])==5) hadFlav=BTagEntry::FLAV_B;

      expEff    = expBtagEff[hadFlav]->Eval(jptForBtag); 
      jetBtagSF = btvsfReaders[hadFlav]->eval_auto_bounds( option, hadFlav, jetaForBtag, jptForBtag);
      jetBtagSF *= expEff>0 ? expBtagEffPy8[hadFlav]->Eval(jptForBtag)/expBtagEff[hadFlav]->Eval(jptForBtag) : 0.;
      
      //updated b-tagging decision with the data/MC scale factor
      myBTagSFUtil->modifyBTagsWithSF(isBTagged, jetBtagSF, expEff);
      //ev.j_btag[k] = isBTagged;
    }
  }
  
  return ev;
}

//details in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> getBTVcalibrationReaders(TString era,BTagEntry::OperatingPoint btagOP)
{
  //start the btag calibration
  TString btagUncUrl(era+"/CSVv2_Moriond17_B_H.csv");
  gSystem->ExpandPathName(btagUncUrl);
  BTagCalibration btvcalib("csvv2", btagUncUrl.Data());

  //start calibration readers for b,c and udsg separately including the up/down variations
  std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> btvCalibReaders;
  btvCalibReaders[BTagEntry::FLAV_B]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders[BTagEntry::FLAV_B]->load(btvcalib,BTagEntry::FLAV_B,"mujets");
  btvCalibReaders[BTagEntry::FLAV_C]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders[BTagEntry::FLAV_C]->load(btvcalib,BTagEntry::FLAV_C,"mujets");
  btvCalibReaders[BTagEntry::FLAV_UDSG]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders[BTagEntry::FLAV_UDSG]->load(btvcalib,BTagEntry::FLAV_UDSG,"incl");

  //all done
  return btvCalibReaders;
}

//the expections are created with the script scripts/saveExpectedBtagEff.py (cf README)
std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> readExpectedBtagEff(TString era,TString btagExpPostFix)
{
  //open up the ROOT file with the expected efficiencies
  TString btagEffExpUrl(era+"/expTageff.root");
  btagEffExpUrl.ReplaceAll(".root",btagExpPostFix+".root");
  gSystem->ExpandPathName(btagEffExpUrl);
  TFile *beffIn=TFile::Open(btagEffExpUrl);
  
  //read efficiency graphs
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> expBtagEff;
  expBtagEff[BTagEntry::FLAV_B]=(TGraphAsymmErrors *)beffIn->Get("b");
  expBtagEff[BTagEntry::FLAV_C]=(TGraphAsymmErrors *)beffIn->Get("c");
  expBtagEff[BTagEntry::FLAV_UDSG]=(TGraphAsymmErrors *)beffIn->Get("udsg");
  beffIn->Close();

  //all done
  return expBtagEff;
}

void applyJetCorrectionUncertainty(TLorentzVector &jp4,JetCorrectionUncertainty *jecUnc,TString direction)
{
    jecUnc->setJetPt(jp4.Pt());
    jecUnc->setJetEta(jp4.Eta());
    double scale = 1.;
    if (direction == "up")
      scale += jecUnc->getUncertainty(true);
    else if (direction == "down")
      scale -= jecUnc->getUncertainty(false);
    
    jp4 *= scale;
}

std::map<TString, std::map<TString, std::vector<double> > > getTrackingEfficiencyMap(TString era) {
  std::map<TString, std::map<TString, std::vector<double> > > trackEffMap;
  
  if(era.Contains("era2016")) {
    trackEffMap["BCDEF"]["binning"] = {-2.4, -1.5, -0.8, 0.8, 1.5, 2.4};
    trackEffMap["BCDEF"]["nominal"] = {0.93, 1.08, 1.01, 1.08, 0.93};
    trackEffMap["BCDEF"]["unc"]     = {0.04, 0.04, 0.03, 0.04, 0.04};
    for (unsigned int i = 0; i < trackEffMap["BCDEF"]["nominal"].size(); i++) {
      trackEffMap["BCDEF"]["up"].push_back(trackEffMap["BCDEF"]["nominal"][i]+trackEffMap["BCDEF"]["unc"][i]);
      trackEffMap["BCDEF"]["down"].push_back(trackEffMap["BCDEF"]["nominal"][i]-trackEffMap["BCDEF"]["unc"][i]);
    }
    trackEffMap["GH"]["binning"] = {-2.4, -1.5, -0.8, 0.8, 1.5, 2.4};
    trackEffMap["GH"]["nominal"] = {1.12, 1.07, 1.04, 1.07, 1.12};
    trackEffMap["GH"]["unc"]     = {0.05, 0.06, 0.03, 0.06, 0.05};
    for (unsigned int i = 0; i < trackEffMap["GH"]["nominal"].size(); i++) {
      trackEffMap["GH"]["up"].push_back(trackEffMap["GH"]["nominal"][i]+trackEffMap["GH"]["unc"][i]);
      trackEffMap["GH"]["down"].push_back(trackEffMap["GH"]["nominal"][i]-trackEffMap["GH"]["unc"][i]);
    }
    /*
    TString effUrl(era+"/pf_tracks.root");
    gSystem->ExpandPathName(effUrl);
    TFile *fIn=TFile::Open(effUrl);
    //lepEffH_["m_singleleptrig"]=(TH2 *)fIn->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA")->Clone();
    pf_eff_=(TH2 *)fIn->Get("eta_pt")->Clone();
    pf_eff_->SetDirectory(0);
    fIn->Close();
    */
  }
  
  return trackEffMap;
}

//Pion tracking efficiencies from https://github.com/pfs/TopLJets2015/blob/80x_rereco_rev/TopAnalysis/src/CorrectionTools.cc#L391-L503
void applyEtaDepTrackingEfficiencySF(MiniEvent_t &ev, std::vector<double> sfs, std::vector<double> etas, int *id) {
  if (sfs.size() != (etas.size() - 1)) std::cout << "applyEtaDepTrackingEfficiencySF error: need one more bin boundary than scale factors: " << sfs.size() << "," << etas.size() << std::endl;
  for (unsigned int i = 0; i < sfs.size(); i++) {
    applyTrackingEfficiencySF(ev, sfs[i], etas[i], etas[i+1], id);
  }
}

//Testing data driven SFs
void applyTrackingEfficiencySF(MiniEvent_t &ev, TH2 *pf_eff_, int *id) {
  if(ev.isData) return;
  
  TRandom* random = new TRandom3(0); // random seed

  for (int k = 0; k < ev.npf; k++) {
    float minEtaForEff( pf_eff_->GetXaxis()->GetXmin() ), maxEtaForEff( pf_eff_->GetXaxis()->GetXmax()-0.01 );
    float etaForEff=TMath::Max(TMath::Min(ev.pf_eta[k],maxEtaForEff),minEtaForEff);
    Int_t etaBinForEff=pf_eff_->GetXaxis()->FindBin(etaForEff);
    
    float minPtForEff( pf_eff_->GetYaxis()->GetXmin() ), maxPtForEff( pf_eff_->GetYaxis()->GetXmax()-0.01 );
    float ptForEff=TMath::Max(TMath::Min(ev.pf_pt[k],maxPtForEff),minPtForEff);
    Int_t ptBinForEff=pf_eff_->GetYaxis()->FindBin(ptForEff);

    float sf=pf_eff_->GetBinContent(etaBinForEff,ptBinForEff);
    
    if (sf <= 1) {
      //if (abs(ev.pf_id[k]) != 211) continue;
      if (random->Rndm() > sf) {
        id[k] = -1;
        continue;
        //make sure that particle does not pass any cuts
        ev.pf_pt[k]  = 1e-20;
        ev.pf_m[k]   = 1e-20;
        ev.pf_eta[k] = 999.;
        ev.pf_c[k]   = 0;
      }
    }
    /* no need to drop since Kalman tracks are matched to PF in Kalman.cc
    for (int k = 0; k < ev.npf; k++) {
      if (abs(ev.k_pf_id[k]) != 211 && abs(ev.k_pf_id[k]) != 321) continue;
      if (ev.k_pf_eta[k] < minEta) continue;
      if (ev.k_pf_eta[k] > maxEta) continue;
      if (random->Rndm() > sf) {
        //make sure that particle does not pass any cuts
        ev.k_pf_pt[k]  = 1e-20;
        ev.k_mass[k]   = 1e-20;
        ev.k_pf_eta[k] = 999.;
      }
    }
    */
  }
  
  delete random;
}

void applyTrackingEfficiencySF(MiniEvent_t &ev, double sf, double minEta, double maxEta, int *id) {
  if(ev.isData) return;
  
  TRandom* random = new TRandom3(0); // random seed

  if (sf <= 1) {
    for (int k = 0; k < ev.npf; k++) {
      if (abs(ev.pf_id[k]) != 211) continue;
      if (ev.pf_eta[k] < minEta) continue;
      if (ev.pf_eta[k] > maxEta) continue;
      if (random->Rndm() > sf) {
        id[k] = -1;
        continue;
        //make sure that particle does not pass any cuts
        ev.pf_pt[k]  = 1e-20;
        ev.pf_m[k]   = 1e-20;
        ev.pf_eta[k] = 999.;
        ev.pf_c[k]   = 0;
      }
    }
    /* no need to drop since Kalman tracks are matched to PF in Kalman.cc
    for (int k = 0; k < ev.npf; k++) {
      if (abs(ev.k_pf_id[k]) != 211 && abs(ev.k_pf_id[k]) != 321) continue;
      if (ev.k_pf_eta[k] < minEta) continue;
      if (ev.k_pf_eta[k] > maxEta) continue;
      if (random->Rndm() > sf) {
        //make sure that particle does not pass any cuts
        ev.k_pf_pt[k]  = 1e-20;
        ev.k_mass[k]   = 1e-20;
        ev.k_pf_eta[k] = 999.;
      }
    }
    */
  }
  else { // sf > 1
    // find charged hadrons that were not reconstructed
    double dRcut = 0.01*0.01; //squared b/c sqrt is expensive!
    //double dRcut = 0.01;
    std::vector<int> chGenNonRecoHadrons;
    int NchGenHadrons = 0;
    for (int g = 0; g < ev.ngpf; g++) {
      if (ev.gpf_pt[g] < 0.9) continue;
      if (ev.gpf_eta[g] < minEta) continue;
      if (ev.gpf_eta[g] > maxEta) continue;
      if (ev.gpf_c[g] == 0) continue;
      if (abs(ev.gpf_id[g]) < 100) continue;
      NchGenHadrons++;
      bool matched = false;
      for (int k = 0; k < ev.npf; k++) {
        if (ev.pf_pt[k] < 0.8) continue;
        if (abs(ev.pf_id[k]) != 211) continue;
        double dEta = ev.gpf_eta[g] - ev.pf_eta[k];
        double dPhi = TVector2::Phi_mpi_pi(ev.gpf_phi[g] - ev.pf_phi[k]);
        double dR = pow(dEta, 2) + pow(dPhi, 2);
        //double dR = sqrt(pow(dEta, 2) + pow(dPhi, 2)); //sqrt is expensive!
        if (dR < dRcut) {
          matched = true;
          break;
        }
      }
      if (!matched) chGenNonRecoHadrons.push_back(g);
    }
    if (chGenNonRecoHadrons.size() == 0) return;
    double promotionProb = TMath::Min(1., NchGenHadrons*(sf-1.)/chGenNonRecoHadrons.size());
    std::vector<int> chGenNonRecoHadronsToPromote;
    for (const int g : chGenNonRecoHadrons) {
      if (random->Rndm() < promotionProb) { //FIXME
        chGenNonRecoHadronsToPromote.push_back(g);
      }
    }
    for (unsigned int i = 0; i < chGenNonRecoHadronsToPromote.size(); i++) {
      int k = ev.npf + i;
      int g = chGenNonRecoHadronsToPromote[i];
      // jet association
      int j = -1;
      double jetR = 0.4*0.4;
      //double jetR = 0.4;
      for (int ij = 0; ij < ev.nj; ij++) {
        double dEta = ev.gpf_eta[g] - ev.j_eta[ij];
        double dPhi = TVector2::Phi_mpi_pi(ev.gpf_phi[g] - ev.j_phi[ij]);
        double dR = pow(dEta, 2) + pow(dPhi, 2);
        //double dR = sqrt(pow(dEta, 2) + pow(dPhi, 2)); //sqrt is expensive!
        if (dR < jetR) {
          j = ij;
          break;
        }
      }
      id[k] = 1;
      ev.pf_j[k]   = j;
      ev.pf_id[k]  = ev.gpf_id[g];
      ev.pf_c[k]   = ev.gpf_c[g];
      ev.pf_pt[k]  = ev.gpf_pt[g];
      ev.pf_eta[k] = ev.gpf_eta[g];
      ev.pf_phi[k] = ev.gpf_phi[g];
      ev.pf_m[k]   = ev.gpf_m[g];
      ev.pf_mother[k] = ev.gpf_mother[g]; //Store mother ID for D mesons (since Kalman filter wasn't run on GEN particles)
      //Update Kalman tracks
      if(ev.gpf_mother[k]==421 || ev.gpf_mother[k]==413) { //don't need J/Psi since this is pion only
      //if(ev.gpf_mother[k]==443 || ev.gpf_mother[k]==421 || ev.gpf_mother[k]==413) {
        int kpf = ++ev.nkpf;
        //ev.njpsi++; //no J/Psi
        ev.nmeson++;
        ev.k_pf_pt[kpf] = ev.gpf_pt[g];
        ev.k_pf_eta[kpf] = ev.gpf_eta[g];
        ev.k_pf_phi[kpf] = ev.gpf_phi[g];
        ev.k_id[kpf] = ev.gpf_mother[g];
        ev.k_chi2[kpf] = 1.;
        ev.k_vtxProb[kpf] = 1.;
        //save pi and K only (211 and 321)
        if(ev.gpf_id[g]==211 || ev.gpf_id[g]==321) ev.k_pf_id[kpf] = 211;
        if(ev.gpf_id[g]==-211 || ev.gpf_id[g]==-321) ev.k_pf_id[kpf] = -211;
        //ev.k_pf_id[kpf] = ev.gpf_id[g];
        ev.k_l3d[kpf] = 1.;
        ev.k_sigmal3d[kpf] = 0.01; //random to pass sig cut > 10
        ev.k_j[kpf] = j;
        //if(ev.gpf_mother[k]==443) ev.k_mass[kpf] = 3.097; //no J/Psi
        if(ev.gpf_mother[k]==421) ev.k_mass[kpf] = 1.864;
        else if(ev.gpf_mother[k]==413) ev.k_mass[kpf] = 2.010;
      }
      ev.pf_dxy[k] = 0.;
      ev.pf_dz[k]  = 0.;
    }
    ev.npf    = ev.npf + chGenNonRecoHadronsToPromote.size();
  }
  
  delete random;
}

std::vector<RunPeriod_t> getRunPeriods(TString era) {
  std::vector<std::pair<TString,float>> runPeriods;
  if(era.Contains("era2016")) {
    runPeriods.push_back(std::pair<TString,float>("BCDEF", 19712.86));
    runPeriods.push_back(std::pair<TString,float>("GH",    16146.178));
  }

  return runPeriods;
}

TString assignRunPeriod(std::vector<RunPeriod_t> &runPeriods, TRandom *rand) {
  float totalLumi(0.);
  for(auto periodLumi : runPeriods) totalLumi += periodLumi.second;

  //randomly pick a run period based on lumi range
  float pickLumi( rand!=0 ? rand->Uniform(totalLumi) : gRandom->Uniform(totalLumi) );
  float testLumi(0);
  int iLumi(0);
  for(auto periodLumi : runPeriods) {
    testLumi += periodLumi.second;
    if(pickLumi < testLumi) break; //random period selected
    else iLumi++;
  }

  return runPeriods[iLumi].first; //return selected period
  //return iLumi;
}

/*
float customSF(float pt, TString runPeriod) {
  if(runPeriod=="BCDEF") {
    if(pt<25) return 0.8139 + 0.0109 * pt;
    else return 0.9175 + 0.00463 * pt;
  }
  else if(runPeriod=="GH") {
    if(pt<25) return 1.062 + 0.00272 * pt;
    else return 1.065 + 0.00208 * pt;
  }
  else return 1.;
}
*/
