#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"

#include "TFile.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>


using namespace std;

//
LeptonEfficiencyWrapper::LeptonEfficiencyWrapper(bool isData,TString era,TString runPeriod, bool debug)
{
  if(isData) return;
  //if(runPeriod.Contains("B") || runPeriod.Contains("C") || runPeriod.Contains("D") || runPeriod.Contains("E") || runPeriod.Contains("F")) runPeriod = "BCDEF";
  if(TString("BCDEF").Contains(runPeriod)) runPeriod_ = "BCDEF";
  else if(TString("GH").Contains(runPeriod)) runPeriod_ = "GH";
  else runPeriod_ = "BCDEF";
  debug_ = debug;
  if(debug_) std::cout << "Initializing efficiencies for " << era << " runPeriod " << runPeriod_ << std::endl;
  init(era,runPeriod_);
}

//
void LeptonEfficiencyWrapper::init(TString era,TString runPeriod)
{
  //2015 dataset
  if(era.Contains("era2015"))
    {
      era_=2015;
      TString lepEffUrl(era+"/muonEfficiencies.root");
      gSystem->ExpandPathName(lepEffUrl);
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH_["m_sel"]=(TH2 *)fIn->Get("m_sel");
      lepEffH_["m_sel"]->SetDirectory(0);
      lepEffH_["m_singleleptrig"]=(TH2 *)fIn->Get("m_trig");
      lepEffH_["m_singleleptrig"]->SetDirectory(0);
      fIn->Close();
      
      lepEffUrl=era+"/electronEfficiencies.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffH_["e_sel"]=(TH2 *)fIn->Get("EGamma_SF2D");
      lepEffH_["e_sel"]->SetDirectory(0);
      fIn->Close();
    }

  //2016 dataset
  else
    {
      era_=2016;
                            //EfficienciesAndSF_RunBCDEF_23SepReReco.root
      TString lepEffUrl(era+"/EfficienciesAndSF_Run"+runPeriod+"_23SepReReco.root");
      gSystem->ExpandPathName(lepEffUrl);
      TFile *fIn=TFile::Open(lepEffUrl);
      //lepEffH_["m_singleleptrig"]=(TH2 *)fIn->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA")->Clone();
      lepEffH_["m_singleleptrig"]=(TH2 *)fIn->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio")->Clone();
      lepEffH_["m_singleleptrig"]->SetDirectory(0);
      fIn->Close();

                    //Tracking_EfficienciesAndSF_RunBCDEF_23SepReReco.root
      lepEffUrl=era+"/Tracking_EfficienciesAndSF_Run"+runPeriod+"_23SepReReco.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffGr_["m_tk_aeta_"+runPeriod]=(TGraphAsymmErrors *)fIn->Get("ratio_eff_aeta_dr030e030_corr");
      lepEffGr_["m_tk_vtx_"+runPeriod]=(TGraphAsymmErrors *)fIn->Get("ratio_eff_vtx_dr030e030_corr");
      fIn->Close();

                    //MuonID_RunBCDEFF_23SepReReco.root
      lepEffUrl=era+"/MuonID_Run"+runPeriod+"_23SepReReco.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);      
      lepEffH_["m_sel"]=(TH2 *)fIn->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio")->Clone();
      lepEffH_["m_sel"]->SetDirectory(0);
      fIn->Close();

                    //MuonIso_RunBCDEFF_23SepReReco.root
      lepEffUrl=era+"/MuonIso_Run"+runPeriod+"_23SepReReco.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);      
      TH2 *isoH=(TH2 *)fIn->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
      for(Int_t xbin=1; xbin<=(lepEffH_["m_sel"])->GetNbinsX(); xbin++)
	for(Int_t ybin=1; ybin<=(lepEffH_["m_sel"])->GetNbinsY(); ybin++)
	  {
	    float sfid(lepEffH_["m_sel"]->GetBinContent(xbin,ybin)), sfiso(isoH->GetBinContent(xbin,ybin));
	    float sfidUnc(lepEffH_["m_sel"]->GetBinError(xbin,ybin)), sfisoUnc(isoH->GetBinError(xbin,ybin));
	    float sf(sfid*sfiso), sfUnc(sqrt(pow(sfid*sfisoUnc,2)+pow(sfidUnc*sfiso,2)));
	    lepEffH_["m_sel"]->SetBinContent(xbin,ybin,sf);
	    lepEffH_["m_sel"]->SetBinError(xbin,ybin,sfUnc);
	  }
      fIn->Close();
      
                            //SingleElectron_TriggerSF_Run2016BCDEF_v1.root
                            //v1 has Ele25_eta2p1 and Ele_27
                            //v2 has Ele32_eta2p1
      lepEffUrl=era+"/SingleElectron_TriggerSF_Run2016_Ele32_eta2p1.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffH_["e_singleleptrig"]=(TH2 *)fIn->Get("SF")->Clone();
      //lepEffH_["e_singleleptrig"]=(TH2 *)fIn->Get("Ele27_WPTight_Gsf")->Clone();
      lepEffH_["e_singleleptrig"]->SetDirectory(0);
      fIn->Close();
                    //ElectronIdTight_egammaEffi_Moriond17_EGM2D.root
      lepEffUrl=era+"/ElectronIdTight_egammaEffi_Moriond17_EGM2D.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffH_["e_sel"]=(TH2 *)fIn->Get("EGamma_SF2D")->Clone();
      lepEffH_["e_sel"]->SetDirectory(0);     
      fIn->Close();
                    //ElectronReco_egammaEffi_Moriond17_EGM2D.root
      lepEffUrl=era+"/ElectronReco_egammaEffi_Moriond17_EGM2D.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffH_["e_reco"]=(TH2 *)fIn->Get("EGamma_SF2D");
      lepEffH_["e_reco"]->SetDirectory(0);
      fIn->Close();

      //Load ee SFs from text file
      TString eeHighUrl=era+"/HLT_DoubleEleLegHigPt.txt";
      gSystem->ExpandPathName(eeHighUrl);
      std::ifstream ifs(eeHighUrl);
      float eta1,eta2,pt1,pt2,sf,sfu;
      while(ifs >> eta1 >> eta2 >> pt1 >> pt2 >> sf >> sfu) {
        std::vector<float> sflist {eta1,eta2,pt1,pt2,sf,sfu};
        eeSFHigh.push_back(sflist);
      }
      TString eeLowUrl=era+"/HLT_DoubleEleLegLowPt.txt";
      gSystem->ExpandPathName(eeLowUrl);
      ifs = ifstream(eeLowUrl);
      while(ifs >> eta1 >> eta2 >> pt1 >> pt2 >> sf >> sfu) {
        std::vector<float> sflist {eta1,eta2,pt1,pt2,sf,sfu};
        eeSFLow.push_back(sflist);
      }
    }
}

//
EffCorrection_t LeptonEfficiencyWrapper::getTriggerCorrection(Leptons leptons)
//EffCorrection_t LeptonEfficiencyWrapper::getTriggerCorrection(std::vector<int> pdgId,std::vector<TLorentzVector> leptons)
{
  EffCorrection_t corr(1.0,0.0);
  if(leptons.size()<1) return corr;
  //if(pdgId.size()<1 || leptons.size()<1) return corr;
  if(era_==2015)
    {
      if(leptons.size()>=2)
	{
	  int cat=abs(leptons[0].getPdgId()*leptons[1].getPdgId());
	  //int cat=abs(pdgId[0]*pdgId[1]);
	  if(cat==13*13)      { corr.first=0.894; corr.second=0.894; }
	  else if(cat==11*13) { corr.first=0.931; corr.second=0.931; }
	  else if(cat==11*11) { corr.first=0.930; corr.second=0.930; }
	}
      else
	{
	  TString hname(abs(leptons[0].getPdgId())==11 ? "e" : "m");
	  //TString hname(abs(pdgId[0])==11 ? "e" : "m");
	  //Correct for electron
	  //  pT vs. eta (no abseta_pt)
	  hname += "_singleleptrig";

	  TH2 *h=lepEffH_[hname];
	  float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
	  float etaForEff=TMath::Max(TMath::Min(float(fabs(leptons[0].Eta())),maxEtaForEff),minEtaForEff);
	  Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);
	  
	  float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
	  float ptForEff=TMath::Max(TMath::Min(float(leptons[0].Pt()),maxPtForEff),minPtForEff);
	  Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);
	  
	  corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
	  corr.second=h->GetBinContent(etaBinForEff,ptBinForEff);
	}
    }
  else
    {
      //cf .https://indico.cern.ch/event/532751/contributions/2170250/subcontributions/196785/attachments/1287403/1915605/Talk_TSG_0806.pdf
      //**************************************
      //*********** HARD CODED SFs ***********
      //**************************************
      if(leptons.size()>=2)
	{
	  float leadPt(TMath::Max(leptons[0].Pt(),leptons[1].Pt())), trailPt(TMath::Min(leptons[0].Pt(),leptons[1].Pt()));
          float leadEta(leadPt==leptons[0].Pt() ? leptons[0].Eta() : leptons[1].Eta()), trailEta(trailPt==leptons[0].Pt() ? leptons[0].Eta() : leptons[1].Eta());
	  int cat=abs(leptons[0].getPdgId()*leptons[1].getPdgId());
	  //int cat=abs(pdgId[0]*pdgId[1]);
	  if(cat==13*13)      
	    { 
	      if(leadPt<40)
		{
		  if(trailPt<40) { corr.first=0.908; corr.second=0.039; }
                  //else           { corr.first=1.000; corr.second=0.000; }
		}
	      else if(leadPt<70)
		{
		  if(trailPt<40)      { corr.first=0.897; corr.second=0.033; }
		  else if(trailPt<70) { corr.first=0.859; corr.second=0.054; }
                  //else                { corr.first=1.000; corr.second=0.000; }
		}
	      else 
		{
		  if(trailPt<40)      { corr.first=0.851; corr.second=0.057; }
		  else if(trailPt<70) { corr.first=0.870; corr.second=0.050; }
		  else                { corr.first=0.769; corr.second=0.174; }
		}
	    }
	  else if(cat==11*13)
	    { 
	      if(abs(leptons[0].getPdgId())==11) { leadPt=leptons[0].Pt(); trailPt=leptons[1].Pt(); }
	      else                  { leadPt=leptons[1].Pt(); trailPt=leptons[0].Pt(); }
	      if(leadPt<40)
                {
                  if(trailPt<40)      { corr.first=0.868; corr.second=0.032; }
                  else if(trailPt<70) { corr.first=0.874; corr.second=0.036; }
                  else                { corr.first=0.933; corr.second=0.060; }
		}
	      else if(leadPt<70)
		{
                  if(trailPt<40)      { corr.first=0.863; corr.second=0.032; }
                  else if(trailPt<70) { corr.first=0.902; corr.second=0.034; }
                  else                { corr.first=0.768; corr.second=0.070; }
                }
	      else
		{
                  if(trailPt<40)      { corr.first=0.901; corr.second=0.042; }
                  else if(trailPt<70) { corr.first=0.838; corr.second=0.052; }
                  else                { corr.first=1.000; corr.second=0.045; }
                }
	    }
          //else if(cat==11*11)         { corr.first=1.00;  corr.second=0.000; } //Using 1 until correct HLT SFs are available
          /*
	  else if(cat==11*11) 
	    {
	      if(leadPt<40)
                {
                  if(trailPt<40) { corr.first=0.853; corr.second=0.056; }
                  //else           { corr.first=1.000; corr.second=0.000; }
		}
	      else if(leadPt<70)
		{
                  if(trailPt<40)      { corr.first=0.954; corr.second=0.026; }
                  else if(trailPt<70) { corr.first=0.944; corr.second=0.042; }
                  //else                { corr.first=1.000; corr.second=0.000; }
                }
	      else
		{
                  if(trailPt<40)      { corr.first=0.918; corr.second=0.046; }
                  else if(trailPt<70) { corr.first=0.984; corr.second=0.037; }
                  else                { corr.first=0.941; corr.second=0.122; }
                }
	    }
          */
          else if(cat==11*11) {
            for(auto &it : eeSFHigh) {
              //[0] is low eta, [1] is high eta
              if(leadEta < it[0] || leadEta > it[1]) continue;
              //[2] is low pT, [3] is high pT
              if(leadPt < it[2] || leadPt > it[3]) continue;
              //[4] is SF, [5] is SF uncertainty
              corr.first=it[4]; corr.second=it[5];
              if(debug_) std::cout << "corr = " << it[4] << std::endl;
              break;
            }
            for(auto &it : eeSFLow) {
              //[0] is low eta, [1] is high eta
              if(trailEta < it[0] || trailEta > it[1]) continue;
              //[2] is low pT, [3] is high pT
              if(trailPt < it[2] || trailPt > it[3]) continue;
              //[4] is SF, [5] is SF uncertainty
              corr.first*=it[4]; corr.second=TMath::Sqrt(corr.second*corr.second + it[5]*it[5]);
              if(debug_) std::cout << "corr = " << it[4] << std::endl;
              break;
            }
          }
	}
      else 
	{
	  TString hname(abs(leptons[0].getPdgId())==11 ? "e" : "m");
	  hname += "_singleleptrig";
	  if( abs(leptons[0].getPdgId())==13 && lepEffH_.find(hname)!=lepEffH_.end() )
	    {
              if(debug_) std::cout << hname << std::endl;
	      TH1 *h=lepEffH_[hname];
	      float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
	      float etaForEff=TMath::Max(TMath::Min(float(fabs(leptons[0].Eta())),maxEtaForEff),minEtaForEff);
	      Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);

	      float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
	      float ptForEff=TMath::Max(TMath::Min(float(leptons[0].Pt()),maxPtForEff),minPtForEff);
	      Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);

	      corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
	      corr.second=h->GetBinError(etaBinForEff,etaBinForEff);
	    }
          //electron histogram has inverted axes and uses eta, not abs(eta)
          //electron histogram uses eta, not abs(eta)
	  else if( abs(leptons[0].getPdgId())==11 && lepEffH_.find(hname)!=lepEffH_.end() ) {
              if(debug_) std::cout << hname << std::endl;
	      TH1 *h=lepEffH_[hname];
              float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(leptons[0].Eta()),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);

              float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(leptons[0].Pt()),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);

	      corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
	      corr.second=h->GetBinError(etaBinForEff,etaBinForEff);
              if(debug_) std::cout << std::endl
                         << minEtaForEff << " " << maxEtaForEff << " " << etaForEff << " "
                         << minPtForEff << " " << maxPtForEff << " " << ptForEff << endl
                         << corr.first << " " << corr.second << endl;
          }
	}
    }

  return corr;
}

//
EffCorrection_t LeptonEfficiencyWrapper::getOfflineCorrection(Particle lep, int nvtx)
//EffCorrection_t LeptonEfficiencyWrapper::getOfflineCorrection(Particle lep, TString runPeriod)
//EffCorrection_t LeptonEfficiencyWrapper::getOfflineCorrection(int pdgId,float pt,float eta, TString runPeriod)
{
  int pdgId = lep.getPdgId();
  float pt(lep.Pt()), eta(lep.Eta());
  EffCorrection_t corr(1.0,0.0);

  //update correction from histo, if found
  TString idstr(abs(pdgId)==11 ? "e" : "m");
  TString hname(idstr);
  hname+="_sel";
  if( lepEffH_.find(hname)!=lepEffH_.end() )
    {
      TH2 *h=lepEffH_[hname];
      float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
      float etaForEff;
      if (minEtaForEff >= 0.) //axis is abseta
        etaForEff=TMath::Max(TMath::Min(float(fabs(eta)),maxEtaForEff),minEtaForEff);
      else //axis is signed eta
        etaForEff=TMath::Max(TMath::Min(float(eta),maxEtaForEff),minEtaForEff);
      Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);

      float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
      float ptForEff=TMath::Max(TMath::Min(pt,maxPtForEff),minPtForEff);
      Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);
      
      corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
      corr.second=h->GetBinError(etaBinForEff,ptBinForEff);

      //tracking efficiency (if available)
      hname=idstr+"_tk_aeta_"+runPeriod_;
      if(lepEffGr_.find(hname)!=lepEffGr_.end() && 0) //No muon tracking SF per POG recommendation https://hypernews.cern.ch/HyperNews/CMS/get/top/2671.html
        {
          if(debug_) std::cout << hname << std::endl;
          Double_t x(0.),xdiff(9999.),y(0.);
          float tkEffSF(1.0),tkEffSFUnc(0);
          for(Int_t ip=0; ip<lepEffGr_[hname]->GetN(); ip++)
            {
              lepEffGr_[hname]->GetPoint(ip,x,y);
              float ixdiff(TMath::Abs(fabs(eta)-x));
              if(ixdiff>xdiff) continue;
              xdiff=ixdiff;
              tkEffSF=y;
              tkEffSFUnc=lepEffGr_[hname]->GetErrorY(ip);
            }
          corr.second = sqrt(pow(tkEffSFUnc*corr.first,2)+pow(tkEffSF*corr.second,2));
          if(debug_) std::cout << "tk eff= " << tkEffSF << std::endl;
          corr.first  = corr.first*tkEffSF;
        }
      /*
      hname=idstr+"_tk_vtx_"+runPeriod_;
      if(lepEffGr_.find(hname)!=lepEffGr_.end())
        {
          if(debug_) std::cout << hname << std::endl;
          Double_t x(0.),xdiff(9999.),y(0.);
          float tkEffSF(1.0),tkEffSFUnc(0);
          for(Int_t ip=0; ip<lepEffGr_[hname]->GetN(); ip++)
            {
              lepEffGr_[hname]->GetPoint(ip,x,y);
              float ixdiff(TMath::Abs(fabs(nvtx)-x));
              if(ixdiff>xdiff) continue;
              xdiff=ixdiff;
              tkEffSF=y;
              tkEffSFUnc=lepEffGr_[hname]->GetErrorY(ip);
            }
          corr.second = sqrt(pow(tkEffSFUnc*corr.first,2)+pow(tkEffSF*corr.second,2));
          if(debug_) std::cout << "tk eff= " << tkEffSF << std::endl;
          corr.first  = corr.first*tkEffSF;
        }
      */
 
      //reco efficiency (if available)
      hname=idstr+"_reco";
      if(lepEffH_.find(hname)!=lepEffH_.end() )
	{
	  TH2 *h=lepEffH_[hname];
	  float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
	  float etaForEff=TMath::Max(TMath::Min(float(fabs(eta)),maxEtaForEff),minEtaForEff);
	  Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);
	  
	  float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
	  float ptForEff=TMath::Max(TMath::Min(pt,maxPtForEff),minPtForEff);
	  Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);
	  
	  corr.second = sqrt(pow(h->GetBinError(etaBinForEff,ptBinForEff)*corr.first,2)+pow(h->GetBinError(etaBinForEff,ptBinForEff)*corr.second,2));
	  corr.first  = corr.first*h->GetBinContent(etaBinForEff,ptBinForEff);
	  
	}
    }

  return corr;

}

EffCorrection_t LeptonEfficiencyWrapper::getOfflineCorrection(pfTrack lep)
{
  int pdgId = lep.getPdgId();
  float pt(lep.Pt()), eta(lep.Eta());
  EffCorrection_t corr(1.0,0.0);

  //update correction from histo, if found
  TString idstr(abs(pdgId)==11 ? "e" : "m");
  TString hname(idstr);
  hname+="_sel";
  if( lepEffH_.find(hname)!=lepEffH_.end() )
    {
      /*
      //for ID only
      TH2 *h=lepEffH_[hname];
      float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
      float etaForEff;
      if (minEtaForEff >= 0.) //axis is abseta
        etaForEff=TMath::Max(TMath::Min(float(fabs(eta)),maxEtaForEff),minEtaForEff);
      else //axis is signed eta
        etaForEff=TMath::Max(TMath::Min(float(eta),maxEtaForEff),minEtaForEff);
      Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);

      float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
      float ptForEff=TMath::Max(TMath::Min(pt,maxPtForEff),minPtForEff);
      Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);
      
      corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
      corr.second=h->GetBinError(etaBinForEff,ptBinForEff);
      */

      //tracking efficiency (if available)
      hname=idstr+"_tk_aeta_"+runPeriod_;
      if(lepEffGr_.find(hname)!=lepEffGr_.end() && 0) //No muon tracking SF per POG recommendation https://hypernews.cern.ch/HyperNews/CMS/get/top/2671.html
        {
          if(debug_) std::cout << hname << std::endl;
          Double_t x(0.),xdiff(9999.),y(0.);
          float tkEffSF(1.0),tkEffSFUnc(0);
          for(Int_t ip=0; ip<lepEffGr_[hname]->GetN(); ip++)
            {
              lepEffGr_[hname]->GetPoint(ip,x,y);
              float ixdiff(TMath::Abs(fabs(eta)-x));
              if(ixdiff>xdiff) continue;
              xdiff=ixdiff;
              tkEffSF=y;
              tkEffSFUnc=lepEffGr_[hname]->GetErrorY(ip);
            }
          corr.second = sqrt(pow(tkEffSFUnc*corr.first,2)+pow(tkEffSF*corr.second,2));
          if(debug_) std::cout << "tk eff= " << tkEffSF << std::endl;
          corr.first  = corr.first*tkEffSF;
        }
      /*
      hname=idstr+"_tk_vtx_"+runPeriod_;
      if(lepEffGr_.find(hname)!=lepEffGr_.end())
        {
          if(debug_) std::cout << hname << std::endl;
          Double_t x(0.),xdiff(9999.),y(0.);
          float tkEffSF(1.0),tkEffSFUnc(0);
          for(Int_t ip=0; ip<lepEffGr_[hname]->GetN(); ip++)
            {
              lepEffGr_[hname]->GetPoint(ip,x,y);
              float ixdiff(TMath::Abs(fabs(nvtx)-x));
              if(ixdiff>xdiff) continue;
              xdiff=ixdiff;
              tkEffSF=y;
              tkEffSFUnc=lepEffGr_[hname]->GetErrorY(ip);
            }
          corr.second = sqrt(pow(tkEffSFUnc*corr.first,2)+pow(tkEffSF*corr.second,2));
          if(debug_) std::cout << "tk eff= " << tkEffSF << std::endl;
          corr.first  = corr.first*tkEffSF;
        }
 
      //reco efficiency (if available)
      hname=idstr+"_reco";
      if(lepEffH_.find(hname)!=lepEffH_.end() )
	{
	  TH2 *h=lepEffH_[hname];
	  float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
	  float etaForEff=TMath::Max(TMath::Min(float(fabs(eta)),maxEtaForEff),minEtaForEff);
	  Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);
	  
	  float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
	  float ptForEff=TMath::Max(TMath::Min(pt,maxPtForEff),minPtForEff);
	  Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);
	  
	  corr.second = sqrt(pow(h->GetBinError(etaBinForEff,ptBinForEff)*corr.first,2)+pow(h->GetBinError(etaBinForEff,ptBinForEff)*corr.second,2));
	  corr.first  = corr.first*h->GetBinContent(etaBinForEff,ptBinForEff);
	  
	}
      */
    }

  return corr;

}

LeptonEfficiencyWrapper::~LeptonEfficiencyWrapper()
{
  //for(auto& it : lepEffH_) it.second->Delete();
}
