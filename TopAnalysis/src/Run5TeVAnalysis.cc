#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"

#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"
#include "TopLJets2015/TopAnalysis/interface/Run5TeVAnalysis.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"

#include <string>
#include <vector>

using namespace std;

//
void Run5TeVAnalysis(TString inFileName,
		     TString outFileName,
		     Int_t channelSelection,
		     Int_t chargeSelection,
                     FlavourSplitting flavourSplitting,
                     TH1F *normH,
                     Bool_t runSysts,
		     TString era)
{
  if(inFileName=="") 
    {
      std::cout << "No inputs specified. return" << std::endl;
      return;
    }

  bool isMC(false);
  if(inFileName.Contains("/MC")) isMC=true;
  std::cout << "Will process " << inFileName << " and save the results in " << outFileName << endl
	    << "Sample will be treated as MC=" << isMC <<  std::endl
	    << "Corrections to be retrieved from era=" << era << std::endl;

  std::map<TString, TGraphAsymmErrors *> expBtagEff;
  BTagSFUtil myBTagSFUtil;
  if(isMC)
    { 
      TString btagEffExpUrl(era+"/expTageff.root");
      TFile *beffIn=TFile::Open(btagEffExpUrl);
      expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEff["others"]=(TGraphAsymmErrors *)beffIn->Get("others");
      beffIn->Close();
    }
 
  std::vector<TString>* inFileNames_p = new std::vector<TString>;
  inFileNames_p->push_back(inFileName);
  const Int_t nFiles = (Int_t)inFileNames_p->size();

  //book histograms
  std::map<TString,TH1 *> histos;
  for(int ij=1; ij<=4; ij++)
    {
      TString pf(Form("%dj",ij));
      histos["lpt_"+pf]  = new TH1F("lpt_"+pf,";Transverse momentum [GeV];Events",20.,0.,200.);
      histos["leta_"+pf] = new TH1F("leta_"+pf,";Pseudo-rapidity [GeV];Events",20.,0.,2.1);
      histos["ht_"+pf]   = new TH1F("ht_"+pf,";H_{T} [GeV];Events",20.,0.,300.);
    }
  histos["njnb"] = new TH1F("njnb",";Category;Events" ,11,0.,11.);
  histos["njnb"]->GetXaxis()->SetBinLabel(1,"1j,0b");
  histos["njnb"]->GetXaxis()->SetBinLabel(2,"1j,1b");
  histos["njnb"]->GetXaxis()->SetBinLabel(3,"2j,0b");
  histos["njnb"]->GetXaxis()->SetBinLabel(4,"2j,1b");
  histos["njnb"]->GetXaxis()->SetBinLabel(5,"2j,2b");
  histos["njnb"]->GetXaxis()->SetBinLabel(6,"3j,0b");
  histos["njnb"]->GetXaxis()->SetBinLabel(7,"3j,1b");
  histos["njnb"]->GetXaxis()->SetBinLabel(8,"3j,2b");
  histos["njnb"]->GetXaxis()->SetBinLabel(9,"4j,0b");
  histos["njnb"]->GetXaxis()->SetBinLabel(10,"4j,1b");
  histos["njnb"]->GetXaxis()->SetBinLabel(11,"4j,2b");
  
  Int_t nSysts(0);
  if(isMC)
    {
      TString expSysts[]={"btagup","btagdn","othertagup","othertagdn","jesup","jesdn","jerup","jerdn","leffup","leffdn"};
      nSysts=sizeof(expSysts)/sizeof(TString);
      histos["njnbshapes_exp"]=new TH2F("njnbshapes_exp",";Category;Systematic uncertainty;Events",11,0,11,nSysts,0,nSysts);
      for(int i=0; i<nSysts; i++)
	histos["njnb_exp"]->GetYaxis()->SetBinLabel(i+1,expSysts[i]);
      histos["njnbshapes_gen"]=new TH2F("njnbshapes_gen",";Category;Systematic uncertainty;Events",11,0,11,500,0,500);
      for(int i=0; i<500;i++)
	histos["njnb_exp"]->GetYaxis()->SetBinLabel(i+1,Form("genUnc%d",i));
      for(int xbin=1; xbin<=histos["njnb"]->GetXaxis()->GetNbins(); xbin++)
	{
	  histos["njnbshapes_exp"]->GetXaxis()->SetBinLabel(xbin,histos["njnb"]->GetXaxis()->GetBinLabel(xbin));
	  histos["njnbshapes_gen"]->GetXaxis()->SetBinLabel(xbin,histos["njnb"]->GetXaxis()->GetBinLabel(xbin));
	}
    }
      
  //prepare histograms
  for(std::map<TString,TH1 *>::iterator it=histos.begin();
      it!=histos.end();
      it++)
    {
      it->second->Sumw2();
      it->second->SetDirectory(0);
    }


  Float_t totalEvtNorm(0);
  for(Int_t fileIter = 0; fileIter < nFiles; fileIter++){

    TString inF(inFileNames_p->at(fileIter));
    if(inF.Contains("/store") && !inF.Contains("root:")) inF="root://eoscms//eos/cms/"+inF;
    TFile* inFile_p = TFile::Open(inF, "READ");
    
    TTree* lepTree_p = (TTree*)inFile_p->Get("ggHiNtuplizer/EventTree");
    TTree* jetTree_p = (TTree*)inFile_p->Get("ak4PFJetAnalyzer/t");
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* hltTree_p = (TTree*)inFile_p->Get("hltanalysis/HltTree");

    //muon variables
    std::vector<float>* muPt_p = 0;
    std::vector<float>* muPhi_p = 0;
    std::vector<float>* muEta_p = 0;
    std::vector<int>* muChg_p = 0;
    std::vector<float>* muChi2NDF_p = 0;
    std::vector<float>* muInnerD0_p = 0;
    std::vector<float>* muInnerDz_p = 0;
    std::vector<int>* muMuonHits_p = 0;
    std::vector<int>* muStations_p = 0;
    std::vector<int>* muTrkLayers_p = 0;
    std::vector<int>* muPixelHits_p = 0;    
    std::vector<float> *muPFChIso_p=0,*muPFPhoIso_p=0,*muPFNeuIso_p=0,*muPFPUIso_p=0;
    lepTree_p->SetBranchStatus("*", 0);
    lepTree_p->SetBranchStatus("muPt", 1);
    lepTree_p->SetBranchStatus("muPhi", 1);
    lepTree_p->SetBranchStatus("muEta", 1);
    lepTree_p->SetBranchStatus("muCharge", 1);
    lepTree_p->SetBranchStatus("muChi2NDF", 1);
    lepTree_p->SetBranchStatus("muInnerD0", 1);
    lepTree_p->SetBranchStatus("muInnerDz", 1);
    lepTree_p->SetBranchStatus("muMuonHits", 1);
    lepTree_p->SetBranchStatus("muStations", 1);
    lepTree_p->SetBranchStatus("muTrkLayers", 1);
    lepTree_p->SetBranchStatus("muPixelHits", 1);    
    lepTree_p->SetBranchStatus("muPFChIso", 1);
    lepTree_p->SetBranchStatus("muPFPhoIso", 1);
    lepTree_p->SetBranchStatus("muPFNeuIso", 1);
    lepTree_p->SetBranchStatus("muPFPUIso", 1);
    lepTree_p->SetBranchAddress("muPt", &muPt_p);
    lepTree_p->SetBranchAddress("muPhi", &muPhi_p);
    lepTree_p->SetBranchAddress("muEta", &muEta_p);
    lepTree_p->SetBranchAddress("muCharge", &muChg_p);
    lepTree_p->SetBranchAddress("muChi2NDF", &muChi2NDF_p);
    lepTree_p->SetBranchAddress("muInnerD0", &muInnerD0_p);
    lepTree_p->SetBranchAddress("muInnerDz", &muInnerDz_p);
    lepTree_p->SetBranchAddress("muMuonHits", &muMuonHits_p);
    lepTree_p->SetBranchAddress("muStations", &muStations_p);
    lepTree_p->SetBranchAddress("muTrkLayers", &muTrkLayers_p);
    lepTree_p->SetBranchAddress("muPixelHits", &muPixelHits_p);    
    lepTree_p->SetBranchAddress("muPFChIso", &muPFChIso_p);
    lepTree_p->SetBranchAddress("muPFPhoIso", &muPFPhoIso_p);
    lepTree_p->SetBranchAddress("muPFNeuIso", &muPFNeuIso_p);
    lepTree_p->SetBranchAddress("muPFPUIso", &muPFPUIso_p);

    //electron variables
    std::vector<float>* elePt_p = 0;
    std::vector<float>* elePhi_p = 0;
    std::vector<float>* eleEta_p = 0;
    std::vector<int>*   eleIDVeto_p=0;
    lepTree_p->SetBranchStatus("elePt", 1);
    lepTree_p->SetBranchStatus("elePhi", 1);
    lepTree_p->SetBranchStatus("eleEta", 1);
    lepTree_p->SetBranchStatus("eleIDVeto", 1);
    lepTree_p->SetBranchAddress("elePt", &elePt_p);
    lepTree_p->SetBranchAddress("elePhi", &elePhi_p);
    lepTree_p->SetBranchAddress("eleEta", &eleEta_p);
    lepTree_p->SetBranchAddress("eleIDVeto", &eleIDVeto_p);
    
    //jet variables
    const int maxJets = 5000;
    Int_t   nref;
    Float_t jtpt[maxJets]; 
    Float_t jteta[maxJets];
    Float_t jtphi[maxJets];
    Float_t jtm[maxJets]; 
    Float_t discr_csvV2[maxJets];
    Float_t refpt[maxJets];
    Float_t refparton_flavor[maxJets];
    jetTree_p->SetBranchStatus("*", 0);
    jetTree_p->SetBranchStatus("nref", 1);
    jetTree_p->SetBranchStatus("jtpt", 1);
    jetTree_p->SetBranchStatus("jtphi", 1);
    jetTree_p->SetBranchStatus("jteta", 1);
    jetTree_p->SetBranchStatus("jtm", 1);
    jetTree_p->SetBranchStatus("discr_csvV2", 1);
    jetTree_p->SetBranchStatus("refpt", 1);
    jetTree_p->SetBranchStatus("refparton_flavor", 1);
    jetTree_p->SetBranchAddress("nref", &nref);
    jetTree_p->SetBranchAddress("jtpt", jtpt);
    jetTree_p->SetBranchAddress("jtphi", jtphi);
    jetTree_p->SetBranchAddress("jteta", jteta);
    jetTree_p->SetBranchAddress("jtm", jtm);
    jetTree_p->SetBranchAddress("discr_csvV2", discr_csvV2);
    jetTree_p->SetBranchAddress("refpt", refpt);
    jetTree_p->SetBranchAddress("refparton_flavorForB", refparton_flavor);

    //event variables
    UInt_t run_, lumi_;
    ULong64_t evt_;
    Int_t hiBin_;
    Float_t vz_;
    std::vector<float> *ttbar_w_p=0;
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("ttbar_w",1);
    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("evt", &evt_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("ttbar_w",&ttbar_w_p);

    //trigger
    int trig = 0;
    std::string triggerName;
    triggerName = isMC ? "HLT_HIL2Mu15ForPPRef_v1" : "HLT_HIL2Mu15_v1"; 
    hltTree_p->SetBranchStatus(triggerName.data(),1);
    hltTree_p->SetBranchAddress(triggerName.data(),&trig);
    
    Int_t nEntries = (Int_t)lepTree_p->GetEntries();
    Int_t entryDiv = ((Int_t)(nEntries/1000));
    
    std::cout << "Analysing " << nEntries << " events isMC=" << isMC
	      << " trigger=" << triggerName << endl;
    
    for(Int_t entry = 0; entry < nEntries; entry++)
      {
	if(entry%entryDiv == 0 && nEntries >= 10000) 
	  std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
      
	//readout this event
	lepTree_p->GetEntry(entry);
	jetTree_p->GetEntry(entry);
	hiTree_p->GetEntry(entry);
	hltTree_p->GetEntry(entry);

	//assign an event weight
	float evWeight(1.0);
	if(isMC && ttbar_w_p->size()) evWeight=ttbar_w_p->at(0);
	totalEvtNorm+=evWeight;

	//select good muons
	//cf. details in https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
	std::vector<TLorentzVector> tightMuons,looseMuons,tightMuonsNonIso;
	std::vector<int> muonCharge;
	const Int_t nMu = (Int_t)muPt_p->size();
	for(Int_t muIter = 0; muIter < nMu; muIter++)
	  {
	    bool passLooseKin( muPt_p->at(muIter) > 15. && TMath::Abs(muEta_p->at(muIter))<2.4);
	    bool passTightKin( muPt_p->at(muIter) > 30. && TMath::Abs(muEta_p->at(muIter))<2.1);
	    bool passLooseId(true);
	    bool passTightId( passLooseId
			      && muChi2NDF_p->at(muIter) < 10
			      && muMuonHits_p->at(muIter) >0 		
			      && muStations_p->at(muIter) >1
			      && TMath::Abs(muInnerD0_p->at(muIter))<0.2
			      && TMath::Abs(muInnerDz_p->at(muIter))<0.5
			      && muPixelHits_p->at(muIter)>0		
			    && muTrkLayers_p->at(muIter)>5);
	    float relIso=(muPFChIso_p->at(muIter)+TMath::Max(muPFPhoIso_p->at(muIter)+muPFNeuIso_p->at(muIter)-0.5*muPFPUIso_p->at(muIter),0.))/muPt_p->at(muIter);
	    bool passTightIso( relIso<0.15);
	    bool passLooseIso( relIso<0.2);
	    
	    //save muon if good
	    TLorentzVector p4(0,0,0,0);
	    p4.SetPtEtaPhiM(muPt_p->at(muIter),muEta_p->at(muIter),muPhi_p->at(muIter), 0.1056583715);
	    if(passTightKin && passTightId && !passTightIso && relIso>0.2) 
	      {
		tightMuonsNonIso.push_back(p4);
		muonCharge.push_back(muChg_p->at(muIter));
	      }
	    if(passTightKin && passTightId && passTightIso)     
	      {
		tightMuons.push_back( p4 );
		muonCharge.push_back(muChg_p->at(muIter));
	      }
	    else if(passLooseKin && passLooseId && passLooseIso) looseMuons.push_back( p4 );
	  }

	//select the muon
	if(channelSelection==1300)
	  {
	    if(tightMuons.size()!=0 || looseMuons.size()!=0) continue;
	    if(tightMuonsNonIso.size()!=1) continue;
	    tightMuons=tightMuonsNonIso;
	  }
	if(channelSelection==13)
	  {
	    if(tightMuons.size()!=1) continue; //=1 tight muon
	    if(looseMuons.size()) continue;    //no extra muons
	    if(chargeSelection!=0)
	      {
		if(muonCharge[0]!=chargeSelection) continue;
	      }
	  }

	//check for extra electrons in the event (veto id is used)
	//see https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns (electron id)
	Int_t nLooseEle(0);
	const Int_t nEle = (Int_t)elePt_p->size();
	for(Int_t eleIter = 0; eleIter < nEle; eleIter++)
	  {
	    
	    if(TMath::Abs(eleEta_p->at(eleIter)) < 20)  continue;
	    if(TMath::Abs(eleEta_p->at(eleIter)) > 2.5) continue;
	    if(eleIDVeto_p->at(eleIter)==0) continue;
	    nLooseEle++;
	  }
	if(nLooseEle>0) continue;
	
	//jet counting
	Int_t njets(0),nBtags(0),htsum(0);
	std::vector<Int_t> njetsVar(4,0),nBtagsVar(4,0);
	for(Int_t jetIter = 0; jetIter < nref; jetIter++)
	  {
	    //cross clean with trigger muon
	    TLorentzVector jp4(0,0,0,0);
	    jp4.SetPtEtaPhiM(jtpt[jetIter],jteta[jetIter],jtphi[jetIter],jtm[jetIter]);
	    if(jp4.DeltaR(tightMuons[0])<0.4) continue;

	    //in tracker region
	    if(TMath::Abs(jp4.Eta())>2.4) continue;

	    if(isMC)
	      {
		std::vector<float> jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),refpt[jetIter]);
		TLorentzVector rawjp4(jp4);
		jp4 *= jerSmear[0];

		//JES varied selections
		if(jp4.Pt()>30./(1+0.028)) njetsVar[0]++;
		if(jp4.Pt()>30./(1-0.028)) njetsVar[1]++;

		//JER variations
		if(rawjp4.Pt()>30./jerSmear[1]) njetsVar[2]++;
		if(rawjp4.Pt()>30./jerSmear[2]) njetsVar[3]++;
	      }
	    
	    //nominal selection
	    if(jp4.Pt()<30) continue;
	    ++njets;
	    htsum+=jp4.Pt();

	    bool passCSVM(discr_csvV2[jetIter]>0.8);

	    if(passCSVM) ++nBtags;
	    
	    //b-tagging variations
	    if(isMC)
	      {
		float jptforBtag(jp4.Pt()>1000. ? 999. : jp4.Pt());
		if(abs(refparton_flavor[jetIter])==5)
		  {
		    float expEff    = expBtagEff["b"]->Eval(jptforBtag); 
		    bool passCSVMUp(passCSVM);
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMUp,1.1,expEff);	
		    nBtagsVar[0]+=passCSVMUp;

		    bool passCSVMDn(passCSVM);
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMDn,0.9,expEff);	
		    nBtagsVar[1]+=passCSVMDn;

		    nBtagsVar[2]+=passCSVM;
		    nBtagsVar[3]+=passCSVM;
		  }
		else
		  {
		    nBtagsVar[0]+=passCSVM;
		    nBtagsVar[1]+=passCSVM;
		    
		    float expEff    = expBtagEff["others"]->Eval(jptforBtag); 
		    bool passCSVMUp(passCSVM);
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMUp,1.3,expEff);	
		    nBtagsVar[2]+=passCSVMUp;
		    
		    bool passCSVMDn(passCSVM);
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMDn,0.7,expEff);	
		    nBtagsVar[3]+=passCSVMDn;
		  }		
	      }
	  }
      
	//control plots for the nominal distribution
	if(njets>0)
	  {
	    TString pf(Form("%dj",TMath::Min(4,njets)));
	    histos["lpt_"+pf]->Fill(tightMuons[0].Pt(),evWeight);
	    histos["leta_"+pf]->Fill(fabs(tightMuons[0].Eta()),evWeight);
	    histos["ht_"+pf]->Fill(htsum,evWeight);
	   
	    Int_t binToFill(0);
	    if(njets==1 && nBtags ==1) binToFill=1;
	    if(njets==2 && nBtags ==0) binToFill=2;
	    if(njets==2 && nBtags ==1) binToFill=3;
	    if(njets==2 && nBtags >=2) binToFill=4;
	    if(njets==3 && nBtags ==0) binToFill=5;
	    if(njets>=3 && nBtags ==1) binToFill=6;
	    if(njets>=3 && nBtags >=2) binToFill=7;
	    if(njets>=4 && nBtags ==0) binToFill=8;
	    if(njets>=4 && nBtags ==1) binToFill=9;
	    if(njets>=4 && nBtags >=2) binToFill=10;
	    histos["njnb"]->Fill(binToFill,evWeight);      
	    
	    if(!isMC || !runSysts) continue;

	    //theory uncertainties (by matrix-element weighting)
	    for(size_t igs=0; igs<ttbar_w_p->size(); igs++)
	      {
		float newWeight( ttbar_w_p->at(igs) );
		((TH2 *)histos["njnb_gen"])->Fill(binToFill,igs,newWeight);
	      }

	    //experimental uncertainties
	    for(Int_t ies=0; ies<nSysts; ies++)
	      {
		float newWeight(evWeight);
		Int_t njets_i(njets), nBtags_i(nBtags);
		if(ies==0) nBtags_i=nBtagsVar[0];
		if(ies==1) nBtags_i=nBtagsVar[1];
		if(ies==2) nBtags_i=nBtagsVar[2];
		if(ies==3) nBtags_i=nBtagsVar[3];
		if(ies==4) njets_i=njetsVar[0];
		if(ies==5) njets_i=njetsVar[1];
		if(ies==6) njets_i=njetsVar[2];
		if(ies==7) njets_i=njetsVar[3];
		if(ies==8) newWeight *= 1.03;
		if(ies==9) newWeight *= 0.97;

		//decide which bin to fill
		Int_t binToFill(0);
		if(njets_i==1 && nBtags_i ==1) binToFill=1;
		if(njets_i==2 && nBtags_i ==0) binToFill=2;
		if(njets_i==2 && nBtags_i ==1) binToFill=3;
		if(njets_i==2 && nBtags_i >=2) binToFill=4;
		if(njets_i==3 && nBtags_i ==0) binToFill=5;
		if(njets_i>=3 && nBtags_i ==1) binToFill=6;
		if(njets_i>=3 && nBtags_i >=2) binToFill=7;
		if(njets_i>=4 && nBtags_i ==0) binToFill=8;
		if(njets_i>=4 && nBtags_i ==1) binToFill=9;
		if(njets_i>=4 && nBtags_i >=2) binToFill=10;
		((TH2 *)histos["njnb_exp"])->Fill(binToFill,ies,newWeight);
	      }

	      
	  }
      }
    
    inFile_p->Close();
    delete inFile_p;    
  }
  
  //dump histograms
  TFile* outFile_p = new TFile(outFileName, "RECREATE");
  outFile_p->cd();
  for(std::map<TString, TH1 *>::iterator it=histos.begin();
      it!=histos.end();
      it++)
    {
      if(isMC && totalEvtNorm>0) it->second->Scale(1./totalEvtNorm);
      it->second->SetDirectory(outFile_p);
      it->second->Write();
    }
  outFile_p->Close();
  delete outFile_p;
  
  return;
}
