#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"

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
  bool isTTJets(false);
  if(inFileName.Contains("/MCTTNominal")) isTTJets=true;

  float totalEvtNorm(1.0);
  if(isMC && normH) totalEvtNorm=normH->GetBinContent(1);
  if(!isMC) runSysts=false;
  std::cout << "Will process " << inFileName << " and save the results in " << outFileName << endl
	    << "Sample will be treated as MC=" << isMC <<  std::endl
	    << "Systematics will be run=" << runSysts << std::endl
	    << "Corrections to be retrieved from era=" << era << std::endl
	    << "Total normalization factor=" << totalEvtNorm << std::endl;

  std::map<TString, TGraphAsymmErrors *> expBtagEff;
  BTagSFUtil myBTagSFUtil;
  if(isMC)
    { 
      TString btagEffExpUrl(era+"/expTageff.root");
      gSystem->ExpandPathName(era+"/expTageff.root");   
      TFile *beffIn=TFile::Open(btagEffExpUrl);
      expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();
      cout << "Read " << expBtagEff.size() << " b-tag efficiency expectations from " << btagEffExpUrl << endl;
    }

  std::vector<TString>* inFileNames_p = new std::vector<TString>;
  inFileNames_p->push_back(inFileName);
  const Int_t nFiles = (Int_t)inFileNames_p->size();

  Int_t nSysts(0);
  TString expSysts[]={"btagup","btagdn","othertagup","othertagdn","jesup","jesdn","jerup","jerdn","leffup","leffdn"};
  
  //book histograms
  std::map<TString,TH1 *> histos;
  histos["lpt"]  = new TH1F("lpt",";Transverse momentum [GeV];Events",20.,0.,200.);
  histos["leta"] = new TH1F("leta",";Pseudo-rapidity;Events",20.,0.,2.1);
  histos["mt"]   = new TH1F("mt",";Transverse Mass [GeV];Events" ,20,0.,200.);
  histos["metpt"]= new TH1F("metpt",";Missing transverse energy [GeV];Events" ,20,0.,200.);

  for(int ij=0; ij<=2; ij++)
    {
      TString pf(Form("%db",ij));
      histos["lpt_"+pf]    = new TH1F("lpt_"+pf,";Transverse momentum [GeV];Events",10.,0.,200.);
      histos["leta_"+pf]   = new TH1F("leta_"+pf,";Pseudo-rapidity;Events",10.,0.,2.1);
      histos["jpt_"+pf]    = new TH1F("jpt_"+pf,";Transverse momentum [GeV];Events",10.,0.,200.);
      histos["jeta_"+pf]   = new TH1F("jeta_"+pf,";Pseudo-rapidity;Events",10.,0.,2.5);
      histos["ht_"+pf]     = new TH1F("ht_"+pf,";H_{T} [GeV];Events",10.,0.,800.);
      histos["metpt_"+pf]  = new TH1F("metpt_"+pf,";Missing transverse energy [GeV];Events" ,10,0.,200.);
      histos["metphi_"+pf] = new TH1F("metphi_" + pf,";MET #phi [rad];Events" ,10,-3.2,3.2);
      histos["mt_"+pf]     = new TH1F("mt_"+pf,";Transverse Mass [GeV];Events" ,10,0.,200.);
      histos["mjj_"+pf]    = new TH1F("mjj_"+pf,";Mass(j,j') [GeV];Events" ,16,0.,200.);
      histos["mlb_"+pf]    = new TH1F("mlb_"+pf,";Mass(l,b) [GeV];Events" ,15,0.,250.);
      histos["njets_"+pf]  = new TH1F("njets_"+pf,";Jet multiplicity;Events" ,6,2.,8.);

      if(isMC && runSysts)
	{
	  nSysts=sizeof(expSysts)/sizeof(TString);
	  histos["mjjshapes_"+pf+"_exp"]=new TH2F("mjjshapes_"+pf+"_exp",";Mass(j,j');Systematic uncertainty;Events",16,0,200,nSysts,0,nSysts);
	  for(int i=0; i<nSysts; i++)
	    histos["mjjshapes_"+pf+"_exp"]->GetYaxis()->SetBinLabel(i+1,expSysts[i]);
	  
	  histos["mjjshapes_"+pf+"_gen"]=new TH2F("mjjshapes_"+pf+"_gen",";Mass(j,j') [GeV];Systematic uncertainty;Events",16,0,200,500,0,500);
	  for(int i=0; i<500;i++)
	    histos["mjjshapes_"+pf+"_gen"]->GetYaxis()->SetBinLabel(i+1,Form("genUnc%d",i));
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

  for(Int_t fileIter = 0; fileIter < nFiles; fileIter++){

    TString inF(inFileNames_p->at(fileIter));
    if(inF.Contains("/store") && !inF.Contains("root:")) inF="root://eoscms//eos/cms/"+inF;
    TFile* inFile_p = TFile::Open(inF, "READ");
    
    TTree* lepTree_p = (TTree*)inFile_p->Get("ggHiNtuplizer/EventTree");
    TTree* jetTree_p = (TTree*)inFile_p->Get("ak4PFJetAnalyzer/t");
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* hltTree_p = (TTree*)inFile_p->Get("hltanalysis/HltTree");
    TTree *pfCand_p  = (TTree *)inFile_p->Get("pfcandAnalyzer/pfTree");
    
    //PF candidates
    std::vector<int> *pfId_p=0;
    std::vector<float> *pfPt_p=0,*pfEta_p=0,*pfPhi_p=0,*pfEnergy_p=0;
    pfCand_p->SetBranchStatus("pfId",     1);
    pfCand_p->SetBranchStatus("pfPt",     1);
    pfCand_p->SetBranchStatus("pfEta",    1);
    pfCand_p->SetBranchStatus("pfPhi",    1);
    pfCand_p->SetBranchStatus("pfEnergy", 1);
    pfCand_p->SetBranchAddress("pfId",     &pfId_p);
    pfCand_p->SetBranchAddress("pfPt",     &pfPt_p);
    pfCand_p->SetBranchAddress("pfEta",    &pfEta_p);
    pfCand_p->SetBranchAddress("pfPhi",    &pfPhi_p);
    pfCand_p->SetBranchAddress("pfEnergy", &pfEnergy_p);

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
    Int_t refparton_flavor[maxJets];
    jetTree_p->SetBranchStatus("*", 0);
    jetTree_p->SetBranchStatus("nref", 1);
    jetTree_p->SetBranchStatus("jtpt", 1);
    jetTree_p->SetBranchStatus("jtphi", 1);
    jetTree_p->SetBranchStatus("jteta", 1);
    jetTree_p->SetBranchStatus("jtm", 1);
    jetTree_p->SetBranchStatus("discr_csvV2", 1);
    jetTree_p->SetBranchStatus("refpt", 1);
    jetTree_p->SetBranchStatus("refparton_flavorForB", 1);
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
    Float_t weight;
    std::vector<float> *ttbar_w_p=0;
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("weight", 1);
    hiTree_p->SetBranchStatus("ttbar_w",1);
    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("evt", &evt_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("weight", &weight);
    hiTree_p->SetBranchAddress("ttbar_w",&ttbar_w_p);
  
    //trigger
    int trig = 0;
    std::string triggerName;
    triggerName = isMC ? "HLT_HIL2Mu15ForPPRef_v1" : "HLT_HIL2Mu15_v1"; 
    hltTree_p->SetBranchStatus(triggerName.data(),1);
    hltTree_p->SetBranchAddress(triggerName.data(),&trig);
    
    Int_t nEntries = (Int_t)lepTree_p->GetEntries();
    
    std::cout << "Analysing " << nEntries << " events isMC=" << isMC
	      << " trigger=" << triggerName << endl;

    for(Int_t entry = 0; entry < nEntries; entry++)
      {
	if(entry%1000==0)
	  {
	    printf("\r [%d/%d] done",entry,nEntries);
	    cout << flush;
	  }

	//readout this event
	lepTree_p->GetEntry(entry);
	jetTree_p->GetEntry(entry);
	hiTree_p->GetEntry(entry);
	hltTree_p->GetEntry(entry);
	pfCand_p->GetEntry(entry);

	//assign an event weight
	float evWeight(1.0);
	if(isMC)
	  {
	    if(ttbar_w_p->size()) evWeight = ttbar_w_p->at(0);
	    evWeight *= totalEvtNorm;
	  }

	//select good muons
	//cf. details in https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
	std::vector<TLorentzVector> tightMuons,looseMuons,tightMuonsNonIso;
	std::vector<int> muonCharge;
	const Int_t nMu = (Int_t)muPt_p->size();
	for(Int_t muIter = 0; muIter < nMu; muIter++)
	  {
	    bool passLooseKin( muPt_p->at(muIter) > 15. && TMath::Abs(muEta_p->at(muIter))<2.4);
	    bool passTightKin( muPt_p->at(muIter) > 25. && TMath::Abs(muEta_p->at(muIter))<2.1);
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
	    bool passLooseIso( relIso<0.25);
	    
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
	    if(tightMuonsNonIso.size()==0) continue;
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
	    
	    if(TMath::Abs(elePt_p->at(eleIter)) < 20)  continue;
	    if(TMath::Abs(eleEta_p->at(eleIter)) > 2.5) continue;
	    if(eleIDVeto_p->at(eleIter)==0) continue;
	    nLooseEle++;
	  }
	if(nLooseEle>0) continue;

	//raw MET (noHF)
	TLorentzVector rawMET(0,0,0,0);	
	for(size_t ipf=0; ipf<pfId_p->size(); ipf++)
	  {
	    Float_t abseta=TMath::Abs(pfEta_p->at(ipf));
	    if(abseta>3.0) continue;
	    rawMET += TLorentzVector(-pfPt_p->at(ipf)*TMath::Cos(pfPhi_p->at(ipf)),
				     -pfPt_p->at(ipf)*TMath::Sin(pfPhi_p->at(ipf)),
				     0,
				     0);
	  }
	
	//transverse mass
	float mt(computeMT(tightMuons[0],rawMET));
	
	//jet counting
	typedef std::vector<TLorentzVector> JetColl_t;
	std::vector<JetColl_t> bJets(9),lightJets(9);
	for(Int_t jetIter = 0; jetIter < nref; jetIter++)
	  {
	    //cross clean with trigger muon
	    TLorentzVector jp4(0,0,0,0);
	    jp4.SetPtEtaPhiM(jtpt[jetIter],jteta[jetIter],jtphi[jetIter],jtm[jetIter]);
	    if(jp4.DeltaR(tightMuons[0])<0.4) continue;
	    
	    //in tracker region
	    if(TMath::Abs(jp4.Eta())>2.4) continue;

	    //systematic variations
	    Int_t jflav(abs(refparton_flavor[jetIter]));	    
	    bool passCSVM(discr_csvV2[jetIter]>0.8),passCSVMUp(passCSVM),passCSVMDn(passCSVM);	    
	    std::vector<float> jerSmear(3,1.0),jesScaleUnc(3,1.0);
	    if(isMC)
	      {
		//jet energy resolution smearing	
		if(refpt[jetIter]>0) jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),refpt[jetIter]);
		TLorentzVector rawjp4(jp4);
		jp4 *= jerSmear[0];

		jesScaleUnc[1]=1.028;
		jesScaleUnc[2]=0.972;

		//b-tagging
		float jptforBtag(jp4.Pt()>1000. ? 999. : jp4.Pt());
		if(jflav==5)
		  {
		    float expEff    = expBtagEff["b"]->Eval(jptforBtag); 
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMUp,1.1,expEff);	
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMDn,0.9,expEff);	
		  }
		else if(jflav==4)
		  {
		    float expEff    = expBtagEff["c"]->Eval(jptforBtag); 
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMUp,1.1,expEff);	
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMDn,0.9,expEff);	
		  }
		else
		  {
		    float expEff    = expBtagEff["udsg"]->Eval(jptforBtag); 
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMUp,1.3,expEff);	
		    myBTagSFUtil.modifyBTagsWithSF(passCSVMDn,0.7,expEff);	
		  }		
	      }

	    if(jp4.Pt()>30)
	      {
		//nominal selection
		if(passCSVM) bJets[0].push_back(jp4);
		else         lightJets[0].push_back(jp4);

		//tag variations affect differently depending on the flavour
		if(jflav==5 || jflav==4)
		  {
		    if(passCSVMUp)
		      {
			bJets[1].push_back(jp4);
		      }
		    else
		      {
			lightJets[1].push_back(jp4);
		      }
		    if(passCSVMDn) 
		      {
			bJets[2].push_back(jp4);
		      }
		    else
		      {
			lightJets[2].push_back(jp4);
		      }
		    if(passCSVM)   
		      {
			bJets[3].push_back(jp4);
			bJets[4].push_back(jp4);
		      }
		    else      
		      {
			lightJets[3].push_back(jp4);
			lightJets[4].push_back(jp4);
		      }
		  }
		else
		  {
		    if(passCSVM)   
		      {
			bJets[1].push_back(jp4);
			bJets[2].push_back(jp4);
		      }
		    else      
		      {
			lightJets[1].push_back(jp4);
			lightJets[2].push_back(jp4);
		      }
		    if(passCSVMUp) 
		      {
			bJets[3].push_back(jp4);
		      }
		    else
		      {
			lightJets[3].push_back(jp4);
		      }
		    if(passCSVMDn) 
		      {
			bJets[4].push_back(jp4);
		      }
		    else
		      {
			lightJets[4].push_back(jp4);
		      }
		  }
	      }
	    
	    for(size_t ivar=0; ivar<2; ivar++)
	      {
		//JES varied selections
		TLorentzVector jesVarP4(jp4); jesVarP4*=jesScaleUnc[ivar+1];
		if(jesVarP4.Pt()>30)
		  {
		    if(passCSVM) bJets[5+ivar].push_back(jesVarP4);
		    else         lightJets[5+ivar].push_back(jesVarP4);
		  }

		//JER varied selections
		TLorentzVector jerVarP4(jp4); jerVarP4*=jerSmear[ivar+1]/jerSmear[0];     
		if(jerVarP4.Pt()>30)
		  {
		    if(passCSVM) bJets[7+ivar].push_back(jerVarP4);
		    else         lightJets[7+ivar].push_back(jerVarP4);
		  }
	      }
	  }
	
	histos["lpt"]->Fill(tightMuons[0].Pt(),evWeight);
	histos["leta"]->Fill(fabs(tightMuons[0].Eta()),evWeight);
	histos["mt"]->Fill(mt,evWeight);
	histos["metpt"]->Fill(rawMET.Pt(),evWeight);

	//
	for(Int_t ivar=0; ivar<=nSysts; ivar++)
	  {
	    Int_t jetIdx(0);
	    if(ivar>=1 && ivar<=8) jetIdx=ivar;

	    //require at least two light jet acompanying the lepton
	    Int_t nljets(lightJets[jetIdx].size());
	    Int_t nbtags(bJets[jetIdx].size());
	    Int_t njets(nljets+nbtags);
	    
	    if(nljets<2) continue;
	    TString pf(Form("%db",TMath::Min(nbtags,2)));
	    
	    //jet-related quantities
	    Float_t mjj( (lightJets[jetIdx][0]+lightJets[jetIdx][1]).M() );
	    Float_t htsum(0);
	    for(Int_t ij=0; ij<nbtags; ij++) htsum += bJets[jetIdx][ij].Pt();
	    for(Int_t ij=0; ij<nljets; ij++) htsum += lightJets[jetIdx][ij].Pt();

	    Float_t mlb( (tightMuons[0]+lightJets[jetIdx][0]).M() );
	    if(nbtags>0)
	      {
		mlb=(tightMuons[0]+bJets[jetIdx][0]).M();
		if(nbtags>1)
		  {
		    mlb=TMath::Min( mlb, Float_t((tightMuons[0]+bJets[jetIdx][1]).M()) );
		  }
	      }

	    //update event weight if needed
	    Float_t iweight(evWeight);
	    if(ivar==9)  iweight*=1.03;
	    if(ivar==10) iweight*=1.03;
	  
	    //fill histos
	    if(ivar==0)
	      {
		histos["lpt_"+pf]->Fill(tightMuons[0].Pt(),iweight);
		histos["leta_"+pf]->Fill(fabs(tightMuons[0].Eta()),iweight);
		histos["ht_"+pf]->Fill(htsum,iweight);
		histos["mjj_"+pf]->Fill(mjj,iweight);
		histos["mlb_"+pf]->Fill(mlb,iweight);
	       
		if(runSysts)
		  {
		    //theory uncertainties (by matrix-element weighting)
		    for(size_t igs=0; igs<ttbar_w_p->size(); igs++)
		      {
			float newWeight( iweight );
			if(isTTJets && normH && normH->GetBinContent(igs+1))
			  {
			    newWeight *= (ttbar_w_p->at(igs)/ttbar_w_p->at(0)) * ( normH->GetBinContent(1)/normH->GetBinContent(igs+1));
			  }
			else
			  {
			    newWeight *= (ttbar_w_p->at(igs)/ttbar_w_p->at(0));
			  }
			((TH2 *)histos["mjjshapes_"+pf+"_gen"])->Fill(mjj,igs,newWeight);
		      }
		  }
		histos["metpt_"+pf]->Fill(rawMET.Pt(),iweight);
		histos["metphi_"+pf]->Fill(rawMET.Phi(),iweight);
		histos["mt_"+pf]->Fill(mt,iweight);    
		if(nbtags)
		  {
		    histos["jpt_"+pf]->Fill(bJets[jetIdx][0].Pt(),iweight);
		    histos["jeta_"+pf]->Fill(fabs(bJets[jetIdx][0].Eta()),iweight);
		  }
		else
		  {
		    histos["jpt_"+pf]->Fill(lightJets[jetIdx][0].Pt(),iweight);
		    histos["jeta_"+pf]->Fill(fabs(lightJets[jetIdx][0].Eta()),iweight);
		  }
		histos["njets_"+pf]->Fill(njets,iweight);
	      }
	    else if (runSysts)
	      {
		((TH2 *)histos["mjjshapes_"+pf+"_exp"])->Fill(mjj,ivar-1,iweight);
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
      it->second->SetDirectory(outFile_p);
      it->second->Write();
    }
  outFile_p->Close();
  delete outFile_p;
  
  return;
}
