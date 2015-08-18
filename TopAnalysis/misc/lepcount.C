void lepcount(){


TString hnames[]={"lepCh_2j_leading","lepCh_3j_leading","lepCh_4j_leading"};
for(size_t i=0; i<sizeof(hnames)/sizeof(TString); i++)
    runComparison(hnames[i]);
}

void runComparison(TString hname)
{
	Float_t wXsec(50100.0);

	TString baseDir = "/afs/cern.ch/user/q/qhassan/public/Analysis/CMSSW_7_0_6_patch1/src/UserCode/TopAnalysis/test/analysis/data/";

	TFile* f0  = TFile::Open(baseDir+"chplots/miniAOD_wjets.root");
	TH1* h1 = (TH1*)f0->Get(hname);
	TH1 *hcutflowf0=(TH1 *) f0->Get("cutflow");
	Float_t norigEventsf0(hcutflowf0->GetBinContent(1));
	h1->Scale(wXsec/norigEventsf0);

  
cout<<"Number of Negative Leptons in W: "<< h1->GetBinContent(2)*100 <<endl <<"Number of Positive Leptons in W: " << h1->GetBinContent(5)*100 << endl;

	Float_t DYXsec(4746.0);

	TFile* f1  = TFile::Open(baseDir+"chplots/miniAOD_DY.root");
	TH1* h2 = (TH1*)f1->Get(hname);
	TH1 *hcutflowf1=(TH1 *) f1->Get("cutflow");
        Float_t norigEventsf1(hcutflowf1->GetBinContent(1));
        h2->Scale(DYXsec/norigEventsf1);

cout<<"Number of Negative Leptons in DY: "<< h2->GetBinContent(2) <<endl <<"Number of Positive Leptons in DY: " << h2->GetBinContent(5) << endl;

	Float_t ttbarXsec(827.05);

	TFile* f2  = TFile::Open(baseDir+"chplots/miniAOD_MC_v2.root");
	TH1* h3 = (TH1*)f2->Get(hname);
	TH1 *hcutflowf2=(TH1 *) f2->Get("cutflow");
        Float_t norigEventsf2(hcutflowf2->GetBinContent(1));
        h3->Scale(ttbarXsec/norigEventsf2);

cout<<"Number of Negative Leptons in MC: "<< h3->GetBinContent(2) <<endl <<"Number of Positive Leptons in MC: " << h3->GetBinContent(5) << endl;

	// QCD	
	 
	Float_t qcd1Xsec(3000114.3);

	TFile* f3  = TFile::Open(baseDir+"chplots/miniAOD_QCD_80_120.root");
	TH1* h4 = (TH1*)f3->Get(hname);
        TH1 *hcutflowf3=(TH1 *) f3->Get("cutflow");
        Float_t norigEventsf3(hcutflowf3->GetBinContent(1));
        h4->Scale(qcd1Xsec/norigEventsf3);

//	h4->Scale(qcd1Xsec/h4->Integral());

//cout<<"Number of Negative Leptons in QCD1: "<< h4->GetBinContent(2) <<endl <<"Number of Positive Leptons in QCD1: " << h4->GetBinContent(5) << endl;

//	cout<<"Number of Negative Leptons in QCD: "<< h4->GetBinContent(2) <<endl <<"Number of Positive Leptons in QCD: " << h4->GetBinContent(5) << endl;
	
	Float_t qcd2Xsec(493200);

	TFile* f4  = TFile::Open(baseDir+"chplots/miniAOD_QCD_120_170.root");
	TH1* h5 = (TH1*)f4->Get(hname);
	TH1 *hcutflowf4=(TH1 *) f4->Get("cutflow");
        Float_t norigEventsf4(hcutflowf4->GetBinContent(1));
        h5->Scale(qcd2Xsec/norigEventsf4);

//	h5->Scale(qcd2Xsec/h5->Integral());	  
 
//cout<<"Number of Negative Leptons in QCD2: "<< h5->GetBinContent(2) <<endl <<"Number of Positive Leptons in QCD2: " << h5->GetBinContent(5) << endl;
/*
	TH1F* h = (TH1*)h4->Clone();
	h->SetName("h");
	h->Reset("ICE");
	h->Add(h5);
*/
	h4->Add(h5);
//cout<<"Number of Negative Leptons in QCDtotal: "<< h4->GetBinContent(2) <<endl <<"Number of Positive Leptons in QCDtotal: " << h4->GetBinContent(5) << endl;

	Float_t qcd3Xsec(120300);

	TFile* f5  = TFile::Open(baseDir+"chplots/miniAOD_QCD_170_300.root");
	TH1* h6 = (TH1*)f5->Get(hname);
	TH1 *hcutflowf5=(TH1 *) f5->Get("cutflow");
        Float_t norigEventsf5(hcutflowf5->GetBinContent(1));
        h6->Scale(qcd3Xsec/norigEventsf5);

//	h6->Scale(qcd3Xsec/h6->Integral());
	  
	h4->Add(h6);

//cout<<"Number of Negative Leptons in QCD3: "<< h6->GetBinContent(2) <<endl <<"Number of Positive Leptons in QCD3: " << h6->GetBinContent(5) << endl;
//cout<<"Number of Negative Leptons in QCDtotal grand: "<< h4->GetBinContent(2) <<endl <<"Number of Positive Leptons in QCDtotal grand: " << h4->GetBinContent(5) << endl;


	Float_t qcd4Xsec(7475);

	TFile* f6  = TFile::Open(baseDir+"chplots/miniAOD_QCD_300_470.root");
	TH1* h7 = (TH1*)f6->Get(hname);
        TH1 *hcutflowf6=(TH1 *) f6->Get("cutflow");
        Float_t norigEventsf6(hcutflowf6->GetBinContent(1));
        h7->Scale(qcd4Xsec/norigEventsf6);

//	h7->Scale(qcd4Xsec/h7->Integral());

	h4->Add(h7);

	Float_t qcd5Xsec(587.1);

	TFile* f7  = TFile::Open(baseDir+"chplots/miniAOD_QCD_470_600.root");
	TH1* h8 = (TH1*)f7->Get(hname);
        TH1 *hcutflowf7=(TH1 *) f7->Get("cutflow");
        Float_t norigEventsf7(hcutflowf7->GetBinContent(1));
        h8->Scale(qcd5Xsec/norigEventsf7);

//	h8->Scale(qcd5Xsec/h8->Integral());

	h4->Add(h8);

	Float_t qcd6Xsec(167);

	TFile* f8  = TFile::Open(baseDir+"chplots/miniAOD_QCD_600_800.root");
	TH1* h9 = (TH1*)f8->Get(hname);
        TH1 *hcutflowf8=(TH1 *) f8->Get("cutflow");
        Float_t norigEventsf8(hcutflowf8->GetBinContent(1));
        h9->Scale(qcd6Xsec/norigEventsf8);
	
//h9->Scale(qcd6Xsec/h9->Integral());

	h4->Add(h9);

	Float_t qcd7Xsec(28.25);

	TFile* f9  = TFile::Open(baseDir+"chplots/miniAOD_QCD_800_1000.root");
	TH1* h10 = (TH1*)f9->Get(hname);
        TH1 *hcutflowf9=(TH1 *) f9->Get("cutflow");
        Float_t norigEventsf9(hcutflowf9->GetBinContent(1));
        h10->Scale(qcd7Xsec/norigEventsf9);

//	h10->Scale(qcd7Xsec/h10->Integral());

	h4->Add(h10);

	Float_t qcd8Xsec(8.195);

	TFile* f10  = TFile::Open(baseDir+"chplots/miniAOD_QCD_1000_1400.root");
	TH1* h11 = (TH1*)f10->Get(hname);
        TH1 *hcutflowf10=(TH1 *) f10->Get("cutflow");
        Float_t norigEventsf10(hcutflowf10->GetBinContent(1));
        h11->Scale(qcd8Xsec/norigEventsf10);
	
//h11->Scale(qcd8Xsec/h11->Integral());

	h4->Add(h11);

	Float_t qcd9Xsec(0.7346);

	TFile* f11  = TFile::Open(baseDir+"chplots/miniAOD_QCD_1400_1800.root");
	TH1* h12 = (TH1*)f11->Get(hname);
        TH1 *hcutflowf11=(TH1 *) f11->Get("cutflow");
        Float_t norigEventsf11(hcutflowf11->GetBinContent(1));
        h12->Scale(qcd9Xsec/norigEventsf11);

//	h12->Scale(qcd9Xsec/h12->Integral());

	h4->Add(h12);

	Float_t qcd10Xsec(0.102);

	TFile* f12  = TFile::Open(baseDir+"chplots/miniAOD_QCD_1800_2400.root");
	TH1* h13 = (TH1*)f12->Get(hname);
        TH1 *hcutflowf12=(TH1 *) f12->Get("cutflow");
        Float_t norigEventsf12(hcutflowf12->GetBinContent(1));
        h13->Scale(qcd10Xsec/norigEventsf12);

//	h13->Scale(qcd10Xsec/h13->Integral());

	h4->Add(h13);

	Float_t qcd11Xsec(0.00644);

	TFile* f13  = TFile::Open(baseDir+"chplots/miniAOD_QCD_2400_3200.root");
	TH1* h14 = (TH1*)f13->Get(hname);
        TH1 *hcutflowf13=(TH1 *) f13->Get("cutflow");
        Float_t norigEventsf13(hcutflowf13->GetBinContent(1));
        h14->Scale(qcd11Xsec/norigEventsf13);

//	h14->Scale(qcd11Xsec/h14->Integral());

	h4->Add(h14);

	Float_t qcd12Xsec(0.000163);

	TFile* f14  = TFile::Open(baseDir+"chplots/miniAOD_QCD_3200.root");
	TH1* h15 = (TH1*)f14->Get(hname);
        TH1 *hcutflowf14=(TH1 *) f14->Get("cutflow");
        Float_t norigEventsf14(hcutflowf14->GetBinContent(1));
        h15->Scale(qcd12Xsec/norigEventsf14);

//	h15->Scale(qcd12Xsec/h15->Integral());

	h4->Add(h15);

cout<<"Number of Negative Leptons in QCD: "<< h4->GetBinContent(2) <<endl <<"Number of Positive Leptons in QCD: " << h4->GetBinContent(5) << endl;

	// Single Top
	 
	Float_t st1Xsec(6.35);

	TFile* f15  = TFile::Open(baseDir+"chplots/miniAOD_TToLeptons_s_channel.root");
	TH1* h16 = (TH1*)f15->Get(hname);
        TH1 *hcutflowf15=(TH1 *) f15->Get("cutflow");
        Float_t norigEventsf15(hcutflowf15->GetBinContent(1));
        h16->Scale(st1Xsec/norigEventsf15);

//	h16->Scale(st1Xsec/h16->Integral());

//cout<<"Number of Negative Leptons in stop1: "<< h16->GetBinContent(2) <<endl <<"Number of Positive Leptons in stop1: " << h16->GetBinContent(5) << endl;

	Float_t st2Xsec(3.97);

	TFile* f16  = TFile::Open(baseDir+"chplots/miniAOD_TBarToLeptons_s_channel.root");
	TH1* h17 = (TH1*)f16->Get(hname);
        TH1 *hcutflowf16=(TH1 *) f16->Get("cutflow");
        Float_t norigEventsf16(hcutflowf16->GetBinContent(1));
        h17->Scale(st2Xsec/norigEventsf16);

//	h17->Scale(st2Xsec/h17->Integral());

	//TH1F* h_stop = (TH1*)h16->Clone();
//	h_stop->SetName("h_stop");
//	h_stop->Reset("ICE");
	h16->Add(h17);

//cout<<"Number of Negative Leptons in stop2: "<< h17->GetBinContent(2) <<endl <<"Number of Positive Leptons in stop2: " << h17->GetBinContent(5) << endl;
//cout<<"Number of Negative Leptons in stop: "<< h16->GetBinContent(2) <<endl <<"Number of Positive Leptons in stop: " << h16->GetBinContent(5) << endl;

	Float_t st3Xsec(103.4);

	TFile* f17  = TFile::Open(baseDir+"chplots/miniAOD_TToLeptons_t_channel.root");
	TH1* h18 = (TH1*)f17->Get(hname);
        TH1 *hcutflowf17=(TH1 *) f17->Get("cutflow");
        Float_t norigEventsf17(hcutflowf17->GetBinContent(1));
        h18->Scale(st3Xsec/norigEventsf17);

//	h18->Scale(st3Xsec/h18->Integral());

	h16->Add(h18);

	Float_t st4Xsec(61.56);

	TFile* f18  = TFile::Open(baseDir+"chplots/miniAOD_TBarToLeptons_t_channel.root");
	TH1* h19 = (TH1*)f18->Get(hname);
        TH1 *hcutflowf18=(TH1 *) f18->Get("cutflow");
        Float_t norigEventsf18(hcutflowf18->GetBinContent(1));
        h19->Scale(st4Xsec/norigEventsf18);

//	h19->Scale(st4Xsec/h19->Integral());

	h16->Add(h19);

	Float_t st5Xsec(35.85);

	TFile* f19  = TFile::Open(baseDir+"chplots/miniAOD_T_tW_channel.root");
	TH1* h20 = (TH1*)f19->Get(hname);
        TH1 *hcutflowf19=(TH1 *) f19->Get("cutflow");
        Float_t norigEventsf19(hcutflowf19->GetBinContent(1));
        h20->Scale(st5Xsec/norigEventsf19);

//	h20->Scale(st5Xsec/h20->Integral());

	h16->Add(h20);

	Float_t st6Xsec(35.85);

	TFile* f20  = TFile::Open(baseDir+"chplots/miniAOD_Tbar_tW_channel.root");
	TH1* h21 = (TH1*)f20->Get(hname);
        TH1 *hcutflowf20=(TH1 *) f20->Get("cutflow");
        Float_t norigEventsf20(hcutflowf20->GetBinContent(1));
        h21->Scale(st6Xsec/norigEventsf20);
	
//h21->Scale(st6Xsec/h21->Integral());

	h16->Add(h21);
	
cout<<"Number of Negative Leptons in S_TOP: "<< h16->GetBinContent(2) <<endl <<"Number of Positive Leptons in S_TOP: " << h16->GetBinContent(5) << endl;

	
	}
