#include<iostream>
#include<iomanip>
#include<sstream>

double pos_ev(0), neg_ev(0), pos_wev(0), neg_wev(0);
float w_ev(0);
int lumi = 100;

    std::string prd(const double x, const int decDigits, const int width) {
    stringstream ss;
//    ss << fixed << right;
    ss << fixed;
    ss.fill(' ');        // fill space around displayed #
    ss.width(width);     // set  width around displayed #
    ss.precision(decDigits); // set # places after decimal
    ss << x;
    return ss.str();
}

    std::string center(const string s, const int w){
    stringstream ss, spaces;
    int padding = w - s.size();                 // count excess room to pad
    for(int i=0; i<padding/2; ++i)
        spaces << " ";
    ss << spaces.str() << s << spaces.str();    // format with padding
    if(padding>0 && padding%2!=0)               // if odd #, add 1 space
        ss << " ";
    return ss.str();
}

std::string align(const string& s, int w, bool left=true){
    int diff=w-s.size();
    if (diff<0){diff=0;}
    std::string padding(diff, ' ');
    stringstream ss;
    if (left){ss << s << padding;}
    else{ss <<padding<<s;}
    return ss.str();
}

//TH2F *error_hist     = new TH2F("error_hist",";error_hist;Ratio" ,10,0.,500.,10,0.,2);
TH1F *error_hist     = new TH1F("error_hist","Ratio" ,9,0,9);

void ele_lepcount(){

cout<<"                   CLOSURE TESTS         "<<endl;								//

TString hnames[]={"lepCh_2j_leading","lepCh_3j_leading","lepCh_4j_leading"};
std::string titles[]={"2j leading","3j leading","4j leading"};

error_hist->GetXaxis()->SetBinLabel(1,"2j (W+DY+t)");
error_hist->GetXaxis()->SetBinLabel(2,"2j (W+DY+t+t#bar{t})");
error_hist->GetXaxis()->SetBinLabel(3,"2j (W+DY+t+t#bar{t}+QCD)");
error_hist->GetXaxis()->SetBinLabel(4,"3j (W+DY+t)");
error_hist->GetXaxis()->SetBinLabel(5,"3j (W+DY+t+t#bar{t})");
error_hist->GetXaxis()->SetBinLabel(6,"3j (W+DY+t+t#bar{t}+QCD)");
error_hist->GetXaxis()->SetBinLabel(7,"4j (W+DY+t)");
error_hist->GetXaxis()->SetBinLabel(8,"4j (W+DY+t+t#bar{t})");
error_hist->GetXaxis()->SetBinLabel(9,"4j (W+DY+t+t#bar{t}+QCD)");

for(size_t i=0; i<sizeof(hnames)/sizeof(TString); i++)
    runComparison(hnames[i], titles[i].c_str(),i*3+1);


error_hist->Draw("E1");

}

void runComparison(TString hname,const char* title, size_t bin_num)
{

	// W+Jets
	Float_t wXsec(50100.0);

//	TString baseDir = "/afs/cern.ch/user/q/qhassan/public/Analysis/CMSSW_7_0_6_patch1/src/UserCode/TopAnalysis/test/analysis/data/";
	TString baseDir = "/afs/cern.ch/user/q/qhassan/public/Analysis/CMSSW_7_0_6_patch1/src/UserCode/TopAnalysis/test/analysis/data/";

	TFile* f0  = TFile::Open(baseDir+"ele_plots/miniAOD_wjets_PU20bx25_new.root");
	TH1* h1 = (TH1*)f0->Get(hname);
float no_scale_negw = h1->GetBinContent(2);
float no_scale_posw = h1->GetBinContent(5);
	TH1 *hcutflowf0=(TH1 *) f0->Get("cutflow");
	Float_t norigEventsf0(hcutflowf0->GetBinContent(1));
	h1->Scale(wXsec/norigEventsf0);

	// DY
	Float_t DYXsec(4746.0);

	TFile* f1  = TFile::Open(baseDir+"ele_plots/miniAOD_DY_PU20bx25_new.root");
	TH1* h2 = (TH1*)f1->Get(hname);
float no_scale_negdy = h2->GetBinContent(2);
float no_scale_posdy = h2->GetBinContent(5);
	TH1 *hcutflowf1=(TH1 *) f1->Get("cutflow");
        Float_t norigEventsf1(hcutflowf1->GetBinContent(1));
        h2->Scale(DYXsec/norigEventsf1);
	
	//tt-bar
	Float_t ttbarXsec(827.05);

	TFile* f2  = TFile::Open(baseDir+"ele_plots/miniAOD_MC_v2_new.root");
	TH1* h3 = (TH1*)f2->Get(hname);
float no_scale_negttbar = h3->GetBinContent(2);
float no_scale_posttbar = h3->GetBinContent(5);
	TH1 *hcutflowf2=(TH1 *) f2->Get("cutflow");
        Float_t norigEventsf2(hcutflowf2->GetBinContent(1));
        h3->Scale(ttbarXsec/norigEventsf2);

	// QCD	
	 
	Float_t qcd1Xsec(3000114.3);

	TFile* f3  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_80_120_Tune4C.root");
	TH1* h4 = (TH1*)f3->Get(hname);
float no_scale_negqcd1 = h4->GetBinContent(2);
float no_scale_posqcd1 = h4->GetBinContent(5);
        TH1 *hcutflowf3=(TH1 *) f3->Get("cutflow");
        Float_t norigEventsf3(hcutflowf3->GetBinContent(1));
        h4->Scale(qcd1Xsec/norigEventsf3);
	
	Float_t qcd2Xsec(493200);

	TFile* f4  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_120_170_Tune4C.root");
	TH1* h5 = (TH1*)f4->Get(hname);
float no_scale_negqcd2 = h5->GetBinContent(2);
float no_scale_posqcd2 = h5->GetBinContent(5);

	TH1 *hcutflowf4=(TH1 *) f4->Get("cutflow");
        Float_t norigEventsf4(hcutflowf4->GetBinContent(1));
        h5->Scale(qcd2Xsec/norigEventsf4);
	h4->Add(h5);

	Float_t qcd3Xsec(120300);

	TFile* f5  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_170_300_Tune4C.root");
	TH1* h6 = (TH1*)f5->Get(hname);
float no_scale_negqcd3 = h6->GetBinContent(2);
float no_scale_posqcd3 = h6->GetBinContent(5);
	TH1 *hcutflowf5=(TH1 *) f5->Get("cutflow");
        Float_t norigEventsf5(hcutflowf5->GetBinContent(1));
        h6->Scale(qcd3Xsec/norigEventsf5);

	h4->Add(h6);

	Float_t qcd4Xsec(7475);

	TFile* f6  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_300_470_Tune4C.root");
	TH1* h7 = (TH1*)f6->Get(hname);
float no_scale_negqcd4 = h7->GetBinContent(2);
float no_scale_posqcd4 = h7->GetBinContent(5);
        TH1 *hcutflowf6=(TH1 *) f6->Get("cutflow");
        Float_t norigEventsf6(hcutflowf6->GetBinContent(1));
        h7->Scale(qcd4Xsec/norigEventsf6);

	h4->Add(h7);

	Float_t qcd5Xsec(587.1);

	TFile* f7  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_470_600_Tune4C.root");
	TH1* h8 = (TH1*)f7->Get(hname);
float no_scale_negqcd5 = h8->GetBinContent(2);
float no_scale_posqcd5 = h8->GetBinContent(5);
        TH1 *hcutflowf7=(TH1 *) f7->Get("cutflow");
        Float_t norigEventsf7(hcutflowf7->GetBinContent(1));
        h8->Scale(qcd5Xsec/norigEventsf7);

	h4->Add(h8);

	Float_t qcd6Xsec(167);

	TFile* f8  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_600_800_Tune4C.root");
	TH1* h9 = (TH1*)f8->Get(hname);
float no_scale_negqcd6 = h9->GetBinContent(2);
float no_scale_posqcd6 = h9->GetBinContent(5);
        TH1 *hcutflowf8=(TH1 *) f8->Get("cutflow");
        Float_t norigEventsf8(hcutflowf8->GetBinContent(1));
        h9->Scale(qcd6Xsec/norigEventsf8);

	h4->Add(h9);

	Float_t qcd7Xsec(28.25);

	TFile* f9  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_800_1000_Tune4C.root");
	TH1* h10 = (TH1*)f9->Get(hname);
float no_scale_negqcd7 = h10->GetBinContent(2);
float no_scale_posqcd7 = h10->GetBinContent(5);
        TH1 *hcutflowf9=(TH1 *) f9->Get("cutflow");
        Float_t norigEventsf9(hcutflowf9->GetBinContent(1));
        h10->Scale(qcd7Xsec/norigEventsf9);

	h4->Add(h10);

	Float_t qcd8Xsec(8.195);

	TFile* f10  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_1000_1400_Tune4C.root");
	TH1* h11 = (TH1*)f10->Get(hname);
float no_scale_negqcd8 = h11->GetBinContent(2);
float no_scale_posqcd8 = h11->GetBinContent(5);
        TH1 *hcutflowf10=(TH1 *) f10->Get("cutflow");
        Float_t norigEventsf10(hcutflowf10->GetBinContent(1));
        h11->Scale(qcd8Xsec/norigEventsf10);
	
	h4->Add(h11);

	Float_t qcd9Xsec(0.7346);

	TFile* f11  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_1400_1800_Tune4C.root");
	TH1* h12 = (TH1*)f11->Get(hname);
float no_scale_negqcd9 = h12->GetBinContent(2);
float no_scale_posqcd9 = h12->GetBinContent(5);
        TH1 *hcutflowf11=(TH1 *) f11->Get("cutflow");
        Float_t norigEventsf11(hcutflowf11->GetBinContent(1));
        h12->Scale(qcd9Xsec/norigEventsf11);


	h4->Add(h12);

	Float_t qcd10Xsec(0.102);

	TFile* f12  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_1800_2400_Tune4C.root");
	TH1* h13 = (TH1*)f12->Get(hname);
float no_scale_negqcd10 = h13->GetBinContent(2);
float no_scale_posqcd10 = h13->GetBinContent(5);
        TH1 *hcutflowf12=(TH1 *) f12->Get("cutflow");
        Float_t norigEventsf12(hcutflowf12->GetBinContent(1));
        h13->Scale(qcd10Xsec/norigEventsf12);


	h4->Add(h13);

	Float_t qcd11Xsec(0.00644);

	TFile* f13  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_2400_3200_Tune4C.root");
	TH1* h14 = (TH1*)f13->Get(hname);
float no_scale_negqcd11 = h14->GetBinContent(2);
float no_scale_posqcd11 = h14->GetBinContent(5);
        TH1 *hcutflowf13=(TH1 *) f13->Get("cutflow");
        Float_t norigEventsf13(hcutflowf13->GetBinContent(1));
        h14->Scale(qcd11Xsec/norigEventsf13);


	h4->Add(h14);

	Float_t qcd12Xsec(0.000163);

	TFile* f14  = TFile::Open(baseDir+"ele_plots/miniAOD_QCD_3200_Tune4C.root");
	TH1* h15 = (TH1*)f14->Get(hname);
float no_scale_negqcd12 = h15->GetBinContent(2);
float no_scale_posqcd12 = h15->GetBinContent(5);
        TH1 *hcutflowf14=(TH1 *) f14->Get("cutflow");
        Float_t norigEventsf14(hcutflowf14->GetBinContent(1));
        h15->Scale(qcd12Xsec/norigEventsf14);

float total_qcd_neg = no_scale_negqcd1 + no_scale_negqcd2 + no_scale_negqcd3 + no_scale_negqcd4 + no_scale_negqcd5 + no_scale_negqcd6 + no_scale_negqcd7 + no_scale_negqcd8 + no_scale_negqcd9 + no_scale_negqcd10 + no_scale_negqcd11 + no_scale_negqcd12;
float total_qcd_pos = no_scale_posqcd1 + no_scale_posqcd2 + no_scale_posqcd3 + no_scale_posqcd4 + no_scale_posqcd5 + no_scale_posqcd6 + no_scale_posqcd7 + no_scale_posqcd8 + no_scale_posqcd9 + no_scale_posqcd10 + no_scale_posqcd11 + no_scale_posqcd12;

	h4->Add(h15);


	// Single Top
	 
	Float_t st1Xsec(6.35);

	TFile* f15  = TFile::Open(baseDir+"ele_plots/miniAOD_TToLeptons_s_channel_new.root");
	TH1* h16 = (TH1*)f15->Get(hname);
float no_scale_negst1 = h16->GetBinContent(2);
float no_scale_posst1 = h16->GetBinContent(5);
        TH1 *hcutflowf15=(TH1 *) f15->Get("cutflow");
        Float_t norigEventsf15(hcutflowf15->GetBinContent(1));
        h16->Scale(st1Xsec/norigEventsf15);
	Float_t st2Xsec(3.97);

	TFile* f16  = TFile::Open(baseDir+"ele_plots/miniAOD_TBarToLeptons_s_channel_new.root");
	TH1* h17 = (TH1*)f16->Get(hname);
float no_scale_negst2 = h17->GetBinContent(2);
float no_scale_posst2 = h17->GetBinContent(5);
        TH1 *hcutflowf16=(TH1 *) f16->Get("cutflow");
        Float_t norigEventsf16(hcutflowf16->GetBinContent(1));
        h17->Scale(st2Xsec/norigEventsf16);
	h16->Add(h17);

	Float_t st3Xsec(103.4);

	TFile* f17  = TFile::Open(baseDir+"ele_plots/miniAOD_TToLeptons_t_channel_new.root");
	TH1* h18 = (TH1*)f17->Get(hname);
float no_scale_negst3 = h18->GetBinContent(2);
float no_scale_posst3 = h18->GetBinContent(5);
        TH1 *hcutflowf17=(TH1 *) f17->Get("cutflow");
        Float_t norigEventsf17(hcutflowf17->GetBinContent(1));
        h18->Scale(st3Xsec/norigEventsf17);
	h16->Add(h18);

	Float_t st4Xsec(61.56);

	TFile* f18  = TFile::Open(baseDir+"ele_plots/miniAOD_TBarToLeptons_t_channel_new.root");
	TH1* h19 = (TH1*)f18->Get(hname);
float no_scale_negst4 = h19->GetBinContent(2);
float no_scale_posst4 = h19->GetBinContent(5);
        TH1 *hcutflowf18=(TH1 *) f18->Get("cutflow");
        Float_t norigEventsf18(hcutflowf18->GetBinContent(1));
        h19->Scale(st4Xsec/norigEventsf18);
	h16->Add(h19);

	Float_t st5Xsec(35.85);

	TFile* f19  = TFile::Open(baseDir+"ele_plots/miniAOD_T_tW_channel_new.root");
	TH1* h20 = (TH1*)f19->Get(hname);
float no_scale_negst5 = h20->GetBinContent(2);
float no_scale_posst5 = h20->GetBinContent(5);
        TH1 *hcutflowf19=(TH1 *) f19->Get("cutflow");
        Float_t norigEventsf19(hcutflowf19->GetBinContent(1));
        h20->Scale(st5Xsec/norigEventsf19);
	h16->Add(h20);

	Float_t st6Xsec(35.85);

	TFile* f20  = TFile::Open(baseDir+"ele_plots/miniAOD_Tbar_tW_channel_new.root");
	TH1* h21 = (TH1*)f20->Get(hname);
float no_scale_negst6 = h21->GetBinContent(2);
float no_scale_posst6 = h21->GetBinContent(5);
        TH1 *hcutflowf20=(TH1 *) f20->Get("cutflow");
        Float_t norigEventsf20(hcutflowf20->GetBinContent(1));
        h21->Scale(st6Xsec/norigEventsf20);
float total_stneg = no_scale_negst1 + no_scale_negst2 + no_scale_negst3 + no_scale_negst4 + no_scale_negst5 + no_scale_negst6;
float total_stpos = no_scale_posst1 + no_scale_posst2 + no_scale_posst3 + no_scale_posst4 + no_scale_posst5 + no_scale_posst6;
	h16->Add(h21);
	
pos_ev = h1->GetBinContent(5)*lumi + h2->GetBinContent(5)*lumi + h3->GetBinContent(5)*lumi + h4->GetBinContent(5)*lumi + h16->GetBinContent(5)*lumi;
neg_ev = h1->GetBinContent(2)*lumi + h2->GetBinContent(2)*lumi + h3->GetBinContent(2)*lumi + h4->GetBinContent(2)*lumi + h16->GetBinContent(2)*lumi;
pos_wev = h1->GetBinContent(5)*lumi;
neg_wev = h1->GetBinContent(2)*lumi;
float pos_mc = h3->GetBinContent(5)*lumi, neg_mc = h3->GetBinContent(2)*lumi, pos_dy = h2->GetBinContent(5)*lumi, neg_dy = h2->GetBinContent(2)*lumi, pos_sintop = h16->GetBinContent(5)*lumi, neg_sintop = h16->GetBinContent(2)*lumi, pos_qcd = h4->GetBinContent(5)*lumi, neg_qcd = h4->GetBinContent(2)*lumi, diff = pos_wev - neg_wev, sum = pos_wev + neg_wev, Rw = diff / sum;

// QCD ERRORS
//
float error_qcd1_pos = ((sqrt(no_scale_posqcd1)*lumi*qcd1Xsec)/norigEventsf3)/Rw, 
error_qcd2_pos = ((sqrt(no_scale_posqcd2)*lumi*qcd2Xsec)/norigEventsf4)/Rw, 
error_qcd3_pos = ((sqrt(no_scale_posqcd3)*lumi*qcd3Xsec)/norigEventsf5)/Rw, 
error_qcd4_pos = ((sqrt(no_scale_posqcd4)*lumi*qcd4Xsec)/norigEventsf6)/Rw, 
error_qcd5_pos = ((sqrt(no_scale_posqcd5)*lumi*qcd5Xsec)/norigEventsf7)/Rw, 
error_qcd6_pos = ((sqrt(no_scale_posqcd6)*lumi*qcd6Xsec)/norigEventsf8)/Rw, 
error_qcd7_pos = ((sqrt(no_scale_posqcd7)*lumi*qcd7Xsec)/norigEventsf9)/Rw, 
error_qcd8_pos = ((sqrt(no_scale_posqcd8)*lumi*qcd8Xsec)/norigEventsf10)/Rw, 
error_qcd9_pos = ((sqrt(no_scale_posqcd9)*lumi*qcd9Xsec)/norigEventsf11)/Rw, 
error_qcd10_pos = ((sqrt(no_scale_posqcd10)*lumi*qcd10Xsec)/norigEventsf12)/Rw, 
error_qcd11_pos = ((sqrt(no_scale_posqcd11)*lumi*qcd11Xsec)/norigEventsf13)/Rw, 
error_qcd12_pos = ((sqrt(no_scale_posqcd12)*lumi*qcd12Xsec)/norigEventsf14)/Rw,
error_qcd1_neg = ((sqrt(no_scale_negqcd1)*lumi*qcd1Xsec)/norigEventsf3)/Rw, 
error_qcd2_neg = ((sqrt(no_scale_negqcd2)*lumi*qcd2Xsec)/norigEventsf4)/Rw, 
error_qcd3_neg = ((sqrt(no_scale_negqcd3)*lumi*qcd3Xsec)/norigEventsf5)/Rw, 
error_qcd4_neg = ((sqrt(no_scale_negqcd4)*lumi*qcd4Xsec)/norigEventsf6)/Rw, 
error_qcd5_neg = ((sqrt(no_scale_negqcd5)*lumi*qcd5Xsec)/norigEventsf7)/Rw, 
error_qcd6_neg = ((sqrt(no_scale_negqcd6)*lumi*qcd6Xsec)/norigEventsf8)/Rw, 
error_qcd7_neg = ((sqrt(no_scale_negqcd7)*lumi*qcd7Xsec)/norigEventsf9)/Rw, 
error_qcd8_neg = ((sqrt(no_scale_negqcd8)*lumi*qcd8Xsec)/norigEventsf10)/Rw, 
error_qcd9_neg = ((sqrt(no_scale_negqcd9)*lumi*qcd9Xsec)/norigEventsf11)/Rw, 
error_qcd10_neg = ((sqrt(no_scale_negqcd10)*lumi*qcd10Xsec)/norigEventsf12)/Rw, 
error_qcd11_neg = ((sqrt(no_scale_negqcd11)*lumi*qcd11Xsec)/norigEventsf13)/Rw, 
error_qcd12_neg = ((sqrt(no_scale_negqcd12)*lumi*qcd12Xsec)/norigEventsf14)/Rw,
error_posqcd = sqrt(pow(error_qcd1_pos,2) + pow(error_qcd2_pos,2) + pow(error_qcd3_pos,2) + pow(error_qcd4_pos,2) + pow(error_qcd5_pos,2) + pow(error_qcd6_pos,2) + pow(error_qcd7_pos,2) + pow(error_qcd8_pos,2) + pow(error_qcd9_pos,2) + pow(error_qcd10_pos,2) + pow(error_qcd11_pos,2) + pow(error_qcd12_pos,2)),
error_negqcd = sqrt(pow(error_qcd1_neg,2) + pow(error_qcd2_neg,2) + pow(error_qcd3_neg,2) + pow(error_qcd4_neg,2) + pow(error_qcd5_neg,2) + pow(error_qcd6_neg,2) + pow(error_qcd7_neg,2) + pow(error_qcd8_neg,2) + pow(error_qcd9_neg,2) + pow(error_qcd10_neg,2) + pow(error_qcd11_neg,2) + pow(error_qcd12_neg,2));
 
// W+JETS and DY
//
float w_dy_pos = pos_wev + pos_dy, 
w_dy_neg = neg_wev + neg_dy, 
error_posw1 = ((sqrt(no_scale_posw)*lumi*wXsec)/norigEventsf0)/Rw, 
error_negw1 = ((sqrt(no_scale_negw)*lumi*wXsec)/norigEventsf0)/Rw, 
error_posdy1 = ((sqrt(no_scale_posdy)*lumi*DYXsec)/norigEventsf1)/Rw, 
error_negdy1 = ((sqrt(no_scale_negdy)*lumi*DYXsec)/norigEventsf1)/Rw, 
pos_error_wdy = sqrt(pow(error_posw1,2) + pow(error_posdy1,2)), 
neg_error_wdy = sqrt(pow(error_negw1,2) + pow(error_negdy1,2)), 
expect_w = (w_dy_pos - w_dy_neg)/Rw, 
corr_w_dyw = ((pos_wev + pos_dy) - (neg_wev + neg_dy))/Rw, 
error_wdy = sqrt(pow(pos_error_wdy,2) + pow(neg_error_wdy,2)), 
stat_error_wdy = (sqrt((w_dy_pos) - (w_dy_neg)))/Rw, 
true_w_wdy = pos_wev + neg_wev;

// W+JETS, DY and SINGLE TOP ONLY
//
float w_dy_stop_pos = pos_wev + pos_dy + pos_sintop, 
w_dy_stop_neg = neg_wev + neg_dy + neg_sintop, 
error_posw2 = ((sqrt(no_scale_posw)*lumi*wXsec)/norigEventsf0)/Rw, 
error_negw2 = ((sqrt(no_scale_negw)*lumi*wXsec)/norigEventsf0)/Rw, 
error_posdy2 = ((sqrt(no_scale_posdy)*lumi*DYXsec)/norigEventsf1)/Rw, 
error_negdy2 = ((sqrt(no_scale_negdy)*lumi*DYXsec)/norigEventsf1)/Rw, 
pos_error_wdystop = sqrt(pow(error_posw2,2) + pow(error_posdy2,2)), 
neg_error_wdystop = sqrt(pow(error_negw2,2) + pow(error_negdy2,2)), 
expect_w = (pos_ev - neg_ev)/Rw, 
corr_w_wdystop = ((pos_wev + pos_dy) - (neg_wev + neg_dy))/Rw, 
error_wdystop = sqrt(pow(pos_error_wdystop,2) + pow(neg_error_wdystop,2)), 
stat_error_wdystop = (sqrt((pos_wev + pos_dy) - (neg_wev + neg_dy)))/Rw, 
true_w_wdystop = pos_wev + neg_wev;

// W+JETS, DY, SINGLE TOP and tt-bar ONLY
//
float w_dy_stop_tt_pos = pos_wev + pos_dy + pos_sintop + pos_mc, 
w_dy_stop_tt_neg = neg_wev + neg_dy + neg_sintop + neg_mc, 
error_posw3 = ((sqrt(no_scale_posw)*lumi*wXsec)/norigEventsf0)/Rw, 
error_negw3 = ((sqrt(no_scale_negw)*lumi*wXsec)/norigEventsf0)/Rw, 
error_posdy3 = ((sqrt(no_scale_posdy)*lumi*DYXsec)/norigEventsf1)/Rw, 
error_negdy3 = ((sqrt(no_scale_negdy)*lumi*DYXsec)/norigEventsf1)/Rw, 
error_postt3 = ((sqrt(no_scale_posttbar)*lumi*ttbarXsec)/norigEventsf2)/Rw, 
error_negtt3 = ((sqrt(no_scale_negttbar)*lumi*ttbarXsec)/norigEventsf2)/Rw, 
pos_error_wdystoptt = sqrt(pow(error_posw3,2) + pow(error_posdy3,2) + pow(error_postt3,2)), 
neg_error_wdystoptt = sqrt(pow(error_negw3,2) + pow(error_negdy3,2) + pow(error_negtt3,2)), 
expect_w = (pos_ev - neg_ev)/Rw, 
corr_w_wdystoptt = ((pos_wev + pos_dy + pos_mc) - (neg_wev + neg_dy + neg_mc))/Rw, 
error_wdystoptt = sqrt(pow(pos_error_wdystoptt,2) + pow(neg_error_wdystoptt,2)), 
stat_error_wdystoptt = (sqrt((pos_wev + pos_dy + pos_mc) - (neg_wev + neg_dy + neg_mc)))/Rw, 
true_w_wdystoptt = pos_wev + neg_wev;

// W+JETS, DY, SINGLE TOP, tt-bar and QCD ONLY
//
float error_posw4 = ((sqrt(no_scale_posw)*lumi*wXsec)/norigEventsf0)/Rw,
error_negw4 = ((sqrt(no_scale_negw)*lumi*wXsec)/norigEventsf0)/Rw,
error_posdy4 = ((sqrt(no_scale_posdy)*lumi*DYXsec)/norigEventsf1)/Rw,
error_negdy4 = ((sqrt(no_scale_negdy)*lumi*DYXsec)/norigEventsf1)/Rw,
error_postt4 = ((sqrt(no_scale_posttbar)*lumi*ttbarXsec)/norigEventsf2)/Rw,
error_negtt4 = ((sqrt(no_scale_negttbar)*lumi*ttbarXsec)/norigEventsf2)/Rw,

pos_error_tot = sqrt(pow(error_posw4,2) + pow(error_posdy4,2) + pow(error_postt4,2) + pow(error_posqcd,2)), 
neg_error_tot = sqrt(pow(error_negw4,2) + pow(error_negdy4,2) + pow(error_negtt4,2) + pow(error_negqcd,2)),
expect_w = (pos_ev - neg_ev)/Rw,
corr_w_tot = ((pos_wev + pos_dy + pos_mc + pos_qcd) - (neg_wev + neg_dy + neg_mc + neg_qcd))/Rw,
error_tot = sqrt(pow(pos_error_tot,2) + pow(neg_error_tot,2)),
stat_error_tot = (sqrt((pos_wev + pos_dy + pos_mc + pos_qcd) - (neg_wev + neg_dy + neg_mc + neg_qcd)))/Rw,
true_w_tot = pos_wev + neg_wev;

float a1 = ( (sqrt(error_wdy)/true_w_wdy)), b1 = (sqrt(error_wdystop)/true_w_wdystop), c1 = (sqrt(error_wdystoptt)/true_w_tot), d1 = (sqrt(error_tot)/true_w_tot);

 error_hist->SetBinContent(bin_num,corr_w_dyw/true_w_wdy);
 error_hist->GetBinError(bin_num,a1);
 error_hist->SetBinContent(bin_num+1,corr_w_wdystoptt/true_w_wdystoptt);
 error_hist->GetBinError(bin_num+1,b1);

 error_hist->SetBinContent(bin_num+2,corr_w_tot/true_w_tot);
 error_hist->GetBinError(bin_num+2,d1);

error_hist->SetMarkerStyle(20);
error_hist->SetMarkerColor(2);
cout<<"a1: "<<a1;
cout<<"\t b1: "<<b1;
cout<<"\t c1: "<<c1;
cout<<"\t d1: "<<d1<<endl;
//cout<<"true: "<<true_w_wdy<<endl;
  
        std::cout << std::string(10*7 + 2*7, '-') << "\n";
        std::cout << center(title,80)        << "\n";
        std::cout << std::string(10*7 + 2*7, '-') << "\n";

        cout << align("Process",10,true)<< " | ";
        cout << align("W+DY",10,true) << " | "; 
        cout << align("W+DY+STOP",10,true) << " | "; 
        cout << align("W+DY+STOP+ttbar",10,true) << " | ";
        cout << align("W+DY+STOP+QCD+ttbar",10,true)<< " | " << "\n";

        cout << align("Corrected",10,true)<< " | ";
        std::cout << prd(corr_w_dyw,3,10) << " | ";
	std::cout << prd(corr_w_wdystop,3,10) << " | ";
	std::cout << prd(corr_w_wdystoptt,3,15) << " | ";
	std::cout << prd(corr_w_tot,3,15) << "     | " << "\n";

	cout << align("Ratio",10,true)<< " | ";
	std::cout << prd(corr_w_dyw/true_w_wdy,3,10) << " | ";
        std::cout << prd(corr_w_wdystop/true_w_wdystop,3,10) << " | ";
        std::cout << prd(corr_w_wdystoptt/true_w_wdystoptt,3,15) << " | ";
        std::cout << prd(corr_w_tot/true_w_tot,3,15) << "     | " << "\n";

	cout << align("Pos_Error",10,true)<< " | ";
        std::cout << prd(pos_error_wdy,3,10) << " | ";
        std::cout << prd(pos_error_wdystop,3,10) << " | ";
        std::cout << prd(pos_error_wdystoptt,3,15) << " | ";
        std::cout << prd(pos_error_tot,3,15) << "     | " << "\n";

        cout << align("Neg_Error",10,true)<< " | ";
        std::cout << prd(neg_error_wdy,3,10) << " | ";
        std::cout << prd(neg_error_wdystop,3,10) << " | ";
        std::cout << prd(neg_error_wdystoptt,3,15) << " | ";
        std::cout << prd(neg_error_tot,3,15) << "     | " << "\n";

        cout << align("Error",10,true)<< " | ";
        std::cout << prd(error_wdy,3,10) << " | ";
        std::cout << prd(error_wdystop,3,10) << " | ";
        std::cout << prd(error_wdystoptt,3,15) << " | ";
        std::cout << prd(error_tot,3,15) << "     | " << "\n";

        cout << align("Stat_Error",10,true)<< " | ";
        std::cout << prd(stat_error_wdy,3,10) << " | ";
        std::cout << prd(stat_error_wdystop,3,10) << " | ";
        std::cout << prd(stat_error_wdystoptt,3,15) << " | ";
        std::cout << prd(stat_error_tot,3,15) << "     | " << "\n";


}	
