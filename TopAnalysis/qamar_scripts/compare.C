/*
	To normalize the Histogram with unit ("1").....
	
*/

void compare()
{

gStyle->SetOptStat(0);
//gStyle->SetOptStat("ne");
//gStyle->SetOptStat("nemr");

	TFile *a = new TFile("/afs/cern.ch/user/q/qhassan/public/Analysis/CMSSW_7_0_6_patch1/src/UserCode/TopAnalysis/test/analysis/data/miniAOD_TT_MC.root");
 
	TH1D *h1=(TH1D*) demo->Get("jetpt");
	TH1D *h2=(TH1D*) demo->Get("bjetpt");
	TH1D *h3=(TH1D*) demo->Get("lightjetpt");
	TCanvas *c1 = new TCanvas("can1","can1",600,1000);

	h1->Draw("hist");
// 	h1->Sumw2();
	h2->Draw("histsame");
//	h2->Sumw2();
	h3->Draw("histsame");
//	h3->Sumw2();
/*		
	double it1 =h1->Integral();
	double it2 =h2->Integral();
	double it3 =h3->Integral();

	cout<< "integral of histogram 1 === " << it1 <<endl;
	cout<< "integral of histogram 2 === " << it2 <<endl;
	cout<< "integral of histogram 3 === " << it3 <<endl;

	h1->Scale(1/it1);	
	h2->Scale(1/it2);
	h3->Scale(1/it3);
*/
	h1->SetLineColor(2);
	h2->SetLineColor(3);
	h3->SetLineColor(1);
	
 Float_t ymax=TMath::Max(TMath::Max(h2->GetMaximum(),h1->GetMaximum()),h3->GetMaximum());
  h1->GetYaxis()->SetRangeUser(0.1,ymax*1.1);

TLegend *legend=new TLegend(0.08,0.93,0.95,0.99);
   legend->SetTextFont(42);
   legend->SetTextSize(0.03);
   legend->SetNColumns(3);
   legend->SetFillStyle(0);
   legend->SetBorderSize(0);
   legend->AddEntry(h1,"bjet","l");
   legend->AddEntry(h2,"bjet","l");
   legend->AddEntry(h3,"lightjet","l");
 legend->Draw();

	c1->Update();
}
