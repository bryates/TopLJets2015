/*
	Compare different distributions	
*/
//void Normalize() {
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TFile.h"

void Normalize() {

CompareAllHistos();
}
void CompareAllHistos()
{
  TString hnames[]={"cutflow","electroneta","electronpt","jetcsv","jeteta","jetpt","jetvtxMass","metphi","metpt","mt_3","mt_41b","mt_42b","mt_4","muondb","muondz","muoneta","muoniso","muonpt","ncsvmjets","nselelectrons","nseljets","nselmuons","nvertices","muonisov2","muonchiso","muonpuchiso","muonneuthadiso","muonphotoniso"};
  for(size_t i=0; i<sizeof(hnames)/sizeof(TString); i++)
    CompareHistos(hnames[i]);
}

void CompareHistos(TString hname)
{	
  TString baseDir("demo");

gStyle->SetOptStat(0);
//gStyle->SetMarkerStyle(0.9);

  TFile *a = new TFile("chargediso_TT_S14_PU40bx50.root");
  TH1 *h1=(TH1*) a->Get(baseDir+"/"+hname);
  h1->Sumw2();
  h1->SetDirectory(0);
//  h1->GetYaxis()->SetRangeUser(0,9000);
  a->Close();
  
  TFile *b = new TFile("chargediso_TT_PU_S14_V5.root");
  TH1 *h2=(TH1*) b->Get(baseDir+"/"+hname);
  h2->Sumw2();
  h2->SetDirectory(0);
//  h2->GetYaxis()->SetRangeUser(0,9000);
  b->Close();

  TFile *c = new TFile("chargediso_TT_PU_S14_V7.root");
  TH1 *h3=(TH1D*) c->Get(baseDir+"/"+hname);
  h3->Sumw2();
  h3->SetDirectory(0);
//  h3->GetYaxis()->SetRangeUser(0,9000);
  c->Close();
if(hname == "jetvtxMass") { 
h1->SetBinContent(1,0); 
h2->SetBinContent(1,0); 
h3->SetBinContent(1,0); 
}
  TCanvas *c1 = new TCanvas("c"+hname,"c"+hname,600,600);
  c1->Range(0,0,1,1);	
  c1->Divide(0,2);
  c1->cd();

  TPad *pad1 = new TPad("pad1", "pad1",0,0.3,1,1);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(2);
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetLeftMargin(0.12);
  pad1->SetBottomMargin(0.005);
  pad1->SetRightMargin(0.05);
  pad1->SetFrameBorderMode(0);
  pad1->Draw();
  pad1->cd();

  h1->Draw("hist");  
  h1->SetLineColor(1);
  h1->SetLineWidth(3);
  h1->SetFillStyle(0);
  
  h2->Draw("histsame");
  h2->SetLineColor(4);
  h2->SetLineWidth(2);
  h2->SetLineStyle(7);
  h2->SetMarkerStyle(20);
  h2->SetMarkerColor(4);
  h2->SetFillStyle(0);

  h3->Draw("histsame");
  h3->SetLineColor(2);
  h3->SetLineWidth(2);
  h3->SetMarkerStyle(24);
  h3->SetMarkerColor(2);
  h3->SetFillStyle(0);

  Float_t ymax=TMath::Max(TMath::Max(h2->GetMaximum(),h1->GetMaximum()),h3->GetMaximum());
  h1->GetYaxis()->SetRangeUser(0.1,ymax*1.1);

TLegend *legend=new TLegend(0.08,0.93,0.95,0.99);
   legend->SetTextFont(42);
   legend->SetTextSize(0.03);
   legend->SetNColumns(3);
   legend->SetFillStyle(0);
   legend->SetBorderSize(0);
   legend->AddEntry(h1,"MG+Pythia (40bx50)","l");
   legend->AddEntry(h2,"NormMix (S14_V5)","l");
   legend->AddEntry(h3,"PreMix (S14_V7)","l");
//   legend->AddEntry(fitFcn,"Global Fit","l");
         legend->Draw();
//
  c1->cd();
  TPad *pad2 = new TPad("pad2", "pad2",0,0,1,0.3);
  pad2->Range(-46.73368,-1.406446,353.2663,0.5209673);
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(0);
  pad2->SetGridy();
  pad2->SetTickx(1);
  pad2->SetTicky(1);
  pad2->SetLeftMargin(0.12);
  pad2->SetBottomMargin(0.3);
  pad2->SetTopMargin(0.005);
  pad2->SetRightMargin(0.05);
  pad2->SetFrameBorderMode(0);  
  pad2->Draw();
  pad2->cd();

  TH1F *h2ratio=(TH1F *) h2->Clone("h2ratio");
  h2ratio->Divide(h1);
  h2ratio->GetYaxis()->SetTitle("Ratio");
  h2ratio->GetYaxis()->SetTitleOffset(0.8);
  h2ratio->GetYaxis()->SetTitleSize(0.08);
  h2ratio->GetYaxis()->SetLabelSize(0.07);
  h2ratio->GetXaxis()->SetTitleSize(0.08);
  h2ratio->GetXaxis()->SetLabelSize(0.07);
  
  TH1F *h3ratio=(TH1F *) h3->Clone("h3ratio");
  h3ratio->Divide(h1);
  
  h2ratio->Draw("e1");
  h3ratio->Draw("e1same");
  h2ratio->GetYaxis()->SetRangeUser(0.34,1.66);

  //  h2ratio->SetMarkerColor(2);
  // h2ratio->SetMarkerStyle(kFullDotLarge);
  // h2ratio->SetMarkerSize(1);
  // h3ratio->SetMarkerColor(4);
  //h3ratio->SetMarkerStyle(kPlus);
  //h3ratio->GetYaxis()->SetRangeUser(1,2);

  c1->cd();
  c1->Modified();
  c1->Update();

  c1->SaveAs(c1->GetName()+TString(".png"));
}
//}

