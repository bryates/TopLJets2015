#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"

#include<iostream>
#include <vector>

enum TTbarSample {AMCATNLOFXFX,AMCATNLOFXFX_HERWIG,M169v5,M175v5};
std::vector<Float_t> acceptanceEst(TTbarSample sample,Int_t njets,TString baseDir="/store/cmst3/user/psilva/LJets2015/5692fcb")
{
  std::vector<TString> files;
  if(sample==AMCATNLOFXFX)
    for(Int_t i=0; i<20; i++) files.push_back(Form("root://eoscms//eos/cms/%s/MC13TeV_TTJets_amcatnloFXFX/MergedMiniEvents_%d.root",baseDir.Data(),i));
  if(sample==AMCATNLOFXFX_HERWIG)
    for(Int_t i=0; i<11; i++) files.push_back(Form("root://eoscms//eos/cms/%s/MC13TeV_TTJets_amcatnloFXFX_herwig/MergedMiniEvents_%d.root",baseDir.Data(),i));
  if(sample==M169v5)
    for(Int_t i=0; i<6; i++) files.push_back(Form("root://eoscms//eos/cms/%s/MC13TeV_TTJets_m169v5/MergedMiniEvents_%d.root",baseDir.Data(),i));
  if(sample==M175v5)
    for(Int_t i=0; i<6; i++) files.push_back(Form("root://eoscms//eos/cms/%s/MC13TeV_TTJets_m175v5/MergedMiniEvents_%d.root",baseDir.Data(),i));
  
  //add results from all files
  TH1 *hcounter=0,*hsystlist=0;
  for(size_t i=0; i<files.size(); i++)
    {
      TFile *inF=TFile::Open(files[i]);
      if(inF==0) continue;
      if(hcounter==0) 
	{
	  hcounter = (TH1*)inF->Get(Form("analysis/fidcounter%d",njets))->Clone("hcounter");
	  hcounter->SetDirectory(0);
	  hsystlist = (TH1*)inF->Get("analysis/generator_initrwgt")->Clone("hsystlist");
	  hsystlist->SetDirectory(0);
	}
      else
	{
	  hcounter->Add((TH1*)inF->Get(Form("analysis/fidcounter%d",njets))->Clone("hcounter") );
	  hsystlist->Add( (TH1*)inF->Get("analysis/generator_initrwgt")->Clone("hsystlist"));
	}
    }

  //report on final result
  float acceptance=hcounter->GetBinContent(2)/hcounter->GetBinContent(1);
  float acceptanceUnc=hcounter->GetBinError(2)/hcounter->GetBinContent(1);
  cout << "Inclusive jet multiplicity >=" << njets << endl
       << "Acc="<< acceptance << "+/-" << acceptanceUnc << endl;

  std::vector<float>toReturn(7,0);
  toReturn[0]=acceptance;
  toReturn[1]=acceptanceUnc;

  TGraph *pdfVariations=new TGraph;
  for(int xbin=3; xbin<=hcounter->GetNbinsX(); xbin++)
    {
      float acceptanceVar=hcounter->GetBinContent(xbin)/hcounter->GetBinContent(1);
      if(acceptanceVar==0) continue;
      float deltaAcceptance=acceptanceVar-acceptance;
      TString varTitle=hsystlist->GetXaxis()->GetBinLabel(xbin);
      if(varTitle.Contains("pdfset") || varTitle.Contains("PDF set"))
	pdfVariations->SetPoint(pdfVariations->GetN(),deltaAcceptance,deltaAcceptance);
      else
	{
	  Float_t muR(1.0), muF(1.0);
	  if(varTitle.Contains("muR=0.2")) muR=2.0;
	  if(varTitle.Contains("muR=0.5")) muR=0.5;
	  if(varTitle.Contains("muF=0.2")) muF=2;
	  if(varTitle.Contains("muF=0.5")) muF=0.5;
	  if(muF==1.0 && muR==0.5) toReturn[2]=deltaAcceptance;
	  if(muF==1.0 && muR==2.0) toReturn[3]=deltaAcceptance;
	  if(muF==0.5 && muR==1.0) toReturn[4]=deltaAcceptance;
	  if(muF==2.0 && muR==1.0) toReturn[5]=deltaAcceptance;
	  cout << "\t " << muR << "\\mu_R " << muF << "\\mu_F : " << deltaAcceptance << endl;	 
	}
    }
  toReturn[6]=pdfVariations->GetRMS();
  cout << "\t PDF : " << toReturn[6] << endl;
  
  return toReturn;
}

void scanAcceptanceEstimateFor(TTbarSample sample,TString baseDir="/store/cmst3/user/psilva/LJets2015/5692fcb")
{
  TGraphErrors *grStat=new TGraphErrors();
  TGraphErrors *grPDF=new TGraphErrors();
  TGraphErrors *grTotal=new TGraphErrors();
  for(Int_t ij=1; ij<5; ij++)
    {
      std::vector<Float_t> result=acceptanceEst(sample,ij,baseDir);
      Int_t np=grStat->GetN();
      grStat->SetPoint(np,ij,result[0]);
      grStat->SetPointError(np,0.5,result[1]);
      grPDF->SetPoint(np,ij,result[0]);
      grPDF->SetPointError(np,0.5,result[6]);
      
      Float_t totalUnc(TMath::Sqrt(pow(result[1],2)
				   +pow(0.5*(fabs(result[2])+fabs(result[3])),2)
				   +pow(0.5*(fabs(result[4])+fabs(result[5])),2)
				   +pow(result[6],2)));

      grTotal->SetPoint(np,ij,result[0]);
      grTotal->SetPointError(np,0.5,totalUnc);
    }

  grStat->SetMarkerStyle(20);
  grStat->SetFillStyle(0);
  grStat->SetTitle("stat");
  grStat->SetName("stat");
  grPDF->SetMarkerStyle(1);
  grPDF->SetFillStyle(1001);
  grPDF->SetFillColor(kGray+1);
  grPDF->SetTitle("pdf");
  grPDF->SetName("pdf");
  grTotal->SetMarkerStyle(1);
  grTotal->SetFillStyle(1001);
  grTotal->SetFillColor(kGray);
  grTotal->SetTitle("total");
  grTotal->SetName("total");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *c=new TCanvas("c","c",500,500);
  grTotal->Draw("a2");
  grPDF->Draw("2");
  grStat->Draw("p");
  grTotal->GetXaxis()->SetTitle("Inclusive jet multiplicity");
  grTotal->GetYaxis()->SetTitle("Acceptance");
  grTotal->GetXaxis()->SetNdivisions(5);
  TLegend *leg=c->BuildLegend();
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetHeader("#bf{CMS} #it{simulation}");
  c->Modified();
  c->Update();
}

