void testjhgoh()
{
  TString baseDir = "/afs/cern.ch/user/q/qhassan/public/Analysis/CMSSW_7_0_6_patch1/src/UserCode/TopAnalysis/test/analysis/data/plots/";
  TFile* f0  = TFile::Open(baseDir+"miniAOD_MC_v2.root");
  TH1* h0jetpt = (TH1*)f0->Get("jetpt_2j_leading");

  Float_t ttbarXsec(827.05);
  TH1 *hcutflowTTbar=(TH1 *) f0->Get("cutflow");
  Float_t norigEventsTTbar(hcutflowTTbar->GetBinContent(1));
  h0jetpt->Scale(ttbarXsec/norigEventsTTbar);

//  f0->Close();

  TFile* f1  = TFile::Open(baseDir+"fminiAOD_MC_v2.root");
  TH1* h1fjetpt = (TH1*)f1->Get("jetpt_2j_leading");

  Float_t ttbarXsec(827.05);
  TH1 *hcutflowTTbar=(TH1 *) f1->Get("cutflow");
  Float_t norigEventsTTbar(hcutflowTTbar->GetBinContent(1));
  h1fjetpt->Scale(ttbarXsec/norigEventsTTbar);

//  f1->Close();

  TH1F* hQCDjetpt = (TH1*)h0jetpt->Clone();
  TH1F* hQCDfjetpt = (TH1*)h1fjetpt->Clone();
  hQCDjetpt->SetName("hQCDjetpt");
  hQCDfjetpt->SetName("hQCDfjetpt");
  hQCDjetpt->Reset();
  hQCDfjetpt->Reset();

  const unsigned int nQCD = 12;
  const char* fileNamesQCD[nQCD] = {
    "80_120", "120_170", "170_300", "300_470", "470_600", "600_800", "800_1000", "1000_1400", "1400_1800", "1800_2400", "2400_3200", "3200",
  };
  float qcdxsec[nQCD]={3000114.3,493200,120300,7475,587.1,167,28.25,8.195,0.7346,0.102,0.00644,0.000163};
  TFile* fQCD[nQCD];
  for ( int i=0; i<nQCD; ++i ) {
    fQCD[i] = TFile::Open(baseDir+Form("/miniAOD_QCD_%s.root", fileNamesQCD[i]));
    TH1 *hcutflow=(TH1 *) fQCD[i]->Get("cutflow");
 TH1* h1jetpt = (TH1*)fQCD[i]->Get("jetpt_2j_leading");
    Float_t norigEvents(hcutflow->GetBinContent(1));
    hQCDjetpt->Add(h1jetpt,qcdxsec[i]/norigEvents);
  }

 for ( int j=0; j<nQCD; ++j ) {
    fQCD[j] = TFile::Open(baseDir+Form("/fminiAOD_QCD_%s.root", fileNamesQCD[j]));
    TH1 *hcutflow=(TH1 *) fQCD[j]->Get("cutflow");
 TH1* h1fjetpt = (TH1*)fQCD[j]->Get("jetpt_2j_leading");
    Float_t norigEvents(hcutflow->GetBinContent(1));
    hQCDfjetpt->Add(h1fjetpt,qcdxsec[j]/norigEvents);
  }

//   TMultiGraph *mg = new TMultiGraph();
//   mg->SetTitle("Exclusion graphs");

  TGraph* grpSoverB = new TGraph;
  for ( int k=0, n=h0jetpt->GetNbinsX(); k<n; ++k )
  {
    const double x = h0jetpt->GetBinCenter(k+1);
    const double y0 = h0jetpt->Integral(k+1, n);
    const double y1 = hQCDjetpt->Integral(k+1, n);

    if ( y1 == 0 ) continue;
    const double sOverB = y0/y1;
    grpSoverB->SetPoint(k, x, sOverB);
  }

  TGraph* secondpt = new TGraph;
  for ( int l=0, n=h1fjetpt->GetNbinsX(); l<n; ++l )
  {
    const double x1 = h1fjetpt->GetBinCenter(l+1);
    const double y01 = h1fjetpt->Integral(l+1, n);
    const double y11 = hQCDfjetpt->Integral(l+1, n);
 if ( y11 == 0 ) continue;
    const double s1OverB1 = y01/y11;
    secondpt->SetPoint(l, x1, s1OverB1);
  }

  TCanvas* c = new TCanvas("c", "c", 500, 500);
  grpSoverB->Draw();
//  secondpt->Draw("[]");

//   mg->Add(grpSoverB);
//   mg->Add(secondpt);
  
//mg->Draw("AC");

  grpSoverB->SetLineColor(2);
  grpSoverB->SetMarkerStyle(20);
  grpSoverB->SetMarkerColor(2);
  
  secondpt->SetLineColor(3);
  secondpt->SetMarkerStyle(20);
  secondpt->SetMarkerColor(3);

  grpSoverB->GetXaxis()->SetTitle("2j_leading");
  grpSoverB->GetYaxis()->SetTitle("S/B");

  legend=new TLegend(0.08,0.93,0.95,0.99);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetNColumns(3);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
//  legend->AddEntry(grpSoverB,"jet pt","ap");
//  legend->AddEntry(secondpt,"charged jet pt","p");
  legend->Draw();
  
  c->Update();
}
