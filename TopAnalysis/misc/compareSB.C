void compareSB()
{
  TString hnames[]={"jetpt_2j_leading"
		    //,"jetpt_2j_nextleading",
		    //"jetpt_3j_leading","jetpt_3j_nextleading",
		    //"jetpt_4j_leading"
		    //,"jetpt_4j_nextleading"
  };
  for(size_t i=0; i<sizeof(hnames)/sizeof(TString); i++)
    runComparison(hnames[i]);
}

void runComparison(TString hname)
{
  Float_t ttbarXsec(827.05);
  TString baseDir = "/afs/cern.ch/user/q/qhassan/public/Analysis/CMSSW_7_0_6_patch1/src/UserCode/TopAnalysis/test/analysis/data/";

  TFile* f0  = TFile::Open(baseDir+"chplots/miniAOD_MC_v2.root");
  TH1* h0jetpt = (TH1*)f0->Get(hname);
  TH1 *hcutflowTTbarf0=(TH1 *) f0->Get("cutflow");
  Float_t norigEventsTTbarf0(hcutflowTTbarf0->GetBinContent(1));
  h0jetpt->Scale(ttbarXsec/norigEventsTTbarf0);
  
  TFile* f1  = TFile::Open(baseDir+"nplots/miniAOD_MC_v2.root");
  TH1* h1fjetpt = (TH1*)f1->Get(hname);
  TH1 *hcutflowTTbarf1=(TH1 *) f1->Get("cutflow");
  Float_t norigEventsTTbarf1(hcutflowTTbarf1->GetBinContent(1));
  h1fjetpt->Scale(ttbarXsec/norigEventsTTbarf1);

 
  //QCD
  TH1F* hQCDjetpt = (TH1*)h0jetpt->Clone("hQCDjetpt");
  hQCDjetpt->Reset("ICE");
  TH1F* hQCDfjetpt = (TH1*)h1fjetpt->Clone("hQCDfjetpt");
  hQCDfjetpt->Reset("ICE");
  const unsigned int nQCD = 12;
  const char* fileNamesQCD[nQCD] = {
    "80_120", "120_170", "170_300", "300_470", "470_600", "600_800", "800_1000", "1000_1400", "1400_1800", "1800_2400", "2400_3200", "3200",
  };
  float qcdxsec[nQCD]={3000114.3,493200,120300,7475,587.1,167,28.25,8.195,0.7346,0.102,0.00644,0.000163};
  for ( int i=0; i<nQCD; ++i ) {
    TFile *inF = TFile::Open(baseDir+Form("chplots/miniAOD_QCD_%s.root", fileNamesQCD[i]));
    TH1 *hcutflow=(TH1 *) inF->Get("cutflow");
    TH1* h1jetpt = (TH1*) inF->Get(hname);
    Float_t norigEvents(hcutflow->GetBinContent(1));
    hQCDjetpt->Add(h1jetpt,qcdxsec[i]/norigEvents);

    inF = TFile::Open(baseDir+Form("nplots/miniAOD_QCD_%s.root", fileNamesQCD[i]));
    hcutflow=(TH1 *) inF->Get("cutflow");
    TH1 *h1fjetpt = (TH1*)inF->Get(hname);
    norigEvents=(hcutflow->GetBinContent(1));
    hQCDfjetpt->Add(h1fjetpt,qcdxsec[i]/norigEvents);
  }

  Int_t nbins=h0jetpt->GetNbinsX();
  cout << "charged TT:" << h0jetpt->Integral(0,nbins+1) << " QCD:" << hQCDjetpt->Integral(0,nbins+1) << endl
       << "neutral TT:" << h1fjetpt->Integral(0,nbins+1) << " QCD:" << hQCDfjetpt->Integral(0,nbins+1) << endl;

  //
  // TCanvas *c=new TCanvas("tmpc","tmpc",500,500);
  // h0jetpt->Draw("hist");
  // h1fjetpt->Draw("e1same");
  // hQCDjetpt->SetLineColor(kRed);
  // hQCDjetpt->Draw("histsame");
  // hQCDfjetpt->SetLineColor(kRed);
  // hQCDfjetpt->Draw("e1same");


  
  TGraphErrors* grpSoverB = new TGraphErrors;
  TGraphErrors* secondpt = new TGraphErrors;
  for ( int k=0, n=h0jetpt->GetNbinsX(); k<n; ++k )
  {

    //charged
    double x = h0jetpt->GetBinCenter(k+1);
    if(x>80) continue;
    double ytt = h0jetpt->Integral(k+1, n+1);
    double yqcd = hQCDjetpt->Integral(k+1, n+1);
    if ( yqcd>0 ){
      double sOverB = ytt/yqcd;
      grpSoverB->SetPoint(k, x, sOverB);
    }
    
    //neutral
    x = h1fjetpt->GetBinCenter(k+1);
    ytt = h1fjetpt->Integral(k+1, n+1);
    yqcd = hQCDfjetpt->Integral(k+1, n+1);
    if ( yqcd >0 ){
      double sOverB = ytt/yqcd;
      secondpt->SetPoint(k, x, sOverB);
    }
  }

  //show
  TCanvas* c = new TCanvas("c"+hname, "c"+hname, 500, 500);
  c->cd();

  grpSoverB->SetLineColor(2);
  grpSoverB->SetMarkerStyle(20);
  grpSoverB->SetMarkerColor(2);

  secondpt->SetLineColor(3);
  secondpt->SetMarkerStyle(20);
  secondpt->SetMarkerColor(3);

  secondpt->GetXaxis()->SetTitle("p_{T} [GeV]");
  secondpt->GetYaxis()->SetTitle("t#bar{t}/Multijets");

  secondpt->Draw("ap");
  grpSoverB->Draw("p"); 


  legend=new TLegend(0.08,0.93,0.95,0.99);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetNColumns(3);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(grpSoverB,"charged jet pt","p");
  legend->AddEntry(secondpt,"neutral jet pt","p");
  legend->Draw();

  //c->Update();
}
