#include <TFile.h>
#include <RooFit.h>
#include <RooWorkspace.h>
#include <RooPlot.h>
#include "convert.h"
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/tdr.C"

TH1F *smoothPlot(TH1F *fit, std::vector<std::pair<float,TString>> range, int scale=1) {
  std::cout << "fitting" << std::endl;
  std::vector<TF1*> fits;
  for(int i = 0; i < range.size(); i++) {
    float low, high=range[i].first;
    if(i==0) low = 0.2;
    else {
      low = range[i-1].first - 0.05;
    }
    std:cout << fit->GetName() << std::endl;
    std::cout << range[i].second << "\t" << low << "\t" << high << std::endl;
    //fit->Fit(range[i].second, "SQR0", "", low, high);
    fit->Fit(range[i].second, "FSMELRQ+", "", low, high);
    fits.push_back((TF1*)fit->GetFunction(range[i].second)->Clone());
  }
  /*
  fit->Fit("pol1", "SQR0", "", 0, 0.6);
  auto f1 = (TF1*)fit->GetFunction("pol1")->Clone();
  fit->Fit("pol2", "SQR0", "", 0.6, 0.85);
  auto f2 = (TF1*)fit->GetFunction("pol2")->Clone();
  fit->Fit("pol1", "SQR0", "", 0.85, 1.0);
  auto f3 = (TF1*)fit->GetFunction("pol1")->Clone();
  std::cout << "fitting done!" << std::endl;
  std::cout << "eval" << std::endl;
  */
  auto gfit = new TGraph();
  for(int i = 1; i <= fit->GetNbinsX(); i++) {
    double x = fit->GetBinCenter(i);
    /*
    if(fit->GetBinCenter(i) < 0.6) {
      std::cout << "fit1" << std::endl;
      gfit->SetPoint(i, x, f1->Eval(x));
    }
    else if(fit->GetBinCenter(i) < 0.85) {
      std::cout << "fit2" << std::endl;
      gfit->SetPoint(i, x, f2->Eval(x));
    }
    else {
      std::cout << "fit3" << std::endl;
      gfit->SetPoint(i, x, f3->Eval(x));
    }
    */
    for(int j = 0; j < range.size(); j++) {
      float high = range[j].first;
      if(fit->GetBinCenter(i) < high) {
        int s = 1;
        if(j==0) s = scale;
        gfit->SetPoint(i, x, s*fits[j]->Eval(x));
        break; //found the fit we needed
      }
    }
    std::cout << fit->GetBinLowEdge(i) << "\t" << fits.back()->Eval(x) << std::endl;
  }
  auto tSpline = new TMVA::TSpline2("spline",gfit);
  auto smooth = new TGraph();
  int num = 1000;
  for(int i = 0; i < num; i++) {
    double x = 2.*i/num;
    double val = tSpline->Eval(x);
    smooth->SetPoint(i, x, val);
  }

  auto smoothHist = (TH1F*)fit->Clone();
  for(int i = 1; i <= fit->GetNbinsX(); i++) {
    smoothHist->SetBinContent(i, smooth->Eval(fit->GetBinCenter(i)));
    std::cout << i << "\t" << fit->GetBinCenter(i) << std::endl;
    //smoothHist->SetBinError(i, smoothup->Eval(fit->GetBinCenter(i)));
  }

  /*
  delete gfit;
  delete tSpline;
  delete smooth;
  */
  return smoothHist;
}

void lerpBin(TH1F *&hist, float x, float xmin, float xmax) {
    int ibind = hist->FindBin(xmin);
    int ibinu = hist->FindBin(xmax);
    int ibin = hist->FindBin(x);
    float v0 = hist->GetBinContent(ibind);
    float v1 = hist->GetBinContent(ibinu);
    //float t =  v0 / v1;
    float t = (x - xmin) / (xmax - xmin);
    float lerp = (1 - t) * v0 + t * v1;
    hist->SetBinContent(ibin, lerp);
}

//FIXME open both epochs and add
void getHist(TString fname, TH1F *&cent, int sample, bool doLerp=true) {
  RooBinning bins(0,1.1);
  //float bin[] = {0, 0.2, 0.6, 0.7, 0.75, 0.8, 0.82, 0.84,0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0};
  std::vector<float> bin;
  if(sample==0) {
    //bin = {0, 0.2, 0.4, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    bin = {0-0.025, 0.2-0.025, 0.4-0.025, 0.55-0.025, 0.6-0.025, 0.65-0.025, 0.7-0.025, 0.75-0.025, 0.8-0.025, 0.85-0.025, 0.9-0.025, 0.95-0.025, 1.0};
    for(int i = 0; i < 15; i++) {
     bins.addBoundary(bin[i]);
    }
  }
  else if (sample==1) {
    bin = {0, 0.2, 0.4, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    for(int i = 0; i < 15; i++) {
     bins.addBoundary(bin[i]);
    }
  }
  else if(sample==2) {
    bin = {0, 0.2, 0.4, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0}; 
    for(int i = 0; i < 17; i++) {
     bins.addBoundary(bin[i]);
    }
  }
  std::cout << fname << std::endl;
  TFile *f = TFile::Open(fname);
  RooPlot *ptfrac;
  if(fname.Contains("_xb")) { // || fname.Contains("")) {
    cent = (TH1F*)f->Get("ptfrac_signal_hist")->Clone("ptfrac_signal");
    cent->SetBinContent(1,0);
    for(int ibin = 5; ibin < 9; ibin++) {
    //int ibin = 8;
    /*
    int ibin(0);
    */
    /*
    if(sample==0)
      ibin = cent->FindBin(0.75);
    */
    /*
    else if(sample==1)
      ibin = cent->FindBin(0.5);
    */
    /*
    int ibinu = ibin+1;
    int ibind = ibin-1;
    if(sample==1) {
      ibind = 4;
      ibinu = 8;
    }       
    float v0 = cent->GetBinContent(ibind);
    float v1 = cent->GetBinContent(ibinu);
    //float t =  v0 / v1;
    float t = (bin[ibin] - bin[ibind]) / (bin[ibinu] - bin[ibind]);
    float lerp = (1 - t) * v0 + t * v1;
    if(fname.Contains("_xb")) {
    std::cout << bin[ibin] << std::endl;
    std::cout << bin[ibind] << std::endl;
    std::cout << bin[ibinu] << std::endl;
    std::cout << cent->GetBinContent(ibin) << std::endl;
    std::cout << t << " " << lerp << std::endl;// << " " << v0 << " " << v1 << std::endl;
    }
    if(fname.Contains("fit"))
    cent->SetBinContent(ibin, lerp);
    */
    }
  }
  else {
  std::cout << "loading workspace" << std::endl;
  RooWorkspace *w = (RooWorkspace*)f->Get("w");
  std::cout << "loading frame" << std::endl;
  auto frame = w->var("ptfrac")->frame();
  std::cout << "loading dataset" << std::endl;
  RooDataSet *sigData = (RooDataSet*)w->data("sigData");
  std::cout << "loading plotting proxy" << std::endl;
  //if(sample!=1)
  ptfrac = sigData->plotOn(frame, RooFit::Binning(bins), RooFit::LineWidth(2));
  //else
  //ptfrac = sigData->plotOn(frame, RooFit::Binning(22), RooFit::LineWidth(2));
  std::cout << "converting histogram" << std::endl;
  if(sample==0)
  cent = (TH1F*)convert(ptfrac,false, bin);
  else if(sample==2)
  cent = (TH1F*)convert(ptfrac,false, bin);
  else
  cent = (TH1F*)convert(ptfrac,false, bin);
  }
  cent->SetDirectory(0);
  //LERP
  //JPsi
  if(doLerp) {
  if(sample==0 && fname.Contains("fit")) lerpBin(cent, 0.75, 0.7, 0.8);
  if(sample==0 && fname.Contains("fit")) lerpBin(cent, 0.85, 0.8, 0.9);
  if(sample==0 && fname.Contains("fit")) {
    int ibin = cent->FindBin(0.2);
    float val = cent->GetBinContent(ibin);
    cent->SetBinContent(ibin, 1.6 * val);
  }
  if(sample==1 && fname.Contains("fit")) {
    int ibin = cent->FindBin(0.3);
    float val = cent->GetBinContent(ibin);
    cent->SetBinContent(ibin, 0.8 * val);
  }
  if(sample==1 && fname.Contains("fit")) lerpBin(cent, 0.4, 0.5, 0.55);
  if(sample==1 && fname.Contains("fit")) lerpBin(cent, 0.55, 0.5, 0.6);
  if(sample==1 && fname.Contains("fit")) lerpBin(cent, 0.6, 0.55, 0.7);
  if(sample==1 && fname.Contains("fit")) lerpBin(cent, 0.75, 0.7, 0.8);
  if(sample==1 && fname.Contains("fit")) lerpBin(cent, 0.8, 0.75, 0.85);
  if(sample==1 && fname.Contains("fit")) lerpBin(cent, 0.85, 0.8, 0.9);
  if(sample==2 && fname.Contains("fit")) lerpBin(cent, 0.65, 0.55, 0.75);
  if(sample==2 && fname.Contains("fit")) lerpBin(cent, 0.45, 0.45, 0.55);
  }
/*
  if(sample==2 && fname.Contains("fit")) {
    int ibin = cent->FindBin(0.3);
    float val = cent->GetBinContent(ibin);
    cent->SetBinContent(ibin, 1.7 * val);
  }
  if(sample==2 && fname.Contains("fit")) lerpBin(cent, 0.55, 0.45, 0.75);
  if(sample==2 && fname.Contains("fit")) lerpBin(cent, 0.65, 0.55, 0.75);
*/
  //if(sample==2 && fname.Contains("fit")) lerpBin(cent, 0.75, 0.65, 0.85);
  /*
  if(sample==2 && fname.Contains("fit")) lerpBin(cent, 0.4, 0., 0.55);
  if(sample==2 && fname.Contains("fit")) lerpBin(cent, 0.3, 0.4, 0.2);
  if(sample==2 && fname.Contains("fit")) {
    int ibin = cent->FindBin(0.9);
    float val = cent->GetBinContent(ibin);
    cent->SetBinContent(ibin, 0.8 * val);
  }
  if(sample==2 && fname.Contains("fit")) {
    int ibin = cent->FindBin(0.95);
    float val = cent->GetBinContent(ibin);
    std::cout << val << std::endl;
    cent->SetBinContent(ibin, 0.8 * val);
  }
  if(sample==2 && fname.Contains("fit")) lerpBin(cent, 0.55, 0.4, 0.65);
  */
  //if(sample==1 && fname.Contains("fit")) lerpBin(cent, 0.4, 0.2, 0.55);
  //if(sample==1 && fname.Contains("fit")) lerpBin(cent, 0.7, 0.65, 0.75);
  //if(sample==1 && fname.Contains("fit")) lerpBin(cent, 0.75, 0.7, 0.8);
  //if(sample==1 && fname.Contains("fit")) lerpBin(cent, 0.55, 0.4, 0.6);
  /*
  */
  /*
  if(sample==0 && fname.Contains("fit")) lerpBin(cent, 0.95, 0.9, 1.0);
  //Manual adjustments to smooth bands
  if(sample==0 && fname.Contains("fit")) cent->SetBinContent(2, 50);
  if(sample==0 && fname.Contains("fit")) cent->SetBinContent(3, 120);
  if(sample==0 && fname.Contains("fit")) cent->SetBinContent(cent->FindBin(0.75), cent->GetBinContent(cent->FindBin(0.75))+15);
  //D0
  if(sample!=1) {
    //int ibin = 8;
    int ibin(0);
    if(sample==0)
      ibin = cent->FindBin(0.75);
    else if(sample==2)
      ibin = cent->FindBin(0.88);
    int ibinu = ibin+1;
    int ibind = ibin-1;
    float v0 = cent->GetBinContent(ibind);
    float v1 = cent->GetBinContent(ibinu);
    std::cout << v0 << " " << v1 << std::endl;
    //float t =  v0 / v1;
    float t = (v1 - v0) / (ibinu - ibind) / (v1 - v0);
    float lerp = (1 - t) * v0 + t * v1;
    if(fname.Contains("fit"))
    cent->SetBinContent(ibin, lerp);
  }
  */
  f->Close();
  std::cout << "file closed" << std::endl;
}

void band(int sample=0, bool fin=false, bool doLerp=false, bool jpT=false) {
TCanvas *c1 = setupCanvas();
setupPad()->cd();

std::vector<TString> samples = { "jpsi", "d0", "d0_mu_tag_mu" };
//TString fname=TString::Format("sPlot/sPlot/TopMass_172v5_Bfrag_genreco_Dz_sPlot_%s.root",samples[sample].Data());
TString fname=TString::Format("sPlot/sPlot/TopMass_172v5_Bfrag_genreco_fit_sPlot_%s.root",samples[sample].Data());
if(sample==2) fname=TString::Format("sPlot/sPlot/TopMass_172v5_Bfrag_genreco_fit_sPlot_%s.root",samples[sample].Data());
//TString fname=TString::Format("sPlot/sPlot/TopMass_172v5_Bfrag_genreco_%s_sPlot_%s.root",samples[sample].Data(),samples[sample].Data());
//fname.ReplaceAll("d0_mu_tag_mu.root","d0_mu_tag_mu.root");
int epoch=0;
if(fname.Contains("1.root")) epoch=1;
else if(fname.Contains("2.root")) epoch=2;
if(sample==2) fname.ReplaceAll(".root",TString::Format(".root"));
else fname.ReplaceAll(".root",TString::Format("_xb.root"));
std::cout << fname << std::endl;
std::map<int, char> epoch_name;
epoch_name[0] = '\0';
epoch_name[1] = 'B';
epoch_name[2] = 'G';
setupCanvas();
setupPad()->cd();

//fname.ReplaceAll(".root", "1.root");
//if(sample==1) fname.ReplaceAll("1.root", "_xb.root");
fname.ReplaceAll(TString::Format("%s",samples[sample].Data()), TString::Format("%s1",samples[sample].Data()));
TH1F *cent, *cent2, *up, *up2, *down, *down2, *data, *data2, *nom, *nom2;
getHist(fname, cent, sample);
if(sample!=1) {
fname.ReplaceAll(TString::Format("%s1",samples[sample].Data()), TString::Format("%s2",samples[sample].Data()));
getHist(fname, cent2, sample);
std::cout << "adding histograms" << std::endl;
cent->Add(cent2);
}
std::cout << "nomalizing histogram" << std::endl;
cent->Scale(1./cent->Integral());

//fname.ReplaceAll("Dz_",TString::Format("Dz_Dzuinc_"));
fname.ReplaceAll("fit_",TString::Format("fitup_"));
//fname.ReplaceAll(TString::Format("%s_",samples[sample].Data()),TString::Format("%sup_",samples[sample].Data()));
if(sample!=1)
fname.ReplaceAll(TString::Format("%s2",samples[sample].Data()), TString::Format("%s1",samples[sample].Data()));
getHist(fname, up, sample);
if(sample!=1) {
fname.ReplaceAll(TString::Format("%s1",samples[sample].Data()), TString::Format("%s2",samples[sample].Data()));
getHist(fname, up2, sample);
up->Add(up2);
}
up->Scale(1./up->Integral());

//fname.ReplaceAll("Dzuinc_",TString::Format("Dzdinc_"));
fname.ReplaceAll("fitup_",TString::Format("fitdown_"));
//fname.ReplaceAll(TString::Format("%sup_",samples[sample].Data()),TString::Format("%sdown_",samples[sample].Data()));
if(sample!=1)
fname.ReplaceAll(TString::Format("%s2",samples[sample].Data()), TString::Format("%s1",samples[sample].Data()));
getHist(fname, down, sample);
if(sample!=1) {
fname.ReplaceAll(TString::Format("%s1",samples[sample].Data()), TString::Format("%s2",samples[sample].Data()));
getHist(fname, down2, sample);
down->Add(down2);
}
down->Scale(1./down->Integral());

//fname.ReplaceAll("172v5_Bfrag_genreco_Dz_Dzdinc",TString::Format("Data"));
fname.ReplaceAll("172v5_Bfrag_genreco_fitdown",TString::Format("Data_Bfrag_genreco"));
//fname.ReplaceAll(TString::Format("172v5_Bfrag_genreco_%sdown",samples[sample].Data()), "Data");
if(sample!=1)
fname.ReplaceAll(TString::Format("%s2",samples[sample].Data()), TString::Format("%s1",samples[sample].Data()));
getHist(fname, data, sample);
if(sample!=1) {
fname.ReplaceAll(TString::Format("%s1",samples[sample].Data()), TString::Format("%s2",samples[sample].Data()));
getHist(fname, data2, sample);
data->Add(data2);
}
data->Scale(1./data->Integral());

/*
fname.ReplaceAll("down_","");
fname.ReplaceAll("172v5_Bfrag_genreco_fitdown","Data");
fname.ReplaceAll("2.root", "1.root");
getHist(fname, nom, sample);
fname.ReplaceAll("1.root", "2.root");
getHist(fname, nom2, sample);
nom->Add(nom2);
nom->Scale(1./nom->Integral());
TFile *f = TFile::Open(fname);
RooWorkspace *w = (RooWorkspace*)f->Get("w");
auto frame = w->var("ptfrac")->frame();
RooDataSet *sigData = (RooDataSet*)w->data("sigData");
RooPlot *ptfrac = sigData->plotOn(frame, RooFit::Binning(22), RooFit::LineWidth(2));
TH1F *cent = (TH1F*)convert(ptfrac,false,0,1.1);
err = sqrt(cent->Integral())/sigData->numEntries();
cent = (TH1F*)convert(ptfrac,true,0,1.1);
cent->SetDirectory(0);
f->Close();

fname.ReplaceAll(TString::Format("%s",samples[sample].Data()),TString::Format("%sup",samples[sample].Data()));
f = TFile::Open(fname);
w = (RooWorkspace*)f->Get("w");
frame = w->var("ptfrac")->frame();
sigData = (RooDataSet*)w->data("sigData");
ptfrac = sigData->plotOn(frame, RooFit::Binning(22), RooFit::LineWidth(2));
TH1F *up = (TH1F*)convert(ptfrac,true,0,1.1);
up->SetDirectory(0);
f->Close();

fname.ReplaceAll(TString::Format("%sup",samples[sample].Data()),TString::Format("%sdown",samples[sample].Data()));
f = TFile::Open(fname);
w = (RooWorkspace*)f->Get("w");
frame = w->var("ptfrac")->frame();
sigData = (RooDataSet*)w->data("sigData");
ptfrac = sigData->plotOn(frame, RooFit::Binning(22), RooFit::LineWidth(2));
TH1F *down = (TH1F*)convert(ptfrac,true,0,1.1);
down->SetDirectory(0);
f->Close();

fname.ReplaceAll("172v5_Bfrag_genreco_fitdown","Data");
f = TFile::Open(fname);
w = (RooWorkspace*)f->Get("w");
frame = w->var("ptfrac")->frame();
sigData = (RooDataSet*)w->data("sigData");
ptfrac = sigData->plotOn(frame, RooFit::Binning(22), RooFit::LineWidth(2));
TH1F *data = (TH1F*)convert(ptfrac,true,0,1.1);
data->SetDirectory(0);
f->Close();
*/

/*
fname.ReplaceAll("Data","172v5_Bfrag_genreco_lepcentral");
f = TFile::Open(fname);
w = (RooWorkspace*)f->Get("w");
frame = w->var("ptfrac")->frame();
sigData = (RooDataSet*)w->data("sigData");
ptfrac = sigData->plotOn(frame, RooFit::Binning(22), RooFit::LineWidth(2));
TH1F *lepcent = (TH1F*)convert(ptfrac,true,0,1.1);
lepcent->SetDirectory(0);
f->Close();

fname.ReplaceAll("central","up");
f = TFile::Open(fname);
w = (RooWorkspace*)f->Get("w");
frame = w->var("ptfrac")->frame();
sigData = (RooDataSet*)w->data("sigData");
ptfrac = sigData->plotOn(frame, RooFit::Binning(22), RooFit::LineWidth(2));
TH1F *lepup = (TH1F*)convert(ptfrac,true,0,1.1);
lepup->SetDirectory(0);
f->Close();

fname.ReplaceAll("up","down");
f = TFile::Open(fname);
w = (RooWorkspace*)f->Get("w");
frame = w->var("ptfrac")->frame();
sigData = (RooDataSet*)w->data("sigData");
ptfrac = sigData->plotOn(frame, RooFit::Binning(22), RooFit::LineWidth(2));
TH1F *lepdown = (TH1F*)convert(ptfrac,true,0,1.1);
lepdown->SetDirectory(0);
f->Close();

fname.ReplaceAll("lepup_","");
fname.ReplaceAll("down_","");
f = TFile::Open(fname);
w = (RooWorkspace*)f->Get("w");
frame = w->var("ptfrac")->frame();
sigData = (RooDataSet*)w->data("sigData");
ptfrac = sigData->plotOn(frame, RooFit::Binning(22), RooFit::LineWidth(2));
TH1F *nom = (TH1F*)convert(ptfrac,true,0,1.1);
//lepdown->SetDirectory(0);
*/

int n = cent->GetNbinsX();
TH1F *fit = (TH1F*)cent->Clone("Fit");
//TH1F *lep = (TH1F*)lepcent->Clone("LEP tune");
for (i=0;i<n;i++) {
   //cent->SetBinError(i, err);
   float u(1),d(1);
   u = abs(max(up->GetBinContent(i)-cent->GetBinContent(i),down->GetBinContent(i)-cent->GetBinContent(i)));
   d = abs(min(up->GetBinContent(i)-cent->GetBinContent(i),down->GetBinContent(i)-cent->GetBinContent(i)));
   fit->SetBinContent(i, cent->GetBinContent(i) + (u-d)/2);
   fit->SetBinError(i, sqrt(pow((u+d)/2,2)));
   fit->SetBinError(i, (u+d)/2);
   if(doLerp) {
   if(fit->GetBinCenter(i) < 0.5 && sample==0)
   fit->SetBinError(i, sqrt(pow((u+d),2)*100));
   if(fit->GetBinCenter(i) < 0.5 && sample==1)
   fit->SetBinError(i, sqrt(pow((u+d),2)*10));
   if(fit->GetBinCenter(i) < 0.6 && sample==2)
   fit->SetBinError(i, sqrt(pow((u+d),2)*10));
   }
   /*
   */
   /*
   //fit->SetBinError(i, sqrt(pow(cent->GetBinError(i),2) + pow((u+d)/2,2)));
   u = abs(max(lepup->GetBinContent(i)-lepcent->GetBinContent(i),lepdown->GetBinContent(i)-lepcent->GetBinContent(i)));
   d = abs(min(lepup->GetBinContent(i)-lepcent->GetBinContent(i),lepdown->GetBinContent(i)-lepcent->GetBinContent(i)));
   lep->SetBinContent(i, lepcent->GetBinContent(i) + (u-d)/2);
   lep->SetBinError(i, sqrt(pow(lepcent->GetBinError(i),2) + pow((u+d)/2,2)));
  */
}

data->SetMarkerStyle(20);
data->SetMarkerColor(kBlack);
data->SetLineColor(kBlack);
data->SetLineWidth(2);
//nom->SetLineWidth(2);
//nom->SetLineColor(kBlue);
if(fname.Contains("mu_tag"))
  data->GetXaxis()->SetTitle("(D^{0}_{#mu} #it{p}_{T} + #mu_{tag} #it{p}_{T})/#Sigma #it{p}_{T}^{ch}");
else if(fname.Contains("d0"))
  data->GetXaxis()->SetTitle("D^{0} #it{p}_{T}/#Sigma #it{p}_{T}^{ch}");
else if(fname.Contains("jpsi"))
  data->GetXaxis()->SetTitle("J/#psi #it{p}_{T}/#Sigma #it{p}_{T}^{ch}");
data->GetYaxis()->SetTitle("1/N dN/d#it{x}_{B}");
data->GetXaxis()->SetRangeUser(0,1.0);//0.975);
if(sample==1) data->GetXaxis()->SetRangeUser(0.2,1.0);
if(sample==1) fit->GetXaxis()->SetRangeUser(0.2,1.0);
data->SetMinimum(0);
//data->SetMaximum(1);
if(sample!=0) data->SetBinContent(1, -1);
data->Draw();
fit->SetFillColor(kGreen+1);
fit->SetMarkerColor(kGreen+1);
/*
lep->SetFillColor(kYellow);
lep->SetMarkerColor(kYellow);
lep->Draw("same e2");
lep->SetTitle("LEP tune");
*/
//fit->DrawCopy("hist");
//fit->Draw("same e3");
std::vector<std::pair<float,TString>> range;
if(sample==0) range = {std::make_pair<float,TString>(0.4, "pol1"), std::make_pair<float,TString>(0.85, "pol2"), std::make_pair<float,TString>(1.0, "pol1")};
if(sample==0) range = {std::make_pair<float,TString>(0.4, "pol1"), std::make_pair<float,TString>(0.85, "pol3"), std::make_pair<float,TString>(1.0, "pol1")};
else if(sample==1) range = {std::make_pair<float,TString>(0.3, "pol1"), std::make_pair<float,TString>(0.7, "pol2"), std::make_pair<float,TString>(1.0, "pol2")};
//else if(sample==1) range = {std::make_pair<float,TString>(0.5, "pol1"), std::make_pair<float,TString>(0.85, "pol2"), std::make_pair<float,TString>(1.0, "pol2")};
else if(sample==2) range = {std::make_pair<float,TString>(0.7, "pol2"), std::make_pair<float,TString>(0.9, "pol1"), std::make_pair<float,TString>(1.0, "pol1")};
    //bin = {0, 0.2, 0.4, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0};
std::cout << "fit" << std::endl;
fit->SetMaximum(.2);
//fit->Draw();
auto fitsm = (TH1F*)smoothPlot(fit, range);
//fitsm = (TH1F*)fit->Clone("fitsm");
fitsm->SetName("fit");
fitsm->SetTitle("fit");
auto fitup = (TH1F*)fit->Clone();
for(int i = 1; i <= fitup->GetNbinsX(); i++) {
  fitup->SetBinContent(i, fit->GetBinError(i));
}
int scale = 1;
if(sample==0) scale = 1;
else if(sample==1) scale = 5;
std::cout << "fit up" << std::endl;
auto fitupsm = (TH1F*)smoothPlot(fitup, range, scale);
for(int i = 1; i <= fitsm->GetNbinsX(); i++) {
  double  val = fitupsm->GetBinContent(i);
  /*
  if(sample==0) {
    if(fitsm->GetBinCenter(i) < 0.4) val *= 10;
    else if(fitsm->GetBinCenter(i) < 0.5) val *= 5;
  }
  */
  fitsm->SetBinError(i, val);
}
if(sample<2) fitsm->Draw("same e3");
else fit->Draw("same e3");
double x[2] = {0.01, 0.2};
auto *fitmask = new TH1F("fitmask", "fitmask", 1, x);
fitmask->SetBinContent(1, 0);
fitmask->SetBinError(1, 0.04);
fitmask->SetFillColor(kWhite);
fitmask->SetMarkerColor(kWhite);
data->SetBinContent(1,-1.);
fitmask->Draw("same e2");
data->Draw("AXIS same");
//g->Draw("same");
//smooth->Draw("same");
//fit->SetTitle("pp: t#bar{t} #rightarrow WWb#bar{b}");
fit->SetTitle("r_{b} fit");
data->Draw("same");
data->SetTitle("Data");
//nom->Draw("same Lhist");
//nom->SetTitle("r_{B}=0.855");
float iniy=0.83;
float dy=0.03;
float ndy=4;
TLegend *leg;
/*
if(sample==1)
leg = new TLegend(0.66, iniy-dy*ndy, 0.94, iniy+0.05);
else
*/
leg = new TLegend(0.21, iniy-dy*ndy, 0.49, iniy+0.05);
leg->SetBorderSize(1);
leg->SetFillStyle(0);      
leg->SetTextFont(43);
leg->SetTextSize(16);
leg->AddEntry( fit, fit->GetTitle(),"f");
leg->AddEntry( data, data->GetTitle(),"lp");
//leg->AddEntry( nom, nom->GetTitle(),"lp");
//leg->AddEntry( lep, lep->GetTitle(),"f");
tdr(data, epoch, fin);
gStyle->SetOptStat(0);

leg->Draw();
/*
*/
if(fin) {
c1->SaveAs(TString::Format("www/meson/tdr/fit_%s_spline_final.pdf",samples[sample].Data()));
c1->SaveAs(TString::Format("www/meson/tdr/fit_%s_spline_final.png",samples[sample].Data()));
}
else {
c1->SaveAs(TString::Format("www/meson/tdr/fit_%s_spline.pdf",samples[sample].Data()));
c1->SaveAs(TString::Format("www/meson/tdr/fit_%s_spline.png",samples[sample].Data()));
}
auto outname = TString::Format("band_%s.root", samples[sample].Data());
std::cout << outname << std::endl;
auto fout = TFile::Open(outname, "RECREATE");
data->SetDirectory(fout);
if(sample<2) fitsm->SetDirectory(fout);
else fit->SetDirectory(fout);
fitmask->SetDirectory(fout);
data->Write();
if(sample<2) fitsm->Write();
else fit->Write();
fitmask->Write();
fout->Close();
}
