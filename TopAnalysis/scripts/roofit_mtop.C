//
// Execute the macro with
// root -l roofit_jpsi.C++
//

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include <RooFit.h>
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooGamma.h"
#include "RooGlobalFunc.h"
#include <RooFitResult.h>
#include "RooWorkspace.h"
#include "RooMinuit.h"
#include <vector>
#include "TopLJets2015/TopAnalysis/interface/CharmEvent.h"
using namespace RooFit;

bool GET_BIT(short x, short b) { return (x & (1<<b) ); }

//void roofit_mtop_BCDEFGH(TString mass="166v5", TString file="jpsi_fit.root") {
void mtop_norm(std::vector<pair<float,float>> &p, TString mass="171.5", short flags=0b00) {
  //TFile *f = new TFile("../BatchJobs/merged.root"); 
  //TFile *f = new TFile("plots/plotter_mtop_BCDEFGH.root");
  mass.ReplaceAll(".","v");
  TFile *f = new TFile("MC13TeV_TTJets_m"+mass+".root");
  TChain *t = new TChain("data");
  t->Add("Chunks/MC13TeV_TTJets_m"+mass+"_*.root");
  //f->ls(); 
  TString name = "massJPsi_l_all_jpsi";
  //TString name = "massJPsi_l_all_jpsi_BCDEF";
  cout << "loading " << name+"_BCDEF mass="+mass << endl;
  TH1F *h1;
  TH1F *h2;
  if(GET_BIT(flags,0)) {
  h1 = (TH1F*)f->Get(name+"_BCDEF"); // hJpsi, hJpsiFit
  h2 = (TH1F*)f->Get(name+"_GH"); // hJpsi, hJpsiFit
  }
  cout << "loading " << name+"_GH mass="+mass << endl;
  TH1F *pu1 = (TH1F*) f->Get("puwgtctr_BCDEF");
  TH1F *pu2 = (TH1F*) f->Get("puwgtctr_GH");
  TH1F *top1 = (TH1F*) f->Get("topptwgt_BCDEF");
  TH1F *top2 = (TH1F*) f->Get("topptwgt_GH");
  float puSF1 = pu1->GetBinContent(1)/pu1->GetBinContent(2);
  float puSF2 = pu2->GetBinContent(1)/pu2->GetBinContent(2);
  float topSF1 = top1->GetBinContent(2)/top1->GetBinContent(1);
  float topSF2 = top2->GetBinContent(2)/top2->GetBinContent(1);
  cout << "PU normalization " << puSF1 << endl;
  cout << "top pT normalization " << topSF1 << endl;
  //h1->Scale(832 * 19716.102 * puSF * topSF);
  cout << "PU normalization " << puSF2 << endl;
  cout << "top pT normalization " << topSF2 << endl;
  //h2->Scale(832 * 16146.178 * puSF * topSF);
  //TH1F *h = (TH1F*)h1->Clone();
  //h->Add(h2);
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  cout << "loaded!" << endl;

  // Declare observable x
  //RooRealVar x("J/#psi+l mass","J/#psi+l mass", 0, 250, "GeV") ;
  RooRealVar jpsi_l_mass("J/#psi+l mass","J/#psi+l mass", 0, 250, "GeV") ;
  //RooRealVar jpsi_mass("jpsi_mass","J/#psi mass", 2.5, 3.4, "GeV") ;
  
  //cout << "creating dataset" << endl;
  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'

  //TTree *t = (TTree*)f->Get("data");
  //TH1F *h1 = new TH1F("h1","h1",50,0,250);
  //TH1F *h2 = new TH1F("h2","h2",50,0,250);
  //h1->Sumw2();
  //h1->SetDirectory(0);
  //h2->Sumw2();
  //h2->SetDirectory(0);
  if(!GET_BIT(flags,0)) {
  t->Draw("jpsi_l_mass>>h1(25,0,250)", "sfs*puwgt*topptwgt*(jpsi_l_mass>0 && jpsi_l_mass<250 && meson_id==443 && jpsi_mass>3.0 && jpsi_mass<3.2 && jpsi_l3d/jpsi_sigmal3d>20 && epoch==1)", "goff");
  t->Draw("jpsi_l_mass>>h2(25,0,250)", "sfs*puwgt*topptwgt*(jpsi_l_mass>0 && jpsi_l_mass<250 && meson_id==443 && jpsi_mass>3.0 && jpsi_mass<3.2  && jpsi_l3d/jpsi_sigmal3d>20 && epoch==2)", "goff");
  h1 = (TH1F*)gDirectory->Get("h1");
  h2 = (TH1F*)gDirectory->Get("h2");
  }
  //RooDataHist dh("dh", "dh", jpsi_l_mass,  h);
  //RooDataSet ds("ds", "ds", RooArgSet(jpsi_l_mass,jpsi_mass), Import(*t), Cut("jpsi_mass>3.0 && jpsi_mass<3.2"));
  //RooDataSet ds("ds", "ds", t, jpsi_l_mass);
  //CharmEvent_t ev;
  //attachToCharmEventTree(t,ev);
  float jpsilm[50],jpsi_mass[50],norm,topptwgt,sfs[50],puwgt[50];
  int epoch[50],meson_id[50];
  t->SetBranchAddress("meson_id", meson_id);
  t->SetBranchAddress("jpsi_l_mass", jpsilm);
  t->SetBranchAddress("jpsi_mass", jpsi_mass);
  t->SetBranchAddress("norm", &norm);
  t->SetBranchAddress("topptwgt", &topptwgt);
  t->SetBranchAddress("sfs", sfs);
  t->SetBranchAddress("puwgt", puwgt);
  t->SetBranchAddress("epoch", epoch);
  RooDataSet ds("ds", "ds", RooArgSet(jpsi_l_mass));
  //TH1F *h1 = new TH1F("h1","h1",100,0,250);
  //TH1F *h2 = new TH1F("h2","h2",100,0,250);
  /*
  for(int i=0; i< t->GetEntries(); i++) {
    t->GetEntry(i);
    for(int j=0; j<2; j++) {
      if(meson_id[j] != 443) continue;
      if(jpsi_mass[j] < 3.0) continue;
      if(jpsi_mass[j] > 3.2) continue;
      if(!(jpsilm[j] > 0)) continue;
      if(jpsilm[j] > 250) continue;
      float scale = 1.;
      scale = sfs[j] * puwgt[j] * topptwgt;// * topSF * puSF;
      //scale = norm * sfs[j] * puwgt[j] * topptwgt * topSF * puSF;
      if(epoch[j]==1)
        scale =  scale * 19716.102 * puSF1 * topSF1;
        //h1->Fill(jpsilm[j], scale);
      else if(epoch[j]==2)
        scale = scale * 16146.178 * puSF2 * topSF2;
        //h2->Fill(jpsilm[j], scale);
      else
        continue;
      jpsi_l_mass = jpsilm[j];
      ds.add(RooArgSet(jpsi_l_mass), scale);

    }
  }
  */
  //RooDataHist dh = *ds.binnedClone();
  h1->Scale(topSF1*puSF1);
  h2->Scale(topSF2*puSF2);
  h1->Scale(832 * 19716.102);
  h2->Scale(832 * 16146.178);
  TH1F *h = (TH1F*)h1->Clone("h");
  h->Add(h2);
  h->Draw();
  RooDataHist dh("dh", "dh", jpsi_l_mass, h);
  /*
  */

  cout << "plotting dataset" << endl;
  // Make plot of binned dataset showing Poisson error bars (RooFit default)
  RooPlot* frame = jpsi_l_mass.frame() ;
  //dh.plotOn(frame);
  //ds.plotOn(frame,Binning(50));
  //frame->Draw();

  cout << "defining variables" << endl;

  RooRealVar a0("a0", "a0", 1);
  RooRealVar a1("a1", "a1", 1);
  RooPolynomial m1("m1", "m1", jpsi_l_mass, RooArgList(a0,a1));

  // Gamma terms
  //RooRealVar g("g","g", 2,1.9,2.1);
  RooRealVar g("g","g", 2.5, 0, 10);
  //RooRealVar g("g","g", 2.5, 2.3, 2.6);
  RooRealVar b("b","b", 32, 30, 100);
  //RooRealVar b("b","b", 35, 30, 40);
  //RooRealVar b("b","b", 36, 35, 37);
  RooRealVar mu("mu","mu", 9, 5, 14);
  //RooRealVar mu("mu","mu", 11, 10, 12);

  // Construct Gaussian PDF for signal
  RooRealVar mean("mean","mean", 70, 60, 90);
  RooRealVar sigma("sigma","sigma", 19, 18, 30);
  //RooRealVar sigma("sigma","sigma", 19, 18.5, 19.5);
  RooRealVar ngsig("ngsig","ngsignal", 100, 0, 10000);
  RooGaussian gauss("gauss","gauss", jpsi_l_mass, mean, sigma);
  //RooGaussian gauss("gauss","gauss", jpsi_l_mass, m1, sigma);

  //  Construct Gamma PDF for signal
  RooRealVar nbsig("nbsig","nbsignal", 100, 0 , 10000);
  RooGamma gamma("gamma","gamma", jpsi_l_mass, g, b, mu);

  cout << "defining model" << endl;
  // Construct a Gaussian+Gamma function to fit the signal component
  RooRealVar alpha("alpha","alpha", 0.45, 0., 1.);
  //RooRealVar alpha("alpha","alpha", 0.45,0.44,0.46);
  //RooRealVar alpha("alpha","alpha", 0.4,0.39,0.41);
  RooAddPdf signalModel("signal model","gauss+gamma",RooArgList(gauss,gamma),RooArgList(alpha));

  // Construct exponential PDF to fit the bkg component
  //RooRealVar lambda("lambda", "slope", -2, -5, 5.);
  //RooExponential expo("expo", "exponential PDF", x, lambda);
  
  // Construct signal + bkg PDF
  //RooRealVar nsig("nsig","#signal events", 4000, 0, 10000) ;
  //RooRealVar nbkg("nbkg","#background events", 4000, 0, 10000) ;
  //RooAddPdf model("model","g+a", RooArgList(signalGauss, expo), RooArgList(nsig,nbkg)) ;
  cout << "fitting model" << endl;
  /*
  signalModel.fitTo(ds);//,Extended());
  */
  signalModel.fitTo(dh,Extended());
  /*
  RooAbsReal *nll = signalModel.createNLL(dh, NumCPU(8), SumW2Error(kTRUE));
  //RooAbsReal *nll = signalModel.createNLL(ds, NumCPU(8), SumW2Error(kTRUE));
  RooMinuit m(*nll);
  m.setPrintLevel(-1); 
  m.setPrintEvalErrors(-1);
  m.migrad();
  m.hesse();
  RooFitResult *r = m.save();
  */
  dh.plotOn(frame);
  //ds.plotOn(frame,Binning(50));
  signalModel.plotOn(frame);
  signalModel.plotOn(frame, Components(gauss),LineStyle(kDashed),LineColor(kRed));
  signalModel.plotOn(frame, Components(gamma),LineStyle(kDashed),LineColor(kBlue));
  //model.plotOn(frame, Components(expo),LineStyle(kDashed));
  //model.paramOn(frame);

  frame->Draw();

  mass.ReplaceAll(".","v");
  mass += "_jpsi";

  c1->SaveAs("MC13TeV_TTJets_m"+mass+".png");
  c1->SaveAs("MC13TeV_TTJets_m"+mass+".pdf");

  //x.setRange("signal model",3.0,3.2);
  //RooAbsReal *intModel = signalGauss.createIntegral(x);
  //cout << ngsig.getVal() << endl;
  //cout << nbkg.getVal() << endl;
  //cout << nsig.getVal() * intModel->getVal() << endl;

  cout << "mean," << mean.getValV() << "," << mean.getError() << endl;
  cout << "sigma," << sigma.getValV() << "," << sigma.getError() << endl;
  cout << "a," << alpha.getValV() << "," << alpha.getError() << endl;
  cout << "g," << g.getValV() << "," << g.getError() << endl;
  cout << "b," << b.getValV() << "," << b.getError() << endl;
  cout << "mu," << mu.getValV() << "," << mu.getError() << endl;

  p.push_back(pair<float,float>(mean.getValV(),mean.getError()));
  p.push_back(pair<float,float>(sigma.getValV(),sigma.getError()));
  p.push_back(pair<float,float>(alpha.getValV(),alpha.getError()));
  p.push_back(pair<float,float>(g.getValV(),g.getError()));
  p.push_back(pair<float,float>(b.getValV(),b.getError()));
  p.push_back(pair<float,float>(mu.getValV(),mu.getError()));

  //x.setRange("sigma",mean.getValV()-sigma.getValV(),mean.getValV()+sigma.getValV());
  //x.setRange("sigma",0,150);
  RooAbsReal *intModel = signalModel.createIntegral(jpsi_l_mass);//,NormSet(x),Range("sigma"));
  RooAbsReal *intgauss = gauss.createIntegral(jpsi_l_mass);//,NormSet(x),Range(0,"sigma"));
  RooAbsReal *intgamma = gamma.createIntegral(jpsi_l_mass);//,NormSet(x));
  cout << "model= " << intModel->getVal() << endl;
  cout << "Gaussian= " << intgauss->getVal() << endl;
  cout << "gamma= " << intgamma->getVal() << endl;
  /*
  RooArgSet *params = signalModel.getVariables();
  params->Print("v");
  */


  return;
}



void roofit_mtop(std::vector<float> &names,
                 std::vector<pair<float,float>> &fit_par,
                 std::vector<pair<float,float>> &fit_err,
                 short flags=0b00, TString file="jpsi_fit.root") {
  //std::vector<float> names = {166.5,171.5,172.5,173.5,175.5,178.5};
  //std::vector<float> names = {166.5,171.5,173.5,175.5,178.5};
  std::vector<std::vector<pair<float,float>>> params;
  //std::vector<pair<float,float>> fit_par;
  //std::vector<pair<float,float>> fit_err;
  for(auto & it : names) {
    std::cout << it << std::endl;
    std::vector<pair<float,float>> param;
    TString tmp_mass = Form("%.1f",it);
    //tmp_mass.ReplaceAll(".","v");
    //if(it == 171.5)
      mtop_norm(param, tmp_mass, flags);
      //return;
    //else
      //mtop(param, tmp_mass);
    params.push_back(param);
  }

  TH1F *mean  = new TH1F("mean","mean;Guassian #mu", 100, -8, 8);
  TH1F *sigma = new TH1F("sigma","sigma;Guassian #sigma", 100, -8, 8);
  TH1F *alpha = new TH1F("alpha","alpha;Guassian #alpha", 100, -8, 8);
  TH1F *gamma = new TH1F("gamma","gamma;Gamma #gamma", 100, -8, 8);
  TH1F *beta  = new TH1F("beta","beta;Gamma #beta", 100, -8, 8);
  TH1F *mu    = new TH1F("mu","mu;Gamma #mu", 100, -8, 8);

  mean->GetYaxis()->SetRangeUser(64,82);
  sigma->GetYaxis()->SetRangeUser(10,30);
  alpha->GetYaxis()->SetRangeUser(0.3,0.6);
  gamma->GetYaxis()->SetRangeUser(1,5);
  beta->GetYaxis()->SetRangeUser(20,50);
  mu->GetYaxis()->SetRangeUser(6,20);
  //int minbin = mean->FindFirstBinAbove(0);
  //cout << minbin << endl;
  //cout << mean->GetBinContent(minbin) << " " << mean->GetBinError(minbin) << endl;
  //cout << mean->GetBinContent(mean->FindFirstBinAbove(0)) << endl;
  //float min = mean->GetBinContent(minbin) - mean->GetBinError(minbin);
  //cout << min << endl;
  //mean->GetYaxis()->SetRange(mean->GetBinContent(minbin) - mean->GetBinError(minbin),78);

  std::vector<TH1F*> hists = {mean, sigma, alpha, gamma, beta, mu};

  for(int in = 0; in < names.size(); in++) {
    cout << names[in] << endl;
    for(int ip = 0; ip < hists.size(); ip++) {
      cout << in << " " << ip << endl;
      cout << hists[ip]->GetName() << endl;
      int binx;
      float &mass = names[in];
      pair<float,float> &p = params[in][ip];
      cout << p.first << " " << p.second << endl;
      binx = hists[ip]->GetXaxis()->FindBin(mass-172.5);
      hists[ip]->SetBinContent(binx, p.first);
      hists[ip]->SetBinError(binx, p.second);
    }
    cout << names[in] << endl;
  }
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  gStyle->SetOptStat(0);
  /*
  int lbin = mean->FindFirstBinAbove(0);
  int ubin = mean->FindLastBinAbove(0);
  lbin = min(lbin, ubin);
  ubin = max(lbin, ubin);
  std::cout << "lbin" << lbin << " " << ubin << std::endl;
  mean->GetYaxis()->SetRangeUser(min((int)(mean->GetBinContent(lbin) - mean->GetBinError(lbin)), (int)(mean->GetBinContent(lbin) + mean->GetBinError(lbin))),
                                 max((int)(mean->GetBinContent(ubin) - mean->GetBinError(ubin))+1, (int)(mean->GetBinContent(ubin) + mean->GetBinError(ubin))+1)+1);
  lbin = sigma->FindFirstBinAbove(0)
  ubin = sigma->FindLastBinAbove(0);
  lbin = min(lbin, ubin);
  ubin = max(lbin, ubin);
  std::cout << "lbin" << lbin << " " << ubin << std::endl;
  sigma->GetYaxis()->SetRangeUser(min((int)(sigma->GetBinContent(lbin) - sigma->GetBinError(lbin)), (int)(sigma->GetBinContent(lbin) + sigma->GetBinError(lbin))),
                                 max((int)(sigma->GetBinContent(ubin) - sigma->GetBinError(ubin))+1, (int)(sigma->GetBinContent(ubin) + sigma->GetBinError(ubin))+1)+1);
  lbin = alpha->FindFirstBinAbove(0);
  ubin = alpha->FindLastBinAbove(0);
  lbin = min(lbin, ubin);
  ubin = max(lbin, ubin);
  std::cout << "lbin" << lbin << " " << ubin << std::endl;
  alpha->GetYaxis()->SetRangeUser(min((int)(alpha->GetBinContent(lbin) - alpha->GetBinError(lbin)), (int)(alpha->GetBinContent(lbin) + alpha->GetBinError(lbin))),
                                 max((int)(alpha->GetBinContent(ubin) - alpha->GetBinError(ubin))+1, (int)(alpha->GetBinContent(ubin) + alpha->GetBinError(ubin))+1)+1);
  lbin = gamma->FindFirstBinAbove(0);
  ubin = gamma->FindLastBinAbove(0);
  lbin = min(lbin, ubin);
  ubin = max(lbin, ubin);
  std::cout << "lbin" << lbin << " " << ubin << std::endl;
  gamma->GetYaxis()->SetRangeUser(min((int)(gamma->GetBinContent(lbin) - gamma->GetBinError(lbin)), (int)(gamma->GetBinContent(lbin) + gamma->GetBinError(lbin))),
                                 max((int)(gamma->GetBinContent(ubin) - gamma->GetBinError(ubin))+1, (int)(gamma->GetBinContent(ubin) + gamma->GetBinError(ubin))+1)+1);
  lbin = beta->FindFirstBinAbove(0);
  ubin = beta->FindLastBinAbove(0);
  lbin = min(lbin, ubin);
  ubin = max(lbin, ubin);
  std::cout << "lbin" << lbin << " " << ubin << std::endl;
  beta->GetYaxis()->SetRangeUser(min((int)(beta->GetBinContent(lbin) - beta->GetBinError(lbin)), (int)(beta->GetBinContent(lbin) + beta->GetBinError(lbin))),
                                 max((int)(beta->GetBinContent(ubin) - beta->GetBinError(ubin))+1, (int)(beta->GetBinContent(ubin) + beta->GetBinError(ubin))+1)+1);
  lbin = mu->FindFirstBinAbove(0);
  ubin = mu->FindLastBinAbove(0);
  lbin = min(lbin, ubin);
  ubin = max(lbin, ubin);
  std::cout << "lbin" << lbin << " " << ubin << std::endl;
  mu->GetYaxis()->SetRangeUser(min((int)(mu->GetBinContent(lbin) - mu->GetBinError(lbin)), (int)(mu->GetBinContent(lbin) + mu->GetBinError(lbin))),
                                 max((int)(mu->GetBinContent(ubin) - mu->GetBinError(ubin))+1, (int)(mu->GetBinContent(ubin) + mu->GetBinError(ubin))+1)+1);
  */
  for(int i = 0; i < hists.size(); i++) {
    hists[i]->Draw();
    TFitResultPtr fit = hists[i]->Fit("pol1","FS");
    Double_t slope = fit->Parameter(1);
    Double_t err = fit->ParError(1);
    fit_par.push_back(std::pair<float,float>(fit->Parameter(0),fit->Parameter(1)));
    fit_err.push_back(std::pair<float,float>(fit->ParError(0),fit->ParError(1)));
    char sign = '+';
    if(fit->Parameter(1)<0) sign = '-';
    TString leg_title = TString::Format("Calibration curve : %0.2f (#pm%0.2f) %c m_{t} %0.2f (#pm%0.2f)",fit->Parameter(0),abs(fit->ParError(0)),sign,abs(slope),abs(err));
    TLegend *leg_calib = new TLegend(0.14,0.75,0.67,0.88,NULL,"brNDC");
    leg_calib->SetBorderSize(0);
    leg_calib->AddEntry(hists[i],leg_title,"lp");
    leg_calib->Draw();
    TString name = hists[i]->GetName();
    name += "_jpsi";
    c1->SaveAs("fit_"+name+".pdf");
    c1->SaveAs("fit_"+name+".png");
  }
}



/*
  // Other PDFs to test
  // Construct Crystal Ball PDF for signal
  RooRealVar cbmean("cbmean", "cbmean" , 3.1, 2.55, 3.5); 
  RooRealVar cbsigma("#sigma", "cbsigma" , 0.02, 0.001, 0.25); 
  RooRealVar ncbsig("ncbsig", "ncbsignal", 4000, 0, 10000); 
  RooRealVar n("n","cbn", 2, 0, 5);
  RooRealVar alpha("#alpha","cbalpha", 1.5, 0, 5);
  RooCBShape cball("cball", "crystal ball", x, mean, cbsigma, alpha, n);

  // Construct Breit-Wigner for signal  
  RooRealVar bwmean("bwmean", "bwmean" , 3.1, 2.55, 3.5); 
  RooRealVar bwsigma("bwsigma", "bwsigma" , 0.02, 0.00001, 10); 
  RooRealVar bwsig("bwsig", "signal", 10, 0, 1000000);
  RooBreitWigner bwgauss("bwgauss","bwgauss", x, bwmean, bwsigma);

  // Construction polynomial PDF for bkg
  RooRealVar a0("a0","a0",1,-10,10);
  RooRealVar a1("a1","a1",1,-10,10);
  RooPolynomial bkg("p1","p1", x, RooArgList(a0,a1),0);

*/
