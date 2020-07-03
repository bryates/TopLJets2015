//
// Execute the macro with
// root -l roofit_d0.C++
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
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/convert.h"
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/tdr.h"
using namespace RooFit;

bool GET_BIT(short x, short b) { return (x & (1<<b) ); }

//void roofit_mtop_BCDEFGH(TString mass="166v5", TString file="d0_fit.root") {
void mtop_norm(std::vector<pair<float,float>> &p, TString mass="172.5", short flags=0b00) {
  //TFile *f = new TFile("../BatchJobs/merged.root"); 
  //TFile *f = new TFile("plots/plotter_mtop_BCDEFGH.root");
  float mtop(mass.Atof());
  mass.ReplaceAll(".","v");
  if(mass != "172v5") mass = "m" + mass;
  TFile *f = new TFile("sPlot/sPlot/TopMass_"+mass+"_sPlot_d01.root");
  TFile *f2 = new TFile("sPlot/sPlot/TopMass_"+mass+"_sPlot_d02.root");
  //f->ls(); 
  TString name = "meson_l_mass_d0_signal";
  //TString name = "massJPsi_l_all_d0_BCDEF";
  cout << "loaded " << mass << endl;
  TCanvas *c1 = setupCanvas();
  setupPad()->cd();
  cout << "loaded!" << endl;

  // Declare observable x
  //RooRealVar x("J/#psi+l mass","J/#psi+l mass", 0, 250, "GeV") ;
  //RooRealVar d0_l_mass("J/#psi+l mass","J/#psi+l mass", 0, 250, "GeV") ;
  RooWorkspace *w = (RooWorkspace*)f->Get("w");
  RooWorkspace *w2 = (RooWorkspace*)f2->Get("w");
  //RooRealVar d0_l_mass = *(RooRealVar*)w->var("J/#Psi+l mass");
  RooRealVar d0_l_mass = *(RooRealVar*)w->var("d0_l_mass");
  //RooRealVar d0_mass("d0_mass","J/#psi mass", 2.5, 3.4, "GeV") ;
  
  //cout << "creating dataset" << endl;
  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'

  cout << "plotting dataset" << endl;
  //RooPlot *frame = (RooPlot*)f->Get("meson_l_mass_d0_signal");
  RooPlot* frame = d0_l_mass.frame() ;
  //auto frame = w->var("J/#Psi+l mass")->frame();
  //cout << "converting" << endl;
  auto meson_l_mass_d0_signal = (RooPlot*)f->Get("meson_l_mass_signal");
  TH1F *l_mass = (TH1F*)convert(meson_l_mass_d0_signal,false,0,250);
  meson_l_mass_d0_signal = (RooPlot*)f2->Get("meson_l_mass_signal");
  TH1F *l_mass2 = (TH1F*)convert(meson_l_mass_d0_signal,false,0,250);
  l_mass->Add(l_mass2);
  //RooDataSet dsn = *(RooDataSet*)w->data("dsSWeights");
  //RooDataSet ds("ds" ,"ds", &dsn, *dsn.get(), 0, "weight");
  RooDataSet ds = *(RooDataSet*)w->data("sigData");
  RooDataSet ds2 = *(RooDataSet*)w2->data("sigData");
  ds.append(ds2);
  //ds.plotOn(frame, Binning(50));
  cout << "making dh" << endl;
  RooDataHist dh("dh", "dh", d0_l_mass, Import(*l_mass));
  //RooDataHist dh("dh", "dh", d0_l_mass, h);
  //RooDataHist dh = *ds.binnedClone();
  cout << "making frame" << endl;

  cout << "defining variables" << endl;

  RooRealVar a0("a0", "a0", 1);
  RooRealVar a1("a1", "a1", 1);
  RooPolynomial m1("m1", "m1", d0_l_mass, RooArgList(a0,a1));

  // Gamma terms
  /*
  //RooRealVar g("g","g", 2,1.9,2.1);
  RooRealVar g("g","g", 2.5, 0, 10);
  //RooRealVar g("g","g", 2.5, 2.3, 2.6);
  //RooRealVar b("b","b", 32, 30, 100);
  RooRealVar b("b","b", 35, 30, 40);
  //RooRealVar b("b","b", 35, 34, 36);
  //RooRealVar mu("mu","mu", 9, 5, 14);
  RooRealVar mu("mu","mu", 12, 11, 13);
  */
  /*
  RooRealVar g("g","g", 2, 0, 20);
  RooRealVar b("b","b", 30, 0, 70);
  RooRealVar mu("mu","mu", 10, 0, 30);
  */

  //constrained versions
  //Variables
  RooConstVar mt("mt","mt", mtop-172.5);
  RooRealVar p0_gamma_bet("p0_gamma_bet","p0_gamma_bet", 27, 0, 100);//27, 0, 100);
  RooRealVar p1_gamma_bet("p1_gamma_bet","p1_gamma_bet", 0.23, 0, 10);// 0, 1);
  RooRealVar p0_gamma_gam("p0_gamma_gam","p0_gamma_gam", 2.1, 0, 100);//0, 5);
  RooRealVar p1_gamma_gam("p1_gamma_gam","p1_gamma_gam", 0.05, 0, 10);//0, 1);
  RooRealVar p0_gamma_mu("p0_gamma_mu","p0_gamma_mu", 0, 5, 50);
  RooRealVar p1_gamma_mu("p1_gamma_mu","p1_gamma_mu", 0.12, 0, 10);
  RooRealVar p0_a("p0_a","p0_a", 0.5, 0, 10);//0.45, 0.6);
  RooRealVar p1_a("p1_a","p1_a", 0.0, 0, 10);//0, 0.005);
  RooRealVar p0_gauss_mu("p0_gauss_mu","p0_gauss_mu", 60, 0, 100);//30, 75);
  RooRealVar p1_gauss_mu("p1_gauss_mu","p1_gauss_mu", 0.15, 0, 10);//0.01, 1);
  RooRealVar p0_gauss_sig("p0_gauss_sig","p0_gauss_sig", 15, 0, 100);//12, 25);
  RooRealVar p1_gauss_sig("p1_gauss_sig","p1_gauss_sig", 0.5, 0, 10);//0, 0.35);

  //Gamma
  RooFormulaVar alpha("alpha","@0+@1*@2", RooArgList(p0_a, p1_a, mt));
  //RooConstVar alpha("alpha","alpha", 0.5);
  RooFormulaVar b("b","@0+@1*@2", RooArgList(p0_gamma_bet, p1_gamma_bet, mt));
  RooFormulaVar g("g","@0+@1*@2", RooArgList(p0_gamma_gam, p1_gamma_gam, mt));
  RooFormulaVar mu("mu","@0+@1*@2", RooArgList(p0_gamma_mu, p1_gamma_mu, mt));

  //Gaussian
  RooFormulaVar mean("mean","@0+@1*@2", RooArgList(p0_gauss_mu, p1_gauss_mu, mt));
  RooFormulaVar sigma("sigma","@0+@1*@2", RooArgList(p0_gauss_sig, p1_gauss_sig, mt));

  RooGamma gamma("gamma","gamma",d0_l_mass,g,b,mu);
  RooGaussian gauss("gauss","gauss",d0_l_mass,mean,sigma);

  /*
  // Construct Gaussian PDF for signal
  //RooRealVar mean("mean","mean", 57, 40, 90);
  RooRealVar mean("mean","mean", 65, 0, 150);
  RooRealVar sigma("sigma","sigma", 16, 8, 30);
  //RooRealVar sigma("sigma","sigma", 19, 18, 30);
  //RooRealVar sigma("sigma","sigma", 19, 18.5, 19.5);
  RooRealVar ngsig("ngsig","ngsignal", 100, 0, 10000);
  RooGaussian gauss("gauss","gauss", d0_l_mass, mean, sigma);
  //RooGaussian gauss("gauss","gauss", d0_l_mass, m1, sigma);

  //  Construct Gamma PDF for signal
  RooRealVar nbsig("nbsig","nbsignal", 100, 0 , 10000);
  RooGamma gamma("gamma","gamma", d0_l_mass, g, b, mu);
  */

  cout << "defining model" << endl;
  // Construct a Gaussian+Gamma function to fit the signal component
  //RooRealVar alpha("alpha","alpha", 0.3, 0., 1.);
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
  //ds.plotOn(frame, Binning(25));
  //signalModel.fitTo(ds, PrintLevel(-1), PrintEvalErrors(-1), Warnings(kFALSE));//,Extended());
  signalModel.fitTo(dh,Extended());
  /*
  RooAbsReal *nll = signalModel.createNLL(dh, NumCPU(8), SumW2Error(kTRUE));
  RooAbsReal *nll = signalModel.createNLL(ds, NumCPU(8), SumW2Error(kTRUE));
  RooMinuit m(*nll);
  m.setPrintLevel(-1); 
  m.setPrintEvalErrors(-1);
  m.migrad();
  m.hesse();
  RooFitResult *r = m.save();
  */
  dh.plotOn(frame);
  signalModel.plotOn(frame);
  signalModel.plotOn(frame, Components(gauss),LineStyle(kDashed),LineColor(kRed));
  signalModel.plotOn(frame, Components(gamma),LineStyle(kDashed),LineColor(kBlue));
  //model.plotOn(frame, Components(expo),LineStyle(kDashed));
  //model.paramOn(frame);

  frame->Draw();
  tdr(frame, 0, false);

  mass.ReplaceAll(".","v");
  mass += "_d0";

  c1->SaveAs("MC13TeV_TTJets_"+mass+".png");
  c1->SaveAs("MC13TeV_TTJets_"+mass+".pdf");

  //x.setRange("signal model",3.0,3.2);
  //RooAbsReal *intModel = signalGauss.createIntegral(x);
  //cout << ngsig.getVal() << endl;
  //cout << nbkg.getVal() << endl;
  //cout << nsig.getVal() * intModel->getVal() << endl;

  /*
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
  */

  float mean_val(p0_gauss_mu.getValV() + p1_gauss_mu.getValV() * mt.getValV());
  float mean_err(sqrt(pow(p0_gauss_mu.getError(),2) + pow(p1_gauss_mu.getError(),2)));
  float sigma_val(p0_gauss_sig.getValV() + p1_gauss_sig.getValV() * mt.getValV());
  float sigma_err(sqrt(pow(p0_gauss_sig.getError(),2) + pow(p1_gauss_sig.getError(),2)));
  float alpha_val(p0_a.getValV() + p1_a.getValV() * mt.getValV());
  float alpha_err(sqrt(pow(p0_a.getError(),2) + pow(p1_a.getError(),2)));
  float gamma_val(p0_gamma_gam.getValV() + p1_gamma_gam.getValV() * mt.getValV());
  float gamma_err(sqrt(pow(p0_gamma_gam.getError(),2) + pow(p1_gamma_gam.getError(),2)));
  float beta_val(p0_gamma_bet.getValV() + p1_gamma_bet.getValV() * mt.getValV());
  float beta_err(sqrt(pow(p0_gamma_bet.getError(),2) + pow(p1_gamma_bet.getError(),2)));
  float mu_val(p0_gamma_mu.getValV() + p1_gamma_mu.getValV() * mt.getValV());
  float mu_err(sqrt(pow(p0_gamma_mu.getError(),2) + pow(p1_gamma_mu.getError(),2)));

  p.push_back(pair<float,float>(mean_val,mean_err));
  p.push_back(pair<float,float>(sigma_val,sigma_err));
  p.push_back(pair<float,float>(alpha_val,alpha_err));
  p.push_back(pair<float,float>(gamma_val,gamma_err));
  p.push_back(pair<float,float>(beta_val,beta_err));
  p.push_back(pair<float,float>(mu_val,mu_err));

  //x.setRange("sigma",mean.getValV()-sigma.getValV(),mean.getValV()+sigma.getValV());
  //x.setRange("sigma",0,150);
  RooAbsReal *intModel = signalModel.createIntegral(d0_l_mass);//,NormSet(x),Range("sigma"));
  RooAbsReal *intgauss = gauss.createIntegral(d0_l_mass);//,NormSet(x),Range(0,"sigma"));
  RooAbsReal *intgamma = gamma.createIntegral(d0_l_mass);//,NormSet(x));
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
                 short flags=0b00, TString file="d0_fit.root") {
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

  mean->GetYaxis()->SetRangeUser(50,80);
  sigma->GetYaxis()->SetRangeUser(20,30);
  alpha->GetYaxis()->SetRangeUser(0.3,0.5);
  gamma->GetYaxis()->SetRangeUser(1.5,3);
  beta->GetYaxis()->SetRangeUser(35,50);
  mu->GetYaxis()->SetRangeUser(6,25);
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
  gStyle->SetOptStat(0);
  /*
  for(auto &it : hists) {
  int lbin = it->FindBin(it->GetMinimumStored());
  int ubin = it->FindBin(it->GetMaximumStored());
  lbin = min(lbin, ubin);
  ubin = max(lbin, ubin);
  std::cout << "lbin" << lbin << " " << ubin << std::endl;
  it->GetYaxis()->SetRangeUser(min((int)(it->GetBinContent(lbin) - it->GetBinError(lbin)), (int)(it->GetBinContent(lbin) + it->GetBinError(lbin))),
                                 max((int)(it->GetBinContent(ubin) - it->GetBinError(ubin))+1, (int)(it->GetBinContent(ubin) + it->GetBinError(ubin))+1)+1);
  it->GetYaxis()->SetRangeUser((int)(it->GetBinContent(lbin) - it->GetBinError(lbin)), (int)(it->GetBinContent(ubin) + it->GetBinError(ubin))+1);
  }
  */
  TCanvas *c1 = setupCanvas();
  setupPad()->cd();
  for(int i = 0; i < hists.size(); i++) {
    hists[i]->Draw();
    tdr(hists[i], 0, false);
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
    name += "_d0";
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
