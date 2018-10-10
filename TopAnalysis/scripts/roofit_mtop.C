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
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/convert.h"
using namespace RooFit;

bool GET_BIT(short x, short b) { return (x & (1<<b) ); }

//void roofit_mtop_BCDEFGH(TString mass="166v5", TString file="jpsi_fit.root") {
void mtop_norm(std::vector<pair<float,float>> &p, TString mass="171.5", short flags=0b00) {
  //TFile *f = new TFile("../BatchJobs/merged.root"); 
  //TFile *f = new TFile("plots/plotter_mtop_BCDEFGH.root");
  mass.ReplaceAll(".","v");
  TFile *f = new TFile("TopMass_"+mass+"_sPlot_jpsi.root");
  //f->ls(); 
  TString name = "meson_l_mass_jpsi_signal";
  //TString name = "massJPsi_l_all_jpsi_BCDEF";
  cout << "loaded " << mass << endl;
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  cout << "loaded!" << endl;

  // Declare observable x
  //RooRealVar x("J/#psi+l mass","J/#psi+l mass", 0, 250, "GeV") ;
  //RooRealVar jpsi_l_mass("J/#psi+l mass","J/#psi+l mass", 0, 250, "GeV") ;
  RooWorkspace *w = (RooWorkspace*)f->Get("w");
  RooRealVar jpsi_l_mass = *(RooRealVar*)w->var("J/#Psi+l mass");
  //RooRealVar jpsi_mass("jpsi_mass","J/#psi mass", 2.5, 3.4, "GeV") ;
  
  //cout << "creating dataset" << endl;
  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'

  cout << "plotting dataset" << endl;
  //RooPlot *frame = (RooPlot*)f->Get("meson_l_mass_jpsi_signal");
  RooPlot* frame = jpsi_l_mass.frame() ;
  //auto frame = w->var("J/#Psi+l mass")->frame();
  //cout << "converting" << endl;
  //TH1F *l_mass = (TH1F*)convert(meson_l_mass_jpsi_signal,true,0,250);
  RooDataSet ds = *(RooDataSet*)w->data("sigData");
  //ds.plotOn(frame, Binning(50));
  cout << "making dh" << endl;
  //RooDataHist dh("dh", "dh", jpsi_l_mass, Import(*l_mass));
  //RooDataHist dh("dh", "dh", jpsi_l_mass, h);
  //RooDataHist dh = *ds.binnedClone();
  cout << "making frame" << endl;

  cout << "defining variables" << endl;

  RooRealVar a0("a0", "a0", 1);
  RooRealVar a1("a1", "a1", 1);
  RooPolynomial m1("m1", "m1", jpsi_l_mass, RooArgList(a0,a1));

  // Gamma terms
  //RooRealVar g("g","g", 2,1.9,2.1);
  RooRealVar g("g","g", 2.5, 0, 10);
  //RooRealVar g("g","g", 2.5, 2.3, 2.6);
  //RooRealVar b("b","b", 32, 30, 100);
  RooRealVar b("b","b", 35, 30, 40);
  //RooRealVar b("b","b", 35, 34, 36);
  //RooRealVar mu("mu","mu", 9, 5, 14);
  RooRealVar mu("mu","mu", 12, 11, 13);

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
  signalModel.fitTo(ds);//,Extended());
  //signalModel.fitTo(dh,Extended());
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
  //dh.plotOn(frame);
  ds.plotOn(frame,Binning(25));
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
  alpha->GetYaxis()->SetRangeUser(0.3,0.7);
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
