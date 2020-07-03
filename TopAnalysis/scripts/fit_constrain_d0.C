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
using namespace RooFit;

//bool GET_BIT(short x, short b) { return (x & (1<<b) ); }
float err(0.);

RooWorkspace create_workspace(bool isData=false) {
  RooWorkspace w("w");

  if(isData)
    w.factory("expr::mt('(a*mtg + b)', a[0.125373], b[154.426], mtg[173,165,183])");
  else
    w.factory("mt[0,-8,8]");
  w.factory("expr::gaus_mean('(a1*mt + a0)', a0[-34.081], a1[0.6297], mt)");//, mt[173,165,180])");
  w.factory("expr::gaus_sigma('(a3*mt + a2)', a2[14.159], a3[0.03039], mt)");
  w.factory("expr::alpha('(a5*mt + a4)', a4[1.062], a5[-0.003668], mt)");
  w.factory("expr::gamma_gamma('(a7*mt + a6)', a7[5.763], a6[-0.02088], mt)");
  w.factory("expr::gamma_beta('(a9*mt + a8)', a9[36.71], a8[0.427], mt)");
  w.factory("expr::gamma_mu('(a11*mt + a10)', a10[53.504], a11[0.3873], mt)");

  w.Print();
  return w;
}

void update_workspace(RooWorkspace *w, std::vector<std::pair<float,float>> &fit_par, std::vector<std::pair<float,float>> &fit_err, bool isData=false) {
  /*
  if(isData)
    w->factory("expr::mt('(a*mtg + b)', a[0.125373], b[154.426], mtg[173,165,183])");
  else
    w->factory("mt[0,-8,8]");
  */
  w->factory("mt[0,-8,8]");
  /*
  w->factory("expr::gaus_mean('(a1*mt + a0)', a0[72.51], a1[0.35], mt)");//, mt[173,165,180])");
  //w->factory("expr::gaus_mean('(a1*mt + a0)', a0[-34.081], a1[0.6297], mt)");//, mt[173,165,180])");
  w->factory("expr::gaus_sigma('(a3*mt + a2)', a2[20.29], a3[0.05], mt)");
  w->factory("expr::alpha('(a5*mt + a4)', a4[0.44], a5[0.003668], mt)");
  //w->factory("expr::alpha('(a5*mt + a4)', a4[1.062], a5[-0.003668], mt)");
  w->factory("expr::gamma_gamma('(a7*mt + a6)', a6[2.47], a7[0.002], mt)");
  //w->factory("expr::gamma_gamma('(a7*mt + a6)', a7[5.763], a6[-0.02088], mt)");
  w->factory("expr::gamma_beta('(a9*mt + a8)', a8[40], a9[0.00], mt)");
  //w->factory("expr::gamma_beta('(a9*mt + a8)', a9[36.71], a8[0.427], mt)");
  w->factory("expr::gamma_mu('(a11*mt + a10)', a10[11.17], a11[0.02], mt)");
  //w->factory("expr::gamma_mu('(a11*mt + a10)', a10[53.504], a11[0.3873], mt)");
  */
  w->factory(TString::Format("expr::gaus_mean('(a1*mt + a0)', a0[%f], a1[%f], mt)",fit_par[0].first, fit_par[0].second));//, mt[173,165,180])");
  w->factory(TString::Format("expr::gaus_sigma('(a3*mt + a2)', a2[%f], a3[%f], mt)",fit_par[1].first, fit_par[1].second));
  w->factory(TString::Format("expr::alpha('(a5*mt + a4)', a4[%f], a5[%f], mt)",fit_par[2].first, fit_par[2].second));
  w->factory(TString::Format("expr::gamma_gamma('(a7*mt + a6)', a6[%f], a7[%f], mt)",fit_par[3].first, fit_par[3].second));
  w->factory(TString::Format("expr::gamma_beta('(a9*mt + a8)', a8[%f], a9[%f], mt)",fit_par[4].first, fit_par[4].second));
  w->factory(TString::Format("expr::gamma_mu('(a11*mt + a10)', a10[%f], a11[%f], mt)",fit_par[5].first, fit_par[5].second));

  w->Print();
}


RooRealVar fit_constrain_d0(RooWorkspace w, std::vector<std::pair<float,float>> &fit_par, std::vector<std::pair<float,float>> &fit_err, TString mass="166v5", short flags=0b10, bool isData=false) {
  bool doBinned = GET_BIT(flags, 0);
  bool allowVary = GET_BIT(flags, 1);
  std::cout << mass << std::endl;
  std::cout << (doBinned ? "binned hist" : "unbinned tree") << std::endl;
  if(mass != "172v5" && mass.Contains("v5")) mass = "m" + mass;
  TString fname = TString("sPlot/sPlot/TopMass_"+mass+"_sPlot_d01.root");
  TFile *f = new TFile(fname);
  fname.ReplaceAll("1.root", "2.root");
  TFile *f2 = new TFile(fname);
  std::cout << fname << std::endl;
  ////TFile *f = new TFile("MC13TeV_TTJets_m"+mass+".root");
  //TChain *t = new TChain("data");
  //t->Add("Chunks/MC13TeV_TTJets_m"+mass+"_*.root");
  //TH1F *pu1 = (TH1F*) f->Get("puwgtctr_BCDEF");
  //TH1F *pu2 = (TH1F*) f->Get("puwgtctr_GH");
  //TH1F *top1 = (TH1F*) f->Get("topptwgt_BCDEF");
  //TH1F *top2 = (TH1F*) f->Get("topptwgt_GH");
  //float puSF1 = pu1->GetBinContent(1)/pu1->GetBinContent(2);
  //float puSF2 = pu2->GetBinContent(1)/pu1->GetBinContent(2);
  //float topSF1 = top1->GetBinContent(2)/top1->GetBinContent(1);
  //float topSF2 = top2->GetBinContent(2)/top1->GetBinContent(1);
  ////TTree *t = (TTree*)f->Get("data");
  //RooRealVar d0_l_mass("meson_l_mass","J/#psi+l mass", 0, 250, "GeV") ;
  //std::cout << "loaded" << std::endl;
  //if(doBinned) {
  //  TH1F *h1 = (TH1F*)gDirectory->Get("massJPsi_l_all_d0_BCDEF");
  //  TH1F *h2 = (TH1F*)gDirectory->Get("massJPsi_l_all_d0_GH");
  //  /*
  //  t->Draw("meson_l_mass>>h1(50,0,250)", "norm*sfs*puwgt*topptwgt*(d0_l_mass>0 && d0_l_mass<250 && d0_mass>3.0 && d0_mass<3.2 && epoch==1)", "goff");
  //  t->Draw("meson_l_mass>>h2(50,0,250)", "norm*sfs*puwgt*topptwgt*(d0_l_mass>0 && d0_l_mass<250 && d0_mass>3.0 && d0_mass<3.2 && epoch==2)", "goff");
  //  TH1F *h1 = (TH1F*)gDirectory->Get("h1");
  //  TH1F *h2 = (TH1F*)gDirectory->Get("h2");
  //  */
  //  h1->Scale(puSF1*topSF1*832*19716.102);
  //  h2->Scale(puSF2*topSF2*832*16146.178);
  //  TH1F *h = (TH1F*)h1->Clone("meson_l_mass");
  //  h->Add(h2);
  //  err = sqrt(h->Integral());
  //  RooDataHist dh("data", "dh", d0_l_mass, h);
  //  w.import(dh);
  //}
  //else {
  //  std::cout << "parsing tree" << std::endl;
  //  float d0lm[50],d0_mass[50],norm[50],topptwgt[50],sfs[50],puwgt[50];
  //  float l3d[50],sl3d[50];
  //  int epoch[50],meson_id[50];
  //  t->SetBranchAddress("meson_id", meson_id);
  //  t->SetBranchAddress("meson_l_mass", d0lm);
  //  t->SetBranchAddress("d0_mass", d0_mass);
  //  t->SetBranchAddress("norm", &norm);
  //  t->SetBranchAddress("topptwgt", &topptwgt);
  //  t->SetBranchAddress("sfs", sfs);
  //  t->SetBranchAddress("puwgt", puwgt);
  //  t->SetBranchAddress("epoch", epoch);
  //  t->SetBranchAddress("meson_l3d", l3d);
  //  t->SetBranchAddress("d0_sigmal3d", sl3d);
  //  RooDataSet ds("data", "ds", RooArgSet(d0_l_mass));
  //  for(int i=0; i< t->GetEntries(); i++) {
  //    t->GetEntry(i);
  //    for(int j=0; j<2; j++) {
  //      if(meson_id[j] != 443) continue;
  //      if(d0_mass[j] < 3.0) continue;
  //      if(d0_mass[j] > 3.2) continue;
  //      if(!(d0lm[j] > 0)) continue;
  //      if(d0lm[j] > 250) continue;
  //      if(l3d[j]/sl3d[j]<20.) continue;
  //      float scale = 1.;
  //      scale = sfs[j] * puwgt[j] * topptwgt[j];// * topSF * puSF;
  //      //scale = norm[j] * sfs[j] * puwgt[j] * topptwgt[j];// * topSF * puSF;
  //      if(epoch[j]==1)
  //        scale =  scale * puSF1 * topSF1 * 19716.102;
  //      else if(epoch[j]==2)
  //        scale = scale * puSF2 * topSF2 * 16146.178;
  //      else
  //        continue;
  //      d0_l_mass = d0lm[j];
  //      ds.add(RooArgSet(d0_l_mass), scale);

  //    }
  //  }
  //  TH1F *h = (TH1F*)f->Get("massJPsi_l_all_d0_BCDEF");
  //  w.import(ds);
  //  std::cout << "parse done!" << std::endl;
  //}
  ////w.import(d0_l_mass);
  //std::cout << "loading fit parameters" << std::endl;
  //int i(0);
  //for(auto & it : fit_par) {
  //  TString par = Form("a%d",(int)i);
  //  w.var(par)->setVal(it.first);
  //  w.var(par)->setConstant(!allowVary);
  //  if(allowVary) {
  //    int j = &it - &fit_par[0];
  //    w.var(par)->setRange(it.first - fit_err[j].first, it.first + fit_err[j].first);
  //  }
  //  i++;
  //  par = Form("a%d",(int)i);
  //  w.var(par)->setVal(it.second);
  //  w.var(par)->setConstant(!allowVary);
  //  if(allowVary) {
  //    int j = &it - &fit_par[0];
  //    w.var(par)->setRange(it.second - fit_err[j].second, it.second + fit_err[j].second);
  //  }
  //  i++;
  //}
  /*
  w.var("a4")->setVal(0.45);
  w.var("a4")->setRange(0.44,0.46);
  w.var("a5")->setVal(0.);
  w.var("a5")->setConstant();
  w.var("a6")->setVal(2.);
  w.var("a6")->setRange(1.9,2.1);
  w.var("a7")->setVal(0.);
  w.var("a7")->setConstant();
  */
  RooWorkspace *u = (RooWorkspace*)f->Get("w");
  RooWorkspace *u2 = (RooWorkspace*)f2->Get("w");
  std::cout << "Current WS" << std::endl;
  u->Print();
  update_workspace(u, fit_par, fit_err, isData);
  RooRealVar *d0_l_mass = (RooRealVar*)u->var("d0_l_mass");
  /*
  RooRealVar g("g","g", 2.5, 0, 10);
  RooRealVar b("b","b", 35, 30, 40);
  RooRealVar mu("mu","mu", 12, 11, 13);

  // Construct Gaussian PDF for signal
  RooRealVar mean("mean","mean", 70, 60, 90);
  RooRealVar sigma("sigma","sigma", 19, 18, 30);
  RooRealVar ngsig("ngsig","ngsignal", 100, 0, 10000);
  RooGaussian gauss("gauss","gauss", *d0_l_mass, mean, sigma);

  //  Construct Gamma PDF for signal
  RooRealVar nbsig("nbsig","nbsignal", 100, 0 , 10000);
  RooGamma gamma("gamma","gamma", *d0_l_mass, g, b, mu);

  RooRealVar alpha("alpha","alpha", 0.45, 0., 1.);
  RooAddPdf signalModel("signal model","gauss+gamma",RooArgList(gauss,gamma),RooArgList(alpha));
  std::cout << "model created" << std::endl;
  */

  /*
  w.factory("Gaussian::gauss(d0_l_mass,gaus_mean,gaus_sigma)");
  w.factory("Gamma::gamma(d0_l_mass,gamma_gamma,gamma_beta,gamma_mu)");
  w.factory("SUM::signalModel(alpha*gauss,gamma)");
  */
  u->factory("Gaussian::gauss(d0_l_mass,gaus_mean,gaus_sigma)");
  u->factory("Gamma::gamma(d0_l_mass,gamma_gamma,gamma_beta,gamma_mu)");
  u->factory("SUM::signalModel(alpha*gauss,gamma)");
  d0_l_mass->Print() ;
  std::cout << "model created" << std::endl;
  u->Print("v");
  /*
  RooPlot* frame = w.var("meson_l_mass")->frame() ;
  */
  d0_l_mass->Print() ;
  RooPlot* frame = d0_l_mass->frame() ;
  std::cout << "frame created" << std::endl;
  /*
  if(doBinned) w.data("data")->plotOn(frame);
  else w.data("data")->plotOn(frame,Binning(25));
  */
  RooDataSet ds = *(RooDataSet*)u->data("sigData");
  RooDataSet ds2 = *(RooDataSet*)u2->data("sigData");
  ds.append(ds2);
  TCanvas *c1 = setupCanvas();
  setupPad()->cd();
  std::cout << "frame plotted" << std::endl;
  //frame->Draw();
  //nll = w.pdf("signalModel")->createNLL(*w.data("data"), NumCPU(8), SumW2Error(kTRUE));
  auto meson_l_mass_d0_signal = (RooPlot*)f->Get("meson_l_mass_signal");
  TH1F *l_mass = (TH1F*)convert(meson_l_mass_d0_signal,false,0,250);
  meson_l_mass_d0_signal = (RooPlot*)f2->Get("meson_l_mass_signal");
  TH1F *l_mass2 = (TH1F*)convert(meson_l_mass_d0_signal,false,0,250);
  l_mass->Add(l_mass2);
  RooDataHist dh("dh", "dh", *d0_l_mass, Import(*l_mass));
  dh.plotOn(frame);
  //ds.plotOn(frame, RooFit::Binning(25));
  /*
  RooAbsReal *nll = u->pdf("signalModel")->createNLL(ds, NumCPU(8), SumW2Error(kTRUE));
  nll->Print();
  //nll = signalModel.createNLL(*w.data("data"), NumCPU(8), SumW2Error(kTRUE));
  std::cout << "NLL created" << std::endl;
  RooMinuit m(*nll);
  m.setPrintLevel(-1); 
  m.setPrintEvalErrors(-1);
  m.migrad();
  m.hesse();
  RooFitResult *r = m.save();
  */
  /*
  w.pdf("signalModel")->plotOn(frame);
  w.pdf("signalModel")->plotOn(frame, Components(*w.pdf("gauss")),LineStyle(kDashed),LineColor(kRed));
  w.pdf("signalModel")->plotOn(frame, Components(*w.pdf("gamma")),LineStyle(kDashed),LineColor(kBlue));
  */
  //u->pdf("signalModel")->fitTo(ds, Extended(), SumW2Error(kFALSE));
  u->pdf("signalModel")->fitTo(dh, Extended(), SumW2Error(kFALSE));
  //u->pdf("signalModel")->fitTo(dh, SumW2Error(kFALSE));
  u->pdf("signalModel")->plotOn(frame);
  u->pdf("signalModel")->plotOn(frame, Components(*u->pdf("gauss")),LineStyle(kDashed),LineColor(kRed));
  u->pdf("signalModel")->plotOn(frame, Components(*u->pdf("gamma")),LineStyle(kDashed),LineColor(kBlue));
  frame->Draw();
  //w.Print();
  //r->Print();
  //w.var("mt")->setError(err);
  if(isData) c1->SaveAs("data_fit_d0.pdf");
  if(isData) c1->SaveAs("data_fit_d0.png");
  else {
   c1->SaveAs(TString::Format("MC_fit_%s_d0.pdf", mass.Data()));
   c1->SaveAs(TString::Format("MC_fit_%s_d0.png", mass.Data()));
  }
  if(isData) {
  RooRealVar mt = *u->var("mt");
  mt.setRange(-4,6);
  RooPlot *likeframe = mt.frame();
  //nll->plotOn(likeframe,ShiftToZero(),LineColor(9));
  likeframe->SetXTitle("M_{t} - 172.5");
  likeframe->SetYTitle("-log(L/L_{max})");
  likeframe->SetMaximum(5.);
  likeframe->Draw();
  TLegend *leg_nll_data = new TLegend(0.58,0.82,0.9,0.9,NULL,"brNDC");
  leg_nll_data->SetHeader(TString::Format("M_{t} = (%3.1f #pm %1.1f) GeV", mt.getVal()+172.5, mt.getError()));
  leg_nll_data->Draw("same");
  c1->SaveAs("data_nll_d0.pdf");
  c1->SaveAs("data_nll_d0.png");
  }
  return *u->var("mt");
  //return *w.var("mt");
}

