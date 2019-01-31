//
// Execute the macro with
// root -l roofit_d0.C++
//

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLegend.h"
#include <RooFit.h>
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooStats/SPlot.h"
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
#include "/afs/cern.ch/user/b/byates/TopAnalysis/interface/CharmEvent.h"
  #include "/afs/cern.ch/user/b/byates/TopAnalysis/src/CharmEvent.cc"
  //#include "convert.h"
//#include "TopAnalysis/interface/CharmEvent.h"
using namespace RooFit;
using namespace RooStats;

bool GET_BIT(short x, short b) { return (x & (1<<b) ); }

//void roofit_mtop_BCDEFGH(TString mass="166v5", TString file="d0_fit.root") {
//void mtop_norm(std::vector<pair<float,float>> &p, TString mass="171.5", short flags=0b00) {
void splot_d0_mu_tag(TString mass="172.5", bool isData=false, TString fragWeight="", int ep=0, bool jpT=false) {
  RooWorkspace w("w",mass);
  float wind(0.045);
/*
void splot_d0_mu(RooWorkspace &w, TString mass="172.5", bool isData=false) {
*/
  //TFile *f = new TFile("../BatchJobs/merged.root"); 
  //TFile *f = new TFile("plots/plotter_mtop_BCDEFGH.root");
  TString syst("");
  //TString dir("/eos/cms/store/user/byates/top18012/Chunks/");
  TString dir("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/Chunks/");
  mass.ReplaceAll(".","v");
  if(mass.Contains("v5") && mass != "172v5") mass = "m" + mass;
  if(mass.Contains("sr") || mass.Contains("erdON") || mass.Contains("Move") || mass.Contains("ue") || mass.Contains("hdamp")) { //ISR,FSR,CR,UE
    syst = mass;
    //mass = "172.5";
    std::cout << "Processing systematics " << TString::Format("MC13TeV_TTJets_%s",syst.Data()) << std::endl;
  }
  else if(mass.Contains("m1")) { //top mass tunes
    syst = mass;
    //mass = "172.5";
    std::cout << "Processing systematics " << TString::Format("MC13TeV_%s",syst.Data()) << std::endl;
  }
  TFile *fin;
  TGraph *g;
  //std::vector<TGraph*> fwgt;
  TH1F *tuneWgt = new TH1F("tuneWgt","tuneWgt",2,0,2);
  TFile *f = new TFile(dir+"../MC13TeV_TTJets_powheg.root"); //open a randome file to get correction histograms
  //TString dir("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/ndau/Chunks/");
  //TFile *f = new TFile(dir+"MC13TeV_TTJets_m"+mass+"_0.root"); //open a randome file to get correction histograms
  //std::vector<RooRealVar> frgWgt;
  TChain *data = new TChain("data");
  if(isData) {
    data->Add(dir+"Data13TeV_*");
    mass="Data";
  }
  //else data->Add("Chunks/MC13TeV_TTJets_m"+mass+"_*.root");
  else {
    if(mass.Contains("TRK") || mass.Contains("TRIGGER") || mass.Contains("LEP") ||mass.Contains("PU") || mass.Contains("PI") || mass.Contains("JER")) {  // SFs
      TString tmp(mass);
      tmp.ReplaceAll("_","/");
      //dir = TString("/eos/cms/store/user/byates/top18012/" + tmp + "/Chunks/");
      dir = TString("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/" + tmp + "/Chunks/");
      std::cout << dir << std::endl;
    }
    //std::vector<TString> mcSamples = { "MC13TeV_TTJets_powheg" };
    std::vector<TString> mcSamples;
    if(syst.Length()>0) mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW",TString::Format("MC13TeV_TTJets_%s",syst.Data()),"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    //Gluon Move and QCD only!
    //if(syst.Length()>0) mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW",TString::Format("MC13TeV_%s",syst.Data()),"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    //else mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_TTJets_m"+mass,"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    else mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_TTJets_powheg","MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    if(mass.Contains("QCD") || mass.Contains("GluonMove")) mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW",TString::Format("MC13TeV_%s*",syst.Data()),"MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
    for(auto & it : mcSamples) {
      //TString mcname = "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/Chunks/";
      TString mcname(dir);
      mcname += it;
      mcname += "_*.root";
      std::cout << "Adding " << it;
      int num = data->Add(mcname);
      std::cout << " (" << num << " files)" << std::endl;
    }
    /*
    std::vector<TString> tunes = {"_down", "_ddown", "_dddown", "", "_cccentral", "_ccentral", "_central", "_uuup", "_uup", "_up" };
    for(auto &it : tunes) {
      fin = TFile::Open("/afs/cern.ch/user/b/byates/TopAnalysis/data/era2016/bfragweights.root");
      fwgt.push_back((TGraph*)fin->Get(it+"Frag"));
      frgWgt.push_back(RooRealVar(it+"Frag", it+"Frag", 1., 0, 2.));
    }
    */
    if(fragWeight.Length() > 0) {
      if(fragWeight.Contains("lep")) {
        fin = TFile::Open("/eos/cms/store/user/byates/top18012/bfragweights_Markus.root");
      }
      else fin = TFile::Open("/eos/cms/store/user/byates/top18012/bfragweights.root");
      mass += "_" + fragWeight;
      fragWeight.ReplaceAll("lep","");
      g = (TGraph*)fin->Get(fragWeight+"Frag");
      std::cout << g->GetName() << std::endl;
    }
  }
  TString fUrl("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_"+mass+"_sPlot_d0_mu_tag_mu.root");
  //TString fUrl("TopMass_"+mass+"_sPlot_d0_mu_tag_mu.root");
  if(ep>0) fUrl.ReplaceAll(".root",TString::Format("%d.root",ep));
  if(jpT) fUrl.ReplaceAll(".root","_jpT.root");
  std::cout << "creating file: "  << fUrl<< std::endl;
  TFile *fout = new TFile(fUrl,"RECREATE");
  //data->Add("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/Chunks/MC13TeV_W*Jets_*.root");
  //disabled top pT reweighting
  //data->Draw("d0_mass>>h(30,1.7,2.0)","norm*sfs*puwgt*(meson_id==42113 && d0_l3d/d0_sigmal3d>10 && HT>180 && d0_sigmal3d>2E-4)");// && j_hadflav[d0_j]==5 && d0_pi_mother==421 && d0_k_mother==421)","goff");
  //data->Draw("d0_mass>>h(30,1.7,2.0)","norm*sfs*puwgt*topptwgt*(meson_id==42113 && d0_l3d/d0_sigmal3d>10 && HT>180 && d0_sigmal3d>2E-4)");// && j_hadflav[d0_j]==5 && d0_pi_mother==421 && d0_k_mother==421)","goff");
  
  CharmEvent_t ev;
  attachToCharmEventTree(data,ev);
  TH1F *h = (TH1F*)gDirectory->Get("h");
  //if(!isData) h->Scale(35864);
  //if(!isData) h->Scale(832*35864);
  //f->ls(); 
  TString name = "massD0_mu_tag_l_all_d0";
  //TString name = "massJPsi_l_all_d0_BCDEF";
  TH1F *pu1 = (TH1F*) f->Get("puwgtctr_BCDEF");
  TH1F *pu2 = (TH1F*) f->Get("puwgtctr_GH");
  TH1F *top1 = (TH1F*) f->Get("topptwgt_BCDEF");
  TH1F *top2 = (TH1F*) f->Get("topptwgt_GH");
  float puSF1 = pu1->GetBinContent(1)/pu1->GetBinContent(2);
  float puSF2 = pu2->GetBinContent(1)/pu2->GetBinContent(2);
  float topSF1 = top1->GetBinContent(2)/top1->GetBinContent(1);
  float topSF2 = top2->GetBinContent(2)/top2->GetBinContent(1);
  pu1->SetDirectory(0);
  pu2->SetDirectory(0);
  top1->SetDirectory(0);
  top2->SetDirectory(0);
  f->Close();
  if(!isData) {
  cout << "PU normalization " << puSF1 << endl;
  cout << "top pT normalization " << topSF1 << endl;
  cout << "PU normalization " << puSF2 << endl;
  cout << "top pT normalization " << topSF2 << endl;
  }
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  cout << "loaded!" << endl;

  // Declare observable x
  RooRealVar d0_mass("d0_mass","M_{K#pi}", 1.7, 2, "GeV") ;
  RooRealVar weight("weight","weight",1.,0.,2.);
  //RooRealVar weight("weight","weight",1.,0.,36000.);
  RooRealVar meson_l_mass("d0_l_mass","D^{0}+l mass", 0, 250, "GeV") ;
  RooRealVar ptfrac("ptfrac","(D^{0} p_{T} + #mu p_{T}) / #Sigma_{ch} p_{T}", 0, 1.1, "") ;
  RooRealVar d0_pt("d0_pt","D^{0}_{#mu} p_{T}", 0, 250, "GeV");
  RooRealVar j_pt_ch("j_pt_ch","#Sigma_{ch} p_{T} charged", 0, 400, "GeV");
  RooRealVar j_pt("j_pt","j p_{T}", 0, 400, "GeV");
  RooRealVar epoch("epoch","epoch",1,2);
  RooRealVar tuneW = RooRealVar("tuneW", "tuneW", 1., 0, 2.);
  /*
  RooArgSet args(d0_mass,ptfrac,meson_l_mass,weight,d0_pt,j_pt_ch,epoch,tuneW);
  for(auto &it : frgWgt)
    args.add(it);
  args.Print();
  */
  
  //cout << "creating dataset" << endl;
  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'

  //RooDataSet dsn("dsn", "dsn", RooArgSet(meson_id,d0_mass,ptfrac,meson_l_mass,weight), Import(*data), Cut("meson_id==42113"));
  RooDataSet dsn("dsn", "dsn", RooArgSet(d0_mass,ptfrac,meson_l_mass,weight,d0_pt,epoch,j_pt_ch,tuneW));
  //RooDataSet dsn("dsn", "dsn", args);
  std::cout << "Total events: " << data->GetEntries() << std::endl;
  for(int i=0; i < data->GetEntries(); i++) {
    ev = {};
    data->GetEntry(i);
    for(int j=0; j<ev.nmeson; j++) {
      //int j(0);
      if(ev.meson_id[j] != 42113) continue;
      if(ev.d0_mass[j] < 1.7) continue;
      if(ev.d0_mass[j] > 2.0) continue;
      if(ep>0 && ev.epoch[j] != ep) continue;
      if(!jpT) {
        //if(ev.j_pt[j]>150) continue;
        if(ev.j_pt_charged[j]>75) continue;
      }
      /*
      else
        if(ev.j_pt[j]>50) continue;
      if(ev.d0_l3d[j]/ev.d0_sigmal3d[j]<10) continue;
      if(ev.d0_sigmal3d[j]<2E-4) continue;
      if(ev.ht<180) continue;
      */
      /*
      if(!ev.isData) {
        if(ev.j_hadflav[ev.d0_j[j]]!=5) continue;
        if(ev.d0_pi_mother[j]!=421) continue;
        if(ev.d0_k_mother[j]!=421) continue;
      }
      */
      //if(j_pt[(int)d0_j[j]]<30) continue;
      /*
      if(!(mesonlm[j] > 0)) continue;
      if(mesonlm[j] > 250) continue;
      */
      float scale = 1.;
      tuneW.setVal(1.);
      if(!isData) {
        scale = ev.norm * ev.xsec * ev.sfs[j] * ev.puwgt[j] * ev.topptwgt;// * topSF * puSF;
        //scale = norm * sfs[j] * puwgt[j] * topptwgt * topSF * puSF;
        if(ev.epoch[j]==1) {
          scale =  scale * 19712.86 * puSF1 * topSF1;
          //h1->Fill(mesonlm[j], scale);
        }
        else if(ev.epoch[j]==2) {
          scale = scale * 16146.178 * puSF2 * topSF2;
          //h2->Fill(mesonlm[j], scale);
        }
        else
          continue;
        tuneWgt->Fill(0.,1.0);
        if(fragWeight.Length() > 0) {
          double frac(0.);
          if(jpT) frac=ev.d0_mu_tag_mu_pt[j]/ev.j_pt[j];
          else frac=ev.d0_mu_tag_mu_pt[j]/ev.j_pt_charged[j];

          if(frac>1.) frac = 1.;
          if(frac<0.1) frac = 0.1;
          /*
          for(auto &it : fwgt)
            args.setRealValue(it->GetName(),scale*it->Eval(frac));
          */
          scale *= g->Eval(frac);
          tuneWgt->Fill(1.,g->Eval(frac));
        }
        else tuneWgt->Fill(1., 1.0); //unit weights to make re-weighting easier to loop over
      }
      epoch.setVal(ev.epoch[j]);
      //d0_mass = ev.d0_mass[j];
      d0_mass.setVal(ev.d0_mass[j]);
      if(ev.d0_mass[j]>(1.864-wind) && ev.d0_mass[j]<(1.864+wind)) { //symmetric around PDG mass = 1.864
      //if(ev.d0_mass[j]>(1.864-0.036) && ev.d0_mass[j]<(1.864+0.036)) { //symmetric around PDG mass = 1.864
      //if(ev.d0_mass[j]>1.83 && ev.d0_mass[j]<1.9 && ev.d0_mu_tag_mu_pt[j]>0) {
      //if(ev.d0_mass[j]>1.8 && ev.d0_mass[j]<1.93) {
      if(!jpT) ptfrac.setVal(ev.d0_mu_tag_mu_pt[j]/ev.j_pt_charged[j]);
      else ptfrac.setVal(ev.d0_mu_tag_mu_pt[j]/ev.j_pt[j]);
      d0_pt.setVal(ev.d0_mu_tag_mu_pt[j]);
      j_pt_ch.setVal(ev.j_pt_charged[j]);
      j_pt.setVal(ev.j_pt[j]);
      //meson_l_mass = ev.d0_l_mass[j];
      meson_l_mass.setVal(ev.d0_l_mass[j]);
      }
      weight.setVal(scale);
      tuneW.setVal(scale);// * tuneW.getVal());
      dsn.add(RooArgSet(d0_mass,meson_l_mass,ptfrac,weight,d0_pt,j_pt_ch,epoch,tuneW));
      dsn.add(RooArgSet(j_pt));
      //dsn.add(RooArgSet(d0_mass,meson_l_mass,ptfrac,weight), scale);
      //within D^0 mass peak (1.864 +/- 0.05)
      //if(ev.d0_l_mass[j] > 1.8 && ev.d0_l_mass[j] < 1.93)
        //dsn.add(RooArgSet(meson_l_mass,ptfrac), scale);

    }
  }
  delete data;
  std::cout << "Done loading tree" << std::endl;

  /*
  tuneW.setConstant();
  if(tuneWgt->GetBinContent(2) > 0) tuneW.setVal(tuneWgt->GetBinContent(2)/tuneWgt->GetBinContent(1));
  RooDataSet tmp("tmp", "tmp", &dsn, *dsn.get(), 0, "weight");
  RooDataSet ds("ds", "ds", &tmp, *tmp.get(), 0, "tuneW");
  */
  /*
  */
  //very inefficient second loop to adjust normalization
  std::cout << "Re-weighting based on norm weight" << std::endl;
  if(!isData) {
  double tuneShape = 1.;
  if(tuneWgt->GetBinContent(2) > 0) tuneShape = tuneWgt->GetBinContent(1)/tuneWgt->GetBinContent(2);
  std::cout << "Norm weight: " << tuneShape << std::endl;
  for(int i = 0; i < dsn.numEntries(); i++) {
    double tmpW = dsn.get(i)->getRealValue("tuneW");
    tuneW.setVal(tmpW * tuneShape);
    weight.setVal(tuneW.getVal());
    //std::cout << tuneW.getVal() << " " << weight.getVal() << std::endl;;
  }
  std::cout << "Re-weighting done!" << std::endl;
  }
  /*
  tmp.Print();
  tmp.addColumn(tuneW, tuneShape);
  tmp.Print();
  //RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, "weight");
  RooDataSet ds("ds", "ds", &tmp, *tmp.get(), 0, "tuneW");
  RooDataSet tmp("tmp", "tmp", &dsn, *dsn.get(), 0, "weight");
  tmp.add(tuneW, tuneShape);
  */
  //RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, "weight");
  //RooRealVar *wgt = (RooRealVar*)dsn.addColumn(tuneW, tuneShape);
  RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, "tuneW");
  //RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, "weight");
  //RooDataSet ds("ds", "ds", &dsn, *dsn.get(), 0, wgt->GetName());
  dsn.Print("v");
  ds.Print("v");
  ds.Print();
  std::cout << ds.weight() << std::endl;
  cout << "plotting dataset" << endl;
  // Make plot of binned dataset showing Poisson error bars (RooFit default)

  cout << "defining variables" << endl;

  // Mean of the J/psi mass peak
  /*
  RooRealVar mean("mean","mean", 1.864, 1.8, 1.93);

  // Construct Gaussian1 PDF for signal
  RooRealVar sigma("sigma","sigma", 0.02, 0.005, 0.06);
  RooRealVar ngsig("ngsig","ngsignal", 5000, 100, 10000000);
  RooGaussian gauss("gauss","gauss", d0_mass, mean, sigma);

  // Construct exponential PDF to fit the bkg component
  RooRealVar lambda("lambda", "slope", -3., -10, 0.);
  //RooRealVar lambda("lambda", "slope", -2.9, -10, 0.);
  RooExponential expo("expo", "exponential PDF", d0_mass, lambda);
  
  // Construct signal + bkg PDF
  RooRealVar nsig("nsig","#signal events", 1300, 0, 10000) ;
  RooRealVar nbkg("nbkg","#background events", 5000, 0, 10000) ;
  RooAddPdf model("model","g+exp", RooArgList(gauss, expo), RooArgList(nsig,nbkg)) ;
  */
  RooRealVar mean("mean","mean", 1.86, 1.864-wind, 1.864+wind);

  // Construct Crystal Ball PDF for signal
  /*
  RooRealVar cbmean("cbmean", "cbmean" , 1.864, 1.85, 1.87); 
  RooRealVar cbsigma("#sigma", "cbsigma" , 0.02, 0.001, 0.25); 
  RooRealVar ncbsig("ncbsig", "ncbsignal", 4000, 0, 10000); 
  RooRealVar n("n","cbn", 5, 0, 10);
  RooRealVar alpha("#alpha","cbalpha", 1, 0, 5);
  RooCBShape cball("cball", "crystal ball", d0_mass, mean, cbsigma, alpha, n);
  */

  // Construct Gaussian1 PDF for signal
  RooRealVar sigma("sigma","sigma", 0.01, 0.005, 0.06);
  RooRealVar ngsig("ngsig","ngsignal", 200, 100, 10000000);
  RooGaussian gauss("gauss","gauss", d0_mass, mean, sigma);

  // Construct exponential PDF to fit the bkg component
  RooRealVar lambda("lambda", "slope", -0.1, -10, 10);
  RooExponential expo("expo", "exponential PDF", d0_mass, lambda);
  
  // Construct signal + bkg PDF
  RooRealVar nsig("nsig","#signal events", 4000, 0, 10000000) ;
  RooRealVar nbkg("nbkg","#background events", 4000, 0, 10000000) ;
  //RooAddPdf model("model","g+exp", RooArgList(cball, expo), RooArgList(nsig,nbkg)) ;
  RooAddPdf model("model","g+exp", RooArgList(gauss, expo), RooArgList(nsig,nbkg)) ;

  cout << "fitting model" << endl;
  meson_l_mass.setConstant();
  ptfrac.setConstant();
  weight.setConstant();
  RooPlot* frame = d0_mass.frame() ;
  model.fitTo(ds, Extended(), SumW2Error(kTRUE));
  ds.plotOn(frame,Binning(30));
  /*
  RooAbsReal *nll = model.createNLL(ds, NumCPU(8), SumW2Error(kTRUE));
  RooMinuit m(*nll);
  m.setPrintLevel(-1); 
  m.setPrintEvalErrors(-1);
  m.migrad();
  m.hesse();
  RooFitResult *r = m.save();
  */

  model.plotOn(frame);
  model.plotOn(frame, Components(expo),LineStyle(kDashed));
  frame->SetTitle("");
  frame->Draw();
  frame->SetName("massD0");
  frame->Write();

  mass = mass + "_d0mu";
  c1->SaveAs("massD0_mu_"+mass+".pdf");
  c1->SaveAs("massD0_mu_"+mass+".png");

  frame = ptfrac.frame();
  if(jpT) ds.plotOn(frame,Binning(22));
  else ds.plotOn(frame,Binning(22));
 
  //model.paramOn(frame);
  /*
  RooDataSet sigData;
  RooDataSet bkgData;
  */

  SPlot sData("sData","An SPlot from mass for sum ch pT", ds, &model, RooArgList(nsig,nbkg));
  RooDataSet sigData = RooDataSet(ds.GetName(), ds.GetTitle(), &ds, *ds.get(), "", "nsig_sw");
  RooDataSet bkgData = RooDataSet(ds.GetName(), ds.GetTitle(), &ds, *ds.get(), "", "nbkg_sw");
  std::cout << sigData.weight() << std::endl;
  std::cout << sData.GetYieldFromSWeight("nsig") << std::endl;
  std::cout << sData.GetYieldFromSWeight("nbkg") << std::endl;
  std::cout << sData.GetYieldFromSWeight("nsig")/sData.GetYieldFromSWeight("nbkg") << std::endl;
  /*
  for(int i = 0; i < 10; i++) {
    std::cout << "weight " << sData.GetSWeight(i,"nsig) << std::endl;
    std::cout << "Total Weight   " << sData.GetSumOfEventSWeight(i) << std::endl;
  }
  */
  /*
  sigData.SetName("sSignal");
  sigData.SetTitle("sSignal");
  bkgData.SetName("sBackground");
  bkgData.SetTitle("sBackground");
  */
  w.import(d0_mass);
  w.import(meson_l_mass);
  w.import(ptfrac);
  w.import(d0_pt);
  w.import(j_pt_ch);
  w.import(j_pt);
  w.import(epoch);
  w.import(weight);
  w.import(tuneW);
  //w.import(bkgData, Rename("bkgData")); //Redundant, just clone dataset and change weight to nsig_bkg
  w.import(sigData, Rename("sigData"));
  w.import(ds, Rename("dsSWeights"));
  /*
  for(auto &it : frgWgt)
    w.import(it);
  */
  /*
  for(int i=0; i < 10; i++) {
    std::cout << "Weight" << sData.GetSWeight(i,"nsig") << std::endl;
    std::cout << "Weight" << sData.GetSWeight(i,"nbkg") << std::endl;
    std::cout << "Total weight" << sData.GetSumOfEventSWeight(i) << std::endl;
  }
  */

  frame->Draw();

  int binning(22);
  if(jpT) binning = 22;

  /*
  frame = d0_mass.frame();
  sigData.plotOn(frame, DataError(RooAbsData::SumW2),
                 RooFit::Name("massD0_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(30));
  frame->SetName("massD0_signal");
  frame->Draw();
  frame->Write();
  c1->SaveAs("massD0_signal_"+mass+".pdf");
  c1->SaveAs("massD0_signal_"+mass+".png");

  frame = d0_mass.frame();
  bkgData.plotOn(frame, DataError(RooAbsData::SumW2),
                 RooFit::Name("massD0_background"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(30));
  frame->SetName("massD0_background");
  frame->Draw();
  frame->Write();
  c1->SaveAs("massD0_background_"+mass+".pdf");
  c1->SaveAs("massD0_background_"+mass+".png");

  RooPlot *frame2 = d0_mass.frame();
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2));
  frame2->Draw();
  */

  RooPlot *frame2 = ptfrac.frame();
  bkgData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("ptfrac_bkg"), RooFit::MarkerColor(419),
                 RooFit::MarkerStyle(24), RooFit::MarkerStyle(24),
                 RooFit::LineWidth(2), RooFit::LineColor(419), Binning(binning));
  frame2->Draw();
  frame2->SetName("ptfrac_bkg");
  TGraph *h_frac = (TGraph*)c1->GetPrimitive("ptfrac_bkg");
  c1->SaveAs("ptfrac_bkg_"+mass+".pdf");
  c1->SaveAs("ptfrac_bkg_"+mass+".png");
  std::cout << "Background: " << frame2->getHist()->Integral() << std::endl;
  frame2 = ptfrac.frame();
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("ptfrac_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(binning));
  frame2->Draw();
  std::cout << "Signal: " << frame2->getHist()->Integral() << std::endl;
  frame2->SetName("ptfrac_mu_tag_signal");
  frame2->SetTitle("");
  frame2->Write();

  if(jpT) mass = mass + "_jpT";
  c1->SaveAs("ptfrac_signal_"+mass+".pdf");
  c1->SaveAs("ptfrac_signal_"+mass+".png");

  frame2 = ptfrac.frame();
  bkgData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("ptfrac_mu_tag_bkg"), RooFit::MarkerColor(419),
                 RooFit::MarkerStyle(24), RooFit::MarkerStyle(24),
                 RooFit::LineWidth(2), RooFit::LineColor(419), Binning(binning));
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("ptfrac_mu_tag_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(binning));
  frame2->SetTitle("");
  frame2->Draw();

  //show a legend
  c1->cd();
  TLegend *leg = new TLegend(0.6,0.75,0.9,0.9,"","brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry("ptfrac_mu_tag_signal", "Signal","p");
  leg->AddEntry("ptfrac_mu_tag_bkg","Background","p");
  leg->Draw();

  c1->SaveAs("ptfrac_"+mass+".pdf");
  c1->SaveAs("ptfrac_"+mass+".png");

  binning = 25;

  frame2 = meson_l_mass.frame();
  bkgData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_mu_tag_bkg"), RooFit::MarkerColor(419),
                 RooFit::MarkerStyle(24), RooFit::MarkerStyle(24),
                 RooFit::LineWidth(2), RooFit::LineColor(419), Binning(binning));
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_mu_tag_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(binning));
  frame2->SetTitle("");
  frame2->Draw();
  delete leg;

  //show a legend
  c1->cd();
  leg = new TLegend(0.6,0.75,0.9,0.9,"","brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry("meson_l_mass_mu_tag_signal", "Signal","p");
  leg->AddEntry("meson_l_mass_mu_tag_bkg","Background","p");
  leg->Draw();

  c1->SaveAs("meson_l_mass_"+mass+".pdf");
  c1->SaveAs("meson_l_mass_"+mass+".png");

  frame2 = meson_l_mass.frame();
  frame2->SetTitle("");
  bkgData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_bkg"), RooFit::MarkerColor(419),
                 RooFit::MarkerStyle(24), RooFit::MarkerStyle(24),
                 RooFit::LineWidth(2), RooFit::LineColor(419), Binning(binning));
  /*
  frame2->Draw();
  c1->SaveAs("meson_l_mass_bkg_"+mass+".pdf");
  c1->SaveAs("meson_l_mass_bkg_"+mass+".png");
  std::cout << "Background: " << frame2->getHist()->Integral() << std::endl;
  frame2 = meson_l_mass.frame();
  sigData.plotOn(frame2, DataError(RooAbsData::SumW2),
                 RooFit::Name("meson_l_mass_signal"), RooFit::MarkerColor(1),
                 RooFit::MarkerStyle(20), RooFit::MarkerStyle(20),
                 RooFit::LineWidth(2), RooFit::LineColor(1), Binning(binning));
  */
  frame2->Draw();
  std::cout << "Signal: " << frame2->getHist()->Integral() << std::endl;
  frame2->SetName("meson_l_mass_mu_tag_signal");
  frame2->SetTitle("");
  frame2->Write();

  c1->SaveAs("meson_l_mass_signal_"+mass+".pdf");
  c1->SaveAs("meson_l_mass_signal_"+mass+".png");

  /*
  TH1F *signalGr=(TH1F*)c1->GetPrimitive("ptfrac_signal");
  signalGr->SaveAs("/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/mtop/ptfrac_signal.pdf");
  signalGr->SaveAs("/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/LJets2015/2016/mtop/ptfrac_signal.png");
  */

  std::cout << "writting workspace to file" << std::endl;
  w.Write();
  std::cout << "writting tuneWgt to file" << std::endl;
  tuneWgt->Write();
  std::cout << "closing file" << std::endl;
  fout->Close();
  std::cout << "DONE!" << std::endl;

  return;
}

