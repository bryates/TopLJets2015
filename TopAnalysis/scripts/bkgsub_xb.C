#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <RooBinning.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooAddPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/param.h"
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/convert.h"
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/tdr.h"

TH1F *nombkgB = nullptr;
TH1F *nombkgG = nullptr;
TH1F *nombkgBW = nullptr;
TH1F *nombkgGW = nullptr;
TH1F *nombkgBc = nullptr;
TH1F *nombkgGc = nullptr;
TH1F *nombkgBo = nullptr;
TH1F *nombkgGo = nullptr;
TH1F *sig_hist = nullptr;
float nsigB = 0;
float nsigG = 0;
std::vector<TString> tune = {"_sdown", "_700", "_725", "_down", "_dddown", "_ddown", "_scentral", "", "_cccentral", "_ccentral", "_925", "_central", "_uuup", "_up"};//, "_1125" };
TString epoch_name[4] = {"_xb", "_BCDEFGH", "_BCDEF", "_GH"};

void bkgsub_xb(TString name, bool isData=false, TString nom = "172v5_Bfrag_genreco") {
  gROOT->ProcessLine("gStyle->SetOptStat(0);");
std::map<TString, float> param;
param["sdown"]=0.655;
param["700"]=0.7;
param["725"]=0.725;
param["down"]=0.755;
param["dddown"]=0.775;
param["ddown"]=0.8;
param["scentral"]=0.825;
param[""]=0.855;
param["cccentral"]=0.875;
param["ccentral"]=0.9;
param["925"]=0.925;
param["central"]=0.955;
param["uuup"]=0.975;
param["uup"]=1.000;
param["up"]=1.055;

  /*
  TH1F *h_tot = nullptr;
  TH1F *bkg_tot = nullptr;
  */
  auto c1 = setupCanvas();
  setupPad()->cd();
  RooBinning bins(0,1.1);
  double *dbin = new double[bin.size()-1];
  for(int i = 0; i < (int)bin.size(); i++) { bins.addBoundary(bin[i]); dbin[i] = bin[i];}
  float nsig(0);
  float nbkg(0);
  float ndiv(0);
  for(int ep = 1; ep <= 2; ep++) {
    TH1F *h = nullptr;
    TH1F *bkg = nullptr;
    TH1F *bkgW = nullptr;
    TH1F *bkgc = nullptr;
    TH1F *bkgo = nullptr;
    TH1F *bkgDss = new TH1F("ptfrac_bkgDss_hist", "ptfrac_bkgDss_hist", bin.size()-1, bin.data());
    TH1F *bkgDssu = new TH1F("ptfrac_bkgDssu_hist", "ptfrac_bkgDssu_hist", bin.size()-1, bin.data());
    std::vector <float> scalesig;
    std::vector <float> scalebkg;
    std::vector <float> scale = {0,0,0};
    if(name.Contains("Data")) {
      std::cout << nombkgB << std::endl;
      std::cout << nombkgG << std::endl;
      if(ep==1) bkg = nombkgB;
      else if(ep==2) bkg = nombkgG;
    }
    if(!isData) {
      std::vector<TString> cuts = {"j_hadflav==5 && mcFile==4", "j_hadflav==4 && mcFile==5", "mcFile<4"};
      std::vector<TString> cutname = {"W", "c", "o"};
      for(int i = 0; i < (int)cuts.size(); i++) {
        RooDataSet *ds = nullptr;
        std::cout << "Processing " << name + cutname[i] << std::endl;
        RooPlot *frame = nullptr;
        for(auto xb : type) {
          if(xb < 0.4) continue;
          std::cout << xb << std::endl;
          TString fname = (TString::Format("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_%s_sPlot_d0%d_%d.root", name.Data(), ep, int(xb*1000) ));
          if(name.Contains("Wup")) fname.ReplaceAll("Wup", "172v5");
          if(name.Contains("Wdown")) fname.ReplaceAll("Wdown", "172v5");
          if(name.Contains("Cup")) fname.ReplaceAll("Cup", "172v5");
          if(name.Contains("Cdown")) fname.ReplaceAll("Cdown", "172v5");
          if(name.Contains("Oup")) fname.ReplaceAll("Oup", "172v5");
          if(name.Contains("Odown")) fname.ReplaceAll("Odown", "172v5");
          std::cout << fname << std::endl;
          TFile *fin = nullptr;
          try {
            fin = new TFile(fname);
          }
          catch(...) {
            std::cout << "There is an error with " << fname << std::endl;
            break;
          }
          if(fin == nullptr) {
            std::cout << fname << " not found!" << std::endl;
            break;
          }
          /*
          */
          std::cout << "here0" << std::endl;
          //auto fin = new TFile(fname);
          std::cout << "here1" << std::endl;
          auto w = (RooWorkspace*)fin->Get("w");
          std::cout << "here2" << std::endl;
          if(ds == nullptr) {
            ds = new RooDataSet("ds", "ds", (RooDataSet*)w->data("dsSWeights"), *w->data("dsSWeights")->get(), 0, "weight");
            frame = w->var("d0_mass")->frame();
          }
          else {
            auto tmpds = new RooDataSet("ds", "ds", (RooDataSet*)w->data("dsSWeights"), *w->data("dsSWeights")->get(), 0, "weight");
            ds->append(*tmpds);
            //delete tmpds;
          }
          std::cout << "here3" << std::endl;
        }
        std::cout << cuts[i] << std::endl;
        std::cout << cutname[i] << std::endl;
        TString fname = (TString::Format("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_172v5_sPlot_d0%d.root", ep ));
        //TString fname = (TString::Format("sPlot/sPlot/TopMass_%s_sPlot_d0%d.root", name.Data(), ep ));
          if(name.Contains("Wup")) fname.ReplaceAll("Wup", "172v5");
          if(name.Contains("Wdown")) fname.ReplaceAll("Wdown", "172v5");
          if(name.Contains("Cup")) fname.ReplaceAll("Cup", "172v5");
          if(name.Contains("Cdown")) fname.ReplaceAll("Cdown", "172v5");
          if(name.Contains("Oup")) fname.ReplaceAll("Oup", "172v5");
          if(name.Contains("Odown")) fname.ReplaceAll("Odown", "172v5");
        auto fin = new TFile(fname);
        auto w = (RooWorkspace*)fin->Get("w");
        auto d = (RooDataSet*)ds->reduce(cuts[i]);
        //auto frame = w->var("d0_mass")->frame();
        auto c1 = setupCanvas();
        c1->cd();
        setupPad()->cd();
        d->plotOn(frame, RooFit::Binning(30), RooFit::DrawOption("AP"));
        auto mass = (TH1F*)convert(frame, false, 1.7, 2.0);
        //auto gauss = (RooBukinPdf*)w->pdf("bukin_novos");
        auto gauss = (RooGaussian*)w->pdf("gauss");
        auto expo = (RooAddPdf*)w->pdf("expo");
        //auto model = new RooAddPdf("model","g+a", RooArgList(*gauss, *expo), RooArgList(*w->var("nsig"),*w->var("nbkg"))) ;
        auto model = (RooAddPdf*)w->pdf("model");
        if(cutname[i] == "c")
          model = (RooAddPdf*)w->pdf("model");
        float totsig = w->var("nsig")->getVal();
        float totbkg = w->var("nbkg")->getVal();
        model->fitTo(*d);
        model->plotOn(frame);
        frame->Draw();
        float chi2 = frame->chiSquare();//"model", "data", 3);
        std::cout << std::endl << cutname[i] << " chi2 " << chi2 << std::endl << std::endl;
        if(name=="172v5") {
        c1->SaveAs(TString::Format("www/meson/morph/mass/massD0_bkgHand_%s%s%s_d0_xb.pdf",name.Data(), cutname[i].Data(), epoch_name[ep+1].Data()));
        c1->SaveAs(TString::Format("www/meson/morph/mass/massD0_bkgHand_%s%s%s_d0_xb.png",name.Data(), cutname[i].Data(), epoch_name[ep+1].Data()));
        }
        //delete d;
        /*
        */
        auto framem = w->var("d0_mass")->frame();
        //w->pdf("model")->plotOn(framem, RooFit::Components("gauss"), RooFit::Binning(60), RooFit::DrawOption("AL"));
        //scalesig[i] = framem->getCurve()->Integral() * w->var("nsig")->getVal();
        /*
        for(int ibin = 0; ibin < 60; ibin++) {
          float imass = mass->GetBinCenter(ibin);
          if(imass < (w->var("mean")->getVal() - 2 * w->var("sigma")->getVal()) || imass > (w->var("mean")->getVal() + 2 * w->var("sigma")->getVal())) continue;
          scalesig[i] += framem->getCurve()->Integral() * mass->GetBinContent(ibin);
        }
        */
        //w->var("d0_mass")->setRange("window", w->var("mean")->getVal() - 2 * w->var("sigp")->getVal(), w->var("mean")->getVal() + 2 * w->var("sigp")->getVal());
        w->var("d0_mass")->setRange("window", w->var("mean")->getVal() - 2 * w->var("sigma")->getVal(), w->var("mean")->getVal() + 2 * w->var("sigma")->getVal());
        //scalesig[i] = w->pdf("gauss")->createIntegral(*w->var("d0_mass"),"window")->getVal() * w->var("nsig")->getVal();
        auto x = *w->var("d0_mass");
        //x.setRange("window",w->var("mean")->getVal() - 2 * w->var("sigp")->getVal(), w->var("mean")->getVal() + 2 * w->var("sigp")->getVal());
        x.setRange("window",w->var("mean")->getVal() - 2 * w->var("sigma")->getVal(), w->var("mean")->getVal() + 2 * w->var("sigma")->getVal());
        //x.setRange("window",1.8,1.93);
        RooAbsReal *intModel = model->createIntegral(x);
        //auto gauss = *(RooGaussian*)w->pdf("gauss");
        RooAbsReal *intsig = gauss->createIntegral(x,x,"window");
        //auto expo = *(RooAddPdf*)w->pdf("bkgModel");
        RooAbsReal *intbkg = expo->createIntegral(x,x,"window");
        auto nsig = *w->var("nsig");
        auto nbkg = *w->var("nbkg");
        float sig = nsig.getVal() * intsig->getVal();
        float sbkg = nbkg.getVal() * intbkg->getVal();
        cout << "Signal= " << sig << endl;
        cout << "Background= " << bkg << endl;
        scalesig.push_back(sig);
        scalebkg.push_back(sbkg);
        cout << "Purity= " << sig / (sig+sbkg) << endl;
        cout << "S/B= " << sig / sbkg << endl;
        cout << "Total sig= " << nsig.getVal() << endl;
        cout << "Total bkg= " << nbkg.getVal() << endl;
        cout << "Purity tot= " << sig / (totsig+totbkg) << endl;
        cout << "Tot sig= " << totsig << endl;
        cout << "Tot bkg= " << totbkg << endl;
        //scalebkg[i] = w->pdf("bkgModel")->createIntegral(*w->var("d0_mass"),"window")->getVal() * w->var("nbkg")->getVal();
        //std::cout << "sig = " << framem->getCurve()->Integral() * w->var("nsig")->getVal() << std::endl;;
        std::cout << "sig = " << scalesig[i] << std::endl;;
        w->pdf("model")->plotOn(framem, RooFit::Components("bkgModel"), RooFit::Binning(60), RooFit::DrawOption("AL"));
        //scalebkg[i] = framem->getCurve()->Integral() * w->var("nbkg")->getVal();
        /*
        for(int ibin = 0; ibin < 60; ibin++) {
          float imass = mass->GetBinCenter(ibin);
          if(imass < (w->var("mean")->getVal() - 2 * w->var("sigma")->getVal()) || imass > (w->var("mean")->getVal() + 2 * w->var("sigma")->getVal())) continue;
          scalebkg[i] += framem->getCurve()->Integral() * mass->GetBinContent(ibin);
        }
        */
        //std::cout << "bkg = " << framem->getCurve()->Integral() * w->var("nbkg")->getVal() << std::endl;;
        std::cout << "bkg = " << scalebkg[i] << std::endl;;
        if((i== 1 && name.EqualTo(nom)) || i!=1) scale.push_back(scalesig[i] / (scalesig[i] + scalebkg[i]));
        //if((i== 1 && name.EqualTo("isr-up")) || i!=1) scale[i] = scalesig[i] / (scalesig[i] + scalebkg[i]);
        std::cout << "sig/bkg = " << scale[i] << std::endl;;

      }
    }
    for(auto b : type) {
      if((&b - &type[0]) == 0) continue;
      if(b < 0.4) continue;
      std::cout << "Runing over xb < " << b << std::endl;
      TString fname = (TString::Format("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_%s_sPlot_d0%d_%d.root", name.Data(), ep, int(b * 1000)));
          if(name.Contains("Wup")) fname.ReplaceAll("Wup", "172v5");
          if(name.Contains("Wdown")) fname.ReplaceAll("Wdown", "172v5");
          if(name.Contains("Cup")) fname.ReplaceAll("Cup", "172v5");
          if(name.Contains("Cdown")) fname.ReplaceAll("Cdown", "172v5");
          if(name.Contains("Oup")) fname.ReplaceAll("Oup", "172v5");
          if(name.Contains("Odown")) fname.ReplaceAll("Odown", "172v5");
      std::cout << "Opening " << fname << std::endl;
      TFile *fin = nullptr;
      fin = new TFile(fname);
      if(fin == nullptr) {
        std::cout << fname << " not found!" << std::endl;
        break;
      }
      std::cout << "Opened!" << std::endl;
      auto w = (RooWorkspace*)fin->Get("w");
      std::cout << w->var("nsig")->getVal() << "\t" << w->var("nbkg")->getVal() << "\t" << w->var("nsig")->getVal() / w->var("nbkg")->getVal() << std::endl;
      nsig += w->var("nsig")->getVal();
      nbkg += w->var("nbkg")->getVal();
      //ndiv += w->var("nsig")->getVal() / w->var("nbkg")->getVal() << std::endl;
      std::cout << "scaling" << std::endl;
      /*
      std::cout << "epoch " << ep << " xb < " << b << "\t"
                << TString::Format("%0.2f", w->var("nsig")->getVal()) << "+/-" << TString::Format("%0.2f", w->var("nsig")->getError()) << "\t"
                << TString::Format("%0.2f", w->var("mean")->getVal()) << "+/-" << TString::Format("%0.2f", w->var("mean")->getError()) << "\t"
                << TString::Format("%0.3f", w->var("sigma")->getVal()) << "+/-" << TString::Format("%0.3f", w->var("sigma")->getError()) << "\t"
                << TString::Format("%0.2f", w->var("lambda")->getVal()) << "+/-" << TString::Format("%0.2f", w->var("lambda")->getError()) << "\t"
                << TString::Format("%0.2f", w->var("ngsigkk")->getVal()) << "+/-" << TString::Format("%0.2f", w->var("ngsigkk")->getError()) << "\t"
                << TString::Format("%0.2f", w->var("meankk")->getVal()) << "+/-" << TString::Format("%0.2f", w->var("meankk")->getError()) << "\t"
                << TString::Format("%0.3f", w->var("sigmakk")->getVal()) << "+/-" << TString::Format("%0.3f", w->var("sigmakk")->getError()) << std::endl;
      */
      auto frame = w->var("ptfrac")->frame();
      TH1F *tmp = nullptr;
      if(!isData) {
        std::cout << "getting bkg" << std::endl;
        auto ds = new RooDataSet("ds", "ds", (RooDataSet*)w->data("sigData"), *w->data("sigData")->get(), 0, "nsig");
        //auto dsn = new RooDataSet("ds", "ds", (RooDataSet*)w->data("dsSWeights"), *w->data("dsSWeights")->get(), 0, "weight");
        //float ratio = 0.319658; //D**
        float ratio_nom = 0.005; //D*
        //float ratio = 0.192; //D*
        float ratio = 0.5; //D*
        std::cout << "Opening fit" << std::endl;
        auto fin = TFile::Open("/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/bfrag/xb_fit.root");
        /*
        auto fin = TFile::Open("/afs/cern.ch/user/b/byates/TopAnalysis/data/era2016/bfragweights.root");
        std::cout << "Loading the sum ch pT morph" << std::endl;
        auto sumptmorph = (TGraph*)fin->Get("Dzrand413Frag")->Clone();
        std::cout << "Loading the sum ch pT morph DONE!" << std::endl;
	RooRealVar tuneW = RooRealVar("tuneW", "tuneW", 1., 0, 2.);
        std::cout << "Morphing ds" << std::endl;
        for(int i = 0; i < dsn->numEntries(); i++) {
          double tmpW = dsn->get(i)->getRealValue("nsig_sw");
          tuneW.setVal(tmpW);// * sumptmorph->Eval(dsn->get(i)->getRealValue("ptfrac")));
          //tuneW.setVal(tmpW * sumptmorph->Eval(dsn->get(i)->getRealValue("ptfrac")));
          //std::cout << tuneW.getVal() << " " << weight.getVal() << std::endl;;
        }
        std::cout << "Morphing ds DONE!" << std::endl;
        auto ds = new RooDataSet("ds", "ds", dsn, *dsn->get(), 0, "tuneW");
        ds->plotOn(frame, RooFit::Binning(bins));
        sig_hist = (TH1F*)convert(frame, false, bin);
        sig_hist->SetDirectory(0);
        */
        /* FIXME use for Dss
        //bkgDss = (TH1F*)fin->Get("bfragAnalysis/xb_DssDz")->Clone();
        bkgDss = (TH1F*)fin->Get("bfragAnalysis/xb_testcharged413")->Clone();
        //bkgDssu = (TH1F*)fin->Get("bfragAnalysis/xb_DssDzu")->Clone();
        //bkgDssu = (TH1F*)fin->Get("bfragAnalysis/xb_DssDzd")->Clone();
        bkgDssu = (TH1F*)fin->Get("bfragAnalysis/xb_testucharged413")->Clone();
        //bkgDss->SetDirectory(0);
        //bkgDssu->SetDirectory(0);
        bkgDss  = (TH1F*)bkgDss->Rebin(bin.size()-1, "", dbin);
        bkgDssu = (TH1F*)bkgDssu->Rebin(bin.size()-1, "", dbin);
        bkgDss->Scale(ratio_nom * w->var("nsig")->getVal() / bkgDss->Integral());
        //bkgDss->Scale(ratio * w->var("nsig")->getVal() / bkgDss->Integral());
        bkgDssu->Scale(ratio * w->var("nsig")->getVal() / bkgDssu->Integral());
        std::cout << "Dss = " << bkgDss->Integral() << "\t" << bkgDssu->Integral() << std::endl;
        //bkgDssu->Scale(1. - 0.20); // 20% up/down for all D**
        //bkgDssu->Scale(1. - 0.035); // +/- 3.5% D* propagated unc.
        */
        std::cout << "nsig=" << w->var("nsig")->getVal() << "\tDss bkg = " << bkgDss->Integral() << std::endl;
        //delete ds;
        ds = new RooDataSet("ds", "ds", (RooDataSet*)w->data("dsSWeights"), *w->data("dsSWeights")->get(), 0, "weight");
        ds->reduce("(j_hadflav==5 && mcFile<=4) || (j_hadflav==4 && mcFile==5)")->plotOn(frame, RooFit::Binning(bins));
        //w->data("sigData")->reduce("(j_hadflav==5 && mcFile<=4) || (j_hadflav==4 && mcFile==5)")->plotOn(frame, RooFit::Binning(bins));
        //w->data("sigData")->reduce("mcFile<=4 || (j_hadflav==4 && mcFile==5)")->plotOn(frame, RooFit::Binning(12));
        std::cout << "getting bkg DONE!" << std::endl;
        std::cout << "converting" << std::endl;
        tmp = (TH1F*)convert(frame, false, bin);
        tmp->SetDirectory(0);
        //tmp->Scale(w->var("nsig")->getVal() / w->var("nbkg")->getVal());
        std::cout << "saving/adding bkg" << std::endl;
        std::cout << bkg << std::endl;
        if(bkg == nullptr) bkg = tmp;
        else { std::cout << "adding" << std::endl; bkg->Add(tmp); }
        std::cout << bkg << std::endl;
        std::cout << bkg->Integral() << std::endl;
        /*
        if(bkg_tot == nullptr) bkg_tot = tmp;
        else bkg_tot->Add(tmp);
        */
        std::cout << "new frame" << std::endl;
        //w->var("d0_mass")->setRange(w->var("mean")->getVal() - 2 * w->var("sigp")->getVal(), w->var("mean")->getVal() + 2 * w->var("sigp")->getVal());
        //w->var("d0_mass")->setRange(w->var("mean")->getVal() - 2 * w->var("sigma")->getVal(), w->var("mean")->getVal() + 2 * w->var("sigma")->getVal());
        auto cut = TString::Format("d0_mass>%0.3f && d0_mass<%0.3f", w->var("mean")->getVal() - 2 * w->var("sigma")->getVal(), w->var("mean")->getVal() + 2 * w->var("sigma")->getVal());
        std::cout << "cut=" << cut << std::endl;
        auto framem = w->var("d0_mass")->frame();
        //ds->reduce("j_hadflav==5 && mcFile==4 && " + cut)->plotOn(framem, RooFit::Binning(30), RooFit::DrawOption("AP"));
        auto dstmp = ds->reduce("j_hadflav==5 && mcFile==4");
        w->var("d0_mass")->setRange(1.7, 2.0);
        dstmp->plotOn(framem, RooFit::Binning(30), RooFit::DrawOption("AP"));
        w->pdf("model")->fitTo(*ds);
        w->pdf("model")->plotOn(framem);
        //w->pdf("model")->plotOn(framem, RooFit::Components("gauss"), RooFit::Binning(30), RooFit::DrawOption("AL"));
        auto bkgModel = (RooAddPdf*)w->pdf("bkgModel");
        auto gauss = (RooGaussian*)w->pdf("gauss");
        //w->pdf("model")->plotOn(framem, RooFit::Components("bkgModel"), RooFit::LineStyle(kDashed));
        //w->pdf("model")->plotOn(framem, RooFit::Components("gauss"), RooFit::LineStyle(kDashed));
        //float scalesig = framem->getCurve()->Integral();
        //w->pdf("model")->plotOn(framem, RooFit::Components("bkgModel"), RooFit::Binning(30), RooFit::DrawOption("AL"));
        framem->Draw();
        c1->SaveAs(TString::Format("www/meson/morph/mass/massD0_bkgHand_%sW%s_d0_%d.pdf",name.Data(), epoch_name[ep+1].Data(), int(b*1000)));
        c1->SaveAs(TString::Format("www/meson/morph/mass/massD0_bkgHand_%sW%s_d0_%d.png",name.Data(), epoch_name[ep+1].Data(), int(b*1000)));
        //float scalebkg = framem->getCurve()->Integral();
        //float scale = scalesig / scalebkg;
        //std::cout << "PDF integral" << scalesig << "\t" << scalebkg << "\t" << scale << std::endl;
        w->var("d0_mass")->setRange(w->var("mean")->getVal() - 2 * w->var("sigma")->getVal(), w->var("mean")->getVal() + 2 * w->var("sigma")->getVal());
        auto framebkg = w->var("ptfrac")->frame();
        std::cout << "plotting W b-jets" << std::endl;
        //ds->reduce("j_hadflav==5 && mcFile<=4")->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        ds->reduce("j_hadflav==5 && mcFile==4 && " + cut)->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        //ds->reduce("j_hadflav==5 && mcFile==4 && d0_mass>1.83 && d0_mass<1.91")->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        //w->data("sigData")->reduce("j_hadflav==5 && mcFile==4")->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        framebkg->Draw();
        auto tmpW = (TH1F*)convert(framebkg, false, bin);
        tmpW->SetDirectory(0);
        tmpW->Scale(w->var("nsig")->getVal() / w->var("nbkg")->getVal());
        //tmpW->Scale(scale[0]);
        if(name.Contains("Wup")) tmpW->Scale(1.25);
        if(name.Contains("Wdown")) tmpW->Scale(0.75);
        if(bkgW == nullptr) bkgW = tmpW;
        else bkgW->Add(tmpW);
        //ds->reduce("j_hadflav==4 && mcFile==5 && " + cut)->plotOn(framem, RooFit::Binning(30), RooFit::DrawOption("AP"));
        w->var("d0_mass")->setRange(1.7, 2.0);
        dstmp = ds->reduce("j_hadflav==4 && mcFile==5");
        dstmp->plotOn(framem, RooFit::Binning(30), RooFit::DrawOption("AP"));
        w->pdf("model")->fitTo(*ds);
        w->pdf("model")->plotOn(framem, RooFit::Components(*bkgModel), RooFit::LineStyle(kDashed));
        w->pdf("model")->plotOn(framem, RooFit::Components(*gauss), RooFit::LineStyle(kDashed));
        //w->pdf("model")->plotOn(framem, RooFit::Components("gauss"), RooFit::Binning(30), RooFit::DrawOption("AL"));
        framem->Draw();
        c1->SaveAs(TString::Format("www/meson/morph/mass/massD0_bkgHand_%sC%s_d0_%d.pdf",name.Data(), epoch_name[ep+1].Data(), int(b*1000)));
        c1->SaveAs(TString::Format("www/meson/morph/mass/massD0_bkgHand_%sC%s_d0_%d.png",name.Data(), epoch_name[ep+1].Data(), int(b*1000)));
        framebkg = w->var("ptfrac")->frame();
        std::cout << "plotting ttbar c-jets" << std::endl;
        //ds->reduce("j_hadflav==4 && mcFile==5")->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        w->var("d0_mass")->setRange(w->var("mean")->getVal() - 2 * w->var("sigma")->getVal(), w->var("mean")->getVal() + 2 * w->var("sigma")->getVal());
        ds->reduce("j_hadflav==4 && mcFile==5 && " + cut)->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        //ds->reduce("j_hadflav==4 && mcFile==5 && d0_mass>1.83 && d0_mass<1.91")->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        //w->data("sigData")->reduce("j_hadflav==4 && mcFile==5")->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        framebkg->Draw();
        auto tmpc = (TH1F*)convert(framebkg, false, bin);
        tmpc->SetDirectory(0);
        //tmpc->Scale(w->var("nsig")->getVal() / w->var("nbkg")->getVal());
        std::cout << "Scaling tmpc" << std::endl;
        std::cout << tmpc->Integral() << std::endl;
        std::cout << tmpc->Integral() << std::endl;
        if(bkgc == nullptr) bkgc = tmpc;
        else bkgc->Add(tmpc);
        //ds->reduce("mcFile<4 && " + cut)->plotOn(framem, RooFit::Binning(30), RooFit::DrawOption("AP"));
        w->var("d0_mass")->setRange(1.7, 2.0);
        dstmp = ds->reduce("mcFile<4");
        dstmp->plotOn(framem, RooFit::Binning(30), RooFit::DrawOption("AP"));
        w->pdf("model")->fitTo(*ds);
        w->pdf("model")->plotOn(framem, RooFit::Components(*bkgModel), RooFit::LineStyle(kDashed));
        w->pdf("model")->plotOn(framem, RooFit::Components(*gauss), RooFit::LineStyle(kDashed));
        //w->pdf("model")->plotOn(framem, RooFit::Components("gauss"), RooFit::Binning(30), RooFit::DrawOption("AL"));
        //w->pdf("model")->plotOn(framem, RooFit::Components("bkgModel"), RooFit::Binning(30), RooFit::DrawOption("AL"));
        framem->Draw();
        c1->SaveAs(TString::Format("www/meson/morph/mass/massD0_bkgHand_%sO%s_d0_%d.pdf",name.Data(), epoch_name[ep+1].Data(), int(b*1000)));
        c1->SaveAs(TString::Format("www/meson/morph/mass/massD0_bkgHand_%sO%s_d0_%d.png",name.Data(), epoch_name[ep+1].Data(), int(b*1000)));
        framebkg->Draw();
        //ds->reduce("mcFile<4 && j_hadflav<5")->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        //ds->reduce("mcFile<4 && j_hadflav<5 && " + cut)->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        w->var("d0_mass")->setRange(w->var("mean")->getVal() - 2 * w->var("sigma")->getVal(), w->var("mean")->getVal() + 2 * w->var("sigma")->getVal());
        ds->reduce("mcFile<4 && " + cut)->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        //w->data("sigData")->reduce("mcFile<4")->plotOn(framebkg, RooFit::Binning(bins), RooFit::DrawOption("AP"));
        framebkg->Draw();
        auto tmpother = (TH1F*)convert(framebkg, false, bin);
        tmpother->SetDirectory(0);
        tmpother->Scale(w->var("nsig")->getVal() / w->var("nbkg")->getVal());
        //tmpother->Scale(scale[2]);
        std::cout << "other=" << tmpother->Integral() << std::endl;
        if(name.Contains("Oup")) tmpother->Scale(1.25);
        if(name.Contains("Odown")) tmpother->Scale(0.75);
        std::cout << "other (scaled)=" << tmpother->Integral() << std::endl;
        if(bkgo == nullptr) bkgo = tmpother;
        else bkgo->Add(tmpother);
        //delete ds;
      }
      //w->data("sigData")->plotOn(frame, RooFit::Binning(12));
      w->data("sigData")->plotOn(frame, RooFit::Binning(bins));
      tmp = (TH1F*)convert(frame, false, bin);
      tmp->SetDirectory(0);
      std::cout << "saving/adding bkg" << std::endl;
      //if(name.Contains("172v5_Bfrag_genreco")) tmp->Scale(((TH1F*)fin->Get("tuneFSR"))->GetBinContent(1) / ((TH1F*)fin->Get("tuneFSR"))->GetBinContent(2));
      //else if(name.Contains("172v5_Bfrag_genreco")) tmp->Scale(((TH1F*)fin->Get("tuneFSR"))->GetBinContent(1) / ((TH1F*)fin->Get("tuneFSR"))->GetBinContent(3));
      if(h == nullptr) h = tmp;
      else h->Add(tmp);
      std::cout << "saving/adding bkg DONE!" << std::endl;
      /*
      if(h_tot == nullptr) h_tot = tmp;
      else h_tot->Add(tmp);
      */
      fin->Close();
      //delete fin;
      //delete w;
    }
    if(!name.EqualTo(nom)) {
    //if(!name.EqualTo("isr-up")) 
      TString fnamenom = (TString::Format("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_%s_sPlot_d0%d_xb.root", nom.Data(), ep));
          if(name.Contains("Wup")) fnamenom.ReplaceAll("Wup", "172v5");
          if(name.Contains("Wdown")) fnamenom.ReplaceAll("Wdown", "172v5");
          if(name.Contains("Cup")) fnamenom.ReplaceAll("Cup", "172v5");
          if(name.Contains("Cdown")) fnamenom.ReplaceAll("Cdown", "172v5");
          if(name.Contains("Oup")) fnamenom.ReplaceAll("Oup", "172v5");
          if(name.Contains("Odown")) fnamenom.ReplaceAll("Odown", "172v5");
      std::cout << "Loading nominal c-jet bkg" << std::endl;
      std::cout << fnamenom << std::endl;
      TFile* finnom = new TFile(fnamenom);
      std::cout << "Getting hist" << std::endl;
      bkgc = (TH1F*)((TH1F*)(finnom->Get("ptfrac_bkgc_hist"))->Clone());
      bkgc->SetDirectory(0);
      std::cout << bkgc << std::endl;
      std::cout << "Getting hist done!" << std::endl;
      std::cout << bkgc->Integral() << std::endl;
      if(name.Contains("Cup")) bkgc->Scale(1.25);
      if(name.Contains("Cdown")) bkgc->Scale(0.75);
      std::cout << scale[1] << std::endl;
      std::cout << "Getting scale" << std::endl;
      TH1F* scale_hist = (TH1F*)((TH1F*)finnom->Get("scale"))->Clone();
      std::cout << scale_hist->Integral() << std::endl;
      scale_hist->SetDirectory(0);
      std::cout << scale_hist << std::endl;
      scale[1] = scale_hist->GetBinContent(2);
      //scale[1] = ((TH1F*)finnom->Get("scale"))->GetBinContent(2);
      std::cout << "Getting scale done!" << std::endl;
      std::cout << "sig/bkg = " << scale[1] << std::endl;;
      std::cout << "Loading nominal other jet bkg" << std::endl;
      bkgo = (TH1F*)((TH1F*)finnom->Get("ptfrac_bkgo_hist"))->Clone();
      bkgo->SetDirectory(0);
      std::cout << bkgo << std::endl;
      if(name.Contains("Oup")) bkgo->Scale(1.25);
      if(name.Contains("Odown")) bkgo->Scale(0.75);
      scale[2] = scale_hist->GetBinContent(3);
      //scale[2] = ((TH1F*)finnom->Get("scale"))->GetBinContent(3);
      std::cout << "Loading nominal W jet bkg" << std::endl;
      bkgW = (TH1F*)((TH1F*)finnom->Get("ptfrac_bkgW_hist"))->Clone();
      bkgW->SetDirectory(0);
      std::cout << bkgW << std::endl;
      if(name.Contains("Wup")) bkgW->Scale(1.25);
      if(name.Contains("Wdown")) bkgW->Scale(0.75);
      scale[0] = scale_hist->GetBinContent(1);
      //scale[0] = ((TH1F*)finnom->Get("scale"))->GetBinContent(1);
      finnom->Close();
      std::cout << "sig/bkg = " << scale[1] << std::endl;;
    }
    /*
    else if(name.EqualTo("172v5_Bfrag_genreco")) {
    //else if(name.EqualTo("isr-up")) 
      bkgc->Scale(scale[1]);
    }
    */
    if(!isData && name == nom) {
      std::cout << "saving bkg" << std::endl;
      if(ep==1) { nombkgB = (TH1F*)bkg->Clone(); nombkgB->SetDirectory(0); nombkgB->Scale(1./h->Integral()); }
      else if(ep==2) { nombkgG = (TH1F*)bkg->Clone(); nombkgG->SetDirectory(0); nombkgG->Scale(1./h->Integral()); }
      if(ep==1) { 
        nombkgBW = (TH1F*)bkgW->Clone(); 
        nombkgBW->SetDirectory(0);
        nombkgBW->Scale(1./h->Integral()); 
        nombkgBc = (TH1F*)bkgc->Clone(); 
        nombkgBc->SetDirectory(0);
        //nombkgBc->Scale(1./h->Integral()); 
        nombkgBo = (TH1F*)bkgo->Clone(); 
        nombkgBo->SetDirectory(0);
        nombkgBo->Scale(1./h->Integral()); 
      }
      else if(ep==2) { 
        nombkgGW = (TH1F*)bkgW->Clone(); 
        nombkgGW->SetDirectory(0);
        nombkgGW->Scale(1./h->Integral()); 
        nombkgGc = (TH1F*)bkgc->Clone(); 
        nombkgGc->SetDirectory(0);
        //nombkgGc->Scale(1./h->Integral()); 
        nombkgGo = (TH1F*)bkgo->Clone(); 
        nombkgGo->SetDirectory(0);
        nombkgGo->Scale(1./h->Integral()); 
      }
    }
    if(isData) {
      if(ep==1) {
        bkg = (TH1F*)nombkgB->Clone();
        bkgW = (TH1F*)nombkgBW->Clone();
        bkgc = (TH1F*)nombkgBc->Clone();
        bkgo = (TH1F*)nombkgBo->Clone();
      }
      else if(ep==2) {
        bkg = (TH1F*)nombkgG->Clone();
        bkgW = (TH1F*)nombkgGW->Clone();
        bkgc = (TH1F*)nombkgGc->Clone();
        bkgo = (TH1F*)nombkgGo->Clone();
      }
      bkg->Scale(h->Integral());
      bkgW->Scale(h->Integral());
      //bkgc->Scale(h->Integral());
      bkgo->Scale(h->Integral());
    }
    bkg = (TH1F*)bkgW->Clone();
    bkg->Add(bkgc);
    bkg->Add(bkgo);
    bkg->SetName("ptfrac_bkg_hist");
    bkg->SetTitle("ptfrac_bkg_hist");
    bkgW->SetName("ptfrac_bkgW_hist");
    bkgW->SetTitle("ptfrac_bkgW_hist");
    bkgc->SetName("ptfrac_bkgc_hist");
    bkgc->SetTitle("ptfrac_bkgc_hist");
    bkgo->SetName("ptfrac_bkgo_hist");
    bkgo->SetTitle("ptfrac_bkgo_hist");
    h->SetName("ptfrac_tot_hist");
    h->SetTitle("ptfrac_tot_hist");
    tdr(h, ep);
    h->Draw();
    tdr(h, ep);
    bkg->Draw("same");
    //bkgc->Draw("same");
    auto sig = (TH1F*)h->Clone("ptfrac_signal_hist");
    //sig = (TH1F*)sig_hist->Clone("ptfrac_signal_hist");
    sig->SetTitle("ptfrac_signal_hist");
    //sig->Add(bkgDss, -1);
    float xb_sig_low = sig->Integral(1,sig->FindBin(0.55)); // D* peaks in lowest bin
    bkgDssu->Scale(0.5 * xb_sig_low / bkgDssu->Integral()); // Scale by 50% of signal
    //bkgDssu->Scale(1. - 0.035); // +/- 3.5% D* propagated unc.
    //sig->Add(bkgDssu);
    std::cout << "adding" << std::endl;
    if(!name.Contains("FSR")) sig->Add(bkg,-1);
    for(int ibin = 1; ibin <= sig->GetNbinsX(); ibin++) {
      std::cout << ibin << ":" << sig->GetBinContent(ibin) << " " << (sig->GetBinContent(ibin) <= 0) << std::endl;
      if(sig->GetBinContent(ibin) <= 0) {
        sig->SetBinContent(ibin,0.);
        std::cout << sig->GetBinContent(ibin) << std::endl;
      }
    }
    sig->SetLineColor(kGreen);
    sig->SetMarkerColor(kGreen);
    bkg->SetLineColor(kRed);
    bkg->SetMarkerColor(kRed);
    sig->Draw("same");
    TString tmpname = TString::Format("%s_%d", TString(name(0, name.First("_"))).Data(), int(1000 * param[TString(name(name.First("_")+1,name.Length())).Data()]));
    if(name.EqualTo("172v5")) tmpname = "172v5_855";
    if(name.EqualTo("172v5_Bfrag_genreco")) tmpname = "172v5_Bfrag_genreco_855";
    if(name.EqualTo("Data")) tmpname = "Data";
    if(name.EqualTo("Data_Bfrag_genreco")) tmpname = "Data_Bfrag_genreco";
    std::cout << tmpname << std::endl;
    if(name.EqualTo("172v5_Bfrag_genreco") || name.EqualTo("Data_Bfrag_genreco")) {
    c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_%s%s_d0_xb.pdf",tmpname.Data(), epoch_name[ep+1].Data()));
    c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_%s%s_d0_xb.png",tmpname.Data(), epoch_name[ep+1].Data()));
    }
    tdr(bkgW,ep);
    bkgW->Draw();
    tdr(bkgW,ep);
        if(name=="172v5") { c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_bkgHand_%sW%s_d0_xb.pdf",name.Data(), epoch_name[ep+1].Data()));
    c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_bkgHand_%sW%s_d0_xb.png",name.Data(), epoch_name[ep+1].Data()));
        }
    tdr(bkgc,ep);
    bkgc->Draw();
    tdr(bkgc,ep);
        if(name=="172v5") {
    c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_bkgHand_%sc%s_d0_xb.pdf",name.Data(), epoch_name[ep+1].Data()));
    c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_bkgHand_%sc%s_d0_xb.png",name.Data(), epoch_name[ep+1].Data()));
    }
    tdr(bkgo,ep);
    bkgo->Draw();
    tdr(bkgo,ep);
        if(name=="172v5") {
    c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_bkgHand_%sother%s_d0_xb.pdf",name.Data(), epoch_name[ep+1].Data()));
    c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_bkgHand_%sother%s_d0_xb.png",name.Data(), epoch_name[ep+1].Data()));
    }
    /*
    if(ep==2) {
      tdr(h_tot);
      h_tot->Draw();
      auto sig_tot = (TH1F*)h_tot->Clone("ptfrac_signal_hist");
      sig_tot->SetTitle("ptfrac_signal_hist");
      sig_tot->Add(bkg_tot,-1);
      sig_tot->SetLineColor(kGreen);
      bkg_tot->SetLineColor(kRed);
      tdr(h_tot);
      c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_%s_BCDEFGH_d0_xb.pdf",name.Data()));
      c1->SaveAs(TString::Format("www/meson/morph/ptfrac/ptfrac_%s_BCDEFGH_d0_xb.png",name.Data()));
    }
    */
    //if(ep==2) begin
    //continue;
    //TString fname = (TString::Format("sPlot/sPlot/TopMass_%s_Dsch50lowint_sPlot_d0%d_xb.root", name.Data(), ep));
    //TString fname = (TString::Format("sPlot/sPlot/TopMass_%s_sumchptrand_sPlot_d0%d_xb.root", name.Data(), ep));
    TString fname = (TString::Format("/eos/cms/store/group/phys_top/byates/sPlot/TopMass_%s_sPlot_d0%d_xb.root", name.Data(), ep));
    std::cout << "Writing " << fname << std::endl;
    auto fout = new TFile(fname, "RECREATE");
    sig->Write();
    bkg->Write();
    bkgW->Write();
    bkgc->Write();
    bkgo->Write();
    bkgDss->Write();
    bkgDssu->Write();
    h->Write();
    TH1F *s = new TH1F("scale", "scale", 3, 0, 3);
    s->Fill(0., scale[0]);
    s->Fill(1., scale[1]);
    s->Fill(2., scale[2]);
    /*
    if(name.EqualTo("172v5")) {
    }
    */
    s->Write();
    fout->Close();
    //end
    /*
    delete fout;
    */
  }
  std::cout << "nsig=" << nsig << std::endl;
  std::cout << "nbkg=" << nbkg << std::endl;
  ndiv = nsig / nbkg;
  std::cout << "ndiv=" << ndiv << std::endl;
}

void bkgsub_xb() {
  gROOT->ProcessLine("gStyle->SetOptStat(0);");

  /*
  for(auto t : tune)
    bkgsub(TString::Format("172v5%s", t.Data()), false);
  */
  bkgsub_xb("172v5_Bfrag_genreco", false, "172v5_Bfrag_genreco");
  bkgsub_xb("Data_Bfrag_genreco", true, "172v5_Bfrag_genreco");
}
