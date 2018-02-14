{
  std::vector<float> names = {166.5,169.5,171.5,172.5,173.5,175.5,178.5};
  std::vector<pair<float,float>> mfits = {};
  std::vector<pair<float,float>> fit_par;
  std::vector<pair<float,float>> fit_err;
  std::vector<RooRealVar> masses;

  //Fit MC and get fit parameters
            //0b vary binned
  short flags(0b11);
  gROOT->ProcessLine(".L splot.C");
  //Fit and plot fitted masses
  gROOT->ProcessLine(".L fit_constrain_d0.C");
  //gROOT->ProcessLine(".L roofit_mtop_d0_unfold.C");
  gROOT->ProcessLine(".L roofit_mtop_d0.C");

  bool isData(false);
  RooWorkspace w = create_workspace(isData);
  roofit_mtop(w,names,fit_par,fit_err,flags);
  return;
  //TH1F *mass = new TH1F("mass","mass;m_{t}^{GEN};m_{t}^{FIT}",100,165,179);//,50,163,180);
  TH1F *mass = new TH1F("mass","mass;m_{t}^{GEN}-172.5 (GeV);m_{t}^{FIT} (GeV)",100,-8,8);//,50,163,180);
  for(auto & it : names) {
    TString tmp_mass = Form("%.1f",it);
    tmp_mass.ReplaceAll(".","v");
    RooRealVar mt = fit_constrain(w,fit_par,fit_err,tmp_mass,flags);
    mt.Print();
    mass->SetBinContent(mass->FindBin(it-172.5),mt.getValV()+172.5);
    mass->SetBinError(mass->FindBin(it-172.5),mt.getError());
    mfits.push_back(pair<float,float>(mt.getValV()+172.5,mt.getError()));
  }

  //Print fit parameters
  for(int i = 0; i < fit_par.size(); i++) {
    std::cout << fit_par[i].first << " + " << fit_par[i].second << " * mt" << std::endl;
    std::cout << fit_err[i].first << "   " << fit_err[i].second << std::endl;
  }
  TFitResultPtr f = mass->Fit("pol1","FS");
  std::pair<float,float> f_err = std::pair<float,float>(f->ParError(0),f->ParError(1));
  std::cout << f->Parameter(0) << " + " << f->Parameter(1) << " * mt" << std::endl;
  std::cout << "(" << f_err.first << ") + (" << f_err.second << ")" << std::endl;
  std::cout << f->Chi2() << " " << f->Ndf() << " " << f->Chi2()/f->Ndf() << std::endl;
  mass->Draw();
  float avg_delta(0.);
  float avg_deltag(0.);
  float avg_var(0.);
  int ndelta(0);
  float mmin(172), mmax(172);
  for(auto & it : names) {
    float itg(it -172.5);
    TString tmp_mass = Form("%.1f",it);
    tmp_mass.ReplaceAll(".","v");
    float calibm = f->Parameter(0) + f->Parameter(1)*itg;
    float deltam = calibm - mfits[&it - &names[0]].first;
    float deltag = calibm - it;
    float deltam_unc = mfits[&it - &names[0]].second;
    //float deltam_unc = min(sqrt(pow(deltam,2) + pow(mfits[&it - &names[0]].second,2)), sqrt(pow(deltam,2) - pow(mfits[&it - &names[0]].second,2)));
    float down = min(abs(deltam - mfits[&it - &names[0]].second), abs(deltam + mfits[&it - &names[0]].second));
    if(abs(deltam) < abs(mfits[&it - &names[0]].second)) deltam = 0;
    else deltam = deltam/abs(deltam) * min(abs(deltam) - abs(mfits[&it - &names[0]].second), abs(deltam) + abs(mfits[&it - &names[0]].second));
    if(deltam != 0) ndelta++;
    //avg_delta += pow(deltam,2);
    avg_delta += deltam;
    //avg_deltag += pow(deltag,2);
    avg_deltag += deltag;
    avg_var += sqrt( pow( deltam_unc, 2) );
    std::cout << it << " GeV -> RooFit mass: " << mfits[&it - &names[0]].first
              << " GeV Calibration line mass: " << calibm << " GeV ("  << deltam << " GeV) +/-" << abs(deltam_unc) << " GeV " << deltag << " GeV" << std::endl;
    mmin = min(mmin, mfits[&it - &names[0]].first);
    mmax = max(mmax, mfits[&it - &names[0]].first);
  }
  //avg_delta /= ndelta;
  avg_delta /= names.size();
  //avg_delta = sqrt(avg_delta);
  avg_deltag /= names.size();
  //avg_deltag = sqrt(avg_deltag);
  avg_var /= names.size();
  for(auto & it : names) {
    float itg(it -172.5);
    TString tmp_mass = Form("%.1f",it);
    tmp_mass.ReplaceAll(".","v");
    float calibm = f->Parameter(0) + f->Parameter(1)*itg;
    std::cout << "Mass: " << it << " GeV -> " << calibm - avg_deltag << " GeV +/- " << sqrt(pow(f->ParError(0),2)+pow(f->ParError(1),2)) << " GeV" << std::endl;
  }
  std::cout << "AVG difference: " << avg_delta << " GeV" << std::endl;
  std::cout << "AVG from GEN: " << avg_deltag << " GeV" << std::endl;
  //mass->GetYaxis()->SetRangeUser(174,177.5);
  int lbin = mass->FindFirstBinAbove(0);
  int ubin = mass->FindLastBinAbove(0);
  /*
  mass->GetYaxis()->SetRangeUser(min((int)(mass->GetBinContent(lbin) - mass->GetBinError(lbin))-1, (int)(mass->GetBinContent(lbin) + mass->GetBinError(lbin))-1),
                                 max((int)(mass->GetBinContent(ubin) - mass->GetBinError(ubin))+5, (int)(mass->GetBinContent(ubin) + mass->GetBinError(ubin))+1)+5);
  */
  mass->GetYaxis()->SetRangeUser((int)mmin-1, (int)mmax+5);
  TString name("");
  if(flags&0x2) name = "_meson_vary";
  else name = "_meson";
  char sign = '+';
  if(f->Parameter(1)<0) sign = '-';
  TString leg_title = TString::Format("Calibration curve : %0.2f (#pm%0.2f) %c m_{t} %0.2f (#pm%0.2f)",f->Parameter(0),abs(f->ParError(0)),sign,f->Parameter(1),abs(f->ParError(1)));
  TLegend *leg_calib = new TLegend(0.14,0.75,0.67,0.88,NULL,"brNDC");
  leg_calib->SetBorderSize(0);
  leg_calib->AddEntry(mass,leg_title,"lp");
  leg_calib->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);  
  latex.SetTextSize(0.025);
  latex.DrawLatex(0.16,0.73,TString::Format("#chi^{2}/n.d.f. = %0.2f",f->Chi2()/f->Ndf()));

  c1->SaveAs("masses"+name+".pdf");
  c1->SaveAs("masses"+name+".png");
}
