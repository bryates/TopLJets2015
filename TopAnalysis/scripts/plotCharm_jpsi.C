#include <TCanvas.h>
#include <TChain.h>
#include <TH1.h>
#include <THStack.h>

void plotCharm_jpsi(int ctausig=10) {
auto c1 = new TCanvas("c1", "c1");
auto data = new TChain("data");
auto datatt = new TChain("data");
auto dataw = new TChain("data");
auto datao = new TChain("data");
auto datad = new TChain("data");
std::vector<TString> mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_TTJets_powheg","MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
std::vector<TString> mcOther = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
for(auto &it : mcSamples)
  data->Add("LJets2015/2016/test/Chunks/"+it+"_*.root");
for(auto &it : mcOther)
  datao->Add("LJets2015/2016/test/Chunks/"+it+"_*.root");
datatt->Add("LJets2015/2016/test/Chunks/MC13TeV_TTJets_powheg_*");
dataw->Add("LJets2015/2016/test/Chunks/MC13TeV_W1Jets_*");
dataw->Add("LJets2015/2016/test/Chunks/MC13TeV_W2Jets_*");
dataw->Add("LJets2015/2016/test/Chunks/MC13TeV_W3Jets_*");
dataw->Add("LJets2015/2016/test/Chunks/MC13TeV_W4Jets_*");

TString val = "jpsi_mass";
int nbins = 30;
float low = 2.8;
float high = 3.4;
std::vector<TString> dataSamples = {"Data13TeV_DoubleEG_2016B","Data13TeV_DoubleEG_2016C","Data13TeV_DoubleEG_2016D","Data13TeV_DoubleEG_2016E","Data13TeV_DoubleEG_2016F","Data13TeV_DoubleMuon_2016B","Data13TeV_DoubleMuon_2016C","Data13TeV_DoubleMuon_2016D","Data13TeV_DoubleMuon_2016E","Data13TeV_DoubleMuon_2016F","Data13TeV_MuonEG_2016B","Data13TeV_MuonEG_2016C","Data13TeV_MuonEG_2016D","Data13TeV_MuonEG_2016E","Data13TeV_MuonEG_2016F","Data13TeV_SingleElectron_2016B","Data13TeV_SingleElectron_2016C","Data13TeV_SingleElectron_2016D","Data13TeV_SingleElectron_2016E","Data13TeV_SingleElectron_2016F","Data13TeV_SingleMuon_2016B","Data13TeV_SingleMuon_2016C","Data13TeV_SingleMuon_2016D","Data13TeV_SingleMuon_2016E","Data13TeV_SingleMuon_2016F"};
for(auto &it : dataSamples)
  datad->Add("LJets2015/2016/test/Chunks/"+it+"_*.root");
TH1F *d = new TH1F("data", "Data", nbins, low, high);
/*
for(auto &it : dataSamples) {
  auto fin = new TFile("LJets2015/2016/test/"+it+".root");
  std::cout << fin->GetName() << std::endl;
  if(TString(fin->GetName()).Contains("G") || TString(fin->GetName()).Contains("H"))
  d->Add((TH1F*)fin->Get("D0_l3d_all_meson_GH")->Clone());
  else
  d->Add((TH1F*)fin->Get("D0_l3d_all_meson_BCDEF")->Clone());
  fin->Close();
  delete fin;
}
*/
TH1F *tot = new TH1F("tot", "tot", nbins, low, high);
//TH1F *tot = new TH1F("tot", "tot", 20, 0, 1);
TH1F *b = new TH1F("b", "b-jets (t#bar{t})", nbins, low, high);
//TH1F *bw = new TH1F("bw", "b (W+jets)", nbins, low, high);
//TH1F *b = new TH1F("c", "c", 20, 0, 1);
TH1F *bw = new TH1F("bw", "b-jets (W+jets)", nbins, low, high);
TH1F *ctt = new TH1F("ctt", "c-jets (t#bar{t})", nbins, low, high);
TH1F *cw = new TH1F("ccw", "c-jets (W+jets)", nbins, low, high);
TH1F *co = new TH1F("co", "c-jets (other)", nbins, low, high);
TH1F *bo = new TH1F("bo", "b-jets (other)", nbins, low, high);
//TH1F *c = new TH1F("c", "c", 20, 0, 1);
TH1F *uds = new TH1F("uds", "uds-jets (all bkg)", nbins, low, high);
//TH1F *uds = new TH1F("uds", "uds", 20, 0, 1);

datad->Draw(TString::Format("%s>>data", val.Data()), TString::Format("meson_id==443 && d0_l3d/d0_sigmal3d>%d", ctausig), "goff");
data->Draw(TString::Format("%s>>tot", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==443 && d0_l3d/d0_sigmal3d>%d)", ctausig), "goff");
//data->Draw(TString::Format("%s>>c", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==443 && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4)", ctausig), "goff");
data->Draw(TString::Format("%s>>uds", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==443 && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)<4)", ctausig), "goff");
datatt->Draw(TString::Format("%s>>ctt", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==443  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4)", ctausig), "goff");
dataw->Draw(TString::Format("%s>>ccw", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==443  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4)", ctausig), "goff");
datao->Draw(TString::Format("%s>>co", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==443  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4)", ctausig), "goff");
datao->Draw(TString::Format("%s>>bo", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==443  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==5)", ctausig), "goff");
dataw->Draw(TString::Format("%s>>bw", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==443  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==5)", ctausig), "goff");
datatt->Draw(TString::Format("%s>>b", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==443 && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==5)", ctausig), "goff");

b->SetFillColor(920);
b->SetLineColor(kBlack);
//bw->SetFillColor(930);
//bw->SetLineColor(kBlack);
ctt->SetFillColor(kRed);
ctt->SetLineColor(kBlack);
bw->SetFillColor(865);
bw->SetLineColor(kBlack);
cw->SetFillColor(kViolet);
cw->SetLineColor(kBlack);
co->SetFillColor(kViolet);
co->SetLineColor(kBlack);
bo->SetFillColor(kRed+2);
bo->SetLineColor(kBlack);
uds->SetFillColor(kBlue);
uds->SetLineColor(kBlack);
d->SetLineColor(kBlack);
d->SetMarkerStyle(20);
b->Scale(35864);
bw->Scale(35864);
ctt->Scale(35864);
cw->Scale(35864);
co->Scale(35864);
bo->Scale(35864);
uds->Scale(35864);

auto hist = new THStack("h", "h");
uds->Add(co);
hist->Add(uds);
//hist->Add(co);
hist->Add(bo);
hist->Add(cw);
hist->Add(bw);
hist->Add(ctt);
hist->Add(b);
hist->Draw("hist");
//d->Draw("same e");
c1->BuildLegend();
hist->SetTitle("");

std::cout << b->Integral() << " " << bw->Integral() << " " << bw->Integral() / (b->Integral() + bw->Integral() + uds->Integral()) << " " << (bw->Integral() / (b->Integral() + bw->Integral() + cw->Integral() + ctt->Integral() + uds->Integral())) * sqrt(1/bw->Integral() + 1/cw->Integral() + 1/ctt->Integral() + (1/b->Integral() + 1/uds->Integral())) << std::endl;
float sum = uds->Integral() + bo->Integral() + cw->Integral() + bw->Integral() + ctt->Integral() + b->Integral();
float bkg = uds->Integral() + cw->Integral() + bw->Integral() + ctt->Integral();
float cbkg = +cw->Integral() + ctt->Integral();
std::cout << "bkg = " << bkg/sum * 100 << " +/- " << bkg/sum * 100 * (1/bkg + 1/sum) << std::endl;
std::cout << "c bkg = " << cbkg/sum * 100 << " +/- " << cbkg/sum * 100 * (1/cbkg + 1/sum) << std::endl;
}
