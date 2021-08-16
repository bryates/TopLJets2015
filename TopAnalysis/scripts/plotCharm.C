#include <TChain.h>
#include <TH1.h>
#include <THStack.h>
#include "/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/2016/mtop/tdr.C"

void plotCharm(int ctausig=10) {
auto c1 = setupCanvas();
c1->cd();
TPad *p1 = new TPad("p1","p1",0.0,0.2,1.0,1.0);
//p1->SetLogy();
p1->SetRightMargin(0.05);
p1->SetLeftMargin(0.12);
p1->SetTopMargin(0.1);
p1->SetBottomMargin(0.12);
p1->Draw();
p1->cd();

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

TString val = "d0_mass";
TString units = "";
//TString val = "d0_mass";
TString title = TString::Format("D^{0} mass%s", units.Data());
//TString title = "D^{0} mass [GeV]";
int nbins = 60;
float low = 1.7;
float high = 2;
std::vector<TString> dataSamples = {"Data13TeV_DoubleEG_2016B","Data13TeV_DoubleEG_2016C","Data13TeV_DoubleEG_2016D","Data13TeV_DoubleEG_2016E","Data13TeV_DoubleEG_2016F","Data13TeV_DoubleMuon_2016B","Data13TeV_DoubleMuon_2016C","Data13TeV_DoubleMuon_2016D","Data13TeV_DoubleMuon_2016E","Data13TeV_DoubleMuon_2016F","Data13TeV_MuonEG_2016B","Data13TeV_MuonEG_2016C","Data13TeV_MuonEG_2016D","Data13TeV_MuonEG_2016E","Data13TeV_MuonEG_2016F","Data13TeV_SingleElectron_2016B","Data13TeV_SingleElectron_2016C","Data13TeV_SingleElectron_2016D","Data13TeV_SingleElectron_2016E","Data13TeV_SingleElectron_2016F","Data13TeV_SingleMuon_2016B","Data13TeV_SingleMuon_2016C","Data13TeV_SingleMuon_2016D","Data13TeV_SingleMuon_2016E","Data13TeV_SingleMuon_2016F","Data13TeV_SingleMuon_2016G","Data13TeV_SingleMuon_2016H_v2","Data13TeV_SingleMuon_2016H_v3","Data13TeV_DoubleEG_2016G","Data13TeV_DoubleEG_2016H_v2","Data13TeV_DoubleEG_2016H_v3","Data13TeV_DoubleMuon_2016G","Data13TeV_DoubleMuon_2016H_v2","Data13TeV_DoubleMuon_2016H_v3","Data13TeV_MuonEG_2016G","Data13TeV_MuonEG_2016H_v2","Data13TeV_MuonEG_2016H_v3","Data13TeV_SingleElectron_2016G","Data13TeV_SingleElectron_2016H_v2","Data13TeV_SingleElectron_2016H_v3"};      

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
TH1F *b = new TH1F("b", "b-jets (t#bar{t})", nbins, low, high);
TH1F *bw = new TH1F("bw", "b-jets (W+jets)", nbins, low, high);
TH1F *ctt = new TH1F("ctt", "c-jets (t#bar{t})", nbins, low, high);
TH1F *cw = new TH1F("ccw", "c-jets (W+jets)", nbins, low, high);
TH1F *co = new TH1F("co", "c-jets (other)", nbins, low, high);
TH1F *uds = new TH1F("uds", "uds-jets (all bkg)", nbins, low, high);
TH1F *totG = new TH1F("totG", "totG", nbins, low, high);
TH1F *bG = new TH1F("bG", "b-jets (t#bar{t})G", nbins, low, high);
TH1F *bwG = new TH1F("bwG", "b-jets (W+jets)G", nbins, low, high);
TH1F *cttG = new TH1F("cttG", "c-jets (t#bar{t})G", nbins, low, high);
TH1F *cwG = new TH1F("ccwG", "c-jets (W+jets)G", nbins, low, high);
TH1F *coG = new TH1F("coG", "c-jets (other)G", nbins, low, high);
TH1F *udsG = new TH1F("udsG", "uds-jets (all bkg)G", nbins, low, high);

datad->Draw(TString::Format("%s>>data", val.Data()), TString::Format("meson_id==421 && d0_l3d/d0_sigmal3d>%d && d0_pt/j_pt_charged>0.2", ctausig), "goff");
data->Draw(TString::Format("%s>>tot", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && d0_l3d/d0_sigmal3d>%d && d0_pt/j_pt_charged>0.2 && epoch==1)", ctausig), "goff");
data->Draw(TString::Format("%s>>totG", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && d0_l3d/d0_sigmal3d>%d && d0_pt/j_pt_charged>0.2 && epoch==2)", ctausig), "goff");
//data->Draw(TString::Format("%s>>c", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4)", ctausig), "goff");
//data->Draw(TString::Format("%s>>uds", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)<4)", ctausig), "goff");
//data->Draw("d0_pt/j_pt_charged>>tot", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421)", ctausig), "goff");
//dataw->Draw(TString::Format("%s>>bw", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==5)", ctausig), "goff");
//data->Draw("d0_pt/j_pt_charged>>b", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==5)", ctausig), "goff");
datatt->Draw(TString::Format("%s>>b", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==5 && d0_pt/j_pt_charged>0.2 && epoch==1)", ctausig), "goff");
datatt->Draw(TString::Format("%s>>bG", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==5 && d0_pt/j_pt_charged>0.2 && epoch==2)", ctausig), "goff");
datatt->Draw(TString::Format("%s>>ctt", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4 && d0_pt/j_pt_charged>0.2 && epoch==1)", ctausig), "goff");
datatt->Draw(TString::Format("%s>>cttG", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4 && d0_pt/j_pt_charged>0.2 && epoch==2)", ctausig), "goff");
dataw->Draw(TString::Format("%s>>ccw", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4 && d0_pt/j_pt_charged>0.2 && epoch==1)", ctausig), "goff");
dataw->Draw(TString::Format("%s>>ccwG", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4 && d0_pt/j_pt_charged>0.2 && epoch==2)", ctausig), "goff");
datao->Draw(TString::Format("%s>>co", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4 && d0_pt/j_pt_charged>0.2 && epoch==1)", ctausig), "goff");
datao->Draw(TString::Format("%s>>coG", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4 && d0_pt/j_pt_charged>0.2 && epoch==2)", ctausig), "goff");
dataw->Draw(TString::Format("%s>>bw", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==5 && d0_pt/j_pt_charged>0.2 && epoch==1)", ctausig), "goff");
dataw->Draw(TString::Format("%s>>bwG", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==5 && d0_pt/j_pt_charged>0.2 && epoch==2)", ctausig), "goff");
//data->Draw("d0_pt/j_pt_charged>>c", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4 && d0_pt/j_pt_charged>0.2 && epoch==1)", ctausig), "goff");
//data->Draw("d0_pt/j_pt_charged>>cG", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==4 && d0_pt/j_pt_charged>0.2 && epoch==2)", ctausig), "goff");
data->Draw(TString::Format("%s>>uds", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==0 && d0_pt/j_pt_charged>0.2 && epoch==1)", ctausig), "goff");
data->Draw(TString::Format("%s>>udsG", val.Data()), TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==0 && d0_pt/j_pt_charged>0.2 && epoch==2)", ctausig), "goff");
//data->Draw("d0_pt/j_pt_charged>>uds", TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421  && d0_l3d/d0_sigmal3d>%d && abs(j_hadflav)==0)", ctausig), "goff");

b->SetFillColor(920);
b->SetLineColor(kBlack);
//bw->SetFillColor(930);
//bw->SetLineColor(kBlack);
ctt->SetFillColor(kRed);
ctt->SetLineColor(kBlack);
bw->SetFillColor(kRed+1);
bw->SetLineColor(kBlack);
cw->SetFillColor(kViolet);
cw->SetLineColor(kBlack);
co->SetFillColor(kViolet);
co->SetLineColor(kBlack);
uds->SetFillColor(kBlue);
uds->SetLineColor(kBlack);
d->SetLineColor(kBlack);
d->SetMarkerStyle(20);
/*
b->Scale(19712.86);
bG->Scale(16146.178);
bw->Scale(19712.86);
bwG->Scale(16146.178);
ctt->Scale(19712.86);
cttG->Scale(16146.178);
cw->Scale(19712.86);
cwG->Scale(16146.178);
co->Scale(19712.86);
coG->Scale(16146.178);
uds->Scale(19712.86);
udsG->Scale(16146.178);
tot->Scale(19712.86);
totG->Scale(16146.178);
tot->Scale(1.06);
totG->Scale(1.13);
*/

b->Add(bG);
bw->Add(bwG);
ctt->Add(cttG);
cw->Add(cwG);
co->Add(coG);
uds->Add(udsG);
tot->Add(totG);

b->Scale(1.27);
bw->Scale(1.27);
ctt->Scale(1.27);
cw->Scale(1.27);
co->Scale(1.27);
uds->Scale(1.27);
tot->Scale(1.165);

float scale = b->Integral() + bw->Integral() + ctt->Integral() + cw->Integral() + co->Integral() + uds->Integral();
std::cout << tot->Integral() / scale << std::endl;
std::cout << d->Integral() / scale << std::endl;

b->Scale(1.165);
bw->Scale(1.165);
ctt->Scale(1.165);
cw->Scale(1.165);
co->Scale(1.165);
uds->Scale(1.165);

auto hist = new THStack("h", "h");
uds->Add(co);
hist->Add(uds);
//hist->Add(co);
hist->Add(cw);
hist->Add(bw);
hist->Add(ctt);
hist->Add(b);
hist->SetMinimum(1);
tot->SetMinimum(1);
tdr(tot,0,0);
tot->Draw();
//hist->Draw("hist");
hist->Draw("same hist");
tdr(tot,0,0);
//tot->Draw("same e2");
d->Draw("same e");
auto leg = p1->BuildLegend();
auto y1 = leg->GetY1();
auto y2 = leg->GetY2();
std::cout << "y1=" << y1 << " y2=" << y2 << std::endl;
float iniy=0.53;
float dy=0.05;
float ndy=5;
leg->SetX1(0.2);
leg->SetX2(0.5);
leg->SetY1(iniy-dy*ndy);
leg->SetY2(iniy+0.05);
hist->SetTitle("");
gStyle->SetOptStat(0);
//gPad->SetLogy();
gPad->Modified();
gPad->Update();
std::cout << hist->GetMaximum() << std::endl;

float sum(0);
for(int i = 0; i < hist->GetNhists(); i++)
  sum += ((TH1F*)((TList*)hist->GetHists())->At(i))->Integral();


std::cout << b->Integral() << " " << bw->Integral() << " " << bw->Integral() / (b->Integral() + bw->Integral() + uds->Integral()) << " " << (bw->Integral() / (b->Integral() + bw->Integral() + cw->Integral() + ctt->Integral() + uds->Integral())) * sqrt(1/bw->Integral() + 1/cw->Integral() + 1/ctt->Integral() + (1/b->Integral() + 1/uds->Integral())) << std::endl;
std::cout << b->Integral() << " " << bw->Integral() << " " << bw->Integral() / (b->Integral() + bw->Integral() + uds->Integral()) << " " << (bw->Integral() / (b->Integral() + bw->Integral() + cw->Integral() + ctt->Integral() + uds->Integral())) * sqrt(1/sum) << std::endl;

std::cout << d->Integral() / tot->Integral() << std::endl;
float bkg = uds->Integral() + cw->Integral() + bw->Integral() + ctt->Integral();
float cbkg = +cw->Integral() + ctt->Integral();
std::cout << "bkg = " << bkg/sum * 100 << " +/- " << bkg/sum * 100 * (1/bkg + 1/sum) << std::endl;
std::cout << "c bkg = " << cbkg/sum * 100 << " +/- " << cbkg/sum * 100 * (1/cbkg + 1/sum) << std::endl;
gStyle->SetOptStat(0);
c1->cd();
auto p2 = new TPad("p2","p2",0.0,0.0,1.0,0.26);
gStyle->SetPaintTextFormat("%0.2f");
p2->SetBottomMargin(0.4);
p2->SetRightMargin(0.05);
p2->SetLeftMargin(0.12);
p2->SetTopMargin(0.1);
p2->SetGridx(false);
p2->SetGridy(false);
p2->Draw();
p2->cd();
auto ratio = (TH1F*)d->Clone("");
ratio->Divide(d, tot, 1, d->Integral() / tot->Integral());
ratio->SetTitle("");
ratio->GetXaxis()->SetTitle(title);
ratio->GetYaxis()->SetTitle(TString::Format("Jets / [%d%s]", nbins, units.Data()));
ratio->GetXaxis()->SetTitleSize(0.1);
ratio->GetXaxis()->SetLabelSize(0.1);
ratio->GetYaxis()->SetTitleSize(0.1);
ratio->GetYaxis()->SetLabelSize(0.1);
ratio->GetYaxis()->SetRangeUser(.8,1.2);
ratio->GetYaxis()->SetNdivisions(4);
ratio->Draw();
}
