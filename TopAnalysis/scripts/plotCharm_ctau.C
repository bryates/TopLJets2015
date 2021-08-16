#include <TCanvas.h>
#include <TChain.h>
#include <TH1.h>
#include <THStack.h>

void plotCharm_ctau(int ctausig=10) {
auto c1 = new TCanvas("c1", "c1");
auto data = new TChain("data");
std::vector<TString> mcSamples = { "MC13TeV_DY10to50","MC13TeV_DY50toInf","MC13TeV_SingleT_t","MC13TeV_SingleT_tW","MC13TeV_SingleTbar_t","MC13TeV_SingleTbar_tW","MC13TeV_TTJets_powheg","MC13TeV_TTWToLNu","MC13TeV_TTWToQQ","MC13TeV_TTZToLLNuNu","MC13TeV_TTZToQQ","MC13TeV_W1Jets","MC13TeV_W2Jets","MC13TeV_W3Jets","MC13TeV_W4Jets","MC13TeV_WWTo2L2Nu","MC13TeV_WWToLNuQQ","MC13TeV_WZ","MC13TeV_ZZ" };
for(auto &it : mcSamples)
  data->Add("LJets2015/2016/test/Chunks/"+it+"_*.root");
//data->Add("LJets2015/2016/test/Chunks/MC13TeV_TTJets_powheg_*");
TH1F *tot = new TH1F("tot", "tot", 200, 0, 1);
TH1F *b = new TH1F("b", "b", 200, 0, 1);
TH1F *c = new TH1F("c", "c", 200, 0, 1);
TH1F *uds = new TH1F("uds", "uds", 200, 0, 1);

data->Draw("d0_l3d>>tot", TString::Format("sfs*topptwgt*puwgt*(meson_id==421 && d0_l3d/d0_sigmal3d<%d)", ctausig), "goff");
//data->Draw("d0_l3d>>tot", TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && d0_l3d/d0_sigmal3d<%d)", ctausig), "goff");
//data->Draw("d0_l3d/d0_sigmal3d>>tot", TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && d0_l3d/d0_sigmal3d>%d)", ctausig), "goff");
data->Draw("d0_l3d>>b", TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && abs(j_hadflav)==5 && d0_l3d/d0_sigmal3d>%d)", ctausig), "goff");
//data->Draw("d0_l3d/d0_sigmal3d>>b", TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && abs(j_hadflav)==5 && d0_l3d/d0_sigmal3d>%d)", ctausig), "goff");
data->Draw("d0_l3d>>c", TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && abs(j_hadflav)==4 && d0_l3d/d0_sigmal3d>%d)", ctausig), "goff");
//data->Draw("d0_l3d/d0_sigmal3d>>c", TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && abs(j_hadflav)==4 && d0_l3d/d0_sigmal3d>%d)", ctausig), "goff");
data->Draw("d0_l3d>>uds", TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && abs(j_hadflav)==0 && d0_l3d/d0_sigmal3d>%d)", ctausig), "goff");
//data->Draw("d0_l3d/d0_sigmal3d>>uds", TString::Format("xsec*norm*sfs*topptwgt*puwgt*(meson_id==421 && abs(j_hadflav)==0 && d0_l3d/d0_sigmal3d>%d)", ctausig), "goff");

b->SetFillColor(920);
b->SetLineColor(kBlack);
c->SetFillColor(kRed);
c->SetLineColor(kBlack);
uds->SetFillColor(kRed+1);
uds->SetLineColor(kBlack);
b->Scale(35864);
c->Scale(35864);
uds->Scale(35864);

gPad->SetLogy();
auto hist = new THStack("h", "h");
hist->Add(uds);
hist->Add(c);
hist->Add(b);
hist->Draw("hist");
c1->BuildLegend();
hist->SetTitle("");

std::cout << b->Integral() << " " << c->Integral() << " " << c->Integral() / (b->Integral() + c->Integral() + uds->Integral()) << std::endl;
std::cout << sqrt(b->Integral()) << " " << sqrt(c->Integral()) << " " << (c->Integral() / (b->Integral() + c->Integral() + uds->Integral())) * sqrt(1/c->Integral() + (1/b->Integral() + 1/uds->Integral())) << std::endl;
std::cout << tot->Integral() << std::endl;
}
