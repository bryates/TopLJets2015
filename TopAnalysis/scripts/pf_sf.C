#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH2F.h>

void pf_sf() {
TString base("/afs/cern.ch/user/b/byates/TopAnalysis/");
TFile *f = new TFile(base + "data/era2016/pf_tracks.root","RECREATE");
float b[] = {-2.4,-1.5,-0.8,-0.4,0,0.4,0.8,1.5,2.4};
float bpt[] = {0,2,4,6,10,20,30,50,100,300};
//float bpt[] = {0,1,2,3,4,5,6,7,8,9,10,20,30,50,100,300};
//float bpt[] = {0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300};
TChain *dataB = new TChain("data");
dataB->Add(base+"LJets2015/2016/eta/Chunks/Data13TeV_*_2016B_*");
dataB->Add(base+"LJets2015/2016/eta/Chunks/Data13TeV_*_2016C_*");
dataB->Add(base+"LJets2015/2016/eta/Chunks/Data13TeV_*_2016D_*");
dataB->Add(base+"LJets2015/2016/eta/Chunks/Data13TeV_*_2016E_*");
dataB->Add(base+"LJets2015/2016/eta/Chunks/Data13TeV_*_2016F_*");

gStyle->SetOptStat(0);
gStyle->SetPaintTextFormat("4.2g");

int nb = sizeof(b)/sizeof(float);
int nbpt = sizeof(bpt)/sizeof(float);
TH2F *hB = new TH2F("hB","hB", nb-1, b, nbpt-1, bpt);
hB->Sumw2();
dataB->Draw("pf_pt:pf_eta>>hB","","colz");
hB->Draw("colz text e");
hB->Write();

TChain *dataG = new TChain("data");
dataG->Add(base+"LJets2015/2016/eta/Chunks/Data13TeV_*_2016G_*");
dataG->Add(base+"LJets2015/2016/eta/Chunks/Data13TeV_*_2016H_v*_*");

TH2F *hG = new TH2F("hG","hG", nb-1, b, nbpt-1, bpt);
hG->Sumw2();
dataG->Draw("pf_pt:pf_eta>>hG","","colz");
hG->Draw("colz text e");
hG->Write();

TH2F *ra = (TH2F*)hB->Clone();
ra->Divide(hB, hG, 1./19712.86, 1./16146.178);
ra->SetName("eta_pt");
ra->SetTitle("");
ra->GetXaxis()->SetTitle("N_{B-F} / N_{GH} 1/(L_{B-F} /L_{GH})PF #eta");
ra->GetYaxis()->SetTitle("N_{B-F} / N_{GH} 1/(L_{B-F} /L_{GH})PF p_{T}");
//Testing errors, ROOT seems to be correct if Sumw2 is OFF
for(int i = ra->GetYaxis()->GetNbins(); i > 0; i--) {
  for(int j = 1; j < ra->GetXaxis()->GetNbins(); j++) {
    int bin = ra->GetBin(j,i);
    float Nb = hB->GetBinContent(bin);
    float Ng = hG->GetBinContent(bin);
    float Nr = ra->GetBinContent(bin);
    float sig = Nr * sqrt(1./Nb + 1/Ng);
    ra->SetBinError(bin, sig);
  }
}
ra->Draw("colz text e");
ra->SetMinimum(0.7);
ra->SetMaximum(1.0);
ra->Write();

hB->Reset();
hG->Reset();

dataB->Draw("d0_pi_pt:d0_pi_eta>>hB","meson_id==421","colz");
hB->SetName("d0piB");
hB->SetTitle("D^{0} #pi");
hB->Draw("colz text e");
hB->Write();

dataG->Draw("d0_pi_pt:d0_pi_eta>>hG","meson_id==421","colz");
hG->SetName("d0piG");
hG->SetTitle("D^{0} #pi");
hG->Draw("colz text e");
hG->Write();
ra->Divide(hB, hG, 1./19712.86, 1./16146.178);
ra->SetName("pi_eta_pt");
ra->SetTitle("");
ra->GetXaxis()->SetTitle("N_{B-F} / N_{GH} 1/(L_{B-F} /L_{GH}) D^{0} #pi #eta");
ra->GetYaxis()->SetTitle("N_{B-F} / N_{GH} 1/(L_{B-F} /L_{GH}) D^{0} #pi p_{T}");
//Testing errors, ROOT seems to be correct if Sumw2 is OFF
for(int i = ra->GetYaxis()->GetNbins(); i > 0; i--) {
  for(int j = 1; j < ra->GetXaxis()->GetNbins(); j++) {
    int bin = ra->GetBin(j,i);
    float Nb = hB->GetBinContent(bin);
    float Ng = hG->GetBinContent(bin);
    float Nr = ra->GetBinContent(bin);
    float sig = Nr * sqrt(1./Nb + 1/Ng);
    ra->SetBinError(bin, sig);
  }
}
ra->Draw("colz text e");
ra->SetMinimum(0.7);
ra->SetMaximum(1.1);
ra->Write();

f->Close();
}
