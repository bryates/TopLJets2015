#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH2F.h>

void pf_sf() {
TFile *f = new TFile("data/era2016/pf_tracks.root","RECREATE");
float b[] = {-2.4,-1.5,-0.8,-0.4,0,0.4,0.8,1.5,2.4};
float bpt[] = {0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300};
TChain *dataB = new TChain("data");
dataB->Add("LJets2015/2016/eta/Data13TeV_*_2016B_*");
dataB->Add("LJets2015/2016/eta/Data13TeV_*_2016C_*");
dataB->Add("LJets2015/2016/eta/Data13TeV_*_2016D_*");
dataB->Add("LJets2015/2016/eta/Data13TeV_*_2016E_*");
dataB->Add("LJets2015/2016/eta/Data13TeV_*_2016F_*");

int nb = sizeof(b)/sizeof(float);
int nbpt = sizeof(bpt)/sizeof(float);
TH2F *hB = new TH2F("hB","hB", nb-1, b, nbpt-1, bpt);
hB->Sumw2();
dataB->Draw("pf_pt:pf_eta>>hB","","colz");
//hB->Draw("colz text e");

TChain *dataG = new TChain("data");
dataG->Add("LJets2015/2016/eta/Data13TeV_*_2016G_*");
dataG->Add("LJets2015/2016/eta/Data13TeV_*_2016H_v*_*");

TH2F *hG = new TH2F("hG","hG", nb-1, b, nbpt-1, bpt);
hG->Sumw2();
dataG->Draw("pf_pt:pf_eta>>hG","","colz");
//hG->Draw("colz text e");

TH2F *ra = (TH2F*)hB->Clone();
ra->Divide(hB, hG, 1./19712.86, 1./16146.178);
ra->SetName("eta_pt");
ra->SetTitle("");
ra->GetXaxis()->SetTitle("N_{B-F} / N_{GH} 1/(L_{B-F} /L_{GH})PF #eta");
ra->GetYaxis()->SetTitle("N_{B-F} / N_{GH} 1/(L_{B-F} /L_{GH})PF p_{T}");
gStyle->SetPaintTextFormat("4.2g");
ra->Draw("colz text e");
ra->Write();
f->Close();
}
