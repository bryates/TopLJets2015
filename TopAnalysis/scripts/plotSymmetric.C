#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include "LJets2015/2016/mtop/tdr.h"

void plotSymmetric() {
  gStyle->SetOptStat(0);
  auto c1 = setupCanvas();
  setupPad()->cd();
  auto fsym = new TFile("LJets2015/2016/test2/MC13TeV_TTJets_symmetric.root");
  auto fasym = new TFile("LJets2015/2016/test/MC13TeV_TTJets_powheg.root");


  auto symB = (TH1F*)fsym->Get("D0oJet_pt_charged_all_meson_BCDEF");
  auto symG = (TH1F*)fsym->Get("D0oJet_pt_charged_all_meson_GH");
  auto asymB = (TH1F*)fasym->Get("D0oJet_pt_charged_all_meson_BCDEF");
  auto asymG = (TH1F*)fasym->Get("D0oJet_pt_charged_all_meson_GH");

  symB->Scale(19712.86);
  asymB->Scale(19712.86);
  symG->Scale(16146.178);
  asymG->Scale(16146.178);
  symB->Add(symG);
  asymB->Add(asymG);

  auto ratio = (TH1F*)asymB->Clone("ratio");
  ratio->Divide(asymB, symB, 1/asymB->Integral(), 1/symB->Integral());
  ratio->GetYaxis()->SetTitle("(1/N Standard) / (1/N Symmetric)");

  tdr(ratio);
  ratio->Draw();
  tdr(ratio);
return;

  asymB->GetYaxis()->SetTitle("1/N dN/d#it{x}_{B} (Jets / 0.05)");
  asymB->SetLineColor(kBlue);
  symB->SetLineColor(kGreen);
  tdr(asymB);
  asymB->DrawNormalized();
  symB->DrawNormalized("same");
  tdr(asymB);

  float iniy=0.83;
  float dy=0.05;
  float ndy=4;
  TLegend *leg = new TLegend(0.21, iniy-dy*ndy, 0.53, iniy+0.05);
  leg->SetBorderSize(1);
  leg->SetFillStyle(0);      
  leg->SetTextFont(43);
  leg->SetTextSize(16);
  leg->AddEntry( asymB, "Standard", "l");
  leg->AddEntry( symB, "Symmetric", "l");
  leg->Draw();
}
