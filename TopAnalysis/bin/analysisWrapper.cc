#include <iostream>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
//#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"
#include "TopLJets2015/TopAnalysis/interface/TOPWidth.h"
#include "TopLJets2015/TopAnalysis/interface/TopSync.h"
//#include "TopLJets2015/TopAnalysis/interface/Run5TeVAnalysis.h"

#include "TH1F.h"
#include "TFile.h"

using namespace std;

//
void printHelp()
{
  cout << "analysisWrapper options are:" << endl
       << "\t --in - input file" << endl
       << "\t --out - output file" << endl
       << "\t --channel - channel to analyze" << endl
       << "\t --charge  - charge selection to apply" << endl
       << "\t --flav    - flavour selection to apply" << endl
       << "\t --runSysts - activate running systematics (1 for up, -1 for down)" << endl
       << "\t --era      - era directory to use for corrections, uncertainties" << endl
       << "\t --runPeriod   - runs to include" << endl
       << "\t --normTag  - normalization tag" << endl
       << "\t --method   - method to run" << endl
       << "\t --verbose   - verbose" << endl
       << "\t --rbFit   - run with fitted rB weights (1 for total fit, 2 for each)" << endl;
}

//
int main(int argc, char* argv[])
{
  //get input arguments
  TString in(""),out(""),era(""),runPeriod("BCDEF"),normTag(""),method(""),weight("genweights.root");
  bool verbose(false);
  int channel(0),charge(0),flav(0);
  short runSysts(0),rbFit(0);
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help") !=string::npos)                     { printHelp(); return -1;} 
    else if(arg.find("--runSysts")!=string::npos && i+1<argc )   { sscanf(argv[i+1],"%hd",&runSysts); i++; }
    else if(arg.find("--channel")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%d",&channel); i++;}
    else if(arg.find("--charge")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&charge); i++;}
    else if(arg.find("--flav")!=string::npos && i+1<argc)     { sscanf(argv[i+1],"%d",&flav); i++;}
    else if(arg.find("--in")!=string::npos && i+1<argc)       { in=argv[i+1]; i++;}
    else if(arg.find("--out")!=string::npos && i+1<argc)      { out=argv[i+1]; i++;}
    else if(arg.find("--normTag")!=string::npos && i+1<argc)  { normTag=argv[i+1]; i++;}
    else if(arg.find("--era")!=string::npos && i+1<argc)      { era=argv[i+1]; i++;}
    else if(arg.find("--runPeriod")!=string::npos && i+1<argc)      { runPeriod=argv[i+1]; i++;}
    else if(arg.find("--method")!=string::npos && i+1<argc)   { method=argv[i+1]; i++;}
    else if(arg.find("--verbose")!=string::npos )             { verbose=true;  }
    else if(arg.find("--rbFit")!=string::npos && i+1<argc )   { sscanf(argv[i+1],"%hd",&rbFit); i++; }
  }

  //open normalization file
  TH1F *normH=0;
  TString normUrl(era+"/"+weight);
  //TString normUrl(era+"/genweights.root");
  if(TString(normTag).Contains("MC13TeV_TTJets_m1")) normUrl.ReplaceAll("genweights","genweights_mass");
  else if(TString(normTag).Contains("MC13TeV_TTJets_fsr")) normUrl.ReplaceAll("genweights","genweights_syst");
  else if(TString(normTag).Contains("MC13TeV_TTJets_isr")) normUrl.ReplaceAll("genweights","genweights_syst");
  else if(TString(normTag).Contains("MC13TeV_TTJets_ue")) normUrl.ReplaceAll("genweights","genweights_syst");
  else if(TString(normTag).Contains("MC13TeV_TTJets_erdON")) normUrl.ReplaceAll("genweights","genweights_syst");
  else if(TString(normTag).Contains("MC13TeV_TTJets_erdON")) normUrl.ReplaceAll("genweights","genweights_syst");
  else if(TString(normTag).Contains("MC13TeV_TTJets_hdamp")) normUrl.ReplaceAll("genweights","genweights_syst");
  else if(TString(normTag).Contains("MC13TeV_TTJets_herwig")) normUrl.ReplaceAll("genweights","genweights_syst");
  else if(TString(normTag).Contains("MC13TeV_GluonMove_erdON")) normUrl.ReplaceAll("genweights","genweights_syst");
  else if(TString(normTag).Contains("MC13TeV_QCD_erdON")) normUrl.ReplaceAll("genweights","genweights_syst");
  //else if(TString(normTag).Contains("MC13TeV_DY50toInf")) normUrl.ReplaceAll("genweights","genweights_ext");
  //else if(TString(normTag).Contains("MC13TeV_DY")) normUrl.ReplaceAll("genweights","genweights_DY_madgraph");
  //else if(TString(normTag).Contains("MC13TeV_DY10to50")) normUrl.ReplaceAll("genweights","genweights_DY_madgraph");
  gSystem->ExpandPathName(normUrl);
  TFile *normF=TFile::Open(normUrl);
  if(normF)
    {
      cout << "Using normalization file " << normUrl << endl;
      normH=(TH1F *)normF->Get(normTag);
      if(normH) normH->SetDirectory(0);
      normF->Close();
    }
  if(normH==0)
    {
      cout << "Check normalization file " << normUrl
	   << " and tag (" << normTag << ")" << endl
	   << "Will run without any" << endl;
      printHelp();
      //return -1;
    }
  
  //check input/output
  if(in=="" || out=="")
    {
      cout << "Check input/output=" << in << "/" << out << endl;
      printHelp();
      return -1;
    }

  //check method to run
  //if(method=="TOP-16-006::RunTop16006")    RunTop16006(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era);
  //else if(method=="TOPWidth::RunTopWidth") RunTopWidth(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era);
  if(method=="TOP::RunTop") RunTop(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era,runPeriod,verbose);
  else if(method=="TOP::RunTopKalman") RunTopKalman(in,out,channel,charge,FlavourSplitting(flav),normH,era,runPeriod,verbose,runSysts,rbFit);
  else if(method=="TOP::RunTopAOD") RunTopAOD(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era,runPeriod,verbose);
  else if(method=="TOP::RunTopSync") RunTopSync(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era,runPeriod,verbose);
  //else if(method=="Run5TeVAnalysis::Run5TeVAnalysis") Run5TeVAnalysis(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era);
  else
    {
      cout << "Check method=" << method <<endl;
      printHelp();
      return -1;
    }

  //all done
  return 0;
}  





