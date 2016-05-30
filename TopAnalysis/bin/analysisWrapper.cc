#include <iostream>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"
#include "TopLJets2015/TopAnalysis/interface/TOPWidth.h"

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
       << "\t --runSysts - activate running systematics" << endl
       << "\t --norm     - ROOT file with normalization weights" << endl
       << "\t --normTag  - normalization tag" << endl
       << "\t --method   - method to run" << endl;
}

//
int main(int argc, char* argv[])
{
  //get input arguments
  TString in(""),out(""),norm(""),normTag(""),method("");
  bool runSysts(false);
  int channel(0),charge(0),flav(0);
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help") !=string::npos)                     { printHelp(); return -1;} 
    else if(arg.find("--runSysts")!=string::npos )            { runSysts=true;  }
    else if(arg.find("--channel")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%d",&channel); i++;}
    else if(arg.find("--charge")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&charge); i++;}
    else if(arg.find("--flav")!=string::npos && i+1<argc)     { sscanf(argv[i+1],"%d",&flav); i++;}
    else if(arg.find("--in")!=string::npos && i+1<argc)       { in=argv[i+1]; i++;}
    else if(arg.find("--out")!=string::npos && i+1<argc)      { out=argv[i+1]; i++;}
    else if(arg.find("--normTag")!=string::npos && i+1<argc)  { normTag=argv[i+1]; i++;}
    else if(arg.find("--norm")!=string::npos && i+1<argc)     { norm=argv[i+1]; i++;}
    else if(arg.find("--method")!=string::npos && i+1<argc)   { method=argv[i+1]; i++;}
  }

  //open normalization file
  TH1F *normH=0;
  TFile *normF=TFile::Open(norm);
  if(normF)
    {
      normH=(TH1F *)normF->Get(normTag);
      if(normH) normH->SetDirectory(0);
      normF->Close();
    }
  if(normH==0)
    {
      cout << "Check normalization file (" << norm << ") and tag (" << normTag << ")" << endl;
      printHelp();
      return -1;
    }
  
  //check input/output
  if(in=="" || out=="")
    {
      cout << "Check input/output=" << in << "/" << out << endl;
      printHelp();
      return -1;
    }

  //check method to run
  if(method=="TOP-16-006::RunTop16006")    RunTop16006(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts);
  else if(method=="TOPWidth::RunTopWidth") RunTopWidth(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts);
  else
    {
      cout << "Check method=" << method <<endl;
      printHelp();
      return -1;
    }

  //all done
  return 0;
}  





