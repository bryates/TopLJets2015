//#include <TFile.h>
//#include <iostream>
//#include <TH1.h>

//#include <vector>

using namespace RooFit;
using namespace RooStats;

//void poison() {}

//gROOT->Reset();

void CountingModel(  int nobs = 11590,           	// number of observed events
		     double b = -11304,           	// number of background events
//                     double xec = 424.5,         	// cross-section
                     double A = 0.088,          	// Acceptance
                     double E = 0.224,         		// sfficiency
                     double E3 = 0.3415,         		// sfficiency
                     double L = 100,               	// Lomonisty
                     double sigmab = 0.3,          	// relative uncertainty in b
//                     double sigmaxec = 0.177,      	// relative uncertainty in xec
                     double sigmaA = 0.0014,      	// relative uncertainty in A
                     double sigmaE = 0.013,       	// relative uncertainty in E
                     double sigmaL = 0.025,         	// relative uncertainty in L
                     double sigmanobs = 0.0093 )        // relative uncertainty in nobs
{
   RooWorkspace w("w");

// make Poisson model * Gaussian constraint
   w.factory("prod:exp(A[0.088,0,100000],E[0.224,0,100000],E3[0.3415,0,100000],L[100,0,100000])"); 
//   w.factory("sum:nexp(xec[23232,0,100000],b[-19496,0,100000])"); 
   w.factory("sum:nexp(s[11590,0,100000],b[-11304,0,100000])"); 

// Poisson of (n | s+b)
   w.factory("Poisson:pdf(nobs[0,1000000],nexp)");
//   w.factory("Gaussian:constraint(b0[19496,100000],b,sigmab[0.3])");
   w.factory("Gaussian:constraint(b0[-11304,100000],b,sigmab[0.3])");
//   w.factory("Gaussian:constraint_xec(xec0[424.5,100000],xec,sigmaxec[0.177])");
   w.factory("Gaussian:constraint_acc(A0[0.088,100000],A,sigmaA[0.0014])");
   w.factory("Gaussian:constraint_eff(E0[0.224,100000],E,sigmaE[0.013])");
   w.factory("Gaussian:constraint_lum(L0[100,100000],L,sigmaL[0.025])");
   w.factory("Gaussian:constraint_obs(nobs0[23232,100000],nobs,sigmanobs[0.0093])");
   w.factory("PROD:model(pdf,constraint,constraint_acc,constraint_eff,constraint_lum,constraint_obs)"); 
//   w.factory("PROD:model(nexp,1/exp)"); 



   w.var("b0")->setVal(b);
   w.var("b0")->setConstant(true); // needed for being treated as global observables
   w.var("sigmab")->setVal(sigmab*b); 

//   w.var("xec0")->setVal(xec);
//   w.var("xec0")->setConstant(true); // needed for being treated as global observables
//   w.var("sigmaxec")->setVal(sigmaxec*xec);

   w.var("A0")->setVal(A);
   w.var("A0")->setConstant(true); // needed for being treated as global observables
   w.var("sigmaA")->setVal(sigmaA*A);

   w.var("E0")->setVal(E);
   w.var("E0")->setConstant(true); // needed for being treated as global observables
   w.var("sigmaE")->setVal(sigmaE*E);

   w.var("L0")->setVal(L);
   w.var("L0")->setConstant(true); // needed for being treated as global observables
   w.var("sigmaL")->setVal(sigmaL*L);

   w.var("nobs0")->setVal(nobs);
   w.var("nobs0")->setConstant(true); // needed for being treated as global observables
   w.var("sigmanobs")->setVal(sigmanobs*nobs);



   ModelConfig mc("ModelConfig",&w);
   mc.SetPdf(*w.pdf("model"));
//   mc.SetParametersOfInterest(*w.var("s"));
   mc.SetParametersOfInterest(*w.var("xec"));
   mc.SetObservables(*w.var("nobs"));
   mc.SetNuisanceParameters(*w.var("b")); 

   mc.SetNuisanceParameters(*w.var("A")); 
   mc.SetNuisanceParameters(*w.var("E")); 
// mc.SetNuisanceParameters(*w.var("b")); 


// these are needed for the hypothesis tests
  
//   mc.SetSnapshot(*w.var("s"));
   mc.SetSnapshot(*w.var("xec"));
   mc.SetGlobalObservables(*w.var("b0"));

//   mc.SetGlobalObservables(*w.var("xec0"));
   mc.SetGlobalObservables(*w.var("A0"));
   mc.SetGlobalObservables(*w.var("E0"));
   mc.SetGlobalObservables(*w.var("L0"));
   mc.SetGlobalObservables(*w.var("L"));
   mc.SetGlobalObservables(*w.var("nobs0"));

   mc.Print(); 

// import model in the workspace 
   w.import(mc); 


// make data set with the namber of observed events
RooDataSet data("data","", *w.var("nobs"));
//RooDataSet data("data","", *w.var("nobs"));
   w.var("nobs")->setVal(nobs);
   data.add(*w.var("nobs") );
// import data set in workspace and save it in a file
   w.import(data);

   w.Print();


   w.pdf("model")->fitTo( *w.data("data") );

   /*
   ProfileLikelihoodCalculator plc(*w.data("data"), mc);
   float confLevel(0.68);
   plc.SetConfidenceLevel(confLevel);
   LikelihoodInterval* interval = plc.GetInterval();
   
   // print out the iterval on the first Parameter of Interest
   RooRealVar* firstPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
   cout << "\n" << confLevel*100 << "% interval on " <<firstPOI->GetName()<<" is : ["<<
     interval->LowerLimit(*firstPOI) << ", "<<
     interval->UpperLimit(*firstPOI) <<"] "<<endl;

   // make a plot
   //   LikelihoodIntervalPlot plot(interval);
   //plot.SetNPoints(50);  // do not use too many points, it could become very slow for some models
   // plot.Draw("");  // use option TF1 if too slow (plot.Draw("tf1"
   */

   TString fileName = "CountingModel.root"; 

// write workspace in the file (recreate file if already existing)  
   w.writeToFile(fileName, true);

} 
