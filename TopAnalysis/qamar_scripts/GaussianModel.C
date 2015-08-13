// make a simple Gaussian model



//#include "RooWorkspace.h"
//#include "RooRealVar.h"
//#include "RooAbsPdf.h"
//#include "RooDataSet.h"
//#include "RooPlot.h"

//#include "TCanvas.h"

using namespace RooFit; 

void GaussianModel(int n = 1000) { 


   RooWorkspace w("w");
// define a Gaussian pdf
w.factory("Gaussian:pdf(x[-10,10],mu[1,-1000,1000],sigma[2,0,1000])");   

   RooAbsPdf * pdf = w.pdf("pdf");   // access object from workspace
   RooRealVar * x = w.var("x");   // access object from workspace
   

// generate n gaussian measurement for a Gaussian(x, 1, 1);
   RooDataSet * data = pdf->generate( *x, n); 

   data->SetName("data");  

// RooFit plotting capabilities    
   RooPlot * pl = x->frame(Title("Gaussian Fit")); 
   data->plotOn(pl); 
   pl->Draw(); 


// now fit the data set  
   pdf->fitTo(*data, Minimizer("Minuit2","Migrad") ); 
 

// plot the pdf on the same RooPlot object we have plotted the data  
  pdf->plotOn(pl);
   pdf->paramOn(pl, Layout(0.6,0.9,0.85));

   pl->Draw();

// import data in workspace (IMPORTANT for saving it ) 
 w.import(*data); 

   w.Print();

// write workspace in the file (recreate file if already existing) 
   w.writeToFile("GaussianModel.root", true);

   cout << "model written to file " << endl;

} 
