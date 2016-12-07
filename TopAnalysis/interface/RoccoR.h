#ifndef ElectroWeakAnalysis_RoccoR
#define ElectroWeakAnalysis_RoccoR
#include "TRandom3.h"
#include "TMath.h"

#include "iostream"
#include "fstream"
#include "sstream"

using namespace std;
using std::cout;
using std::endl;
using std::vector;

struct CrystalBall{
    double pi;
    double SPiO2;
    double S2;

    double m;
    double s;
    double a;
    double n;

    double B;
    double C;
    double D;
    double N;

    double NA;
    double Ns;
    double NC;
    double F;
    double G;
    double k;

    double cdfMa;
    double cdfPa;

  CrystalBall();
  CrystalBall(double m_, double s_, double a_, double n_);
  void init(double m_, double s_, double a_, double n_);
  double pdf(double x);
  double cdf(double x);
  double invcdf(double u);
};

class RocRes{
    private:
	static const int NMAXETA=12;
	static const int NMAXTRK=12;

	int NETA;
	int NTRK;
	int NMIN;

	double BETA[NMAXETA+1];
	double ntrk[NMAXETA][NMAXTRK+1];
	double dtrk[NMAXETA][NMAXTRK+1];

	double width[NMAXETA][NMAXTRK];
	double alpha[NMAXETA][NMAXTRK];
	double power[NMAXETA][NMAXTRK];

	double rmsA[NMAXETA][NMAXTRK];
	double rmsB[NMAXETA][NMAXTRK];
	double rmsC[NMAXETA][NMAXTRK];

	double kDat[NMAXETA];
	double kRes[NMAXETA];

	TRandom3 random; //only used while deriving corrections now
	
	int getBin(double x, const int NN, const double *b);


    public:
	enum TYPE {MC, Data, Extra};

	CrystalBall  cb[NMAXETA][NMAXTRK]; 

	RocRes();
	int getEtaBin(double feta);
	int getNBinDT(double v, int H);
	int getNBinMC(double v, int H);
	void dumpParams();
	void init(string filename);

	~RocRes(){}

	double Sigma(double pt, int H, int F);
	double kSpreadDet(double gpt, double rpt, double eta, int nlayers, double w);
	double kSmearDet(double pt, double eta, TYPE type, double v, double u);
	double kExtraDet(double pt, double eta, int nlayers, double u, double w);
	void fillFitData(int &H, int &F, int &D, double &xmc, double &xdt, double &Rmc, double &Rdt, double pt, double eta);
};


class RocOne{
    private:
	enum TYPE{MC, DT};
	static const int NMAXETA=24;
	static const int NMAXPHI=16;
	static const double MPHI;

	int NETA;
	int NPHI;

	double BETA[NMAXETA+1];
	double DPHI;

	double M[2][NMAXETA][NMAXPHI];
	double A[2][NMAXETA][NMAXPHI];
	double D[2][NMAXETA];

	RocRes RR;

	int getBin(double x, const int NN, const double *b);
	int getBin(double x, const int nmax, const double xmin, const double dx);

    public:
	RocOne();
	~RocOne(){}
	RocOne(string filename, int iTYPE=0, int iSYS=0, int iMEM=0);
	bool checkSYS(int iSYS, int iMEM, int kSYS=0, int kMEM=0);
	bool checkTIGHT(int iTYPE, int iSYS, int iMEM, int kTYPE=0, int kSYS=0, int kMEM=0);
	void init(string filename, int iTYPE=0, int iSYS=0, int iMEM=0);

	double kScaleDT(int Q, double pt, double eta, double phi);
	double kScaleMC(int Q, double pt, double eta, double phi, double kSMR=1);
	double kScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers, double u, double w);
	double kScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers, double w);
	double kGenSmearDet(double pt, double eta, double v, double u);
};


class RoccoR{
    public:
	static const int Nset=5;
	enum TYPE{Default, Stat, ModPt, CorDM, FitDM};
	static const int Nmem[Nset];

	RoccoR(); //temporary, will change 
	RoccoR(std::string dirname); //temporary, will change 
	~RoccoR();
	
	double kScaleDT(int Q, double pt, double eta, double phi, TYPE T, int m);

	double kScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers);	      //only for default, non-reproducible
	double kScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers); //only for default, non-reproducible

	double kScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers, double u, double w, TYPE T, int m);  
	double kScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers, double w, TYPE T, int m); 

	std::vector<std::vector<double> > kkScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers);
	std::vector<std::vector<double> > kkScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers);

    private:
	TRandom3 random;
	RocOne *RC[Nset][100];

};


#endif

