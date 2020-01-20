/*
  ctime_cor.cc
  This is a macro to optimize correction parameters for ctime
  
  Toshiyuki Gogami, Sep 25, 2019
*/

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TObjArray.h"
#include "TApplication.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TVector3.h"

#include <fstream>
#include <iostream>
using namespace std;

const double me = 0.000511;
const double mk = 0.493677;
const double mp = 0.938272;
const double mL = 1.115683;
const double mS = 1.192642;
const double HTkin_mom_scale = 2.218/2.100;

extern double CalcMM(double ee, double* pvec_ep, double* pvec_k, double mt);

extern double calcf2t_zt(double* P, 
			 double xf, double xpf,
			 double yf, double ypf);
extern double calcf2t_3rd(double*, 
			  double, double, 
			  double, double, 
			  double, int);
extern double calcf2t_4th(double*,
			  double, double, 
			  double, double,
			  double, int);
extern double calcf2t_5th(double*,
			  double, double, 
			  double, double,
			  double, int);
extern double tune(double* pa, int j);
extern void fcn(int &nPar, double* /*grad*/, 
		double &fval, double* param, int /*iflag*/);

// =================================================== //
// ==== Offset and scaling factors for matrices ====== //
// =================================================== //
const double  XFPm=-0.7,  XpFPm=-0.15; 
const double  YFPm=-0.05, YpFPm=-0.18;
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74; 
const double  XFPr=1.3,   XpFPr=0.27; 
const double  YFPr=0.1,   YpFPr=0.10; 
const double  Xptr=0.15,  Yptr=0.08, Momr=0.18; 
const double  PLm = 25.4, PLr=0.7; 
const double  Ztm = -0.15,Ztr=0.35; 


//const int npeak = 1;
double pcent = 0.0;
double pcent_real = 0.0;
double selection_width = 1.0;
//double selection_width_2[npeak2] = {5.0};
//int nL1, nS, nL2;

//const int nParamT  = 56;  // Number of parameters; 3rd
const int nParamT  = 112;  // Number of parameters (56 x 2 ); 3rd
//const int nParamT = 126;  // Number of parameters
//const int nParamT  = 252;  // Number of parameters (126 x 2); 4th
//const int nParamT2 = 504;  // Number of parameters (252 x 2); 5th
const int nmax = 10000;    // Number of events used for tuning

double x[nmax],   y[nmax]; // Reft HRS
double xp[nmax], yp[nmax];
//double tgang_xp[nmax], tgang_yp[nmax];
double x2[nmax], y2[nmax]; // Light HRS
double xp2[nmax],yp2[nmax];
double coincoin[nmax];
//double tgang_xp2[nmax], tgang_yp2[nmax];
double avz[nmax];          // averated vertex z
double beam_mom[nmax];
int tune_id[nmax];

int peak_flag[nmax];
int ntune_event = 0;
double OptPar[nParamT]; 

const double hrs_ang = 13.2 * 3.14159 / 180.;

int main(int argc, char** argv){
  TApplication* app = new TApplication("app", &argc, argv);
  int nite = 0;  // The number of tuning iteration
  char RootFileName1[500];
  char RootFileName2[500];
  if(argc==1){
    cout << " nite=0; no tuning" << endl;
    sprintf(RootFileName1,"h2_20190919.root");
    nite = 0;
  }
  else if(argc==2){
    nite = atoi(argv[1]);
    sprintf(RootFileName1,"h2_20190919.root");
  }
  else {
    cout << " Please type as follows " << endl;
    cout << " ./ctime_cor  (NITE)  (ROOTFILE)" << endl;
    return 99;
  }

    //const int nite = 1;  // The number of tuning iteration
  
  // =================================== //
  // ======== General conditions ======= //
  // =================================== //
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  // ======================================== //
  // ======= Opening a ROOT file ============ //
  // ======================================== //
  TFile* f1 = new TFile(RootFileName1);
  TTree* t1 = (TTree*)f1->Get("tree");
  //TFile* f2 = new TFile(RootFileName2);
  //TTree* t2 = (TTree*)f2->Get("tree");
  //Double_t trig1;
  double ent = t1->GetEntries();
  //ent = 50000; // for test
  //if(nite>0) ent = 400000;
  
  cout << endl;
  cout << " The number of events read from a ROOT file: " << ent << endl;
  int evshift = 0;
  //evshift = 0;
  
  const int max = 100;
  double Xpt[max], Ypt[max];
  double Xt, Yt;
  double Zt[max];
  double XFP[max],  YFP[max];
  double XpFP[max], YpFP[max];
  
  //double rtime_s0[max], ltime_s0[max];
  //double rtime_s2[max], ltime_s2[max];
  //double rtime[max], ltime[max];
  //double rpathl[max], lpathl[max];
  //double rpathl_s2[max], lpathl_s2[max];
  double a1, a2;
  double a1_2, a2_2;
  double mom1[max], mom2[max];
  //const int f1n = 64;
  //double rf1tdc[f1n];
  //double lf1tdc[f1n];
  double rvz[max], lvz[max];
  double rvz_2[max], lvz_2[max];
  double th1[max], ph1[max];
  double th2[max], ph2[max];
  double th1_2[max], ph1_2[max];
  double th2_2[max], ph2_2[max];
  Int_t runnum, runnum_2;
  double hallap;
  double hallap_2;
  //double r_s2_t_pads[max];
  //double l_s2_t_pads[max];
  //double r_s2_nthit;
  //double l_s2_nthit;
  double r_th_fp[max];
  double l_th_fp[max];
  double r_ph_fp[max];
  double l_ph_fp[max];
  double l_x_fp[max];
  double r_x_fp[max];
  double l_y_fp[max];
  double r_y_fp[max];
  double r_th_fp_2[max];
  double l_th_fp_2[max];
  double r_ph_fp_2[max];
  double l_ph_fp_2[max];
  double l_x_fp_2[max];
  double r_x_fp_2[max];
  double l_y_fp_2[max];
  double r_y_fp_2[max];
  //const int n = 16;
  //double r_s2_la_c[n];
  //double r_s2_ra_c[n];
  //double l_s2_la_c[n];
  //double l_s2_ra_c[n];
  //double rbeta[max];
  //double lbeta[max];
  //double nhit, nhit_R;
  double ps_asum;
  double ps_asum_2;
  //double a1_tdc[24];
  //double a2_tdc[26];
  double ctime[max];
  double ctime_2[max];
  double vz_mean[max];
  double vz_mean_2[max];
  double l_cer;
  double l_cer_2;
  

  t1->SetBranchAddress("fEvtHdr.fRun", &runnum    );
  t1->SetBranchAddress("HALLA_p", &hallap );
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);
  t1->SetBranchAddress("L.cer.asum_c", &l_cer);
  t1->SetBranchAddress("vz_mean", &vz_mean);
  t1->SetBranchAddress("R.tr.p", &mom1);
  t1->SetBranchAddress("L.tr.p", &mom2);
  t1->SetBranchAddress("R.tr.vz", &rvz);
  t1->SetBranchAddress("L.tr.vz", &lvz);
  t1->SetBranchAddress("R.tr.tg_th", &th1);
  t1->SetBranchAddress("R.tr.tg_ph", &ph1);
  t1->SetBranchAddress("L.tr.tg_th", &th2);
  t1->SetBranchAddress("L.tr.tg_ph", &ph2);
  t1->SetBranchAddress("R.tr.x",   &r_x_fp);
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("R.tr.y",   &r_y_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("R.tr.th",  &r_th_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("R.tr.ph",  &r_ph_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);
  t1->SetBranchAddress("R.ps.asum_c", &ps_asum);
  t1->SetBranchAddress("ctime", &ctime);
  
  
  TFile* fnew = new TFile("check_tune.root","recreate");
  //TTree* tnew = new TTree("tree","check for LHRS SS");
  //tnew->Branch("Xt",   &Xt,   "Xt/D");
  //tnew->Branch("Yt",   &Yt,   "Yt/D");
  //tnew->Branch("Zt",   &Zt,   "Zt[100]/D");
  //tnew->Branch("Xpt",  &Xpt,  "Xpt[100]/D");
  //tnew->Branch("Ypt",  &Ypt,  "Ypt[100]/D");
  //tnew->Branch("XFP",  &XFP,  "XFP[100]/D");
  //tnew->Branch("YFP",  &YFP,  "YFP[100]/D");
  //tnew->Branch("XpFP", &XpFP, "XpFP[100]/D");
  //tnew->Branch("YpFP", &YpFP, "YpFP[100]/D");
  

  char name_MctimecorL[500];
  char name_MctimecorR[500];
  sprintf(name_MctimecorL,"iniL.dat"); 
  sprintf(name_MctimecorR,"iniR.dat"); 
  ifstream MctimecorL(name_MctimecorL);
  ifstream MctimecorR(name_MctimecorR);
  for (int i=0 ; i<nParamT/2. ; i++){
    double par=0.;
    int p=0;
    MctimecorL >> par >> p >> p >> p >> p >> p;
    OptPar[i] = par;
    MctimecorR >> par >> p >> p >> p >> p >> p;
    OptPar[i+56] = par;
  }
  MctimecorL.close();
  MctimecorR.close();

  ntune_event = 0;
  
  TH1F* h1 = new TH1F("h1","",200,-10.0,10.0);
  h1->GetXaxis()->SetTitle("ctime (ns)");
  //h1->GetXaxis()->CenterTitle();
  //h1->GetXaxis()->SetRangeUser(-0.2,0.2);
  //h1->GetYaxis()->SetTitle("x");
  TH1F* h1_select = (TH1F*)h1->Clone("h1_select");
  
  TH1F* h2 = new TH1F("h2","",200,-50.,150.0);
  TH1F* h2_L = (TH1F*)h2->Clone("h2_L");
  TH1F* h2_S = (TH1F*)h2->Clone("h2_S");
  
  TH1F* h1_2 = new TH1F("h1_2","",200,-10.0,10.0);
  h1_2->GetXaxis()->SetTitle("ctime (ns)");
  TH1F* h2_2 = new TH1F("h2_2","",200,-50.,150.0);
  TH1F* h2_L_2 = (TH1F*)h2_2->Clone("h2_L_2");
  //TH1F* h2_S_2 = (TH1F*)h2_2->Clone("h2_S_h22");
  
  char tempc[500];
  char tempc2[500];
  
  for(int i=0 ; i<nmax ; i++){
    x[i]    = -2222.0; // x at FP
    y[i]    = -2222.0; // y at FP
    xp[i]   = -2222.0; // x' at FP
    yp[i]   = -2222.0; // y' at FP
    
    //tgang_xp[nmax] = -2222.0;
    //tgang_yp[nmax] = -2222.0;
    
    x2[i]    = -2222.0; // x at FP
    y2[i]    = -2222.0; // y at FP
    xp2[i]   = -2222.0; // x' at FP
    yp2[i]   = -2222.0; // y' at FP
    //tgang_xp2[nmax] = -2222.0;
    //tgang_yp2[nmax] = -2222.0;

    avz[i]   = -2222.0;
    
    //peak_flag[i] = -1;
    //tune_id[i]   = -1;
    
    //z_recon[ntune_event] = -2222.0;
    //holegroup[ntune_event] = -1;
  }
  
  for (int i=0 ; i< ent ; i++){
    // ----- Initialization ------- //
    for(int j=0 ; j<max ; j++){
      mom1[j]   = -2222.0;
      mom2[j]   = -2222.0;
      th1[j]   = -2222.0;
      ph1[j]   = -2222.0;
      th2[j]   = -2222.0;
      ph2[j]   = -2222.0;
      ctime[j] = -2222.0;
      vz_mean[j] = -2222.0;
      rvz[j]   = -2222.0;
      lvz[j]   = -2222.0;
    }
    a1 = -2222.0;
    a2 = -2222.0;
    ps_asum = -2222.0;
    l_cer = -2222.0;
    
    t1->GetEntry(i);
    
    //if(i+evshift<ent) t1->GetEntry(i+evshift); 
    //else t1->GetEntry(i-ent+evshift);
    double kcenter = 0.0;
    ctime[0] = ctime[0] - kcenter;

    if(fabs(rvz[0]-lvz[0])<0.025
       && fabs(vz_mean[0])<0.1
       && a1 < 2.0
       && a2 > 2.0
       && a2 < 18.0
       && l_cer>1500.0
       && fabs(ctime[0])<5.0 ){
      
      //h1->Fill(ctime[0]);
      double cor1, cor2;
      cor1 = calcf2t_3rd(OptPar, 
			 (l_x_fp[0]-XFPm)/XFPr, 
			 (l_th_fp[0]-XpFPm)/XpFPr, 
			 (l_y_fp[0]-YFPm)/YFPr, 
			 (l_ph_fp[0]-YpFPm)/YpFPr, 
			 (vz_mean[0]-Ztm)/Ztr,
			 1);
      cor2 = calcf2t_3rd(OptPar, 
			 (r_x_fp[0]-XFPm)/XFPr, 
			 (r_th_fp[0]-XpFPm)/XpFPr, 
			 (r_y_fp[0]-YFPm)/YFPr, 
			 (r_ph_fp[0]-YpFPm)/YpFPr, 
			 (vz_mean[0]-Ztm)/Ztr,
			 2);
      double ctime_wcor = ctime[0] + cor1 + cor2;
      
      h1->Fill(ctime_wcor);

      if(pcent-selection_width < ctime_wcor
	 && ctime_wcor < pcent+selection_width){
	
	if( ntune_event<nmax ){
	  x[ntune_event]  = (l_x_fp[0]-XFPm)/XFPr;    // scaled 
	  y[ntune_event]  = (l_y_fp[0]-YFPm)/YFPr;    // scaled
	  xp[ntune_event] = (l_th_fp[0]-XpFPm)/XpFPr; // scaled
	  yp[ntune_event] = (l_ph_fp[0]-YpFPm)/YpFPr; // scaled
	  x2[ntune_event]  = (r_x_fp[0]-XFPm)/XFPr;   // scaled 
	  y2[ntune_event]  = (r_y_fp[0]-YFPm)/YFPr;   // scaled
	  xp2[ntune_event] = (r_th_fp[0]-XpFPm)/XpFPr;// scaled
	  yp2[ntune_event] = (r_ph_fp[0]-YpFPm)/YpFPr;// scaled
	  avz[ntune_event] = (vz_mean[0]-Ztm)/Ztr;    // scaled
	  
	  coincoin[ntune_event] = ctime[0];
	  h1_select->Fill(ctime[0]);
	  
	  ntune_event++;
	}// if (ntune_event)
      }// if (ctime +/- selection_width)
    } // if (kaon selection etc.)
  }// for (ent)
  

  cout << " The number of events selected to be used for tuning: "
       << ntune_event << endl;
  cout << endl;

  
  double chi_sq[nite];
  double x[nite];
  if (nite>0) cout << " Tuning started: " << endl;
  
  for(int i=0 ; i<nite ; i++){
    // --------------------------- //
    // ---- Parameter tuning ----- //
    // --------------------------- //
    x[i] = i+1;

    chi_sq[i] = tune(OptPar,i); 
    
    cout << " Tuning# = " << i << ": chisq = "
	 << chi_sq[i] << endl;
    cout << endl;
    
    sprintf(tempc,  "./newpar/newpar_ctimeR_%d.dat",i); 
    sprintf(tempc2, "./newpar/newpar_ctimeL_%d.dat",i); 

    ofstream * ofs1 = new ofstream(tempc);
    ofstream * ofs2 = new ofstream(tempc2);
    int nppp = 0;
    const int nn = 3; // 3rd order matrix using xf, xpf, y, ypf, and zt
    for(int i=0 ; i<nn+1 ; i++){
      for(int e=0 ; e<nn+1 ; e++){
	for(int d=0 ; d<nn+1 ; d++){
	  for(int c=0 ; c<nn+1 ; c++){
	    for(int b=0 ; b<nn+1 ; b++){
	      for(int a=0 ; a<nn+1 ; a++){  
		if(a+b+c+d+e==i){
		  *ofs1 << OptPar[nppp+56] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  *ofs2 << OptPar[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  nppp++;
		  //	      cout << Plen_opt[nppp] 
		  //		   << " " << a 
		  //		   << " " << b
		  //		   << " " << c
		  //		   << " " << d << endl;
		}
	      }
	    }
	  }
	}
      }
    }
    ofs1->close();
    ofs1->clear();
    ofs2->close();
    ofs2->clear();

  }
  
  // =================================== //
  // ======== Draw histograms ========== //
  // =================================== //
  TCanvas* c1 = new TCanvas("c1","c1");
  h1->Draw();
  h1_select->SetLineColor(2);
  h1_select->Draw("same");

  //TCanvas* c2 = new TCanvas("c2","c2");
  //h2->Draw();
  //h2_L->SetLineColor(2);
  //h2_S->SetLineColor(9);
  //h2_L->Draw("same");
  //h2_S->Draw("same");
  //
  //TCanvas* c3 = new TCanvas("c3","c3");
  //h1_2->Draw();
  //TCanvas* c4 = new TCanvas("c4","c4");
  //h2_2->Draw();
  //h2_L_2->SetLineColor(2);
  //h2_L_2->Draw("same");
  
  TObjArray h(1);
  h.Add(c1);
  //h.Add(c2);
  h.Add(h1);
  h.Add(h1_select);
  //h.Add(h2);
  //h.Add(h2_L);
  //h.Add(h2_S);
  
  //h.Add(h1_2); 
  //h.Add(h2_2);
  //h.Add(h2_L_2);
  
  //for(int i=0 ; i<npeak ; i++){}
  
  if(nite>0){
    TGraph * gr1 = new TGraph(nite,x,chi_sq);
    gr1->SetName("gr1");
    TCanvas * c_chi2 = new TCanvas("c_chi2","c_chi2",600,600);
    c_chi2->cd(1); gr1->Draw("*la");
    h.Add(gr1);
    h.Add(c_chi2);
  }
  h.Write();

  app->Run();
  
  return 0;
} 

//////////////////////////////////////////////////
double calcf2t_zt(double* P, double xf, double xpf, 
		  double yf, double ypf)
//////////////////////////////////////////////////
{
  // -----3rd order ----- 
  // This is the third order claculation byb  using 35 parameter
  const int nMatT=3;  // These are for the RHRS use same for the LHRS
  const int nXf=3;
  const int nXpf=3;
  const int nYf=3;
  const int nYpf=3;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0;
  
  for (int n=0;n<nMatT+1;n++){ 
    for (d=0;d<n+1;d++){
      for (c=0;c<n+1;c++){ 
	for (b=0;b<n+1;b++){
	  for (a=0;a<n+1;a++){ 
	    
	    if (a+b+c+d==n){
	      if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf){
		x = pow(xf,double(a))*pow(xpf,double(b))*
		  pow(yf,double(c))*pow(ypf,double(d));
	      }
	      else{
		x = 0.;
	      }
	      Y += x*P[npar]; 
	      npar++;
	    }
	    
	  }
	}
      }    
    }
  }
  
  return Y; 
  
}

//////////////////////////////////////////////////
double calcf2t_3rd(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt,
		   int flag=0)
//////////////////////////////////////////////////
{
  // ------------------------------------------------ //
  // ----- 3rd order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //
  const int nMatT=3;  
  const int nXf=3;
  const int nXpf=3;
  const int nYf=3;
  const int nYpf=3;
  const int nZt=3;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		if(flag==2) Y += x*P[npar+56]; // RHRS
		else Y += x*P[npar];    // LHRS
	      npar++;
	      }
	      
	    }
	  }
	}
      }    
    }
  }
  
  return Y; 
  
}

//////////////////////////////////////////////////
double calcf2t_4th(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt,
		   int flag=0)
//////////////////////////////////////////////////
{
  // ------------------------------------------------ //
  // ----- 4rd order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=4;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		if(flag==2) Y += x*P[npar+126]; // LHRS
		else Y += x*P[npar];    // RHRS
	      npar++;
	      }
	      
	    }
	  }
	}
      }    
    }
  }
  
  return Y; 
  
}




//////////////////////////////////////////////////
double calcf2t_5th(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt,
		   int flag=0)
//////////////////////////////////////////////////
{
  // ------------------------------------------------ //
  // ----- 5th order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //
  const int nMatT=5;  
  const int nXf=5;
  const int nXpf=5;
  const int nYf=5;
  const int nYpf=5;
  const int nZt=5;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		if(flag==2) Y += x*P[npar+252]; // LHRS
		else Y += x*P[npar];    // RHRS
	      npar++;
	      }
	      
	    }
	  }
	}
      }    
    }
  }
  
  return Y; 
  
}


// #############################################################
double tune(double* pa, int j) 
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10]; 
  int ierflg = 0;
  int allparam = nParamT;
  //cout << allparam << endl;
  TMinuit* minuit = new TMinuit(allparam);
  minuit->SetFCN(fcn);

  double start[allparam];
  double step[allparam];
  const int nMatT=3;  
  const int nXf=3;
  const int nXpf=3;
  const int nYf=3;
  const int nYpf=3;
  const int nZt=3; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  start[npar]    = pa[npar];
		  step[npar]     = pa[npar] * 5.0e-3; 
		  start[npar+56] = pa[npar+56];
		  step[npar+56]  = pa[npar+56] * 5.0e-3;  
		  
		  if(step[npar]==0) step[npar]=0.05;
		  if(step[npar+56]==0) step[npar+56]=0.05;
		}
		else{
		  // ----- No tuning ------ //
		  start[npar]     = 0.0;
		  step[npar]      = 0.0;
		  start[npar+56] = 0.0;
		  step[npar+56]  = 0.0;
		}
		npar++;
	      }
	    }
	  }
	}
      }    
    }
  }
  
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1);
  
  double LLim[allparam];// Lower limit for all of the parameter
  double ULim[allparam];// Upper limit for all of the parameter
  char pname[500];
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    //start[i] = pa[i]; 
    //step[i] = 1.0e-3;
    
    //LLim[i] = pa[i] - pa[i]*0.8;
    //ULim[i] = pa[i] + pa[i]*0.8;
    LLim[i] = pa[i] - 5.0; // temp
    ULim[i] = pa[i] + 5.0; // temp
    
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }
  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0; // original
  //arglist[0] = 1.0; // test
  //arglist[0] = 0.0;   // test
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  
  // ~~~~ Migrad + Simplex  ~~~~ 
  arglist[0] = 50000;
  //arglist[1] = 0.01; // origial
  arglist[1] = 0.02; // 
  //arglist[1] = 0.1; // 
  minuit -> mnexcm("MINImize",arglist,2,ierflg); // Chi-square minimization
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double er;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);

  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){
    minuit -> GetParameter(i,OptPar[i],er);
  }
  
  return chi2;
}

// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  
  const double sigma = 1.0; // MeV/c2
  double residual = 0.0;
  double chi2 = 0.0;
  double ctime_wcor;
  double cor1, cor2;
  
  chi2 = 0.0;
  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;
    cor1 = calcf2t_3rd(param,
		       x[i], xp[i], y[i], yp[i],
		       avz[i], 1);
    cor2 = calcf2t_3rd(param,
		       x2[i], xp2[i], y2[i], yp2[i],
		       avz[i], 2);
    ctime_wcor = coincoin[i] + cor1 + cor2;
    
    residual = ctime_wcor - pcent_real;
    chi2 = chi2 + pow(residual/sigma, 2.0);
  }
  
  fval = chi2;
    
}

double CalcMM(double ee, double* pvec_ep, double* pvec_k, double mt){
  
  double pe = ee;
  double Ee = sqrt(me*me + pe*pe);
  TVector3 vec_e (0.0, 0.0, pe);
  
  double pep  = pvec_ep[0];
  double xpep = pvec_ep[1];
  double ypep = pvec_ep[2];
  double px_ep, py_ep, pz_ep;
  pz_ep = pep / sqrt(1.0 + xpep*xpep + ypep*ypep);
  px_ep = xpep * pz_ep;
  py_ep = ypep * pz_ep;
  TVector3 vec_ep (px_ep, py_ep, pz_ep);
  vec_ep.RotateY(hrs_ang);
  //double Eep = sqrt(vec_ep * vec_ep);
  double Eep = sqrt(pep*pep + me*me);
  
  double pk  = pvec_k[0];
  double xpk = pvec_k[1];
  double ypk = pvec_k[2];
  double px_k, py_k, pz_k;
  pz_k = pk / sqrt(1.0 + xpk*xpk + ypk*ypk);
  px_k = xpk * pz_k;
  py_k = ypk * pz_k;
  TVector3 vec_k (px_k, py_k, pz_k);
  vec_k.RotateY(-hrs_ang);
  //double Ek = sqrt(vec_k * vec_k);
  double Ek = sqrt(pk*pk + mk*mk);
  
  double missingE2 = 0.0, missingP2 = 0.0, missingM2 = 0.0;
  missingE2 = pow(Ee + mt - Ek - Eep, 2.0);
  missingP2 = (vec_e - vec_ep - vec_k) * (vec_e - vec_ep - vec_k);
  missingM2 = missingE2 - missingP2;
  
  double MissingMass = 0.0;
  MissingMass = sqrt(missingM2);

  return MissingMass;
  
}





