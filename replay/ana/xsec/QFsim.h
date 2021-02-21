#ifndef QFsim_h
#define QFsim_h 1
#include <TLorentzVector.h>
#include "define.h"
#include "Param.h"
#include <complex.h>
const int nmax =200;
const int nmax2 =200;

class QFsim{

public:
  QFsim();
  ~QFsim();
  void SetRootList(string ifrname);
  void SetHist();
  void FitQFL();
  void FitQFnnL();
  void FitQFnnL_new();
  void Write();
  void CalcPS(TH1D* hexp);
  void CalcPS_c(TH1D* hexp, TH1D* hsim, int model);
  
  TFile* fH_sim;
  TFile* fT_sim;
  TFile* fexp;
  TFile* fT_sim_10keV;
  TFile* fH_sim_10keV;
  TFile* fexp_10keV;
  TFile* ofr;
  TTree* Tnew;
  
  int bin_mm, xbin_size;
  double min_mm,max_mm;
  double NL_bg;

  TGraphErrors* gchi_L;
  TGraphErrors* gchi_nnL;
  TGraphErrors* gchi_fsi[100];
  TGraphErrors* gdiff_L;
  TGraphErrors* gdiff_fsi[100];
  
  //=== CalcPS =====//
  TGraph* g_s[100];
  TGraph* g_n[100];
  TGraph* g_val[100];
  TGraph* g_ps[100];  

  TGraph* g_s_orig;
  TGraph* g_n_orig;
  TGraph* g_val_orig;
  TGraph* g_ps_orig;    

  
  TGraph* g_ps_sig[100];
  TGraph* g_val_sig[100];
  TGraph* g_mm_sig[100];

  TH2D* hPS[100];
  TH2D* hPval[100];

  TH2D* hPS_orig;
  TH2D* hPval_orig;
  
  // --- CalcPS Paramters----//
  int Nsig=10;
  double min_sigma = 1.0;
  double max_sigma = 2.0;
  
  const int nvp=4;
  TH1D* hexp;
  TH1D* hexp_c;
  TH1D* hexp_nnL;
  TH1D* hexp_peak_L;
  TH1D* hexp_L;
  TH1D* hexp_peak;
  TH1D* hexp_acc;
  TH1D* hexp_10keV;
  TH1D* hexp_10keV_s;
  TH1D* hexp_10keV_c;
  TH1D* hexp_acc_10keV;
  TH1D* hmm_10keV;
  TH1D* hmm_10keV_s;
  TH1D* hmm_10keV_c;

  TH1D* hH_nnLsim_10keV;
  TH1D* hH_nnLsim_10keV_s;  
  TH1D* hmm_fsi1_10keV;
  TH1D* hmm_fsi2_10keV;
  TH1D* hmm_fsi3_10keV;
  TH1D* hmm_fsi1_10keV_s;
  TH1D* hmm_fsi2_10keV_s;
  TH1D* hmm_fsi3_10keV_s;  
  TH1D* hmm_fsi1_10keV_c;
  TH1D* hmm_fsi2_10keV_c;
  TH1D* hmm_fsi3_10keV_c;



  
  TH1D* hT_Lsim;
  TH1D* hH_Lsim;
  TH1D* hH_nnLsim;
  TH1D* hH_nnLsim_s;
  TH1D* hT_nnLsim;
  TH1D* hT_nnLsim_s;  
  TH1D* hT_Lsim_s;
  TH1D* hH_Lsim_s;
  TH1D* hnnL_sim;
  TH1D* hnnL_sim_s;

  TH1D* hmm;
  TH1D* hmm_fsi1;
  TH1D* hmm_fsi2;
  TH1D* hmm_fsi3;
  TH1D* hmm_s;
  TH1D* hmm_fsi1_s;
  TH1D* hmm_fsi2_s;
  TH1D* hmm_fsi3_s;
  TH1D* hmm_fsi1_c;
  TH1D* hmm_fsi2_c;
  TH1D* hmm_fsi3_c;  
  TH1D* hmm_c;
  double  wmin_L;
  double  wmin_nnL;
};


QFsim::QFsim(){};
QFsim::~QFsim(){};

#endif
