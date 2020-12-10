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
  
  TFile* fH_sim;
  TFile* fT_sim;
  TFile* fexp;
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
  const int nvp=4;
  TH1D* hexp;
  TH1D* hexp_c;
  TH1D* hexp_nnL;
  TH1D* hexp_peak_L;
  TH1D* hexp_L;
  TH1D* hexp_peak;
  TH1D* hexp_acc;
  
  TH1F* hT_Lsim;
  TH1F* hH_Lsim;
  TH1F* hH_nnLsim;
  TH1F* hH_nnLsim_s;
  TH1F* hT_nnLsim;
  TH1F* hT_nnLsim_s;  
  TH1F* hT_Lsim_s;
  TH1F* hH_Lsim_s;
  TH1F* hnnL_sim;
  TH1F* hnnL_sim_s;

  TH1F* hmm;
  TH1F* hmm_fsi1;
  TH1F* hmm_fsi2;
  TH1F* hmm_fsi3;
  TH1F* hmm_s;
  TH1F* hmm_fsi1_s;
  TH1F* hmm_fsi2_s;
  TH1F* hmm_fsi3_s;
  TH1F* hmm_fsi1_c;
  TH1F* hmm_fsi2_c;
  TH1F* hmm_fsi3_c;  
  TH1F* hmm_c;
  double  wmin_L;
  double  wmin_nnL;
};


QFsim::QFsim(){};
QFsim::~QFsim(){};

#endif
