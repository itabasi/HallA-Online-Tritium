#ifndef t0calib_h
#define t0calib_h 1
using namespace std;
#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"
#include "Tree.h"
#include "ParamMan.h"
#include "define.h"


class t0calib: public Tree{

public:
  t0calib();
  ~t0calib();
private:
  Setting *set;
  ParamMan *param;
  coincalib *coinc;
public:
  void NewRoot(string ofname);
  void MakeHist();
  void Fill();
  
  //===== Histoguram ======//
  const ns2seg=16;
  TH1D* hcoinRs2T[ns2seg];
  TH1D* hcoinRs2B[ns2seg];
  TH1D* hcoinLs2T[ns2seg];
  TH1D* hcoinLs2T[ns2seg];
  TH1D* hs0Rs2T[ns2seg];
  TH1D* hs0Rs2B[ns2seg];
  TH1D* hs0Ls2T[ns2seg];
  TH1D* hs0Ls2T[ns2seg];
  double min_ct,max_ct;
  int bin_ct;
  double min_tof,max_tof;
  int bin_tof;
  
};

t0calib::t0calib(){
  set = new Setting();
  set->Initialize();
  coinc = new coincalib(); 
};

t0calib::~t0calib(){};
