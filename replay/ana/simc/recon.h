#ifndef recon_h
#define recon_h 1
#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"
#include "define.h"

class recon{

public:
  recon();
  ~recon();
  
  void SetRoot(string ifname);
  void SetBranch();
  void Calib();
  void Fill();
  void MTParam_R();
  void MTParam_L();
  void MTP_mom();  
  void SetMatrix(string mtparam);
  void NewRoot(string ofnameo);
  double Eloss(double yp, double z, char* arm);
  void Write();
  TFile* ifp;
  TChain* T;
  //------ NewROOT -----------//
  TFile* ofp;
  TTree* tnew;


  float R_tr_x,R_tr_th,R_tr_y,R_tr_ph;
  float L_tr_x,L_tr_th,L_tr_y,L_tr_ph;
  float R_tr_tg_th,R_tr_tg_ph;
  float L_tr_tg_th,L_tr_tg_ph;
  float R_Ras_x, L_Ras_x;
  float mm;
  float zposi,R_tr_vz,L_tr_vz;  

  float R_p,Rp_x,Rp_y,Rp_z;
  float L_p,Lp_x,Lp_y,Lp_z;
  float Ee,Ee_,Ek;
  double B_p;

  double mm_c;
  int ENum;
  string param[10];
};

  recon::recon(){T=new TChain("SNT");};
  recon::~recon(){};


  
#endif
