#ifndef coincalib_h
#define coincalib_h 1
using namespace std;
#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"
#include "Tree.h"
#include "ParamMan.h"
#include "define.h"

  bool single=false;
  int nite=0;
  int MODE=0;
class coincalib : public Tree
{
  public:
  coincalib();
  ~coincalib();
  
  private:
  Setting *set;
  ParamMan *param;
  public:
  //  double CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit);
  //  double CoinCalc_c(int RS2_seg, int LS2_seg, int rhit, int lhit);
  void CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit);
  void PathCalc(int rhit, int lhit);
  void SetRoot(string ifname);
  void SetRunList(string ifname);
  void NewRoot(string ofname);
  void ReadParam(string name);
  void SetMT(bool rhrs);

  double SetOff();
  void CoinTuning(string ofname, int MODE);
  void MTParam(string mtparam, int MODE);
  void MTParam_R(string mtparam);
  void MTParam_L(string mtparam);
  void EventSelect();
  double tune(double*pa, int MODE);
  void Fill();
  void Write();
  void Close();

  //=== Parameters ===//
  TFile* ofr;
  TTree* tnew;
  double coin_offset;
  double coint,coint_c;
  double Rp,Rp_c;
  bool z_flag;
  bool pid_flag;
  bool ct_flag;
  bool kaon_flag;
  bool pion_flag;
  bool fp_flag;
  double tdc_time;
  double Rtof,Ltof;
  int R_s2pad,L_s2pad;
  double chi_sq1[100],chi_sq2[100],chi_sq[100];
  TH1D* hcoin_cut;
  double Rtof_c,Ltof_c;
  double Rs2_t,Ls2_t;
  double beta_r,beta_l;
  double R_pathl,L_pathl;
  double R_pathl_c,L_pathl_c;
  double Beta_R,Beta_L;
  
};
  
  coincalib::coincalib(){
  set = new Setting();
  mode=MODE;
  cout<<"Analysis mode "<<mode<<endl;
  };

coincalib::~coincalib(){};




#endif
