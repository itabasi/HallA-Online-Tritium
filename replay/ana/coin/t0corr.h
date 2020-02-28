#ifndef t0corr_h
#define t0corr_h 1
using namespace std;
#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"
#include "Tree.h"
#include "ParamMan.h"
#include "define.h"
#include "coincalib.h"


const  int ns2seg=16;
const  int ncanvas=19;
class t0corr : public Tree
{

 public:
  t0corr();
  ~t0corr();

  void SetRoot(string ifname);
  void SetRunList(string ifname);
  void NewRoot(string ofname);
  void ReadParam(string name);
  void MakeHist();
  void GetParam(string ofname);
  void CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit);
  void Fill();
  void Fit();
  void Print(string ofname);
  void Draw();
  void PathCalc(int rhit, int lhit);
  void Write();

  void Close();
  
  private:
  Setting *set;
  ParamMan *param;


 public:


  ofstream * ofp;
  TFile* ofr;
  TTree* tnew;


  TH1D* hRs2t[ns2seg];
  TH1D* hRs2b[ns2seg];
  TH1D* hLs2t[ns2seg];
  TH1D* hLs2b[ns2seg];
  TH1D* hcoin;
  TH2D* hcoin_r;
  TH2D* hcoin_l;
  TH2D* hcoin_rk;
  TH2D* hcoin_lk;  
  TH1D* hcoin_k;
  TF1* fRs2t[ns2seg];
  TF1* fRs2b[ns2seg];
  TF1* fLs2t[ns2seg];
  TF1* fLs2b[ns2seg];
  TF1* fcoin;
  TF1* fcoin_k;
  TGraphErrors* gRct;
  TGraphErrors* gLct;

  TCanvas* c[ncanvas];


  bool z_flag, pid_flag;

  double min_ct, max_ct;
  int bin_ct;

  double min_pi, max_pi;
  int bin_pi;

  double min_k, max_k;
  int bin_k;
 
  double min_fit[ns2seg], max_fit[ns2seg];
  int bin_fit[ns2seg];

  double Rs2t_posi[ns2seg], Rs2b_posi[ns2seg], Ls2t_posi[ns2seg], Ls2b_posi[ns2seg];
  double Rs2t_posi_ch[ns2seg], Rs2b_posi_ch[ns2seg], Ls2t_posi_ch[ns2seg], Ls2b_posi_ch[ns2seg];
  double ref_ch;
  double mean_pi, mean_k;

  double Rtof, Rtof_c,Ltof,Ltof_c,R_pathl, L_pathl,Rp, tdc_time,Rs2_t,Ls2_t,Beta_R,Beta_L,R_pathl_c, L_pathl_c,Beta_K,Rtof_K,coin_k;
  int R_s2_pad,L_s2_pad; 
  double LS2T_off[ns2seg],LS2B_off[ns2seg],RS2T_off[ns2seg],RS2B_off[ns2seg];
  double Offset[4][ns2seg];
  int z_cut,pid_cut;
  int runnum;
 



};

t0corr::t0corr(){
  
  set = new Setting();



  


}

t0corr::~t0corr(){}


#endif
