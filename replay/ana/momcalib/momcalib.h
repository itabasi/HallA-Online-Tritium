#ifndef momalib_h
#define momcalib_h 1
#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"
#include "define.h"
#include "tree.h"
//using namespace std;





//=================================================//
//============= MomCalib Class ====================//
//=================================================//

class momcalib : public tree
{

 public:
  momcalib();
  ~momcalib();

 public:
  //  void SetBranch(string ifname, bool rarm);
  //  void NewBranch(string ofname, bool rarm);
  void EventSelection(double ww);
  void SetRoot(string ifname);
  void SingleRoot(string ifname);
  void MakeHist();
  void SetInput(string infile);
  string SetVal(string line, string name);
  void NewRoot(string ofname);
  void MTParam(string mtparam);
  void MTParam_R();
  void MTParam_L();
  void MTP_mom();
  void MTP_weight();
  void ParamCorr();
  void SetAlEvents(double weight, double mean,double width);
  void CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit);
  void zCorr(bool rarm, bool larm);
  void rasCorr(bool rarm, bool larm );
  void angCorr(bool rth, bool rph, bool lth, bool lph);
  void plossCorr(bool PLoss);
  void momCorr(bool rarm, bool larm);
  void momCalc();
  void nmatrix(int n);
  double AC_npe(int nac,int seg,double adc);
  void GetACParam();
  double SSx(double z, double th);
  double SSy(double z, double ph);  
  int mode(string arm,char* Target, int F1tdc);
  void MomTuning(string ofname);
  double Eloss(double yp,double z,  char* arm);
  double tune(double* pa, int j, int angflag);
  double tune_Tkine(double* pa, int j, int angflag);   
  double* MixiedEvent(double * coin_acc);
  void MixedEvent();
  void Fill();
  void Close();


  Setting *set = new Setting();
  
  
  //------mode -------------//
  char* target;
  int tdc_mode=0;
  
  //------ MTParam --------//
  string param[100];
  bool MT[10];  
  bool ploss;
  bool ploss_f;
  //------ NewROOT -----------//
  TFile* ofp;
  TTree* tnew;
  double coin_time;
  double Dpe,Dpe_,Dpk;
  double RP,LP;
  //----- SSx && SSy ------//
  const double l0 = 100.3; // length target to Sieve slit
  double rssx,rssy,lssx,lssy,rssx_c,rssy_c,lssx_c,lssy_c;

  //----- MakeHist ---------//  
  //      Fill    //
  TH1F* hLz;
  TH1F* hLth;
  TH1F* hLph;
  TH1F* hLz_c;
  TH1F* hLz_rc;
  TH1F* hLth_c;
  TH1F* hLph_c;
  TH1F* hRp;
  TH1F* hRp_c;
  TH1F* hRp_ec;  
  TH1F* hLp;
  TH1F* hLp_c;
  TH1F* hLp_ec;
  TH2F* hLss;
  TH2F* hLss_c;
  TH2F* hRss;
  TH2F* hRss_c;
  TH1F* hz_Al;
  //   EventSelection //
  TH1F* hRz_es;
  TH1F* hRth_es;
  TH1F* hRph_es;
  TH1F* hRz_es_c;
  TH1F* hRz_es_rc;
  TH1F* hRth_es_c;
  TH1F* hRph_es_c;
  TH1F* hLz_es;
  TH1F* hLth_es;
  TH1F* hLph_es;
  TH1F* hLz_es_c;
  TH1F* hLz_es_rc;
  TH1F* hLth_es_c;
  TH1F* hLph_es_c;
  TH1F* hmm_Al;
  TH2F* hmm_Al_nnL;
  TH2F* hmm_Al_nnL_zp;
  TH2F* hmm_Al_nnL_zn;
  TH1F* hmm_Al_zp;
  TH1F* hmm_Al_zn;
  TH1F* hmm_Al_zp_c;
  TH1F* hmm_Al_zn_c;  
  TH1F* hmm_Al_mom;
  TH1F* hmm_Al_acc;
  TH1F* hmm_Al_select;
  TH1F* hmm_Al_cut;
  TH1F* hmm_Al_cut_acc;
  TH1F* hmm_Al_select_c;
  TH1F* hRp_es;
  TH1F* hRp_es_c;
  TH1F* hRp_es_ec;  
  TH1F* hLp_es;
  TH1F* hLp_es_c;
  TH1F* hLp_es_ec;
  TH2F* hLss_es;
  TH2F* hLss_es_c;
  TH2F* hRss_es;
  TH2F* hRss_es_c;    


  TH1F* hct;
  TH1F* hct2; 
  TH1F* hLz_cut;
  TH1F* hRz_cut;  
  TH1F* hct_cut;
  TH1F* hac1_cut;
  TH1F* hac2_cut;
  TH1F* hcer_cut;

  double dpe,dpe_,dpk;
  //----- tuning check -----//
  TH1F* hmm_L_t[100];
  TH1F* hmm_S_t[100];
  TH1F* hmm_t[100];  
  TH1F* hmm_L_tc[100];
  TH1F* hmm_S_tc[100];
  TH1F* hmm_tc[100];  

  //------------------------//  
  TH1F* hRz;
  TH1F* hRth;
  TH1F* hRph;
  TH1F* hRz_c;
  TH1F* hRz_rc;  
  TH1F* hRth_c;
  TH1F* hRph_c;    
  TH1F* hmm_select;
  TH1F* hmm_L;
  TH1F* hmm_S;
  TH1F* hmm_nnL;
  TH1F* hmm_nnL_acc;
  TH1F* hmm_cut;

  TH1F* hmm_nnL_select;
  TH1F* hmm_nnL_cut;
  TH1F* hmm_nnL_cut_acc;

  
  
  TH1F* hmm_b;
  TH1F* hmm_cut_H;
  TH1F* hmm_cut_T;
  TH1F* hLam;
  TH1F* hSig;


  double coin_offset; // H1 mode
  
  double min_z,max_z;
  int bin_z;
  double min_th,max_th;
  int bin_th;
  double min_ph,max_ph;
  int bin_ph;
  double min_lp,max_lp;
  int bin_lp;
  double min_rp,max_rp;
  int bin_rp;
  double min_bp,max_bp;
  int bin_bp;
  double min_ssx,max_ssx,min_ssy,max_ssy;
  int bin_ssx,bin_ssy;
  double min_mm,max_mm;
  int bin_mm;

  double min_ct,max_ct;
  int bin_ct;
  double tdc_time=0.056;//[ns]
  
  //------ MTtuing -----------//


  //------ EventSelection ---//
  int TNum;
  div_t d;
  double chi_Al=0.0;
  double chi_L=0.0;
  double chi_S=0.0;
  double chi2_init=0.0;
  //--- Event Selection mass cut ---//
  //  double mmL_range = 0.005; // [GeV/c^2]
  //  double mmS_range = 0.005; // [GeV/c^2]
  double mmL_range = 0.003; // [GeV/c^2]
  double mmS_range = 0.003; // [GeV/c^2]
  //  double mmL_range = 0.002; // [GeV/c^2]
  //  double mmS_range = 0.002; // [GeV/c^2]
  
    int nlam=0;
    int nsig=0;
    
 // Cut Parameters //
 const double coin_cutmin=-248;
 const double coin_cutmax=-244; 
 const double rpathl_cutmin=28.7;
 const double rpathl_cutmax=29.4;
 const double lpathl_cutmin=28.6;
 const double lpathl_cutmax=29.2;
 const double rbeta_cutmin=0.0;
 const double rbeta_cutmax=1.0;
 const double lbeta_cutmin=0.9;
 const double lbeta_cutmax=1.0;
 const double Rvz_cutmin=-0.1;
 const double Rvz_cutmax= 0.1;
 const double Lvz_cutmin=-0.1;
 const double Lvz_cutmax= 0.1;
 const double Rx_cutmin= -0.4;
 const double Rx_cutmax= 0.4;
 // Event Flag //
  bool nnL_run;
  bool Rs0_flag;
  bool Ls0_flag;
  bool Rs2_flag;
  bool Ls2_flag;
  bool Scinti_flag;
  bool ac1_flag;
  bool ac2_flag;
  bool gs_flag;
  bool sh_flag;
  bool ps_flag;
  bool trig_flag;
  bool track_flag;
  bool Rpathl_flag;
  bool Lpathl_flag;
  bool coin_flag;
  bool x_flag;
  bool y_flag;
  bool z_flag;
  bool Al_flag;
  bool Lam_run;
  bool Sig_run;
  bool RPID_flag;
  bool LPID_flag;
  bool mom_select;
  int tuned_num;
  double mm;
  double mm_L,mm_L_b;
  double mm_L_c;
  double mm_L_pz;
  double mm_acc;
  double mm_pi;
  double MM_nnL,MM_Al,MM_Al_nnL;
  double MM_nnL_b,MM_Al_b;
  double MM_L_b,MM_L;
  double mm_Al,mm_Al_acc,mm_Al_c,mm_Al_acc_c;
  double mm_nnL,mm_nnL_acc,mm_nnL_c,mm_nnL_acc_c;
  double mm_Al_nnL, mm_X;
  double ct;
  double ctime;
  double Ee,Ee_,Ek,Epi;
  double Ee_2,Ek2;
  double pe,pe_,pk,ppi;
  double coin_t,coin_tc;
  double rbeta,rbeta_k,lbeta;
  double Rs2_off,Ls2_off; 
  double Rs2_tcorr,Ls2_tcorr;
  int Ls2pads,Rs2pads;
  bool cut_ac1,cut_ac2,cut_beta;
  int nac1,nac2,nac3,n;
  double tof_r,tof_l; 
  double rpathl,lpathl;
  double corr_R,corr_L;
  double rpath_corr,lpath_corr;
  int tflag;
  int tev;
  double Lp_z,Lp_x,Lp_y;
  double Rp_z,Rp_x,Rp_y;
  double LP_z,LP_x,LP_y;
  double RP_z,RP_x,RP_y;  
  bool cut, pid_cut,ct_cut,z_cut;
  
    //---- MomTuning ----//
    double chi_sq[1000],chi_sq1[100],chi_sq2[100];
    TGraphErrors* gchi_p  = new TGraphErrors();
    TGraphErrors* gchi_Rp = new TGraphErrors();
    TGraphErrors* gchi_Lp = new TGraphErrors();
    TGraphErrors* gchi_w  = new TGraphErrors();
    //----- ELoss ------//
    
    
    //----  ParamCorr -----//
    int nPC=0;
    bool PC_flag=false;
    bool Lp_scale=false; 
};

//=====================================//
//======== Momcalib ==================//
//=====================================//

momcalib::momcalib(){

 set->Initialize();
 
};
momcalib::~momcalib(){};

#endif
