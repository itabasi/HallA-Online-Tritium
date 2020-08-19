#ifndef momalib_h
#define momcalib_h 1
#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"
#include "define.h"
#include "tree.h"
//using namespace std;

bool RHRS;
bool LHRS;
bool single  = false;
bool MT_f[10];
bool Initial = false;
//bool Al_mode = true;
bool Al_mode = false;
//bool nnL_mode =true;
bool nnL_mode  = false;
bool form_mode = false;
bool Al_check  = false;
//bool form_mode = true;
//bool Initial=true;
const int MAX=1000;
int nite=0;
int MODE=5;
double weight=0.0;
double weightT=0.0;
double weightAl=0.0;
double weightAl_def=0.0;
double weight_def =0.0;
double weightnnL=0.0;
double range_MgL=0.003; 
double Al_range =0.0;
double MgL_weight;
double MgL_mean;
double MgL_width;
bool Goga=false;

extern double s2f1_off(int i,char* ARM,char* MODE, int KINE);
extern void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
extern void fcn_new(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
extern double tuning(double* pa, int j, int MODE); 
extern  double Calc_ras(double a,double b,double c){return  a *b + c;};  
extern double calcf2t_ang(double* P,double xf, double xpf, double yf, double fpf,double z);
extern double calcf2t_zt(double* P, double xf, double xpf, double yf, double ypf);
extern double calcf2t_mom(double* P, double xf, double xpf, double yf, double ypf, double zt);
extern double calcf2t_mom_RL( double* RP, double Rxf, double Rxpf, double Ryf, double Rypf, double Rzt,
			      double* LP, double Lxf, double Lxpf, double Lyf, double Lypf, double Lzt	     );
extern double expgaus(double *x, double *par);
extern double expgaus_mean(double *x, double *par);
extern double expgaus_sigma(double *x, double *par);
extern double Get_expgaus_mean(double * par);
extern double Get_expgaus_sigma(double * par, double mean);


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


void momcalib::nmatrix(int n){

  if(n==2){nnp=2; nParamTp=21;
  }else if(n==3){nnp=3; nParamTp=56;
  }else if(n==4){nnp=4; nParamTp=126;
  }else if(n==5){nnp=5; nParamTp=252;
  }else if(n==6){nnp=6; nParamTp=462;
  }else{cout<<"Matrix parameter's Errror "<<endl;exit(1);}


};



//////////////////////////////////////////////////////////////

void momcalib::SetAlEvents(double weight, double mean,double width){

  MgL_width = width;
  MgL_mean  = mean;
  MgL_weight  = weight;

  cout<<"=======< Al Tuning Parameters >=========="<<endl;
  cout<<" weight : "<<MgL_weight<<endl;
  cout<<" mean : "<<MgL_mean<<" MeV"<<endl;
  cout<<" width : "<<MgL_width<<" MeV"<<endl;

  // MeV to GeV
  MgL_width *=1./1000.;
  MgL_mean  *=1./1000.;
}

////////////////////////////////////////////////////////////

int momcalib::mode(string arm, char* Target, int F1tdc){
  cout<<endl;
  cout<<"======================================="<<endl;
  if(arm=="R"){
    MODE=-1;
    cout<<"====== RHRS momentum tuning ==========="<<endl;
    
  }else if(arm=="L"){
    MODE=1;
    cout<<"====== LHRS momentum tuning ==========="<<endl;    
  } else if(arm=="C" || arm=="I"){
    MODE=0;
    cout<<"====== R & L momentum tuning ==========="<<endl;    
    if(arm=="I"){
      cout<<"======    Initial matrix     ==========="<<endl;
      Initial=true;}
    }else if(arm=="Al"){
      cout<<"====== w/ Al momentum tuning ==========="<<endl;   
      Al_mode=true; 
      MODE=0;
      
    }else{cout<<"Please Select Tuning Mode !"<<endl; exit(1);}
    
    cout<<"======================================="<<endl;
    
    target=Target;

  // tdc_time tdc_time=0.56;//[ns]
    //    cout<<"F1tdc_mode : "<<typeid(F1tdc).name()<<endl;
  if(F1tdc==1){tdc_time=0.056;// [ns]
    tdc_mode=1;
    Lp_scale=false;
    //    coin_offset=464.13; // H1 mode
    coin_offset=464.73; // H1 mode
  }else if(F1tdc==2){tdc_time=0.058; //[ns]
    tdc_mode=2;
    //    coin_offset=470.13; // H2 mode
    coin_offset=470.63; // H2 mode
    //    Lp_scale=true;
  }else if(F1tdc==3){
    tdc_time=0.058; //[ns]
    tdc_mode=2;
    //    coin_offset=470.13; // H2 mode
    coin_offset=470.63; // H2 mode
    Lp_scale=true;    

  }else{cout<<"Error F1tdc modes "<<F1tdc<<endl;exit(1);}


  
  cout<<"tdc resolution [ns]"<<tdc_time<<endl;
  cout<<"Target : "<<target<<endl;
  cout<<"Momentum order : "<<nnp<< " #Param: "<<nParamTp*2<<endl;
  cout<<"Lp scale mode "<<Lp_scale<<endl;
  cout<<"Initial matrix mode : "<<Initial<<endl;
  cout<<"Al mode "<<Al_mode<<endl;
  return MODE;
}



//-----------------------------//
//--------- SetRun  -----------//
//-----------------------------//

void momcalib::SingleRoot(string ifname){

  SetRun(ifname);
  SetBranch();
   //  T->SetBranchStatus("coin",1);
  //  T->SetBranchAddress("coin",&coin_time);
  T->SetBranchStatus("ct",1);
  T->SetBranchAddress("ct",&coin_time);  

};



//-----------------------------//
//--------- SetRoot -----------//
//-----------------------------//

void momcalib::SetRoot(string ifname){

  ChainTree(ifname);
  SetBranch();
};

//-----------------------------//
//--------- NewRoot -----------//
//-----------------------------//

void momcalib::NewRoot(string ofname){

  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  tnew=new TTree("T","Momentum matrix tuning");
  tnew =T->CloneTree(0);
  tnew->Branch("mm", &mm,"mm/D");  
  tnew->Branch("mm_acc", &mm_acc,"mm_acc/D");    
  tnew->Branch("mm_L", &mm_L_b,"mm_L/D");
  tnew->Branch("mm_L_c", &mm_L,"mm_L_c/D");
  tnew->Branch("mm_nnL", &mm_nnL,"mm_nnL/D");
  tnew->Branch("Al_tuning", &Al_check,"Al_tuning/B");    
  tnew->Branch("mm_Al", &mm_Al,"mm_Al/D");
  tnew->Branch("mm_X", &mm_X,"mm_X/D");
  tnew->Branch("mm_Al_nnL", &mm_Al_nnL,"mm_Al_nnL/D");
  tnew->Branch("mm_Al_acc", &mm_Al_acc,"mm_Al_acc/D");  
  tnew->Branch("mm_Al_c", &mm_Al_c,"mm_Al_c/D");
  tnew->Branch("mm_Al_acc_c", &mm_Al_acc_c,"mm_Al_acc_c/D");

  tnew->Branch("mm_nnL", &mm_nnL,"mm_nnL/D");
  tnew->Branch("mm_nnL_acc", &mm_nnL_acc,"mm_nnL_acc/D");  
  tnew->Branch("mm_nnL_c", &mm_nnL_c,"mm_nnL_c/D");
  tnew->Branch("mm_nnL_acc_c", &mm_nnL_acc_c,"mm_nnL_acc_c/D");  
  

  tnew->Branch("mm_pi", &mm_pi,"mm_pi/D");
  tnew->Branch("coin", &coin_t,"coin_t/D");  
  tnew->Branch("coin_acc", &ctime,"coin_acc/D");  
  tnew->Branch("tflag", &tflag,"tflag/I");  
  tnew->Branch("tev", &tev,"tev/I");
  tnew->Branch("dpe",&Dpe,"dpe/D");
  tnew->Branch("dpe_",&Dpe_,"dpe_/D");
  tnew->Branch("dpk",&Dpk,"dpk/D");
};

//=====================================================================//

void momcalib::MakeHist(){

  
  min_z=-0.2;
  max_z=0.2;
  bin_z=400;
  hLz=new TH1F("hLz","",bin_z,min_z,max_z);
  set->SetTH1(hLz,"LHRS z-vertex w/ correction","z [m]","Counts");
  hLz->SetLineColor(4);
  hLz->SetFillColor(4);  
  hLz->SetFillStyle(3002);
  hLz_c=new TH1F("hLz_c","",bin_z,min_z,max_z);  
  set->SetTH1(hLz_c,"LHRS z-vertex z parameter tuing","z [m]","Counts");
  hLz_c->SetLineColor(8);
  hLz_c->SetFillColor(8);  
  hLz_c->SetFillStyle(3002);
  hLz_rc=new TH1F("hLz_rc","",bin_z,min_z,max_z);  
  set->SetTH1(hLz_rc,"LHRS z-vertex w/ raster & parameters correction","z [m]","Counts");
  hLz_rc->SetLineColor(2);
  hLz_rc->SetFillColor(2);  
  hLz_rc->SetFillStyle(3002);
  
  hRz=new TH1F("hRz","",bin_z,min_z,max_z);
  set->SetTH1(hRz,"RHRS z-vertex ","z [m]","Counts");
  hRz->SetLineColor(4);
  hRz->SetFillColor(4);  
  hRz->SetFillStyle(3002);
  hRz_c=new TH1F("hRz_c","",bin_z,min_z,max_z);  
  set->SetTH1(hRz_c,"RHRS z-vertex z parameter tuing","z [m]","Counts");
  hRz_c->SetLineColor(8);
  hRz_c->SetFillColor(8);  
  hRz_c->SetFillStyle(3002);
  hRz_rc=new TH1F("hRz_rc","",bin_z,min_z,max_z);  
  set->SetTH1(hRz_rc,"RHRS z-vertex w/ raster & parameters correction","z [m]","Counts");
  hRz_rc->SetLineColor(2);
  hRz_rc->SetFillColor(2);  
  hRz_rc->SetFillStyle(3002);


  //===== Event Selection ==============//
  hLz_es=new TH1F("hLz_es","",bin_z,min_z,max_z);
  set->SetTH1(hLz_es,"LHRS z-vertex w/ correction","z [m]","Counts");
  hLz_es->SetLineColor(4);
  hLz_es->SetFillColor(4);  
  hLz_es->SetFillStyle(3002);
  hLz_es_c=new TH1F("hLz_es_c","",bin_z,min_z,max_z);  
  set->SetTH1(hLz_es_c,"LHRS z-vertex z parameter tuing","z [m]","Counts");
  hLz_es_c->SetLineColor(8);
  hLz_es_c->SetFillColor(8);  
  hLz_es_c->SetFillStyle(3002);
  hLz_es_rc=new TH1F("hLz_es_rc","",bin_z,min_z,max_z);  
  set->SetTH1(hLz_es_rc,"LHRS z-vertex w/ raster & parameters correction","z [m]","Counts");
  hLz_es_rc->SetLineColor(2);
  hLz_es_rc->SetFillColor(2);  
  hLz_es_rc->SetFillStyle(3002);
  
  hRz_es=new TH1F("hRz_es","",bin_z,min_z,max_z);
  set->SetTH1(hRz_es,"RHRS z-vertex ","z [m]","Counts");
  hRz_es->SetLineColor(4);
  hRz_es->SetFillColor(4);  
  hRz_es->SetFillStyle(3002);
  hRz_es_c=new TH1F("hRz_es_c","",bin_z,min_z,max_z);  
  set->SetTH1(hRz_es_c,"RHRS z-vertex z parameter tuing","z [m]","Counts");
  hRz_es_c->SetLineColor(8);
  hRz_es_c->SetFillColor(8);  
  hRz_es_c->SetFillStyle(3002);
  hRz_es_rc=new TH1F("hRz_es_rc","",bin_z,min_z,max_z);  
  set->SetTH1(hRz_es_rc,"RHRS z-vertex w/ raster & parameters correction","z [m]","Counts");
  hRz_es_rc->SetLineColor(2);
  hRz_es_rc->SetFillColor(2);  
  hRz_es_rc->SetFillStyle(3002);

  
  hz_Al=new TH1F("hz_Al","",bin_z,min_z,max_z);
  set->SetTH1(hz_Al," z-vertex w/ correction","z [m]","Counts");
  hz_Al->SetLineColor(4);
  hz_Al->SetFillColor(4);  
  hz_Al->SetFillStyle(3002);

  
  min_th=-0.1;
  max_th=0.1;
  bin_th=1000;
  
  hLth=new TH1F("hLth","",bin_th,min_th,max_th);
  set->SetTH1(hLth,"LHRS th-vertex ","theta [rad]","Counts");
  hLth->SetLineColor(4);
  hLth->SetFillColor(4);  
  hLth->SetFillStyle(3002);
  hLth_c=new TH1F("hLth_c","",bin_th,min_th,max_th);  
  set->SetTH1(hLth_c,"LHRS theta-vertex w/ parameter tuning","theta [rad]","Counts");
  hLth_c->SetLineColor(2);
  hLth_c->SetFillColor(2);  
  hLth_c->SetFillStyle(3002);
  hRth=new TH1F("hRth","",bin_th,min_th,max_th);
  set->SetTH1(hRth,"RHRS th-vertex ","theta [rad]","Counts");
  hRth->SetLineColor(4);
  hRth->SetFillColor(4);  
  hRth->SetFillStyle(3002);  
  hRth_c=new TH1F("hRth_c","",bin_th,min_th,max_th);  
  set->SetTH1(hRth_c,"RHRS th-vertex w/ paramater tuning","theta [rad]","Counts");  
  hRth_c->SetLineColor(2);
  hRth_c->SetFillColor(2);  
  hRth_c->SetFillStyle(3002);


  //====== event selection ========//

  hLth_es=new TH1F("hLth_es","",bin_th,min_th,max_th);
  set->SetTH1(hLth_es,"LHRS th-vertex ","theta [rad]","Counts");
  hLth_es->SetLineColor(4);
  hLth_es->SetFillColor(4);  
  hLth_es->SetFillStyle(3002);
  hLth_es_c=new TH1F("hLth_es_c","",bin_th,min_th,max_th);  
  set->SetTH1(hLth_es_c,"LHRS theta-vertex w/ parameter tuning","theta [rad]","Counts");
  hLth_es_c->SetLineColor(2);
  hLth_es_c->SetFillColor(2);  
  hLth_es_c->SetFillStyle(3002);
  hRth_es=new TH1F("hRth_es","",bin_th,min_th,max_th);
  set->SetTH1(hRth_es,"RHRS th-vertex ","theta [rad]","Counts");
  hRth_es->SetLineColor(4);
  hRth_es->SetFillColor(4);  
  hRth_es->SetFillStyle(3002);  
  hRth_es_c=new TH1F("hRth_es_c","",bin_th,min_th,max_th);  
  set->SetTH1(hRth_es_c,"RHRS th-vertex w/ paramater tuning","theta [rad]","Counts");  
  hRth_es_c->SetLineColor(2);
  hRth_es_c->SetFillColor(2);  
  hRth_es_c->SetFillStyle(3002);

  
  
  
  min_ph=-0.1;
  max_ph=0.1;
  bin_ph=1000;

  hLph=new TH1F("hLph","",bin_ph,min_ph,max_ph);
  set->SetTH1(hLph,"LHRS ph-vertex ","phi [rad]","Counts");
  hLph->SetLineColor(4);
  hLph->SetFillColor(4);  
  hLph->SetFillStyle(3002);  
  hLph_c=new TH1F("hLph_c","",bin_ph,min_ph,max_ph);  
  set->SetTH1(hLph_c,"LHRS ph-vertex w/ parameter tuning","phi [rad]","Counts");
  hLph_c->SetLineColor(2);
  hLph_c->SetFillColor(2);  
  hLph_c->SetFillStyle(3002);
  hRph=new TH1F("hRph","",bin_ph,min_ph,max_ph);
  set->SetTH1(hRph,"RHRS ph-vertex ","phi [rad]","Counts");
  hRph->SetLineColor(4);
  hRph->SetFillColor(4);  
  hRph->SetFillStyle(3002);  
  hRph_c=new TH1F("hRph_c","",bin_ph,min_ph,max_ph);  
  set->SetTH1(hRph_c,"RHRS ph-vertex parameter tuning","phi [rad]","Counts");  
  hRph_c->SetLineColor(2);
  hRph_c->SetFillColor(2);  
  hRph_c->SetFillStyle(3002);


  //======= EventSelection ===============//

  hLph_es=new TH1F("hLph_es","",bin_ph,min_ph,max_ph);
  set->SetTH1(hLph_es,"LHRS ph-vertex ","phi [rad]","Counts");
  hLph_es->SetLineColor(4);
  hLph_es->SetFillColor(4);  
  hLph_es->SetFillStyle(3002);  
  hLph_es_c=new TH1F("hLph_es_c","",bin_ph,min_ph,max_ph);  
  set->SetTH1(hLph_es_c,"LHRS ph-vertex w/ parameter tuning","phi [rad]","Counts");
  hLph_es_c->SetLineColor(2);
  hLph_es_c->SetFillColor(2);  
  hLph_es_c->SetFillStyle(3002);
  hRph_es=new TH1F("hRph_es","",bin_ph,min_ph,max_ph);
  set->SetTH1(hRph_es,"RHRS ph-vertex ","phi [rad]","Counts");
  hRph_es->SetLineColor(4);
  hRph_es->SetFillColor(4);  
  hRph_es->SetFillStyle(3002);  
  hRph_es_c=new TH1F("hRph_es_c","",bin_ph,min_ph,max_ph);  
  set->SetTH1(hRph_es_c,"RHRS ph-vertex parameter tuning","phi [rad]","Counts");  
  hRph_es_c->SetLineColor(2);
  hRph_es_c->SetFillColor(2);  
  hRph_es_c->SetFillStyle(3002);


  
  
  min_ssx = -8.0; //[cm]
  max_ssx =  8.0; //[cm]
  min_ssy =  -5.0; //[cm]
  max_ssy =  5.0; //[cm]
  bin_ssx =(int)((max_ssx - min_ssx)/0.05);
  bin_ssy =(int)((max_ssy - min_ssy)/0.05);

  hLss=new TH2F("hLss","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hLss,"Sieve Slit pattern in LHRS w/o parameter tuning","ssy [cm]","ssx [cm]");
  hLss_c=new TH2F("hLss_c","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hLss_c,"Sieve Slit pattern in LHRS w/  parameter tuning","ssy [cm]","ssx [cm]");

  hRss=new TH2F("hRss","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hRss,"Sieve Slit pattern in RHRS w/o parameter tuning","ssy [cm]","ssx [cm]");
  hRss_c=new TH2F("hRss_c","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hRss_c,"Sieve Slit pattern in RHRS w/  parameter tuning","ssy [cm]","ssx [cm]");


  //====== Event Slection hist =====================//
  hLss_es=new TH2F("hLss_es","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hLss_es,"Sieve Slit pattern in LHRS w/o parameter tuning","ssy [cm]","ssx [cm]");
  hLss_es_c=new TH2F("hLss_es_c","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hLss_es_c,"Sieve Slit pattern in LHRS w/  parameter tuning","ssy [cm]","ssx [cm]");

  hRss_es=new TH2F("hRss_es","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hRss_es,"Sieve Slit pattern in RHRS w/o parameter tuning","ssy [cm]","ssx [cm]");
  hRss_es_c=new TH2F("hRss_es_c","",bin_ssy,min_ssy,max_ssy,bin_ssx,min_ssx,max_ssx);
  set->SetTH2(hRss_es_c,"Sieve Slit pattern in RHRS w/  parameter tuning","ssy [cm]","ssx [cm]");

  
  
  min_lp = 1.8;
  max_lp = 2.5;
  bin_lp =(int)((max_lp - min_lp)*10000.);

  hLp=new TH1F("hLp","",bin_lp,min_lp,max_lp);
  set->SetTH1(hLp,"LHRS momentum w/o parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hLp->SetLineColor(4);
  hLp->SetFillColor(4);  
  hLp->SetFillStyle(3002);
  hLp_c=new TH1F("hLp_c","",bin_lp,min_lp,max_lp);
  set->SetTH1(hLph_c,"LHRS momentum w/ parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hLp_c->SetLineColor(8);
  hLp_c->SetFillColor(8);  
  hLp_c->SetFillStyle(3002);
  hLp_ec=new TH1F("hLp_ec","",bin_lp,min_lp,max_lp);
  set->SetTH1(hLp_ec,"LHRS momentum w/ parameter tuning & Energy loss correction"," p [GeV/c]","Counts/ 100 keV ");
  hLp_ec->SetLineColor(2);
  hLp_ec->SetFillColor(2);  
  hLp_ec->SetFillStyle(3002);
  
  min_rp = 1.5;
  max_rp = 2.3;
  bin_rp =(int)((max_rp - min_rp)*1000);

  hRp=new TH1F("hRp","",bin_rp,min_rp,max_rp);
  set->SetTH1(hRph,"RHRS momentum w/o parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hRp->SetLineColor(4);
  hRp->SetFillColor(4);  
  hRp->SetFillStyle(3002);
  hRp_c=new TH1F("hRp_c","",bin_rp,min_rp,max_rp);
  set->SetTH1(hRph_c,"RHRS momentum w/ parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hRp_c->SetLineColor(8);
  hRp_c->SetFillColor(8);  
  hRp_c->SetFillStyle(3002);
  hRp_ec=new TH1F("hRp_ec","",bin_rp,min_rp,max_rp);
  set->SetTH1(hRp_ec,"RHRS momentum w/ parameter tuning & Energy loss correction"," p [GeV/c]","Counts/ 100 keV ");
  hRp_ec->SetLineColor(2);
  hRp_ec->SetFillColor(2);  
  hRp_ec->SetFillStyle(3002);  


  //====== Event Slection hist =====================//

  hLp_es=new TH1F("hLp_es","",bin_lp,min_lp,max_lp);
  set->SetTH1(hLp_es,"LHRS momentum w/o parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hLp_es->SetLineColor(4);
  hLp_es->SetFillColor(4);  
  hLp_es->SetFillStyle(3002);
  hLp_es_c=new TH1F("hLp_es_c","",bin_lp,min_lp,max_lp);
  set->SetTH1(hLp_es_c,"LHRS momentum w/ parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hLp_es_c->SetLineColor(8);
  hLp_es_c->SetFillColor(8);  
  hLp_es_c->SetFillStyle(3002);
  hLp_es_ec=new TH1F("hLp_es_ec","",bin_lp,min_lp,max_lp);
  set->SetTH1(hLp_es_ec,"LHRS momentum w/ parameter tuning & Energy loss correction"," p [GeV/c]","Counts/ 100 keV ");
  hLp_es_ec->SetLineColor(2);
  hLp_es_ec->SetFillColor(2);  
  hLp_es_ec->SetFillStyle(3002);

  hRp_es=new TH1F("hRp_es","",bin_rp,min_rp,max_rp);
  set->SetTH1(hRp_es,"RHRS momentum w/o parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hRp_es->SetLineColor(4);
  hRp_es->SetFillColor(4);  
  hRp_es->SetFillStyle(3002);
  hRp_es_c=new TH1F("hRp_es_c","",bin_rp,min_rp,max_rp);
  set->SetTH1(hRp_es_c,"RHRS momentum w/ parameter tuning "," p [GeV/c]","Counts/ 100 keV ");
  hRp_es_c->SetLineColor(8);
  hRp_es_c->SetFillColor(8);  
  hRp_es_c->SetFillStyle(3002);
  hRp_es_ec=new TH1F("hRp_es_ec","",bin_rp,min_rp,max_rp);
  set->SetTH1(hRp_es_ec,"RHRS momentum w/ parameter tuning & Energy loss correction"," p [GeV/c]","Counts/ 100 keV ");
  hRp_es_ec->SetLineColor(2);
  hRp_es_ec->SetFillColor(2);  
  hRp_es_ec->SetFillStyle(3002);  

  
  
  min_ct=-50.0;
  max_ct=50.0;
  //  bin_ct=(double)((max_ct-min_ct)/tdc_time);
  bin_ct=4000;

  hct=new TH1F("hct","",bin_ct,min_ct,max_ct);
  set->SetTH1(hct,"Coincidence time","coin-time [ns]","Counts");  
  
  hct2=new TH1F("hct2","",bin_ct,min_ct,max_ct);
  set->SetTH1(hct2,"Coincidence time","coin-time [ns]","Counts");  

  
  //----- Cut hist --------//
  hLz_cut=new TH1F("hLz_cut","",bin_z,min_z,max_z);
  set->SetTH1(hLz_cut,"LHRS z cut event","z [m]","Counts");
  hLz_cut->SetLineColor(2);
  hLz_cut->SetFillColor(2);  
  hLz_cut->SetFillStyle(3002);

  hRz_cut=new TH1F("hRz_cut","",bin_z,min_z,max_z);
  set->SetTH1(hRz_cut,"RHRS z cut event","z [m]","Counts");
  hRz_cut->SetLineColor(2);
  hRz_cut->SetFillColor(2);  
  hRz_cut->SetFillStyle(3002);    

 
  hct_cut=new TH1F("hct_cut","",bin_ct,min_ct,max_ct);
  set->SetTH1(hct,"Coincidence time w/ cut","coin-time [ns]","Counts");  

  hct_cut->SetLineColor(2);
  hct_cut->SetFillColor(2);  
  hct_cut->SetFillStyle(3002);  


  
  min_mm = -100;// [MeV]
  max_mm = +200; // [MeV]
  int mm_width=1; // MeV
  double mmbin_width=(double)mm_width; //[MeV]
  bin_mm=(int)((max_mm-min_mm)/mmbin_width);

  hmm_select=new TH1F("hmm","",bin_mm,min_mm,max_mm);
  set->SetTH1(hmm_select,"Missing mass Event Selction","M_{X}-M_{#Lambda} [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));
  hmm_cut=new TH1F("hmm_cut","",bin_mm,min_mm,max_mm);
  set->SetTH1(hmm_cut,"Missing mass  w/ cut hist ","M_{X}-M_{#Lambda} [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));

  hmm_cut->SetLineColor(1);
  hmm_cut->SetFillColor(6);  
  hmm_cut->SetFillStyle(3002);


  hmm_b=new TH1F("hmm_b","",bin_mm,min_mm,max_mm);
  set->SetTH1(hmm_b,"Missing mass  w/ cut hist ","M_{X}-M_{#Lambda} [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));

  hmm_b->SetLineColor(1);
  hmm_b->SetFillColor(5);  
  hmm_b->SetFillStyle(3002);

  hmm_cut_H=new TH1F("hmm_cut_H","",bin_mm,min_mm,max_mm);
  set->SetTH1(hmm_cut_H,"Missing mass  w/ cut hist ","M_{X}-M_{#Lambda} [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));

  hmm_cut_H->SetLineColor(1);
  hmm_cut_H->SetFillColor(2);  
  hmm_cut_H->SetFillStyle(3002);

  hmm_cut_T=new TH1F("hmm_cut_T","",bin_mm,min_mm,max_mm);
  set->SetTH1(hmm_cut_T,"Missing mass  w/ cut hist ","M_{X}-M_{#Lambda} [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));

  hmm_cut_T->SetLineColor(1);
  hmm_cut_T->SetFillColor(3);  
  hmm_cut_T->SetFillStyle(3002);




  //------ Event Selection hist ----------------//
  hmm_L=new TH1F("hmm_L","",bin_mm,min_mm,max_mm);
  set->SetTH1(hmm_L,"Lambda missing mass Event Selction","M_{X}-M_{#Lambda}  [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));  
  hmm_L->SetLineColor(2);
  hmm_L->SetFillColor(2);  
  hmm_L->SetFillStyle(3002);

  hmm_S=new TH1F("hmm_S","",bin_mm,min_mm,max_mm);    
  set->SetTH1(hmm_S,"Sigma missing mass Event Selction","M_{X}-M_{#Lambda} [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));
  hmm_S->SetLineColor(kBlue);
  hmm_S->SetFillColor(kBlue);  
  hmm_S->SetFillStyle(3002);

  hmm_Al=new TH1F("hmm_Al","",600,-300,300);
  set->SetTH1(hmm_Al,"Missing mass Event Selction","M_{X}-M_{core+#Lambda} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));

  hmm_Al_acc=new TH1F("hmm_Al_acc","",600,-300,300);
  set->SetTH1(hmm_Al_acc,"Missing mass Event Selction","M_{X}-M_{core+#Lambda} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));


  hmm_Al_select=new TH1F("hmm_Al_select","",600,-300,300);
  set->SetTH1(hmm_Al_select,"Missing mass Event Selction","M_{X}-M_{core+#Lambda} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));
  hmm_Al_select->SetLineColor(kBlue);
  hmm_Al_select->SetFillColor(kBlue);  
  hmm_Al_select->SetFillStyle(3002);


  hmm_Al_select_c=new TH1F("hmm_Al_select_c","",600,-300,300);
  set->SetTH1(hmm_Al_select_c,"Missing mass Event Selction","M_{X}-M_{core+#Lambda} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));
  hmm_Al_select_c->SetLineColor(8);
  hmm_Al_select_c->SetFillColor(8);  
  hmm_Al_select_c->SetFillStyle(3002);

  hmm_Al_cut=new TH1F("hmm_Al_cut","",600,-300,300);
  set->SetTH1(hmm_Al_cut,"Missing mass Event Selction","M_{X}-M_{core+#Lambda} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));
  hmm_Al_cut->SetLineColor(kRed);
  hmm_Al_cut->SetFillColor(kRed);  
  hmm_Al_cut->SetFillStyle(3002);

  hmm_Al_cut_acc=new TH1F("hmm_Al_cut_acc","",600,-300,+300);
  set->SetTH1(hmm_Al_cut_acc,"Missing mass Event Selction","M_{X}-M_{core+#Lambda} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));
  hmm_Al_cut_acc->SetLineColor(kRed);
  hmm_Al_cut_acc->SetFillColor(kRed);  
  hmm_Al_cut_acc->SetFillStyle(3002);



  hmm_Al_zp=new TH1F("hmm_Al_zp","",600,-300,300);
  set->SetTH1(hmm_Al_zp,"Al missing mass z>0 events","M_{X}-M_{core+#Lambda} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));

  hmm_Al_zn=new TH1F("hmm_Al_zn","",600,-300,300);
  set->SetTH1(hmm_Al_zn,"Al missing mass z<0 events","M_{X}-M_{core+#Lambda} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));


  hmm_Al_zp_c=new TH1F("hmm_Al_zp_c","",600,-300,300);
  set->SetTH1(hmm_Al_zp_c,"Al missing mass z>0 events after Al tuning","M_{X}-M_{core+#Lambda} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));

  hmm_Al_zn_c=new TH1F("hmm_Al_zn_c","",600,-300,300);
  set->SetTH1(hmm_Al_zn_c,"Al missing mass z<0 events afeter Al tuning","M_{X}-M_{core+#Lambda} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));  
  

  
  hmm_Al_mom=new TH1F("hmm_Al_mom","",300,-0300,300);
  set->SetTH1(hmm_Al_mom,"Missing mass Event Selction","M_{X}-M_{core+#Lambda} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));



  hmm_Al_nnL=new TH2F("hmm_Al_nnL","",100,-50,50,100,-50,50);
  set->SetTH2(hmm_Al_nnL,"Missing mass Al vs nnL ","Al missing mass[MeV/c^{2}]","nnL missing mass[MeV/c^{2}]");


  hmm_Al_nnL_zp=new TH2F("hmm_Al_nnL_zp","",100,-50,50,100,-50,50);
  set->SetTH2(hmm_Al_nnL_zp,"Missing mass Al vs nnL ","Al missing mass[MeV/c^{2}]","nnL missing mass[MeV/c^{2}]");

  hmm_Al_nnL_zn=new TH2F("hmm_Al_nnL_zn","",100,-50,50,100,-50,50);
  set->SetTH2(hmm_Al_nnL_zn,"Missing mass Al vs nnL ","Al missing mass[MeV/c^{2}]","nnL missing mass[MeV/c^{2}]");  
  

  
  hmm_nnL=new TH1F("hmm_nnL","",600,-300,300);
  set->SetTH1(hmm_nnL,"nnL missing mass Event Selction","M_{X}-M_{nnL} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));  

  //  hmm_nnL->SetFillStyle(3002);

  hmm_nnL_acc=new TH1F("hmm_nnL_acc","",600,-300,300);
  set->SetTH1(hmm_nnL_acc,"nnL_acc missing mass Event Selction","M_{X}-M_{nnL}  [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));  
  hmm_nnL->SetLineColor(1);
  //  hmm_nnL->SetFillColor(2);  
  //  hmm_nnL->SetFillStyle(3002);

  hmm_nnL_select=new TH1F("hmm_nnL_select","",600,-300,300);
  set->SetTH1(hmm_nnL_select,"nnL_select missing mass Event Selction","M_{X}-M_{nnL_select} [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));  
  hmm_nnL_select->SetLineColor(kBlue);
  hmm_nnL_select->SetFillColor(kBlue);  
  hmm_nnL_select->SetFillStyle(3002);



  hmm_nnL_cut=new TH1F("hmm_nnL_cut","",600,-300,300);
  set->SetTH1(hmm_nnL_cut,"nnL missing mass Event Selction","M_{X}-M_{nnL}  [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));  

  hmm_nnL_cut->SetLineColor(4);
  hmm_nnL_cut->SetFillColor(4);  
  hmm_nnL_cut->SetFillStyle(3002);


  hmm_nnL_cut_acc=new TH1F("hmm_nnL_cut_acc","",600,-300,300);
  set->SetTH1(hmm_nnL_cut_acc,"nnL missing mass Event Selction","M_{X}-M_{nnL}  [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));  

  hmm_nnL_cut_acc->SetLineColor(4);
  hmm_nnL_cut_acc->SetFillColor(4);  
  hmm_nnL_cut_acc->SetFillStyle(3002);



  hmm_nnL_select=new TH1F("hmm_nnL_select","",600,-300,300);
  set->SetTH1(hmm_nnL_select,"nnL missing mass Event Selction","M_{X}-M_{nnL}  [MeV/c^{2}]",Form("Counts/%d MeV",mm_width));  

  hmm_nnL_select->SetLineColor(6);
  hmm_nnL_select->SetFillColor(6);  
  hmm_nnL_select->SetFillStyle(3002);




 //------- tuning check ----------------------//
  for(int i=0;i<nite;i++){

    hmm_L_t[i]=new TH1F(Form("hmm_L_t_%d",i),"",bin_mm,min_mm,max_mm);
    set->SetTH1(hmm_L_t[i],Form("Lambda tuning hist # %d",i),"mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));
    hmm_L_tc[i]=new TH1F(Form("hmm_L_tc_%d",i),"",bin_mm,min_mm,max_mm);
    set->SetTH1(hmm_L_tc[i],Form("Lambda tuning hist # %d",i),"mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));
    hmm_S_t[i]=new TH1F(Form("hmm_S_t_%d",i),"",bin_mm,min_mm,max_mm);
    set->SetTH1(hmm_S_t[i],Form("Lambda tuning hist # %d",i),"mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));
    hmm_S_tc[i]=new TH1F(Form("hmm_S_tc_%d",i),"",bin_mm,min_mm,max_mm);
    set->SetTH1(hmm_S_tc[i],Form("Lambda tuning hist # %d",i),"mass [GeV/c^{2}]",Form("Counts/%d MeV",mm_width));        
    
    }
  

  
  
}

//==========================================================

void momcalib::MTParam(string mtparam){

  cout<<endl;
  cout<<"==============================="<<endl;
  cout<<"=== Input Matrix Parameters ==="<<endl;
  cout<<"==============================="<<endl;

  string buf;
  int s=0;
  ifstream ifp(Form("%s",mtparam.c_str()),ios::in);
  if (ifp.fail()){ cerr << "failed open files" <<mtparam.c_str()<<endl; exit(1);}
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >>param[s];
    cout<<param[s]<<endl;
    s++;
  }

  if(MODE==0 || MODE==-1) {MTParam_R();cout<<" Input RHRS Matrix parameter "<<endl;}
  if(MODE==0 || MODE==1) {MTParam_L();cout<<" Input LHRS Matrix parameter "<<endl;}
  MTP_mom();

  for(int i=0;i<10;i++){MT[i]=false; MT_f[i]=false;}

  
  //======= Tuning selection flag =====================//
  //--------- RHRS ------------------------//

  if(Initial){
  MT[0] = true;  // RHRS z correction
  MT[1] = true;  // RHRS raster correction
  MT[2] = true;  // RHRS theta correction
  MT[3] = true;  // RHRS phi correction
  //--------- LHRS -----------------------//
  MT[4] = true;  // LHRS z correction
  MT[5] = true;  // LHRS raster correction
  MT[6] = true;  // LHRS theta correction
  MT[7] = true;  // LHRS phi correction
  //-------- momentum calibration ---------//
  MT[8] = true; // RHRS momentum correction  
  MT[9] = true; // LHRS momentum correction  

  //================================================//

  MT_f[0] = true;  // RHRS z correction
  MT_f[1] = true;  // RHRS raster correction
  MT_f[2] = true;  // RHRS theta correction
  MT_f[3] = true;  // RHRS phi correction
  //--------- LHRS -----------------------//
  MT_f[4] = true;  // LHRS z correction
  MT_f[5] = true;  // LHRS raster correction
  MT_f[6] = true;  // LHRS theta correction
  MT_f[7] = true;  // LHRS phi correction

  MT[8] = false; // RHRS momentum correction  
  MT[9] = false; // LHRS momentum correction  
  ploss = false;  // Energy Loss 
 
  //-------- momentum calibration ---------//
  MT_f[8] = false; // RHRS momentum correction  
  MT_f[9] = false; // LHRS momentum correction  
  ploss_f = false;  // Energy Loss

  }else{
  
    MT[8] = true; // RHRS momentum correction  
    MT[9] = true; // LHRS momentum correction  
    ploss = true;  // Energy Loss 
    
    //-------- momentum calibration ---------//
    MT_f[8] = true; // RHRS momentum correction  
    MT_f[9] = true; // LHRS momentum correction  
    ploss_f = true;  // Energy Loss
    
  }
  //================================================//

   
  
  cout<<endl;
  cout<<"======== Correction Parameters ========="<<endl;
  cout<<"          < Event Selection >                                < Fill >"<<endl;;
  if(MT[0])cout<<" RHRS z      correction                  ";
  else     cout<<" RHRS z                    no correction ";
if(MT_f[0])cout<<"             correction                  "<<endl;  
 else      cout<<"                           no correction "<<endl;  
  if(MT[1])cout<<" RHRS raster correction                  ";
  else     cout<<" RHRS raster               no correction ";
if(MT_f[1])cout<<"             correction                  "<<endl;  
 else      cout<<"                           no correction "<<endl;  
  
  if(MT[2])cout<<" RHRS theta  correction                  ";
  else     cout<<" RHRS theta                no correction ";
if(MT_f[2])cout<<"             correction                  "<<endl;  
 else      cout<<"                           no correction "<<endl;  
  if(MT[3])cout<<" RHRS phi    correction                  ";
  else     cout<<" RHRS phi                  no correction ";
if(MT_f[3])cout<<"             correction                  "<<endl;  
 else      cout<<"                           no correction "<<endl;
  if(MT[4])cout<<" LHRS z      correction                  ";
  else     cout<<" LHRS z                    no correction ";
if(MT_f[4])cout<<"             correction                  "<<endl;  
 else      cout<<"                           no correction "<<endl;    
  if(MT[5])cout<<" LHRS raster correction                  ";
  else     cout<<" LHRS raster               no correction ";
  if(MT_f[5])cout<<"             correction                  "<<endl;  
 else      cout<<"                           no correction "<<endl;  
  if(MT[6])cout<<" LHRS theta  correction                  ";
  else     cout<<" LHRS theta                no correction ";
if(MT_f[6])cout<<"             correction                  "<<endl;  
 else      cout<<"                           no correction "<<endl;  
  if(MT[7])cout<<" LHRS phi    correction                  ";
  else     cout<<" LHRS phi                  no correction ";
if(MT_f[7])cout<<"             correction                  "<<endl;  
 else      cout<<"                           no correction "<<endl;  
  if(MT[8])cout<<" RHRS mom    correction                  ";
  else     cout<<" RHRS mom                  no correction ";
if(MT_f[8])cout<<"             correction                  "<<endl;  
 else      cout<<"                           no correction "<<endl;  
  if(MT[9])cout<<" LHRS mom    correction                  ";
  else     cout<<" LHRS mom                  no correction ";
if(MT_f[9])cout<<"             correction                  "<<endl;  
 else      cout<<"                           no correction "<<endl;  
  if(ploss)cout<<" Energy Los  correction                  ";
  else     cout<<" Energy Los                no correction ";
if(ploss_f)cout<<"             correction                  "<<endl;  
 else      cout<<"                           no correction "<<endl;  
  cout<<endl;
}



//=====================================================================//

void momcalib::MTParam_R(){

  //=================//
  //==== RHRS =======//
  //=================//

  
  //====== RHRS z parameters ======//

    char name_Mzt[500];
    sprintf(name_Mzt, param[0].c_str()); // optimized
    ifstream Mzt(name_Mzt);

   for(int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt >> par >> p >> p >> p >> p;
    Pzt[i]=par;
    //    cout<<"R Mzt : "<<Pzt[i]<<endl;
   }
  Mzt.close();

  
  //====== RHRS raster paramters =======//
    char name_Mras[500];
    sprintf(name_Mras, param[1].c_str()); // optimized
    //    cout<<"RHRS Raster parameters file: "<<name_Mras<<endl;
  ifstream Mras(name_Mras);

  for (int i=0;i<nParamT_ras;i++){

    Mras >> Pras[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) Pras[i] = 0.0;
  }


  Mras.close();    

  
  //===== RHRS theta parameters ======// 
    char name_Mxpt[500];
    sprintf(name_Mxpt, param[2].c_str()); // optimized
  ifstream Mxpt(name_Mxpt);

    for(int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt >> par >> p >> p >> p >> p >> p;
    Pxpt[i]  = par;
  }
  Mxpt.close();  
 //===== RHRS phi parameters ======//
  char name_Mypt[500];
    sprintf(name_Mypt, param[3].c_str()); // optimized  
  ifstream Mypt(name_Mypt);

  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mypt >> par >> p >> p >> p >> p >> p;
    Pypt[i]  = par;
  }
  Mypt.close();    



  
};
//////////////////////////////////////////////////////////////

void momcalib::MTParam_L(){

  //=================//
  //===== LHRS ======//
  //=================//

  
  //====== LHRS z parameters ======//  
  char name_Mzt_L[500];
  sprintf(name_Mzt_L,param[4].c_str()); // optimized
  ifstream Mzt_L(name_Mzt_L); 
  for (int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt_L >> par >> p >> p >> p >> p;
    Pzt_L[i]=par;
  }
 Mzt_L.close();

  //====== LHRS raster paramters =======//
    char name_Mras_L[500];
    sprintf(name_Mras_L, param[5].c_str()); // optimized
    //    cout<<"LHRS Raster parameters file: "<<name_Mras_L<<endl;
  ifstream Mras_L(name_Mras_L);

  for (int i=0;i<nParamT_ras;i++){

    Mras_L >> Pras_L[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) Pras_L[i] = 0.0;
  }
  
  Mras_L.close();    

 
 //===== LHRS theta parameters ======// 
  char name_Mxpt_L[500];
    sprintf(name_Mxpt_L, param[6].c_str()); // optimized
  ifstream Mxpt_L(name_Mxpt_L);
  //  cout<<"LHRS theta parameters file: "<<name_Mxpt_L<<endl;  
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt_L >> par >> p >> p >> p >> p >> p;
    //   cout<<"LHRS theta : "<<par<<endl;
    Pxpt_L[i]  = par;
  }
  Mxpt_L.close();

  
 //===== LHRS phi parameters ===x==//
  char name_Mypt_L[500];
    sprintf(name_Mypt_L, param[7].c_str()); // optimized
  ifstream Mypt_L(name_Mypt_L);

  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mypt_L >> par >> p >> p >> p >> p >> p;
    Pypt_L[i]  = par;    
  }
  Mypt_L.close();    


}

//========================================================//


void momcalib::MTP_mom(){

  //====== RHRS Momentum parameters ========//
    char name_Mpt[500];
    sprintf(name_Mpt, param[8].c_str()); // optimized
    ifstream Mpt(name_Mpt);
    if (Mpt.fail() && Initial==0){ cerr << "failed open files" <<name_Mpt<<endl; exit(1);}
   for(int i=0;i<nParamTp;i++){
    double par=0.;
    int p=0;
    Mpt >> par >> p >> p >> p >> p >> p;
    Prp[i]=par;
    Opt_par_R[i]=par;
    Opt_par[i]=par;
   }

 
   Mpt.close();

   

  //====== LHRS Momentum parameters ========//
    char name_Mpt_L[500];
    sprintf(name_Mpt_L, param[9].c_str()); // optimized
    ifstream Mpt_L(name_Mpt_L);
    if (Mpt_L.fail()  && Initial==0){ cerr << "failed open files" <<name_Mpt_L<<endl; exit(1);}
   for(int i=0;i<nParamTp;i++){
    double par=0.;
    int p=0;
    Mpt_L >> par >> p >> p >> p >> p >> p;
    Plp[i]=par;
    Opt_par_L[i]=par;
    Opt_par[i+nParamTp]=par;  // Both momentum paramters
    //cout<<Form("Opt_par[%d]",i)<<Opt_par[i+nParamTp]<<endl;
   }
  Mpt_L.close();



}

//===============================================================================//

void momcalib::zCorr(bool rarm, bool larm){


  //======== Nomalization =============//
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Rz[0]   = (Rz[0]-Ztm)/Ztr;
    Lx_fp[0]  = (Lx_fp[0]-XFPm)/XFPr; 
    Lth_fp[0] = (Lth_fp[0]-XpFPm)/XpFPr;
    Ly_fp[0]  = (Ly_fp[0]-YFPm)/YFPr;
    Lph_fp[0] = (Lph_fp[0]-YpFPm)/YpFPr;
    Lz[0]   = (Lz[0]-Ztm)/Ztr;    
  //===================================//    

    
    if(rarm) Rz[0]   = calcf2t_zt(Pzt, Rx_fp[0], Rth_fp[0], Ry_fp[0], Rph_fp[0]); // nomalized
    if(larm) Lz[0]   = calcf2t_zt(Pzt_L, Lx_fp[0], Lth_fp[0], Ly_fp[0], Lph_fp[0]); //nomalized
    

    //========== Scaled =================//
    Rx_fp[0]  = Rx_fp[0]  * XFPr + XFPm;
    Rth_fp[0] = Rth_fp[0] * XpFPr + XpFPm;
    Ry_fp[0]  = Ry_fp[0]  * YFPr + YFPm;
    Rph_fp[0] = Rph_fp[0] * YpFPr + YpFPm;
    Lx_fp[0]  = Lx_fp[0]  * XFPr + XFPm;
    Lth_fp[0] = Lth_fp[0] * XpFPr + XpFPm;
    Ly_fp[0]  = Ly_fp[0]  * YFPr + YFPm;
    Lph_fp[0] = Lph_fp[0] * YpFPr + YpFPm;
    //===================================//
    
    Rz[0]   = Rz[0] * Ztr + Ztm;     // scaled
    Lz[0]   = Lz[0] * Ztr + Ztm;     // scaled    

}

//=================================================================================//

void momcalib::rasCorr(bool rarm, bool larm){


  //======== Nomalization =============//
  /*
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Rz[0]   = (Rz[0]-Ztm)/Ztr;

    Lx_fp[0]  = (Lx_fp[0]-XFPm)/XFPr; 
    Lth_fp[0] = (Lth_fp[0]-XpFPm)/XpFPr;
    Ly_fp[0]  = (Ly_fp[0]-YFPm)/YFPr;
    Lph_fp[0] = (Lph_fp[0]-YpFPm)/YpFPr;
    Lz[0]   = (Lz[0]-Ztm)/Ztr;
  */
  //===================================//    



    //======== Raster Correction ==========================//    

    if(rarm){
    RasterCor= Calc_ras(R_Ras_x, Pras[2], Pras[0]);
    RasterCor = RasterCor/tan(hrs_ang);
    //    Rz[0]  = Rz[0]*Ztr +Ztm; // scaled     
    Rz[0] = Rz[0] + RasterCor; // correction
    //    Rz[0]   = (Rz[0]-Ztm)/Ztr;    // nomalization
    }
    
    if(larm){
    RasterCor_L = Calc_ras(L_Ras_x, Pras_L[2], Pras_L[0]);
    RasterCor_L = RasterCor_L/tan(hrs_ang);
    //    Lz[0]   = Lz[0]*Ztr +Ztm;     // scaled
    Lz[0]   = Lz[0] + RasterCor_L;
    // cout<<"Lz "<<Lz[0]<<" RasL "<<RasterCor_L<<" Par2 "<<Pras_L[2]<<" Par0 "<<Pras_L[0]<<" rasx "<<L_Ras_x<<endl;
    //    Lz[0]    =  (Lz[0]  -  Ztm)/Ztr;    // nomalization
    }
    //====================================================//

    //========== Scaled ==================//
    /*
    Rx_fp[0]  = Rx_fp[0]  * XFPr + XFPm;
    Rth_fp[0] = Rth_fp[0] * XpFPr + XpFPm;
    Ry_fp[0]  = Ry_fp[0]  * YFPr + YFPm;
    Rph_fp[0] = Rph_fp[0] * YpFPr + YpFPm;

    Lx_fp[0]  = Lx_fp[0]  * XFPr + XFPm;
    Lth_fp[0] = Lth_fp[0] * XpFPr + XpFPm;
    Ly_fp[0]  = Ly_fp[0]  * YFPr + YFPm;
    Lph_fp[0] = Lph_fp[0] * YpFPr + YpFPm;

    Rz[0]  = Rz[0]*Ztr +Ztm; // scaled
    Lz[0]  = Lz[0]*Ztr +Ztm; // scaled

    */
    
    //===================================//

    
}

//===============================================================================//

void momcalib::angCorr(bool rth, bool rph,bool lth, bool lph ){


  //======== Nomalization =============//
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Rz[0]     = (Rz[0]-Ztm)/Ztr;
    Rth[0]    = (Rth[0] - Xptm)/Xptr;
    Rph[0]    = (Rph[0] - Yptm)/Yptr;    
    Lx_fp[0]  = (Lx_fp[0]-XFPm)/XFPr; 
    Lth_fp[0] = (Lth_fp[0]-XpFPm)/XpFPr;
    Ly_fp[0]  = (Ly_fp[0]-YFPm)/YFPr;
    Lph_fp[0] = (Lph_fp[0]-YpFPm)/YpFPr;
    Lz[0]     = (Lz[0]-Ztm)/Ztr;
    Lth[0]    = (Lth[0] - Xptm)/Xptr; 
    Lph[0]    = (Lph[0] - Yptm)/Yptr;     
  //===================================//    


    if(rth)Rth[0]  = calcf2t_ang(Pxpt,Rx_fp[0], Rth_fp[0],Ry_fp[0], Rph_fp[0],Rz[0]); // nomalized
    if(rph)Rph[0]  = calcf2t_ang(Pypt,Rx_fp[0], Rth_fp[0],Ry_fp[0], Rph_fp[0],Rz[0]); // nomalized    
    if(lth)Lth[0]  = calcf2t_ang(Pxpt_L, Lx_fp[0], Lth_fp[0],Ly_fp[0], Lph_fp[0], Lz[0]); // nomalized 
    if(lph)Lph[0]  = calcf2t_ang(Pypt_L, Lx_fp[0], Lth_fp[0],Ly_fp[0], Lph_fp[0], Lz[0]); // nomalized   


    //========== Scaled ==================//
    Rx_fp[0]  = Rx_fp[0]  * XFPr + XFPm;
    Rth_fp[0] = Rth_fp[0] * XpFPr + XpFPm;
    Ry_fp[0]  = Ry_fp[0]  * YFPr + YFPm;
    Rph_fp[0] = Rph_fp[0] * YpFPr + YpFPm;

    Lx_fp[0]  = Lx_fp[0]  * XFPr + XFPm;
    Lth_fp[0] = Lth_fp[0] * XpFPr + XpFPm;
    Ly_fp[0]  = Ly_fp[0]  * YFPr + YFPm;
    Lph_fp[0] = Lph_fp[0] * YpFPr + YpFPm;
    
    //===================================//

   
    Rph[0]  = Rph[0] * Yptr + Yptm; // scaled
    Rth[0]  = Rth[0] * Xptr + Xptm; // scaled
    Rz[0]   = Rz[0]  *  Ztr + Ztm; // scaled
    Lth[0]  = Lth[0] * Xptr + Xptm;  // scaled    
    Lph[0]  = Lph[0] * Yptr + Yptm;  // scaled    
    Lz[0]   = Lz[0]  * Ztr  + Ztm; // scaled
    
}

//================================================================================//

void momcalib::plossCorr(bool PLoss){

 
  
    if(PLoss){
     
    dpe=0.0; dpe_=0.0; dpk=0.0;
    dpe   = Eloss(0.0,0.0,"B");
    dpe_  = Eloss(Lph[0],Lz[0],"L");
    dpk   = Eloss(Rph[0],Rz[0],"R");
    pe    = pe    - dpe;
    Lp[0] = Lp[0] + dpe_;
    Rp[0] = Rp[0] + dpk;
    }// Energy loss


}

//===============================================================================//
void momcalib::momCorr(bool rarm, bool larm){

  Lp_x=-100.;
  Lp_y=-100.;
  Lp_z=-100.;
  Rp_x=-100.;
  Rp_y=-100.;
  Rp_z=-100.;

  
 
  //======== Nomalization =============//
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Lx_fp[0]  = (Lx_fp[0]-XFPm)/XFPr; 
    Lth_fp[0] = (Lth_fp[0]-XpFPm)/XpFPr;
    Ly_fp[0]  = (Ly_fp[0]-YFPm)/YFPr;
    Lph_fp[0] = (Lph_fp[0]-YpFPm)/YpFPr;

    Rz[0]   = (Rz[0]-Ztm)/Ztr;
    Lz[0]   = (Lz[0]-Ztm)/Ztr;
    
    Rp[0]   = (Rp[0]-PRm)/PRr;    
    Lp[0]   = (Lp[0]-PLm)/PLr;    

    //===================================//    

    
    if(rarm)Rp[0] = calcf2t_mom(Opt_par_R, Rx_fp[0], Rth_fp[0], Ry_fp[0], Rph_fp[0],Rz[0]);
    if(larm)Lp[0] = calcf2t_mom(Opt_par_L, Lx_fp[0], Lth_fp[0], Ly_fp[0], Lph_fp[0],Lz[0]);    

      
    //========== Scaled ==================//
    Rx_fp[0]  = Rx_fp[0]  * XFPr + XFPm;
    Rth_fp[0] = Rth_fp[0] * XpFPr + XpFPm;
    Ry_fp[0]  = Ry_fp[0]  * YFPr + YFPm;
    Rph_fp[0] = Rph_fp[0] * YpFPr + YpFPm;

    Lx_fp[0]  = Lx_fp[0]  * XFPr + XFPm;
    Lth_fp[0] = Lth_fp[0] * XpFPr + XpFPm;
    Ly_fp[0]  = Ly_fp[0]  * YFPr + YFPm;
    Lph_fp[0] = Lph_fp[0] * YpFPr + YpFPm;  

    //===================================//

    
    Rz[0]  = Rz[0] * Ztr + Ztm; // scaled
    Rp[0]  = Rp[0] * PRr + PRm; // scaled
    Lz[0]  = Lz[0] * Ztr + Ztm; // scaled
    Lp[0]  = Lp[0] * PLr + PLm; // scaled    


    //====== Scaled Lp =========//


    //    if(runnum>=111552){
    //      Lp[0]=2.21807/2.1*Lp[0];
    //    }


    if(runnum<=111220 || (111480<=runnum && 111542>=runnum))Lp[0]=Lp[0];
    else if(larm)Lp[0]=2.21807/2.1*Lp[0];

    //    plossCorr(ploss);
      plossCorr(ploss_f);
    
    //==== Right Hand Coorinate ====//

    Lp_z = Lp[0]/sqrt(1.0*1.0 + Lth[0]*Lth[0] + Lph[0]*Lph[0]);
    Rp_z = Rp[0]/sqrt(1.0*1.0 + Rth[0]*Rth[0] + Rph[0]*Rph[0]);

    Lp_x = Lp_z*     Lth[0] ;
    Lp_y = Lp_z*     Lph[0] ;
    Rp_x = Rp_z*     Rth[0] ;
    Rp_y = Rp_z*     Rph[0] ; 


    
    
}

//===============================================================================//

void momcalib::momCalc(){

  LP_x=-100.;
  LP_y=-100.;
  LP_z=-100.;
  RP_x=-100.;
  RP_y=-100.;
  Rp_z=-100.;

  //==== Right Hand Coorinate ====//

    LP_z = LP/sqrt(1.0*1.0 + Lth[0] * Lth[0] + Lph[0] * Lph[0] );
    RP_z = RP/sqrt(1.0*1.0 + Rth[0] * Rth[0] + Rph[0] * Rph[0] );


   
    /*
    LP_x = LP_z*    Lth[0] ;
    LP_y = LP_z*    Lph[0] ;
    RP_x = RP_z*    Rth[0] ;
    RP_y = RP_z*    Rph[0] ; 
    */

    LP_x = LP_z*   -Lth[0] ;
    LP_y = LP_z*   -Lph[0] ;
    RP_x = RP_z*   -Rth[0] ;
    RP_y = RP_z*   -Rph[0] ;     
    
}

//===============================================================================//

double momcalib::SSx(double z,double th){

  double ssx=-100.0;
      for(int j=0 ; j<nfoil ; j++){
	l[j]=0.0;
	l[j] = sqrt(pow(l0,2.0) + pow(fcent_real[j]*100.,2.0) -2.0*l0*fcent_real[j]*100.*cos(hrs_ang));
	dth[j] = asin(l0/l[j]*sin(hrs_ang)) - hrs_ang;
	projectf[j] = cos( dth[j] );
	if(fcent[j]-selection_width<z
	   && z<fcent[j]+selection_width){
	  ssx = - th * l[j] * projectf[j];
	}
      }
      return ssx;
}

//==============================================================================//

double momcalib::SSy(double z,double ph){

  double ssy=-100.0;
      for(int j=0 ; j<nfoil ; j++){
	l[j]=0.0;
	l[j] = sqrt(pow(l0,2.0) + pow(fcent_real[j]*100.,2.0) -2.0*l0*fcent_real[j]*100.*cos(hrs_ang));
	dth[j] = asin(l0/l[j]*sin(hrs_ang)) - hrs_ang;
	projectf[j] = cos( dth[j] );
	
	if(fcent[j]-selection_width<z
	   && z<fcent[j]+selection_width){
	  ssy = - ph * l[j] * projectf[j];
	}
      }
      return ssy;
}


//==============================================================================//

void momcalib::ParamCorr(){


    //===================================================================//
    //============================ Tuning value =========================//
    //===================================================================//



  
 
 // RHRS //

  

    //   
    //    hRth->Fill(Rth[0]);    
    //    hRph->Fill(Rph[0]);





  //======== Nomalization =============//
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Rth[0]    = (Rth[0] - Xptm)/Xptr;
    Rph[0]    = (Rph[0] - Yptm)/Yptr;    
    Rz[0]   = (Rz[0]-Ztm)/Ztr;
  //===================================//    



    
    if(MT[0])Rz[0] = calcf2t_zt(Pzt, Rx_fp[0], Rth_fp[0], Ry_fp[0], Rph_fp[0]);

    
    if(MT[1]){
    RasterCor= Calc_ras(R_Ras_x, Pras[2], Pras[0]);
    RasterCor = RasterCor/tan(hrs_ang);
    Rz[0]  = Rz[0]*Ztr +Ztm; // scaled     
    Rz[0] = Rz[0] - RasterCor; // correction    
    hRz_c->Fill(Rz[0]);    
    Rz[0]   = (Rz[0]-Ztm)/Ztr;    // nomalization
    }
    
    
    if(MT[2])Rth[0]  = calcf2t_ang(Pxpt,Rx_fp[0], Rth_fp[0],Ry_fp[0], Rph_fp[0],Rz[0]);
    if(MT[3])Rph[0] = calcf2t_ang(Pypt,Rx_fp[0], Rth_fp[0],Ry_fp[0], Rph_fp[0],Rz[0]);
    if(MT[8])Rp[0] = calcf2t_mom(Opt_par_R, Rx_fp[0], Rth_fp[0],Ry_fp[0], Rph_fp[0],Rz[0]);

    
    Rph[0]  = Rph[0]*Yptr +Yptm;
    Rth[0]  = Rth[0]*Xptr +Xptm;
    Rz[0]  = Rz[0]*Ztr +Ztm;


    Rx_fp[0]  = Rx_fp[0]  * XFPr + XFPm;
    Rth_fp[0] = Rth_fp[0] * XpFPr + XpFPm;
    Ry_fp[0]  = Ry_fp[0]  * YFPr + YFPm;
    Rph_fp[0] = Rph_fp[0] * YpFPr + YpFPm;


    //    hRth_c->Fill(Rth[0]);        
    //    hRph_c->Fill(Rph[0]);    


    
 // LHRS //



    //   
    //    hLth->Fill(Lth[0]);    
    //    hLph->Fill(Lph[0]);


    //=========  nomalization  ================//
    Lx_fp[0]  = (Lx_fp[0]-XFPm)/XFPr; 
    Lth_fp[0] = (Lth_fp[0]-XpFPm)/XpFPr;
    Ly_fp[0]  = (Ly_fp[0]-YFPm)/YFPr;
    Lph_fp[0] = (Lph_fp[0]-YpFPm)/YpFPr;
    Lz[0]    =  (Lz[0]  -  Ztm)/Ztr; 
    Lth[0]    = (Lth[0] - Xptm)/Xptr; 
    Lph[0]    = (Lph[0] - Yptm)/Yptr; 
    //========================================//
    
    


    if(MT[4])  Lz[0]   = calcf2t_zt(Pzt_L, Lx_fp[0], Lth_fp[0], Ly_fp[0], Lph_fp[0]);

    //======== Raster Correction ==========================//
    if(MT[5]){
    RasterCor_L = Calc_ras(L_Ras_x, Pras_L[2], Pras_L[0]);
    RasterCor_L = RasterCor_L/tan(hrs_ang);
    Lz[0]   = Lz[0]*Ztr +Ztm;     // scaled
    Lz[0] = Lz[0]+RasterCor_L;
    //    hLz_c->Fill(Lz[0]);    
    Lz[0]    =  (Lz[0]  -  Ztm)/Ztr;    // nomalization

    }
    //====================================================//
    
    if(MT[6])Lth[0]  = calcf2t_ang(Pxpt_L, Lx_fp[0], Lth_fp[0],Ly_fp[0], Lph_fp[0], Lz[0]);
    if(MT[7])Lph[0]  = calcf2t_ang(Pypt_L, Lx_fp[0], Lth_fp[0],Ly_fp[0], Lph_fp[0], Lz[0]);
    if(MT[9])Lp[0] = calcf2t_mom(Opt_par_L, Lx_fp[0], Lth_fp[0],Ly_fp[0], Lph_fp[0],Lz[0]);

    
    Lth[0]  = Lth[0]*Xptr +Xptm;  // scaled    
    Lph[0]  = Lph[0]*Yptr +Yptm;  // scaled
    Lz[0]   = Lz[0]*Ztr +Ztm;     // scaled

    Lx_fp[0]  = Lx_fp[0]  * XFPr + XFPm;
    Lth_fp[0] = Lth_fp[0] * XpFPr + XpFPm;
    Ly_fp[0]  = Ly_fp[0]  * YFPr + YFPm;
    Lph_fp[0] = Lph_fp[0] * YpFPr + YpFPm;



    //    hLth_c->Fill(Lth[0]);        
    //    hLph_c->Fill(Lph[0]);    


    
}


//================================================================================//
void momcalib::EventSelection(double ww){

  cout<<endl;

  cout<<"=========================================="<<endl;
  cout<<"========= Event Selection ================"<<endl;
  cout<<"=========================================="<<endl;  

  cout<<"Al weight : "<<ww<<endl;
  int nLam=0;
  int nSig=0;
  int nLamT=0;
  int nAl=0;
  int nnnL=0;
  int NeedTNum;
  int NeedTNum_par=0;
  cout<<"mode "<<MODE<<endl;
  if(MODE==0){NeedTNum=nParamTp*2;}
  else if(MODE==-1 || MODE==1){NeedTNum=nParamTp;}
  else{ cout<<"Error mode setting !! ("<<MODE<<")"<<endl;exit(1);}

   TNum=T->GetEntries();

   int break_num=nmax;
     //(int)NeedTNum*2;

   //   TNum=1000000;   
  cout<<"Select Events : "<<TNum<<endl;
  d=div(TNum,100);
  int nd=0;
  tuned_num=0;  
  for (int i=0;i<TNum;i++){
    
    // ===== Initialization ====== //    
    for(int j=0;j<Max;j++){
      Lp[j]=-2222.0;
      Rp[j]=-2222.0;
      Ls2tpads[j]=-2222.0;
      Rs2tpads[j]=-2222.0;
      rtrpathl[j]=-2222.0;
      rs2pathl[j]=-2222.0;
      Rx_fp[j]=-2222.0;
      Ry_fp[j]=-2222.0;
      Rth_fp[j]=-2222.0;
      Rph_fp[j]=-2222.0;
      Rx[j]=-2222.0;
      Ry[j]=-2222.0;
      Rz[j]=-2222.0;
      Rth[j]=-2222.0;
      Rph[j]=-2222.0;
      Lx_fp[j]=-2222.0;
      Ly_fp[j]=-2222.0;
      Lth_fp[j]=-2222.0;
      Lph_fp[j]=-2222.0;
      Lx[j]=-2222.0;
      Ly[j]=-2222.0;
      Lz[j]=-2222.0;
      Lth[j]=-2222.0;
      Lph[j]=-2222.0;
    }
    rssx=-2222.0;
    rssy=-2222.0;
    lssx=-2222.0;
    lssy=-2222.0;    
    rssx_c=-2222.0;
    rssy_c=-2222.0;
    lssx_c=-2222.0;
    lssy_c=-2222.0;    

    mm    =-2222.0;
    mm_L  =-2222.0;
    mm_pi =-2222.0;
    mm_Al =-2222.0;
    Ra1sum=-2222.0;
    Ra2sum=-2222.0;
    Rgssum=-2222.0;
    Lpssum=-2222.0;    
    coin_time=-2222.0;


    Lshsum=-2222.0;
    R_Ras_x=-222222.0;
    L_Ras_x=-222222.0;
    coin_t=-2222.0;
    runnum=0;
    //---- Events Flag -----//
    Rs2_flag=false;
    Ls2_flag=false;
    Rs0_flag=false;
    Ls0_flag=false;
    Scinti_flag=false;
    ac1_flag=false;
    ac2_flag=false;
    gs_flag=false;
    sh_flag=false;
    ps_flag=false;
    trig_flag=false;
    track_flag=false;
    Rpathl_flag=false;
    Lpathl_flag=false;
    coin_flag=false;
    x_flag=false;
    y_flag=false;
    z_flag=false;
    Al_flag=false;
    RPID_flag=false;
    LPID_flag=false;
    nnL_run=false;
    Lam_run=false;
    mom_select=false;
    T->GetEntry(i);


    
    
    pe=hallap/1000.;
    double Bp_b=pe;
    double Rp_b=Rp[0];
    double Lp_b=Lp[0];
    if(Ltrn!=1 || Rtrn!=1)continue;
    
    hLz_es  ->Fill(Lz[0]);
    hLth_es ->Fill(Lth[0]);
    hLph_es ->Fill(Lph[0]);    
    hLp_es  ->Fill(Lp[0]);

    lssx=SSx(Lz[0],Lth[0]);    
    lssy=SSy(Lz[0],Lph[0]);        
    hLss_es->Fill(lssy,lssx);
    rssx=SSx(Rz[0],Rth[0]);    
    rssy=SSy(Rz[0],Rph[0]);        
    hRss_es->Fill(rssy,rssx);
    hRz_es  ->Fill(Rz[0]);
    hRth_es ->Fill(Rth[0]);
    hRph_es ->Fill(Rph[0]);    
    hRp_es  ->Fill(Rp[0]);
    
    
 //===== Parameter Tuning ============//
    zCorr(MT[0],MT[4]);    
    hLz_es_c->Fill(Lz[0]);
    hRz_es_c->Fill(Rz[0]);    
    rasCorr(MT[1],MT[5]);
    hLz_es_rc->Fill(Lz[0]);
    hRz_es_rc->Fill(Rz[0]);

    angCorr(MT[2],MT[3],MT[6],MT[7]);
    hLth_es_c->Fill(Lth[0]);
    hRth_es_c->Fill(Rth[0]);
    hLph_es_c->Fill(Lph[0]);
    hRph_es_c->Fill(Rph[0]);
    lssx_c=SSx(Lz[0],Lth[0]);    
    lssy_c=SSy(Lz[0],Lph[0]);        
    hLss_es_c->Fill(lssy_c,lssx_c);
    rssx_c=SSx(Rz[0],Rth[0]);    
    rssy_c=SSy(Rz[0],Rph[0]);        
    hRss_es_c->Fill(rssy_c,rssx_c);



    momCorr(MT[8],MT[9]);
    hRp_es_c->Fill(Rp[0]);
    hLp_es_c->Fill(Lp[0]);
    
    
 //====================================//





    
 //------ Set Physics value --------//
    Ee=sqrt(pow(pe,2)+pow(Me,2));
    Ee_=sqrt(pow(Lp[0],2)+pow(Me,2));
    Ek=sqrt(pow(Rp[0],2)+pow(MK,2));
    Ls2pads=(int)Ls2tpads[0];
    Rs2pads=(int)Rs2tpads[0];
    rpathl=rtrpathl[0]+rs2pathl[0]; // R-HRS path length S2 -RF
    lpathl=ltrpathl[0]+ls2pathl[0]; // L-HRS path length S2 -RF
    rbeta=Rp[0]/Ek; 
    rpath_corr=rpathl/rbeta/c;
    lbeta=Lp[0]/Ee_; 
    lpath_corr=lpathl/lbeta/c;


    
    TVector3 B_v(0.0,0.0,sqrt(Ee*Ee-Me*Me));
    TVector3 R_v(Rp_x, Rp_y, Rp_z);
    TVector3 L_v(Lp_x, Lp_y, Lp_z);

    R_v.RotateX(   13.2/180.*3.14);
    L_v.RotateX( - 13.2/180.*3.14);

        

    mm = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
	       - (B_v - L_v - R_v)*(B_v - L_v - R_v) );

    mm_Al = sqrt( (Ee + MAl - Ee_ - Ek)*(Ee + MAl - Ee_ - Ek)
		  - (B_v - L_v - R_v)*(B_v - L_v - R_v) );

    mm_nnL = sqrt( (Ee + MTr - Ee_ - Ek)*(Ee + MTr - Ee_ - Ek)
	       - (B_v - L_v - R_v)*(B_v - L_v - R_v) );



   //--- Set Coincidence time ---------// 

    if(runnum<=111368){tdc_time=0.056; tdc_mode=1;}
    else {tdc_time=0.058; tdc_mode=2;}

 Rs2_off= s2f1_off(Rs2pads,"R",target,tdc_mode);
 Ls2_off= s2f1_off(Ls2pads,"L",target,tdc_mode);
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 coin_t=tof_r-tof_l+rpath_corr-lpath_corr-pathl_off-coin_offset; //  coin Path & Offset  correction   

 if(single)coin_t=coin_time;

 //---- Event Selection ----------//
  if(coin_t<1.0 && -1.0<coin_t)coin_flag=true;
  //  if((Rvz_cutmin<Rz[0] && Rz[0]<Rvz_cutmax) &&  ( Lvz_cutmin<Lz[0] && Lz[0]<Lvz_cutmax)) z_flag=true;
  if(fabs(Rz[0]-Lz[0])<0.025 && Rz[0]>-100 && Lz[0]>-100 && fabs(Rz[0]+Lz[0])/2.0 < 0.1 ) z_flag=true;
  //  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads)
  track_flag=true;
  if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax) Rpathl_flag=true;
  if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax) Lpathl_flag=true;
  if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0) Rs0_flag=true;
  if(-LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0) Ls0_flag=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)   Rs2_flag=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)      Ls2_flag=true;
  if(DRevtype==5)                              trig_flag=true;
  //  if(Ra1sum_p<a1_th)      ac1_flag=true;
  //  if(Ra2sum_p>a2_th)      ac2_flag=true;
  //===== AC NPE SUM CUT =====//
  if(AC1_npe_sum<a1_th)      ac1_flag=true;
  if(AC2_npe_sum>a2_th)      ac2_flag=true;
  
  //  if(fabs(Rz[0]-Lz[0])<0.025 && Rz[0]>-100 && Lz[0]>-100 && fabs(Rz[0]+Lz[0])/2.0 > 0.1) Al_flag=true;

  if(fabs(Rz[0]-Lz[0])<0.025 && Rz[0]>-100 && Lz[0]>-100 && (Rz[0]+Lz[0])/2.0 > 0.1) Al_flag=true; // test z >0 tuning
  //if(fabs(Rz[0]-Lz[0])<0.025 && Rz[0]>-100 && Lz[0]>-100 && (Rz[0]+Lz[0])/2.0 <-0.1) Al_flag=true; // test z <0 tuning

  gs_flag=true;
  ps_flag=true;
  sh_flag=true;
  
  if(ac1_flag && ac2_flag && gs_flag) RPID_flag=true;
  if(Lcersum>2000.)LPID_flag=true;  
  if(Rs2_flag && Rs0_flag && Ls2_flag && Ls0_flag)      Scinti_flag=true;
  Scinti_flag=true;
  
  if((runnum>111220 && runnum<111480) || (runnum>111576 && runnum<111698) 
     || (runnum>111738 && runnum<111830) ) nnL_run=true;
  else Lam_run=true;
  
  //  cout<<"runnum "<<runnum<<" Lam_run "<<Lam_run<<" nnL_run "<<nnL_run<<endl;


     //      Lp[0]=2.21807/2.1*Lp[0];
  
    // ------------------------------------------- //
    // ------- Event selection for tuning -------- //
    // ------------------------------------------- //







  
  if(trig_flag && Scinti_flag && track_flag && Rpathl_flag && Lpathl_flag){


    if(RPID_flag && LPID_flag && z_flag){
	hLz_cut->Fill(Lz[0]);
	hRz_cut->Fill(Rz[0]);}  
    

 
      if(trig_flag && RPID_flag && LPID_flag && z_flag){
        hct->Fill(coin_t);
	if(coin_flag)hct_cut->Fill(coin_t);
      }



      //======== Production Events ============//


      //====== Al Events Selection =============//
      bool zp_flag= false;
      if(fabs(Rz[0]-Lz[0])<0.025 && Rz[0]>-100 && Lz[0]>-100 && (Rz[0]+Lz[0])/2.0 > 0.1) zp_flag=true; // z >0
      bool zn_flag= false;
      if(fabs(Rz[0]-Lz[0])<0.025 && Rz[0]>-100 && Lz[0]>-100 && (Rz[0]+Lz[0])/2.0 < -0.1) zn_flag=true; // z <0      



      ////    z >0 Alminum events ////
     if(trig_flag && coin_flag && RPID_flag && LPID_flag  && zp_flag){
	hmm_Al_zp->Fill((mm_Al-MMgL)*1000.); 
      }

      ////    z <0 Alminum events ////
     if(trig_flag && coin_flag && RPID_flag && LPID_flag && zn_flag){
	hmm_Al_zn->Fill((mm_Al-MMgL)*1000.); 
      }
      
      
      
      if( trig_flag &&  RPID_flag && LPID_flag && Al_flag &&( (-45<coin_t && coin_t <-15) || (15<coin_t && coin_t<45) )   )
	hmm_Al_acc->Fill((mm_Al-MMgL)*1000.);
      
    //=====================================//
      
      //      Al_check =false;
      //      if(i%2==0)Al_check =true;
      Al_check =true;


      
      if(trig_flag && coin_flag && RPID_flag && LPID_flag && Al_flag && Al_check){
      hz_Al->Fill((Rz[0]+Lz[0])/2.0);
      hmm_Al->Fill((mm_Al-MMgL)*1000.);
      if(mom_select)      hmm_Al_mom->Fill((mm_Al-MMgL)*1000);
      bool MgL_flag = false;
      //      bool MgL_flag2= false;
      //double range_MgL=0.002; 

      //      double MgL_mean[4]={-4.0e-3, 6.0e-3, 20.0e-3, 30.0e-3}; //Bishnu parameters
      //      double MgL_mean[4]={0.00,0.07,0.18,0.2};
      //      double MgL_mean[4]={-100,-100,-100,0.168}; // GeV
      //      double MgL_mean[4]={-100,-100,-100,0.0475}; // GeV
      //      double MgL_mean[4]={-100,-100,-100,0.01}; // GeV
      //      double MgL_mean[4]={0.018,0.010,0.025,-100}; // GeV
      //double MgL_mean[4]={0.015,-100,-100,-100}; // GeV
      // double MgL_mean[4]={0.025,-100,-100,-100}; // GeV
      //      double MgL_mean[4]={-0.007,-100,-100,-100}; // GeV
      //      double MgL_mean[4]={-100,-0.020,-100,-100}; // Al_posi
      //      double MgL_mean[4]={-100., -17.6e-3, -7.0e-3, 1.8e-3}; //28 Al L peak
      //      double MgL_mean[4]={-100., -100, -7.0e-3, 100}; //28 Al L peak
      //     double MgL_mean[4]={-100., -100, -17.60e-3, 100}; //28 Al L peak
      //      double MgL_mean[4]={-100., -17.6e-3, -100.,-100.}; //28 Al L peak
      //      double MgL_mean[4]={-100., -100, -7.0e-3,-100.}; //28 Al L peak
      /*
      double MgL_mean[4]={-100., -100, 1.8e-3,-100.}; //28 Al L peak
      int MgL_peak=-1;
      
      for(int m=0;m<4;m++)
	if(fabs(mm_Al-MMgL -MgL_mean[m])<range_MgL){
	  MgL_flag=true;
	  MgL_peak=m;
	  break;
	}
*/
      if(fabs(mm_Al-MMgL -MgL_mean)<MgL_width)MgL_flag =true;
      //      if(fabs(mm_Al-MMgL)+Al_th<0.01 && Al_mode)MgL_flag=true;
      /*
      if( -range_MgL<mm_Al-MMgL - Al_th && mm_Al-MMgL - Al_th < range_MgL  && Al_mode)MgL_flag=true;  
      if( -range_MgL<mm_Al-MMgL - Al_th2 &&mm _Al-MMgL - Al_th2 < range_MgL  && Al_mode)MgL_flag2=true;  
      */

      //    if(MgL_flag && Al_mode){

    if(MgL_flag && Al_mode ){

      mass_flag[ntune_event]=2;
      //      mass_ref[ntune_event] = MMgL + MgL_mean[MgL_peak];
      mass_ref[ntune_event] = MMgL + MgL_mean;
      chi_Al+=(mm_Al-mass_ref[ntune_event])*(mm_Al-mass_ref[ntune_event])/0.002/0.002;

      //      cout<<"mm_MgL "<<mm_Al<<" MgL_mean "<<MgL_mean[MgL_peak]*1000.<<" mass_ref "<<mass_ref[ntune_event]<<" peak num "<<MgL_peak<<" diff "<<mm_Al-MMgL-MgL_mean[MgL_peak]<<endl;
      //      if(MgL_flag)mass_ref[ntune_event] = MMgL+Al_th;
      //      if(MgL_flag2)mass_ref[ntune_event] = MMgL+Al_th2;
      MM[ntune_event] = mm_Al;
      nrun[ntune_event] = runnum;
      rx_fp[ntune_event]  = Rx_fp[0];
      rth_fp[ntune_event] = Rth_fp[0];
      ry_fp[ntune_event]  = Ry_fp[0];
      rph_fp[ntune_event] = Rph_fp[0];
      rx[ntune_event]  = Rx[0];
      rth[ntune_event] = Rth[0];
      ry[ntune_event]  = Ry[0];
      rph[ntune_event] = Rph[0];
      rz[ntune_event] = Rz[0];      
      lx_fp[ntune_event]  = Lx_fp[0];
      lth_fp[ntune_event] = Lth_fp[0];
      ly_fp[ntune_event]  = Ly_fp[0];
      lph_fp[ntune_event] = Lph_fp[0];
      lx[ntune_event]  = Lx[0];
      lth[ntune_event] = Lth[0];
      ly[ntune_event]  = Ly[0];
      lph[ntune_event] = Lph[0];
      lz[ntune_event] = Lz[0];
      Rras_x[ntune_event] =R_Ras_x;
      Lras_x[ntune_event] =L_Ras_x;

      drp[ntune_event] =dpk;
      dlp[ntune_event] =dpe_;      
      dbp[ntune_event] =dpe;

      beam_p[ntune_event]=Bp_b;
      rp[ntune_event] =Rp_b;
      lp[ntune_event] =Lp_b;      
      bp[ntune_event] =Bp_b;

      rpc[ntune_event] =Rp[0];
      lpc[ntune_event] =Lp[0];      
      bpc[ntune_event] =pe;

      //      hz_Al->Fill((Rz[0]+Lz[0])/2.0);
      hmm_Al_select->Fill((mm_Al-MMgL)*1000.);
      tuned_num=i;
      nAl++;

      tevent[ntune_event]=i;
      ntune_event++;



      }
  

      }





  //======== ACCidencel B.G. ==============///
      if(nnL_run && trig_flag && RPID_flag && LPID_flag && z_flag &&( (-45<coin_t && coin_t <-15) || (15<coin_t && coin_t<45)))
	hmm_nnL_acc->Fill((mm_nnL-MnnL)*1000.);



    if(trig_flag && coin_flag && RPID_flag && LPID_flag && z_flag){
      mm_L=mm;
      if(Lam_run)      hmm_select->Fill((mm-ML)*1000.);
      if(nnL_run)      hmm_nnL->Fill((mm_nnL-MnnL)*1000.);



     
      //====== Tuning Events Selction =========//
      bool Lam_flag=false;
      bool Sig_flag=false;
      bool nnL_flag=false;


      if(Initial){
	if(1.09 < mm && mm < 1.115 && ntune_event<nmax ) Lam_flag=true;
	if(1.16 < mm && mm < 1.20 && ntune_event<nmax ) Sig_flag=true;
	//      if(1.09 < mm && mm < 1.12 && ntune_event<nmax ) Lam_flag=true;
	//      if(1.16 < mm && mm < 1.20 && ntune_event<nmax ) Sig_flag=true;
      
      }else if(form_mode){

	if(ML-mmL_range + paramL[3] + 0.002< mm && mm < ML+mmL_range  + paramL[3]+0.002 && ntune_event<nmax && ( ( MT[8] && MT[9] ) || single) )   Lam_flag=true;
	if(MS0-mmS_range + paramS[3] +0.002 < mm && mm < MS0 + mmS_range + paramS[3] +0.002 && ntune_event<nmax &&  ( ( MT[8] && MT[9] ) || single) ) Sig_flag=true;


      }else{
	
	if(ML-mmL_range < mm && mm < ML+mmL_range && ntune_event<nmax && ( ( MT[8] && MT[9] ) || single) )   Lam_flag=true;
	if(1.12 < mm && mm < 1.15 && ntune_event<nmax && ( MT[8]==0 || MT[9]==0 ) && single==0 )         Lam_flag=true;
	if(MS0-mmS_range < mm && mm < MS0+mmS_range && ntune_event<nmax &&  ( ( MT[8] && MT[9] ) || single) ) Sig_flag=true;
	if(1.20 < mm && mm < 1.23 && ntune_event<nmax && ( MT[8]==0 ||  MT[9]==0 ) && single==0)        Sig_flag=true;
      }
   

      if(nnL_run && ( fabs(mm_nnL-MnnL)<0.003  )&& nnL_mode)nnL_flag=true;


 

      // as a test //
      //      Lam_flag=false;
      //      Sig_flag=false;


      if(Lam_flag && Lam_run){

      nrun[ntune_event] = runnum;
      mass_flag[ntune_event]=0;
      mass_ref[ntune_event] = ML;
      MM[ntune_event]=mm;
      rx_fp[ntune_event]  = Rx_fp[0];
      rth_fp[ntune_event] = Rth_fp[0];
      ry_fp[ntune_event]  = Ry_fp[0];
      rph_fp[ntune_event] = Rph_fp[0];
      rx[ntune_event]  = Rx[0];
      rth[ntune_event] = Rth[0];
      ry[ntune_event]  = Ry[0];
      rph[ntune_event] = Rph[0];
      rz[ntune_event] = Rz[0];      
      lx_fp[ntune_event]  = Lx_fp[0];
      lth_fp[ntune_event] = Lth_fp[0];
      ly_fp[ntune_event]  = Ly_fp[0];
      lph_fp[ntune_event] = Lph_fp[0];
      lx[ntune_event]  = Lx[0];
      lth[ntune_event] = Lth[0];
      ly[ntune_event]  = Ly[0];
      lph[ntune_event] = Lph[0];
      lz[ntune_event] = Lz[0];
      Rras_x[ntune_event] =R_Ras_x;
      Lras_x[ntune_event] =L_Ras_x;

      drp[ntune_event] =dpk;
      dlp[ntune_event] =dpe_;      
      dbp[ntune_event] =dpe;

      beam_p[ntune_event]=Bp_b;
      rp[ntune_event] =Rp_b;
      lp[ntune_event] =Lp_b;      
      bp[ntune_event] =Bp_b;

      rpc[ntune_event] =Rp[0];
      lpc[ntune_event] =Lp[0];      
      bpc[ntune_event] =pe;

      hmm_L->Fill((mm-ML)*1000.);
      tuned_num=i;
      if(runnum<111555) nLam++;
      else nLamT++;

      chi_L+=(mm-ML)*(mm-ML)/0.002/0.002;

      tevent[ntune_event]=i;
      ntune_event++;
      }
      

      if(Sig_flag && Lam_run){

      nrun[ntune_event] = runnum;
      mass_flag[ntune_event]=1;
      mass_ref[ntune_event] = MS0;
      //      beam_p[ntune_event]=pe;

      MM[ntune_event]=mm;
      rx_fp[ntune_event]  = Rx_fp[0];
      rth_fp[ntune_event] = Rth_fp[0];
      ry_fp[ntune_event]  = Ry_fp[0];
      rph_fp[ntune_event] = Rph_fp[0];
      rx[ntune_event]  = Rx[0];
      rth[ntune_event] = Rth[0];
      ry[ntune_event]  = Ry[0];
      rph[ntune_event] = Rph[0];
      rz[ntune_event] = Rz[0];      
      lx_fp[ntune_event]  = Lx_fp[0];
      lth_fp[ntune_event] = Lth_fp[0];
      ly_fp[ntune_event]  = Ly_fp[0];
      lph_fp[ntune_event] = Lph_fp[0];
      lx[ntune_event]  = Lx[0];
      lth[ntune_event] = Lth[0];
      ly[ntune_event]  = Ly[0];
      lph[ntune_event] = Lph[0];
      lz[ntune_event] = Lz[0];
      Rras_x[ntune_event] =R_Ras_x;
      Lras_x[ntune_event] =L_Ras_x;            
      //      rp[ntune_event] =Rp[0];
      //      lp[ntune_event] =Lp[0];      
      beam_p[ntune_event]=Bp_b;
      rp[ntune_event] =Rp_b;
      lp[ntune_event] =Lp_b;      
      bp[ntune_event] =Bp_b;

      drp[ntune_event] =dpk;
      dlp[ntune_event] =dpe_;      
      dbp[ntune_event] =dpe;

      rpc[ntune_event] =Rp[0];
      lpc[ntune_event] =Lp[0];      
      bpc[ntune_event] =pe;
      hmm_S->Fill((mm-ML)*1000.);
      chi_S+=(mm-MS0)*(mm-MS0)/0.002/0.002;
      nSig++;
      tevent[ntune_event]=i;
      ntune_event++;
      tuned_num=i;
      }



      if(nnL_flag && nnL_run){

      nrun[ntune_event] = runnum;
      mass_flag[ntune_event]=3;
      mass_ref[ntune_event] = MnnL;
      MM[ntune_event]=mm_nnL;
      rx_fp[ntune_event]  = Rx_fp[0];
      rth_fp[ntune_event] = Rth_fp[0];
      ry_fp[ntune_event]  = Ry_fp[0];
      rph_fp[ntune_event] = Rph_fp[0];
      rx[ntune_event]  = Rx[0];
      rth[ntune_event] = Rth[0];
      ry[ntune_event]  = Ry[0];
      rph[ntune_event] = Rph[0];
      rz[ntune_event] = Rz[0];      
      lx_fp[ntune_event]  = Lx_fp[0];
      lth_fp[ntune_event] = Lth_fp[0];
      ly_fp[ntune_event]  = Ly_fp[0];
      lph_fp[ntune_event] = Lph_fp[0];
      lx[ntune_event]  = Lx[0];
      lth[ntune_event] = Lth[0];
      ly[ntune_event]  = Ly[0];
      lph[ntune_event] = Lph[0];
      lz[ntune_event] = Lz[0];
      Rras_x[ntune_event] =R_Ras_x;
      Lras_x[ntune_event] =L_Ras_x;            
      beam_p[ntune_event]=Bp_b;
      rp[ntune_event] =Rp_b;
      lp[ntune_event] =Lp_b;      
      bp[ntune_event] =Bp_b;

      drp[ntune_event] =dpk;
      dlp[ntune_event] =dpe_;      
      dbp[ntune_event] =dpe;

      rpc[ntune_event] =Rp[0];
      lpc[ntune_event] =Lp[0];      
      bpc[ntune_event] =pe;
      hmm_nnL_select->Fill((mm_nnL-MnnL)*1000.);
      nnnL++;
      tevent[ntune_event]=i;
      ntune_event++;
      tuned_num=i;
      }



      
      //      if(runnum<=111220 || (111480<=runnum && 111542>=runnum)) scale[ntune_event-1]=false; // pe=2.1 GeV
      //      else scale[ntune_event-1]=true;  // pe_ = 2.2 GeV

    }
  } // event selection

  //=============== COMMENT OUT ========================================//
  if(i==d.quot * 10*nd){cout<<10*nd<<" % Filled ("<<i<<" / "<<TNum<<")"<<endl; nd++;}
  if(ntune_event==NeedTNum/10*NeedTNum_par){cout<<"          "<<10*NeedTNum_par<<" % Tuning Events "<<ntune_event<<" / "<<NeedTNum<<endl;NeedTNum_par++;}
  if(ntune_event>=break_num || ntune_event==nmax){cout<<"Get enough tuning events "<<ntune_event<<endl;tuned_num=i;break;}
  //====================================================================//

  }//End Fill


  weight  = (double)nLam/(double)nSig;
  weightT = (double)nLam/(double)nLamT;
  weightAl= (double)nLam/ (double)nAl;
  weightnnL= (double)nLam/ (double)nnnL;
  if(nSig==0)  weight  = 1.0;
  if(nLamT==0) weightT = 1.0;
  if(nAl==0 || nLam==0)   weightAl= 1.0;
  if(nnnL==0)   weightnnL= 1.0;
  weightAl=weightAl*MgL_weight;
  weightAl_def=weightAl;
  weight_def = weight;
  Al_range = range_MgL;
  //  weight *=2.0;
  // chi square //
  chi_Al = chi_Al/nAl;
  chi_L  = chi_L/nLam;
  chi_S  = chi_S/nSig;

  chi2_init= chi_Al+chi_L+chi_S;
  


  //  weightnnL=10.;
  //  weightAl=0.0;
  
  cout<<"Tuning Events: "<<tuned_num<<endl;
  cout<<"Select Events: "<<ntune_event<<endl;
  cout<<"Select Lambda events : "<<nLam<<endl;
  cout<<"Select Sigma  events : "<<nSig<<endl;
  cout<<"Select Lam(T) events : "<<nLamT<<endl;
  cout<<"Select Al     events : "<<nAl<<endl;
  cout<<"Select nnL    events : "<<nnnL<<endl;
  cout<<"Weight Al : "<<MgL_weight<<endl;
  cout<<"Ratio Lam/Sig :"<<weight<<endl;
  cout<<"Ratio LamH/LamT :"<<weightT<<endl;
  cout<<"Ratio LamH/Al :"<<weightAl<<endl;
  cout<<"Ratio LamH/T :"<<weightnnL<<endl;
  cout<<" Chi2 Lam    :"<<chi_L<<endl;
  cout<<" Chi2 Sigma  :"<<chi_S<<endl;
  cout<<" Chi2 Al     :"<<chi_Al<<endl;
  cout<<" Chi2 sum    :"<<chi2_init<<endl;
}

//========================================================//




void momcalib::MomTuning(string ofname){


  if (nite>0){
    cout << "======================================================" <<endl;
    cout << "=================  Tuning started ====================" << endl;
    cout << "======================================================" <<endl;}
    char tempc[500],tempc2[500];
    const  char* new_tempc=ofname.c_str();
    cout<<"new marix file: "<<new_tempc<<endl;
    
    ofstream * ofs1;
    ofstream * ofs2;  

    //====== Set mean & sigma ========//

    //    /*
    //    TF1* fmean = new TF1("fmean","expgaus_mean",-0.300,0.300,4);
    //    fmean->SetParameters(parL[0],parL[1],parL[2],parL[3]);
    //    Lam_mean = fmean-> Integral(-0.05,0.1);
    //    TF1* fsigma = new TF1("fsigma","expgaus_sigma",-0.300,0.300,5);
    //    fsigma->SetParameters(parL[0],parL[1],parL[2],parL[3],Lam_mean);
    //    Lam_sigma = fsigma-> Integral(-0.05,0.1);

    //    */

    cout<<" Lambda  mean : "<<Lam_mean*1000.<< " MeV"<<endl;
    cout<<" Lambda  sigma : "<<Lam_sigma*1000<<" MeV"<<endl;


    for(int i=0 ; i<nite ; i++){

    cout<<"tuning i: "<<i+1<<" /"<<nite<<endl;
    // --------------------------- //
    // ---- Parameter tuning ----- //
    // --------------------------- //
    
    if(MODE==-1){
    cout<<"------- Rp tuning -----"<<endl;
    chi_sq1[i]=0.0;
    chi_sq1[i] = tune(Opt_par_R,i,-1);   // Rp
    }else if(MODE==1){
    cout<<"------- Lp tuning -----"<<endl;    
    chi_sq2[i]=0.0;
    chi_sq2[i] = tune(Opt_par_L,i,1);  // Lp
    }else if(MODE==0){
    cout<<"---------pk & pe tuning ------"<<endl;
    chi_sq[i]=0.0;
    weightAl=weightAl_def*(i);
    //    Al_range = range_MgL - 0.001*i;
    //    if(Al_range<=0)range_MgL - 0.001*(i-1) -0.0001*i;
    //    weight=weight_def*(i+1);
    cout<<"weightAl "<<weightAl<<endl;
    cout<<"weight "<<weight<<endl;
    chi_sq[i] = tune(Opt_par,i,0);

    }

    cout << " Tuning# = " << i+1 << ": chisq = ";

    if(MODE==-1){
      cout  << chi_sq1[i] << endl;
      gchi_Rp->SetPoint(i,i,chi_sq1[i]);          
    }else if(MODE== 1){
      cout  << chi_sq2[i] << endl;
      gchi_Lp->SetPoint(i,i,chi_sq2[i]);     
    }else if(MODE== 0){
      cout  << chi_sq[i]  << endl;
      gchi_p ->SetPoint(i,i,chi_sq[i]);    
      if(Al_mode && i==0)gchi_w ->SetPoint(0,0,chi2_init);
      if(Al_mode)gchi_w ->SetPoint(i+1,weightAl,chi_sq[i]);   

    }

    
    if(MODE==-1 || MODE==0){
    sprintf(tempc,  "%s_Rp_%d_MODE%d.dat",new_tempc,i,MODE); 
    cout<<"new matrix Rp: "<<tempc<<endl;
    ofs1 = new ofstream(tempc);}

    if(MODE==1 ||MODE==0){
    sprintf(tempc2, "%s_Lp_%d_MODE%d.dat",new_tempc,i,MODE);
    cout<<"new matrix Lp: "<<tempc2<<endl;    
    ofs2 = new ofstream(tempc2);}


    int nppp = 0;
    
    for(int i=0 ; i<nnp+1 ; i++){
      for(int e=0 ; e<nnp+1 ; e++){
	for(int d=0 ; d<nnp+1 ; d++){
	  for(int c=0 ; c<nnp+1 ; c++){
	    for(int b=0 ; b<nnp+1 ; b++){
	      for(int a=0 ; a<nnp+1 ; a++){  
		if(a+b+c+d+e==i){
		  if(MODE==-1){		  
		  *ofs1 << Opt_par_R[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  }else if(MODE==1){
		  *ofs2 << Opt_par_L[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  }else if(MODE==0){

		    *ofs1 << Opt_par[nppp] 
			  << " " << a 
			  << " " << b
			  << " " << c
			  << " " << d
			  << " " << e << endl;
		    *ofs2 << Opt_par[nParamTp+nppp] 
			  << " " << a 
			  << " " << b
			  << " " << c
			  << " " << d
			  << " " << e << endl;		    
		    

		    
		  }

		  nppp++;

		}
	      }
	    }
	  }
	}
      }
    }

    if(MODE==-1 || MODE==0){
      ofs1->close();
      ofs1->clear();
    }else if(MODE==0 || MODE==1){
      ofs2->close();
      ofs2->clear();}


    
  }

  
  cout<<"========== Tuning is done ============="<<endl;
  
  

}

//========================================================//

double momcalib::Eloss(double yp, double z,char* arm){


  
  double x;

  
  //  if(arm) x = - tan(hrs_ang-yp); //yp : phi [rad] RHRS
  //  else    x = - tan(hrs_ang+yp); //yp : phi [rad] LHRS

  //----- Original coordinate  -------//
  // Definition by K. Suzuki  (fixed Oct. 23rd, 2019)       //
  // R-HRS : right hand coordinate (Unticlockwise rotation) //
  // L-HRS : left  hand coordinate (    Clockwise rotation) //


  
  //  if(arm=="R")        x = - hrs_ang - yp; //yp : phi [rad] RHRS
  //  else if(arm=="L")   x = - hrs_ang + yp; //yp : phi [rad] LHRS
  //  if(arm=="R")        x = - hrs_ang - yp; //yp : phi [rad] RHRS
  //  else if(arm=="L")   x = - hrs_ang + yp; //yp : phi [rad] LHRS
  //  elase x=0.0;


  // Gogami Model //
  if(arm=="R")        x = - hrs_ang + yp; //yp : phi [rad] RHRS
  else if(arm=="L")   x = - hrs_ang - yp; //yp : phi [rad] LHRS
  else x=0.0;

  
  double ph[3],pl[2];
  double dEloss;
  bool high;
  
  if(z>0.08)high=false; //low : low  energy Loss (z> 8 cm)  pol1 function
  else high=true;       //high  : high  energy Loss (z< 8 cm) sin function
  
    //==== thickness 0.400 mm ========//

  if(arm=="R"){
    ph[0] = -1.31749;
    ph[1] = -4.61513;
    ph[2] = 2.03687;
    pl[0] = 3.158e-2;
    pl[1] = 4.05819e-1;
    
  }else if(arm=="L"){
    ph[0] = -1.3576;
    ph[1] = -4.5957;
    ph[2] = 2.0909;
    pl[0] = 6.2341e-3;
    pl[1] = 4.0336e-1;
    
  }

  double dEloss_h = ph[0]*sin(ph[1]*x)+ph[2];
  double dEloss_l = pl[0]*x +pl[1];

  
  if(high)dEloss = dEloss_h;
  else dEloss = dEloss_l;
  
  
  
  //==== thickness 0.4 mm in beam energy loss ======//
  if(arm=="B")dEloss=0.1843; //[MeV/c]



  //======== Al Flame Energy Loss ========//
  // Upstream && Downsream //
  // thickness 0.4 mm //
  if(z<-0.1){
    if(arm=="B")dEloss  = 0.1345;
    if(arm=="R")dEloss += 0.803;
    if(arm=="L")dEloss += 0.809;
  }
  else if(z>0.1){
    if(arm=="B")dEloss  = 0.264;
    if(arm=="R")dEloss += 0.313;
    if(arm=="L")dEloss += 0.313;
    
  }

  dEloss=dEloss/1000.; // [GeV/c]  
  return dEloss;
  
}

//========================================================//

void momcalib::Fill(){

  cout<<endl;
  cout<<"======================================"<<endl;
  cout<<"============= Event Fill ============="<<endl;
  cout<<"======================================"<<endl;  

  int t=0;

  ENum= T->GetEntries();
  cout<<"Events :"<<ENum<<endl;
     if(ENum<10000)d=div(ENum,1000);
   else   d=div(ENum,10000);

     bool fill_flag;
     bool fill_Al;
     bool fill_nnL;
     //     ENum=100000;
  for (int i=0;i<ENum;i++){

    // ===== Initialization ====== //    

    for(int j=0;j<Max;j++){
      Lp[j]=-2222.0;
      Rp[j]=-2222.0;
      Ls2tpads[j]=-2222.0;
      Rs2tpads[j]=-2222.0;
      rtrpathl[j]=-2222.0;
      rs2pathl[j]=-2222.0;
      Rx_fp[j]=-2222.0;
      Ry_fp[j]=-2222.0;
      Rth_fp[j]=-2222.0;
      Rph_fp[j]=-2222.0;
      Rx[j]=-2222.0;
      Ry[j]=-2222.0;
      Rth[j]=-2222.0;
      Rph[j]=-2222.0;
      Lx_fp[j]=-2222.0;
      Ly_fp[j]=-2222.0;
      Lth_fp[j]=-2222.0;
      Lph_fp[j]=-2222.0;
      Lx[j]=-2222.0;
      Ly[j]=-2222.0;
      Lth[j]=-2222.0;
      Lph[j]=-2222.0;
    }

    
    mm     =-2222.0;
    mm_L   =-2222.0;
    mm_L_pz=-2222.0;    
    mm_pi  =-2222.0;
    mm_acc =-2222.0;
    mm_Al  =-2222.0;
    MM_nnL = -2222.0;
    MM_Al  = -2222.0;
    mm_Al_acc =-2222.0;
    mm_Al_c = -2222.0;
    mm_Al_acc_c = -2222.0;
    MM_L=-2222.0;
    mm_X=-2222.0;
    MM_L_b=-2222.0;
    mm_nnL =-2222.0;
    mm_nnL_acc=-2222.0;
    mm_nnL_c=-2222.0;
    mm_nnL_acc_c=-2222.0;
    mm_Al_nnL=-2222.0;
    Ra1sum =-2222.0;
    Ra2sum =-2222.0;
    Rgssum =-2222.0;
    Lpssum =-2222.0;    
    Lshsum =-2222.0;
    coin_time=-2222.0;
    R_Ras_x=-100010.;
    L_Ras_x=-100000.;
    Dpe=-10;
    Dpe_=-10;
    Dpk=-10;
    RP=-100.;
    LP=-100.;
    mm_L_b=-2222.0;
    mm_L=-2222.0;
    //---- Events Flag -----//
    Rs2_flag=false;
    Ls2_flag=false;
    Rs0_flag=false;
    Ls0_flag=false;
    Scinti_flag=false;
    ac1_flag=false;
    ac2_flag=false;
    gs_flag=false;
    sh_flag=false;
    ps_flag=false;
    trig_flag=false;
    track_flag=false;
    Rpathl_flag=false;
    Lpathl_flag=false;
    coin_flag=false;
    x_flag=false;
    y_flag=false;
    z_flag=false;
    RPID_flag=false;
    LPID_flag=false;
    Al_flag=false;
    nnL_run=false;
    tflag=-1;
    tev=-1;
    fill_flag=false;

    T->GetEntry(i);

    RP=Rp[0];
    LP=Lp[0];
    if(Ltrn!=1 || Rtrn!=1)continue;
    
    
    if(tevent[t]==i){
      tev=1;
      t++;}
    
    if(i<tuned_num)tflag=1;
    else if(i>tuned_num)tflag=0;

    hLz  ->Fill(Lz[0]);
    hLth ->Fill(Lth[0]);
    hLph ->Fill(Lph[0]);    
    hLp  ->Fill(Lp[0]);

    lssx=SSx(Lz[0],Lth[0]);    
    lssy=SSy(Lz[0],Lph[0]);        
    hLss->Fill(lssy,lssx);
    rssx=SSx(Rz[0],Rth[0]);    
    rssy=SSy(Rz[0],Rph[0]);        
    hRss->Fill(rssy,rssx);

   
    hRz  ->Fill(Rz[0]);
    hRth ->Fill(Rth[0]);
    hRph ->Fill(Rph[0]);    
    hRp  ->Fill(Rp[0]);

    momCalc();
    pe=hallap/1000.;
    
 //===== Parameter Tuning ============//


    zCorr(MT_f[0],MT_f[4]);
    hLz_c->Fill(Lz[0]);
    hRz_c->Fill(Rz[0]);        

    rasCorr(MT_f[1],MT_f[5]);
    hLz_rc->Fill(Lz[0]);
    hRz_rc->Fill(Rz[0]);


    angCorr(MT_f[2],MT_f[3],MT_f[6],MT_f[7]);
    hLth_c->Fill(Lth[0]);
    hRth_c->Fill(Rth[0]);
    hLph_c->Fill(Lph[0]);
    hRph_c->Fill(Rph[0]);
    lssx_c=SSx(Lz[0],Lth[0]);    
    lssy_c=SSy(Lz[0],Lph[0]);        
    hLss_c->Fill(lssy_c,lssx_c);
    rssx_c=SSx(Rz[0],Rth[0]);    
    rssy_c=SSy(Rz[0],Rph[0]);        
    hRss_c->Fill(rssy_c,rssx_c);


    momCorr(MT_f[8],MT_f[9]);
    hRp_c->Fill(Rp[0]);
    hLp_c->Fill(Lp[0]);

    //    plossCorr(ploss_f);   
    //    hRp_ec->Fill(Rp[0]);
    //    hLp_ec->Fill(Lp[0]);

    //    Dpe = pe;
    //    Dpe_= Lp[0];
    //    Dpk = Rp[0];
    Dpe = dpe*1000.;  // MeV
    Dpe_= dpe_*1000.; //MeV      
    Dpk = dpk*1000.;  // MeV   
 //====================================//

    
 //------ Set Physics value --------//
    Ee=sqrt(pow(pe,2)+pow(Me,2));
    Ee_=sqrt(pow(Lp[0],2)+pow(Me,2));
    Ek=sqrt(pow(Rp[0],2)+pow(MK,2));
    Ls2pads=(int)Ls2tpads[0];
    Rs2pads=(int)Rs2tpads[0];
    rpathl=rtrpathl[0]+rs2pathl[0]; // R-HRS path length S2 -RF
    lpathl=ltrpathl[0]+ls2pathl[0]; // L-HRS path length S2 -RF
    rbeta=Rp[0]/Ek; 
    rpath_corr=rpathl/rbeta/c;
    lbeta=Lp[0]/Ee_; 
    lpath_corr=lpathl/lbeta/c;
    

    // before tuning parameters //


    Ee_2 = sqrt(pow(LP,2)+pow(Me,2));
    Ek2  = sqrt(pow(RP,2)+pow(MK,2));    
    TVector3 B_vb(0.0,0.0,sqrt(Ee*Ee-Me*Me));
    TVector3 R_vb(RP_x, RP_y, RP_z);
    TVector3 L_vb(LP_x, LP_y, LP_z);
    R_vb.RotateX(   13.2/180.*3.14);
    L_vb.RotateX( - 13.2/180.*3.14);


   MM_Al_b = sqrt( (Ee + MAl - Ee_2 - Ek2)*(Ee + MAl - Ee_2 - Ek2)
		 - (B_vb - L_vb - R_vb)*(B_vb - L_vb - R_vb) );
   MM_Al_b = (MM_Al_b - MMgL)*1000.;

   MM_nnL_b = sqrt( (Ee + MTr - Ee_2 - Ek2)*(Ee + MTr - Ee_2 - Ek2)
		  - (B_vb - L_vb - R_vb)*(B_vb - L_vb - R_vb) );
   
   MM_nnL_b = (MM_nnL_b - MnnL)*1000.;
   
   MM_L_b = sqrt( (Ee + Mp - Ee_2 - Ek2)*(Ee + Mp - Ee_2 - Ek2)
		 - (B_vb - L_vb - R_vb)*(B_vb - L_vb - R_vb) );
   MM_L_b = (MM_L_b - ML)*1000.;

    // after tuning parameters //
    
    TVector3 B_v(0.0,0.0,sqrt(Ee*Ee-Me*Me));
    TVector3 R_v(Rp_x, Rp_y, Rp_z);
    TVector3 L_v(Lp_x, Lp_y, Lp_z);
    R_v.RotateX(   13.2/180.*3.14);
    L_v.RotateX( - 13.2/180.*3.14);

    
    
    mm = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
	       - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
    
 
    TVector3 L_vz, R_vz, B_vz;
   
   double mm_pz;
   mm_pz = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
		 - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
   
   
   MM_Al = sqrt( (Ee + MAl - Ee_ - Ek)*(Ee + MAl - Ee_ - Ek)
		 - (B_v - L_v - R_v)*(B_v - L_v - R_v) );

   MM_Al = (MM_Al - MMgL)*1000.;

   MM_nnL = sqrt( (Ee + MTr - Ee_ - Ek)*(Ee + MTr - Ee_ - Ek)
		  - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
   
   MM_nnL = (MM_nnL - MnnL)*1000.;
   
   MM_L = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
		  - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
   
   MM_L = (MM_L - ML)*1000.;

   mm_X = sqrt( (Ee + 100. - Ee_ - Ek)*(Ee + 100. - Ee_ - Ek)
		  - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
   
   mm_X = (mm_X - 100.)*1000.;   


   //   if(runnum==111157 && nev==240177)cout<<"Ee "<<Ee<<" Ee_ "<<Ee_<<" Ek "<<Ek<<" A "<<(Ee + Mp - Ee_ - Ek)<<endl;
   
     //cout<<"A "<< (Ee + Mp - Ee_ - Ek)<<" B "<<sqrt((B_v - L_v - R_v)*(B_v - L_v - R_v))<<" mm "<<MM_L<<endl;
 //--- Set Coincidence time ---------// 
 // Rs2_off=RS2_off_H1[Rs2pads];
 // Ls2_off=LS2_off_H1[Ls2pads];

    if(runnum<=111368){
      tdc_time=0.056; 
      tdc_mode=1;
      coin_offset=464.73; // H1 mode
    }else {
      tdc_time=0.058; 
      tdc_mode=2;
      coin_offset=470.63; // H2 mode
}

 


 Rs2_off= s2f1_off(Rs2pads,"R",target,tdc_mode);
 Ls2_off= s2f1_off(Ls2pads,"L",target,tdc_mode);
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 coin_t=tof_r-tof_l+rpath_corr-lpath_corr-pathl_off-coin_offset; //  coin Path & Offset  correction

 // coin_time=coin_t;
 // if(single)
   coin_t=coin_time;

 
 //---- Event Selection ----------//
  if(coin_t<1.0 && -1.0<coin_t)coin_flag=true; //changed by itabashi 2020/06/22
  //  if(Rvz_cutmin<Rz[0] && Rz[0]<Rvz_cutmax && Lvz_cutmin<Lz[0] && Lz[0]<Lvz_cutmax) z_flag=true; 
  if(fabs(Rz[0]-Lz[0])<0.025 && Rz[0]>-100 && Lz[0]>-100 && fabs(Rz[0]+Lz[0])/2.0 < 0.1 ) z_flag=true;
  //  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads)   track_flag=true; 
  track_flag=true;  //changed by itabashi 2020/06/22
  //  if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax) Rpathl_flag=true;
  //  if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax) Lpathl_flag=true;
  
  Rpathl_flag=true; // as a test 
  Lpathl_flag=true; // as a test
  if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0) Rs0_flag=true;
  if(-LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0) Ls0_flag=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)   Rs2_flag=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)      Ls2_flag=true;
  if(DRevtype==5)                              trig_flag=true;
  //  if(Ra1sum_p<a1_th)      ac1_flag=true;
  //  if(Ra2sum_p>a2_th)      ac2_flag=true;
  //===== AC NPE SUM CUT =====//
  if(AC1_npe_sum<a1_th)      ac1_flag=true;
  if(AC2_npe_sum>a2_th)      ac2_flag=true;
  
  gs_flag=true;
  ps_flag=true;
  sh_flag=true;



  
  if(ac1_flag && ac2_flag && gs_flag) RPID_flag=true;
  if(Lcersum>2000.)LPID_flag=true;  
  //  if(Rs2_flag && Rs0_flag && Ls2_flag && Ls0_flag) Scinti_flag=true;
  Scinti_flag=true;
   
   if(trig_flag && Scinti_flag && track_flag && Rpathl_flag && Lpathl_flag){

     // if(trig_flag ){
     if(trig_flag && RPID_flag && LPID_flag && z_flag){hct2->Fill(coin_t);}
     if(fabs(Rz[0]-Lz[0])<0.025 && Rz[0]>-100 && Lz[0]>-100 && fabs(Rz[0]+Lz[0])/2.0 > 0.1) Al_flag=true;
     if((runnum>111220 && runnum<111480) || (runnum>111576 && runnum<111698) 
	|| (runnum>111738 && runnum<111830) ) nnL_run=true;


     
  //======== Production Events ============//

      Al_check =false;
      Al_check =true;

      
      bool zp_flag= false;
      if(fabs(Rz[0]-Lz[0])<0.025 && Rz[0]>-100 && Lz[0]>-100 && (Rz[0]+Lz[0])/2.0 > 0.1) zp_flag=true; // z >0
      bool zn_flag= false;
      if(fabs(Rz[0]-Lz[0])<0.025 && Rz[0]>-100 && Lz[0]>-100 && (Rz[0]+Lz[0])/2.0 < -0.1) zn_flag=true; // z <0
      
      ////    z >0 Alminum events ////
     if(trig_flag && coin_flag && RPID_flag && LPID_flag && Al_flag && zp_flag){
       hmm_Al_zp_c->Fill(MM_Al);
       hmm_Al_nnL_zp->Fill(MM_Al,MM_nnL);
     }
     
     ////    z <0 Alminum events ////
     if(trig_flag && coin_flag && RPID_flag && LPID_flag && Al_flag && zn_flag){
	hmm_Al_zn_c->Fill(MM_Al);
	hmm_Al_nnL_zn->Fill(MM_Al,MM_nnL);
     }
     




	
      
      if(trig_flag && coin_flag && RPID_flag && LPID_flag && Al_flag){
	if(Al_check)hmm_Al_cut->Fill(MM_Al);
	mm_Al = MM_Al_b;
	mm_Al_c = MM_Al;
	mm_Al_nnL = MM_nnL;
	    hmm_Al_nnL->Fill(MM_Al,MM_nnL);
	for(int j=0;j<ntune_event;j++)
	  if(mass_flag[j]==2 && tevent[j]==i){
	    hmm_Al_select_c->Fill(MM_Al);

	  }
	fill_flag=true;
      }

   
   
    if(coin_flag && RPID_flag && LPID_flag && z_flag){

      /*
  if(fabs(Rz[0]-Lz[0])<0.025 && Rz[0]>-100 && Lz[0]>-100 && fabs(Rz[0]+Lz[0])/2.0 < 0.1 ) z_flag=true;
  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads) track_flag=true;
  if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax) Rpathl_flag=true;
  if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax) Lpathl_flag=true;
  if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0) Rs0_flag=true;
  if(-LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0) Ls0_flag=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)   Rs2_flag=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)      Ls2_flag=true;
  if(DRevtype==5)                              trig_flag=true;a
       */
      
      fill_flag=true;
      if(nnL_run){hmm_nnL_cut->Fill(MM_nnL);
      mm_nnL = MM_nnL_b;
      mm_nnL_c = MM_nnL;
      
      }else { // hydrogen run
      mm_L=MM_L;
      mm_L_b=MM_L_b;
      mm_L_pz=mm_pz;
      hmm_cut->Fill((mm-ML)*1000.);

      if(runnum>=111552)hmm_cut_T->Fill(MM_L);
      else hmm_cut_H->Fill(MM_L);
      }
    }

  //======== ACCidencel B.G. ==============///
    if(trig_flag && RPID_flag && LPID_flag && z_flag &&
       ((-45<coin_t && coin_t <-15) || (15<coin_t && coin_t<45))){
      fill_flag=true;
      ct=coin_t;
      mm_acc=mm;
      mm_nnL_acc = MM_nnL_b;
      mm_nnL_acc_c = MM_nnL;
      if(nnL_run)hmm_nnL_cut_acc->Fill(MM_nnL);
              while(1){
	       if(-3.0<ct && ct<3.0){ctime=ct; break;}
	       else if(ct<-3.0) ct=ct+6;
	       else if(3.0<ct)  ct=ct-6;
	      }
    }


    if(trig_flag && RPID_flag && LPID_flag && Al_flag &&
       ( (-45<coin_t && coin_t <-15) || (15<coin_t && coin_t<45) )){
      fill_flag=true;
      hmm_Al_cut_acc->Fill(MM_Al);
      mm_Al_acc = MM_Al_b;
      mm_Al_acc_c = MM_Al;
    }



    //=====================================//

    
   }


  
   //   if(Initial && fill_flag)
   

   
   if(fill_flag)   tnew->Fill();
   if(i % (d.quot * 1000) == 0)cout<<i<<" / "<<ENum<<endl;
  }

  
}

//========================================================//

double* momcalib::MixiedEvent(double * coin_acc){


  


}
//=======================================================//

void momcalib::Close(){

   tnew->Write();


   hmm_Al_cut_acc->Scale(2./60.);
   hmm_nnL_cut_acc->Scale(2./60.);
   hmm_Al_acc->Scale(2./60.);
   hmm_nnL_acc->Scale(2./60.);

   hLz->Write();
   hz_Al->Write();
   hLth->Write();
   hLph->Write();
   hLz_c->Write();
   hLz_rc->Write();   
   hLth_c->Write();
   hLph_c->Write();
   hLss->Write();
   hLss_c->Write();
   hLp->Write();
   hLp_c->Write();   
   hLp_ec->Write();
   hRz->Write();
   hRth->Write();
   hRph->Write();
   hRz_c->Write();
   hRz_rc->Write();   
   hRth_c->Write();
   hRph_c->Write();
   hRss->Write();
   hRss_c->Write();
   hRp->Write();
   hRp_c->Write();   
   hRp_ec->Write();

   //Event Selection hist //
   hLz_es->Write();
   hLth_es->Write();
   hLph_es->Write();
   hLz_es_c->Write();
   hLz_es_rc->Write();   
   hLth_es_c->Write();
   hLph_es_c->Write();
   hLss_es->Write();
   hLss_es_c->Write();
   hLp_es->Write();
   hLp_es_c->Write();   
   hLp_es_ec->Write();
   hRz_es->Write();
   hRth_es->Write();
   hRph_es->Write();
   hRz_es_c->Write();
   hRz_es_rc->Write();   
   hRth_es_c->Write();
   hRph_es_c->Write();
   hRss_es->Write();
   hRss_es_c->Write();
   hRp_es->Write();
   hRp_es_c->Write();   
   hRp_es_ec->Write();   
   
   hLz_cut->Write();
   hRz_cut->Write();
   hct->Write();
   hct2->Write();   
   hct_cut->Write();

   hmm_Al_zp->Write();
   hmm_Al_zn->Write();
   hmm_Al_zp_c->Write();
   hmm_Al_zn_c->Write();
   hmm_select->Write();
   hmm_L->Write();
   hmm_nnL->Write();
   hmm_nnL_select->Write();
   hmm_nnL_acc->Write();
   hmm_nnL_cut_acc->Write();
   hmm_nnL_cut->Write();
   hmm_S->Write();
   hmm_cut->Write();
   hmm_cut_H->Write();
   hmm_cut_T->Write();
   hmm_Al->Write();
   hmm_Al_acc->Write();
   hmm_Al_mom->Write();
   hmm_Al_cut->Write();
   hmm_Al_cut_acc->Write();
   hmm_Al_select->Write();
   hmm_Al_select_c->Write();
   hmm_Al_nnL->Write();
   hmm_Al_nnL_zp->Write();
   hmm_Al_nnL_zn->Write();
   gchi_p->SetName("gchi");
   gchi_p->Write();
   gchi_w->SetName("gchi_w");
   gchi_w->SetMarkerStyle(20);
   gchi_w->SetMarkerColor(2);
   gchi_w->Write();
   ofp->Close();

}

/////////////////////////////////////////////////////////////

/*
void momcalib::GetACParam(){


  cout<<"==================================="<<endl;
  cout<<"========== GetACParam ============="<<endl;
  cout<<"==================================="<<endl;
  // taken by /ac/param/offset_ac.dat 
  string pname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/ac/param/offset_ac.dat";
  ifstream ifp(pname.c_str(),ios::in);
  if (ifp.fail()){ cerr << "failed open files" <<pname.c_str()<<endl; exit(1);}
  cout<<" Param file : "<<pname.c_str()<<endl;
  
  string buf;
  int AC,Seg;
  double off,pe;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> AC >> Seg >> off >> pe;

    if(AC==1){
      ac1_off[Seg]=off;
      ac1_1pe[Seg]=pe;
    }else if(AC==2){
      ac2_off[Seg]=off;
      ac2_1pe[Seg]=pe;
    }else{
      cout<<"Error :"<<endl; exit(1);
    }

  }
}

*/

//////////////////////////////////////////////////////////

/*
double momcalib::AC_npe(int nac,int seg,double adc){

  

  double npe,ac_off,ac_1pe;

  if(nac==1){
    ac_off=ac1_off[seg];
    ac_1pe=ac1_1pe[seg];
  }else if(nac==2){
    ac_off=ac2_off[seg];
    ac_1pe=ac2_1pe[seg];
  }else {
    cout<<"Error : falid Get AC parameters "<<endl; exit(1);}

  //  npe=(adc)/(ac_1pe - ac_off); // Just correct gain
  //  if(nac==1)npe=(adc)/(ac_1pe - ac_off)*2.0;     // Gogami AC DB was changed gain 400 -> 200
  if(nac==1)npe=(adc)/(ac_1pe - ac_off);     //  gain 400
  else if(nac==2)  npe=(adc)/(ac_1pe - ac_off);
    // in this case, we need scale gain parameter 2 times
  return npe;  
}
*/

//=========================================================//
//================= Function ==============================//
//=========================================================//

// ####################################################
double s2f1_off(int i,char* ARM,char* MODE, int KINE){
// ####################################################

  double RS2_offset[16],LS2_offset[16];
  if(MODE=="H" && KINE==2){

    double RS2_R_off[16]={-8361.42, -8395.25, -8414.89, -8419.06, -8362.64, -8381.55, -8370.53, -8392.66, -8389.77, -8393.96, -8388.11, -8381.73, -8333.95, -8348.93, -8363.93, -8360.30};
    double RS2_L_off[16]={-8473.92, -8470.25, -8489.89, -8494.06, -8512.64, -8494.05, -8520.53, -8505.16, -8502.27, -8468.96, -8500.61, -8494.23, -8521.45, -8498.93, -8476.43, -8472.80};
    double LS2_R_off[16]={-12441.14, -12490.70, -12579.43, -12601.39, -12471.56, -12471.38, -12658.08, -12656.28, -12690.65, -12489.77, -12701.56, -12675.30, -12696.36, -12491.35, -12709.36, -12539.99 };
    double LS2_L_off[16]={-12141.14, -12190.70, -12091.93, -12076.39, -12209.06, -12208.88, -12058.08, -12056.28, -12015.65, -12227.27, -12026.56, -12000.30, -11983.86, -12191.35, -11996.86, -12239.99 };

    double  RS2_off_H2[16];
    double  LS2_off_H2[16];
    for(int l=0;l<16;l++){
      RS2_off_H2[l]=RS2_R_off[l] + RS2_L_off[l];
      LS2_off_H2[l]=LS2_R_off[l] + LS2_L_off[l];
    }
     // double  RS2_off_H2[16]={-16911.4,-16864.3,-16900,-16897,-16873.8,-16868.4,-16901.1,-16876.8,-16895.4,-16860.9,-16893.1,-16884.4,-16847.3,-16842.7,-16836.9,-16882.6};
      // double  LS2_off_H2[16]={-25336.9,-25386.6,-25367.5,-25392.3,-25391.1,-25386.2,-25422,-25428.9,-25417.3,-25426.8,-25438.7,-25383.4,-25396,-25418.5,-25436.4,-26082.1};

      
  LS2_offset[i]=LS2_off_H2[i];
  RS2_offset[i]=RS2_off_H2[i];
  }


  if(MODE=="H" && KINE==1){

double  RS2_off_H1[16]={-16828.7,-16863,-16894,-16893.3,-16870.9,-16867.2,-16900.3,-16876.8,17554.1,-16861.6,-16895,-16890.7,-16854.6,-16852.9,-16850.5,-16861.9};
double  LS2_off_H1[16]={-25335,-25385.6,-25367,-25392.1,-25391.7,-25386.4,-25422.1,-25428.9,-25414.9,-25424.7,-25436.9, -25381.2,-25390,-25413.4,-25428.7,-26640.8};
  LS2_offset[i]=LS2_off_H1[i];
  RS2_offset[i]=RS2_off_H1[i];
  }

 double s2f1_offset; 
 if(ARM=="R")s2f1_offset=RS2_offset[i];
 else  if(ARM=="L")s2f1_offset=LS2_offset[i];
 else {cout<<"false read out !!"<<endl;}

  return s2f1_offset;

  }





// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  
  const double sigma = 0.002;//[GeV/c^2]
  double ztR      = 0.0;
  double refpos   = 0.0;
  double residual = 0.0;
  double ang      = 0.0;
  double sspos    = 0.0;
  double total_chi2 = 0.0;
  double mass,mass_MgL,mass_nnL;
  double Ee,Ee_,Ek;
  double rbeta,lbeta;
  double ref_rp=0.0;
  double ref_lp=0.0;
  double rp_ref=1.8;
  double lp_ref=2.1;
  //  double rp,lp;
  double Chi2=0.0;
  double chi2=0.0;
  double rx_c,ry_c,rth_c,rph_c,rz_c,lx_c,ly_c,lth_c,lph_c;
  double par_R[nParamTp],par_L[nParamTp];
  double parR[nParamTp],parL[nParamTp];
  double lp_x;
  double lp_y;
  double lp_z;
  double rp_x;
  double rp_y;
  double rp_z;
  double Rp,Lp,Bp;
  double MgL_weight;
  //  momcalib* mom=new momcalib();
  double dpk, dpe_,dpe;
  TVector3 L_v, R_v, B_v;
  //======= Coincidence analysis ========//
  
  for(int i=0;i<nParamTp;i++){
    if(MODE==0){
      par_R[i]=param[i];
      par_L[i]=param[i+nParamTp];
      parR[i]=par_R[i];
      parL[i]=par_L[i];
    }else if(MODE==-1){
      par_R[i]=param[i];
      par_L[i]=Opt_par_L[i];
    }else if(MODE==1){
      par_R[i]=Opt_par[i];
      par_L[i]=param[i];
      //      cout<<"param "<<i<<" : "<<par_L[i]<<endl;
    }else {
      par_R[i]=Opt_par[i];
      par_L[i]=Opt_par_L[i];
    }
  }


  //======================================//

  
  for(int i=0 ; i<ntune_event ; i++){

    //====== Initialization ========//

    residual =0.0;
    ref_rp=0.0;
    ref_lp=0.0;
    ref_rp=rp[i];
    ref_lp=lp[i];
    mass  = 0.0;
    mass_MgL = 0.0;
    mass_nnL = 0.0;
    dpe_=0.0;
    dpk=0.0;
    lp_x=0.0;
    lp_y=0.0;
    lp_z=0.0;
    rp_x=0.0;
    rp_y=0.0;
    rp_z=0.0;
    Rp=0.0;
    Lp=0.0;
    Bp=0.0;


    Rp = rp[i];
    Lp = lp[i];
    Bp = bp[i];

    dpk  = drp[i];
    dpe_ = dlp[i];
    dpe  = dbp[i];



    if(nrun[i]<=111220 || (111480<=nrun[i] && 111542>=nrun[i])) scale[i]=false; // pe=2.1 GeV
    else scale[i]=true;


    if(scale[i])lp_ref=2.21807;
    else lp_ref=2.1;


    // ===== Al range Flag ======//
    //    if(mass_flag[i]==2 && Al_range<fabs(MM[i]-mass_ref[i]))continue;

    
    // RHRS //

    //======== Scaled paramters ===============//

    rx_fp[i]  = (rx_fp[i]-XFPm)/XFPr;
    rth_fp[i] = (rth_fp[i]-XpFPm)/XpFPr;
    ry_fp[i]  = (ry_fp[i]-YFPm)/YFPr;
    rph_fp[i] = (rph_fp[i]-YpFPm)/YpFPr;
    rz[i]     = (rz[i]-Ztm)/Ztr;
    Rp        = (Rp - PRm)/PRr;

									     
    //==== momentum tuning ==========//

    if(MT_f[8])Rp = calcf2t_mom(par_R, rx_fp[i], rth_fp[i], ry_fp[i], rph_fp[i],rz[i]);


    //==== Scaled to Nomal =======//    
    rx_fp[i]  = rx_fp[i]  * XFPr + XFPm;
    rth_fp[i] = rth_fp[i] * XpFPr + XpFPm;
    ry_fp[i]  = ry_fp[i]  * YFPr + YFPm;
    rph_fp[i] = rph_fp[i] * YpFPr + YpFPm;
    rz[i]     = rz[i] * Ztr +Ztm;
    Rp        = Rp * PRr +PRm;    
    //=============================//
   

    // LHRS //
    
    //======== Scaled paramters ===============//

    lx_fp[i]  = (lx_fp[i]-XFPm)/XFPr;
    lth_fp[i] = (lth_fp[i]-XpFPm)/XpFPr;
    ly_fp[i]  = (ly_fp[i]-YFPm)/YFPr;
    lph_fp[i] = (lph_fp[i]-YpFPm)/YpFPr;
    lz[i]     = (lz[i]-Ztm)/Ztr;
    Lp        = (Lp - PLm)/PLr;
    //==== momentum tuning ==========//

    if(MT_f[9])Lp = calcf2t_mom(par_L, lx_fp[i], lth_fp[i], ly_fp[i], lph_fp[i], lz[i]);

    //==== Scaled to Nomal =======//
    
    lx_fp[i]  = lx_fp[i]  * XFPr + XFPm;
    lth_fp[i] = lth_fp[i] * XpFPr + XpFPm;
    ly_fp[i]  = ly_fp[i]  * YFPr + YFPm;
    lph_fp[i] = lph_fp[i] * YpFPr + YpFPm;
    lz[i]     = lz[i]*Ztr + Ztm;
    Lp        = Lp   *PLr + PLm;
    
    //=============================//


    if(scale[i] && MT_f[9])Lp = Lp *2.21807/2.1; 


    //==== E loss Calib =====//

    Bp = Bp - dpe;
    Rp = Rp + dpk;
    Lp = Lp + dpe_;
        
    //------ Set Physics value --------//
    Ee=sqrt(pow(Bp,2)+pow(Me,2));
    Ee_=sqrt(pow(Lp,2)+pow(Me,2));
    Ek=sqrt(pow(Rp,2)+pow(MK,2));
    
    rbeta=Rp/Ek; 
    lbeta=Lp/Ee_; 

    //==== Right Hand coordinate =====//
    //    rp_z = Rp/sqrt(1.0*1.0 + pow( rth[i]) ,2.0) + pow(  rph[i] ),2.0) );
    //    lp_z = Lp/sqrt(1.0*1.0 + pow( lth[i]) ,2.0) + pow(  lph[i] ),2.0) );
    rp_z = Rp/sqrt(1.0*1.0 + rth[i]*rth[i] + rph[i]*rph[i]);
    lp_z = Lp/sqrt(1.0*1.0 + lth[i]*lth[i] + lph[i]*lph[i] );

    lp_x = lp_z *    lth[i] ;
    lp_y = lp_z *    lph[i] ; 
    rp_x = rp_z *    rth[i] ;
    rp_y = rp_z *    rph[i] ;


    B_v.SetXYZ(0.0,0.0,sqrt(Ee*Ee-Me*Me));
    R_v.SetXYZ(rp_x, rp_y, rp_z);
    L_v.SetXYZ(lp_x, lp_y, lp_z);


    R_v.RotateX(  13.2/180.*3.14);
    L_v.RotateX( -13.2/180.*3.14);


    mass = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
		 - (B_v - L_v - R_v)*(B_v - L_v - R_v) );

    mass_MgL = sqrt( (Ee + MAl - Ee_ - Ek)*(Ee + MAl - Ee_ - Ek)
		     - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
    
    mass_nnL = sqrt( (Ee + MTr - Ee_ - Ek)*(Ee + MTr - Ee_ - Ek)
		     - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
    

    //======================//
    //====== Residual ======//
    //======================//




    if(mass_flag[i]==0 && scale[i]==0)      residual = (mass - mass_ref[i]); // Lambda H kinematics
    else if(mass_flag[i]==1)                residual = (mass - mass_ref[i])*weight;  // weight Sigma0:
    else if(scale[i] && mass_flag[i]==0)    residual = (mass - mass_ref[i])*weightT;  // weight Lambda T kinematics  
    else if(mass_flag[i]==3 )               residual = (mass_nnL - mass_ref[i])*weightnnL;  // weight nnL
    else if(mass_flag[i]==2 )               residual = (mass_MgL - mass_ref[i])*weightAl;  // MgL events


    
    // conditon //
    //    if((fabs( Rp -rp_ref )>0.1 || fabs( Lp - lp_ref) >0.1 ))residual=residual*1.5; //rp_ref=1.8 GeV, lp_ref=2.1 GeV


    //    if(mass_flag[i]==2)cout<<"i "<<i<<" runnum "<<nrun[i]<<" mass "<<mass_MgL<< " mass dif "<<(mass_MgL-mass_ref[i])*1000.<<" mass_ref "<<mass_ref[i]<<" weightAl "<<weightAl<<" residual "<<residual<<" Rp "<<rp[i]<<" Rp_c "<<Rp<<" Lp "<<lp[i]<<" Lp_c "<<Lp<<" scale "<<scale[i]<<endl;



    //    if(pow(residual/sigma,2.0)/ntune_event>100)
    //cout<<"i "<<i<<" mass_flag "<<mass_flag[i]<<" mass "<<mass_MgL<<" mass_ref "<<mass_ref[i]<<" residual "<<residual<<" chi "<<pow(residual/sigma,2.0)/ntune_event<<" chi_sum "<<chi2<<" weight "<<weightAl<<endl;

 
     chi2=chi2 + pow(residual/sigma,2.0)/ntune_event;
     
     //     if(mass_ref[i]>MMgL-1. && mass_ref[i]!=MnnL)cout<<"i "<<i<<" mm_MgL "<<mass_MgL<<" mass_ref "<<mass_ref[i]<<" residual "<<residual<<" weight "<<weightAl<<" chi "<<chi2<<endl;
    
  }//end for 

  
  //  cout<<" chi2 "<<chi2 <<endl;
  
  
  fval=chi2;
  //  delete mom;  
  
  
}//end fcn



// #############################################################
void fcn_new(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
    cout<<"test new"<<endl;
  double sigma = 0.002;//[GeV/c^2]
  double mean  =0.0;
  double ztR      = 0.0;
  double refpos   = 0.0;
  double residual = 0.0;
  double ang      = 0.0;
  double sspos    = 0.0;
  double total_chi2 = 0.0;
  double mass,mass_MgL,mass_nnL;
  double Ee,Ee_,Ek;
  double rbeta,lbeta;
  double ref_rp=0.0;
  double ref_lp=0.0;
  double rp_ref=1.8;
  double lp_ref=2.1;
  //  double rp,lp;
  double Chi2=0.0;
  double chi2=0.0;
  double rx_c,ry_c,rth_c,rph_c,rz_c,lx_c,ly_c,lth_c,lph_c;
  double par_R[nParamTp],par_L[nParamTp];
  double parR[nParamTp],parL[nParamTp];
  double lp_x;
  double lp_y;
  double lp_z;
  double rp_x;
  double rp_y;
  double rp_z;
  double Rp,Lp,Bp;
  double MgL_weight;
  double x[1];
  double xmax[1];
  double PD;
  //  momcalib* mom=new momcalib();
  double dpk, dpe_,dpe;
  TVector3 L_v, R_v, B_v;
  //======= Coincidence analysis ========//

  for(int i=0;i<nParamTp;i++){
    if(MODE==0){
      par_R[i]=param[i];
      par_L[i]=param[i+nParamTp];
      parR[i]=par_R[i];
      parL[i]=par_L[i];
    }else if(MODE==-1){
      par_R[i]=param[i];
      par_L[i]=Opt_par_L[i];
    }else if(MODE==1){
      par_R[i]=Opt_par[i];
      par_L[i]=param[i];
      //      cout<<"param "<<i<<" : "<<par_L[i]<<endl;
    }else {
      par_R[i]=Opt_par[i];
      par_L[i]=Opt_par_L[i];
    }
  }



  //======================================//

  
  for(int i=0 ; i<ntune_event ; i++){

    //====== Initialization ========//

    residual =0.0;
    ref_rp=0.0;
    ref_lp=0.0;
    ref_rp=rp[i];
    ref_lp=lp[i];
    mass  = 0.0;
    mass_MgL = 0.0;
    mass_nnL = 0.0;
    dpe_=0.0;
    dpk=0.0;
    lp_x=0.0;
    lp_y=0.0;
    lp_z=0.0;
    rp_x=0.0;
    rp_y=0.0;
    rp_z=0.0;
    Rp=0.0;
    Lp=0.0;
    Bp=0.0;
    PD=0.0;

    Rp = rp[i];
    Lp = lp[i];
    Bp = bp[i];

    dpk  = drp[i];
    dpe_ = dlp[i];
    dpe  = dbp[i];
    x[0]=100.0;
    xmax[0]=0.0;
    sigma= 100.0;
    mean=100.0;
    if(nrun[i]<=111220 || (111480<=nrun[i] && 111542>=nrun[i])) scale[i]=false; // pe=2.1 GeV
    else scale[i]=true;


    if(scale[i])lp_ref=2.21807;
    else lp_ref=2.1;


    // RHRS //

    //======== Scaled paramters ===============//

    rx_fp[i]  = (rx_fp[i]-XFPm)/XFPr;
    rth_fp[i] = (rth_fp[i]-XpFPm)/XpFPr;
    ry_fp[i]  = (ry_fp[i]-YFPm)/YFPr;
    rph_fp[i] = (rph_fp[i]-YpFPm)/YpFPr;
    rz[i]     = (rz[i]-Ztm)/Ztr;
    Rp        = (Rp - PRm)/PRr;

									     
    //==== momentum tuning ==========//

    if(MT_f[8])Rp = calcf2t_mom(par_R, rx_fp[i], rth_fp[i], ry_fp[i], rph_fp[i],rz[i]);


    //==== Scaled to Nomal =======//    
    rx_fp[i]  = rx_fp[i]  * XFPr + XFPm;
    rth_fp[i] = rth_fp[i] * XpFPr + XpFPm;
    ry_fp[i]  = ry_fp[i]  * YFPr + YFPm;
    rph_fp[i] = rph_fp[i] * YpFPr + YpFPm;
    rz[i]     = rz[i] * Ztr +Ztm;
    Rp        = Rp * PRr +PRm;    
    //=============================//
   

    // LHRS //
    
    //======== Scaled paramters ===============//

    lx_fp[i]  = (lx_fp[i]-XFPm)/XFPr;
    lth_fp[i] = (lth_fp[i]-XpFPm)/XpFPr;
    ly_fp[i]  = (ly_fp[i]-YFPm)/YFPr;
    lph_fp[i] = (lph_fp[i]-YpFPm)/YpFPr;
    lz[i]     = (lz[i]-Ztm)/Ztr;
    Lp        = (Lp - PLm)/PLr;
    //==== momentum tuning ==========//

    if(MT_f[9])Lp = calcf2t_mom(par_L, lx_fp[i], lth_fp[i], ly_fp[i], lph_fp[i], lz[i]);

    //==== Scaled to Nomal =======//
    
    lx_fp[i]  = lx_fp[i]  * XFPr + XFPm;
    lth_fp[i] = lth_fp[i] * XpFPr + XpFPm;
    ly_fp[i]  = ly_fp[i]  * YFPr + YFPm;
    lph_fp[i] = lph_fp[i] * YpFPr + YpFPm;
    lz[i]     = lz[i]*Ztr + Ztm;
    Lp        = Lp   *PLr + PLm;
    
    //=============================//


    if(scale[i] && MT_f[9])Lp = Lp *2.21807/2.1; 


    //==== E loss Calib =====//


    Bp = Bp - dpe;
    Rp = Rp + dpk;
    Lp = Lp + dpe_;
        
    //------ Set Physics value --------//
    Ee=sqrt(pow(Bp,2)+pow(Me,2));
    Ee_=sqrt(pow(Lp,2)+pow(Me,2));
    Ek=sqrt(pow(Rp,2)+pow(MK,2));
    
    rbeta=Rp/Ek; 
    lbeta=Lp/Ee_; 

    //==== Right Hand coordinate =====//
    /*
    rp_z = Rp/sqrt(1.0*1.0 + rth[i]*rth[i] + rph[i]*rph[i] );
    lp_z = Lp/sqrt(1.0*1.0 + lth[i]*lth[i] + lph[i]*lph[i] );
    lp_x = lp_z *   lth[i] ;
    lp_y = lp_z *   lph[i] ; 
    rp_x = rp_z *   rth[i] ;
    rp_y = rp_z *   rph[i] ;
    */
    rp_z = Rp/sqrt(1.0*1.0 + rth[i]*rth[i] + rph[i]*rph[i] );
    lp_z = Lp/sqrt(1.0*1.0 + lth[i]*lth[i] + lph[i]*lph[i] );
    lp_x = lp_z *  lth[i] ;
    lp_y = lp_z *  lph[i] ; 
    rp_x = rp_z *  rth[i] ;
    rp_y = rp_z *  rph[i] ;    

    B_v.SetXYZ(0.0,0.0,sqrt(Ee*Ee-Me*Me));
    R_v.SetXYZ(rp_x, rp_y, rp_z);
    L_v.SetXYZ(lp_x, lp_y, lp_z);


    R_v.RotateX(  13.2/180.*3.14);
    L_v.RotateX( -13.2/180.*3.14);
    

    mass = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
		 - (B_v - L_v - R_v)*(B_v - L_v - R_v) );

    mass_MgL = sqrt( (Ee + MAl - Ee_ - Ek)*(Ee + MAl - Ee_ - Ek)
		     - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
    
    mass_nnL = sqrt( (Ee + MTr - Ee_ - Ek)*(Ee + MTr - Ee_ - Ek)
		     - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
    


    //======================//
    //====== Residual ======//
    //======================//


    if(mass_flag[i]==0 || mass_flag[i]==1)      residual = mass - mass_ref[i]; // Lambda H kinematics
    else if(mass_flag[i]==2 )                    residual = mass_MgL - mass_ref[i];  // MgL events
    else if(mass_flag[i]==3 )                    residual = mass_nnL - mass_ref[i];  // weigth nnL


    if(mass_flag[i]==0 || mass_flag[i]==1 ){

      if(mass_flag[i]==0){       xmax[0] = paramL[3]; mean = Lam_mean; sigma = Lam_sigma; 
      }else if(mass_flag[i]==1){ xmax[0] = paramS[3]; mean = Sig_mean; sigma = Sig_sigma;}
       
      x[0] = mass - mass_ref[i] - mean;
      
      if(mass_flag[i]==0)PD = (expgaus(x,paramL)/expgaus(xmax,paramL))*(expgaus(x,paramL)/expgaus(xmax,paramL));
      else if(mass_flag[i]==1) PD = (expgaus(x,paramS)/expgaus(xmax,paramS))*(expgaus(x,paramS)/expgaus(xmax,paramS));
      residual = residual/ PD;

    }else {sigma=0.002; mean=0.0;}
    //====== weight =======//

    if(mass_ref[i]==MS0)                    residual *= weight;  // weigth Sigma0:
    else if(scale[i] && mass_ref[i]==ML)    residual *= weightT;  // weigth Lambda T kinematics  
    else if(mass_flag[i]==2 )               residual *= weightAl;  // MgL events
    else if(mass_flag[i]==3 )               residual *= weightnnL;  // weigth nnL
    if((fabs( Rp -rp_ref )>0.1 || fabs( Lp - lp_ref) >0.1 ))residual=residual*3.0; //rp_ref=1.8 GeV, lp_ref=2.1 GeV



    chi2=chi2 + pow(residual/sigma,2.0)/ntune_event;


    //    cout<<"i "<<i<<" sigma "<<sigma<<" residual "<<residual<<" pow "<<pow(residual/sigma,2.0)/ntune_event<<endl;

  }//end for 
  
  
  //  cout<<" chi2 "<<chi2 <<endl;
  
  
  fval=chi2;
  //  delete mom;  
  
  
}//end fcn




// #############################################################
double momcalib::tune(double* pa, int j, int MODE) 
// #############################################################
{

  
  double chi = 0.0;
  double arglist[10]; 
  int ierflg = 0;
  int allparam = nParamTp;
  int arm=1;

  if(MODE==0){
    allparam=nParamTp*2;
    arm=2;  }



  cout<<"mode "<<MODE<<" allParam "<<allparam<<endl;


  TMinuit* minuit= new TMinuit(allparam);

  //  if(form_mode)  minuit->SetFCN(fcn_new); // fcn Chi-square function as a test
  //  else
    minuit->SetFCN(fcn); // fcn Chi-square function
  
  double start[allparam];
  double step[allparam];


  const int nMatT =nnp;  
  const int nXf   =nnp;
  const int nXpf  =nnp;
  const int nYf   =nnp;
  const int nYpf  =nnp;
  const int nZt   =nnp;
  const int nnn   =nnp;
  int nfix=0;
  //  nfix=3;

  // The number of order is reduced for test (4-->2)


  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;


  for(int f=0;f<arm;f++){
    for (int n=0;n<nMatT+1;n++){
      for(e=0;e<n+1;e++){
	for (d=0;d<n+1;d++){
	  for (c=0;c<n+1;c++){ 
	    for (b=0;b<n+1;b++){
	      for (a=0;a<n+1;a++){ 
		if (a+b+c+d+e==n){
		  if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt && a+b+c+d+e<=nnn){

		    if(n<nfix){
		    start[npar] = pa[npar];
		    step[npar] = 0.0;

		    }else{

		    start[npar] = pa[npar];
		    step[npar] = 1.0e-3;
		    }

		  }
		  else{
		    start[npar] = 0.0;
		    step[npar] = 0.0;

		  }
		  npar++;

		}
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
   
    LLim[i] = pa[i] - 5.0; // temp
    ULim[i] = pa[i] + 5.0; // temp
    //    step[i] = 10.0;  // test
    //    LLim[i] = pa[i]; // temp
    //    ULim[i] = pa[i]; // temp

    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
	  
  }




  
  // ~~~~ Strategy ~~~~
  //  arglist[0] = 2.0; // original
  arglist[0] = 1.0; // test
  //arglist[0] = 0.0;   // test
  minuit->mnexcm("SET STR",arglist,1,ierflg);
 
  // ~~~~ Migrad + Simplex  ~~~~ 
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg); // Chi-square minimization
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double er;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);

  if(amin>0) chi=amin;

  
  for(int i=0 ; i<allparam ; i++){
    if(MODE==0){
      minuit -> GetParameter(i,Opt_par[i],er);
      if(i<nParamTp){minuit -> GetParameter(i,Opt_par_R[i],er);}// RHRS momentum parameter
      else if(nParamTp<=i){minuit -> GetParameter(i,Opt_par_L[i-nParamTp],er);}// RHRS momentum parameters
      
      
    }else if(MODE==-1){minuit -> GetParameter(i,Opt_par_R[i],er);// RHRS momentum parameters
    }else if(MODE==1){minuit -> GetParameter(i,Opt_par_L[i],er);// LHRS momentum parameters
    }
  }

  
  return chi;
}

 
// ###################################################
double calcf2t_mom(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt)
// ###################################################
{
  // ----------------------------------------------------------------- //
  // ------ 4rd order using xf, xpf, yf, ypf, zt, xt, xpt, yt, ytp --- //
  // ----------------------------------------------------------------- //



  
  const int nMatT =nnp;  
  const int nXf   =nnp;
  const int nXpf  =nnp;
  const int nYf   =nnp;
  const int nYpf  =nnp;
  const int nZt   =nnp;

  
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
		
		Y += x*P[npar]; 
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



// ###################################################
double calcf2t_ang(double* P, double xf, double xpf, 
		     double yf, double ypf, double zt)
// ####################################################
{
  // ------------------------------------------------ //
  // ----- 4rd order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //
  
  const int nMatT=nn;  
  const int nXf=nn;
  const int nXpf=nn;
  const int nYf=nn;
  const int nYpf=nn;
  const int nZt=nn;
  
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
		Y += x*P[npar]; 
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
      	

// #################################################
double calcf2t_zt(double* P, double xf, double xpf, 
                 double yf, double ypf){
// ###############################################


  const int nMatT=nnz; 
  const int nXf=nnz;
  const int nXpf=nnz;
  const int nYf=nnz;
  const int nYpf=nnz;

  
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


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus(double *x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
  //  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus_mean(double *x, double *par){
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
  //  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));
  val = val * x[0];
  
  return val;

}

double expgaus_sigma(double *x, double *par){
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
  //  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));
  val = val *(x[0] - par[4])* (x[0] - par[4]);
  
  return val;

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double Get_expgaus_mean(double *par){

  TF1* fmean = new TF1("fmean","expgaus_mean",-0.300,0.300,4);
  fmean->SetParameters(par[0],par[1],par[2],par[3]);
  double  mean= fmean-> Integral(-0.05,0.1);
  
  return mean;

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double Get_expgaus_sigma(double *par, double mean){
 

  TF1* fsigma = new TF1("fsigma","expgaus_sigma",-0.300,0.300,5);
  fsigma->SetParameters(par[0],par[1],par[2],par[3], mean);
  double sigma = fsigma-> Integral(-0.05,0.1);
   
  return sigma;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#endif
