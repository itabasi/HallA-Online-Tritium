#ifndef t0calib_h
#define t0calib_h 1

using namespace std;

#include "Setting.h"
#include "Tree.h"
#include "ParamMan.h"
#include "define.h"
#include "Param.h"
//#include "ana_Lambda.h"


struct EventSelection{

  int RS2_seg;
  int LS2_seg;
  double Rpathl,Lpathl;
  double RS2_t,RS2_b,LS2_t,LS2_b;
  double RS2_t_off[16],RS2_b_off[16],LS2_t_off[16],LS2_b_off[16];
  double RS2_ref,LS2_ref;
  double ct_g, ct,ct_b;
  double RF1_t, RF1_b,LF1_t,LF1_b,RF1_ref,LF1_ref;
  double RS2_off[16],LS2_off[16];
  double Rp,Lp;
};

struct EventSelection ev;


struct Fill{

  int RS2_seg;
  int LS2_seg;
  double Rpathl,Lpathl;
  double RS2_t,RS2_b,LS2_t,LS2_b;
  //  double RS2_t_off[16],RS2_b_off[16],LS2_t_off[16],LS2_b_off[16];
  double RS2_ref,LS2_ref;
  double ct_g, ct;
  double RF1_t, RF1_b,LF1_t,LF1_b,RF1_ref,LF1_ref;
  double RS2_off[16],LS2_off[16];
  double Rp,Lp;
  double ct_b;
};

struct Fill tr;



class t0calib : public Tree
{

 public:
  t0calib();
  ~t0calib();
 private:
  Setting* set;
  ParamMan *param;
  
 public:
  void SetRoot(string ifname);
  void EventSelection();
  double tune(double *pa, int j);
  void Tuning(string ofname);
  void SetHist();
  void Fill();
  void NewRoot(string ofrname);
  //--------------------------------//
  //--------< Parameters > ---------//
  //--------------------------------//

  string ofroot;
  TFile* ofp;
  TFile* ofr;
  TTree* T;
  
  double RS2_off[16],LS2_off[16],coin_offset;
  double momR,momL;
  double Rpathl, Lpathl;
  int Rs2_pad, Ls2_pad;
  //  int niter;
  double ct;
  double RF1[100],LF1[100];

  //  int nS2 =16;
  TH1D* hct_select;
  TH1D* hct_ev;
  TH1D* hct_select_c;
  TH1D* hct_c;

  TH1D* hct_RS2[nS2];
  TH1D* hct_LS2[nS2];
  TF1*  fct_RS2[nS2];
  TF1*  fct_LS2[nS2];
  
};

t0calib::t0calib(){
  set = new Setting();
  ofp = new TFile(ofroot.c_str());
  

};

t0calib::~t0calib(){};


#endif
