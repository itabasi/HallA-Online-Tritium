#ifndef xsec_h
#define xsec_h 1
#include <TLorentzVector.h>
#include "define.h"
#include "Param.h"

const int nmax =200;

class xsec{


  
public:
  xsec();
  ~xsec();
  void SetParam(string ifpname);
  //  bool SetVal(string line,string name, string val);
  string SetVal(string line,string name);
  void Calc_XS();
  void SetLAccept(string ifpname);
  void SetLAccept_new(string ifpname);
  void SetRAccept(string ifpname);
  double PhaseShift();
  TVector3 HallACoordinate(TVector3 V, bool rhrs);
  double VP_Flux(TLorentzVector B_v, TLorentzVector L_v);
  double NGamma(double Qc);
  double NHyp(double Nhyp);
  double dOmega_L(double Ee_, double theta);
  double dOmega_R(double Ek  , double theta);
  double NTarget(string mode);
  double NKaon(double Ek);
  //---------- SetParam ----------//
  string RHRS_accept_file;
  string LHRS_accept_file;
  string mode;
  double Nhyp,Nhyp2;  
  double Bp_mean,Lp_mean,Rp_mean;
  double Lpx,Lpy,Lpz,Bp,Ee,Ee_,Lp,Rp,Rpx,Rpy,Rpz,Ek;
  double theta_L_max,theta_R_max;
  double hrs_angle = 13.2*PI/180.;
  double gen_theta_accept =  6.7272080607340551e-2;
  double Qc;
  double kaon_survival_ratio;
  // ----- SetXS ----//
  int nEe_=0; int nEk=0;
  int nth   = 100;
  int nph   = 100;
  int nth_R = 100;
  int nph_R = 100;
  double LE[nmax],LOmega[nmax][nmax],LTheta[nmax],LEMax,LEMin;
  double RE[nmax],ROmega[nmax][nmax],RTheta[nmax],REMax,REMin;
  // ------ NGamma ------//
  TLorentzVector B_v,L_v,R_v;
  // ------ NHyp --------//
  
};

xsec::xsec(){};
xsec::~xsec(){};

#endif
