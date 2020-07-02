#ifndef xsec_h
#define xsec_h 1
#include <TLorentzVector.h>
#include "define.h"
#include "Param.h"

const int nmax =1000;

class xsec{


  
public:
  xsec();
  ~xsec();
  void SetParam();
  void Calc_XS();
  void SetLAccept(string ifpname);
  void SetRAccept(string ifpname);
  TVector3 HallACoordinate(TVector3 V, bool rhrs);
  double VP_Flux(TLorentzVector B_v, TLorentzVector L_v);
  double NGamma(double Qc);
  double dOmega_L(double Ee_);
  double dOmega_R(double Ek);
  double NTarget(string mode);

  double Lpx,Lpy,Lpz,Bp,Ee,Ee_,Lp,Lp_mean;
  double hrs_angle = 13.2*PI/180.;
  double gen_theta_accep =  6.7272080607340551e-2;
  double Qc;
  double kaon_survival_ratio;
  // ----- SetXS ----//
  int nEe_=0; int nEk=0;
  int nth =100;
  int nph =100;
  double LE[nmax],LOmega[nmax],LEMax,LEMin;
  double RE[nmax],ROmega[nmax];
  // ------ NGamma ------//
  TLorentzVector B_v,L_v;
  
};

xsec::xsec(){};
xsec::~xsec(){};

#endif
