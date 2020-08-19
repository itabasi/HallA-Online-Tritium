#ifndef xsec_h
#define xsec_h 1
#include <TLorentzVector.h>
#include "define.h"
#include "Param.h"
#include <complex.h>
const int nmax =200;
const int nmax2 =200;
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
  complex<double> vlkkp(double skp,double skpp,int lcall);
  complex<double> vofq(double q2);
  //  void DGauss(double* Y,double* WY, int N);
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


  //------- PhaseShift ------//
  double ampj,ampj2,amtag,amtag2,fmu,fmunn,amu,fnucl,fnucl2,hbarc2;
  double skz,skz2,eon,eoff,rho,ekpj,ektag,lpjmax,skp,skpp,skp2,skpp2,vv1,uofq,fthecm,q2;
  double POL[nmax2],deltal[60];
  double skk[nmax2],x[nmax2],w[nmax2],wt[nmax2],pol[nmax2+1],pols[nmax2+1],xvq[nmax2],wvq[nmax2];
  complex<double> gren[nmax2+1][nmax2+1],ul[nmax2+1][nmax2+1],gren2[nmax2+1][nmax2+1],v[nmax2],ton[60],tborn[60];
  complex<double> fthe[nmax2],fthe1[nmax2],Rl[nmax2+1],V[nmax2+1];
  complex<double> xi, vq,v2,vv,delta,det,ron,fborn,fthed;
  double va,b_a,vr,b_r,b_r2,b_a2,vqa,vqr,ronr,tonre,tonim,tbornre,tbornim,perre,perim,delrad,deldeg,deldega,tcross,tcross1;
  
  const int n1 =40;
  const int n2 =40;
  //  const int n1 =4;
  //  const int n2 =4;  
  //  const int n1 =6;
  //  const int n2 =6;
  int model =1;
  const int L_value =3;
  double xxx, totalc;
  int npot =n1+n2;
  int nrp1,lmax,lmax1;
  double tol =0.01;
  bool nr_mode = true;
  //  bool nr_mode=false;
};

xsec::xsec(){};
xsec::~xsec(){};

#endif
