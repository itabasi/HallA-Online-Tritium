/*
  CalcKinematics_eek.hh
  
  Toshiyuki Gogami, October 16, 2017
*/

#ifndef CalcKinematics_H
#define CalcKinematics_H 1

#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include "const.hh"
using namespace std;


class CalcKinematics_eek{
public:
  CalcKinematics_eek();
  ~CalcKinematics_eek();
  
private:
  double theta_spec_cent;
  int flag;    // Switch of hypernuclei (hyperons)
  double mt, mhyp, AMN_t;
  double q;    // Momentum of virtual photon (MeV/c)
  double Qsq;  // Qsquare (MeV/c)2^{2}
  double theta_gam; // Theta between e and photon
  double phi_gam;   // Phi
  double theta_gk;  // Theta between photon and K+
  double phi_gk;    // Phi
  double Eep;   // Energy of e beam (MeV)
  double pbeam; // Momentum of e beam (MeV/c)
  double pepri; // Momentum of e' (MeV/c)
  double pkaon; // Momentum of Kaon+ (MeV/c)
  double pLamb; // Momentum transfer to Lambda (MeV/c)
  double omega; // Ee - Eep (MeV)
  double xp_eep, yp_eep;
  double xp_eg, yp_eg;
  double xp_gk, yp_gk;
  double xp_k, yp_k;
  double xt,yt,zt; // Production point used for further simualtions
  double evid;     // Event ID
  
  double mmcheck;
  double mmcheck2;
  
public:
  //void Calc(double,double,double,double);
  int Calc(double,double,double,double);
  double CheckMM();
  double CheckMM2();

public:
  void SetFlag(int dragon_flag){ flag=dragon_flag;};
  TVector3 GetEpMom(){return TVector3(pepri,xp_eep,yp_eep);};
  TVector3 GetKMom() {return TVector3(pkaon,xp_k,yp_k);};
  void SetGenPos(double x,double y,double z){
    xt = x;
    yt = y;
    zt = -z;
  };
  TVector3 GetGenPos()     {return TVector3(xt,yt,zt);};
  void     SetEventID(int id){evid = id;};
  int      GetEventID()    {return evid;};
  double   GetOmega()      {return omega;};
  double   GetLambdaMom()  {return pLamb;};
  double   GetQsquare()    {return Qsq;};
  double   GetThetaEGamma(){return theta_gam;};
  double   GetThetaGammaK(){return theta_gk;};
  double   GetCheckedMM()  {return mmcheck;};
  double   GetAssumedMhyp(){return mhyp;};
  int GetFlag(){return flag;};
  
};

#endif
