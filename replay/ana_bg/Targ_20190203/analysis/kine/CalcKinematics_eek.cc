/*
  CalcKinematics_eek.cc
  
  Toshiyuki Gogami, October 16, 2017
*/

#include "CalcKinematics_eek.hh"
#include "const.hh"
#include <TVector3.h>
#include <TRandom3.h>
#include <iostream>
#include <fstream>
using namespace std;

CalcKinematics_eek::CalcKinematics_eek()
  :mt(0.0), mhyp(0.0), AMN_t(1.0),
   q(0.0), Qsq(0.0), 
   theta_gam(0.0), phi_gam(0.0),
   theta_gk(0.0), phi_gk(0.0),
   Eep(0.0), pbeam(0.0),
   pepri(0.0), pkaon(0.0), pLamb(0.0), 
   omega(0.0), 
   xp_eep(0.0), yp_eep(0.0),
   xp_eg(0.0), yp_eg(0.0),
   xp_gk(0.0), yp_gk(0.0),
   xp_k(0.0), yp_k(0.0),
   xt(0.0), yt(0.0), zt(0.0),
   evid(0)   
{}

CalcKinematics_eek::~CalcKinematics_eek()
{}

int CalcKinematics_eek::Calc(double pe, double pep, 
			     double theta_ep, double phi_ep)
{
  double Ebeam = sqrt(pe*pe + me*me); // Energy of beam electron
  pepri = pep; // Momentum of scattered electron
  mt = mp;     // Mass of a target nucleus
  mhyp = mL;   // Mass of a hypernucleus
  AMN_t = 1.0; // Atomic number
  
  if(flag==1){ // p(e,e'K+)L
    mt = mp, mhyp = mL, AMN_t=1;
  }
  else if(flag==300){ // 3H(e,e'K+)nnL
    mt = mH3, mhyp = mnnL, AMN_t=3;
  }
  else if(flag==400){ // 4He(e,e'K+)4LH
    mt = mHe4, mhyp = mH3 + mL - 2.0, AMN_t=4;
  }
  else if(flag==1200){ // 12C(e,e'K+)12LB
    mt = mC12, mhyp = mB12L, AMN_t=12;
  }
  else if(flag==4000){ // 40Ca(e,e'K+)40LK
    mt = mCa40, mhyp = mK39+mL-20.0, AMN_t=40;
  }
  else if(flag==4800){ // 48Ca(e,e'K+)48LK
    mt = mCa48, mhyp = mK47+mL-20.0, AMN_t=48;
  }
  else {
    flag = 11111;
    mt = mp, mhyp = mL, AMN_t=1;
  }
  
  double msca = mk; // kaon+
  Eep   = sqrt(pep*pep+me*me);
  omega = Ebeam - Eep;

  pbeam = sqrt(Ebeam*Ebeam - me*me);
  TVector3 pvec_e (0.0, 0.0, pbeam);
  TVector3 pvec_ep (pep*sin(theta_ep)*cos(phi_ep), 
		    pep*sin(theta_ep)*sin(phi_ep),
		    pep*cos(theta_ep));

  q = sqrt((pvec_e-pvec_ep)*(pvec_e-pvec_ep));
  Qsq = -(omega*omega - q*q); // GeV^{2}
  theta_gam = asin(pep/q*sin(theta_ep));
  phi_gam = phi_ep + 3.14159;
  TVector3 pvec_gam (q*sin(theta_gam)*cos(phi_gam), 
		     q*sin(theta_gam)*sin(phi_gam),
		     q*cos(theta_gam));
  
  // ========================================== //
  // ==== Random vlues for angles between ===== //
  // ==== gamma and K+.                   ===== //
  // ========================================== //
  TRandom3* rnd1 = new TRandom3();
  rnd1->SetSeed();
  theta_gk = acos(rnd1->Uniform(0.99,1.00));
  rnd1->SetSeed();
  phi_gk = rnd1->Uniform(0.0,3.14159*2.0);
  
  xp_eep = tan(theta_ep) * cos(phi_ep);
  yp_eep = tan(theta_ep) * sin(phi_ep);
  
  xp_eg = tan(theta_gam) * cos(phi_gam);
  yp_eg = tan(theta_gam) * sin(phi_gam);
  
  xp_gk = tan(theta_gk) * cos(phi_gk);
  yp_gk = tan(theta_gk) * sin(phi_gk);

  theta_spec_cent = 12.5*deg2rad;
  
  xp_k = tan(atan(xp_gk) + atan(xp_eg));
  yp_k = tan(atan(yp_gk) + atan(yp_eg));

  double eps, epsL;
  double ssq;
  eps = 1./(1.0+2.0*fabs(q*q)/(Qsq)*pow(tan(theta_ep/2.0),2.0));
  epsL = Qsq/omega/omega * eps;
  ssq = sqrt( -Qsq  + mt*(2*omega+mt) );
  
  double A=0.0, B=0.0, C=0.0, D=0.0;
  A = mhyp*mhyp + Qsq - mt*mt - 2.0*omega*mt;
  B = -(A-msca*msca)/2.0;
  C = q * cos(theta_gk);
  D = (mt+omega)*(mt+omega);
  double theta = theta_gk; // theta_gamma_K
  
  double pk1,pk2;
  double AA, BB, CC;

  AA = C*C - D;
  BB = 2.0*C*B;
  CC = B*B - D*msca*msca;
  pk1 = (-BB + sqrt(BB*BB-4.0*AA*CC))/(2*AA); // 
  pk2 = (-BB - sqrt(BB*BB-4.0*AA*CC))/(2*AA); // 

  pkaon = pk2;
  
  double pL;
  pL = sqrt(q*q + pk2*pk2 - 2.0*q*pk2*cos(theta));
  pLamb = pL;

  if(pkaon>0) return 1;
  else return 0;
  
}


double CalcKinematics_eek::CheckMM(){
  double e1,e2,e3;
  
  // --- Beam --- // 
  TVector3 vec1(0.0,0.0,pbeam);
  e1 = sqrt(me*me + vec1*vec1);
  
  // --- Scattered electron --- // 
  double tempx, tempy, tempz;
  tempz = pepri/sqrt(1.0+xp_eep*xp_eep+yp_eep*yp_eep);
  tempx = tempz * xp_eep;
  tempy = tempz * yp_eep;
  TVector3 vec2(tempx, tempy, tempz);
  e2 = sqrt(me*me + vec2*vec2);
  
  // --- Kaon+ --- // 
  tempz = pkaon/sqrt(1.0+xp_k*xp_k+yp_k*yp_k);
  tempx = tempz * xp_k;
  tempy = tempz * yp_k;
  TVector3 vec3(tempx, tempy, tempz);
  e3 = sqrt(mk*mk + vec3*vec3);
  
  mmcheck = pow(e1+mt-e2-e3,2.0) 
    - (vec1-vec2-vec3)*(vec1-vec2-vec3);
  mmcheck = sqrt(mmcheck);

  return mmcheck;
}

double CalcKinematics_eek::CheckMM2(){
  double e1,e2;
  
  // --- virtual photon --- //
  TVector3 vec1(0.0,0.0,q);
  e1 = omega;
  
  // --- Kaon+ --- // 
  double tempx, tempy, tempz;
  tempz = pkaon/sqrt(1.0+xp_gk*xp_gk+yp_gk*yp_gk);
  tempx = tempz * xp_gk;
  tempy = tempz * yp_gk;
  TVector3 vec2(tempx, tempy, tempz);
  e2 = sqrt(mk*mk + vec2*vec2);
  
  mmcheck2 = pow(e1+mt-e2,2.0) 
    - (vec1-vec2)*(vec1-vec2);
  mmcheck2 = sqrt(mmcheck2);

  return mmcheck2;
}
