#ifndef fsi_h
#define fsi_h 1
#include <TLorentzVector.h>
#include "define.h"
#include "Param.h"
#include <complex.h>
#include <TRandom.h>

const int nmax =200;
const int nmax2 =200;
const int nmax3 =4000;
class fsi{

  
public:
  fsi();
  ~fsi();
  void SetBranch(string ifrname);
  void SetParam(string ifpname);
  void SetMom();
  void SetHist();
  void JostParam();
  double Chi2(TH1D* hmm, TH1D* hmmFSI);
  void SetFermiMom(string ifpname);
  TLorentzVector CalcPn(TLorentzVector P2n);
  TLorentzVector CalcPn2(TLorentzVector P2n);
  TLorentzVector CalcPn_deu(TLorentzVector P2n);
  TLorentzVector CalcPn_test(TLorentzVector P2n);
  double PhaseShift(double q, int LL, int potential);
  
  double CalcQ(TLorentzVector Pn, TLorentzVector Pp, TLorentzVector Pv, TLorentzVector PK);
  double CalcQ2(TLorentzVector Pn, TLorentzVector Pp, TLorentzVector Pv, TLorentzVector PK);
  complex<double> vlkkp(double skp,double skpp,int lcall);
  complex<double> vofq(double q2);
  void NewRoot(string ofname);
  void InfluenceFactor(string pname, int l);
  void InfluenceFactorVscale(string pname, int l);
  double InfluenceFactor2(double kreol, int potential);
  TRandom random;
  //------- NewRoot    -------//
  TFile * ofr;
  TTree* Tnew;
  float fsiw;

  double fac0,fac1,fac2,fac3,fac1t,fac2t,fac3t,fac1s,fac2s,fac3s,ton1s,ton1t,ton1,ton2,ton2s,ton2t,ton3,ton3s,ton3t,toff1,toff1s,toff1t,toff2,toff2s,toff2t,toff3,toff3s,toff3t;

  double era1,era1s,era1t,era2,era2s,era2t,era3,era3s,era3t;
  
  double fac1_2,fac2_2,fac3_2,fac1_2t,fac2_2t,fac3_2t,fac1_2s,fac2_2s,fac3_2s;
  double fac1_test,fac2_test,fac3_test;
  double Prel,Prel2,Prel_test;
  double mm;
  double pL,pn,pn2;
  double ranth,ranph;
  double pLx,pLy,pLz,pnx,pny,pnz,pnx2,pny2,pnz2;
  double theta_L;
  
  //-------- SetFermiMom ------//
  double  pval[nmax3],mprob[nmax3],dmprob[nmax3];
  int nump;
  //-------- SetParam ---------//
  double krel1[nmax3];
  double krel2[nmax3];
  double krel3[nmax3];
  double w1[nmax3];
  double w2[nmax3];
  double w3[nmax3];

  double w1t[nmax3];
  double w2t[nmax3];
  double w3t[nmax3];
  
  int np1,np2,np3;
  
  //-------- SetHist ---------//

  TGraphErrors* gf0;
  TGraphErrors* gf1;
  TGraphErrors* gf2;
  TGraphErrors* gf3;
  TGraphErrors* gf1t;
  TGraphErrors* gf2t;
  TGraphErrors* gf3t;
  TGraphErrors* gf1_all;
  TGraphErrors* gf2_all;
  TGraphErrors* gf3_all;  
  TGraphErrors* gxf1;
  TGraphErrors* gxf2;
  TGraphErrors* gxf3;
  TGraphErrors* gx1;
  TGraphErrors* gx2;
  TGraphErrors* gx3;

  TGraphErrors* gI0;
  TGraphErrors* gI1;
  TGraphErrors* gI2; 
  TGraphErrors* gI3; 
  TGraphErrors* gI[100];
  TGraphErrors* gIs[100];
  TGraphErrors* gIt[100];
  
  TGraphErrors* gI00;
  TGraphErrors* gI01;
  TGraphErrors* gI02; 
  TGraphErrors* gI03;   


  TGraphErrors* gERA_I0;
  TGraphErrors* gERA_I1;
  TGraphErrors* gERA_I2;
  TGraphErrors* gERA_I3;

  TGraphErrors* gJost_I0;
  TGraphErrors* gJost_I1;
  TGraphErrors* gJost_I2;
  TGraphErrors* gJost_I3;
  
  
  TGraphErrors* grad1;
  TGraphErrors* grad2; 
  TGraphErrors* grad3;     

  TGraphErrors* grad00;
  TGraphErrors* grad01;
  TGraphErrors* grad02; 
  TGraphErrors* grad03;

  TGraphErrors* gradI1;
  TGraphErrors* gradI2; 
  TGraphErrors* gradI3;
  
  
  TGraphErrors* gradl1[20];
  TGraphErrors* gradl2[20];
  TGraphErrors* gradl3[20];

  TGraphErrors* gIl1[20];
  TGraphErrors* gIl2[20];
  TGraphErrors* gIl3[20];

  TGraphErrors* gIl01[20];
  TGraphErrors* gIl02[20];
  TGraphErrors* gIl03[20];  
  
  TGraphErrors* gfr1[20];
  TGraphErrors* gfr2[20];
  TGraphErrors* gfr3[20];

  TGraphErrors* gfi1[20];
  TGraphErrors* gfi2[20];
  TGraphErrors* gfi3[20];
  

  TGraphErrors* gfth1;
  TGraphErrors* gfth2;
  TGraphErrors* gfth3;



  
  // T matrix calc
  
  TGraphErrors* gIton1s;
  TGraphErrors* gIton2s;
  TGraphErrors* gIton3s;

  TGraphErrors* gIton1t;
  TGraphErrors* gIton2t;
  TGraphErrors* gIton3t;

  TGraphErrors* gIton1;
  TGraphErrors* gIton2;
  TGraphErrors* gIton3;  


  TGraphErrors* gItoff1s;
  TGraphErrors* gItoff2s;
  TGraphErrors* gItoff3s;

  TGraphErrors* gItoff1t;
  TGraphErrors* gItoff2t;
  TGraphErrors* gItoff3t;

  TGraphErrors* gItoff1;
  TGraphErrors* gItoff2;
  TGraphErrors* gItoff3;    


  // ERA calc.
  
  TGraphErrors* gIera1;
  TGraphErrors* gIera2;
  TGraphErrors* gIera3;  
  
  TGraphErrors* gIera1s;
  TGraphErrors* gIera2s;
  TGraphErrors* gIera3s;  
  
  TGraphErrors* gIera1t;
  TGraphErrors* gIera2t;
  TGraphErrors* gIera3t;  
  
  
  
  TF1* fJl[100];
  TF1* fJls[100];
  TF1* fJlt[100];
  TH1F* hmm1;
  TH1F* hmm_fsi1;
  TH1F* hmm_fsi2;
  TH1F* hmm_fsi3;

  TH1F* hmm_fsi1t;
  TH1F* hmm_fsi2t;
  TH1F* hmm_fsi3t;

  TH1F* hmm_fsi1s;
  TH1F* hmm_fsi2s;
  TH1F* hmm_fsi3s;
  
  
  TH1F* hmm_fsi1_2;
  TH1F* hmm_fsi2_2;
  TH1F* hmm_fsi3_2;
  
  TH1F* hmm_fsi1_2s;
  TH1F* hmm_fsi2_2s;
  TH1F* hmm_fsi3_2s;

  TH1F* hmm_fsi1_2t;
  TH1F* hmm_fsi2_2t;
  TH1F* hmm_fsi3_2t;

  TH1F* hmm;
  TH1F* hmm_L;
  TH1F* hmm_nnL;


  TH1F* hmm_Jl[100];
  TH1F* hmm_Jls[100];
  TH1F* hmm_Jlt[100];

  TH1F* hmm_Jl_2[100];
  TH1F* hmm_Jl_2s[100];
  TH1F* hmm_Jl_2t[100];  
  //-------- JostParam ---------//
  string vname[100];
  double r_s[100],a_s[100],r_t[100],a_t[100];
  int vmax=10;
  //-------- SetMom   --------//

  

  int ENum;
  float vertex_p_P,vertex_p_xptar,vertex_p_yptar;
  float vertex_e_P,vertex_e_xptar,vertex_e_yptar;
  float vertex_p_Px,vertex_p_Py,vertex_p_Pz;
  float vertex_e_Px,vertex_e_Py,vertex_e_Pz;
  float Lth_gen, Lph_gen, Lp_gen;
  float Rth_gen, Rph_gen, Rp_gen;
  float vertex_uq_x,vertex_uq_y,vertex_uq_z,vertex_q,vertex_nu;
  float pferx,pfery,pferz,pfer,Mrec,efer;
  float mmnuc,Em,mm_L,mm_nnL;
  float Trec;
  //  TFile * ifr;
  TChain * T;
  
  //------- Influence Factor  --//
  ofstream* ofp0;
  ofstream* ofp1;
  ofstream* ofp2;
  ofstream* ofp3;
  ofstream* ofp1t;
  ofstream* ofp2t;
  ofstream* ofp3t;
  ofstream* ofp[100];
  ofstream* ofps[100];
  ofstream* ofpt[100];
  //------- PhaseShift ------//
  double ampj,ampj2,amtag,amtag2,fmu,fmunn,amu,fnucl,fnucl2,hbarc2;
  double skz,skz2,eon,eoff,rho,ekpj,ektag,lpjmax,skp,skpp,skp2,skpp2,vv1,uofq,fthecm,q2;
  double POL[nmax2],deltal[60];
  double skk[nmax2],x[nmax2],w[nmax2],wt[nmax2],pol[nmax2+1],pols[nmax2+1],xvq[nmax2],wvq[nmax2];
  complex<double> gren[nmax2+1][nmax2+1],ul[nmax2+1][nmax2+1],gren2[nmax2+1][nmax2+1],v[nmax2],v0[nmax2],ton[60],ton0[60],tborn[60],delta_rad[60],jon[60],toff[nmax3];
  complex<double> fthe[nmax2],fthe1[nmax2],fthe0[nmax2],Rl[nmax2+1],Rl0[nmax2+1],V[nmax2+1],fthel[20],ftheL;
  complex<double>Psi[180],Psi0[180],psi[180],psi0[180];
  complex<double> xi, vq,v2,vv,delta,det,ron,ron0,fborn,fborn0,fthed;
  double va,b_a,vr,b_r,b_r2,b_a2,vqa,vqr,ronr,tonre,tonim,tbornre,tbornim,perre,perim,delrad,delrad0,deldeg,deldega,tcross,cross[180],tcross1,tcross0,tcrossl;
  complex<double> expL[20];
  double I0,I1,Il[20],Il0[20],IJ,IERA,Iton,Itoff;
  double f2[180];
  TF1* fI0;
  const int n1 =40;
  const int n2 =40;
  int model =0;
  const int L_value =10;
  double xxx, totalc;
  int npot =n1+n2;
  int nrp1,lmax,lmax1;
  double tol =0.01;
  bool nr_mode = true;

};

fsi::fsi(){};
fsi::~fsi(){};

#endif
