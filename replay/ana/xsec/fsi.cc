#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
using namespace std;
#include <TApplication.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>
#include <TCut.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>
#include <TText.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TColor.h>
#include <TPaveText.h>
#include <TRandom.h>
#include <TMinuit.h>
#include "fsi.h"
#include "define.h"
#include "Param.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_specfunc.h>
using Eigen::MatrixXd;
using namespace Eigen;
extern void PL(double X, int L ,double * POL);
extern void PL2(double X, int L ,double * POL);
extern double Bessel(int l, double x);
extern void DGauss(double* Y ,double *WY, int N);
extern complex<double>LEQ1(complex<double>* A, complex<double>* B,int N,int LA);
extern complex<double>LEQ2(complex<double>* A, complex<double>* B, complex<double>* R,int N,int LA);
extern void LEQ3(complex<double>* A, complex<double>* B, complex<double>* R, int N);
extern double Jl0(double p0,double a0, double r0);
extern double fJ(double* x, double* par);
extern double fJ2(double* x, double* par);
extern double rad_ERA(double*x, double* par);
bool SIMC    = true;
bool TClone  = false;
//bool E09     = false;
bool Deu     = false;
//bool Deu     = true;
bool E09     = true;
bool Lam     = true;
bool single  = true;
//bool Vscale  = true;
bool Vscale  = false;
//bool Cha_mode =true;
bool Cha_mode =false;
bool Ton_mode =false;
/////////////////////////////////////////////////////////////////////////

void fsi::NewRoot(string ofname){

  ofr = new TFile(ofname.c_str(),"recreate");

  Tnew = new TTree("T","FSI calculation");
  
  if(TClone){
    Tnew = T->CloneTree(0);
    Tnew -> Branch("mm",&mm,"mm/D");
    Tnew -> Branch("fac1",&fac1,"fac1/D");
    Tnew -> Branch("fac2",&fac2,"fac2/D");
    Tnew -> Branch("fac3",&fac3,"fac3/D");
    Tnew -> Branch("fac1s",&fac1s,"fac1s/D");
    Tnew -> Branch("fac2s",&fac2s,"fac2s/D");
    Tnew -> Branch("fac3s",&fac3s,"fac3s/D");        
    Tnew -> Branch("fac1t",&fac1t,"fac1t/D");
    Tnew -> Branch("fac2t",&fac2t,"fac2t/D");
    Tnew -> Branch("fac3t",&fac3t,"fac3t/D");    
    Tnew -> Branch("fac1_2",&fac1_2,"fac1_2/D");
    Tnew -> Branch("fac2_2",&fac2_2,"fac2_2/D");
    Tnew -> Branch("fac3_2",&fac3_2,"fac3_2/D");
    Tnew -> Branch("fac1_2s",&fac1_2s,"fac1_2s/D");
    Tnew -> Branch("fac2_2s",&fac2_2s,"fac2_2s/D");
    Tnew -> Branch("fac3_2s",&fac3_2s,"fac3_2s/D");    
    Tnew -> Branch("fac1_2t",&fac1_2t,"fac1_2t/D");
    Tnew -> Branch("fac2_2t",&fac2_2t,"fac2_2t/D");
    Tnew -> Branch("fac3_2t",&fac3_2t,"fac3_2t/D");
    
    Tnew -> Branch("Perl",&Prel,"Prel/D");
    Tnew -> Branch("Perl2",&Prel2,"Prel2/D");
    Tnew -> Branch("pL",&pL,"pL/D");
    Tnew -> Branch("pn",&pn,"pn/D");
    Tnew -> Branch("pn2",&pn2,"pn2/D");
    Tnew -> Branch("ranth",&ranth,"ranth/D");
    Tnew -> Branch("ranph",&ranph,"ranph/D");
    Tnew -> Branch("pLx",&pLx,"pLx/D");
    Tnew -> Branch("pLy",&pLy,"pLy/D");
    Tnew -> Branch("pLz",&pLz,"pLz/D");
    Tnew -> Branch("pnx",&pnx,"pnx/D");
    Tnew -> Branch("pny",&pny,"pny/D");
    Tnew -> Branch("pnz",&pnz,"pnz/D");
    Tnew -> Branch("pnx2",&pnx2,"pnx2/D");
    Tnew -> Branch("pny2",&pny2,"pny2/D");
    Tnew -> Branch("pnz2",&pnz2,"pnz2/D");    
    Tnew -> Branch("theta_L",&theta_L,"theta_L/D");    
					   
    
    
  }
  
}



/////////////////////////////////////////////////////////////////////////

void fsi::SetParam(string ifpname){

  cout<<"=================================="<<endl;
  cout<<"=== Set Influence Param =========="<<endl;
  cout<<"=================================="<<endl;

  string buf;
  int s=0;
  string vnames[10];
  ifstream ifp(ifpname.c_str(),ios::in);
  if(ifp.fail()){cout<<"Could not find Files "<<ifpname<<endl;exit(1);}

  while(!ifp.eof()){
    getline(ifp,buf);
    if(buf[0]=='#')continue;
    stringstream sbuf(buf);
    sbuf >> vnames[s];
    cout<<" file : "<<vnames[s]<<endl;
    s++;
  }

  
  ifstream ifp1(vnames[0].c_str(),ios::in);
  ifstream ifp2(vnames[1].c_str(),ios::in);
  ifstream ifp3(vnames[2].c_str(),ios::in);


  int i1=0;
  
  while(!ifp1.eof()){
    getline(ifp1,buf);
    if(buf[0]=='#' || buf.length()==0 )continue;
    stringstream sbuf(buf);
    sbuf >> krel1[i1] >> w1[i1] >> w1t[i1];
    i1++;
  }


  np1=i1;
  
  int i2=0;
  while(!ifp2.eof()){
    getline(ifp2,buf);
    stringstream sbuf(buf);
    if(buf[0]=='#' ||  buf.length()==0 )continue;
    sbuf >> krel2[i2] >> w2[i2] >> w2t[i2];
    //    cout<<"i "<<i2<<" krel "<<krel2[i2]<<" w2 "<<w2[i2]<<endl;
    i2++;
  }

  np2 = i2;
  
  int i3=0;
  while(!ifp3.eof()){
    getline(ifp3,buf);
    if(buf[0]=='#' || buf.length()==0 )continue;
    stringstream sbuf(buf);
    sbuf >> krel3[i3] >> w3[i3] >> w3t[i3];
    i3++;
  }
  
  np3 =i3;


  
  
}

/////////////////////////////////////////////////////////////////////////

void fsi::JostParam(){

  vname[0] = "E09";     a_s[0] = -2.68; r_s[0] = 2.91; a_t[0] = -1.66; r_t[0] = 3.33;
  vname[1] = "Verma";   a_s[1] = -2.29; r_s[1] = 3.15; a_t[1] = -1.77; r_t[1] = 3.25;
  vname[2] = "Jue_A";   a_s[2] = -1.60; r_s[2] = 1.33; a_t[2] = -1.6; r_t[2] = 3.15;
  vname[3] = "Jue_B";   a_s[3] = -0.57; r_s[3] = 7.65; a_t[3] = -1.94; r_t[3] = 2.42;
  vname[4] = "NSC97a";  a_s[4] = -0.77; r_s[4] = 6.09; a_t[4] = -2.15; r_t[4] = 2.71;
  vname[5] = "NSC97b";  a_s[5] = -0.97; r_s[5] = 5.09; a_t[5] = -2.09; r_t[5] = 2.80;
  vname[6] = "NSC97c";  a_s[6] = -1.28; r_s[6] = 4.22; a_t[6] = -2.07; r_t[6] = 2.86;
  vname[7] = "NSC97d";  a_s[7] = -1.82; r_s[7] = 3.52; a_t[7] = -1.94; r_t[7] = 3.01;
  vname[8] = "NSC97e";  a_s[8] = -2.24; r_s[8] = 3.24; a_t[8] = -1.83; r_t[8] = 3.14;
  vname[9] = "NSC97f";  a_s[9] = -2.68; r_s[9] = 3.07; a_t[9] = -1.67; r_t[9] = 3.34;
  

  vmax=10;
  
}

/////////////////////////////////////////////////////////////////////////

void fsi::SetFermiMom(string ifpname){

  ifstream ifp(ifpname.c_str(),ios::in);
  if(ifp.fail())cout<<"Not found files "<<ifpname<<endl;
  string buf;

  //Initialization
  for(int  i=0;i<nmax3;i++){
    pval[i] = 0.0;
    mprob[i] = 0.0;
    dmprob[i] = 0.0;
  }
  

  // ==== Get pval & mprob ========//
  int ii=0;
  double order ;
  string p,pp;
  
  while(!ifp.eof()){
  getline(ifp,buf);
  if(buf.length()==0)break;
  stringstream sbuf(buf);
  sbuf >> pval[ii] >> dmprob[ii] >> mprob[ii];
  pval[ii] /=1000.; // MeV to GeV
  //  cout<<"i "<<ii<<" pval "<<pval[ii]<<" dmprob "<<dmprob[ii]<<" mprob "<<mprob[ii]<<endl;
  nump=ii;
  ii++;
  }
  
}


/////////////////////////////////////////////////////////////////////////

void fsi::InfluenceFactor(string pname, int ll){
  
  int ii=0;

  //  fI0 = new TF1("fI0","[0]/x*",0,);
  bool test =true;
  //  test =false;
  string pname0 = pname + "_3He.dat";
  string pname1 = pname + "_Verma.dat";
  string pname2 = pname + "_JulichA.dat";
  string pname3 = pname + "_JulichB.dat";
  
  //  string pname1t = pname + "_Verma_T.dat";
  //  string pname2t = pname + "_JulichA_T.dat";
  //  string pname3t = pname + "_JulichB_T.dat";

  ofp0 = new ofstream(pname0.c_str());
  ofp1 = new ofstream(pname1.c_str());
  ofp2 = new ofstream(pname2.c_str());
  ofp3 = new ofstream(pname3.c_str());

  //  ofp1t = new ofstream(pname1t.c_str());
  //  ofp2t = new ofstream(pname2t.c_str());
  //  ofp3t = new ofstream(pname3t.c_str());  


  *ofp0 << "### FSI tabel with 3He Potential (2-Gauss potential) #####"<<endl;

  
  *ofp1 << "### FSI tabel with Verma Potential (2-Gauss potential) #####"<<endl;
  *ofp1<<"# krel [MeV] # Influence Factor(S) # Influence Factor(T)"<<endl;

  *ofp2 << "### FSI tabel with Julich A Potential (2-Gauss potential) #####"<<endl;
  *ofp2<<"# krel [MeV] # Influence Factor(S) # Influence Factor(T)"<<endl;
  //  *ofp2<<"# Weight # krel [MeV] "<<endl;

  *ofp3 << "### FSI tabel with Julich B Potential (2-Gauss potential) #####"<<endl;
  *ofp3<<"# krel [MeV] # Influence Factor(S) # Influence Factor(T)"<<endl;
  //  *ofp3<<"# Weight # krel [MeV] "<<endl;


  

  //  *ofp1t << "### FSI tabel with Verma Potential Triplet (2-Gauss potential) #####"<<endl;
  //  *ofp1t<<"# Weight # krel [MeV] "<<endl;

  //  *ofp2t << "### FSI tabel with Julich A Potential Triplet (2-Gauss potential) #####"<<endl;
  //  *ofp2t<<"# Weight # krel [MeV] "<<endl;

  //  *ofp3t << "### FSI tabel with Julich B Potential Triplet (2-Gauss potential) #####"<<endl;
  //  *ofp3t<<"# Weight # krel [MeV] "<<endl;  

  
  double qmin = 0.0;
  double qmax = 2000;
  int imax = 500;
  //  double qmax = 300;
  //  int imax = 300;
  if(test){

    qmax =250;
    imax = 25;
  }
  
  for(int i=1; i<imax;i++){
  
    double qi = (qmax - qmin)/(double)imax * (double)i +qmin;

    
    fac0   = PhaseShift(qi, ll, 0);
    gI00   -> SetPoint(ii,qi,I0);
    gI0    -> SetPoint(ii,qi,I1);
    grad00 -> SetPoint(ii,qi,delrad0);
    fac0   = I0;
    
    fac1   = PhaseShift(qi, ll, 1);
    gx1    -> SetPoint(ii,qi,tcross0);
    gxf1   -> SetPoint(ii,qi,tcross);
    gI01   -> SetPoint(ii,qi,I0);
    gI1    -> SetPoint(ii,qi,I1);
    grad1  -> SetPoint(ii,qi,delrad);
    grad01 -> SetPoint(ii,qi,delrad0);
    gradI1 -> SetPoint(ii,delrad0,fac1);

    ton1s  = Iton;
    toff1s = Itoff;
    era1s  = IERA;
    
    for(int l=0;l<=ll;l++){
      gradl1[l] -> SetPoint(ii,qi,deltal[l]);
      gIl1[l]   -> SetPoint(ii,qi,Il[l]);
      gIl01[l]  -> SetPoint(ii,qi,Il0[l]);
      gfr1[l]   -> SetPoint(ii,qi,real(fthel[l]));
      gfi1[l]   -> SetPoint(ii,qi,imag(fthel[l]));
    }


    
    if(ii==0){
      for(int th =0;th<180;th++)
	gfth1->SetPoint(th,th,cross[th]);
    }


    fac2 = PhaseShift(qi, ll, 2);
    gx2    ->SetPoint(ii,qi,tcross0);
    gxf2   ->SetPoint(ii,qi,tcross);
    gI02   ->SetPoint(ii,qi,I0);
    gI2    ->SetPoint(ii,qi,I1);
    grad2  ->SetPoint(ii,qi,delrad);
    grad02 ->SetPoint(ii,qi,delrad0);
    gradI2 ->SetPoint(ii,delrad0,fac2);

    ton2s = Iton;
    toff2s = Itoff;
    era2s  = IERA;
  
    for(int l=0;l<=ll;l++){
      gradl2[l] -> SetPoint(ii,qi,deltal[l]);
      gIl2[l]   -> SetPoint(ii,qi,Il[l]);
      gIl02[l]  -> SetPoint(ii,qi,Il0[l]);
      gfr2[l]   -> SetPoint(ii,qi,real(fthel[l]));
      gfi2[l]   -> SetPoint(ii,qi,imag(fthel[l]));
    }

    if(ii==0){
      for(int th =0;th<180;th++)
	gfth2->SetPoint(th,th,cross[th]);
    }

	
    fac3 = PhaseShift(qi, ll, 3);
    gx3    -> SetPoint(ii,qi,tcross0);
    gxf3   -> SetPoint(ii,qi,tcross);    
    gI03   -> SetPoint(ii,qi,I0);
    gI3    -> SetPoint(ii,qi,I1);
    grad3  -> SetPoint(ii,qi,delrad);
    grad03 -> SetPoint(ii,qi,delrad0);    
    gradI3 -> SetPoint(ii,delrad0,fac3);

    ton3s  = Iton;
    toff3s = Itoff;
    era3s  = IERA;


    
    for(int l=0;l<=ll;l++){
      gradl3[l] -> SetPoint(ii,qi,deltal[l]);
      gIl3[l]   -> SetPoint(ii,qi,Il[l]);
      gIl03[l]  -> SetPoint(ii,qi,Il0[l]);
      gfr3[l]   -> SetPoint(ii,qi,real(fthel[l]));
      gfi3[l]   -> SetPoint(ii,qi,imag(fthel[l]));
    }

    if(ii==0){
      for(int th =0;th<180;th++)
	gfth3->SetPoint(th,th,cross[th]);
    }
      //      gfth3->SetPoint(th,th,real( fthe[th]*conj(fthe[th]) ));
    
    if(fac1<0 || fac1>100.) fac1 = 1.0;
    if(fac2<0 || fac2>100.) fac2 = 1.0;
    if(fac3<0 || fac3>100.) fac3 = 1.0;


    fac1t  = PhaseShift(qi, ll, -1);
    ton1t  = Iton;
    toff1t = Itoff;
    era1t  = IERA;
    
    fac2t  = PhaseShift(qi, ll, -2);
    ton2t  = Iton;
    toff2t = Itoff;
    era2t  = IERA;
    
    fac3t  = PhaseShift(qi, ll, -3);
    ton3t  = Iton;
    toff3t = Itoff;
    era3t  = IERA;

    
    //    fac1t = (fac1 + 3*fac1t)/4.;
    //    fac2t = (fac2 + 3*fac2t)/4.;
    //    fac3t = (fac3 + 3*fac3t)/4.;	
    
    ton1 = (ton1s + 3.0*ton1t)/4.0;
    ton2 = (ton2s + 3.0*ton2t)/4.0;
    ton3 = (ton3s + 3.0*ton3t)/4.0;

    toff1 = (toff1s + 3.0*toff1t)/4.0;
    toff2 = (toff2s + 3.0*toff2t)/4.0;
    toff3 = (toff3s + 3.0*toff3t)/4.0;

    era1 = (era1s + 3.0*era1t)/4.0;
    era2 = (era2s + 3.0*era2t)/4.0;
    era3 = (era3s + 3.0*era3t)/4.0;


    //==== T-operater on-shell =======//
    
    gIton1s-> SetPoint(ii,qi,ton1s);
    gIton2s-> SetPoint(ii,qi,ton2s);
    gIton3s-> SetPoint(ii,qi,ton3s);
    
    gIton1t-> SetPoint(ii,qi,ton1t);
    gIton2t-> SetPoint(ii,qi,ton2t);
    gIton3t-> SetPoint(ii,qi,ton3t);

    gIton1-> SetPoint(ii,qi,ton1);
    gIton2-> SetPoint(ii,qi,ton2);
    gIton3-> SetPoint(ii,qi,ton3);


    //==== T-operater half off-shell =======//    

    gItoff1s-> SetPoint(ii,qi,toff1s);
    gItoff2s-> SetPoint(ii,qi,toff2s);
    gItoff3s-> SetPoint(ii,qi,toff3s);

    gItoff1t-> SetPoint(ii,qi,toff1t);
    gItoff2t-> SetPoint(ii,qi,toff2t);
    gItoff3t-> SetPoint(ii,qi,toff3t);

    gItoff1 -> SetPoint(ii,qi,toff1);
    gItoff2 -> SetPoint(ii,qi,toff2);
    gItoff3 -> SetPoint(ii,qi,toff3);


    //==== ERA caluclaiton =======//    

    gIera1s-> SetPoint(ii,qi,era1s);
    gIera2s-> SetPoint(ii,qi,era2s);
    gIera3s-> SetPoint(ii,qi,era3s);

    gIera1t-> SetPoint(ii,qi,era1t);
    gIera2t-> SetPoint(ii,qi,era2t);
    gIera3t-> SetPoint(ii,qi,era3t);

    gIera1 -> SetPoint(ii,qi,era1);
    gIera2 -> SetPoint(ii,qi,era2);
    gIera3 -> SetPoint(ii,qi,era3);

    
    //====== Influence Factor ===========//
    
    gf1 -> SetPoint(ii,qi,fac1);
    gf2 -> SetPoint(ii,qi,fac2);
    gf3 -> SetPoint(ii,qi,fac3);

    gf1t -> SetPoint(ii,qi,fac1t);
    gf2t -> SetPoint(ii,qi,fac2t);
    gf3t -> SetPoint(ii,qi,fac3t);    

    gf1_all -> SetPoint(ii,qi,(fac1 + 3.0* fac1t)/4.0);
    gf2_all -> SetPoint(ii,qi,(fac2 + 3.0* fac1t)/4.0);
    gf3_all -> SetPoint(ii,qi,(fac3 + 3.0* fac3t)/4.0);        

    
    *ofp0 << qi << " " << fac0 << endl;
    *ofp1 << qi << " " << fac1 <<" "<< fac1t <<endl;
    *ofp2 << qi << " " << fac2 <<" "<< fac2t <<endl;
    *ofp3 << qi << " " << fac3 <<" "<< fac3t <<endl;

    
    
    //   *ofp0 << qi << " " << fac0 <<endl;
    //   *ofp1 << qi << " " << fac1 <<endl;
    //   *ofp2 << qi << " " << fac2 <<endl;
    //   *ofp3 << qi << " " << fac3 <<endl;
    //   *ofp1t << qi << " " << fac1t <<endl;
    //   *ofp2t << qi << " " << fac2t <<endl;
    //   *ofp3t << qi << " " << fac3t <<endl;

    
    ii++;
  }


  
  for(int i=0;i<vmax;i++){
    fJl[i]->Write();
    fJls[i]->Write();
    fJlt[i]->Write();
  }
  
  for(int l=0;l<=ll;l++){
    gradl1[l]->Write();
    gradl2[l]->Write();
    gradl3[l]->Write();
    gIl1[l]  -> Write();
    gIl2[l]  -> Write();
    gIl3[l]  -> Write();
    gIl01[l] -> Write();
    gIl02[l] -> Write();
    gIl03[l] -> Write();    
    gfr1[l]  -> Write();
    gfi1[l]  -> Write();
    gfr2[l]  -> Write();
    gfi2[l]  -> Write();
    gfr3[l]  -> Write();
    gfi3[l]  -> Write();    
  }


  grad_Cha_1s->Write();
  grad_Cha_2s->Write();
  grad_Cha_3s->Write();
  frad_ERA_1s->Write();
  frad_ERA_2s->Write();
  frad_ERA_3s->Write();
  
  gf0->Write();
  gf1->Write();
  gf2->Write();
  gf3->Write();
  gf1t->Write();
  gf2t->Write();
  gf3t->Write();
  gf1_all->Write();
  gf2_all->Write();
  gf3_all->Write();    
  //  gxf1->Write();
  //  gxf2->Write();
  //  gxf3->Write();
  //  gx1->Write();
  //  gx2->Write();
  //  gx3->Write();
  gI0->Write();
  gI1->Write();
  gI2->Write();
  gI3->Write();
  gI00->Write();
  gI01->Write();
  gI02->Write();
  gI03->Write();
  grad1->Write();
  grad2->Write();
  grad3->Write();
  grad00->Write();
  grad01->Write();
  grad02->Write();
  grad03->Write();  
  gradI1->Write();
  gradI2->Write();
  gradI3->Write();
  gfth1->Write();
  gfth2->Write();
  gfth3->Write();
  gIton1s->Write();
  gIton1t->Write();
  gIton1 ->Write();
  gIton2s->Write();
  gIton2t->Write();
  gIton2 ->Write();
  gIton3s->Write();
  gIton3t->Write();
  gIton3 ->Write();
  gItoff1s->Write();
  gItoff1t->Write();
  gItoff1 ->Write();
  gItoff2s->Write();
  gItoff2t->Write();
  gItoff2 ->Write();
  gItoff3s->Write();
  gItoff3t->Write();
  gItoff3 ->Write();  

  gIera1s->Write();
  gIera1t->Write();
  gIera1 ->Write();
  gIera2s->Write();
  gIera2t->Write();
  gIera2 ->Write();
  gIera3s->Write();
  gIera3t->Write();
  gIera3 ->Write();


  
  ofp0->close();
  ofp1->close();
  ofp2->close();
  ofp3->close();

  //  ofp1t->close();
  //  ofp2t->close();
  //  ofp3t->close();  
  
  
}


////////////////////////////////////////////////////////////

void fsi::InfluenceFactorVscale(string pname, int ll){


  cout<<"===================================================="<<endl;
  cout<<"======< InfluenceFacroe Vscale Calculation >======="<<endl;
  cout<<"===================================================="<<endl;
  
  int ii=0;

  bool test =true;
  //  test =false;

  string pnames[100];
  int nv =30;
  double ss;

   
  double qmin = 0.0;
  double qmax = 2000;
  int imax = 500;
  if(test){
    qmax =250;
    imax =125;
  }
  

  double va_s = -167.34; // Verma Va Singlet
  double va_t = -132.42; // Verma Va Triplet
  double b_as = 1.1;     // Verma beta Singlet
  double b_at = 1.1;     // Verma beta Triplet
  double re  = 2.91 ;   // Verma effective range
  double fac,fac_s,fac_t;
  
  double a_s = -2.29; // Verma scattering range Singlet
  double r_s =  3.15; // Verma effective  range Singlet
  double a_t = -1.77; // Verma scattering range Triplet
  double r_t =  3.25; // Verma effective  range Singlet
  
  for(int n =0;n<nv;n++){

    ss = (1.0 + (double)n/(double)nv);

    cout<<" n : "<<n<<" scale fac "<<ss<<endl;
    pnames[n] = pname + Form("_scale%d.dat",n);
    ofp[n] = new ofstream(pnames[n].c_str());
    *ofp[n] << Form("### FSI tabel with V scale fac(%4f) (2-Gauss potential) #####",ss)<<endl;
    gI[n]= new TGraphErrors();
    gIs[n]= new TGraphErrors();
    gIt[n]= new TGraphErrors();
    gI[n]->SetName(Form("gI_%d",n));
    gIs[n]->SetName(Form("gIqs_%d",n));
    gIt[n]->SetName(Form("gIt_%d",n));

    
  for(int i=1; i<imax;i++){
    double qi = (qmax - qmin)/(double)imax * (double)i +qmin;

    // Singlet  Calculation //
    v_a = va_s * ss;   beta_a = b_as; 
    a_0 = a_s;         r_0 = r_s;
    fac_s   = PhaseShift(qi, ll, 0);
    
    if(fac_s<0 || fac_s>100.) fac_s = 1.0;
    gIs[n]    -> SetPoint(i-1,qi,fac_s);

    // Triplet Calculatuion //

    v_a = va_t * ss;    beta_a = b_at; 
    a_0 = a_t;         r_0 = r_t;
    fac_t   = PhaseShift(qi, ll, 0);
   
    if(fac_s<0 || fac_s>100.) fac_s = 1.0;
    gIt[n]    -> SetPoint(i-1,qi,fac_t);    
  
    // 0.25 * 1S0 + 0.75 * 3S1 //
    fac = (fac_s + 3.0 * fac_t)/4.0;
    gI[n]    -> SetPoint(i-1,qi,fac);    

    *ofp[n] << qi << " " << fac_s <<" "<< fac_t <<endl;
    
  } // end for i

  gI[n] ->Write();
  gIs[n]->Write();
  gIt[n]->Write();
  gIT[n]->Write();
  
  ofp[n]->close();

  
  } // end for n
  



    

  
  
}



///////////////////////////////////////////////////////////////////////

void fsi::SetHist(){

 
  gf0 = new TGraphErrors();
  gf0->SetName("gf0");
  gf0->SetTitle("Verma Potentail FSI Influence Factor ; P_{#Lambda-n} [MeV] ; FSI Influence factor");
  gf0->SetMarkerStyle(31);
  gf0->SetMarkerColor(8);
  
  gf1 = new TGraphErrors();
  gf1->SetName("gf1");
  gf1->SetTitle("Verma Potentail FSI Influence Factor ; P_{#Lambda-n} [MeV] ; FSI Influence factor");
  gf1->SetMarkerStyle(31);
  gf1->SetMarkerColor(1);
  gf2 = new TGraphErrors();
  gf2->SetName("gf2");
  gf2->SetTitle("Julich A Potentail FSI Influence Factor ; P_{#Lambda-n} [MeV] ; FSI Influence factor");
  gf2->SetMarkerStyle(31);
  gf2->SetMarkerColor(2);
  gf3 = new TGraphErrors();
  gf3->SetName("gf3");
  gf3->SetTitle("Julich B Potentail FSI Influence Factor ; P_{#Lambda-n} [MeV] ; FSI Influence factor");
  gf3->SetMarkerStyle(31);
  gf3->SetMarkerColor(4);



  gf1t = new TGraphErrors();
  gf1t->SetName("gf1t");
  gf1t->SetTitle("Verma Potentail FSI Influence Factor ; P_{#Lambda-n} [MeV] ; FSI Influence factor");
  gf1t->SetMarkerStyle(31);
  gf1t->SetMarkerColor(1);
  gf2t = new TGraphErrors();
  gf2t->SetName("gf2t");
  gf2t->SetTitle("Julich A Potentail FSI Influence Factor ; P_{#Lambda-n} [MeV] ; FSI Influence factor");
  gf2t->SetMarkerStyle(31);
  gf2t->SetMarkerColor(2);
  gf3t = new TGraphErrors();
  gf3t->SetName("gf3t");
  gf3t->SetTitle("Julich B Potentail FSI Influence Factor ; P_{#Lambda-n} [MeV] ; FSI Influence factor");
  gf3t->SetMarkerStyle(31);
  gf3t->SetMarkerColor(4);



  gf1_all = new TGraphErrors();
  gf1_all->SetName("gf1_all");
  gf1_all->SetTitle("Verma Potentail FSI Influence Factor ; P_{#Lambda-n} [MeV] ; FSI Influence factor");
  gf1_all->SetMarkerStyle(31);
  gf1_all->SetMarkerColor(1);
  gf2_all = new TGraphErrors();
  gf2_all->SetName("gf2_all");
  gf2_all->SetTitle("Julich A Potentail FSI Influence Factor ; P_{#Lambda-n} [MeV] ; FSI Influence factor");
  gf2_all->SetMarkerStyle(31);
  gf2_all->SetMarkerColor(2);
  gf3_all = new TGraphErrors();
  gf3_all->SetName("gf3_all");
  gf3_all->SetTitle("Julich B Potentail FSI Influence Factor ; P_{#Lambda-n} [MeV] ; FSI Influence factor");
  gf3_all->SetMarkerStyle(31);
  gf3_all->SetMarkerColor(4);  




  
  gxf1 = new TGraphErrors();
  gxf1->SetName("gxf1");
  gxf1->SetTitle("Verma Potentail w FSI Cross Section; P_{#Lambda-n} [MeV] ; Cross Section");
  gxf1->SetMarkerStyle(2);
  gxf1->SetMarkerColor(1);

  gxf2 = new TGraphErrors();
  gxf2->SetName("gxf2");
  gxf2->SetTitle("Julich A Potential w FSI Cross Section; P_{#Lambda-n} [MeV] ; Cross Section");
  gxf2->SetMarkerStyle(2);
  gxf2->SetMarkerColor(2);


  gxf3 = new TGraphErrors();
  gxf3->SetName("gxf3");
  gxf3->SetTitle("Julich B Potential w FSI Cross Section; P_{#Lambda-n} [MeV] ; Cross Section");
  gxf3->SetMarkerStyle(2);
  gxf3->SetMarkerColor(4);
  


  gx1 = new TGraphErrors();
  gx1->SetName("gx1");
  gx1->SetTitle("Verma Potentail w FSI Cross Section; P_{#Lambda-n} [MeV] ; Cross Section");
  gx1->SetMarkerStyle(4);
  gx1->SetMarkerColor(1);

  gx2 = new TGraphErrors();
  gx2->SetName("gx2");
  gx2->SetTitle("Julich A Potential w FSI Cross Section; P_{#Lambda-n} [MeV] ; Cross Section");
  gx2->SetMarkerStyle(4);
  gx2->SetMarkerColor(2);


  gx3 = new TGraphErrors();
  gx3->SetName("gx3");
  gx3->SetTitle("Julich B Potential w FSI Cross Section; P_{#Lambda-n} [MeV] ; Cross Section");
  gx3->SetMarkerStyle(4);
  gx3->SetMarkerColor(4);
  

  gI0 = new TGraphErrors();
  gI0->SetName("gI0");
  gI0->SetTitle("3He w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor");
  gI0->SetMarkerStyle(2);
  gI0->SetMarkerColor(8);  
  
  
  gI1 = new TGraphErrors();
  gI1->SetName("gI1");
  gI1->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor");
  gI1->SetMarkerStyle(2);
  gI1->SetMarkerColor(1);

  gI2 = new TGraphErrors();
  gI2->SetName("gI2");
  gI2->SetTitle("Julich A w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor");
  gI2->SetMarkerStyle(2);
  gI2->SetMarkerColor(2);

  gI3 = new TGraphErrors();
  gI3->SetName("gI3");
  gI3->SetTitle("Julich B w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor");
  gI3->SetMarkerStyle(2);
  gI3->SetMarkerColor(4);
 

  
  gI00 = new TGraphErrors();
  gI00->SetName("gI00");
  gI00->SetTitle("3He potential w FSI Influence Factor (Jost Func.); P_{#Lambda-n} [MeV] ; Influence Factor");
  gI00->SetMarkerStyle(4);
  gI00->SetMarkerColor(3);
  
  gI01 = new TGraphErrors();
  gI01->SetName("gI01");
  gI01->SetTitle("Verma Potentail w FSI Influence Factor (Jost Func.); P_{#Lambda-n} [MeV] ; Influence Factor");
  gI01->SetMarkerStyle(4);
  gI01->SetMarkerColor(1);

  gI02 = new TGraphErrors();
  gI02->SetName("gI02");
  gI02->SetTitle("Julich A w FSI Influence Factor (Jost Func.); P_{#Lambda-n} [MeV] ; Influence Factor");
  gI02->SetMarkerStyle(4);
  gI02->SetMarkerColor(2);

  gI03 = new TGraphErrors();
  gI03->SetName("gI03");
  gI03->SetTitle("Julich B w FSI Influence Factor (Jost Func.); P_{#Lambda-n} [MeV] ; Influence Factor");
  gI03->SetMarkerStyle(4);
  gI03->SetMarkerColor(4);


  

  grad1 = new TGraphErrors();
  grad1->SetName("grad1");
  grad1->SetTitle("Vema Potential w FSI Phase Shift; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
  grad1->SetMarkerStyle(3);
  grad1->SetMarkerColor(1);

  grad2 = new TGraphErrors();
  grad2->SetName("grad2");
  grad2->SetTitle("Julich A Potential w FSI Phase Shift; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
  grad2->SetMarkerStyle(3);
  grad2->SetMarkerColor(2);  

  grad3 = new TGraphErrors();
  grad3->SetName("grad3");
  grad3->SetTitle("Julich B Potential w FSI Phase Shift; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
  grad3->SetMarkerStyle(3);
  grad3->SetMarkerColor(4);    


  grad00 = new TGraphErrors();
  grad00->SetName("grad00");
  grad00->SetTitle("3He potential w FSI Phase Shift; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
  grad00->SetMarkerStyle(5);
  grad00->SetMarkerColor(3);

  
  grad01 = new TGraphErrors();
  grad01->SetName("grad01");
  grad01->SetTitle("Vema Potential w FSI Phase Shift; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
  grad01->SetMarkerStyle(5);
  grad01->SetMarkerColor(1);

  grad02 = new TGraphErrors();
  grad02->SetName("grad02");
  grad02->SetTitle("Julich A Potential w FSI Phase Shift; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
  grad2->SetMarkerStyle(5);
  grad2->SetMarkerColor(2);  

  grad03 = new TGraphErrors();
  grad03->SetName("grad03");
  grad03->SetTitle("Julich B Potential w FSI Phase Shift; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
  grad03->SetMarkerStyle(5);
  grad03->SetMarkerColor(4);    


  double prel[12]={10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0,120.0,140.0};
  double del_cha_1s[12]={0.11656,022035,0.31112,0.38399,0.43910,0.47819,0.50366,0.51790,0.52093,0.50058,0.46577};
  double del_cha_2s[12]={0.8084,0.15907,0.23309,0.30124,0.36257,0.41669,0.46362,0.50368,0.56531,0.60631,0.63124};
  double del_cha_3s[12]={0.2951,0.5728,0.8598,0.10854,0.12950,0.14601,0.16874,0,17445,0.16700,0.14875};
  
  grad_Cha_1s = new TGraphErrors();
  grad_Cha_1s->SetName("grad_Cha_1s");
  grad_Cha_1s->SetTitle("Vema Potential w FSI Phase Shift (Singlet) ; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
  grad_Cha_1s->SetMarkerStyle(8);
  grad_Cha_1s->SetMarkerColor(1);

  grad_Cha_2s = new TGraphErrors();
  grad_Cha_2s->SetName("grad_Cha_2s");
  grad_Cha_2s->SetTitle("Julich A Potential w FSI Phase Shift (Singlet) ; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
  grad_Cha_2s->SetMarkerStyle(8);
  grad_Cha_2s->SetMarkerColor(2);  

  grad_Cha_3s = new TGraphErrors();
  grad_Cha_3s->SetName("grad_Cha_3s");
  grad_Cha_3s->SetTitle("Julich B Potential w FSI Phase Shift (Singlet) ; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
  grad_Cha_3s->SetMarkerStyle(8);
  grad_Cha_3s->SetMarkerColor(4);  

  frad_ERA_1s =new TF1("frad_ERA_1s",rad_ERA,0,250,2);
  frad_ERA_1s->SetParameters(-2.29,3.15);
  frad_ERA_2s =new TF1("frad_ERA_2s",rad_ERA,0,250,2);
  frad_ERA_2s->SetParameters(-1.60,1.33);  
  frad_ERA_3s =new TF1("frad_ERA_3s",rad_ERA,0,250,2);
  frad_ERA_3s->SetParameters(-0.57,7.65);  

  
  
  for(int i=0;i<12;i++){
    grad_Cha_1s->SetPoint(i,prel[i],del_cha_1s[i]);
    grad_Cha_2s->SetPoint(i,prel[i],del_cha_2s[i]);
    grad_Cha_3s->SetPoint(i,prel[i],del_cha_3s[i]);
  }
  
  

  gradI1 = new TGraphErrors();
  gradI1->SetName("gradI1");
  gradI1->SetTitle("Vema Potential w FSI Phase Shift; phase shift [rad]; Influence Factor ");
  gradI1->SetMarkerStyle(3);
  gradI1->SetMarkerColor(1);

  gradI2 = new TGraphErrors();
  gradI2->SetName("gradI2");
  gradI2->SetTitle("Vema Potential w FSI Phase Shift; phase shift [rad]; Influence Factor ");  
  gradI2->SetMarkerStyle(3);
  gradI2->SetMarkerColor(2);  

  gradI3 = new TGraphErrors();
  gradI3->SetName("gradI3");
  gradI3->SetTitle("Julich A Potential w FSI Phase Shift; P_{#Lambda-n} [MeV] ; phase shift [rad] ");  

  gradI3->SetMarkerStyle(3);
  gradI3->SetMarkerColor(4);    


  gfth1 = new TGraphErrors();
  gfth1->SetName("gfth1");
  gfth1->SetTitle("Verma Potentail w FSI Influence Factor; Scattering Angle [rad] ; Influence Factor");
  gfth1->SetMarkerStyle(2);
  gfth1->SetMarkerColor(1);

  gfth2 = new TGraphErrors();
  gfth2->SetName("gfth2");
  gfth2->SetTitle("Verma Potentail w FSI Influence Factor; Scattering Angle [rad]  ; Influence Factor");
  gfth2->SetMarkerStyle(2);
  gfth2->SetMarkerColor(2);

  gfth3 = new TGraphErrors();
  gfth3->SetName("gfth3");
  gfth3->SetTitle("Verma Potentail w FSI Influence Factor;  Scattering Angle [rad]  ; Influence Factor");
  gfth3->SetMarkerStyle(2);
  gfth3->SetMarkerColor(4);
  

  gIton1 = new TGraphErrors();
  gIton1->SetName("gIton1");
  gIton1->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (0.75 ^{3}S_{1} + 0.25 ^{1}S_{0})");
  gIton1->SetMarkerStyle(7);
  gIton1->SetMarkerColor(1);

  gIton1s = new TGraphErrors();
  gIton1s->SetName("gIton1s");
  gIton1s->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{1}S_{0})");
  gIton1s->SetMarkerStyle(2);
  gIton1s->SetMarkerColor(1);

  gIton1t = new TGraphErrors();
  gIton1t->SetName("gIton1t");
  gIton1t->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{3}S_{1})");
  gIton1t->SetMarkerStyle(2);
  gIton1t->SetMarkerColor(1);    





  gIton2 = new TGraphErrors();
  gIton2->SetName("gIton2");
  gIton2->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (0.75 ^{3}S_{1} + 0.25 ^{1}S_{0})");
  gIton2->SetMarkerStyle(7);
  gIton2->SetMarkerColor(2);

  gIton2s = new TGraphErrors();
  gIton2s->SetName("gIton2s");
  gIton2s->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{1}S_{0})");
  gIton2s->SetMarkerStyle(2);
  gIton2s->SetMarkerColor(2);

  gIton2t = new TGraphErrors();
  gIton2t->SetName("gIton2t");
  gIton2t->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{3}S_{1})");
  gIton2t->SetMarkerStyle(2);
  gIton2t->SetMarkerColor(2);    

  



  gIton3 = new TGraphErrors();
  gIton3->SetName("gIton3");
  gIton3->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (0.75 ^{3}S_{1} + 0.25 ^{1}S_{0})");
  gIton3->SetMarkerStyle(7);
  gIton3->SetMarkerColor(3);

  gIton3s = new TGraphErrors();
  gIton3s->SetName("gIton3s");
  gIton3s->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{1}S_{0})");
  gIton3s->SetMarkerStyle(3);
  gIton3s->SetMarkerColor(3);

  gIton3t = new TGraphErrors();
  gIton3t->SetName("gIton3t");
  gIton3t->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{3}S_{1})");
  gIton3t->SetMarkerStyle(3);
  gIton3t->SetMarkerColor(3);    





  // half off-shell
  
  gItoff1 = new TGraphErrors();
  gItoff1->SetName("gItoff1");
  gItoff1->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (0.75 ^{3}S_{1} + 0.25 ^{1}S_{0})");
  gItoff1->SetMarkerStyle(3);
  gItoff1->SetMarkerColor(1);

  gItoff1s = new TGraphErrors();
  gItoff1s->SetName("gItoff1s");
  gItoff1s->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{1}S_{0})");
  gItoff1s->SetMarkerStyle(3);
  gItoff1s->SetMarkerColor(1);

  gItoff1t = new TGraphErrors();
  gItoff1t->SetName("gItoff1t");
  gItoff1t->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{3}S_{1})");
  gItoff1t->SetMarkerStyle(3);
  gItoff1t->SetMarkerColor(1);    


  gItoff2 = new TGraphErrors();
  gItoff2->SetName("gItoff2");
  gItoff2->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (0.75 ^{3}S_{1} + 0.25 ^{1}S_{0})");
  gItoff2->SetMarkerStyle(2);
  gItoff2->SetMarkerColor(2);  

  gItoff2s = new TGraphErrors();
  gItoff2s->SetName("gItoff2s");
  gItoff2s->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{1}S_{0})");
  gItoff2s->SetMarkerStyle(2);
  gItoff2s->SetMarkerColor(2);

  gItoff2t = new TGraphErrors();
  gItoff2t->SetName("gItoff2t");
  gItoff2t->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{3}S_{1})");
  gItoff2t->SetMarkerStyle(2);
  gItoff2t->SetMarkerColor(2);    

  
  gItoff3 = new TGraphErrors();
  gItoff3->SetName("gItoff3");
  gItoff3->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (0.75 ^{3}S_{1} + 0.25 ^{1}S_{0})");
  gItoff3->SetMarkerStyle(7);
  gItoff3->SetMarkerColor(3);

  gItoff3s = new TGraphErrors();
  gItoff3s->SetName("gItoff3s");
  gItoff3s->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{1}S_{0})");
  gItoff3s->SetMarkerStyle(3);
  gItoff3s->SetMarkerColor(3);

  gItoff3t = new TGraphErrors();
  gItoff3t->SetName("gItoff3t");
  gItoff3t->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{3}S_{1})");
  gItoff3t->SetMarkerStyle(3);
  gItoff3t->SetMarkerColor(3);    


  //======== Influence Factor ERA caluration ==========//


  gIera1 = new TGraphErrors();
  gIera1->SetName("gIera1");
  gIera1->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (0.75 ^{3}S_{1} + 0.25 ^{1}S_{0})");
  gIera1->SetMarkerStyle(7);
  gIera1->SetMarkerColor(1);

  gIera1s = new TGraphErrors();
  gIera1s->SetName("gIera1s");
  gIera1s->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{1}S_{0})");
  gIera1s->SetMarkerStyle(2);
  gIera1s->SetMarkerColor(1);

  gIera1t = new TGraphErrors();
  gIera1t->SetName("gIera1t");
  gIera1t->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{3}S_{1})");
  gIera1t->SetMarkerStyle(2);
  gIera1t->SetMarkerColor(1);    





  gIera2 = new TGraphErrors();
  gIera2->SetName("gIera2");
  gIera2->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (0.75 ^{3}S_{1} + 0.25 ^{1}S_{0})");
  gIera2->SetMarkerStyle(7);
  gIera2->SetMarkerColor(2);

  gIera2s = new TGraphErrors();
  gIera2s->SetName("gIera2s");
  gIera2s->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{1}S_{0})");
  gIera2s->SetMarkerStyle(2);
  gIera2s->SetMarkerColor(2);

  gIera2t = new TGraphErrors();
  gIera2t->SetName("gIera2t");
  gIera2t->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{3}S_{1})");
  gIera2t->SetMarkerStyle(2);
  gIera2t->SetMarkerColor(2);    

  



  gIera3 = new TGraphErrors();
  gIera3->SetName("gIera3");
  gIera3->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (0.75 ^{3}S_{1} + 0.25 ^{1}S_{0})");
  gIera3->SetMarkerStyle(7);
  gIera3->SetMarkerColor(3);

  gIera3s = new TGraphErrors();
  gIera3s->SetName("gIera3s");
  gIera3s->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{1}S_{0})");
  gIera3s->SetMarkerStyle(3);
  gIera3s->SetMarkerColor(3);

  gIera3t = new TGraphErrors();
  gIera3t->SetName("gIera3t");
  gIera3t->SetTitle("Verma Potentail w FSI Influence Factor (ERA); P_{#Lambda-n} [MeV] ; Influence Factor (^{3}S_{1})");
  gIera3t->SetMarkerStyle(3);
  gIera3t->SetMarkerColor(3);    


  for(int l=0;l<20;l++){
    gradl1[l] = new TGraphErrors();
    gradl1[l]->SetName(Form("gradl1_%d",l));
    gradl1[l]->SetTitle("Vema Potential w FSI Phase Shift; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
    gradl1[l]->SetMarkerStyle(2);
    gradl1[l]->SetMarkerColor(l+1);
    
    gradl2[l] = new TGraphErrors();
    gradl2[l]->SetName(Form("gradl2_%d",l));
    gradl2[l]->SetTitle("Vema Potential w FSI Phase Shift; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
    gradl2[l]->SetMarkerStyle(3);
    gradl2[l]->SetMarkerColor(l+1);

    gradl3[l] = new TGraphErrors();
    gradl3[l]->SetName(Form("gradl3_%d",l));
    gradl3[l]->SetTitle("Vema Potential w FSI Phase Shift; P_{#Lambda-n} [MeV] ; phase shift [rad] ");
    gradl3[l]->SetMarkerStyle(5);
    gradl3[l]->SetMarkerColor(l+1);


  gIl1[l] = new TGraphErrors();
  gIl1[l]->SetName(Form("gIl1_%d",l));
  gIl1[l]->SetTitle("Vema Potential w FSI Phase Shift; phase shift [rad]; Influence Factor ");
  gIl1[l]->SetMarkerStyle(2);
  gIl1[l]->SetMarkerColor(l+1);

  gIl2[l] = new TGraphErrors();
  gIl2[l]->SetName(Form("gIl2_%d",l));
  gIl2[l]->SetTitle("Vema Potential w FSI Phase Shift; phase shift [rad]; Influence Factor ");
  gIl2[l]->SetMarkerStyle(3);
  gIl2[l]->SetMarkerColor(l+1);
  
  gIl3[l] = new TGraphErrors();
  gIl3[l]->SetName(Form("gIl3_%d",l));
  gIl3[l]->SetTitle("Vema Potential w FSI Phase Shift; phase shift [rad]; Influence Factor ");
  gIl3[l]->SetMarkerStyle(5);
  gIl3[l]->SetMarkerColor(l+1);
  

  gIl01[l] = new TGraphErrors();
  gIl01[l]->SetName(Form("gIl01_%d",l));
  gIl01[l]->SetTitle("Vema Potential w FSI Phase Shift; phase shift [rad]; Influence Factor ");
  gIl01[l]->SetMarkerStyle(2);
  gIl01[l]->SetMarkerColor(l+1);

  gIl02[l] = new TGraphErrors();
  gIl02[l]->SetName(Form("gIl02_%d",l));
  gIl02[l]->SetTitle("Vema Potential w FSI Phase Shift; phase shift [rad]; Influence Factor ");
  gIl02[l]->SetMarkerStyle(3);
  gIl02[l]->SetMarkerColor(l+1);  

  gIl03[l] = new TGraphErrors();
  gIl03[l]->SetName(Form("gIl03_%d",l));
  gIl03[l]->SetTitle("Vema Potential w FSI Phase Shift; phase shift [rad]; Influence Factor ");
  gIl03[l]->SetMarkerStyle(5);
  gIl03[l]->SetMarkerColor(l+1);

  

  gfi1[l] = new TGraphErrors();
  gfi1[l]->SetName(Form("gfi1_%d",l));
  gfi1[l]->SetTitle("Vema Potential w FSI Scatering Amplitude; krel [MeV];Im[fthe] ");
  gfi1[l]->SetMarkerStyle(2);
  gfi1[l]->SetMarkerColor(l+1);

  gfi2[l] = new TGraphErrors();
  gfi2[l]->SetName(Form("gfi2_%d",l));
  gfi2[l]->SetTitle("Vema Potential w FSI Scatering Amplitude; krel [MeV];Im[fthe] ");
  gfi2[l]->SetMarkerStyle(3);
  gfi2[l]->SetMarkerColor(l+1);  

  gfi3[l] = new TGraphErrors();
  gfi3[l]->SetName(Form("gfi3_%d",l));
  gfi3[l]->SetTitle("Vema Potential w FSI Scatering Amplitude; krel [MeV];Im[fthe] ");
  gfi3[l]->SetMarkerStyle(2);
  gfi3[l]->SetMarkerColor(l+1);  



  gfr1[l] = new TGraphErrors();
  gfr1[l]->SetName(Form("gfr1_%d",l));
  gfr1[l]->SetTitle("Vema Potential w FSI Scatering Amplitude; krel [MeV];Im[fthe] ");
  gfr1[l]->SetMarkerStyle(2);
  gfr1[l]->SetMarkerColor(l+1);

  gfr2[l] = new TGraphErrors();
  gfr2[l]->SetName(Form("gfr2_%d",l));
  gfr2[l]->SetTitle("Vema Potential w FSI Scatering Amplitude; krel [MeV];Im[fthe] ");
  gfr2[l]->SetMarkerStyle(3);
  gfr2[l]->SetMarkerColor(l+1);  

  gfr3[l] = new TGraphErrors();
  gfr3[l]->SetName(Form("gfr3_%d",l));
  gfr3[l]->SetTitle("Vema Potential w FSI Scatering Amplitude; krel [MeV];Im[fthe] ");
  gfr3[l]->SetMarkerStyle(2);
  gfr3[l]->SetMarkerColor(l+1);  
  
  
  


  
  }






  
  
  double min_mm = -300.;
  double max_mm =  200.;
  int bin_mm = (int)(max_mm-min_mm)/2;  // 2 MeV bin
  //  bin_mm = (int)(max_mm-min_mm)/1;  // 1 MeV bin
  // bin_mm = (int)(max_mm-min_mm)/5; // 5 MeV bin
  int bin_size = (int)(max_mm - min_mm)/bin_mm;
  hmm = new TH1F("hmm",Form("Missing Mass w/o FSI ; -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);

  hmm ->SetLineColor(8);
  hmm ->SetFillColor(8);
  hmm ->SetFillStyle(3002);
  
  hmm_fsi1 = new TH1F("hmm_fsi1",Form("Missing Mass w/ FSI (Verma potenrial); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi1 ->SetLineColor(1);
  
  hmm_fsi2 = new TH1F("hmm_fsi2",Form("Missing Mass w/ FSI (Julich A potenrial); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi2 ->SetLineColor(2);
  hmm_fsi3 = new TH1F("hmm_fsi3",Form("Missing Mass w/ FSI (Julich B potenrial); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi3 ->SetLineColor(4);


  hmm_fsi1_2 = new TH1F("hmm_fsi1_2",Form("Missing Mass (Calc B) w/ FSI (Verma potenrial); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi1_2 ->SetLineColor(1);
  hmm_fsi1_2 ->SetLineStyle(10);
  
  hmm_fsi2_2 = new TH1F("hmm_fsi2_2",Form("Missing Mass (Calc B) w/ FSI (Julich A potenrial); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi2_2 ->SetLineColor(2);
  hmm_fsi2_2 ->SetLineStyle(10);  

  hmm_fsi3_2 = new TH1F("hmm_fsi3_2",Form("Missing Mass (Calc B) w/ FSI (Julich B potenrial); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi3_2 ->SetLineColor(4);
  hmm_fsi3_2 ->SetLineStyle(10);


  // Singlet

  hmm_fsi1s = new TH1F("hmm_fsi1s",Form("Missing Mass w/ FSI (Verma potenrial Singlet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi1s ->SetLineColor(1);
  
  hmm_fsi2s = new TH1F("hmm_fsi2s",Form("Missing Mass w/ FSI (Julich A potenrial Singlet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi2s ->SetLineColor(2);
  hmm_fsi3s = new TH1F("hmm_fsi3s",Form("Missing Mass w/ FSI (Julich B potenrial Singlet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi3s ->SetLineColor(4);


  hmm_fsi1_2s = new TH1F("hmm_fsi1_2s",Form("Missing Mass (Calc B) w/ FSI (Verma potenrial Singlet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi1_2s ->SetLineColor(1);
  hmm_fsi1_2s ->SetLineStyle(10);
  
  hmm_fsi2_2s = new TH1F("hmm_fsi2_2s",Form("Missing Mass (Calc B) w/ FSI (Julich A potenrial Singlet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi2_2s ->SetLineColor(2);
  hmm_fsi2_2s ->SetLineStyle(10);  

  hmm_fsi3_2s = new TH1F("hmm_fsi3_2s",Form("Missing Mass (Calc B) w/ FSI (Julich B potenrial Singlet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi3_2s ->SetLineColor(4);
  hmm_fsi3_2s ->SetLineStyle(10);    


  // Triplet

  hmm_fsi1t = new TH1F("hmm_fsi1t",Form("Missing Mass w/ FSI (Verma potenrial Triplet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi1t ->SetLineColor(6);
  
  hmm_fsi2t = new TH1F("hmm_fsi2t",Form("Missing Mass w/ FSI (Julich A potenrial Triplet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi2t ->SetLineColor(2);
  hmm_fsi3t = new TH1F("hmm_fsi3t",Form("Missing Mass w/ FSI (Julich B potenrial Triplet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi3t ->SetLineColor(4);

  hmm_fsi1_2t = new TH1F("hmm_fsi1_2t",Form("Missing Mass (Calc B) w/ FSI (Verma potenrial Triplet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi1_2t ->SetLineColor(6);
  hmm_fsi1_2t ->SetLineStyle(10);
  
  hmm_fsi2_2t = new TH1F("hmm_fsi2_2t",Form("Missing Mass (Calc B) w/ FSI (Julich A potenrial Triplet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi2_2t ->SetLineColor(2);
  hmm_fsi2_2t ->SetLineStyle(10);  

  hmm_fsi3_2t = new TH1F("hmm_fsi3_2t",Form("Missing Mass (Calc B) w/ FSI (Julich B potenrial Triplet); -B_{#Lambda} [MeV] ; Counts/%d MeV",bin_size),bin_mm,min_mm,max_mm);
  hmm_fsi3_2t ->SetLineColor(4);
  hmm_fsi3_2t ->SetLineStyle(10);    

  

  for(int i=0;i<vmax;i++){

    fJls[i] = new TF1(Form("fJls_%d",i),fJ,0,2000,2);
    fJls[i]->SetTitle(Form("%s_s ; P_{rel} [MeV] ;Influence Factor (^{1}S_{0})",vname[i].c_str()));
    fJls[i]->SetParameter(0,a_s[i]);
    fJls[i]->SetParameter(1,r_s[i]);
    fJls[i]->SetNpx(2000);
    fJlt[i] = new TF1(Form("fJlt_%d",i),fJ,0,2000,2);
    fJlt[i]->SetTitle(Form("%s_t ; P_{rel} [MeV] ;Influence Factor (^{3}S_{1})",vname[i].c_str()));
    fJlt[i]->SetParameter(0,a_t[i]);
    fJlt[i]->SetParameter(1,r_t[i]);    
    fJlt[i]->SetNpx(2000);
    fJl[i] = new TF1(Form("fJl_%d",i),fJ2,0,2000,4);
    fJl[i]->SetTitle(Form("%s ; P_{rel} [MeV] ;Influence Factor (0.25 ^{1}S_{0} + 0.75 ^{3}S_{1})",vname[i].c_str()));
    fJl[i]->SetParameter(0,a_s[i]);
    fJl[i]->SetParameter(1,r_s[i]);        
    fJl[i]->SetParameter(2,a_t[i]);
    fJl[i]->SetParameter(3,r_t[i]);        
    fJl[i]->SetNpx(2000);
    //    cout<<"test "<<i<<endl;
    
    hmm_Jl[i]= new TH1F(Form("hmm_Jl_%d",i),Form("w/ FSI (Jost Function ) %s Potential (0.25 ^{1}S_{0} + 0.75 ^{3}S_{1})",vname[i].c_str()),bin_mm,min_mm,max_mm);
    hmm_Jls[i]= new TH1F(Form("hmm_Jls_%d",i),Form("w/ FSI (Jost Function ) %s Potential (^{1}S_{0})",vname[i].c_str()),bin_mm,min_mm,max_mm);
    hmm_Jlt[i]= new TH1F(Form("hmm_Jlt_%d",i),Form("w/ FSI (Jost Function ) %s Potential (^{3}S_{1})",vname[i].c_str()),bin_mm,min_mm,max_mm);        

    hmm_Jl_2[i]= new TH1F(Form("hmm_Jl2_%d",i),Form("w/ FSI (Jost Function ) %s Potential (0.25 ^{1}S_{0} + 0.75 ^{3}S_{1})",vname[i].c_str()),bin_mm,min_mm,max_mm);
    hmm_Jl_2s[i]= new TH1F(Form("hmm_Jl_2s_%d",i),Form("w/ FSI (Jost Function ) %s Potential (^{1}S_{0})",vname[i].c_str()),bin_mm,min_mm,max_mm);
    hmm_Jl_2t[i]= new TH1F(Form("hmm_Jl_2t_%d",i),Form("w/ FSI (Jost Function ) %s Potential (^{3}S_{1})",vname[i].c_str()),bin_mm,min_mm,max_mm);


    hmm_Jl[i]->SetLineColor(i+1);
    hmm_Jls[i]->SetLineColor(i+1);
    hmm_Jlt[i]->SetLineColor(i+1);
    hmm_Jl_2[i]->SetLineColor(i+1);
    hmm_Jl_2s[i]->SetLineColor(i+1);
    hmm_Jl_2t[i]->SetLineColor(i+1);
    fJl[i]->SetLineColor(i+1);
    fJls[i]->SetLineColor(i+1);
    fJlt[i]->SetLineColor(i+1);
    if(i==9){
    hmm_Jl[i]->SetLineColor(i+2);
    hmm_Jls[i]->SetLineColor(i+2);
    hmm_Jlt[i]->SetLineColor(i+2);
    hmm_Jl_2[i]->SetLineColor(i+2);
    hmm_Jl_2s[i]->SetLineColor(i+2);
    hmm_Jl_2t[i]->SetLineColor(i+2);
    fJl[i]->SetLineColor(i+2);
    fJls[i]->SetLineColor(i+2);
    fJlt[i]->SetLineColor(i+2);

    }
  }


 
  hmm_L = new TH1F("hmm_L","H mass Missing Mass w/o FSI ; -B_{#Lambda} [MeV] ; Counts/1MeV",bin_mm,min_mm,max_mm);
  hmm_nnL = new TH1F("hmm_nnL","T mass Missing Mass w/o FSI ; -B_{#Lambda} [MeV] ; Counts/1MeV",bin_mm,min_mm,max_mm);    

}

///////////////////////////////////////////////////////////////////////


double Chi2(TH1D* hmm, TH1D* hmmFSI){

  int nbins = hmm->GetXaxis()->GetNbins();
  int nmin  = hmm->GetXaxis()->GetXmin();
  int nmax  = hmm->GetXaxis()->GetXmax();

  double ymax0 = hmm->GetYaxis()->GetXmax();
  double ymax1 = hmmFSI->GetYaxis()->GetXmax();

  double y0[nbins],y1[nbins];
  int bin_x0=0;
  double x0;
  for(int ibin = 0;ibin<nbins;ibin++){
    x0 = hmm -> GetXaxis()->GetBinCenter(ibin);
    if( x0 < 0.0 )bin_x0 = ibin;
    y0[ibin]  =  hmm    -> GetBinContent(ibin);
    y1[ibin]  =  hmmFSI -> GetBinContent(ibin);
  }


  double chi2,chi2_min,w,wmin;
  int wmax=100;
  double width = 0.1;
  chi2_min=100.;
  for(int wi=0; wi<wmax;wi++){
    chi2=0.0;
    w = ymax0/ymax1 + (double)(-wmax/2. + wi)*width; 
    for(int ibin = bin_x0;ibin<nbins;ibin++){
      chi2 +=  pow((y0[ibin] -w*y1[ibin])/y0[ibin],2.0);  
    }
    if(chi2_min > chi2){
      chi2_min = chi2;
      wmin = w;
    }
  } // for w

  return wmin;
  

}

//////////////////////////////////////////////////////////////////////

void fsi::SetBranch(string ifrname){



  if(SIMC){
    cout<<"SetBranch ; "<<ifrname<<endl;
    T = new TChain("SNT");  
    T->Add(ifrname.c_str());
    T->SetBranchAddress("Rp_gen",& vertex_p_P);
    T->SetBranchAddress("Rth_gen",& vertex_p_xptar);
    T->SetBranchAddress("Rph_gen",& vertex_p_yptar);
    T->SetBranchAddress("Lp_gen",& vertex_e_P);
    T->SetBranchAddress("Lth_gen",& vertex_e_xptar);
    T->SetBranchAddress("Lph_gen",& vertex_e_yptar);
    T->SetBranchAddress("uq_x_gen",& vertex_uq_x);
    T->SetBranchAddress("uq_y_gen",& vertex_uq_y);
    T->SetBranchAddress("uq_z_gen",& vertex_uq_z);
    T->SetBranchAddress("q_gen",& vertex_q);
    T->SetBranchAddress("nu_gen",& vertex_nu);
    T->SetBranchAddress("pfer",& pfer);  //MeV
    T->SetBranchAddress("efer",& efer);  //MeV
    T->SetBranchAddress("pferx",& pferx);  //MeV
    T->SetBranchAddress("pfery",& pfery);  //MeV
    T->SetBranchAddress("pferz",& pferz);  //MeV
    T->SetBranchAddress("Mrec",&  Mrec);  //MeV
    T->SetBranchAddress("mmnuc",& mmnuc);  //MeV
    T->SetBranchAddress("Em_v",& Em);  //MeV
    T->SetBranchAddress("mm_nnL",& mm_nnL);  //MeV
    T->SetBranchAddress("mm_L",& mm_L);  //MeV
    T->SetBranchAddress("Trec",& Trec);  //MeV
    
    
    
  }


  
}

///////////////////////////////////////////////////////////////////////


void fsi::SetMom(){

  cout<<"================================"<<endl;
  cout<<"====== SetMom  ================="<<endl;
  cout<<"================================"<<endl;
  

  TLorentzVector P2n,P2n_2, Pn,Pn_2, Pp, Pv,PK,Pe_,Pe;
  double P2nx,P2ny,P2nz,E2n,Ppx,Ppy,Ppz,Ep;
  double Pvx,Pvy,Pvz,Ev,Pkx,Pky,Pkz,Ek;
  double Ee, Ee_,Pe_x,Pe_y,Pe_z;

  ENum=  T->GetEntries();
  //  double Ee = 4.3185;
  //  ENum =10;
  cout<<" ENum : "<<ENum<<endl;
  for(int n=0;n<ENum;n++){


    Prel  = 0.0;
    Prel2 = 0.0;
    fac1  = 0.0;
    fac2  = 0.0;
    fac3  = 0.0;
    fac1t  = 0.0;
    fac2t  = 0.0;
    fac3t  = 0.0;
    fac1s  = 0.0;
    fac2s  = 0.0;
    fac3s  = 0.0;    
    fac1_2  = 0.0;
    fac2_2  = 0.0;
    fac3_2  = 0.0;
    fac1_2s  = 0.0;
    fac2_2s  = 0.0;
    fac3_2s  = 0.0;
    fac1_2t  = 0.0;
    fac2_2t  = 0.0;
    fac3_2t  = 0.0;
    
    T->GetEntry(n);

    // Set unit GeV 
    pfer       = pfer/1000.;
    efer       = efer/1000.;
    vertex_p_P = vertex_p_P/1000.;
    vertex_e_P = vertex_e_P/1000.;
    vertex_nu  = vertex_nu/1000.;
    Mrec       = Mrec/1000.;
    Em         = Em/1000.;
    // Set Proton Kinematics   //
    Ppx = pfer * pferx; // GeV
    Ppy = pfer * pfery; // GeV
    Ppz = pfer * pferz; // GeV
    //    Ep  = sqrt(pfer*pfer + Mp*Mp);
    Ep = efer;
    Pp.SetPxPyPzE(Ppx,Ppy,Ppz,Ep);

    // Set Virtual Photon Kinematics //
    
    Ev  = vertex_nu; //  Ee - Ee_ GeV
    Ee = Ev + vertex_e_P;
    Pe_z = vertex_e_P/sqrt(1.0 + vertex_e_xptar*vertex_e_xptar + vertex_e_yptar* vertex_e_yptar);
    Pe_x = Pe_z * vertex_e_xptar;
    Pe_y = Pe_z * vertex_e_yptar;
    Ee_ = sqrt(vertex_e_P*vertex_e_P + Me*Me);    
    Pe_.SetPxPyPzE(Pe_x,Pe_y,Pe_z,Ee_);
    Pe_.RotateX(-13.2*PI/180.);
    Pe.SetPxPyPzE(0.0,0.0,Ee,Ee);
    Pv = Pe - Pe_;


    // Set Kaon Kinematics //
    
    Pkz = vertex_p_P/sqrt(1.0 + vertex_p_xptar*vertex_p_xptar + vertex_p_yptar* vertex_p_yptar);
    Pkx = Pkz * vertex_p_xptar;
    Pky = Pkz * vertex_p_yptar;
    Ek = sqrt(vertex_p_P*vertex_p_P + MK*MK);
    PK.SetPxPyPzE(Pkx,Pky,Pkz,Ek);
    PK.RotateX(13.2*PI/180.);
    //    cout<<" Ev "<<Pv.E()<<" Pv "<<Pv.P()<<" Mag "<<Pv.Mag()<<endl;

    // Set Nuetron Kinematics // 


    
    P2nx = - Ppx;
    P2ny = - Ppy;
    P2nz = - Ppz;
    E2n  = sqrt(pfer*pfer + Mrec*Mrec);
    P2n.SetPxPyPzE(P2nx,P2ny,P2nz,E2n);
    P2n_2= P2n;

    Pn = CalcPn(P2n);
    pnx = Pn.Px();
    pny = Pn.Py();
    pnz = Pn.Pz();
    pn  = Pn.P();
    
    Pn_2 = CalcPn2(P2n_2);    
    pnx2 = Pn_2.Px();
    pny2 = Pn_2.Py();
    pnz2 = Pn_2.Pz();
    pn2  = Pn_2.P();

    if(Deu){
      //      double En_d = sqrt(Pp.Px()*Pp.Px() +Pp.Py()*Pp.Py()+Pp.Pz()*Pp.Pz() + Mn*Mn);
      double En_d = MD - Ep + Em;
      Pn.SetPxPyPzE(-Pp.Px(), -Pp.Py(), -Pp.Pz(), En_d);
      Pn_2.SetPxPyPzE(-Pp.Px(), -Pp.Py(), -Pp.Pz(), En_d);      
      //      Pn.SetPxPyPzE(-Pp.Px(), -Pp.Py(), -Pp.Pz(), Trec);
      //      Pn_2.SetPxPyPzE(-Pp.Px(), -Pp.Py(), -Pp.Pz(), Trec);
      //      cout<<"Deu mode "<<Deu<<" En_d "<<En_d<<" Trec "<<Trec<<endl;
    }
    //    Pn.SetPxPyPzE(Pnx,Pny,Pnz,En);
    //    Prel = CalcQ(Pn,Pp,APv,PK); // MeV


    
    TLorentzVector Pp_2,Pv_2,PK_2;
    Pp_2 =Pp;
    Pv_2 =Pv;
    PK_2 =PK;
    
    Prel = CalcQ(Pn,Pp,Pv,PK); // MeV
    Prel2 = CalcQ2(Pn_2,Pp_2,Pv_2,PK_2); // MeV  New version 


 
    
    for(int i=0;i<np1;i++){
      if(Prel < krel1[i]){
	if(i==0 || i==np1-1){
	  fac1s = w1[i];
	  fac1t = w1t[i];
	  fac1  = (fac1s + 3.0*fac1t)/4.0;	  
	}else {
	  double a = (w1[i+1] -w1[i])/(krel1[i+1] - krel1[i]);
	  double b = w1[i] - a*krel1[i];
	  double at = (w1t[i+1] -w1t[i])/(krel1[i+1] - krel1[i]);
	  double bt = w1t[i] - at*krel1[i];
	  
	  fac1s   =  a*Prel  + b;
	  fac1t   =  at*Prel  + bt;
	  fac1    =(fac1s + 3.0*fac1t)/4.0;

	}
	
	break;
      }else if(i==np1-1){
	fac1s = w1[i];
	fac1t = w1t[i];
	fac1  = (fac1s + 3.0*fac1t)/4.0;	  
      }
    } // end fac1
    
    
    for(int i=0;i<np2;i++){
      if(Prel < krel2[i]){
	if(i==0 || i==np2-1){
 	  fac2s = w2[i];
	  fac2t = w2t[i];
	  fac2  = (fac2s + 3.0*fac2t)/4.0;	  	  
	}else {
	  double a = (w2[i+1] -w2[i])/(krel2[i+1] - krel2[i]);
	  double b = w2[i] - a*krel2[i];

	  double at = (w2t[i+1] -w2t[i])/(krel2[i+1] - krel2[i]);
	  double bt = w2t[i] - at*krel2[i];
	  
	  fac2s =  a*Prel   + b;
	  fac2t =  at*Prel  + bt;
	  fac2  = (fac2s + 3.0*fac2t)/4.0;	  
	}
	break;
      }else if(i==np2-1){
	fac2s = w2[i];
	fac2t = w2t[i];
	fac2  = (fac2s + 3.0*fac2t)/4.0;	  	  
      }
    } // end fac2
    
     
    for(int i=0;i<np3;i++){
      if(Prel < krel3[i]){
	if(i==0 || i==np3-1){
	  fac3s   = w3[i];
	  fac3t  = w3t[i];
	  fac3   = (fac3s + 3.0*fac3t)/4.0;	  	  	  
	}else {
	  double a = (w3[i+1] -w3[i])/(krel3[i+1] - krel3[i]);
	  double b = w3[i] - a*krel3[i];
	  double at = (w3t[i+1] -w3t[i])/(krel3[i+1] - krel3[i]);
	  double bt = w3t[i] - at*krel3[i];
	  
	  fac3s =  a*Prel   + b;
	  fac3t =  at*Prel  + bt;
	  fac3  = (fac3s + 3.0*fac3t)/4.0;	  	  
	}
	
	break;
      }else if(i==np3-1){
	fac3   = w3[i];
	fac3t  = w3t[i];
	fac3   = (fac3s + 3.0*fac3t)/4.0;	  	  	  
      }
    } // end fac3
    



    //  Prel2  //


    
    for(int i=0;i<np1;i++){
      if(Prel2 < krel1[i]){
	if(i==0 || i==np1-1){
	  fac1_2s = w1[i];
	  fac1_2t = w1t[i];
	  fac1_2  = (fac1_2s + 3.0*fac1_2t)/4.0;	  
	}else {
	  double a = (w1[i+1] -w1[i])/(krel1[i+1] - krel1[i]);
	  double b = w1[i] - a*krel1[i];
	  double at = (w1t[i+1] -w1t[i])/(krel1[i+1] - krel1[i]);
	  double bt = w1t[i] - at*krel1[i];
	  
	  fac1_2s   =  a*Prel2  + b;
	  fac1_2t   =  at*Prel2  + bt;
	  fac1_2    =(fac1_2s + 3.0*fac1_2t)/4.0;

	}
	
	break;
      }else if(i==np1-1){
	fac1_2s = w1[i];
	fac1_2t = w1t[i];
	fac1_2  = (fac1_2s + 3.0*fac1_2t)/4.0;	  
      }
    } // end fac1_2
    


    for(int i=0;i<np2;i++){
      if(Prel2 < krel2[i]){
	if(i==0 || i==np2-1){
 	  fac2_2s = w2[i];
	  fac2_2t = w2t[i];
	  fac2_2  = (fac2_2s + 3.0*fac2_2t)/4.0;	  	  
	}else {
	  double a = (w2[i+1] -w2[i])/(krel2[i+1] - krel2[i]);
	  double b = w2[i] - a*krel2[i];

	  double at = (w2t[i+1] -w2t[i])/(krel2[i+1] - krel2[i]);
	  double bt = w2t[i] - at*krel2[i];
	  
	  fac2_2s =  a*Prel2   + b;
	  fac2_2t =  at*Prel2  + bt;
	  fac2_2  = (fac2_2s + 3.0*fac2_2t)/4.0;	  
	}
	break;
      }else if(i==np2-1){
	fac2_2s = w2[i];
	fac2_2t = w2t[i];
	fac2_2  = (fac2_2s + 3.0*fac2_2t)/4.0;	  	  
      }
    } // end fac2_2
    
     
    for(int i=0;i<np3;i++){
      if(Prel2 < krel3[i]){
	if(i==0 || i==np3-1){
	  fac3_2s   = w3[i];
	  fac3_2t  = w3t[i];
	  fac3_2   = (fac3_2s + 3.0*fac3_2t)/4.0;	  	  	  
	}else {
	  double a = (w3[i+1] -w3[i])/(krel3[i+1] - krel3[i]);
	  double b = w3[i] - a*krel3[i];
	  double at = (w3t[i+1] -w3t[i])/(krel3[i+1] - krel3[i]);
	  double bt = w3t[i] - at*krel3[i];
	  
	  fac3_2s =  a*Prel2   + b;
	  fac3_2t =  at*Prel2  + bt;
	  fac3_2  = (fac3_2s + 3.0*fac3_2t)/4.0;	  	  
	}
	
	break;
      }else if(i==np3-1){
	fac3_2s  = w3[i];
	fac3_2t  = w3t[i];
	fac3_2   = (fac3_2s + 3.0*fac3_2t)/4.0;	  	  	  
      }
    } // end fac3_2
    


    

    

    if(!E09){
    fac1 = (fac1-1.0)*2.0 +1.0;
    fac2 = (fac2-1.0)*2.0 +1.0;
    fac3 = (fac3-1.0)*2.0 +1.0;
    fac1_2 = (fac1_2-1.0)*2.0 +1.0;
    fac2_2 = (fac2_2-1.0)*2.0 +1.0;
    fac3_2 = (fac3_2-1.0)*2.0 +1.0;

    fac1s = (fac1s-1.0)*2.0 +1.0;
    fac2s = (fac2s-1.0)*2.0 +1.0;
    fac3s = (fac3s-1.0)*2.0 +1.0;
    fac1_2s = (fac1_2s-1.0)*2.0 +1.0;
    fac2_2s = (fac2_2s-1.0)*2.0 +1.0;
    fac3_2s = (fac3_2s-1.0)*2.0 +1.0;

    fac1t = (fac1t-1.0)*2.0 +1.0;
    fac2t = (fac2t-1.0)*2.0 +1.0;
    fac3t = (fac3t-1.0)*2.0 +1.0;
    fac1_2t = (fac1_2t-1.0)*2.0 +1.0;
    fac2_2t = (fac2_2t-1.0)*2.0 +1.0;
    fac3_2t = (fac3_2t-1.0)*2.0 +1.0;        
    
    }
    


    /*    
    //test //
    if(Prel>400){
    fac1 = 1.0;
    fac2 = 1.0;
    fac3 = 1.0;
    }
    if(Prel2>400){
    fac1_2 = 1.0;
    fac2_2 = 1.0;
    fac3_2 = 1.0;
    }
    */

    

              mm =(mmnuc - MnnL )*1000.; // -BL [MeV]
    if(E09)   mm =(mmnuc - MH3L )*1000.; // -BL [MeV]
    if(Deu)   mm =(mmnuc - (ML+Mn) )*1000.; // -BL [MeV]


    hmm -> Fill(mm);

    // Single //
    hmm_fsi1s ->Fill(mm , fac1s);
    hmm_fsi2s ->Fill(mm , fac2s);
    hmm_fsi3s ->Fill(mm , fac3s);
    hmm_fsi1_2s ->Fill(mm , fac1_2s);
    hmm_fsi2_2s ->Fill(mm , fac2_2s);
    hmm_fsi3_2s ->Fill(mm , fac3_2s);


    for(int i=0;i<vmax;i++){
      double s = Jl0(Prel,a_s[i],r_s[i]);
      double t = Jl0(Prel,a_t[i],r_t[i]);
      double s2 = Jl0(Prel2,a_s[i],r_s[i]);
      double t2 = Jl0(Prel2,a_t[i],r_t[i]);

      // Cut off Parameter
      double cut_off = 100000.;
      if(Prel  > cut_off){  s= 1.0; t =1.0;}
      if(Prel2 > cut_off){  s2= 1.0; t2 =1.0;}
    if(!E09){
    s = (s-1.0)*2.0 +1.0;
    t = (t-1.0)*2.0 +1.0;
    s2 = (s2-1.0)*2.0 +1.0;
    t2 = (t2-1.0)*2.0 +1.0;    
    }
    
      hmm_Jls[i] ->Fill(mm,s );
      hmm_Jlt[i] ->Fill(mm,t );
      hmm_Jl[i] ->Fill(mm,(s+3.0*t)/4.0 );

      hmm_Jl_2s[i] ->Fill(mm,s2 );
      hmm_Jl_2t[i] ->Fill(mm,t2 );
      hmm_Jl_2[i] ->Fill(mm,(s2+3.0*t2)/4.0 );
      
    }
    

    // Triplet //
    hmm_fsi1t ->Fill(mm , fac1t);
    hmm_fsi2t ->Fill(mm , fac2t);
    hmm_fsi3t ->Fill(mm , fac3t);
    hmm_fsi1_2t ->Fill(mm , fac1_2t);
    hmm_fsi2_2t ->Fill(mm , fac2_2t);
    hmm_fsi3_2t ->Fill(mm , fac3_2t);        
    // Single +3. Triple /4.0 
    hmm_fsi1 ->Fill(mm , fac1);
    hmm_fsi2 ->Fill(mm , fac2);
    hmm_fsi3 ->Fill(mm , fac3);
    hmm_fsi1_2 ->Fill(mm , fac1_2);
    hmm_fsi2_2 ->Fill(mm , fac2_2);
    hmm_fsi3_2 ->Fill(mm , fac3_2);    

    hmm_L->Fill(mm_L);
    hmm_nnL->Fill(mm_nnL);
    
    Tnew->Fill();    
    
    
  }
  


  for(int i=0;i<vmax;i++){
    hmm_Jl[i]->Write();
    hmm_Jl_2[i]->Write();
    hmm_Jls[i]->Write();
    hmm_Jl_2s[i]->Write();
    hmm_Jlt[i]->Write();
    hmm_Jl_2t[i]->Write();    
  }
  
  
  hmm->Write();
  hmm_L->Write();
  hmm_nnL->Write();
  hmm_fsi1->Write();
  hmm_fsi2->Write();
  hmm_fsi3->Write();
  hmm_fsi1_2->Write();
  hmm_fsi2_2->Write();
  hmm_fsi3_2->Write();  
  hmm_fsi1s->Write();
  hmm_fsi2s->Write();
  hmm_fsi3s->Write();
  hmm_fsi1_2s->Write();
  hmm_fsi2_2s->Write();
  hmm_fsi3_2s->Write();    
  hmm_fsi1t->Write();
  hmm_fsi2t->Write();
  hmm_fsi3t->Write();
  hmm_fsi1_2t->Write();
  hmm_fsi2_2t->Write();
  hmm_fsi3_2t->Write();  

  
}

///////////////////////////////////////////////////////////////////////


TLorentzVector fsi::CalcPn(TLorentzVector P2n){


  // Neutron 4-vector Calculation code //
  // P2n 2 neutron (recoil ) 4-momentum vector [GeV]
  //
  // Assumption : pn = p2n/2.0

  TLorentzVector Pn1,Pn2;
  // n1,n2 are number of nuetrons
  // C.M coordinate 
  TVector3 B; // Lorentz Boost
  B = P2n.BoostVector();
  P2n.Boost(-B);



  
  // cout<<"E2n .cm "<<P2n.E()<<endl;
  // Generate rondom n1 3-vector in C.M.
  double E2n = P2n.E();
  double En1,En2, pn1,pn2;
  En1 = E2n/2.0;
  En2 = E2n/2.0;
  if(En1> Mn ) pn1 = sqrt(En1*En1 - Mn*Mn);
  else       pn1 = 0.0;


  double px1,py1,pz1,px2,py2,pz2;
  ranth = acos( random.Uniform()*2.0 -1.0);
  ranph = random.Uniform()* 2.0 *PI;

  px1 = pn1 * sin(ranth) * cos(ranph);
  py1 = pn1 * sin(ranth) * sin(ranph);
  pz1 = pn1 * cos(ranth);
  Pn1.SetPxPyPzE(px1, py1, pz1 , En1);


  
  Pn2 = P2n - Pn1;
  
  // C.M to Lab System
  P2n.Boost(B);
  Pn1.Boost(B);
  Pn2.Boost(B);


  
  return Pn2;
  
  
  
}

///////////////////////////////////////////////////////////////////////



TLorentzVector fsi::CalcPn_test(TLorentzVector P2n){


  // Neutron 4-vector Calculation code //
  // P2n 2 neutron (recoil ) 4-momentum vector [GeV]
  //
  // Assumption : pn = p2n/2.0

  TLorentzVector Pn1,Pn2;
  // n1,n2 are number of nuetrons
  // C.M coordinate 
  TVector3 B; // Lorentz Boost
  B = P2n.BoostVector();
  P2n.Boost(-B);



  
  // cout<<"E2n .cm "<<P2n.E()<<endl;
  // Generate rondom n1 3-vector in C.M.
  double E2n = P2n.E();
  double En1,En2, pn1,pn2;
  En1 = E2n/2.0;
  En2 = E2n/2.0;
  if(En1> Mn ) pn1 = sqrt(En1*En1 - Mn*Mn);
  else       pn1 = 0.0;


  double px1,py1,pz1,px2,py2,pz2;
  ranth = acos( random.Uniform()*2.0 -1.0);
  ranph = random.Uniform()* 2.0 *PI;

  px1 = pn1 * sin(ranth) * cos(ranph);
  py1 = pn1 * sin(ranth) * sin(ranph);
  pz1 = pn1 * cos(ranth);
  Pn1.SetPxPyPzE(px1, py1, pz1 , En1);


  
  Pn2 = P2n - Pn1;
  
  // C.M to Lab System
  P2n.Boost(B);
  Pn1.Boost(B);
  Pn2.Boost(B);


  
  return Pn2;
  
  
  
}

///////////////////////////////////////////////////////////////////////


TLorentzVector fsi::CalcPn2(TLorentzVector P2n){


  // Neutron 4-vector Calculation code //
  // P2n 2 neutron (recoil ) 4-momentum vector [GeV]
  //
  // Assumption : pn = p2n/2.0

  TLorentzVector Pn1,Pn2;
  // n1,n2 are number of nuetrons
  // C.M coordinate 
  TVector3 B; // Lorentz Boost
  B = P2n.BoostVector();
  P2n.Boost(-B);
  // Generate rondom n1 3-vector in C.M.
  double E2n = P2n.E();
  double En1,En2, pn1,pn2, Mn1,Mn2;
  if(E09) Mn1 =Mp;
  else    Mn1 =Mn;  // nuetron mass
  Mn2 =Mn;  // nuetron mass
  // Generate Fermi Momentum in removal system
  double pferhi,pferlo,pfer2;
  double ranprob = random.Uniform();
  int ii =0;

  while(ranprob >= mprob[ii] && ii <= nump)ii++;
  if(ii == nump-1){
    pferlo = (pval[ii-1] + pval[ii] )/2.;
    pferhi = pval[nump-1]; 
  }else if(ii==0){
    pferlo = 0.0;
    pferhi = (pval[ii] + pval[ii+1] )/2.;
  }else {
    pferlo = (pval[ii-1] + pval[ii] )/2.;
    pferhi = (pval[ii] + pval[ii+1] )/2.;
  }
  

  pn1 = pferlo + (pferhi - pferlo) * random.Uniform();
  pn2 = pn1;
  double Emiss  = Mn1 + Mn2 - P2n.Mag();
  double m2  = (Pn2.Mag())  - Mn1 + Emiss;
  
  En2 = sqrt(m2*m2 + pn2*pn2);
  En1 = P2n.Mag() - En2;

  //=============================================
  

  double px1,py1,pz1,px2,py2,pz2;
  ranth = acos( random.Uniform()*2.0 -1.0);
  ranph = random.Uniform()* 2.0 *PI;

  px1 = pn1 * sin(ranth) * cos(ranph);
  py1 = pn1 * sin(ranth) * sin(ranph);
  pz1 = pn1 * cos(ranth);
  Pn1.SetPxPyPzE(px1, py1, pz1 , En1);
  
  Pn2 = P2n - Pn1;
  
  // C.M to Lab System
  P2n.Boost(B);
  Pn1.Boost(B);
  Pn2.Boost(B);
  
  return Pn2;
  
  
  
}



///////////////////////////////////////////////////////////////////////



double fsi::CalcQ(TLorentzVector Pn, TLorentzVector Pp, TLorentzVector Pv, TLorentzVector PK){

  // Lab system to C.M system TVector(E,Px,Py,Pz)
  // Pp : 4-Momentum of proton  [GeV] 
  // Pn : 4-Momentum of nuetron [GeV/]
  // pv : 4-Momentum of vertual photon [GeV]
  // pK : 4-Momentum of Kaon     [GeV]
  // Out put



  // relative momentum (Lambda-n) Center of mass system //
  // Loretnz Boost in C.M. system // 
  TLorentzVector Pcm;
  Pcm = Pv + Pp;
  TVector3 B; // Loretnz boost
  B = Pcm.BoostVector();
  TLorentzVector Ppi;
  Ppi = Pp;
  Ppi.Boost(-B);
  // final nuetron 4-momentum
  TLorentzVector Pnf;
  Pnf = Pn;
  double Ecm = Pcm.E();
  
  // final kaon 4-momentum in C.M. system
  TLorentzVector PKf;
  PKf = PK;
  PKf.Boost(-B);

  // final Lambda 4- momentum in C.M. system
  TLorentzVector PLf;
  PLf.SetPxPyPzE(-PKf.Px(), -PKf.Py(), -PKf.Pz(),sqrt(PKf.P()*PKf.P() +ML*ML  ) );

  
  PLf.Boost(B); // C.M. to Lab system
  
  TLorentzVector Pvi;
  Pvi = Pv;
  
  TLorentzVector Q;
  double q;
  
  Q = Pnf + PLf;
  TVector3 B2;
  B2 =Q.BoostVector();

  Pnf.Boost(-B2);
  PLf.Boost(-B2);
  
  q =  Pnf.P() + PLf.P() ;
  
  
  // Retrun GeV to MeV
  return q*1000.;


}

////////////////////////////////////////////////////////////////

double fsi::CalcQ2(TLorentzVector Pn, TLorentzVector Pp, TLorentzVector Pv, TLorentzVector PK){

  // Lab system to C.M system TVector(E,Px,Py,Pz)
  // Pp : 4-Momentum of proton  [GeV] 
  // Pn : 4-Momentum of nuetron [GeV/]
  // pv : 4-Momentum of vertual photon [GeV]
  // pK : 4-Momentum of Kaon     [GeV]
  // Out put of relative momentum (Lambda-n) Center of mass system //



  TLorentzVector PL;
  PL = Pp + Pv - PK;


  pL = PL.P();
  pLx = PL.Px();
  pLy = PL.Py();
  pLz = PL.Pz();

  TVector3 PLv;
  PLv.SetXYZ(pLx/pL,pLy/pL,pLz/pL);
  TVector3 Rv;
  Rv.SetXYZ(0.0,0.0,1.0);
  //  Rv.RotateX(13.2*PI/180.);

  theta_L = acos(Rv * PLv );
  theta_L *= 180./PI;
  
  
  // C.M. coordinate
  
  TLorentzVector Pcm;
  Pcm = PL + Pn;
  //  pn = Pn.P();

  
  TVector3 B;
  B = Pcm.BoostVector();

  PL.Boost(-B);
  Pn.Boost(-B);

  TLorentzVector Pr;
  Pr = PL -Pn;
  double q = Pr.P();


  
  // Retrun GeV to MeV
  return q*1000.;


}



/////////////////////////////////////////////////////////////////////////
double fsi::PhaseShift(double qq, int LL, int potential){

  
  
  //======= Energy Order is MeV ========//
  // fmu : dueteron mode
  amu    = Mu*1000.;
  fnucl  = Mp*1000.;
  fnucl2 = fnucl*fnucl;
  hbarc2 = hbarc*hbarc;
  ampj   =  (ML*1000.);
  ampj2  = ampj*ampj;
  amtag  = Mn*1000.;
  amtag2 = amtag*amtag;
  fmu = ampj*amtag/(ampj+amtag);
  fmunn = ampj*fnucl/(ampj+fnucl);
  xi=complex<double>(0.0,1.0);
 
  model = potential;
  skz = qq;
  nrp1 = n1+n2;
  ftheL = complex<double>(0.0,0.0);

  if(Cha_mode && model==2 && skz >100)model=1;
    
    double a0,r0;
         if(model==0){  a0 = -2.68; r0 = 2.91;}
    else if(model==1 ){ a0 = -2.29; r0 = 3.15;}
    else if(model==2 ){ a0 = -1.60; r0 = 1.33;}
    else if(model==3 ){ a0 = -0.57; r0 = 7.65;}
    else if(model==-1){ a0 = -1.77; r0 = 3.25;}
    else if(model==-2){ a0 = -1.6;  r0 = 3.15;}
    else if(model==-3){ a0 = -1.94; r0 = 2.42;}
    
	 if(Vscale){
	   a0 = a_0;
	   r0 = r_0;
	 }

	 double R0 = 3.15; //test
	 //test
	 if(Cha_mode)r0 =3.15;


	 
  // Kinematics
    if(nr_mode){
      skz = skz/hbarc; // [1/fm]
      skz2 = skz*skz;  
      eon = hbarc2*skz2/2.0/fmu; // Kinematics Energy [MeV]
      rho = fmu*skz/hbarc2; // Desity of states/energy [1/(MeV*fm^3)]

    }else {
      skz2 = skz*skz;
      ekpj  = sqrt(ampj2   + skz2);
      ektag = sqrt( amtag2 + skz2);
      eon = ekpj + ektag;
      rho = skz*ekpj*ektag/eon;
    }

    // Get Gauss pts and wtss for integral v(q)Pl(x)dx (x=cos(theta))
    DGauss(xvq,wvq,npot);
    DGauss(x,w,n1);
    for( int i =0; i<n1;i++){ // DO 20
      skk[i] = skz  * (x[i]+1.0);
      wt[i]  = w[i] * skz;

    } // end 20
    DGauss(x,w,n2);
    for( int i =0; i<n2;i++){ // DO 21
      int  ii = i +n1;
      double coss = cos(PI*(x[i]+1.0)/4.0);
      skk[ii] = 2.0 * skz +2.0 *skz * tan(PI*(x[i]+1.0 )/4.0);
      wt[ii]  = 2.0* PI * 0.25 *skz/coss/coss*w[i];
    } // end 21

    // Set Last Segment
    skk[nrp1] = skz;

    // start the l-loop Do 30
    // LL = the physical L . (0,1,2, ....)

    lpjmax =LL;
    
    for(int ll=0;ll<=lpjmax;ll++){ // DO 30
      //      int l = ll+1;
      int l = ll;
      // start K',K'' loop to 40 & 50

      for(int kp=0;kp <= nrp1;kp++){ // DO 40
	skp = skk[kp];
	for(int kpp=0;kpp <= nrp1;kpp++){ // DO 50
	  skpp  = skk[kpp];
	  skpp2 = skpp*skpp;
	  // Initialization
	  gren[kp][kpp] = 0.0;
	  ul[kp][kpp]   = 0.0;
	  // ============== //
	  v2 = vlkkp(skp,skpp,ll);
	  ul[kp][kpp]  = v2;
	  eoff = hbarc2*skpp2/2.0/fmu;
	  if(!nr_mode)eoff=sqrt(ampj2+skpp2)+ sqrt(amtag2+skpp2);
	  if(kpp != nrp1)gren[kp][kpp] = wt[kpp]*skpp2*ul[kp][kpp]/(eon - eoff); 
	  if(kp==kpp)  delta = complex<double>(1.0,0.0);
	  else 	       delta = complex<double>(0.0,0.0);

	  gren[kp][kpp] = delta - gren[kp][kpp];

	  // gren : I - V * Go  
	} // END 50
      } // END 40


      // Set the constant vector i.e. the B column matrix
      // in thet matrix equation AX = B

      
      for(int i=0;i<nmax2;i++){
	v[i]  = ul[i][nrp1];  // DO 60
	v0[i] = ul[i][nrp1]; 
      }
      

      int nmax3=(nrp1+1)*(nrp1+1);
      complex<double>A[nmax3],A0[nmax3];
      for(int i=0;i<=nrp1;i++){
	for(int j=0;j<=nrp1;j++){

	  //	  if(i!=j) A0[i*(nrp1+1)+j] = complex<double>(0.0,0.0);
	  //	  else A0[i*(nrp1+1)+j]=gren[i][j];
	  A[i*(nrp1+1)+j] = gren[i][j];

	  if(i==j) A0[i*(nrp1+1)+j] = gren[i][j];
		     // complex<double>(1.0,0.0);
	  else A0[i*(nrp1+1)+j] = complex<double>(0.0,0.0);
	  
	}
      }
      
      // AX =B matrix : R(E-VGo)=V   (C.11) 
      // A : gren[i][j] = E -V * Go (C.11)
      // B : V
      // X : R =(E - V * Go)^-1 * V (A) = A^-1 * B -> retrun X :  R | k> half-off shell 

      

      ron  = complex<double>(0.0,0.0);
      ron0 = complex<double>(0.0,0.0);
      LEQ3(A,v,Rl,nrp1);
      LEQ3(A0,v0,Rl0,nrp1);



      //==============================================//
      //=======< T matrix calculation >===============//
      //==============================================//

      
      // calc T matrix : T = R - i*pi*R*delta(E-Ho)T (C.12)
      // delta(E-Ho) = rho*delta(k-k')
      // <k|T = 1/(1+ i*rho*<k|R|k>)<k|R    (5.60)
      // T = 1.0
      complex<double>T[nmax3];
      toff[l] = complex<double>(0.0,0.0);

      
      for(int i=0;i<=nrp1;i++){

	T[i] = complex<double>(0.0,0.0);
	T[i] = 1.0/(1.0 + xi* rho * Rl[nrp1] * PI )* Rl[i];
	toff[l] += T[i];
      }


      toff[l] = toff[l]/double(nrp1+1.0);


      // Calc Influence Factor
      // I =  <k| (1 + TGo)
      






      
      
      ron       = Rl[nrp1];
      ron0      = Rl0[nrp1];
      ronr      = real(ron);
      delrad0   = atan(skz/(-1./a0 + r0*skz*skz/2.0));
      delrad    = -atan(PI*rho*ronr);
      deltal[l] = delrad;

      //      ton[l]    = ron/(1.0 +xi*PI*rho*ron);
      
      ton[l]    = T[nrp1];
      ton0[l]   = ron0/(1.0 +xi*PI*rho*ron0);
      jon[l]    = ton[l]/(1.0 - xi*PI*rho*ton[l]);
      vv        = vlkkp(skz,skz,ll);
      tborn[l]  = vv;
      tonre     = real(ton[l]);
      tonim     = imag(ton[l]);
      tbornre   = real(tborn[l]);
      tbornim   = imag(tborn[l]);
      perre     = abs(tbornre -tonre)/abs(tonre);
      perim     = abs(tbornim -tonim)/abs(tonim);
      deldeg    = delrad*180./PI;
      deldega   = abs(deldeg);
      double xx = cos(0.0);
      PL(xx, l, pol);
      
      
      fthel[l] = -(2.0*(double)l +1.0)*ton[l]*pol[l]/4.0/PI*4.0*PI*PI*rho/skz;
      ftheL   += fthel[l];
      //      fthel[l] = ftheL;
      expL[l]  = complex<double>(cos(deltal[l]), sin(deltal[l]));

      
    } // END 30


    double totalc =0.0;
    fborn =  complex<double>(0.0,0.0);
    
    for(int ithe=0;ithe<180;ithe++){ // DO 80
      double the = double(ithe);
      double ther = PI* the/180.;
      double q  = 2.0 * skz * sin(ther/2.0);
      double q2 = q*q;
      xxx = cos(ther);
      vq  = vofq(q2);
      PL(xxx, lpjmax, pol);
      fthe[ithe]   = complex<double>(0.0,0.0);
      fthe0[ithe]  = complex<double>(0.0,0.0);
      fthe1[ithe]  = complex<double>(0.0,0.0);      
      fthed        = complex<double>(0.0,0.0);
      Psi[ithe]    = complex<double>(0.0,0.0);
      Psi0[ithe]   = complex<double>(0.0,0.0);
      psi[ithe]    = complex<double>(0.0,0.0);
      psi0[ithe]   = complex<double>(0.0,0.0);
      
      for(int l =0;l <= lpjmax;l++){ // DO 90
	double ll=(double)l;
	fthe[ithe]   =  fthe[ithe]   + (2.0*ll +1.0)*ton[l]*pol[l]/4.0/PI;
	fthe0[ithe]  =  fthe0[ithe]  + (2.0*ll +1.0)*ton0[l]*pol[l]/4.0/PI;
	fthe1[ithe]  =  fthe1[ithe]  + (2.0*ll +1.0)*(ton[l] -tborn[l])*pol[l]/4.0/PI;
	

	
	double xxx0 = skz*r0;
	double xxx  = skz*r0+deltal[l];

	psi0[ithe] =  psi0[ithe]   + (2.0*ll +1.0)*pow(xi,l)*sin(xxx0-1./2.*ll*PI)/(xxx0)*pol[l];
	psi[ithe]  =  psi[ithe]    + (2.0*ll +1.0)*pow(xi,l)*expL[l]*sin(xxx -1./2.*ll*PI )/(xxx0)*pol[l];

	
	Psi0[ithe] =  Psi0[ithe]   + (2.0*ll +1.0)*pow(xi,l)*jn(l,xxx0)*pol[l];
	Psi[ithe]  =  Psi[ithe]    + (2.0*ll +1.0)*pow(xi,l)*expL[l]*(xxx)*jn(l,xxx)/(xxx0)*pol[l];

      } // END 90


      fthe[ithe]   = - fthe[ithe]*4.0*PI*PI*rho/skz;
      fthe0[ithe]  = - fthe0[ithe]*4.0*PI*PI*rho/skz;
      fthe1[ithe]  = -(fthe1[ithe]+vq)*4.0*PI*PI*rho/skz;
      fborn        += -vq*4.0*PI*PI*rho/skz;

      totalc = totalc + fabs(fthe[ithe])*fabs(fthe[ithe])*sin(ther)*PI/180.*2*PI;
      cross[ithe] = real((fthe[ithe])* conj(fthe[ithe]) )+1.0;


      
    } // END 80


    tcross   = 4.0 * PI/skz *  imag(fthe[0]);
    tcross0  = 4.0 * PI/skz *  imag(fthe0[0]);
    tcross1  = 4.0 * PI/skz *  imag(fthe1[0]);
    tcrossl  = 4.0 * PI/skz *  imag(ftheL);

  


    // Jost Fanction
    
    double alpha =  ( 1. + sqrt(1. - 2.*r0/a0)) / r0;
    double beta  =  (-1. + sqrt(1. - 2.*r0/a0)) / r0;
    
    complex<double> Jl = (skz - xi*beta)/(skz - xi*alpha);
    
    complex<double> exp =complex<double>(cos(skz*r0),sin(skz*r0));
    complex<double> phi   =  exp + ftheL/r0*exp;
    complex<double> phi_  =  exp - ftheL/r0*exp; 
    complex<double> phi0  =  exp;


    I0  = 1./( real(Jl*conj(Jl)) ); // Jost Fanction
    if(Cha_mode)I1 = pow( sin(skz*hbarc*r0 + deltal[LL])/sin(skz*hbarc*r0) ,2.0); // ERA sin function (Cha's D thisis)
    else I1 = pow( sin(skz*r0 + deltal[LL])/sin(skz*r0) ,2.0); // ERA sin function
    





    
    
    //    IERA = pow( sin(skz*r0 + deltal[LL])/sin(skz*r0) ,2.0); // ERA
    IERA = pow( sin(skz*R0 + deltal[LL])/sin(skz*R0) ,2.0); // ERA test
    
    // T-operator matrix : I = <k'|T
    // if( k'==k) Iton (on-shell)
    // Iton  : <k|T|k> on-shell
    // Itoff : Sum(k'=0 ... nrp1+1 ) <k|T|k'>/(nrp1+1)  half-off-shell
    
    Iton  = sqrt( real( ton[LL]  * conj(ton[LL] )   )  );
    Itoff = sqrt( real( toff[LL] * conj(toff[LL])   )  );
    

    
    double I =real( phi/phi0 *conj(phi/phi0) );
    
    for(int l=0;l<=LL;l++){
      Il[l]  = real(fthel[l]* conj(fthel[l])) + 1.0;
      Il0[l] = imag(fthel[l])/skz+ 1.0;
    }


    //    cout<<"qq "<<qq<<" skz "<<skz<<" skz*hbarc "<<skz*hbarc<<endl;
    
    if( I0 >100.  ) I0 = 0.0;
    if( I1 >100.  ) I1 = 0.0;
    if(Cha_mode) I = I1;
    
    if(tcross0>0 || I <100.) return I;
    else     return 1.0;


}

///////////////////////////////////////////////////////////////////////


double fsi::CalcIfac(int LL, int potential, double Prel){

  
  //======= Energy Order is MeV ========//
  // fmu : dueteron mode
  amu    = Mu*1000.;
  fnucl  = Mp*1000.;
  fnucl2 = fnucl*fnucl;
  hbarc2 = hbarc*hbarc;
  ampj   =  (ML*1000.);
  ampj2  = ampj*ampj;
  amtag  = Mn*1000.;
  amtag2 = amtag*amtag;
  fmu = ampj*amtag/(ampj+amtag);
  fmunn = ampj*fnucl/(ampj+fnucl);
  xi=complex<double>(0.0,1.0);


  
  model = potential;
  //  cout<<" Potential "<<model<<endl;
  nrp1 = n1+n2;
  ftheL = complex<double>(0.0,0.0);

  
  //  skz = 100.0; // MeV test 
  skz = Prel; // MeV test 
    // Kinematics

 
    if(nr_mode){
      skz = skz/hbarc; // [1/fm]
      skz2 = skz*skz;  
      eon = hbarc2*skz2/2.0/fmu; // Kinematics Energy [MeV]
      rho = fmu*skz/hbarc2; // Desity of states/energy [1/(MeV*fm^3)]

    }else {
      skz2 = skz*skz;
      ekpj  = sqrt(ampj2   + skz2);
      ektag = sqrt( amtag2 + skz2);
      eon = ekpj + ektag;
      rho = skz*ekpj*ektag/eon;
    }
    
    // Get Gauss pts and wtss for integral v(q)Pl(x)dx (x=cos(theta))
  
    DGauss(xvq,wvq,npot);
    DGauss(x,w,n1);
    for( int i =0; i<n1;i++){ // DO 20
      skk[i] = skz  * (x[i]+1.0);
      wt[i]  = w[i] * skz;

    } // end 20
    DGauss(x,w,n2);
    for( int i =0; i<n2;i++){ // DO 21
      int  ii = i +n1;
      double coss = cos(PI*(x[i]+1.0)/4.0);
      skk[ii] = 2.0 * skz +2.0 *skz * tan(PI*(x[i]+1.0 )/4.0);
      wt[ii]  = 2.0* PI * 0.25 *skz/coss/coss*w[i];
      //      cout<<"i "<<i<<" skk "<<skk[i]<<" skz "<<skz*hbarc<<endl;
    } // end 21

    // Set Last Segment
    skk[nrp1] = skz;

    // start the l-loop Do 30
    // LL = the physical L . (0,1,2, ....)

    lpjmax =LL;


    

    //==============================================//
    //============= Loop Calc T-operator ===========//
    //==============================================//



    
    for(int ll=0;ll<=lpjmax;ll++){ // DO 30
      //      int l = ll+1;
      int l = ll;
      // start K',K'' loop to 40 & 50

      for(int ki =0; ki<=nrp1;ki++){  // Set skz


	skz = skk[ki]; // [1/fm]
	skz2 = skz*skz;  
	eon = hbarc2*skz2/2.0/fmu; // Kinematics Energy [MeV]
	rho = fmu*skz/hbarc2; // Desity of states/energy [1/(MeV*fm^3)]

	
	for(int kp=0;kp <= nrp1;kp++){ // DO 40
	  skp = skk[kp];
	  for(int kpp=0;kpp <= nrp1;kpp++){ // DO 50
	    skpp  = skk[kpp];
	    skpp2 = skpp*skpp;
	    // Initialization
	    gren[kp][kpp] = complex<double>(0.0,0.0);
	    ul[kp][kpp]   = complex<double>(0.0,0.0);
	    // ============== //
	    v2 = vlkkp(skp,skpp,ll);
	    ul[kp][kpp]  = v2;
	    eoff = hbarc2*skpp2/2.0/fmu;
	    if(!nr_mode)eoff=sqrt(ampj2+skpp2)+ sqrt(amtag2+skpp2);
	    if(kpp != nrp1)gren[kp][kpp] = wt[kpp]*skpp2*ul[kp][kpp]/(eon - eoff); 
	    if(kp==kpp)  delta = complex<double>(1.0,0.0);
	    else 	       delta = complex<double>(0.0,0.0);
	    
	    gren[kp][kpp] = delta - gren[kp][kpp];
	    //	    cout<<"kp "<<kp<<" kpp "<<kpp<<" gren "<<gren[kp][kpp]<<endl;
	    // gren : I - V * Go  
	  } // END 50
	} // END 40
	
	
	

	// Set the constant vector i.e. the B column matrix
	// in thet matrix equation AX = B
	
	
	for(int i=0;i<nmax2;i++){
	  v[i]  = ul[i][nrp1];  // DO 60
	  v0[i] = ul[i][nrp1]; 
	}
	
	
      int nmax3=(nrp1+1)*(nrp1+1);
      complex<double>A[nmax3],A0[nmax3];
      for(int i=0;i<=nrp1;i++){
	for(int j=0;j<=nrp1;j++){
	  
	  //	  if(i!=j) A0[i*(nrp1+1)+j] = complex<double>(0.0,0.0);
	  //	  else A0[i*(nrp1+1)+j]=gren[i][j];
	  A[i*(nrp1+1)+j] = gren[i][j];
	  //	  cout<<"i "<<i<<" j "<<j<<" A "<<gren[i][j]<<endl;
	}
      }
      
      // AX =B matrix : R(E-VGo)=V   (C.11) 
      // A : gren[i][j] = E -V * Go (C.11)
      // B : V
      // X : R =(E - V * Go)^-1 * V (A) = A^-1 * B -> retrun X :  R | k> half-off shell 

      
      
      LEQ3(A,v,Rl,nrp1);

      
      
      
      //==============================================//
      //=======< T matrix calculation >===============//
      //==============================================//
      
      
      // calc T matrix : T = R - i*pi*R*delta(E-Ho)T (C.12)
      // delta(E-Ho) = rho*delta(k-k')
      // <k|T = 1/(1+ i*rho*<k|R|k>)<k|R    (5.60)
      complex<double>T[nmax3]; // T-operator in 2-body scattering 
     
      for(int i=0;i<=nrp1;i++){
	
	T[i] = complex<double>(0.0,0.0);
	T[i] = 1.0/(1.0 + xi* rho * Rl[nrp1] * PI )* Rl[i];
	To[ki][i][ll] = T[i];
	//	cout<<"ll "<<ll<<" ki "<<ki<<" i "<<i<<" To "<<To[ki][i][ll]
	//	    <<" Rl "<<Rl[i]<<" Rl[nrp1 ]"<<Rl[nrp1]<<endl;
      }



      }
      



      
      // Calc Influence Factor
      // I =  <k| (1 + TGo)

      for(int kp=0;kp <= nrp1;kp++){ // DO 40
	skp = skk[kp];
	Ifac[kp][ll] = 0.0;
	for(int kpp=0;kpp <= nrp1;kpp++){ // DO 50
	  skpp  = skk[kpp];
	  skpp2 = skpp*skpp;
	  // Initialization
	  Ts[kp][kpp]    = complex<double>(0.0,0.0);
	  v2 = vlkkp(skp,skpp,ll);
	  eoff = hbarc2*skpp2/2.0/fmu;
	  if(!nr_mode)eoff=sqrt(ampj2+skpp2)+ sqrt(amtag2+skpp2);
	  if(kpp != nrp1)Ts[kp][kpp] = wt[kpp]*skpp2*To[kp][kpp][ll]/(eon - eoff); 
	  if(kp==kpp)  delta = complex<double>(1.0,0.0);
	  else 	       delta = complex<double>(0.0,0.0);

	  Ts[kp][kpp] = delta + Ts[kp][kpp];
	  Ifac[kp][ll]   += real( sqrt( real(Ts[kp][kpp])* conj(Ts[kp][kpp] )) );

	  // Ts : E - To (T-operator) * Go  
	} // END 50



	I0 =Ifac[kp][ll]; // Ifac = Ts*Ts

      } // END 40
    }  // end loop L




    
    return I0;
    
}

////////////////////////////////////////////////////////////////////////

/*
double fsi::CalcIfac_new(int LL, int potential, double Prel){

  cout<<"======================================="<<endl;
  cout<<"===  start Calc Influence Factor  ====="<<endl;
  cout<<"======================================="<<endl;
  
  //======= Energy Order is MeV ========//
  // fmu : dueteron mode
  amu    = Mu*1000.;
  fnucl  = Mp*1000.;
  fnucl2 = fnucl*fnucl;
  hbarc2 = hbarc*hbarc;
  ampj   =  (ML*1000.);
  ampj2  = ampj*ampj;
  amtag  = Mn*1000.;
  amtag2 = amtag*amtag;
  fmu = ampj*amtag/(ampj+amtag);
  fmunn = ampj*fnucl/(ampj+fnucl);
  xi=complex<double>(0.0,1.0);


  
  model = potential;
  nrp1 = n1+n2;
  ftheL = complex<double>(0.0,0.0);
  skz = Prel; // MeV test 


 
    if(nr_mode){
      skz = skz/hbarc; // [1/fm]
      skz2 = skz*skz;  
      eon = hbarc2*skz2/2.0/fmu; // Kinematics Energy [MeV]
      rho = fmu*skz/hbarc2; // Desity of states/energy [1/(MeV*fm^3)]

    }else {
      skz2 = skz*skz;
      ekpj  = sqrt(ampj2   + skz2);
      ektag = sqrt( amtag2 + skz2);
      eon = ekpj + ektag;
      rho = skz*ekpj*ektag/eon;
    }
    
    // Get Gauss pts and wtss for integral v(q)Pl(x)dx (x=cos(theta))
  
    DGauss(xvq,wvq,npot);
    DGauss(x,w,n1);
    for( int i =0; i<n1;i++){ // DO 20
      skk[i] = skz  * (x[i]+1.0);
      wt[i]  = w[i] * skz;

    } // end 20
    DGauss(x,w,n2);
    for( int i =0; i<n2;i++){ // DO 21
      int  ii = i +n1;
      double coss = cos(PI*(x[i]+1.0)/4.0);
      skk[ii] = 2.0 * skz +2.0 *skz * tan(PI*(x[i]+1.0 )/4.0);
      wt[ii]  = 2.0* PI * 0.25 *skz/coss/coss*w[i];
      //      cout<<"i "<<i<<" skk "<<skk[i]<<" skz "<<skz*hbarc<<endl;
    } // end 21

    // Set Last Segment
    skk[nrp1] = skz;

    // start the l-loop Do 30
    // LL = the physical L . (0,1,2, ....)

    lpjmax =LL;
    int ll = LL;

    


    
    //==============================================//
    //============= Loop Calc T-operator ===========//
    //==============================================//



      for(int ki =0; ki<=nrp1;ki++){  // Set skz


	skz = skk[ki]; // [1/fm]
	skz2 = skz*skz;  
	eon = hbarc2*skz2/2.0/fmu; // Kinematics Energy [MeV]
	rho = fmu*skz/hbarc2; // Desity of states/energy [1/(MeV*fm^3)]

      

	skz = skk[ki]; // [1/fm]
	skz2 = skz*skz;  
	eon = hbarc2*skz2/2.0/fmu; // Kinematics Energy [MeV]
	rho = fmu*skz/hbarc2; // Desity of states/energy [1/(MeV*fm^3)]

	
	for(int kp=0;kp <= nrp1;kp++){ // DO 40
	  skp = skk[kp];
	  for(int kpp=0;kpp <= nrp1;kpp++){ // DO 50
	    skpp  = skk[kpp];
	    skpp2 = skpp*skpp;
	    // Initialization
	    gren[kp][kpp] = complex<double>(0.0,0.0);
	    ul[kp][kpp]   = complex<double>(0.0,0.0);
	    // ============== //
	    v2 = vlkkp(skp,skpp,ll);
	    ul[kp][kpp]  = v2;
	    eoff = hbarc2*skpp2/2.0/fmu;
	    if(!nr_mode)eoff=sqrt(ampj2+skpp2)+ sqrt(amtag2+skpp2);
	    if(kpp != nrp1)gren[kp][kpp] = wt[kpp]*skpp2*ul[kp][kpp]/(eon - eoff); 
	    if(kp==kpp)  delta = complex<double>(1.0,0.0);
	    else 	       delta = complex<double>(0.0,0.0);
	    
	    gren[kp][kpp] = delta - gren[kp][kpp];
	    //	    cout<<"kp "<<kp<<" kpp "<<kpp<<" gren "<<gren[kp][kpp]<<endl;
	    // gren : I - V * Go  
	  } // END 50
	} // END 40
	
	
	

	// Set the constant vector i.e. the B column matrix
	// in thet matrix equation AX = B
	
	
	for(int i=0;i<nmax2;i++){
	  v[i]  = ul[i][nrp1];  // DO 60
	  v0[i] = ul[i][nrp1]; 
	}
	
	
      int nmax3=(nrp1+1)*(nrp1+1);
      complex<double>A[nmax3],A0[nmax3];
      for(int i=0;i<=nrp1;i++){
	for(int j=0;j<=nrp1;j++){
	  
	  //	  if(i!=j) A0[i*(nrp1+1)+j] = complex<double>(0.0,0.0);
	  //	  else A0[i*(nrp1+1)+j]=gren[i][j];
	  A[i*(nrp1+1)+j] = gren[i][j];
	  //	  cout<<"i "<<i<<" j "<<j<<" A "<<gren[i][j]<<endl;
	}
      }
      
      // AX =B matrix : R(E-VGo)=V   (C.11) 
      // A : gren[i][j] = E -V * Go (C.11)
      // B : V
      // X : R =(E - V * Go)^-1 * V (A) = A^-1 * B -> retrun X :  R | k> half-off shell 

      
      
      LEQ3(A,v,Rl,nrp1);

      
      
      
      //==============================================//
      //=======< T matrix calculation >===============//
      //==============================================//
      
      
      // calc T matrix : T = R - i*pi*R*delta(E-Ho)T (C.12)
      // delta(E-Ho) = rho*delta(k-k')
      // <k|T = 1/(1+ i*rho*<k|R|k>)<k|R    (5.60)
      complex<double>T[nmax3]; // T-operator in 2-body scattering 
     
      for(int i=0;i<=nrp1;i++){
	
	T[i] = complex<double>(0.0,0.0);
	T[i] = 1.0/(1.0 + xi* rho * Rl[nrp1] * PI )* Rl[i];
	To[ki][i][LL] = T[i];

      }


}
      



      
      // Calc Influence Factor
      // I =  <k| (1 + TGo)

      for(int kp=0;kp <= nrp1;kp++){ // DO 40
	skp = skk[kp];
	Ifac[kp][ll] = 0.0;
	for(int kpp=0;kpp <= nrp1;kpp++){ // DO 50
	  skpp  = skk[kpp];
	  skpp2 = skpp*skpp;
	  // Initialization
	  Ts[kp][kpp]    = complex<double>(0.0,0.0);
	  v2 = vlkkp(skp,skpp,ll);
	  eoff = hbarc2*skpp2/2.0/fmu;
	  if(!nr_mode)eoff=sqrt(ampj2+skpp2)+ sqrt(amtag2+skpp2);
	  if(kpp != nrp1)Ts[kp][kpp] = wt[kpp]*skpp2*To[kp][kpp][ll]/(eon - eoff); 
	  if(kp==kpp)  delta = complex<double>(1.0,0.0);
	  else 	       delta = complex<double>(0.0,0.0);

	  Ts[kp][kpp] = delta + Ts[kp][kpp];
	  Ifac[kp][ll]   += real( sqrt( real(Ts[kp][kpp])* conj(Ts[kp][kpp] )) );

	  // Ts : E - To (T-operator) * Go  
	} // END 50
	//	gIT[ll]->SetPoint(kp,skp*hbarc,Ifac[kp][ll]);
	//      gIT[ll]->Write();	
	//	cout<<"kp "<<kp<<" Ifac "<<Ifac[kp][ll]<<endl;
	I0 =Ifac[kp][ll];

      } // END 40





    
    return I0;
    
}

*/

///////////////////////////////////////////////////////////////////////

void fsi::GetIfac(int LL, int potentail, string pname){

  cout<<"======================================="<<endl;
  cout<<"===  start Calc Influence Factor  ====="<<endl;
  cout<<"======================================="<<endl;




  string pname1 = pname + "_Verma.dat";
  string pname2 = pname + "_JulichA.dat";
  string pname3 = pname + "_JulichB.dat";

  ofp[1] = new ofstream(pname1.c_str());
  ofp[2] = new ofstream(pname2.c_str());
  ofp[3] = new ofstream(pname3.c_str());


  
  *ofp[1] << "### FSI tabel with Verma Potential (2-Gauss potential) calculated by T-operator (on-shell approximation) #####"<<endl;
  *ofp[1]<<"# krel [MeV] # Influence Factor(S) # Influence Factor(T)"<<endl;

  *ofp[2] << "### FSI tabel with Julich A Potential (2-Gauss potential)  calculated by T-operator (on-shell approximation)#####"<<endl;
  *ofp[2]<<"# krel [MeV] # Influence Factor(S) # Influence Factor(T)"<<endl;

  *ofp[3] << "### FSI tabel with Julich B Potential (2-Gauss potential)  calculated by T-operator (on-shell approximation)#####"<<endl;
  *ofp[3]<<"# krel [MeV] # Influence Factor(S) # Influence Factor(T)"<<endl;
  

  bool test=true;
  test =false;
  int npmax =50;
  double Prel=0.0;          
  for(int i=1;i<=3;i++){ // Potential
    gIT[i]= new TGraphErrors();
    gIT[i]->SetName(Form("gIT_%d",i));    
    gIT[i]->SetMarkerStyle(7);
    gIT[i]->SetMarkerColor(i);
    gITs[i]= new TGraphErrors();
    gITs[i]->SetName(Form("gITs_%d",i));    
    gITs[i]->SetMarkerStyle(3);
    gITs[i]->SetMarkerColor(i);
    gITt[i]= new TGraphErrors();
    gITt[i]->SetName(Form("gITt_%d",i));    
    gITt[i]->SetMarkerStyle(3);
    gITt[i]->SetMarkerColor(i);
    Prel=0.0;
    if(test)npmax=5;  
    for(int p=0;p<npmax;p++){
      Prel += 3.*(double)(p+1);
      double Is = CalcIfac(0,i,Prel);
      double It = CalcIfac(0,-i,Prel);
      double I  = (Is+3.0*It)/4.0;
      gITs[i]->SetPoint(p,Prel,Is);
      gITt[i]->SetPoint(p,Prel,It);
      gIT[i] ->SetPoint(p,Prel,I );
      cout<<"potential "<<i<<" i "<<p<<" / "<<npmax<<" Prel "<<Prel<<" I "<<I<<endl;
      
      *ofp[i] << Prel << " " << Is <<" "<< It <<endl;
      
      
    }
    
    gITs[i]->Write();
    gITt[i]->Write();
    gIT[i] ->Write();
    ofp[i] ->close();      
  }
    
  
  
}


////////////////////////////////////////////////////////////////////////

complex<double> fsi::vlkkp(double skp,double skpp, int lcall){

  //  int l = lcall +1;
  int l = lcall;
  double xx=0.0;
  complex<double>vv(0.0,0.0);
  for(int i =0;i<npot;i++){
      
    q2 = skp*skp+skpp*skpp -2.0*skp*skpp*xvq[i]; //(skp -skpp)^2
    vq = vofq(q2);
    xx = xvq[i];
    PL(xx,l,pol);
    vv =  vv + vq *wvq[i]*pol[l];
  }
  
  vv = 2.0*PI*vv +xi*0.0;
  return vv;

};

////////////////////////////////////////////////////////////////////////

complex<double> fsi::vofq(double q2){
  
  complex<double> vqa,vqr,vq;

  vr   = 246.80; // MeV
  b_r  = 0.82;   // fm

  
 if(model==1){ // Verma model
    va    = -167.34; // MeV
    b_a   = 1.1;  //fm  
  }else if(model ==2){// julich A
    va   = -373.94; //  MeV 
    b_a  =  0.79;   //  fm
  }else if(model ==3){// julich B 
    va   = -131.49; //  MeV 
    b_a  =  1.095;   //  fm
 }else if(model==-1){
   va   = -132.42;
   b_a  = 1.1;
 }else if(model==-2){
   va   = -144.14;
   b_a  = 1.059;
 }else if(model==-3){
   va   = -189.60;
   b_a  = 0.964;
   
 }


 
 // test 

 if(Vscale){

   va  =   v_a;
   b_a =   beta_a;
   
 }
 
 /*

 if(model==1){ // Verma model
    va    = -167.34; // MeV
    b_a   = 1.1;  //fm  
  }else if(model ==2){// julich A
    va    = -167.34*1.5; // MeV
    b_a   = 1.1;  //fm  
  }else if(model ==3){// julich B
    va    = -167.34*2.0; // MeV
    b_a   = 1.1;  //fm     
 }else if(model==-1){
   va   = -132.42;
   b_a  = 1.1;
 }else if(model==-2){
   va   = -132.42*1.5;
   b_a  = 1.1;   
 }else if(model==-3){
   va   = -132.42*2.0;
   b_a  = 1.1;   
   
}
 
 */

  
  if(!nr_mode){
    b_a = b_a/hbarc;
    b_r = b_r/hbarc;
  }

  b_a2 = b_a*b_a;
  b_r2 = b_r*b_r;
  

  // Vq : Integral exp(-iqx)*V(x)
  // -> Fourier transform V(q)= sqrt(pi/a)^3*exp(-q2/4a)  
  
  vqa = va*pow(b_a,3)/(8.0*PI*PI)*sqrt(PI)/exp(q2*b_a2/4.0); // Attractive
  vqr = vr*pow(b_r,3)/(8.0*PI*PI)*sqrt(PI)/exp(q2*b_r2/4.0); // Repulsive
  vq = vqa + vqr  + xi*0.0;
  

  //  if(skz*hbarc==100 && model==1)cout<<"skz "<<skz*hbarc<<" q2 "<<q2<<" vq "<<vq<<" vqa "<<vqa<<" vqr "<<vqr<<endl;

  return vq;
};


///////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////
///////////////////////// Main func /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){


  int ch; char* mode="C";
  string ifname = "../rootfiles/simc/nnL_simc.root";
  string ofname = "./test.root";
  string ofPname ="./param/test";
  string pname ="./param/potential.list";
  bool write_mode = false;
  bool root_mode  = false;
  bool param_mode = false;
  int l=0;
  string L;
  extern char *optarg;

  while((ch=getopt(argc,argv,"h:p:s:r:w:f:s:l:TDHeISbcop"))!=-1){
    switch(ch){
            
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      root_mode =true;
      break;

    case 's':
      ifname = optarg;
      cout<<"input root filename : "<<ofname<<endl;
      root_mode =true;
      break;

    case 'S':
      cout<<" Vscale mode ON "<<endl;
      Vscale =true;
      break;

      
    case 'D':

      Deu =true;
      E09 =false;
      cout<<"Deuteron mode : "<<Deu<<endl;      
      break;

    case 'H':

      E09 =true;
      Deu =false;
      cout<<"3He  mode : "<<E09<<endl;      
      break;                  

    case 'T':

      E09 =false;
      Deu =false;
      cout<<"Tritium  mode : "<<endl;      
      break;                        
    case 'e':
      break;
    case 'p':
      pname = optarg;
      cout<<"input Param filename : "<<pname<<endl;
      param_mode = true;
      break;
      
    case 'r':
      ofname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;

      
    case 'w':
      ofPname  = optarg;
      cout<<"output new parameter filename : "<<ofPname<<endl;
      write_mode =true;
      break;


    case 'l':
      L = optarg;
      cout<<"L max : "<<L<<endl;
      l = atoi(L.c_str());
      
   case 'b':
      cout<<"BACH MODE!"<<endl;
      break;
      

    case 'h':
      cout<<"-f : input root  filename"<<endl;
      cout<<"-m : input matrix filename"<<endl;      
      cout<<"-w : output matrix filename"<<endl;
      cout<<"-r : output root filename"<<endl;
      cout<<"-o : output root & tuning file name "<<endl;

      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;

    }
  }




  
  string fermi_mom="./param/2h_AV18_K.dat";
  TClone = root_mode;
  cout<<"TClone "<<TClone<<endl;
  TApplication *theApp =new TApplication("App",&argc,argv);
  gSystem->Load("libMinuit");
  fsi* FSI =new fsi();
  FSI->JostParam();
  if(root_mode)  FSI -> SetBranch(ifname);
  FSI -> NewRoot(ofname);
  cout<<"param_mode "<<param_mode<<endl;
  if(param_mode) FSI -> SetParam(pname);
  FSI -> SetHist();
  FSI -> SetFermiMom(fermi_mom);
  if(write_mode && Ton_mode)  FSI ->GetIfac(0,1, ofPname); // T-operator calc.
  if(write_mode && !Vscale) FSI -> InfluenceFactor(ofPname,l);
  if(write_mode &&  Vscale) FSI -> InfluenceFactorVscale(ofPname,l);
  if(param_mode) FSI -> SetMom();

  cout<<endl;
  cout<<"======================================= "<<endl;
  if(root_mode)  cout<<"Input root file   : "<<ifname<<endl;
  if(param_mode) cout<<"Input param file  : "<<pname<<endl;
  if(write_mode) cout<<"OS;utput param file : "<<ofPname<<endl;
  cout<<"Output root  file : "<<ofname<<endl;;

  gSystem->Exit(1);
  theApp->Run(); 
  return 0;
  
}//end main



///////////////////////////////////////////////////////////////////
///////////// Function ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////





//####################################################

void DGauss(double* Y, double* WY, int N){

  
  double V,Z,Z1,P,P1,DP,DZ;
  //long long V,Z,Z1,P,P1,DP,DZ;
  double x,w,y,wy;
  int NMAX=200;
  double X[NMAX], W[NMAX],POL[NMAX+1];
  int NP1 = N+1;
  int M =(int)N/2;
  double polw[nmax2];
  int KMAX =100;
  Z = PI/(2.0*(double)N);
  //cout<<"DGAUS N "<<N<<" M "<<M<<endl;
  //  X[0] = 0.0;
  //  Y[0] = 0.0;
  //  PL(Y[0],N,polw);
  //  Y[0] = 2.0*(1.0-Y[0]*Y[0])/fabs((double(N+1)*polw[N]*double(N+1)*polw[N]));


      
  for(int L=1; L<=M;L++){    // DO 20
    for(int K=1;KMAX;K++){ // DO 5
      PL(Z,N,POL);
      P = POL[N];
      //  DP = (double)N* (Z*POL[NP1] - POL[N])/(Z*Z -1.0);
      DP = (double)N* (Z*POL[N] - POL[N-1])/(Z*Z -1.0);
      Z = Z- P/DP;
      if(abs(P) < 1.0e-12)break;
      if(K==KMAX-1)return;
    } // for K END 5
    
    X[L] = Z;
    //==========< calc weight > ====================//
    V = 2.0/( (1.0 -Z*Z)*DP*DP ); // original
    W[L] =V;
    //==============================================//

   
    
    if(L == M )break; // GOTO 30    
    DZ =Z;
    if(L >= 1)DZ = (X[L] - X[L-1])*0.5;
    for(int K=0;K<nmax2;K++){  // DO 17
      Z = Z+DZ;
      PL(Z,N,POL);
      P = POL[N];
      Z1 = Z  + DZ;
      PL(Z1,N,POL);
      P1 = POL[N];
      if(P*P1 < 0.0 )break; // TO 18
    } // END 17
    
    
    Z = (P1*Z -P*Z1)/(P1-P);
    
    
  } // for L END 20
  


  for(int NEG =1 ;NEG <N+1;NEG++){  // DO 40
    if(NEG <= M)Y[NEG-1] = - X[M+1-NEG];
    if(NEG >  M)Y[NEG-1] =   X[NEG-M];
    if(NEG <= M)WY[NEG-1] =  W[M+1-NEG];
    if(NEG >  M)WY[NEG-1] =  W[NEG-M];
    //    cout<<" NEG "<<NEG<<" Y "<<Y[NEG]<<" YW "<<WY[NEG]<<endl;
  } // END 40

  return ;
  
};

//#############################################################


void PL(double X, int L, double * POL){

  double S;
  POL[0] = 1.0;
  POL[1] = X;
  if(L <=1)return;
  for(int jj=1;jj<L;jj++){
    S=(double)(jj+1);
    POL[jj+1] =( (2.0*S -1.0)*X*POL[jj] -(S-1.0)*POL[jj-1] )/S;
  }
  //  cout<<"POL L "<<POL[L]<<" POL L+1 "<<POL[L+1]<<endl;
  return;
};


//################################################################


void PL2(double X, int L, double * POL){

  double S;
  double pol_L,pol_L1;
  //  if(X>1 )X =1.0;
  POL[0] = 1.0;
  POL[1] = X;
  if(L <=1)return;
  for(int jj=1;jj<=L;jj++){
    S=(double)(jj+1);

    POL[jj+1] =( (2.0*S -1.0)*X*POL[jj] -(S-1.0)*POL[jj-1] )/S;
    if(jj==L-1)   pol_L = POL[jj+1]*1.0e15;
    if(jj==L)   pol_L1 = POL[jj+1]*1.0e15;
  }
  cout<<"POL[L] "<<pol_L<<" L+1 "<<pol_L1<<" diff "<<pol_L-pol_L1<<endl;
  return;
};


//################################################################


complex<double>LEQ1(complex<double>*A, complex<double>* B,int N,int LA){


  // LEQ1 solves the matrix equation AX =B
  // A = Cofficient matrix, solution stored in B
  // B = constant matrix
  // N = number of equation s and unknowns
  // LA = dimension of the first subscript of A
  // **NOTE** the matrix A[i][j] must be symmetric
  // This routine was originally used in WIZARD

  
  complex<double> DET,R,T;
  double CABS,T1,T2;
  complex<double> t1,t2;
  int LNA = LA*N-LA;


  // ==================//

  for(int i=1;i<N;i++){ // DO 50
    int ii = i-1;
    int LJ = -LA;
    for(int j=0;j<=ii;j++){ // DO 51
      LJ = LJ + LA;
      int LIJ = i + LJ;
      int LJJ = j + LJ;
      t1 = A[LIJ]*conj(A[LIJ]);
      t2 = A[LJJ]*conj(A[LJJ]);
      T1 = real(t1);
      T2 = real(t2);
      //      cout<<"i "<<i<<" j "<<j<<" Aij "<<A[LIJ]<<" T2 "<<T2<<" T1 "<<T1<<endl;
      
      //      cout<<" j "<<j<<" LJ "<<LJ<<endl;
      if(T1 == 0.0) continue;
      if(T2 < T1){ // GO 10

	int LNAO = LNA+1;
	//	for(int LKO=0;LKO<LNAO;LKO++){ // DO 20
	for(int LKO=1;LKO<=LNAO;LKO = LKO + LA){ // DO 20
	  
	  int LK  =  LKO -1;
	  int LJK =  j + LK;
	  int LIK =  i + LK;
	  T       =  A[LJK];
	  A[LJK]  =  A[LIK];
	  A[LIK]  =  -T;
	  //	  cout<<" DO 20 : LKO "<<LKO<<" / "<<LNAO<<endl;
	} // END 20
	
	T    = B[j];
	B[j] = B[i];
	B[i] = -T;
	
      } // END 10
      // GO TO 30
      R = A[LIJ]/A[LJJ];
      int LLJ = LJ +LA;
      for(int LK = LLJ; LK<LNA; LK += LA){ // DO 40
	int	LIK = i +LK;
	int 	LJK = j +LK;
	A[LIK] = A[LIK] - R* A[LJK];
      }
	B[i] = B[i] - R*B[j];
	//	cout<<" i "<<i<<" j "<<j<<" Bi "<<B[i]<<endl; 
    } // END 51
  } // END 50

  
  int LNN = N + LNA;
  //int LNN = LNA;
  B[N] = B[N]/A[LNN];

  //  cout<<" B[N] "<<B[N]<<" A[LNN] "<<A[LNN]<<endl;
  //  cout<<" LNN "<<LNN<<endl;
  int LI = LNA;
  int LII =0;
  
  for(int ii=1;ii<N;ii++){ // DO 70
    LI = LI - LA;
    int i = N -ii +1;
    int KK = i+1;
    T =0.0;
    int LK = LA* i - LA;
    for(int k =KK;k<N;k++){ // DO 60

      LK = LK + LA;
      int LIK = i +LK;
      T =  T + A[LIK] * B[k];
      //      cout<<"Do 60 : k "<<k<<" / "<<N<<" T "<<T<<" A "<<A[LIK]<<" B "<<B[k]<<endl;
    }
    LII = i + LI;
    if(fabs( real(T) )<1.0e10)    B[i] = (B[i] - T )/A[LII];
    else B[i] =0.0;
    //    cout<<"i "<<i<<" B "<<B[i]<<" T "<<T<<" A "<<A[LII]<<endl;
  } // END 70
  DET = (1.0,0.0);
  
  for(int  j=1;j<N;j++){ // DO 80
    int LJJ = j + ( j-1)*LA;
    DET = DET*A[LJJ];
  }
  

  return DET;

}


//=============================================================================

complex<double>LEQ2(complex<double>*A, complex<double>* B, complex<double>* R,int N,int LA){


  // LEQ2 solves the matrix equation AX =B
  // A = Cofficient matrix, solution stored in B
  // B = constant matrix
  // N = number of equation s and unknowns
  // LA = dimension of the first subscript of A
  // **NOTE** the matrix A[i][j] must be symmetric
  // This routine was originally used in WIZARD

  complex<double> ron=0.0;

  // get index
  int nrow = N;
  int ncol = N;
  int ielem=0;
  MatrixXcd AA(ncol,nrow);
  VectorXcd BB(ncol);
  for(int icol =0; icol<ncol; icol++){
    for(int irow =0; irow<nrow; irow++){
      ielem = icol*nrow + irow;
      if(icol<N && irow<N)AA(icol,irow) = A[ielem];
      else AA(icol,irow) = 0.0;
    }
  }

  
  for(int icol =0;icol<N;icol++){
    if(icol<N) BB(icol) = B[icol];
    else  BB(icol) = 0.0;
  }
  VectorXcd X = AA.colPivHouseholderQr().solve(BB);
  //  cout<<" Solved X is "<<X<<endl;
  for(int icol =0;icol<N;icol++){
    //	B[icol] = X(icol);
    R[icol] = X(icol);
    //    cout<<" icol "<<icol<<" X "<<X(icol)<<" B "<<B[icol]<<endl;
	//    cout<<" icol "<<icol<<" B "<<B[icol]<<endl;
  }


  //  cout<<"N "<<N<<" N/4 "<<N/4<<" ron "<<X(N/4)<<" B "<<B[N/4]<<endl;
  //  ron = X(N/4);

  // on shell analysis //
  // Ron = Von + Vhalf *Go * Rhalf
  // ron = V(nrp1-1) + V(k',nrp1-1) * Go(k',nrp1-1) *X(k')
  //  ron = BB(nrp1-1) + AA(,nrp1-1)*X
 
  //  for(int i=0;i<N;i++)ron += AA(i,N/4)*X(i);
  //  ron += BB(N/4);
  //    cout<<" ron "<<ron<<" X "<<X(N-1)<<endl;

  return ron;

};


//##################################################################



void LEQ3(complex<double>*A, complex<double>* B, complex<double>* R,int N){


  // LEQ3 solves the matrix equation AX =B
  // A = Cofficient matrix, solution stored in B
  // B = constant matrix
  // N = number of equation s and unknowns
  // LA = dimension of the first subscript of A
  // **NOTE** the matrix A[i][j] must be symmetric
  // This routine was originally used in WIZARD

  complex<double> ron=0.0;

  // get index
  int nrow = N+1;
  int ncol = N+1;
  int ielem=0;
  MatrixXcd AA(ncol,nrow);
  VectorXcd BB(ncol);
  for(int icol =0; icol < ncol; icol++){
    for(int irow =0; irow < nrow; irow++){
      ielem = icol*nrow + irow;
      AA(icol,irow) = A[ielem];
      //      if(irow==N || icol==N)
      //      cout<<"icol "<<icol<<" irow "<<irow<<" A "<<A[ielem]<<endl;
    }
  }

    for(int icol =0;icol < ncol;icol++){
      BB(icol) = B[icol];
  }
  
  VectorXcd X = AA.colPivHouseholderQr().solve(BB);
  for(int icol =0;icol <= N;icol++){
    R[icol] = X(icol);
    
    //    cout<<" i "<<icol<<" Rl "<<R[icol]<<" V "<<B[icol]<<endl;
  }



};

//////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////

double Jl0(double p0,double a0, double r0){


  
  double  krel = p0/197.3269718;
  double alpha =  ( 1. + sqrt(1. - 2.*r0/a0)) / r0;
  double beta  =  (-1. + sqrt(1. - 2.*r0/a0)) / r0;
  complex<double> xi(0.0,1.0);
  complex<double> Jl = (krel - xi*beta)/(krel - xi*alpha);
  double I = 1./real(Jl* conj(Jl));
  
  return I;


};

double fJ(double *x ,double *par){

  double  krel = x[0]/197.3269718;
  double a0 = par[0];
  double r0 = par[1];
  double alpha =  ( 1. + sqrt(1. - 2.*r0/a0)) / r0;
  double beta  =  (-1. + sqrt(1. - 2.*r0/a0)) / r0;
  complex<double> xi(0.0,1.0);
  complex<double> Jl = (krel - xi*beta)/(krel - xi*alpha);
  double I = 1./real(Jl* conj(Jl));
  return I;  

};

double fJ2(double *x ,double *par){

  double  krel = x[0]/197.3269718;
  double as = par[0];
  double rs = par[1];
  double at = par[2];
  double rt = par[3];
  double alpha =  ( 1. + sqrt(1. - 2.*rs/as)) / rs;
  double beta  =  (-1. + sqrt(1. - 2.*rs/as)) / rs;
  complex<double> xi(0.0,1.0);
  complex<double> Jls = (krel - xi*beta)/(krel - xi*alpha);

  double alpha_t =  ( 1. + sqrt(1. - 2.*rt/at)) / rt;
  double beta_t  =  (-1. + sqrt(1. - 2.*rt/at)) / rt;
  complex<double> Jlt = (krel - xi*beta_t)/(krel - xi*alpha_t);  
  
  double I = 0.25/real(Jls* conj(Jls)) + 0.75/real(Jlt* conj(Jlt));

  //  cout<<"test "<<endl;
  return I;  

};


double rad_ERA(double *x ,double *par){

  double krel = x[0]/197.3269718;
  double a0   = par[0]; // scattering length
  double r0   = par[1]; // effective range


  // ERA :  k/cos(delta) = 1/a +r0k^2/2

  double delta = atan(krel/(-1./a0 +r0*krel*krel/2.));
   return delta;
  
};
