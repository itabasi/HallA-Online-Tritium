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
#include <TVector3.h>
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
#include "xsec.h"
#include "define.h"
#include "Param.h"
#include <TLorentzVector.h>
#include <TVector3.h>

///////////////////////////////////////////////////////////////////////

void xsec::SetParam(){


  Bp = 4.318; // Beam Energy [GeV]
  Lp = 2.100; // Scattered electron mom [GeV/c]
  Lp_mean = 2.100; // Scattered electron mom [GeV/c]
  Ee  = sqrt(Bp*Bp + Me*Me);
  Ee_ = sqrt(Lp*Lp + Me*Me);
  kaon_survival_ratio = 0.171;
  Qc =4.7; //[C]
  
}

//////////////////////////////////////////////////////////////////////

//double xsec::VP_Flux(double Q,double Ee, double Ee_,double theta_ee_){
double xsec::VP_Flux(TLorentzVector B_v, TLorentzVector L_v){

  // calcration of virtual photon flux //
  // Input parameters
  // B_v Beam 4-vector (px,py,pz,E)
  // L_v Scattered electron  4-vector (px,py,pz,E)
  // Out put vertual photon flux
  if(L_v.E()<=0.0)return 0.0;
  TLorentzVector VF_v; // Vertual Photon 4-vector
  VF_v = B_v - L_v;
  double nu =   VF_v.E(); // Ee - Ee_ [Gev]
  double Q2 = - VF_v.M2();  
  double q2 =   VF_v.P()*VF_v.P();  
  double theta_ee_ = B_v.Angle(L_v.Vect());
  double Ee  = B_v.E(); // Beam Energy [GeV]
  double Ee_ = L_v.E(); // Scattered electron Energy [GeV]
  //  double Ev  = (Ee - Ee_) + q2/(2.0*Mp);
  double Ev  = (Ee - Ee_) + Q2/(2.0*Mp);
  double epsilon = 1.0/(1.0 + 2.0*q2/Q2*tan(theta_ee_/2.0)*tan(theta_ee_/2.0));
  double Gamma = alpha /(2.0*PI*PI *Q2) * Ev/(1.0 - epsilon)*Ee_/Ee;
  //cout<<"Q2 "<<Q2<<" Ev "<<Ev<<" epsilon "<<epsilon<<" Ee "<<Ee<<" Ee_ "<<Ee_<<endl;
  //cout<<"Gamma "<<Gamma<<endl;
  //  Gamma = 3.2e-9;//test
  
  return Gamma;

}

/////////////////////////////////////////////////////////////////////


TVector3 xsec::HallACoordinate(TVector3 V, bool rhrs){

  // change coorinate spherical to Hall A
  // input sperical coordinate
  // V.X() : right array
  // V.Y() : gravity array
  // V.Z() : beam arrray
  // output Hall A coordinate
  // V.X() : gracity array
  // V.Y() : left array
  // V.Z() : beam arrray

  V.RotateZ(PI/2.0);
  if(rhrs) V.RotateX(+ hrs_angle);
  else     V.RotateX(- hrs_angle);
  
  return V;

}

//////////////////////////////////////////////////////////////////////
double xsec::NGamma(double Qc){

  // calculation of number of vertual photon
  // Qc : Total beam charge [C]
  double theta,phi;
  B_v.SetPxPyPzE(0.0,0.0,Bp,Ee);
  double Ngamma = 0.0;
  double width_Ee_   = fabs(LEMax-LEMin)/(double)nEe_;
  double width_theta = gen_theta_accep/(double)nth;
  double srad,rad;
  double VP=0.0;
  // check flag ////
  //defolt sum_theta : true , test : false
  bool sum_theta = true ;
  bool test      = false;
  //////////////////
  double Theta;
  double VPF=0.0;
  double width_phi = 2.0*PI/double(nph);
  TVector3 LP_v;
  
  for(int iEe_=0; iEe_<nEe_;iEe_++){
    if(sum_theta){
      srad = dOmega_L(LE[iEe_]);
      rad = acos(1.0-srad/(2.0*PI));
      width_theta =rad/(double)nth;
    }else       nth=1;
    for(int ith=0;ith<nth;ith++){
      for(int iph=0;iph<nph;iph++){
	if(sum_theta) theta = rad*(double)ith/(double)nth;
	else  theta =0.0;

	phi   = 2.0*PI*double(iph)/double(nph);
	Lp = sqrt(LE[iEe_]*LE[iEe_] - Me*Me);
	Lpx = Lp*sin(theta)*cos(phi);
	Lpy = Lp*sin(theta)*sin(phi);
	Lpz = Lp*cos(theta);
	LP_v.SetXYZ(Lpx,Lpy,Lpz);
	LP_v = HallACoordinate(LP_v,false);
       	L_v.SetPxPyPzE(LP_v.X(),LP_v.Y(),LP_v.Z(),LE[iEe_]);
	if(sum_theta)VPF += VP_Flux(B_v,L_v) *width_Ee_*width_theta*sin(theta)*width_phi;
	else  VPF += VP_Flux(B_v,L_v) *width_Ee_* dOmega_L(LE[iEe_]);

      }      
    }
  }

  

    theta = 0.0;
    phi =0.0;
    Lp  = Lp_mean;
    Lpx = Lp*sin(theta)*cos(phi);
    Lpy = Lp*sin(theta)*sin(phi);
    Lpz = Lp*cos(theta);
    LP_v.SetXYZ(Lpx,Lpy,Lpz);
    LP_v = HallACoordinate(LP_v,false);
    L_v.SetPxPyPzE(L_v.X(),L_v.Y(),L_v.Z(),Lp);
    if(test){
    VPF =0.0;
    VPF += VP_Flux(B_v,L_v) * 0.006*Lp_mean*(0.045*2.0);
    }

    Ngamma =VPF*(Qc/e);    

    cout<<"Virtual Photon Flux ; "<<VP_Flux(B_v,L_v) <<" [/GeV/sr]"<<endl;
    cout<<"Integral VP Flux : "<<VPF<<endl;

    
  return Ngamma;
  
}

////////////////////////////////////////////////////////////////////////

double xsec::dOmega_L(double Ee_){
  double iEe_;
  double dOmega;
  for(int i=0; i<nEe_;i++){    
    iEe_ = LE[i];
    if(iEe_ > Ee_){
      if(i== nEe_ -1)dOmega = LOmega[i];
      else if(i == 0)dOmega = LOmega[i];
      else dOmega = (LOmega[i] + LOmega[i-1])/2.0;
      break;  }

  }
  return dOmega;
  
}
/////////////////////////////////////////////////////////////////

double xsec::dOmega_R(double Ek){
  double iEk;
  double dOmega;
  for(int i=0; i<nEk;i++){    
    iEk = RE[i];
    if(iEk > RE[i]){
      if(i== nEk -1)dOmega = ROmega[i];
      else if(iEk == 0)dOmega = ROmega[i];
      else dOmega = (ROmega[i] + ROmega[i-1])/2.0;
      break;  }
  }
  return dOmega;
  
}

///////////////////////////////////////////////////////////////////////

double xsec::NTarget(string mode){

  double NT,rho,thickness; double A;
  if(mode=="T"){thickness = thickness_3H; A=3.0; 
  }else if(mode=="H"){thickness =thickness_1H; A=1.0;
  }else if(mode=="He"){thickness = thickness_3He; A=3.0;
  }
  
  NT = Na * thickness / A;

  return NT;

  


}

////////////////////////////////////////////////////////////////////////

void xsec::SetLAccept(string ifpname){

  cout<<endl;
  cout<<"========================================"<<endl;
  cout<<"===== Input LHRS Acceptance Param ======"<<endl;
  cout<<"========================================"<<endl;
  cout<<" input file : "<<ifpname<<endl;
  string buf;
  ifstream ifp(ifpname.c_str(),ios::in);
  nEe_ =0;
  LEMax=0.0;
  LEMin=1000.;
  if(ifp.fail()){cerr<<" Could not open file "<<endl; exit(1);}
  while(!ifp.eof() || nEe_ > nmax){
    getline(ifp,buf);
    if(buf[0]=='#')continue;
    stringstream sbuf(buf);
    sbuf >> LE[nEe_] >> LOmega[nEe_];
    if(LEMax<LE[nEe_])LEMax=LE[nEe_];
    if(LEMin>LE[nEe_] && LE[nEe_]>0)LEMin=LE[nEe_];
    //    cout<<" Ee_  "<<LE[nEe_]<<" LAccept "<<LOmega[nEe_]<<endl;
    nEe_++;
  }
  
}

////////////////////////////////////////////////////////////////////////

void xsec::SetRAccept(string ifpname){

  cout<<endl;
  cout<<"========================================"<<endl;
  cout<<"===== Input RHRS Acceptance Param ======"<<endl;
  cout<<"========================================"<<endl;  
  cout<<" input file : "<<ifpname<<endl;
  string buf;
  ifstream ifp(ifpname.c_str(),ios::in);
  nEk =0;
  if(ifp.fail()){cerr<<" Could not open file "<<endl; exit(1);}
  while(!ifp.eof() || nEk > nmax){
    getline(ifp,buf);
    if(buf[0]=='#')continue;
    stringstream sbuf(buf);
    sbuf >> RE[nEk] >> ROmega[nEk];
    //    cout<<" Ee_  "<<LE[nEe_]<<" LAccept "<<LOmega[nEe_]<<endl;
    nEk++;
  }
  
}

////////////////////////////////////////////////////////////////////////


void xsec::Calc_XS(){
  cout<<"Charge : "<<Qc<<" [C] "<<endl;
  double Ngamma = NGamma(Qc);
  cout<<"Ngamma : "<<Ngamma<<endl;
  double Nt = NTarget("H");
  cout<<"Nt : "<<Nt<<endl;


  double NL = 2000.;

  double  XS = NL/Nt/Ngamma/kaon_survival_ratio/0.006 *cm_to_barn*1.0e9;
  cout<<" Cross Section : "<<XS<<" [nb/sr]"<<endl;
  
}

////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////
///////////////////////// Main func /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){


  int ch; char* mode="C";
  string ifname = "../rootfiles/momcalib/test.root";
  string ofname = "../rootfiles/momcalib/test.root";
  string ofMTPname ="";
  extern char *optarg;

  while((ch=getopt(argc,argv,"h:s:w:f:s:bcop"))!=-1){
    switch(ch){
      
      
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;

    case 's':

      ifname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;      
      
    case 'r':

      ofname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;

      
    case 'w':
      ofMTPname  = optarg;
      cout<<"output new parameter filename : "<<ofMTPname<<endl;      
      break;

      
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



  string name_LAcc ="./param/LHRS_accept.param";
  string name_RAcc ="./param/test.param";
  
  TApplication *theApp =new TApplication("App",&argc,argv);
  gSystem->Load("libMinuit");

  xsec* XS =new xsec();
  XS->SetParam();
  XS->SetLAccept(name_LAcc);
  XS->SetRAccept(name_RAcc);
  XS->Calc_XS();
  
   gSystem->Exit(1);
   theApp->Run();

 
  return 0;
  
}//end main

