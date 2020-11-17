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
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;
using namespace Eigen;

///////////////////////////////////////////////////////////////////////

void xsec::SetParam(string ifpname){

 
  Bp = 4.318; // Beam Energy [GeV]
  Lp = 2.100; // Scattered electron mom [GeV/c]
  Lp_mean = 2.100; // Scattered electron mom [GeV/c]
  Rp = 1.823;
  kaon_survival_ratio = 0.171;
  kaon_eff = 0.94;
  Qc =4.7; //[C]

  ifstream ifp(ifpname.c_str());
  if(ifp.fail()){cout<<"Could not find file : "<<ifpname<<endl;exit(1);}

  
  string line;
  string name;
  string param;
  int i=0;
  while(1){
    getline(ifp,line);
    if(line[0]=='#')continue;
    if(ifp.eof())break;
    name = "mode:";          param =SetVal(line,name); if(param!="false"){ mode    = param; }
    name = "Beam momentum:"; param =SetVal(line,name);if(param!="false"){ Bp_mean = stod(param); }
    name = "RHRS momentum:"; param =SetVal(line,name);if(param!="false"){ Rp_mean = stod(param); }
    name = "LHRS momentum:"; param =SetVal(line,name);if(param!="false"){ Lp_mean = stod(param);  }
    name = "Beam Charge:";   param =SetVal(line,name);if(param!="false"){ Qc      = stod(param);  }
    name = "RHRS Acceptance file:";   param =SetVal(line,name);if(param!="false"){ RHRS_accept_file = param; }
    name = "LHRS Acceptance file:";   param =SetVal(line,name);if(param!="false"){ LHRS_accept_file = param; }
    name = "Number of Hypernuclei:";  param =SetVal(line,name);if(param!="false"){ Nhyp = stod(param);  }
    name = "Number of Hypernuclei2:"; param =SetVal(line,name);if(param!="false"){ Nhyp2 = stod(param);  }

  }

  Lp = Lp_mean;
  Rp = Rp_mean;
  Bp = Bp_mean;
  Ee  = sqrt(Bp*Bp + Me*Me);
  Ee_ = sqrt(Lp*Lp + Me*Me);
  Ek  = sqrt(Rp*Rp + MK*MK);
  SetRAccept(RHRS_accept_file);
  SetLAccept(LHRS_accept_file);

}

//////////////////////////////////////////////////////////////////

string xsec::SetVal(string line, string name){

  string val;
  if (line.compare(0,name.size(),name) !=0) return "false"; 
  
  string p = line.substr(name.size());
  string unit;
  stringstream sbuf(p);
  sbuf >> val >>unit;
  cout<<name<<" "<<val<<" "<<unit<<endl;
  return val;
  
}

//////////////////////////////////////////////////////////////////////


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

////////////////////////////////////////////////////////////////////

double xsec::NKaon(double Ek){
  

}

/////////////////////////////////////////////////////////////////////

double xsec::NHyp(double Nhyp){

  double srad,rad;
  double theta,phi;
  double width_Ek =(REMax -REMin)/double(nEk);
  double width_theta =theta_R_max/(double)nth_R;
  double width_phi =2.0*PI/(double)nph_R;
  double NHYP=0.0;
  double accept=0.0;


  //======== Uniform Nhyp(Ek) ==========//
  Nhyp = Nhyp / double(nEk);
  
  for(int iEk=0;iEk<nEk;iEk++){
    accept = 0.0;
    for(int ith=0;ith<nth_R;ith++){
      theta = theta_R_max*(double)ith/(double)nth_R;
      for(int iph=0;iph<nph_R;iph++){
	phi   = 2.0*PI*double(iph)/double(nph_R);
	accept += (dOmega_R(RE[iEk],theta)* sin(theta)*width_theta*width_phi);	
      }
    }
    
    if(accept>0)NHYP += Nhyp/accept;

  } // End Ek
  cout<<"width_phi "<<width_phi<<endl;
  cout<<"Nhyp "<<Nhyp<<" NHYP "<<NHYP<<endl;
  return NHYP;
  
}

/////////////////////////////////////////////////////////////////////

double xsec::Lab_to_CM(TLorentzVector R_v, TLorentzVector VF_v){


  TLorentzVector PT,Pcm,R_cm_v;
  double Mtar;

  if(mode=="H")Mtar = Mp;
  else if(mode=="T")Mtar = MT;
  else if(mode=="he")Mtar = MHe3;
  
  PT.SetPxPyPzE(0.0,0.0,0.0,Mtar);

  Pcm = VF_v + PT;
  TVector3 B;
  
  B = Pcm.BoostVector();
  double u  = R_v.P()/R_v.E();

  R_cm_v =R_v;
  R_cm_v.Boost(-B);


  double beta = B.Mag();
  double gamma = 1./sqrt(1.0 - beta*beta);
  double accept  = 6.0e-3; 
  double cos = 1.0 -accept/2./PI;
  double u_ = R_cm_v.P()/R_cm_v.E();
  double cos_ = (u*cos- beta)/(1.0 - beta*u*cos)/u_;
  double sin_ = sqrt(1.0 - cos_*cos_);
  double accept_cm = 2.0*PI*(1.0 - cos_);

  double XS_Lab_to_cm = gamma*(1.0 + cos_*beta/u_)/pow( sin_*sin_ + gamma*gamma*pow(cos_ + beta/u_ ,2.0)  ,3./2.);


  //  cout<<"VF "<<VF_v.P()<<" B "<<B.Mag()<<" Mtar "<<Mtar<<endl;
  //  cout<<"u "<<u<<" u_ "<<u_<<" cos "<<cos<<" cos_ "<<cos_<<endl;
  //  cout<<"R_v "<<R_v.P()<<" R_cm_v "<<R_cm_v.P()<<" beta "<<beta<<" fac "<<XS_Lab_to_cm<<endl;
  
  return XS_Lab_to_cm;

  
  
}


//////////////////////////////////////////////////////////////////////
double xsec::NGamma(double Qc){

  // calculation of number of vertual photon
  // Qc : Total beam charge [C]
  double theta,phi;
  B_v.SetPxPyPzE(0.0,0.0,Bp,Ee);
  double Ngamma = 0.0;
  double width_Ee_   = fabs(LEMax-LEMin)/(double)nEe_;
  //  double width_theta = gen_theta_accept/(double)nth;
  double width_theta = theta_L_max/(double)nth;
  double width_phi = 2.0*PI/double(nph);
  
  double srad,rad;  double VP=0.0;
  double Theta;
  double dOmega;
  double VPF=0.0;
  TVector3 LP_v;
  
  for(int iEe_=0; iEe_<nEe_;iEe_++){
    for(int ith=0;ith<nth;ith++){
      theta = theta_L_max*(double)ith/(double)nth;
      //      theta = gen_theta_accept*(double)ith/(double)nth;
      for(int iph=0;iph<nph;iph++){
	phi   = 2.0*PI*double(iph)/double(nph);
	Lp = sqrt(LE[iEe_]*LE[iEe_] - Me*Me);
	Lpx = Lp*sin(theta)*cos(phi);
	Lpy = Lp*sin(theta)*sin(phi);
	Lpz = Lp*cos(theta);
	LP_v.SetXYZ(Lpx,Lpy,Lpz);
	LP_v = HallACoordinate(LP_v,false);
       	L_v.SetPxPyPzE(LP_v.X(),LP_v.Y(),LP_v.Z(),LE[iEe_]);
	VPF += VP_Flux(B_v,L_v)*dOmega_L(LE[iEe_],theta)*sin(theta)*width_Ee_*width_theta*width_phi;
      }
    }
  }


  Ngamma =VPF*(Qc/e);    
  cout<<"Virtual Photon Flux ; "<<VP_Flux(B_v,L_v) <<" [/GeV/sr]"<<endl;
  cout<<"Integral VP Flux : "<<VPF<<endl;

    
  return Ngamma;
  
}

////////////////////////////////////////////////////////////////////////

double xsec::dOmega_L(double Ee_, double theta){
  double iEe_,ith;
  double dOmega=0.0;

  for(int i=0; i<nEe_;i++){    
    iEe_ = LE[i];
    if(iEe_ >= Ee_){
      if(theta > theta_L_max){
       	dOmega = LOmega[i][nth-1];
	return dOmega;
      }
      
      for(int j=0;j<nth;j++){
	ith = LTheta[j];
	if(ith >= theta){
	  if(i== nEe_-1)dOmega = LOmega[i][j];
	  else if(i == 0)dOmega = LOmega[i][j];
	  else if(LOmega[i][j] >0 && LOmega[i-1][j]>0)dOmega = (LOmega[i][j] + LOmega[i-1][j])/2.0;
	  else if(LOmega[i][j] >0)dOmega = LOmega[i][j];
	  return dOmega;
	}
      }
    }
  }



  
  return 0.0;
  
}

/////////////////////////////////////////////////////////////////


double xsec::dOmega_R(double Ek, double theta){
  double iEk,ith;
  double dOmega=0.0;

  for(int i=0; i<nEk;i++){    
    iEk = RE[i];
    if(iEk >= Ek){
      if(theta > theta_R_max){
       	dOmega = ROmega[i][nth-1];
	return dOmega;
      }
      
      for(int j=0;j<nth;j++){
	ith = RTheta[j];
	if(ith >= theta){
	  if(i== nEk-1)dOmega = ROmega[i][j];
	  else if(i == 0)dOmega = ROmega[i][j];
	  else if(ROmega[i][j] >0 && ROmega[i-1][j]>0)dOmega = (ROmega[i][j] + ROmega[i-1][j])/2.0;
	  else if(ROmega[i][j] >0)dOmega = ROmega[i][j];

	  //	  cout<<" Ek "<<RE[i]<<" theta "<<theta<<" dOmega "<<dOmega<<endl;
	  return dOmega;
	}
      }
    }
  }

  
  return 0.0;
  
}

///////////////////////////////////////////////////////////////////////

double xsec::NTarget(string mode){

  double NT,rho,thickness; double A;
  if(mode=="T"){thickness = thickness_3H; A=3.0; 
  }else if(mode=="H"){thickness =thickness_1H; A=1.0;
  }else if(mode=="3He"){thickness = thickness_3He; A=3.0;
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
  nth  =0;
  LEMax=0.0;
  LEMin=1000.;
  double theta,mom,acc;

  //  initilization
    for(int i=0;i<nmax;i++){
      LE[i]=0.0;
      LTheta[i]=0.0;
      for(int j=0;j<nmax;j++)LOmega[i][j]=0.0;

    }
  
  if(ifp.fail()){cerr<<" Could not open file "<<endl; exit(1);}
  while(!ifp.eof() || nEe_ > nmax){
    getline(ifp,buf);
    if(buf[0]=='#')continue;
    stringstream sbuf(buf);
    //    sbuf >> LE[nEe_] >> LTheta[nEe_] >> LOmega[nEe_][nth];
    sbuf >> mom >> theta >> acc;
    
    if(nEe_ ==0 && LE[nEe_]==0){
      LE[nEe_] = mom;
      nth=0;
    }else if(LE[nEe_] < mom){
      LE[nEe_+1] = mom;
      nEe_++;
      nth=0;
    }
    

    LTheta[nth]= theta;
    LOmega[nEe_][nth]  = acc;
    if(LEMax<LE[nEe_])LEMax=LE[nEe_];
    if(LEMin>LE[nEe_] && LE[nEe_]>0)LEMin=LE[nEe_];
    if(theta_L_max < LTheta[nth])theta_L_max = LTheta[nth];
    nth++; 
  }
  nth= nth -1;
  nEe_ = nEe_ +1;

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
  nth_R  =0;
  REMax=0.0;
  REMin=1000.;
  double theta,mom,acc;

  //  initilization
    for(int i=0;i<nmax;i++){
      RE[i]=0.0;
      RTheta[i]=0.0;
      for(int j=0;j<nmax;j++)ROmega[i][j]=0.0;

    }
  
  if(ifp.fail()){cerr<<" Could not open file "<<endl; exit(1);}
  while(!ifp.eof() || nEk > nmax){
    getline(ifp,buf);
    if(buf[0]=='#')continue;
    stringstream sbuf(buf);
    //    sbuf >> RE[nEk] >> RTheta[nEk] >> ROmega[nEk][nth_R];
    sbuf >> mom >> theta >> acc;
    
    if(nEk ==0 && RE[nEk]==0){
      RE[nEk] = mom;
      nth_R=0;
    }else if(RE[nEk] < mom){
      RE[nEk+1] = mom;
      nEk++;
      nth_R=0;
    }
    

    RTheta[nth_R]= theta;
    ROmega[nEk][nth_R]  = acc;
    if(REMax<RE[nEk])REMax=RE[nEk];
    if(REMin>RE[nEk] && RE[nEk]>0)REMin=RE[nEk];
    if(theta_R_max < RTheta[nth_R])theta_R_max = RTheta[nth_R];
    nth_R++; 
  }

  nth_R= nth_R -1;
  nEk = nEk+1;
}



////////////////////////////////////////////////////////////////////////


void xsec::Calc_XS(){
  cout<<"Charge : "<<Qc<<" [C] "<<endl;
  double Ngamma = NGamma(Qc);
  cout<<"Ngamma : "<<Ngamma<<endl;
  double Nt = NTarget(mode);
  cout<<"Nt : "<<Nt<<endl;
  cout<<"Nhyp "<<Nhyp<<endl;
  double  XS,XS2;
  double XS_CM;
  XS = NHyp(Nhyp)/Nt/Ngamma/kaon_survival_ratio/kaon_eff*cm_to_barn*1.0e9;
  

  // XS calc CM coordinate
  TLorentzVector Rv, Lv, Bv, VF_v;
  Bv.SetPxPyPzE(0.0,0.0,Bp,sqrt(Bp*Bp + Me*Me));
  Lv.SetPxPyPzE(0.0,0.0,Lp_mean,sqrt(Lp_mean*Lp_mean + Me*Me));
  Lv.RotateX( - hrs_angle);
  Rv.SetPxPyPzE(0.0,0.0,Rp_mean,sqrt(Rp_mean*Rp_mean + MK*MK));
  Rv.RotateX( + hrs_angle);
  VF_v = Bv - Lv;
  XS_CM = XS*Lab_to_CM(Rv,VF_v);


  if(mode=="H")  XS2  = NHyp(Nhyp2)/Nt/Ngamma/kaon_survival_ratio *cm_to_barn*1.0e9;
  cout<<"Cross Section : "<<XS<<" [nb/sr]"<<endl;
  if(mode=="H")  cout<<" Cross Section2 : "<<XS2<<" [nb/sr]"<<endl;

  cout<<"Cross Section (CM): "<<XS_CM<<" [nb/sr]"<<endl;

  
}

/////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////
///////////////////////// Main func /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){


  int ch; char* mode="C";
  string ifname = "../rootfiles/momcalib/test.root";
  string ofname = "../rootfiles/momcalib/test.root";
  string ofMTPname ="";
  string pname ="./param/Lam_data.in";
  extern char *optarg;

  while((ch=getopt(argc,argv,"h:s:w:f:s:bcop"))!=-1){
    switch(ch){
      
      
    case 'f':
      pname = optarg;
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
  XS->SetParam(pname);
  cout<<"input param name : "<<pname<<endl;
  //  XS->SetLAccept(name_LAcc);
  //  XS->SetRAccept(name_RAcc);
  XS->Calc_XS();
  // double test=  XS->PhaseShift();
  gSystem->Exit(1);
  theApp->Run();

 
  return 0;
  
}//end main



///////////////////////////////////////////////////////////////////
///////////// Function ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


