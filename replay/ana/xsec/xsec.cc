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
extern double calc_kaon_survival_ratio(TLorentzVector R_v, double pathL);
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


  cout<<"Charge : "<<Qc<<" [C] "<<endl;
  Ngamma = NGamma(Qc);
  cout<<"Ngamma : "<<Ngamma<<endl;
  Nt = NTarget(mode);
  cout<<"Nt : "<<Nt<<endl;
  cout<<"Nhyp "<<Nhyp<<endl;

  Mtar,Mhyp;
  if(mode=="H"){Mtar = Mp; Mhyp = ML;}
  else if(mode=="T"){Mtar = MT; Mhyp = MnnL;}
  else if(mode=="he"){Mtar = MHe3; Mhyp = MH3L;}

  
  
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
  //  Nhyp = Nhyp / double(nEk);
  /*  
  for(int iEk=0;iEk<nEk;iEk++){
    accept = 0.0;
    for(int ith=0;ith<nth_R;ith++){
      theta = theta_R_max*(double)ith/(double)nth_R;
      for(int iph=0;iph<nph_R;iph++){
	phi   = 2.0*PI*double(iph)/double(nph_R);
	accept += (dOmega_R(RE[iEk],theta)* sin(theta)*width_theta*width_phi);	
      }
    }
  */

  accept = 0.0;
  for(int iEk=0;iEk<nEk;iEk++){
    for(int ith=0;ith<nth_R;ith++)
      accept += dOmega_R(RE[iEk],theta);	
    if(accept>0)NHYP += Nhyp/accept;
		   
		   } // End Ek
  cout<<"width_phi "<<width_phi<<endl;
  cout<<"Nhyp "<<Nhyp<<" NHYP "<<NHYP<<endl;
  return NHYP;
  
}


double xsec::GetAccept_R(double Rp, double Rth){

  double accept =0.0;
  double srad,rad;
  double theta,phi;
  double width_Ek =(REMax -REMin)/double(nEk);
  double width_theta =theta_R_max/(double)nth_R;
  double width_phi =2.0*PI/(double)nph_R;
  int iEk, ith;
  double accept_total=0.0;
  for(int i=0;i<nEk;i++)
    if(Rp<RE[i] || i ==nEk-1){ iEk =i; break;}
  for(int i=0;i<nth_R;i++)
    accept_total += ROmega[iEk][i];

  
    //    if(Rth	< RTheta[i] || i==nth_R-1){
    //	  ith = i; width_theta = RTheta[i+1]- RTheta[i]; break;}
    //  accept = ROmega[iEk][ith];

  


  

  return accept_total;
  

  //  if(accept_total <=0.0) return 0.0;
  //  return accept/accept_total;  // [sr]



  
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
      Lp = sqrt(LE[iEe_]*LE[iEe_] - Me*Me);
      Lpx = Lp*sin(theta)*cos(phi);
      Lpy = Lp*sin(theta)*sin(phi);
      Lpz = Lp*cos(theta);
      LP_v.SetXYZ(Lpx,Lpy,Lpz);
      LP_v = HallACoordinate(LP_v,false);
      L_v.SetPxPyPzE(LP_v.X(),LP_v.Y(),LP_v.Z(),LE[iEe_]);
      VPF += VP_Flux(B_v,L_v)*dOmega_L(LE[iEe_],theta)*width_Ee_;

    }
  }
  

  Ngamma =VPF*(Qc/e);    
  //  cout<<"Virtual Photon Flux ; "<<VP_Flux(B_v,L_v) <<" [/GeV/sr]"<<endl;
  //  cout<<"Integral VP Flux : "<<VPF<<endl;

    
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
    if(iEk >= Ek || i==nEk-1){
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
	  return dOmega;
	}
      }
    }
  }


  
  return -1.0;
  
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

double xsec::Efficiency(string mode, TLorentzVector Rv, double RpathL){

  double eff_total=1.0;
  //  kaon_survival_ratio = 0.171;
  kaon_survival_ratio = calc_kaon_survival_ratio(Rv, RpathL);
  //  cout<<" survival ratio "<< kaon_survival_ratio<<endl;
  eff_total *=kaon_survival_ratio;
  //  kaon_eff = 0.9;
  eff_ac    = 0.59;  eff_total *= eff_ac;
  eff_vz    = 0.83;  eff_total *= eff_vz;
  eff_coin  = 0.97;  eff_total *= eff_coin;
  eff_track = 0.98;  eff_total *= eff_track;
  eff_chi   = 1.0;   eff_total *= eff_chi;
  eff_daq   = 0.95;  eff_total *= eff_daq;

  
  return eff_total;

}

///////////////////////////////////////////////////////////////////////

void xsec::Calc_XS(){
  double XS,XS2;
  double XS_CM;
  TLorentzVector R_v;
  R_v.SetPxPyPzE(0.0,0.0,Rp,sqrt(Rp*Rp+MK*MK));
  double RpathL=23.;
  double eff_total = Efficiency("H",R_v,23.);
  cout<<" Efficency : "<<eff_total<<endl;
  XS = NHyp(Nhyp)/Nt/Ngamma/eff_total*cm_to_barn*1.0e9;
  

  // XS calc CM coordinate
  TLorentzVector Rv, Lv, Bv, VF_v;
  Bv.SetPxPyPzE(0.0,0.0,Bp,sqrt(Bp*Bp + Me*Me));
  Lv.SetPxPyPzE(0.0,0.0,Lp_mean,sqrt(Lp_mean*Lp_mean + Me*Me));
  Lv.RotateX( - hrs_angle);
  Rv.SetPxPyPzE(0.0,0.0,Rp_mean,sqrt(Rp_mean*Rp_mean + MK*MK));
  Rv.RotateX( + hrs_angle);
  VF_v = Bv - Lv;
  XS_CM = XS*Lab_to_CM(Rv,VF_v);

  if(mode=="H")  XS2  = NHyp(Nhyp2)/Nt/Ngamma/eff_total *cm_to_barn*1.0e9;
  cout<<"Cross Section : "<<XS<<" [nb/sr]"<<endl;
  if(mode=="H")  cout<<"Cross Section(Sigma) : "<<XS2<<" [nb/sr]"<<endl;

  cout<<"Cross Section (CM): "<<XS_CM<<" [nb/sr]"<<endl;

  
}

/////////////////////////////////////////////////////////////////////////

double xsec::GetXS(TLorentzVector B_v, TLorentzVector L_v, TLorentzVector R_v){


  double XS;
  double XS_CM;
  double accept_R = 2.0*PI*(1.0 - cos(theta_R_max));
  // B_v : 4-vector in Beam 
  // L_v : 4-vector in LHRS LHRS coordinate
  // R_v : 4-vector in RHRS RHRS coordinate
  //  double Nt = NTarget(mode);
  //  double Ngamma = NGamma(Qc);
  double Bp,Bpx,Bpy,Bpz;
  double Rp,Rpx,Rpy,Rpz, Rth,Rph,Rz;
  double Lp,Lpx,Lpy,Lpz,Lth,Lph,Lz;
  double ER,EL,EB;
  Rp = R_v.P(); Rpx = R_v.Px(); Rpy = R_v.Py(); Rpz = R_v.Pz();
  Lp = L_v.P(); Lpx = L_v.Px(); Lpy = L_v.Py(); Lpz = L_v.Pz();
  Bp = B_v.P(); Bpx = B_v.Px(); Bpy = B_v.Py(); Bpz = B_v.Pz();
  ER = R_v.E(); EL = L_v.E();   EB = B_v.E();
  //  Rth = Rpz/Rpx; Rph = Rpz/Rpy;
  //  Lth = Lpz/Lpx; Lph = Lpz/Lpy;
  Rth = R_v.Theta(); Rph = R_v.Phi();  // [rad]
  Lth = L_v.Theta(); Lph = L_v.Phi();  // [rad]
  
  accept_R = GetAccept_R(Rp, Rth);
  if(accept_R<=0) return 0.0;
  double eff_total = Efficiency(mode, R_v,tr.RpathL);
  XS    = 1./Nt/Ngamma/accept_R/eff_total*cm_to_barn;

  //  cout<<"XS "<<XS<<" eff "<<eff_total<<" accept "<<accept_R<<" Ngamma "<<Ngamma<<endl;

  return XS; // [b/sr]
  
  
}



/////////////////////////////////////////////////////////////////////////


void xsec::SetBranch(string ifrname){

  add_tree(ifrname);
  pack_tree();
  readtreeHRSR();
  readtreeHRSL();

  tree->SetBranchStatus("Bp_c"     ,1  );
  tree->SetBranchStatus("Rp"     ,1  );
  tree->SetBranchStatus("Lp"     ,1  );
  tree->SetBranchStatus("Rth"      ,1  );
  tree->SetBranchStatus("Rph"      ,1  );
  tree->SetBranchStatus("Lth"      ,1  );
  tree->SetBranchStatus("Lph"      ,1  );
  tree ->SetBranchStatus("ntr_r",  1);
  tree ->SetBranchStatus("ntr_l",  1);
  tree ->SetBranchStatus("mm_nnL", 1);
  tree ->SetBranchStatus("pid_cut", 1);
  tree ->SetBranchStatus("z_cut",   1);
  tree ->SetBranchStatus("ct_cut",  1);
  tree ->SetBranchStatus("dpe",   1);
  tree ->SetBranchStatus("dpe_",  1);
  tree ->SetBranchStatus("dpk",   1);
  
  tree->SetBranchAddress("Bp_c"          ,&Bp_c );
  tree->SetBranchAddress("Lp"          ,&Lp_c );
  tree->SetBranchAddress("Rp"          ,&Rp_c );
  tree ->SetBranchAddress("Rth"          ,&Rth_c);
  tree ->SetBranchAddress("Lth"          ,&Lth_c);
  tree ->SetBranchAddress("Rph"          ,&Rph_c);
  tree ->SetBranchAddress("Lph"          ,&Lph_c);
  tree ->SetBranchAddress("ntr_r",&ntr_r);
  tree ->SetBranchAddress("ntr_l",&ntr_l);
  tree ->SetBranchAddress("mm_nnL",&mm_nnL);
  tree ->SetBranchAddress("pid_cut", &pid_cut);
  tree ->SetBranchAddress("z_cut",   &z_cut);
  tree ->SetBranchAddress("ct_cut",  &ct_cut);
  tree ->SetBranchAddress("dpe",     &dpe);
  tree ->SetBranchAddress("dpe_",    &dpe_);
  tree ->SetBranchAddress("dpk",     &dpk);
  
  
}


////////////////////////////////////////////////////////////////////////

void xsec::NewRoot(string ofrname){

  ofr  = new TFile(ofrname.c_str(),"recreate");
  Tnew = new TTree("T","T");
  Tnew = tree->CloneTree(0);
  Tnew ->Branch("mm",&tr.mm,"mm/D");
  Tnew ->Branch("xsec",&tr.xsec,"xsec/D");
  Tnew ->Branch("xsec_cm",&tr.xsec_cm,"xsec_cm/D");
  Tnew ->Branch("dOmage_R",&tr.dOmega_RHRS,"dOmega_R/D");


  //======= Make Hist ==========//
  int bin_mm;
  double min_mm,max_mm;
  min_mm =-100.;
  max_mm = 200;
  bin_mm=(int)(max_mm - min_mm); // 1 MeV Bin
  hmm =new TH1D("hmm","Missing Mass ; -B_{#Lambda} [MeV] ; Counts/2MeV",bin_mm,min_mm,max_mm);
  hmm_xsec =new TH1D("hmm_xsec","Missing Mass ; -B_{#Lambda} [MeV] ; d#sigma/d#Omega_{K} [(nb/sr) / 1MeV]",bin_mm,min_mm,max_mm);
  hmm_xsec_cm =new TH1D("hmm_xsec_cm","Missing Mass ; -B_{#Lambda} [MeV] ; d#sigma/d#Omega_{K}^{C.M.} [(nb/sr) / 1MeV]",bin_mm,min_mm,max_mm);



  /*
  bin_Rp =(int)(REMax -REMin)*1000; // Bin size MeV
  bin_Rth =(int)(theta_R_max -0.0)*1000; // bin size mrad

  hRp  = new TH1D("hRp","",bin_Rp,REMin,REMax);
  hRth = new TH1D("hRth","",bin_Rth,0.0,theta_R_max);
  hRp_Rth = new TH2D("hRp_Rth",bin_Rp,REMin,REMax,bin_Rth,0.0,theta_R_max);
  
  */
  
  
}



////////////////////////////////////////////////////////////////////////

void xsec::Loop(){

  cout<<"===================================="<<endl;
  cout<<"=============< Loop >==============="<<endl;
  cout<<"===================================="<<endl;



  
  int ENum = tree->GetEntries();
  
  for(int iev=0;iev<ENum;iev++){

   
    tree->GetEntry(iev);


	//===== Initialization =========//
	tr.mm   =-1000.;
	tr.xsec =-1000.;
	tr.RpathL =0.0;
		
	//==============================//
	Ee   = sqrt(Bp_c*Bp_c + Me*Me);
	L_E  = sqrt(Lp_c*Lp_c + Me*Me);
	R_E  = sqrt(Rp_c*Rp_c + MK*MK);

	
	if(ntr_r==1)
	  tr.RpathL = R_tr_pathl[0]+R_s2_trpath[0];
	else tr.RpathL =29.0;

	
	
	//===== Right Hand Coordinate ====//
	double R_pz = Rp_c/sqrt(1.0*1.0 + Rth_c* Rth_c + Rph_c* Rph_c );
	double R_px = R_pz * Rth_c ;
	double R_py = R_pz * Rph_c ;
	double L_pz = Lp_c/sqrt(1.0*1.0 +  Lth_c* Lth_c +  Lph_c* Lph_c);
	double L_px = L_pz * Lth_c ;
	double L_py = L_pz * Lph_c ;
	
	TVector3 L_v3,R_v3,B_v3;
	B_v3.SetXYZ(0.0,0.0,Bp_c);
	L_v3.SetXYZ(L_px,L_py,L_pz);
	R_v3.SetXYZ(R_px,R_py,R_pz);
	
	TLorentzVector L_v,R_v,B_v, T_v,MM_v;
	B_v.SetPxPyPzE(0.0,0.0,Bp,Ee);
	T_v.SetPxPyPzE(0.0,0.0,0.0,Mtar);
	L_v.SetPxPyPzE(L_px,L_py,L_pz,L_E);
	R_v.SetPxPyPzE(R_px,R_py,R_pz,R_E);
	
  
	// === Get XS =======//
	// HRS coordinate //
	tr.xsec = GetXS(B_v,L_v,R_v);
	tr.xsec *=1.0e9; // [b/sr] -> [nb/sr]
	tr.dOmega_RHRS = dOmega_RHRS;
	//===== Calc MM ======//
	// Rotate X HRS to Hall A coordinate
	R_v.RotateX( 13.2/180.*PI);
	L_v.RotateX(-13.2/180.*PI);
	MM_v = B_v + T_v - L_v -R_v;
	tr.mm = (MM_v.Mag() - Mhyp)*1000.; // -BL [MeV]
	tr.xsec_cm = Lab_to_CM(R_v, B_v - L_v); // CS in CM

	
	if(z_cut>0 && ct_cut>0 && pid_cut){
	hmm->Fill(tr.mm);
	hmm_xsec->Fill(tr.mm,tr.xsec);
       	hmm_xsec_cm->Fill(tr.mm, tr.xsec * tr.xsec_cm);
	}

	
       	Tnew->Fill();
	if(iev%(ENum/10)==0) cout<<"Filled "<<iev<<" / "<<ENum<<endl;
  } // END LOOP


  //======= Write =======//

  Tnew        -> Write();
  hmm         -> Write();
  hmm_xsec    -> Write();
  hmm_xsec_cm -> Write();
  
}


/////////////////////////////////////////////////////////////////////////
///////////////////////// Main func /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){


  int ch; char* mode="C";
  string ifname = "../rootfiles/momcalib/test.root";
  string ofname = "../rootfiles/momcalib/test.root";
  string ofMTPname ="";
  string pname ="./input/accept/Lam_data.in";
  extern char *optarg;

  while((ch=getopt(argc,argv,"h:r:w:f:s:bcop"))!=-1){
    switch(ch){
      
      
    case 'f':
      pname = optarg;
      cout<<"input filename : "<<pname<<endl;
      break;

    case 's':
      ifname = optarg;
      cout<<"input root filename : "<<ifname<<endl;      
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



  string name_LAcc ="./param/accept/LHRS_accept.param";
  string name_RAcc ="./param/accept/RHRS_accept.param";
  
  TApplication *theApp =new TApplication("App",&argc,argv);
  gSystem->Load("libMinuit");

  xsec* XS =new xsec();
  XS->SetParam(pname);
  cout<<"input param name : "<<pname<<endl;
  //  XS->SetLAccept(name_LAcc);
  //  XS->SetRAccept(name_RAcc);
  XS->Calc_XS();
  // double test=  XS->PhaseShift();
  XS -> SetBranch(ifname);
  XS -> NewRoot(ofname);
  XS -> Loop();
  gSystem->Exit(1);
  theApp->Run();

 
  return 0;
  
}//end main



///////////////////////////////////////////////////////////////////
///////////// Function ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////



double calc_kaon_survival_ratio(TLorentzVector R_v, double pathL){

  // calculation of kon survial ratio

  
  double kaon_lifetime = 1.238e-8; // [sec]
  double c             = 2.99792458e8;  // [sec] speed of light
  double ctau          = kaon_lifetime*c; // [m]
  double gamma_K = R_v.Gamma();  
  double survival_ratio = exp(-pathL/ctau/gamma_K) ;

  return survival_ratio ; 
  
};
