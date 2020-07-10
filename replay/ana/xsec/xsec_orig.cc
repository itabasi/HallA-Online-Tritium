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

extern void PL(double X, int L ,double * POL);
extern void PL2(double X, int L ,double * POL);
extern void DGauss(double* Y ,double *WY, int N);
extern complex<double>LEQ1(complex<double>* A, complex<double>* B,int N,int LA);
///////////////////////////////////////////////////////////////////////

void xsec::SetParam(string ifpname){

 
  Bp = 4.318; // Beam Energy [GeV]
  Lp = 2.100; // Scattered electron mom [GeV/c]
  Lp_mean = 2.100; // Scattered electron mom [GeV/c]
  Rp = 1.823;
  kaon_survival_ratio = 0.171;
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

  XS = NHyp(Nhyp)/Nt/Ngamma/kaon_survival_ratio*cm_to_barn*1.0e9;
  
  if(mode=="H")  XS2  = NHyp(Nhyp2)/Nt/Ngamma/kaon_survival_ratio *cm_to_barn*1.0e9;
  cout<<"Cross Section : "<<XS<<" [nb/sr]"<<endl;
  if(mode=="H")  cout<<" Cross Section2 : "<<XS2<<" [nb/sr]"<<endl;
  
}

/////////////////////////////////////////////////////////////////////////
/*
complex<double>vofq(double q2){


  if(model==1){
    va = -167.34; // MeV
    b_a = 1.1;    // fm
  }else(model==2){
      va  = -373,94;
      b_a = 0.79;
    }
  
  vr  = 246.80;
  b_r = 0.82;
  
  if(nr_mode){
    b_a = b_a/hbarc;
    b_r = b_r/hbarc; }
  

      b_a2 = b_a * b_a ;
      b_r2 = b_r * b_r ;

      nv =80;
      vqa = va*b_a*b_a*b_a/(8.0*PI*PI);
      vqr = vr*b_r*b_r*b_r/(8.0*PI*PI);
      vq  = vqa + vqr + xi*0.0;
      
      return vq ;
      
  
}

*/
/////////////////////////////////////////////////////////////////////////
double xsec::PhaseShift(){

  cout<<endl;
  cout<<" ========================================= "<<endl;
  cout<<" ========== Start PhaseShift ============== "<<endl;
  cout<<" ========================================== "<<endl;
  
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
  xi=(0.0,1.0);
  for(int ntst =1; ntst<=15;ntst++){

    skz = double(ntst)*10.;
    nrp1 = n1+n2;

    // Kinematics
 
    if(nr_mode){
      skz = skz/hbarc;
      skz2 = skz*skz;
      eon = hbarc*hbarc*skz2/2.0/fmu;
      rho = fmu*skz/hbarc2;
    }else {
      skz2 = skz*skz;
      ekpj  = sqrt(ampj2   + skz2);
      ektag = sqrt( amtag2 + skz2);
      eon = ekpj + ektag;
      rho = skz*ekpj*ektag/eon;
    }

    // Get Gauss pts and wtss for integral v(q)Pl(x)dx (x=cos(theta))

    //    npot =80;
    
    DGauss(xvq,wvq,npot);
    DGauss(x,w,n1);
    for( int i =0; i<n1;i++){ // DO 20
      skk[i] = skz  * (x[i]+1.0);
      wt[i]  = w[i] * w[i] * skz;
      //      cout<<"i "<<i<<" skk "<<skk[i]<<" skz "<<skz<<" x "<<x[i]<<" w "<<wt[i]<<endl;
    } // end 20

    DGauss(x,w,n2);
    for( int i =0; i<n2;i++){ // DO 21
      int  ii = i +n1;
      double coss = cos(PI*(x[i]+1.0)/4.0);
      skk[ii] = 2.0 * skz +2.0 *skz * tan(PI*(x[i]+1.0 )/4.0);
      wt[ii]  = 2.0* PI * 0.25 *skz/coss/coss*w[i]; 
      //     cout<<"i "<<i<<" skk "<<skk[i]<<" skz "<<skz<<" x "<<x[i]<<" w "<<wt[i]<<endl;
    } // end 21

    // Set Last Segment
    skk[nrp1] = skz;

    // start the l-loop Do 30
    // LL = the physical L . (0,1,2, ....)

    lpjmax =L_value;

    for(int ll=0;ll<=lpjmax;ll++){ // DO 30
      int l = ll+1;
      // start K',K'' loop to 40 & 50
      for(int kp=0;kp<nrp1;kp++){ // DO 40

	skp = skk[kp];
	gren[kp][nrp1] = 0.0;
	
	for(int kpp=0;kpp<nrp1;kpp++){ // DO 50
	  // Initialization
	  gren[kp][kpp] = 0.0;
	  ul[kp][kpp]   = 0.0;
	  // ============== //
	  v2 = vlkkp(skp,skp,ll);
	  ul[kp][kpp]  = v2;
	  skpp2 = skpp*skpp;
	  eoff = hbarc2*skpp2/2.0/fmu;
	  if(nr_mode)eoff=sqrt(ampj2+skpp2)+ sqrt(amtag2+skpp2);

	  delta =(0.0,0.0);
	  if(kpp == nrp1){
	    if(kp == kpp)delta=(1.0,0.0);
	    gren[kp][kpp] = delta - gren[kp][kpp];
	  }else{
	    gren[kp][kpp] = wt[kpp]*skpp2*ul[kp][kpp]/(eon - eoff);
	  }
	} // END 50
      } // END 40


      // Set the constant vector i.e. the B column matrix
      // in thet matrix equation AX = B


      for(int i=0;i<nrp1;i++)v[i] = ul[i][nrp1];  // DO 60



  // added by itabashi
  int nmax3=(nmax2+1)*(nmax2+1);
  complex<double>A[nmax3];
  for(int i=0;i<nmax2+1;i++){
    for(int j=0;j<nmax2+1;j++){
      A[i*(nmax2+1)+j]=gren[i][j];
      //      if(i<nrp1 && j<nrp1)      cout<<"i "<<i<<" j "<<j<<" A "<<gren[i][j]<<endl;
    }
  }

      
      det = LEQ1(A,v,nrp1,nmax2+1);
      //      det =1.0; // as a test
      //      cout<<"det "<<det<<endl;
      ron = v[nrp1];
      ronr = real(ron);
      delrad = tan(PI*rho*ronr);
      delta1[0] = delrad +xi*0.0;
      ton[0]    = ron/(1.0 +xi*PI*rho*ron);
      vv = vlkkp(skz,skz,ll);
      tborn[0] = vv;
      tonre = real(ton[0]);
      tonim = imag(ton[0]);
      tbornre = real(tborn[0]);
      tbornim = imag(tborn[0]);
      perre = abs(tbornre -tonre)/abs(tonre);
      perim = abs(tbornim -tonim)/abs(tonim);
      deldeg = delrad*180./PI;
      deldega = abs(deldeg);

      // test if tborn and ton are same within given tolerance
      // deffined at the top of program as tol
      if(perre <= tol && perim <= tol){
	lmax =ll;
	lmax1 = lmax -1;
	break;
      }
      
    } // END 30

    totalc = 0.0;

    // Now do the sum for the scattering amplitude
    // by using Born approximation, it essentially includes all partical waves

    double totalc =0.0;

    for(int ithe=0;ithe<180;ithe++){ // DO 80

      double the = double(ithe -1.0);
      double ther = PI* the/180.;
      double q  = 2.0 * skz * sin(ther/2.0);
      double q2 = q*q;
      xxx = cos(ther);
      PL(xxx, lmax1, pol);
      fthe[ithe]=(0.0,0.0);
      fthe1[ithe]=(0.0,0.0);      
      fthed =(0.0,0.0);

      for(int ll =0;ll<=lmax1;ll++){ // DO 90
	int l = ll+1;
	fthe[ithe]  = fthe[ithe] + (2.0*ll +1.0)*ton[l]*pol[l]/4.0/PI;
	fthe1[ithe] =  fthe1[ithe]+(2.0* ll + 1.0)*(ton[l] -tborn[l])*pol[l]/4.0/PI; 
      } // END 90

      
      fthe[ithe] = - fthe[ithe]*4.0*PI*PI*rho/skz;
      totalc = totalc + abs(fthe[ithe])*abs(fthe[ithe])*sin(ther)*PI/180.*2*PI;
      fthe1[ithe] = -(fthe1[ithe]+vq)*4.0*PI*PI*rho/skz;
      fborn = -vq*4.0*PI*PI*rho/skz;
      
    } // END 80

    tcross   = 4.0 * PI/skz * imag(fthe[0]);
    tcross1  = 4.0 * PI/skz*  imag(fthe1[0]);
   
  }// end ntst



  

}




////////////////////////////////////////////////////////////////////////

complex<double> xsec::vlkkp(double skp,double skpp, int lcall){

  
  int l = lcall +1;
  double xx=0.0;
  complex<double>vv(0.0,0.0);
  for(int i =0;i<npot;i++){

    q2 = skp*skp*skpp*skpp -2.0 * skp*skpp*xvq[i];
    vq = vofq(q2);
    xx = xvq[i];
    PL(xx,l,pol);
    vv =  vv + vq *wvq[i]*pol[i];
  }

  vv = vv*vv*PI +xi*0.0;

  return vv;

};

////////////////////////////////////////////////////////////////////////

complex<double> xsec::vofq(double q2){

  complex<double> vqa;
  
  if(model==1){ // Verma model
    va    = -167.34; // MeV
    b_a   = 1.1;  //fm  
  }else if(model ==2){// Hulich A
    va   = -373.94; //  MeV 
    b_a  =  0.79;   //  fm
  }

  vr   = 246.80; // MeV
  b_r  = 0.82;   // fm

  if(nr_mode){
    b_a = b_a/hbarc;
    b_r = b_r/hbarc;
  }

  b_a2 = b_a*b_a;
  b_r2 = b_r*b_r;

  int   nv = npot;
  vqa = pow(va*b_a,3)/(8.0*PI*PI)*sqrt(PI)/exp(q2*b_a2/4.0); // Attractive
  vqr = pow(vr*b_r,3)/(8.0*PI*PI)*sqrt(PI)/exp(q2*b_r2/4.0); // Repulsive
  vq = vqa + vqr  + xi*0.0;
  return vq;
};


///////////////////////////////////////////////////////////////////////


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
  //  XS->SetLAccept(name_LAcc);
  //  XS->SetRAccept(name_RAcc);
  XS->Calc_XS();
  XS->PhaseShift();
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
  X[0] = 0.0;
  Y[0] = 0.0;
  PL(Y[0],N,polw);
  Y[0] = 2.0*(1.0-Y[0]*Y[0])/fabs((double(N+1)*polw[N]*double(N+1)*polw[N]));


      
  for(int L=1; L<=M;L++){    // DO 20
    for(int K=1;KMAX;K++){ // DO 5
      PL(Z,N,POL);
      P = POL[NP1];
      //  DP = (double)N* (Z*POL[NP1] - POL[N])/(Z*Z -1.0);
      DP = (double)N* (Z*POL[NP1-1] - POL[N-1])/(Z*Z -1.0);
      Z = Z- P/DP;
      if(abs(P) < 1.0e-12){break;}
      if(K==KMAX-1){cout<<" faled DGaus "<<endl;return;}
    } // for K END 5

    cout<<" after break "<<" L "<<L<<" Z "<<Z<<" N "<<N<<" DP "<<DP<<" P "<<" DP-P "<<DP-P<<endl;
    
    X[L] = Z;
    //==========< calc weight > ====================//
    PL(Z,N,polw);
    V = 2.0/( (1.0 -Z*Z)*DP*DP ); // original
    //V = 2.0*(1.0 - Z*Z )/fabs((double(N+1)*polw[N]*double(N+1)*polw[N])) ;
    W[L] =V;
    cout<<"after calc w  Z "<<Z<<" N "<<N<<endl;
    //==============================================//

   
    
    if(L == M )break; // GOTO 30    
    DZ =Z;
    if(L >= 1)DZ = (X[L] - X[L-1])*0.5;
    //      cout<<"L "<<L<<" XL "<<X[L]<<" XL-1 "<<X[L-1]<<" DZ "<<DZ<<endl; 


    for(int K=0;K<nmax2;K++){  // DO 17
      Z = Z+DZ;
      PL(Z,N,POL);
      //      cout<<"k "<<K<<" POL[NP1] "<<POL[NP1]<<endl; P =
      P = POL[NP1];
      //      P = POL[NP1-1]; // itabashix
      //      P = POL[N];
      Z1 = Z  + DZ;
      PL(Z1,N,POL);
      P1 = POL[NP1];
      //      P1 = POL[NP1-1]; // itabashi
      //      P1 = POL[N-1];
      if(P*P1 < 0.0 ){     // cout<<"P "<<P<<" P1 "<<P1<<endl;
       	cout<<"K "<<K<<" Z "<<Z<<" Z1 "<<Z1<<" DZ "<<DZ<<" P "<<P<<" P1 "<<P1<<endl;
	break; }// TO 18

    } // END 17

    
      Z = (P1*Z -P*Z1)/(P1-P);

    
  } // for L END 20



  for(int NEG =1 ;NEG <N+1;NEG++){  // DO 40
    if(NEG <= M)Y[NEG] = - X[M+1-NEG];
    if(NEG >  M)Y[NEG] =   X[NEG-M];
    if(NEG <= M)WY[NEG] =  W[M+1-NEG];
    if(NEG >  M)WY[NEG] =  W[NEG-M];
    //    cout<<"  i "<<NEG<<" Y "<<Y[NEG]<<NEG<<" WY "<<WY[NEG]<<endl;
  } // END 40

  cout<<endl<<" END DGAUSS "<<endl;
  for(int NEG =0 ;NEG <N+1;NEG++)cout<<"  i "<<NEG<<" Y "<<Y[NEG]<<NEG<<" WY "<<WY[NEG]<<endl;
  cout<<endl;
  return ;
  
};

//#############################################################


void PL(double X, int L, double * POL){

  double S;
  //  if(X>1 )X =1.0;
  POL[0] = 1.0;
  POL[1] = X;
  if(L <=1)return;
  for(int jj=1;jj<=L;jj++){
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
    for(int j=0;j<ii;j++){ // DO 51

      LJ = LJ + LA;
      int LIJ = i + LJ;
      int LJJ = j + LJ;
      t1 = A[LIJ]*conj(A[LIJ]);
      t2 = A[LJJ]*conj(A[LJJ]);
      T1 = real(t1);
      T2 = real(t2);

      if(T1 == 0.0) continue;
      if(T2 < T1){ // GO 10

	int LNAO = LNA+1;

	//	for(int LKO=0;LKO<LNAO;LKO++){ // DO 20
	for(int LKO=0;LKO<LNAO;LKO = LKO + LA){ // DO 20
	  
	  int LK  =  LKO -1;
	  int LJK =  j + LK;
	  int LIK =  i + LK;
	  T       =  A[LJK];
	  A[LJK]  =  A[LIK];
	  A[LIK]  =  -T;
	  
	} // END 20
	
	T    = B[j];
	B[j] = B[i];
	B[i] = -T;
	
      } // END 10
      // GO TO 30
      R = A[LIJ]/A[LJJ];
      int LLJ = LJ +LA;
      for(int LK = LLJ; LK<LNA; LK += LA) // DO 40
	B[i] = B[i] - R*B[j];
    } // END 51
  } // END 50

  int LNN = N + LNA;
  B[N] = B[N]/A[LNN];
  int LI = LNA;
  int LII =0;
  
  for(int ii=1;ii<N;ii++){ // DO 70
    LI = LI - LA;
    int i = N -ii +1;
    int KK = i+1;
    T =0.0;
    int LK = LA* i - LA;
    LII =0; // added by itabashi
    for(int k =KK;k<N;k++){ // DO 60

      LK = LK + LA;
      int LIK = i +LK;
      T =  T + A[LIK] * B[k];
      LII = i + LI; 
    }
    B[i] = (B[i] - T )/A[LII];
  } // END 70

  for(int  j=1;j<N;j++){ // DO 80
    int LJJ = j + ( j-1)*LA;
    DET = DET*A[LJJ];
  }
  return DET;

};
