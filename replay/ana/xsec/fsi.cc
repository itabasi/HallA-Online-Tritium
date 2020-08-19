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
using Eigen::MatrixXd;
using namespace Eigen;
extern void PL(double X, int L ,double * POL);
extern void PL2(double X, int L ,double * POL);
extern void DGauss(double* Y ,double *WY, int N);
extern complex<double>LEQ1(complex<double>* A, complex<double>* B,int N,int LA);
extern complex<double>LEQ2(complex<double>* A, complex<double>* B, complex<double>* R,int N,int LA);
extern void LEQ3(complex<double>* A, complex<double>* B, complex<double>* R, int N);

bool SIMC    = true;
bool TClone  = false;
//bool E09     = false;
bool E09     = true;
bool Lam     = true;
bool single  = true;
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
    Tnew -> Branch("fac1_2",&fac1_2,"fac1_2/D");
    Tnew -> Branch("fac2_2",&fac2_2,"fac2_2/D");
    Tnew -> Branch("fac3_2",&fac3_2,"fac3_2/D");    
    Tnew -> Branch("Perl",&Prel,"Prel/D");
    Tnew -> Branch("Perl2",&Prel2,"Prel2/D");
    Tnew -> Branch("pL",&pL,"pL/D");
    Tnew -> Branch("pn",&pn,"pn/D");
    Tnew -> Branch("ranth",&ranth,"ranth/D");
    Tnew -> Branch("ranph",&ranph,"ranph/D");
    Tnew -> Branch("pLx",&pLx,"pLx/D");
    Tnew -> Branch("pLy",&pLy,"pLy/D");
    Tnew -> Branch("pLz",&pLz,"pLz/D");
    Tnew -> Branch("pnx",&pnx,"pnx/D");
    Tnew -> Branch("pny",&pny,"pny/D");
    Tnew -> Branch("pnz",&pnz,"pnz/D");
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
  string vname[10];
  ifstream ifp(ifpname.c_str(),ios::in);
  if(ifp.fail()){cout<<"Could not find Files "<<ifpname<<endl;exit(1);}

  while(!ifp.eof()){
    getline(ifp,buf);
    if(buf[0]=='#')continue;
    stringstream sbuf(buf);
    sbuf >> vname[s];
    cout<<" file : "<<vname[s]<<endl;
    s++;
  }

  
  ifstream ifp1(vname[0].c_str(),ios::in);
  ifstream ifp2(vname[1].c_str(),ios::in);
  ifstream ifp3(vname[2].c_str(),ios::in);


  int i1=0;
  
  while(!ifp1.eof()){
    getline(ifp1,buf);
    if(buf[0]=='#' || buf.length()==0 )continue;
    stringstream sbuf(buf);
    sbuf >> krel1[i1] >> w1[i1];
    i1++;
  }


  np1=i1;
  
  int i2=0;
  while(!ifp2.eof()){
    getline(ifp2,buf);
    stringstream sbuf(buf);
    if(buf[0]=='#' ||  buf.length()==0 )continue;
    sbuf >> krel2[i2] >> w2[i2];
    //    cout<<"i "<<i2<<" krel "<<krel2[i2]<<" w2 "<<w2[i2]<<endl;
    i2++;
  }

  np2 = i2;
  
  int i3=0;
  while(!ifp3.eof()){
    getline(ifp3,buf);
    if(buf[0]=='#' || buf.length()==0 )continue;
    stringstream sbuf(buf);
    sbuf >> krel3[i3] >> w3[i3];
    i3++;
  }
  
  np3 =i3;
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
  
  string pname1t = pname + "_Verma_T.dat";
  string pname2t = pname + "_JulichA_T.dat";
  string pname3t = pname + "_JulichB_T.dat";

  ofp0 = new ofstream(pname0.c_str());
  ofp1 = new ofstream(pname1.c_str());
  ofp2 = new ofstream(pname2.c_str());
  ofp3 = new ofstream(pname3.c_str());

  ofp1t = new ofstream(pname1t.c_str());
  ofp2t = new ofstream(pname2t.c_str());
  ofp3t = new ofstream(pname3t.c_str());  


  *ofp0 << "### FSI tabel with 3He Potential (2-Gauss potential) #####"<<endl;
  *ofp0<<"# Weight # krel [MeV] "<<endl;
  
  *ofp1 << "### FSI tabel with Verma Potential (2-Gauss potential) #####"<<endl;
  *ofp1<<"# Weight # krel [MeV] "<<endl;

  *ofp2 << "### FSI tabel with Julich A Potential (2-Gauss potential) #####"<<endl;
  *ofp2<<"# Weight # krel [MeV] "<<endl;

  *ofp3 << "### FSI tabel with Julich B Potential (2-Gauss potential) #####"<<endl;
  *ofp3<<"# Weight # krel [MeV] "<<endl;


  *ofp1t << "### FSI tabel with Verma Potential Triplet (2-Gauss potential) #####"<<endl;
  *ofp1t<<"# Weight # krel [MeV] "<<endl;

  *ofp2t << "### FSI tabel with Julich A Potential Triplet (2-Gauss potential) #####"<<endl;
  *ofp2t<<"# Weight # krel [MeV] "<<endl;

  *ofp3t << "### FSI tabel with Julich B Potential Triplet (2-Gauss potential) #####"<<endl;
  *ofp3t<<"# Weight # krel [MeV] "<<endl;  

  
  double qmin = 0.0;
  double qmax = 2000;
  int imax = 1000;

  if(test){

    qmax =1000;
    imax = 50;
  }
  
  for(int i=1; i<imax;i++){
  
    double qi = (qmax - qmin)/(double)imax * (double)i +qmin;

    
    fac0 = PhaseShift(qi, ll, 0);
    gI00   -> SetPoint(ii,qi,I0);
    gI0    -> SetPoint(ii,qi,I1);
    grad00 -> SetPoint(ii,qi,delrad0);
    fac0 = I0;
    
    fac1 = PhaseShift(qi, ll, 1);
    gx1    -> SetPoint(ii,qi,tcross0);
    gxf1   -> SetPoint(ii,qi,tcross);
    gI01   -> SetPoint(ii,qi,I0);
    gI1    -> SetPoint(ii,qi,I1);
    grad1  -> SetPoint(ii,qi,delrad);
    grad01 -> SetPoint(ii,qi,delrad0);
    gradI1 -> SetPoint(ii,delrad0,fac1);

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
      //      gfth1->SetPoint(th,th,real( fthe[th]*conj(fthe[th]) ));
    //real(ftheL* conj(ftheL) )+1.0;
    
    fac2 = PhaseShift(qi, ll, 2);
    gx2    ->SetPoint(ii,qi,tcross0);
    gxf2   ->SetPoint(ii,qi,tcross);
    gI02   ->SetPoint(ii,qi,I0);
    gI2    ->SetPoint(ii,qi,I1);
    grad2  ->SetPoint(ii,qi,delrad);
    grad02 ->SetPoint(ii,qi,delrad0);
    gradI2 ->SetPoint(ii,delrad0,fac2);

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
      //      gfth2->SetPoint(th,th,real( fthe[th]*conj(fthe[th]) ));
	
    fac3 = PhaseShift(qi, ll, 3);
    gx3    -> SetPoint(ii,qi,tcross0);
    gxf3   -> SetPoint(ii,qi,tcross);    
    gI03   -> SetPoint(ii,qi,I0);
    gI3    -> SetPoint(ii,qi,I1);
    grad3  -> SetPoint(ii,qi,delrad);
    grad03 -> SetPoint(ii,qi,delrad0);    
    gradI3 -> SetPoint(ii,delrad0,fac3);

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


    fac1t = PhaseShift(qi, ll, -1);
    fac2t = PhaseShift(qi, ll, -2);
    fac3t = PhaseShift(qi, ll, -3);

    fac1t = (fac1 + 3*fac1t)/4.;
    fac2t = (fac2 + 3*fac2t)/4.;
    fac3t = (fac3 + 3*fac3t)/4.;	
    
    gf1 -> SetPoint(ii,qi,fac1);
    gf2 -> SetPoint(ii,qi,fac2);
    gf3 -> SetPoint(ii,qi,fac3);

    gf1t -> SetPoint(ii,qi,fac1t);
    gf2t -> SetPoint(ii,qi,fac2t);
    gf3t -> SetPoint(ii,qi,fac3t);    
    
   *ofp0 << qi << " " << fac0 <<endl;
   *ofp1 << qi << " " << fac1 <<endl;
   *ofp2 << qi << " " << fac2 <<endl;
   *ofp3 << qi << " " << fac3 <<endl;

   *ofp1t << qi << " " << fac1t <<endl;
   *ofp2t << qi << " " << fac2t <<endl;
   *ofp3t << qi << " " << fac3t <<endl;
   
    ii++;
  }
  
  gf0->Write();
  gf1->Write();
  gf2->Write();
  gf3->Write();
  gf1t->Write();
  gf2t->Write();
  gf3t->Write();  
  gxf1->Write();
  gxf2->Write();
  gxf3->Write();
  gx1->Write();
  gx2->Write();
  gx3->Write();
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

  ofp0->close();
  ofp1->close();
  ofp2->close();
  ofp3->close();

  ofp1t->close();
  ofp2t->close();
  ofp3t->close();  
  
  
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
  gI0->SetTitle("3He w FSI Influence Factor; P_{#Lambda-n} [MeV] ; Influence Factor");
  gI0->SetMarkerStyle(2);
  gI0->SetMarkerColor(8);  
  
  
  gI1 = new TGraphErrors();
  gI1->SetName("gI1");
  gI1->SetTitle("Verma Potentail w FSI Influence Factor; P_{#Lambda-n} [MeV] ; Influence Factor");
  gI1->SetMarkerStyle(2);
  gI1->SetMarkerColor(1);

  gI2 = new TGraphErrors();
  gI2->SetName("gI2");
  gI2->SetTitle("Julich A w FSI Influence Factor; P_{#Lambda-n} [MeV] ; Influence Factor");
  gI2->SetMarkerStyle(2);
  gI2->SetMarkerColor(2);

  gI3 = new TGraphErrors();
  gI3->SetName("gI3");
  gI3->SetTitle("Julich B w FSI Influence Factor; P_{#Lambda-n} [MeV] ; Influence Factor");
  gI3->SetMarkerStyle(2);
  gI3->SetMarkerColor(4);
 

  
  gI00 = new TGraphErrors();
  gI00->SetName("gI00");
  gI00->SetTitle("3He potential w FSI Influence Factor; P_{#Lambda-n} [MeV] ; Influence Factor");
  gI00->SetMarkerStyle(4);
  gI00->SetMarkerColor(3);
  
  gI01 = new TGraphErrors();
  gI01->SetName("gI01");
  gI01->SetTitle("Verma Potentail w FSI Influence Factor; P_{#Lambda-n} [MeV] ; Influence Factor");
  gI01->SetMarkerStyle(4);
  gI01->SetMarkerColor(1);

  gI02 = new TGraphErrors();
  gI02->SetName("gI02");
  gI02->SetTitle("Julich A w FSI Influence Factor; P_{#Lambda-n} [MeV] ; Influence Factor");
  gI02->SetMarkerStyle(4);
  gI02->SetMarkerColor(2);

  gI03 = new TGraphErrors();
  gI03->SetName("gI03");
  gI03->SetTitle("Julich B w FSI Influence Factor; P_{#Lambda-n} [MeV] ; Influence Factor");
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
  int bin_mm = (int)(max_mm-min_mm)/2;
  
  hmm = new TH1F("hmm","Missing Mass w/o FSI ; -B_{#Lambda} [MeV] ; Counts/1MeV",bin_mm,min_mm,max_mm);

  hmm ->SetLineColor(1);
  
  hmm_fsi1 = new TH1F("hmm_fsi1","Missing Mass w/ FSI (Verma potenrial); -B_{#Lambda} [MeV] ; Counts/1MeV",bin_mm,min_mm,max_mm);
  hmm_fsi1 ->SetLineColor(6);
  
  hmm_fsi2 = new TH1F("hmm_fsi2","Missing Mass w/ FSI (Verma potenrial); -B_{#Lambda} [MeV] ; Counts/1MeV",bin_mm,min_mm,max_mm);
  hmm_fsi2 ->SetLineColor(2);
  hmm_fsi3 = new TH1F("hmm_fsi3","Missing Mass w/ FSI (Verma potenrial); -B_{#Lambda} [MeV] ; Counts/1MeV",bin_mm,min_mm,max_mm);
  hmm_fsi3 ->SetLineColor(4);


  hmm_fsi1_2 = new TH1F("hmm_fsi1_2","Missing Mass (Calc B) w/ FSI (Verma potenrial); -B_{#Lambda} [MeV] ; Counts/1MeV",bin_mm,min_mm,max_mm);
  hmm_fsi1_2 ->SetLineColor(6);
  hmm_fsi1_2 ->SetLineStyle(10);
  
  hmm_fsi2_2 = new TH1F("hmm_fsi2_2","Missing Mass (Calc B) w/ FSI (Verma potenrial); -B_{#Lambda} [MeV] ; Counts/1MeV",bin_mm,min_mm,max_mm);
  hmm_fsi2_2 ->SetLineColor(2);
  hmm_fsi2_2 ->SetLineStyle(10);  

  hmm_fsi3_2 = new TH1F("hmm_fsi3_2","Missing Mass (Calc B) w/ FSI (Verma potenrial); -B_{#Lambda} [MeV] ; Counts/1MeV",bin_mm,min_mm,max_mm);
  hmm_fsi3_2 ->SetLineColor(4);
  hmm_fsi3_2 ->SetLineStyle(10);  
  
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
    fac1_2  = 0.0;
    fac2_2  = 0.0;
    fac3_2  = 0.0;
    
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

    Pn_2 = CalcPn2(P2n_2);    
    
    //    Pn.SetPxPyPzE(Pnx,Pny,Pnz,En);
    //    Prel = CalcQ(Pn,Pp,APv,PK); // MeV

    TLorentzVector Pp_2,Pv_2,PK_2;
    Pp_2 =Pp;
    Pv_2 =Pv;
    PK_2 =PK;
    
    Prel = CalcQ(Pn,Pp,Pv,PK); // MeV
    Prel2 = CalcQ2(Pn_2,Pp_2,Pv_2,PK_2); // MeV  New version 

    //    cout<<" Pn "<<Pn.P()<<" Pn_2 "<<Pn_2.P()<<" Perl "<<Prel<<" Prel2 "<<Prel2<<endl;
    
    
    for(int i=0;i<np1;i++){
      //      cout<<" i "<<i<<" Prel "<<Prel<<" krel1 "<<krel1[i]<<endl;
      if(Prel < krel1[i]){
	if(i==0 || i==np1-1)fac1 = w1[i];
	else {
	  double a = (w1[i+1] -w1[i])/(krel1[i+1] - krel1[i]);
	  double b = w1[i] - a*krel1[i];
	  fac1   =  a*Prel  + b;}
	break;
      }else if(i==np1-1)fac1 = w1[i];
    } // end fac1
    
    
    for(int i=0;i<np2;i++){
      if(Prel < krel2[i]){
	if(i==0 || i==np2-1)fac2 = w2[i];
	else {
	  double a = (w2[i+1] -w2[i])/(krel2[i+1] - krel2[i]);
	  double b = w2[i] - a*krel2[i];
	  fac2 = a*Prel + b;
	  //	  cout<<" n "<<n<<" fac2 "<<fac2<<endl;
	}
	break;
      }else if(i==np2-1)fac2 = w2[i];
    } // end fac2
    
     
    for(int i=0;i<np3;i++){
      if(Prel < krel3[i]){
	if(i==0 || i==np3-1)fac3 = w3[i];
	else {
	  double a = (w3[i+1] -w3[i])/(krel3[i+1] - krel3[i]);
	  double b = w3[i] - a*krel3[i];
	  fac3 = a*Prel + b;  }
	break;
      }else if(i==np3-1)fac3 = w3[i];
    } // end fac1
    
    


    //  Prel2  //

    for(int i=0;i<np1;i++){
      //      cout<<" i "<<i<<" Prel2 "<<Prel2<<" krel1 "<<krel1[i]<<endl;
      if(Prel2 < krel1[i]){
	if(i==0 || i==np1-1)fac1_2 = w1[i];
	else {
	  double a = (w1[i+1] -w1[i])/(krel1[i+1] - krel1[i]);
	  double b = w1[i] - a*krel1[i];
	  fac1_2   =  a*Prel2  + b;}
	break;
      }else if(i==np1-1)fac1_2 = w1[i];
    } // end fac1_2
    
    
    for(int i=0;i<np2;i++){
      if(Prel2 < krel2[i]){
	if(i==0 || i==np2-1)fac2_2 = w2[i];
	else {
	  double a = (w2[i+1] -w2[i])/(krel2[i+1] - krel2[i]);
	  double b = w2[i] - a*krel2[i];
	  fac2_2 = a*Prel2 + b;
	  //	  cout<<" n "<<n<<" fac2_2 "<<fac2_2<<endl;
	}
	break;
      }else if(i==np2-1)fac2_2 = w2[i];
    } // end fac2_2
    



     
    for(int i=0;i<np3;i++){
      if(Prel2 < krel3[i]){
	if(i==0 || i==np3-1)fac3_2 = w3[i];
	else {
	  double a = (w3[i+1] -w3[i])/(krel3[i+1] - krel3[i]);
	  double b = w3[i] - a*krel3[i];
	  fac3_2 = a*Prel2 + b;  }
	break;
      }else if(i==np3-1)fac3_2 = w3[i];
    } // end fac3
    
    

    // I(nn) = I(n)*2


    if(E09){
    fac1 = (fac1-1.0)     +1.0;
    fac2 = (fac2-1.0)     +1.0;
    fac3 = (fac3-1.0)     +1.0;
    fac1_2 = (fac1_2-1.0) +1.0;
    fac2_2 = (fac2_2-1.0) +1.0;
    fac3_2 = (fac3_2-1.0) +1.0;    
    }else{
    fac1 = (fac1-1.0)*2.0 +1.0;
    fac2 = (fac2-1.0)*2.0 +1.0;
    fac3 = (fac3-1.0)*2.0 +1.0;
    fac1_2 = (fac1_2-1.0)*2.0 +1.0;
    fac2_2 = (fac2_2-1.0)*2.0 +1.0;
    fac3_2 = (fac3_2-1.0)*2.0 +1.0;    
    }
    



    
    if(E09)   mm =(mmnuc - MH3L )*1000.; // -BL [MeV]
    else      mm =(mmnuc - MnnL )*1000.; // -BL [MeV]
    hmm -> Fill(mm);
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
  


  

  hmm->Write();
  hmm_L->Write();
  hmm_nnL->Write();
  hmm_fsi1->Write();
  hmm_fsi2->Write();
  hmm_fsi3->Write();

  hmm_fsi1_2->Write();
  hmm_fsi2_2->Write();
  hmm_fsi3_2->Write();  

  
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
  pn = Pn.P();

  
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


    
    double a0,r0;
         if(model==0){  a0 = -2.68; r0 = 2.91;}
    else if(model==1 ){ a0 = -2.29; r0 = 3.15;}
    else if(model==2 ){ a0 = -1.60; r0 = 1.33;}
    else if(model==3 ){ a0 = -0.57; r0 = 7.65;}
    else if(model==-1){ a0 = -1.77; r0 = 3.25;}
    else if(model==-2){ a0 = -1.6;  r0 = 3.15;}
    else if(model==-3){ a0 = -1.94; r0 = 2.42;}

  // Kinematics
 
    if(nr_mode){
      skz = skz/hbarc;
      skz2 = skz*skz;
      eon = hbarc2*skz2/2.0/fmu;
      rho = fmu*skz/hbarc2;
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
	  if(kp==kpp) delta = complex<double>(1.0,0.0);
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
      

      // A : gren[i][j] = E -V * Go (C.11)
      // B : V
      // R : =(E - V * Go)^-1 * V (A) = A^-1 * B -> set B 


      ron  = complex<double>(0.0,0.0);
      ron0 = complex<double>(0.0,0.0);
      LEQ3(A,v,Rl,nrp1);
      LEQ3(A0,v0,Rl0,nrp1);


      ron       = Rl[nrp1];
      ron0      = Rl0[nrp1];
      ronr      = real(ron);
      delrad0   = atan(skz/(-1./a0 + r0*skz*skz/2.0));
      delrad    = -atan(PI*rho*ronr);
      deltal[l] = delrad;
      ton[l]    = ron/(1.0 +xi*PI*rho*ron);
      ton0[l]   = ron0/(1.0 +xi*PI*rho*ron0);
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

      //      cout<<"perre "<<1./perre<<" perim "<<1./perim<<endl;
      //      cout<<"ton "<<ton[l]<<" xi "<<xi<<" ron "<<ron<<" rho "<<rho<<" A "<<(1.0 +xi*PI*rho*ron)<<endl;
      
      
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

      for(int l =0;l <= lpjmax;l++){ // DO 90
	double ll=(double)l;
	fthe[ithe]   =  fthe[ithe]   + (2.0*ll +1.0)*ton[l]*pol[l]/4.0/PI;
	fthe0[ithe]  =  fthe0[ithe]  + (2.0*ll +1.0)*ton0[l]*pol[l]/4.0/PI;
	fthe1[ithe] =  fthe1[ithe] + (2.0*ll +1.0)*(ton[l] -tborn[l])*pol[l]/4.0/PI;

      } // END 90


      fthe[ithe]   = - fthe[ithe]*4.0*PI*PI*rho/skz;
      fthe0[ithe]  = - fthe0[ithe]*4.0*PI*PI*rho/skz;
      fthe1[ithe]  = -(fthe1[ithe]+vq)*4.0*PI*PI*rho/skz;
      fborn        += -vq*4.0*PI*PI*rho/skz;
      //     fborn0      = 4.0*PI*PI*rho/skz;
      totalc = totalc + fabs(fthe[ithe])*fabs(fthe[ithe])*sin(ther)*PI/180.*2*PI;

      cross[ithe] = real((fthe[ithe])* conj(fthe[ithe]) )+1.0;


      
    } // END 80


    //    cout<<"fthe "<<fthe[0]<<" fthe0 "<<fthe0[0]<<" ratio "<<pow((fthe[0]+fthe0[0])/fthe0[0],2.0) <<endl;
    
    tcross   = 4.0 * PI/skz *  imag(fthe[0]);
    tcross0  = 4.0 * PI/skz *  imag(fthe0[0]);
    tcross1  = 4.0 * PI/skz *  imag(fthe1[0]);
    tcrossl  = 4.0 * PI/skz *  imag(ftheL);

    
    double alpha =  sqrt(pow(1./r0,2.0) - 2./(a0*r0)) ;
    double beta  =  -2.0/r0 + alpha;
    complex<double> Jl = (skz - xi*beta)/(skz - xi*alpha);

    



    //    I0  = pow( sin(skz*r0 + deltal[LL])/sin(skz*r0) ,2.0);

    I0  = 1./( real(Jl*conj(Jl)) )/4.0/PI +1.0;
    I1  = tcrossl/4.0/PI +1.0;

    if(model==0)I1 = pow( sin(skz*r0 + deltal[LL])/sin(skz*r0) ,2.0) +1.0;
    
    double I = real(ftheL* conj(ftheL) )+1.0;

    //    cout<<"Jl "<<Jl<<" Jl* "<<conj(Jl)<<" I0 "<<Jl*conj(Jl)<<" Jl^2 "<<Jl*Jl<<endl;     
    
    for(int l=0;l<=LL;l++){
      Il[l]  = real(fthel[l]* conj(fthel[l])) + 1.0;
      Il0[l] = imag(fthel[l])/skz+ 1.0;
    }
    
    if( I0 >100.  ) I0 = 0.0;
    if( I1 >100.  ) I1 = 0.0;



    
    if(tcross0>0 || I <100.) return I;
    else     return 1.0;


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

  //test
  //  vr =0.0;
  //  b_r=1.0;
  
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
  

 //test
 // va   = 0.0; // MeV
 // b_a  = 0.0;   // fm

  
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
  string ifname = "../rootfiles/simc/3H_test.root";
  string ofname = "./test.root";
  string ofPname ="./param/test";
  string pname ="./param/potential.list";
  bool write_mode = false;
  bool root_mode  = false;
  bool param_mode = false;
  int l=0;
  string L;
  extern char *optarg;

  while((ch=getopt(argc,argv,"h:p:s:r:w:f:s:l:Ibcop"))!=-1){
    switch(ch){
            
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      root_mode =true;
      break;

    case 's':
      ifname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
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
  if(root_mode)  FSI -> SetBranch(ifname);
  FSI -> NewRoot(ofname);
  if(param_mode) FSI -> SetParam(pname);
  FSI -> SetHist();
  FSI -> SetFermiMom(fermi_mom);
  if(write_mode) FSI -> InfluenceFactor(ofPname,l);
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
