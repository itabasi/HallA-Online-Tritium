double corr_R_adc(int i);
double corr_R_x(int i);
double corr_R_th(int i);
double corr_R_alig(int i);
double corr_L_adc(int i);
double corr_L_x(int i);
double corr_L_th(int i);
double corr_L_alig(int i);
double s2f1_off(int i,char* ARM,char* MODE,int KINE);
const double c=299792458e-9;// [m/ns]
const double mk=493.7e-3;// Kaon mass [GeV/c^2]
const double me=0.511e-3;// electron mass [GeV/c^2] 
const double ml=1115.7e-3;//Lambda mass [GeV/c^2]
const double mn=939.6e-3; // neutron mass [GeV/c^2]
const double mpi=139.6e-3;// pion mass [GeV/c^2]
const  int nth=3; //th num
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
using namespace std;
#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "Setting.h"

//===================================================================//
//============================= Main ================================//
//===================================================================//

int main(int argc, char** argv){

  int ch; char* mode;
  int kine=1;// 1: hydrogen kinematics 2:tritium kinematics
  double tdc_time=58.0e-3;//[ns]
  string ifname ="/home/itabashi/jlab_nnL/ita_scripts/rootfiles/Lambda_small1_mm_ana1230.root";
  string ofname = "/pdf/hydro1_AC_eff_test.pdf";
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  bool test_flag=false;
  bool Pion_flag=false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"h:f:w:n:bcop:GHTt12P"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 'w':
      output_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
  
    case 'G':
    mode="G";
      break;
  
    case 'H':
    mode="H";
      break;

    case 'T':
      mode="T";    
	break;

    case 't':
      test_flag=true;    
	break;


    case 'P':
      Pion_flag=true;    
	break;

  case '1':
    tdc_time=56.23e-3;//[ns]
    kine=1;
      break;

  case '2':
    tdc_time=58e-3;//[ns]
    kine=2;
      break;


    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
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


  TApplication *theApp =new TApplication("App",&argc,argv);
 if(draw_flag==0)gROOT->SetBatch(1);


 Setting *set = new Setting();
 set->Initialize();




 
  cout<<"mode :"<<mode<<endl;
  cout<<"tdc_time[ns]: "<<tdc_time<<endl;

 
 

 //============= Set Branch Status ==================//

  int max=1000; 
  double RF1[max],LF1[max];
  double Rs0r_ac[max],Rs0l_ac[max],Ls0r_ac[max],Ls0l_ac[max];
  double Rs2r_ac[max],Rs2l_ac[max],Ls2r_ac[max],Ls2l_ac[max];
  double Rs0r_tc[max],Rs0l_tc[max],Ls0r_tc[max],Ls0l_tc[max];
  double Rs2r_tc[max],Rs2l_tc[max],Ls2r_tc[max],Ls2l_tc[max];
  double Ra1t[max],Ra1a[max],Ra1a_p[max],Ra1a_c[max],Ra1sum;
  double Ra2t[max],Ra2a[max],Ra2a_p[max],Ra2a_c[max],Ra2sum;
  double La1t[max],La1a[max],La1a_p[max],La1a_c[max],La1sum;
  double La2t[max],La2a[max],La2a_p[max],La2a_c[max],La2sum;
  double Rp[max],Rpx[max],Rpy[max],Lp[max],Lpx[max],Lpy[max];
  double Rth[max],Rph[max],Rx[max],Rvz[max],Lth[max],Lph[max],Lx[10],Ly[10],Lvz[100];
  double Rbeta[max],Lbeta[max];
  double rs2pathl[max],rs0pathl[max],rtrpathl[max];
  double ls2pathl[max],ls0pathl[max],ltrpathl[max];
  double trigger[100];
  double hallap;
  double Rs2tpads[100],Ls2tpads[100];
  double Rs2trpad[100],Ls2trpad[100];
  //---- Gogami root ---------//
  double ctime[1000];
  double DRT5;
  double Rs2ra[100];
  double Rs2la[100];
  double Ls2ra[100];
  double Ls2la[100];



  //=== z analysis parameters ====//
 double min_al=0.1;
 double max_al=0.2;
 double p0[1000];
 double p1[1000];
 double p2[1000];
 double p0_err[1000];
 double p1_err[1000];
 double p2_err[1000];
 double min_z,max_z;
 min_z=-0.2;
 max_z=0.2;
 int bin_z=500;
 double Rz;
 TGraphErrors* gAl_sig=new TGraphErrors();
 set->SetGr(gAl_sig,"Al peak Sigma Run dependence","run num","width [m]");
 TGraphErrors* gAl_mean=new TGraphErrors();
 set->SetGr(gAl_mean,"Al peak mean Run dependence","run num","mean [m]");
 TH1F* hz[1000];
 double top_x[1000],top_y[1000];
 char* r="r";
 
 int l=0; // run counts

 int x; // run number
 int y; // sub run number
 int z; // run number small 3 order
 
 //=== Tohoku nnL =====//
 x=28; 
 y=35;
 z=31;

 //===== ifarm ========//
 x=64;
 y=71;
 z=67;

  //===========================================================//
  //============== RUN STUDY ==================================//
  //===========================================================//

 int runnum,sub_num,i;
  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname, run_num,sub,run_i;

  cout<<"start analysis !"<<endl;
  while(1){

     //TChain //
  TChain* T=new TChain("T");

  
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> runname;
    run_num=buf.substr(x,6);
    runnum=atoi(run_num.c_str());
    sub=buf.substr(y,1);
    runnum=atoi(run_num.c_str());
    run_i=buf.substr(z,3);
    i=atoi(run_i.c_str());
    // cout<<"run num: "<<runnum<<endl;
    //    cout<<"sub num: "<<sub<<endl;

    if(sub!=r){T->Add(runname.c_str());}
    else if(sub==r){
           T->Add(runname.c_str());
	   cout<<"run num: "<<runnum<<endl;
    int evnt=T->GetEntries();
    if(test_flag){evnt=100000;}
    cout<<"Get Entries: "<<evnt<<endl;

  T->SetBranchStatus("*",0);  
 //------ Right Arm -------------//

  // T->SetBranchStatus("RTDC.F1FirstHit",1);
  // T->SetBranchAddress("RTDC.F1FirstHit",RF1); 
  // T->SetBranchStatus("R.s2.t_pads",1);
  // T->SetBranchAddress("R.s2.t_pads",Rs2tpads);
  // T->SetBranchStatus("R.s2.trpad",1);
  //T->SetBranchAddress("R.s2.trpad",Rs2trpad);
  // T->SetBranchStatus("R.s2.ra",1);
  //T->SetBranchAddress("R.s2.ra",Rs2ra);
  // T->SetBranchStatus("R.s2.la",1);
  // T->SetBranchAddress("R.s2.la",Rs2la);
  // T->SetBranchStatus("R.a1.asum_c",1);
  // T->SetBranchAddress("R.a1.asum_c",&Ra1sum);
  // T->SetBranchStatus("R.a2.asum_c",1);
  //T->SetBranchAddress("R.a2.asum_c",&Ra2sum);

  // path length//
  //T->SetBranchStatus("R.s2.trpath",1); 
  // T->SetBranchAddress("R.s2.trpath",rs2pathl); 
  //T->SetBranchStatus("R.tr.pathl",1);  
  //T->SetBranchAddress("R.tr.pathl",rtrpathl);
 // Target positon information //
  // T->SetBranchStatus("R.tr.p",1);
  // T->SetBranchAddress("R.tr.p",Rp);
  //T->SetBranchStatus("R.tr.vz",1);    
 T->SetBranchAddress("R.tr.vz",Rvz); 

 //------ Left Arm ---------------//
 /*
 T->SetBranchStatus("LTDC.F1FirstHit",1);
 T->SetBranchAddress("LTDC.F1FirstHit",LF1); 
 T->SetBranchStatus("L.s2.t_pads",1);
 T->SetBranchAddress("L.s2.t_pads",Ls2tpads);
 T->SetBranchStatus("L.s2.trpad",1);
 T->SetBranchAddress("L.s2.trpad",Ls2trpad);
 T->SetBranchStatus("L.s2.ra",1);
 T->SetBranchAddress("L.s2.ra",Ls2ra);
 T->SetBranchStatus("L.s2.la",1);
 T->SetBranchAddress("L.s2.la",Ls2la);
  // path length//
 T->SetBranchStatus("L.s2.trpath",1); 
 T->SetBranchAddress("L.s2.trpath",ls2pathl); 
 T->SetBranchStatus("L.tr.pathl",1);   
 T->SetBranchAddress("L.tr.pathl",ltrpathl);
 T->SetBranchStatus("L.tr.p",1);
 T->SetBranchAddress("L.tr.p",Lp);  
 T->SetBranchStatus("L.tr.vz",1);    
 T->SetBranchAddress("L.tr.vz",Lvz);
 T->SetBranchStatus("L.tr.x",1);    
 T->SetBranchAddress("L.tr.x",Lx);
 T->SetBranchStatus("L.tr.y",1);    
 T->SetBranchAddress("L.tr.y",Ly);
 T->SetBranchStatus("L.tr.ph",1);    
 T->SetBranchAddress("L.tr.ph",Lph);
 T->SetBranchStatus("L.tr.th",1);    
 T->SetBranchAddress("L.tr.th",Lth);

 */


 //========= Definition of Hist ========================// 
 

 hz[i]=new TH1F(Form("hz[%d]",i),"",bin_z,min_z,max_z);
 set->SetTH1(hz[i],"z hist","z [m]","Counts");



 
 for(int k=0;k<evnt;k++){
   T->GetEntry(k);

   Rz=-100.;
   Rz=Rvz[0];
   hz[i]->Fill(Rz);

 

 }//end fill

 
 //====== Fitting =========//


 TF1* fAl=new TF1("fAl","gausn(0)",min_z,max_z);
 hz[i]->Fit("fAl","Rq","",min_al,max_al);
 p0[i]=fAl->GetParameter(0);
 p1[i]=fAl->GetParameter(1);
 p2[i]=fAl->GetParameter(2); 
 p0_err[i]=fAl->GetParError(0);
 p1_err[i]=fAl->GetParError(1);
 p2_err[i]=fAl->GetParError(2); 

 gAl_sig->SetPoint(l,i,p2[i]);
 gAl_mean->SetPoint(l,i,p1[i]);
 gAl_sig->SetPointError(l,0,p2_err[i]);
 gAl_mean->SetPointError(l,0,p1_err[i]);



 l=l+1; 
    }//end TChain

    delete T;
  }// end while

  TCanvas* c0=new TCanvas("c0","c0");
  TCanvas* c1=new TCanvas("c1","c1");

  c0->cd();
  gAl_mean->GetYaxis()->SetRangeUser(0.13,0.14);
  gAl_mean->Draw("AP");
  c1->cd();
  gAl_sig->GetYaxis()->SetRangeUser(0.01,0.02);
  gAl_sig->Draw("AP");



 TString name;
 if(output_flag){
 name.Form(ofname.c_str());
 c0->Print(name+"[","pdf");
 c0->Print(name,"pdf");
 c1->Print(name,"pdf");
 c1->Print(name +"]","pdf");
}


    
 if(draw_flag==0)gSystem->Exit(1);
 theApp->Run();
 return 0;



}//end Main


//==============================================//
//========== Defined Function ==================//
//=============================================//

double s2f1_off(int i,char* ARM,char* MODE, int KINE){


  double RS2_offset[16],LS2_offset[16];
  if(MODE=="H" && KINE==2){
 
 double  RS2_off_H2[16]={-16911.4,-16864.3,-16900,-16897,-16873.8,-16868.4,-16901.1,-16876.8,-16895.4,-16860.9,-16893.1,-16884.4,-16847.3,-16842.7,-16836.9,-16882.6};
 double  LS2_off_H2[16]={-25336.9,-25386.6,-25367.5,-25392.3,-25391.1,-25386.2,-25422,-25428.9,-25417.3,-25426.8,-25438.7,-25383.4,-25396,-25418.5,-25436.4,-26082.1};
 
  LS2_offset[i]=LS2_off_H2[i];
  RS2_offset[i]=RS2_off_H2[i];
  }


  if(MODE=="H" && KINE==1){
    
    //double  RS2_off_H1[16]={-16911.4,-16864.9,-16900,-16897.6,-16874.8,-16869.3,-16901.1,-16876.8,-16895.6,-16860.3,-16892.6,-16885,-16847.3,-16843.3,-16838.4,-16882.6};
    //double  LS2_off_H1[16]={-25336.9,-25385.7,-25367,-25392.2,-25391,-25386.3,-25422,-25428.9,-25415.2,-25425,-25438,-25381,-25394.4,-25417.5,-25432.8,-26082.1};

double  RS2_off_H1[16]={-16828.7,-16863,-16894,-16893.3,-16870.9,-16867.2,-16900.3,-16876.8,-16876.8,-16861.6,-16895,-16890.7,-16854.6,-16852.9,-16850.5,-16861.9};
double  LS2_off_H1[16]={-25335,-25385.6,-25367,-25392.1,-25391.7,-25386.4,-25422.1,-25428.9,-25414.9,-25424.7,-25436.9, -25381.2,-25390,-25413.4,-25428.7,-26640.8};
  LS2_offset[i]=LS2_off_H1[i];
  RS2_offset[i]=RS2_off_H1[i];
  }

 double s2f1_offset; 
 if(ARM=="R")s2f1_offset=RS2_offset[i];
 else  if(ARM=="L")s2f1_offset=LS2_offset[i];
 else {cout<<"false read out !!"<<endl;}

  return s2f1_offset;

}
