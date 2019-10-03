#ifndef VDCt0_h
#define VDCt0_h 1
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "TApplication.h"
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
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
#include "Setting.h"

const   int nmax=400;
const   int nwire=368;

bool tuned_flag;


class VDCt0{

 public:
  VDCt0();
  ~VDCt0();

  
  
 public:
  void SetRunList(string ifname);
  void SetRun(int runnum);
  void SetBranch();
  void NewRoot(string ofname);
  void MakeHist();
  void Fill();
  void GetOffset(string paraname);
  double Findt0(bool rarm,char* plane, int wire);
  double Findt0_time(bool rarm, char* plane, int wire);
  void Deft0(string ifname);
  void SetRawt0();
  void Sett0();
  void Write(string ofname);
  void MakeRoot(string ofname);
  void MakeRoot_c(string ofname);  
  void Draw();
  void Print(string ofname);
  void Print_c(string ofname);
  
  Setting* set;
  //== SetRunList ====//
    int ENum;
    TChain* T;
    double off[nwire][8];
    
  //=== SetBranch ====//
  double evtype, HallA_p;  
  double Ru1_nhit[nmax],Ru2_nhit[nmax],Rv1_nhit[nmax],Rv2_nhit[nmax];
  int NRu1_nhit,NRu2_nhit,NRv1_nhit,NRv2_nhit;
  double Ru1_wire[nmax],Ru2_wire[nmax],Rv1_wire[nmax],Rv2_wire[nmax];
  int NRu1_wire,NRu2_wire,NRv1_wire,NRv2_wire;
  double Ru1_rtime[nmax],Ru2_rtime[nmax],Rv1_rtime[nmax],Rv2_rtime[nmax];
  int NRu1_rtime,NRu2_rtime,NRv1_rtime,NRv2_rtime;
  double Ru1_time[nmax],Ru2_time[nmax],Rv1_time[nmax],Rv2_time[nmax];
  int NRu1_time,NRu2_time,NRv1_time,NRv2_time;

  double Lu1_nhit[nmax],Lu2_nhit[nmax],Lv1_nhit[nmax],Lv2_nhit[nmax];
  int NLu1_nhit,NLu2_nhit,NLv1_nhit,NLv2_nhit;
  double Lu1_wire[nmax],Lu2_wire[nmax],Lv1_wire[nmax],Lv2_wire[nmax];
  int NLu1_wire,NLu2_wire,NLv1_wire,NLv2_wire;
  double Lu1_rtime[nmax],Lu2_rtime[nmax],Lv1_rtime[nmax],Lv2_rtime[nmax];
  int NLu1_rtime,NLu2_rtime,NLv1_rtime,NLv2_rtime;
 double Lu1_time[nmax],Lu2_time[nmax],Lv1_time[nmax],Lv2_time[nmax];
  int NLu1_time,NLu2_time,NLv1_time,NLv2_time;      

  //===== NewRoot =====//
  TFile* fnew;
  TTree* tnew;
  double Ru1_rt_p[500],Ru2_rt_p[500],Rv1_rt_p[500],Rv2_rt_p[500];
  double Lu1_rt_p[500],Lu2_rt_p[500],Lv1_rt_p[500],Lv2_rt_p[500];


  //===== MakeHist =====//
  TH2F* hRu1;
  TH2F* hRu2;
  TH2F* hRv1;
  TH2F* hRv2;
  TH2F* hLu1;
  TH2F* hLu2;
  TH2F* hLv1;
  TH2F* hLv2;

  TH2F* hRu1_c;
  TH2F* hRu2_c;
  TH2F* hRv1_c;
  TH2F* hRv2_c;
  TH2F* hLu1_c;
  TH2F* hLu2_c;
  TH2F* hLv1_c;
  TH2F* hLv2_c;  
  
  TGraphErrors* gRu1;
  TGraphErrors* gRu2;
  TGraphErrors* gRv1;
  TGraphErrors* gRv2;
  TGraphErrors* gLu1;
  TGraphErrors* gLu2;
  TGraphErrors* gLv1;
  TGraphErrors* gLv2;

  TGraphErrors* gRu1_c;
  TGraphErrors* gRu2_c;
  TGraphErrors* gRv1_c;
  TGraphErrors* gRv2_c;
  TGraphErrors* gLu1_c;
  TGraphErrors* gLu2_c;
  TGraphErrors* gLv1_c;
  TGraphErrors* gLv2_c;  
  
  TH1D* hRu1_rtime[nwire];
  TH1D* hRu2_rtime[nwire];
  TH1D* hRv1_rtime[nwire];
  TH1D* hRv2_rtime[nwire];
  TH1D* hLu1_rtime[nwire];
  TH1D* hLu2_rtime[nwire];
  TH1D* hLv1_rtime[nwire];
  TH1D* hLv2_rtime[nwire];

  TH1D* hRu1_time[nwire];
  TH1D* hRu2_time[nwire];
  TH1D* hRv1_time[nwire];
  TH1D* hRv2_time[nwire];
  TH1D* hLu1_time[nwire];
  TH1D* hLu2_time[nwire];
  TH1D* hLv1_time[nwire];
  TH1D* hLv2_time[nwire];    


  TF1* fLu1_t0[nwire];
  TF1* fLu2_t0[nwire];
  TF1* fLv1_t0[nwire];
  TF1* fLv2_t0[nwire];  
  TF1* fRu1_t0[nwire];
  TF1* fRu2_t0[nwire];
  TF1* fRv1_t0[nwire];
  TF1* fRv2_t0[nwire];

  TF1* fLu1_rt0[nwire];
  TF1* fLu2_rt0[nwire];
  TF1* fLv1_rt0[nwire];
  TF1* fLv2_rt0[nwire];  
  TF1* fRu1_rt0[nwire];
  TF1* fRu2_rt0[nwire];
  TF1* fRv1_rt0[nwire];
  TF1* fRv2_rt0[nwire];    

  
  double min_rtime,max_rtime;
  int bin_rtime;
  double min_time,max_time;
  int bin_time;

  //==== FIndT0 ======//
  //    double dy,y,yb,ya;
    double a_min,a_max,slope_1,slope_2;
    double t0_min,t0_max,t0_1,t0_2;  
    double dt0;
    double t0offset;    
  //==== Draw =====//

  TCanvas* c0[11];// = new TCanvas("c0","c0");
  //c0->Divide(4,8);
  TCanvas* c1[11];// = new TCanvas("c1","c1");
  //  c1->Divide(4,8);
  TCanvas* c2[11];// = new TCanvas("c2","c2");
  //  c2->Divide(4,8);  
  TCanvas* c3[11];// = new TCanvas("c3","c3");
  //  c3->Divide(4,8);
   TCanvas* c4[11];// = new TCanvas("c0","c0");
  //c0->Divide(4,8);
  TCanvas* c5[11];// = new TCanvas("c1","c1");
  //  c1->Divide(4,8);
  TCanvas* c6[11];// = new TCanvas("c2","c2");
  //  c2->Divide(4,8);  
  TCanvas* c7[11];// = new TCanvas("c3","c3");
  //  c3->Divide(4,8);
  TCanvas* c10;
  /*
  TCanvas* c4 = new TCanvas("c4","c4");
  c4->Divide(4,8);  
  TCanvas* c5 = new TCanvas("c5","c5");
  c5->Divide(4,8);  
  TCanvas* c6 = new TCanvas("c6","c6");
  c6->Divide(4,8);  
  TCanvas* c7 = new TCanvas("c7","c7");
  c7->Divide(4,8);  
  TCanvas* c8 = new TCanvas("c8","c8");
  c8->Divide(4,8);  
  TCanvas* c9 = new TCanvas("c9","c9");
  c9->Divide(4,8);
  TCanvas* c10= new TCanvas("c10","c10");
  c10->Divide(4,8);  
  TCanvas* c11= new TCanvas("c11","c11");
  c11->Divide(4,8);
  */

  TLine* line_time[nwire];
  TLine* line_rtime[nwire];   
  //  TCanvas* c12= new TCanvas("c12","c12");
  //  TCanvas* c13= new TCanvas("c13","c13");
  //  TCanvas* c14= new TCanvas("c14","c14");
  //  TCanvas* c15= new TCanvas("c15","c15");  


  
  //=== Def param ====//
  double  Ru1t0_def[nwire],Ru2t0_def[nwire],Rv1t0_def[nwire],Rv2t0_def[nwire];
  double  Lu1t0_def[nwire],Lu2t0_def[nwire],Lv1t0_def[nwire],Lv2t0_def[nwire];  
  double T0_def[nwire][4];
  //==== Write ====//
  double Ru1t0[nwire],Ru2t0[nwire],Rv1t0[nwire],Rv2t0[nwire];
  double Lu1t0[nwire],Lu2t0[nwire],Lv1t0[nwire],Lv2t0[nwire];
  double Ru1t0_err[nwire],Ru2t0_err[nwire],Rv1t0_err[nwire],Rv2t0_err[nwire];
  double Lu1t0_err[nwire],Lu2t0_err[nwire],Lv1t0_err[nwire],Lv2t0_err[nwire];
  double Ru1t0_c[nwire],Ru2t0_c[nwire],Rv1t0_c[nwire],Rv2t0_c[nwire];
  double Lu1t0_c[nwire],Lu2t0_c[nwire],Lv1t0_c[nwire],Lv2t0_c[nwire];
  double Ru1t0_err_c[nwire],Ru2t0_err_c[nwire],Rv1t0_err_c[nwire],Rv2t0_err_c[nwire];
  double Lu1t0_err_c[nwire],Lu2t0_err_c[nwire],Lv1t0_err_c[nwire],Lv2t0_err_c[nwire];

  
  double Ru1t0_min[nwire],Ru1t0_max[nwire],Ru2t0_min[nwire],Ru2t0_max[nwire],Rv1t0_min[nwire],Rv1t0_max[nwire],Rv2t0_min[nwire],Rv2t0_max[nwire];
  double Lu1t0_min[nwire],Lu1t0_max[nwire],Lu2t0_min[nwire],Lu2t0_max[nwire],Lv1t0_min[nwire],Lv1t0_max[nwire],Lv2t0_min[nwire],Lv2t0_max[nwire];
  

  //==== MakeRoot =====//

  //==== Findt0 =====//
  double  Lu1_p0[nwire],Lu1_p1[nwire],Lu2_p0[nwire],Lu2_p1[nwire],Lv1_p0[nwire],Lv1_p1[nwire],Lv2_p0[nwire],Lv2_p1[nwire];
 double  Ru1_p0[nwire],Ru1_p1[nwire],Ru2_p0[nwire],Ru2_p1[nwire],Rv1_p0[nwire],Rv1_p1[nwire],Rv2_p0[nwire],Rv2_p1[nwire];

  
};


VDCt0::VDCt0(){
  cout<<"start VDC T0 tuning "<<endl;
  set=new Setting();
  set->Initialize();
  for(int i=0;i<nwire;i++)
    for(int j=0;j<8;j++)off[i][j]=0.0;
  
}
VDCt0::~VDCt0(){}


/////////////////////////////////////////////////////////////////////////////////



void VDCt0::SetRunList(string ifname){

  T=new TChain("T");
  
    ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname,runnum,runnum_group;
  int num;
  int run_group;
  int k=0;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;

    for(int i=0;i<(int)runname.length();i++){
      string   st=runname.substr(i,1);
      if(st=="1"){
	   string run_test0=runname.substr(i,6);
	  int run_test=atoi(run_test0.c_str());
	  if(run_test>111111)num=run_test;
      }
    }


    
    
    //    runnum = runname.substr(35,6);
    //    num =atoi(runnum.c_str());
    //    runnum_group = runname.substr(35,5);    
    //    run_group=atoi(runnum_group.c_str());
    
    T->Add(runname.c_str());
  }
  ENum=T->GetEntries();
  cout<<"Events: "<<ENum<<endl; 
  
  
}

//////////////////////////////////////////////////////////

void VDCt0::SetRun(int runnum){
  T=new TChain("T");
  //####################
  int sum_run=3;
  //#####################
  //  cout<<"TChain run number : "<<runnum<<" - "<<runnum+sum_run-1<<endl;

  //  const string ROOTfilePath="/data/opt_small/VDC/initial";
  const string ROOTfilePath="/data2/opt_small/VDC/initial/";  
  const string root =".root";
  ostringstream str;
  int run;
  for(int i=0;i<sum_run;i++){
    //==== Initialization of string =======//
    str.str("");
    
    str.clear(ostringstream::goodbit);
    //=====================================//
  run=runnum+i;
  str<< ROOTfilePath <<"tritium_" <<run;
  string basename = str.str().c_str();
  string rootfile = basename + root ;
   //======== SUB RUN ===================//
   int sub=1;
   while ( !gSystem->AccessPathName(rootfile.c_str()) ) {
        T->Add(rootfile.c_str());
        cout << "ROOT file " << rootfile << " added to TChain." << endl;
        rootfile = basename + "_" + sub + ".root";
	sub++;
   }
  }//end for

  ENum=T->GetEntries();
  cout<<"Events: "<<ENum<<endl; 

  
}

/////////////////////////////////////////////////////

void VDCt0::GetOffset(string paraname){

  ifstream ifparam(paraname.c_str());  
  string buf;
  if (ifparam.fail()){ cerr << "failed open files" <<paraname<<endl; exit(1);}

  int plane=-1;
  int i=0;

  while( getline(ifparam,buf) ){

    if( buf[0]=='#' ){
      i=0;
      plane++;}
    if( ifparam.eof() || i+1 > nwire){continue;}

    ifparam >> off[i][plane] >> off[i+1][plane] >> off[i+2][plane] >> off[i+3][plane]
	    >> off[i+4][plane] >> off[i+5][plane] >> off[i+6][plane] >> off[i+7][plane];
    i=i+8;
  }//end while dat file

  
}


/////////////////////////////////////////////////////

void VDCt0::SetBranch(){

  T->SetBranchStatus("*",0);
  //  T->SetBranchStatus("fEvtHdr.fRun"   ,1);
  //  T->SetBranchStatus("fEvtHdr.fEvtNum" ,1);
  T->SetBranchStatus("HALLA_p"        ,1);          T->SetBranchAddress("HALLA_p"        ,&HallA_p); 
  T->SetBranchStatus("DR.evtypebits"        ,1);    T->SetBranchAddress("DR.evtypebits"        ,&evtype);

  
  //==========================//
  //=========== VDC ==========//
  //==========================//

  
  //========== RHRS VDC ==========//
  
   T->SetBranchStatus("R.vdc.u1.nhit"          ,1);            T->SetBranchAddress("R.vdc.u1.nhit"          ,Ru1_nhit);
   T->SetBranchStatus("R.vdc.u2.nhit"          ,1);            T->SetBranchAddress("R.vdc.u2.nhit"          ,Ru2_nhit);
   T->SetBranchStatus("R.vdc.v1.nhit"          ,1);            T->SetBranchAddress("R.vdc.v1.nhit"          ,Rv1_nhit);
   T->SetBranchStatus("R.vdc.v2.nhit"          ,1);            T->SetBranchAddress("R.vdc.v2.nhit"          ,Rv2_nhit);
  //  T->SetBranchStatus("Ndata.R.vdc.u1.nhit"          ,1);      T->SetBranchAddress("Ndata.R.vdc.u1.nhit"          ,&NRu1_nhit);
  //  T->SetBranchStatus("Ndata.R.vdc.u2.nhit"          ,1);      T->SetBranchAddress("Ndata.R.vdc.u2.nhit"          ,&NRu2_nhit);
  //  T->SetBranchStatus("Ndata.R.vdc.v1.nhit"          ,1);      T->SetBranchAddress("Ndata.R.vdc.v1.nhit"          ,&NRv1_nhit);
  //  T->SetBranchStatus("Ndata.R.vdc.v2.nhit"          ,1);      T->SetBranchAddress("Ndata.R.vdc.v2.nhit"          ,&NRv2_nhit);
  T->SetBranchStatus("R.vdc.u1.wire"          ,1);            T->SetBranchAddress("R.vdc.u1.wire"          ,Ru1_wire);
  T->SetBranchStatus("R.vdc.u2.wire"          ,1);            T->SetBranchAddress("R.vdc.u2.wire"          ,Ru2_wire);
  T->SetBranchStatus("R.vdc.v1.wire"          ,1);            T->SetBranchAddress("R.vdc.v1.wire"          ,Rv1_wire);
  T->SetBranchStatus("R.vdc.v2.wire"          ,1);            T->SetBranchAddress("R.vdc.v2.wire"          ,Rv2_wire);
  T->SetBranchStatus("Ndata.R.vdc.u1.wire"          ,1);      T->SetBranchAddress("Ndata.R.vdc.u1.wire"          ,&NRu1_wire);
  T->SetBranchStatus("Ndata.R.vdc.u2.wire"          ,1);      T->SetBranchAddress("Ndata.R.vdc.u2.wire"          ,&NRu2_wire);
  T->SetBranchStatus("Ndata.R.vdc.v1.wire"          ,1);      T->SetBranchAddress("Ndata.R.vdc.v1.wire"          ,&NRv1_wire);
  T->SetBranchStatus("Ndata.R.vdc.v2.wire"          ,1);      T->SetBranchAddress("Ndata.R.vdc.v2.wire"          ,&NRv2_wire);
  T->SetBranchStatus("R.vdc.u1.rawtime"          ,1);         T->SetBranchAddress("R.vdc.u1.rawtime"          ,Ru1_rtime);         
  T->SetBranchStatus("R.vdc.u2.rawtime"          ,1);         T->SetBranchAddress("R.vdc.u2.rawtime"          ,Ru2_rtime);         
  T->SetBranchStatus("R.vdc.v1.rawtime"          ,1);         T->SetBranchAddress("R.vdc.v1.rawtime"          ,Rv1_rtime);         
  T->SetBranchStatus("R.vdc.v2.rawtime"          ,1);         T->SetBranchAddress("R.vdc.v2.rawtime"          ,Rv2_rtime);         
  T->SetBranchStatus("Ndata.R.vdc.u1.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.u1.rawtime"          ,&NRu1_rtime);         
  T->SetBranchStatus("Ndata.R.vdc.u2.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.u2.rawtime"          ,&NRu2_rtime);         
  T->SetBranchStatus("Ndata.R.vdc.v1.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.v1.rawtime"          ,&NRv1_rtime);         
  T->SetBranchStatus("Ndata.R.vdc.v2.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.v2.rawtime"          ,&NRv2_rtime);         

  T->SetBranchStatus("R.vdc.u1.time"          ,1);         T->SetBranchAddress("R.vdc.u1.time"          ,Ru1_time);         
  T->SetBranchStatus("R.vdc.u2.time"          ,1);         T->SetBranchAddress("R.vdc.u2.time"          ,Ru2_time);         
  T->SetBranchStatus("R.vdc.v1.time"          ,1);         T->SetBranchAddress("R.vdc.v1.time"          ,Rv1_time);         
  T->SetBranchStatus("R.vdc.v2.time"          ,1);         T->SetBranchAddress("R.vdc.v2.time"          ,Rv2_time);         


  /*
  T->SetBranchStatus("R.vdc.u1.time"          ,1);
  T->SetBranchStatus("R.vdc.u2.time"          ,1);
  T->SetBranchStatus("R.vdc.v1.time"          ,1);
  T->SetBranchStatus("R.vdc.v2.time"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.u1.time"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.u2.time"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.v1.time"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.v2.time"          ,1);  
  T->SetBranchStatus("R.vdc.u1.dist"          ,1);
  T->SetBranchStatus("R.vdc.u2.dist"          ,1);
  T->SetBranchStatus("R.vdc.v1.dist"          ,1);
  T->SetBranchStatus("R.vdc.v2.dist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.u1.dist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.u2.dist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.v1.dist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.v2.dist"          ,1);  
  T->SetBranchStatus("R.vdc.u1.ddist"          ,1);
  T->SetBranchStatus("R.vdc.u2.ddist"          ,1);
  T->SetBranchStatus("R.vdc.v1.ddist"          ,1);
  T->SetBranchStatus("R.vdc.v2.ddist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.u1.ddist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.u2.ddist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.v1.ddist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.v2.ddist"          ,1);    
  T->SetBranchStatus("R.vdc.u1.trdist"          ,1);
  T->SetBranchStatus("R.vdc.u2.trdist"          ,1);
  T->SetBranchStatus("R.vdc.v1.trdist"          ,1);
  T->SetBranchStatus("R.vdc.v2.trdist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.u1.trdist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.u2.trdist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.v1.trdist"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.v2.trdist"          ,1);  
  */
  
  //========== LHRS VDC ===============//

  T->SetBranchStatus("L.vdc.u1.nhit"          ,1);            T->SetBranchAddress("L.vdc.u1.nhit"          ,Lu1_nhit);
  T->SetBranchStatus("L.vdc.u2.nhit"          ,1);            T->SetBranchAddress("L.vdc.u2.nhit"          ,Lu2_nhit);
  T->SetBranchStatus("L.vdc.v1.nhit"          ,1);            T->SetBranchAddress("L.vdc.v1.nhit"          ,Lv1_nhit);
  T->SetBranchStatus("L.vdc.v2.nhit"          ,1);            T->SetBranchAddress("L.vdc.v2.nhit"          ,Lv2_nhit);
  //  T->SetBranchStatus("Ndata.L.vdc.u1.nhit"          ,1);      T->SetBranchAddress("Ndata.L.vdc.u1.nhit"          ,&NLu1_nhit);
  //  T->SetBranchStatus("Ndata.L.vdc.u2.nhit"          ,1);      T->SetBranchAddress("Ndata.L.vdc.u2.nhit"          ,&NLu2_nhit);
  //  T->SetBranchStatus("Ndata.L.vdc.v1.nhit"          ,1);      T->SetBranchAddress("Ndata.L.vdc.v1.nhit"          ,&NLv1_nhit);
  //  T->SetBranchStatus("Ndata.L.vdc.v2.nhit"          ,1);      T->SetBranchAddress("Ndata.L.vdc.v2.nhit"          ,&NLv2_nhit);
  T->SetBranchStatus("L.vdc.u1.wire"          ,1);            T->SetBranchAddress("L.vdc.u1.wire"          ,Lu1_wire);
  T->SetBranchStatus("L.vdc.u2.wire"          ,1);            T->SetBranchAddress("L.vdc.u2.wire"          ,Lu2_wire);
  T->SetBranchStatus("L.vdc.v1.wire"          ,1);            T->SetBranchAddress("L.vdc.v1.wire"          ,Lv1_wire);
  T->SetBranchStatus("L.vdc.v2.wire"          ,1);            T->SetBranchAddress("L.vdc.v2.wire"          ,Lv2_wire);
  T->SetBranchStatus("Ndata.L.vdc.u1.wire"          ,1);      T->SetBranchAddress("Ndata.L.vdc.u1.wire"          ,&NLu1_wire);
  T->SetBranchStatus("Ndata.L.vdc.u2.wire"          ,1);      T->SetBranchAddress("Ndata.L.vdc.u2.wire"          ,&NLu2_wire);
  T->SetBranchStatus("Ndata.L.vdc.v1.wire"          ,1);      T->SetBranchAddress("Ndata.L.vdc.v1.wire"          ,&NLv1_wire);
  T->SetBranchStatus("Ndata.L.vdc.v2.wire"          ,1);      T->SetBranchAddress("Ndata.L.vdc.v2.wire"          ,&NLv2_wire);
  T->SetBranchStatus("L.vdc.u1.rawtime"          ,1);         T->SetBranchAddress("L.vdc.u1.rawtime"          ,Lu1_rtime);         
  T->SetBranchStatus("L.vdc.u2.rawtime"          ,1);         T->SetBranchAddress("L.vdc.u2.rawtime"          ,Lu2_rtime);         
  T->SetBranchStatus("L.vdc.v1.rawtime"          ,1);         T->SetBranchAddress("L.vdc.v1.rawtime"          ,Lv1_rtime);         
  T->SetBranchStatus("L.vdc.v2.rawtime"          ,1);         T->SetBranchAddress("L.vdc.v2.rawtime"          ,Lv2_rtime);         
  T->SetBranchStatus("Ndata.L.vdc.u1.rawtime"          ,1);   T->SetBranchAddress("Ndata.L.vdc.u1.rawtime"          ,&NLu1_rtime);         
  T->SetBranchStatus("Ndata.L.vdc.u2.rawtime"          ,1);   T->SetBranchAddress("Ndata.L.vdc.u2.rawtime"          ,&NLu2_rtime);         
  T->SetBranchStatus("Ndata.L.vdc.v1.rawtime"          ,1);   T->SetBranchAddress("Ndata.L.vdc.v1.rawtime"          ,&NLv1_rtime);         
  T->SetBranchStatus("Ndata.L.vdc.v2.rawtime"          ,1);   T->SetBranchAddress("Ndata.L.vdc.v2.rawtime"          ,&NLv2_rtime);         
  T->SetBranchStatus("L.vdc.u1.time"          ,1);         T->SetBranchAddress("L.vdc.u1.time"          ,Lu1_time);         
  T->SetBranchStatus("L.vdc.u2.time"          ,1);         T->SetBranchAddress("L.vdc.u2.time"          ,Lu2_time);         
  T->SetBranchStatus("L.vdc.v1.time"          ,1);         T->SetBranchAddress("L.vdc.v1.time"          ,Lv1_time);         
  T->SetBranchStatus("L.vdc.v2.time"          ,1);         T->SetBranchAddress("L.vdc.v2.time"          ,Lv2_time);         


  

  /*
  T->SetBranchStatus("L.vdc.u1.time"          ,1);
  T->SetBranchStatus("L.vdc.u2.time"          ,1);
  T->SetBranchStatus("L.vdc.v1.time"          ,1);
  T->SetBranchStatus("L.vdc.v2.time"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.u1.time"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.u2.time"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.v1.time"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.v2.time"          ,1);  
  T->SetBranchStatus("L.vdc.u1.dist"          ,1);
  T->SetBranchStatus("L.vdc.u2.dist"          ,1);
  T->SetBranchStatus("L.vdc.v1.dist"          ,1);
  T->SetBranchStatus("L.vdc.v2.dist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.u1.dist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.u2.dist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.v1.dist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.v2.dist"          ,1);  
  T->SetBranchStatus("L.vdc.u1.ddist"          ,1);
  T->SetBranchStatus("L.vdc.u2.ddist"          ,1);
  T->SetBranchStatus("L.vdc.v1.ddist"          ,1);
  T->SetBranchStatus("L.vdc.v2.ddist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.u1.ddist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.u2.ddist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.v1.ddist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.v2.ddist"          ,1);    
  T->SetBranchStatus("L.vdc.u1.trdist"          ,1);
  T->SetBranchStatus("L.vdc.u2.trdist"          ,1);
  T->SetBranchStatus("L.vdc.v1.trdist"          ,1);
  T->SetBranchStatus("L.vdc.v2.trdist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.u1.trdist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.u2.trdist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.v1.trdist"          ,1);
  T->SetBranchStatus("Ndata.L.vdc.v2.trdist"          ,1);  
    
  */
  
}

///////////////////////////////////////////////////

void VDCt0::NewRoot(string ofname){

  
  fnew = new TFile(Form("%s",ofname.c_str()),"recreate");
  tnew =new TTree("T",ofname.c_str());

  tnew->Branch("Ru1_rt_p",Ru1_rt_p,"Ru1_rt_p[500]/D");
  tnew->Branch("Ru2_rt_p",Ru2_rt_p,"Ru2_rt_p[500]/D");
  tnew->Branch("Rv1_rt_p",Rv1_rt_p,"Rv1_rt_p[500]/D");
  tnew->Branch("Rv2_rt_p",Rv2_rt_p,"Rv2_rt_p[500]/D");
  tnew->Branch("Lu1_rt_p",Lu1_rt_p,"Lu1_rt_p[500]/D");
  tnew->Branch("Lu2_rt_p",Lu2_rt_p,"Lu2_rt_p[500]/D");
  tnew->Branch("Lv1_rt_p",Lv1_rt_p,"Lv1_rt_p[500]/D");
  tnew->Branch("Lv2_rt_p",Lv2_rt_p,"Lv2_rt_p[500]/D");

 
}

//////////////////////////////////////////////////////

void VDCt0::MakeHist(){

  
  min_rtime=1500.;
  max_rtime=3000.;
  bin_rtime=1500;
  double bin_rtdc=1.0;
  bin_rtime=(int )((max_rtime - min_rtime)/bin_rtdc);

  min_time=-100.0;
  max_time=500.0;
  double conv_tdc=0.5; //[ns/bin]
  bin_time=(int)((max_time-min_time)/conv_tdc); //0.5 [ns/bin]
  //  bin_time=(int)((max_time-min_time)/(conv_tdc*2)); //1.0 [ns/bin]

  
  double min_wire=0;
  double max_wire=(double)nwire;
  int bin_wire=nwire;

  hRu1=new TH2F("hRu1","",bin_wire,min_wire,max_wire,bin_rtime,min_rtime,max_rtime);
  set->SetTH2(hRu1,"RHRS U1 Wire vs Raw TDC hist","#wire","raw time [ch]");  
  hRu2=new TH2F("hRu2","",nwire+1,0,nwire,bin_rtime,min_rtime,max_rtime);
  set->SetTH2(hRu2,"RHRS U2 Wire vs Raw TDC hist","#wire","raw time [ch]");    
  hRv1=new TH2F("hRv1","",nwire+1,0,nwire,bin_rtime,min_rtime,max_rtime);
  set->SetTH2(hRv1,"RHRS V1 Wire vs Raw TDC hist","#wire","raw time [ch]");      
  hRv2=new TH2F("hRv2","",nwire+1,0,nwire,bin_rtime,min_rtime,max_rtime);
  set->SetTH2(hRv2,"RHRS V2 Wire vs Raw TDC hist","#wire","raw time [ch]");        
  hLu1=new TH2F("hLu1","",nwire+1,0,nwire,bin_rtime,min_rtime,max_rtime);
  set->SetTH2(hLu1,"LHRS U1 Wire vs Raw TDC hist","#wire","raw time [ch]");  
  hLu2=new TH2F("hLu2","",nwire+1,0,nwire,bin_rtime,min_rtime,max_rtime);
  set->SetTH2(hLu2,"LHRS U2 Wire vs Raw TDC hist","#wire","raw time [ch]");  
  hLv1=new TH2F("hLv1","",nwire+1,0,nwire,bin_rtime,min_rtime,max_rtime);
  set->SetTH2(hLv1,"LHRS V1 Wire vs Raw TDC hist","#wire","raw time [ch]");  
  hLv2=new TH2F("hLv2","",nwire+1,0,nwire,bin_rtime,min_rtime,max_rtime);  
  set->SetTH2(hLv2,"LHRS V2 Wire vs Raw TDC hist","#wire","raw time [ch]");



  hRu1_c=new TH2F("hRu1_c","",bin_wire,min_wire,max_wire,bin_time,min_time,max_time);
  set->SetTH2(hRu1_c,"RHRS U1 Wire vs Raw TDC hist","#wire","time [ns]");  
  hRu2_c=new TH2F("hRu2_c","",nwire+1,0,nwire,bin_time,min_time,max_time);
  set->SetTH2(hRu2_c,"RHRS U2 Wire vs Raw TDC hist","#wire"," time [ns]");    
  hRv1_c=new TH2F("hRv1_c","",nwire+1,0,nwire,bin_time,min_time,max_time);
  set->SetTH2(hRv1_c,"RHRS V1 Wire vs Raw TDC hist","#wire","time [ns]");      
  hRv2_c=new TH2F("hRv2_c","",nwire+1,0,nwire,bin_time,min_time,max_time);
  set->SetTH2(hRv2_c,"RHRS V2 Wire vs Raw TDC hist","#wire","time [ns]");        
  hLu1_c=new TH2F("hLu1_c","",nwire+1,0,nwire,bin_time,min_time,max_time);
  set->SetTH2(hLu1_c,"LHRS U1 Wire vs Raw TDC hist","#wire","time [ns]");  
  hLu2_c=new TH2F("hLu2_c","",nwire+1,0,nwire,bin_time,min_time,max_time);
  set->SetTH2(hLu2_c,"LHRS U2 Wire vs Raw TDC hist","#wire","time [ns]");  
  hLv1_c=new TH2F("hLv1_c","",nwire+1,0,nwire,bin_time,min_time,max_time);
  set->SetTH2(hLv1_c,"LHRS V1 Wire vs Raw TDC hist","#wire","time [ns]");  
  hLv2_c=new TH2F("hLv2_c","",nwire+1,0,nwire,bin_time,min_time,max_time);  
  set->SetTH2(hLv2_c,"LHRS V2 Wire vs Raw TDC hist","#wire"," time [ns]");



  
  for(int i=0;i<nwire;i++){
    //======== raw time hist =======================//
    hRu1_rtime[i]=new TH1D(Form("hRu1_rtime_%d",i),"",bin_rtime,min_rtime,max_rtime);
    set->SetTH1(hRu1_rtime[i],Form("RVDC U1 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");
    hRu2_rtime[i]=new TH1D(Form("hRu2_rtime_%d",i),"",bin_rtime,min_rtime,max_rtime);
    set->SetTH1(hRu2_rtime[i],Form("RVDC U2 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");
    hRv1_rtime[i]=new TH1D(Form("hRv1_rtime_%d",i),"",bin_rtime,min_rtime,max_rtime);
    set->SetTH1(hRv1_rtime[i],Form("RVDC V1 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");    
    hRv2_rtime[i]=new TH1D(Form("hRv2_rtime_%d",i),"",bin_rtime,min_rtime,max_rtime);
    set->SetTH1(hRv2_rtime[i],Form("RVDC V2 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");    
    hLu1_rtime[i]=new TH1D(Form("hLu1_rtime_%d",i),"",bin_rtime,min_rtime,max_rtime);
    set->SetTH1(hLu1_rtime[i],Form("LVDC U1 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");
    hLu2_rtime[i]=new TH1D(Form("hLu2_rtime_%d",i),"",bin_rtime,min_rtime,max_rtime);
    set->SetTH1(hLu2_rtime[i],Form("LVDC U2 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");    
    hLv1_rtime[i]=new TH1D(Form("hLv1_rtime_%d",i),"",bin_rtime,min_rtime,max_rtime);
    set->SetTH1(hLv1_rtime[i],Form("LVDC V1 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");
    hLv2_rtime[i]=new TH1D(Form("hLv2_rtime_%d",i),"",bin_rtime,min_rtime,max_rtime);
    set->SetTH1(hLv2_rtime[i],Form("LVDC V2 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");

    //======== time hist =======================//
    hRu1_time[i]=new TH1D(Form("hRu1_time_%d",i),"",bin_time,min_time,max_time);
    set->SetTH1(hRu1_time[i],Form("RVDC U1 wire-%d  tdc hist ",i),"time [ns]","Counts");
    hRu2_time[i]=new TH1D(Form("hRu2_time_%d",i),"",bin_time,min_time,max_time);
    set->SetTH1(hRu2_time[i],Form("RVDC U2 wire-%d  tdc hist ",i),"time [ns]","Counts");
    hRv1_time[i]=new TH1D(Form("hRv1_time_%d",i),"",bin_time,min_time,max_time);
    set->SetTH1(hRv1_time[i],Form("RVDC V1 wire-%d  tdc hist ",i),"time [ns]","Counts");    
    hRv2_time[i]=new TH1D(Form("hRv2_time_%d",i),"",bin_time,min_time,max_time);
    set->SetTH1(hRv2_time[i],Form("RVDC V2 wire-%d  tdc hist ",i),"time [ns]","Counts");    
    hLu1_time[i]=new TH1D(Form("hLu1_time_%d",i),"",bin_time,min_time,max_time);
    set->SetTH1(hLu1_time[i],Form("LVDC U1 wire-%d  tdc hist ",i),"time [ns]","Counts");
    hLu2_time[i]=new TH1D(Form("hLu2_time_%d",i),"",bin_time,min_time,max_time);
    set->SetTH1(hLu2_time[i],Form("LVDC U2 wire-%d  tdc hist ",i),"time [ns]","Counts");    
    hLv1_time[i]=new TH1D(Form("hLv1_time_%d",i),"",bin_time,min_time,max_time);
    set->SetTH1(hLv1_time[i],Form("LVDC V1 wire-%d  tdc hist ",i),"time [ns]","Counts");
    hLv2_time[i]=new TH1D(Form("hLv2_time_%d",i),"",bin_time,min_time,max_time);
    set->SetTH1(hLv2_time[i],Form("LVDC V2 wire-%d  tdc hist ",i),"time [ns]","Counts");    
    
  }


  
  //=============== raw time hist ==================//
  gRu1=new TGraphErrors();
  gRu1->SetTitle("RVDC U1 raw time Offset; #wire ; T0 [ch ]");
  gRu1->SetMarkerSize(1.0);
  gRu1->SetMarkerColor(2);
  gRu1->SetMarkerStyle(20);  
  gRu2=new TGraphErrors();
  gRu2->SetTitle("RVDC U2 raw time Offset; #wire ; T0 [ch ]");  
  gRu2->SetMarkerSize(1.0);
  gRu2->SetMarkerColor(2);
  gRu2->SetMarkerStyle(20);
  gRv1=new TGraphErrors();
  gRv1->SetTitle("RVDC V1 raw time Offset; #wire ; T0 [ch ]");  
  gRv1->SetMarkerSize(1.0);
  gRv1->SetMarkerColor(2);
  gRv1->SetMarkerStyle(20);
  gRv2=new TGraphErrors();
  gRv2->SetTitle("RVDC V2 raw time Offset; #wire ; T0 [ch ]");
  gRv2->SetMarkerSize(1.0);
  gRv2->SetMarkerColor(2);
  gRv2->SetMarkerStyle(20);    
  gLu1=new TGraphErrors();
  gLu1->SetTitle("LVDC U1 raw time Offset; #wire ; T0 [ch ]");  
  gLu1->SetMarkerSize(1.0);
  gLu1->SetMarkerColor(2);
  gLu1->SetMarkerStyle(20);    
  gLu2=new TGraphErrors();
  gLu2->SetTitle("LVDC U2 raw time Offset; #wire ; T0 [ch ]");    
  gLu2->SetMarkerSize(1.0);
  gLu2->SetMarkerColor(2);
  gLu2->SetMarkerStyle(20);      
  gLv1=new TGraphErrors();
  gLv1->SetTitle("LVDC V1 raw time Offset; #wire ; T0 [ch ]");  
  gLv1->SetMarkerSize(1.0);
  gLv1->SetMarkerColor(2);
  gLv1->SetMarkerStyle(20);      
  gLv2=new TGraphErrors();
  gLv2->SetTitle("LVDC V2 raw time Offset; #wire ; T0 [ch ]");  
  gLv2->SetMarkerSize(1.0);
  gLv2->SetMarkerColor(2);
  gLv2->SetMarkerStyle(20);

  //=========== time ====================//
  
  gRu1_c=new TGraphErrors();
  gRu1_c->SetTitle("RVDC U1 Offset; #wire ; T0 [ns]");
  gRu1_c->SetMarkerSize(1.0);
  gRu1_c->SetMarkerColor(2);
  gRu1_c->SetMarkerStyle(20);  
  gRu2_c=new TGraphErrors();
  gRu2_c->SetTitle("RVDC U2 Offset; #wire ; T0 [ns]");  
  gRu2_c->SetMarkerSize(1.0);
  gRu2_c->SetMarkerColor(2);
  gRu2_c->SetMarkerStyle(20);
  gRv1_c=new TGraphErrors();
  gRv1_c->SetTitle("RVDC V1 Offset; #wire ; T0 [ns]");  
  gRv1_c->SetMarkerSize(1.0);
  gRv1_c->SetMarkerColor(2);
  gRv1_c->SetMarkerStyle(20);
  gRv2_c=new TGraphErrors();
  gRv2_c->SetTitle("RVDC V2 Offset; #wire ; T0 [ns]");
  gRv2_c->SetMarkerSize(1.0);
  gRv2_c->SetMarkerColor(2);
  gRv2_c->SetMarkerStyle(20);    
  gLu1_c=new TGraphErrors();
  gLu1_c->SetTitle("LVDC U1 Offset; #wire ; T0 [ns]");  
  gLu1_c->SetMarkerSize(1.0);
  gLu1_c->SetMarkerColor(2);
  gLu1_c->SetMarkerStyle(20);    
  gLu2_c=new TGraphErrors();
  gLu2_c->SetTitle("LVDC U2 Offset; #wire ; T0 [ns]");    
  gLu2_c->SetMarkerSize(1.0);
  gLu2_c->SetMarkerColor(2);
  gLu2_c->SetMarkerStyle(20);      
  gLv1_c=new TGraphErrors();
  gLv1_c->SetTitle("LVDC V1 Offset; #wire ; T0 [ns]");  
  gLv1_c->SetMarkerSize(1.0);
  gLv1_c->SetMarkerColor(2);
  gLv1_c->SetMarkerStyle(20);      
  gLv2_c=new TGraphErrors();
  gLv2_c->SetTitle("LVDC V2 Offset; #wire ; T0 [ns]");  
  gLv2_c->SetMarkerSize(1.0);
  gLv2_c->SetMarkerColor(2);
  gLv2_c->SetMarkerStyle(20);  
  

  //------ T0 fitting function ----//


  for(int i=0; i<nwire;i++){

    fLu1_t0[i]=new TF1(Form("fLu1_t0_%d",i),"[0]*x+[1]",min_time,max_time);
    fLu1_t0[i]->SetLineColor(2);
    fLu1_t0[i]->SetNpx(2000);    
    fLu2_t0[i]=new TF1(Form("fLu2_t0_%d",i),"[0]*x+[1]",min_time,max_time);
    fLu2_t0[i]->SetLineColor(2);
    fLu2_t0[i]->SetNpx(2000);        
    fLv1_t0[i]=new TF1(Form("fLv1_t0_%d",i),"[0]*x+[1]",min_time,max_time);
    fLv1_t0[i]->SetLineColor(2);
    fLv1_t0[i]->SetNpx(2000);        
    fLv2_t0[i]=new TF1(Form("fLv2_t0_%d",i),"[0]*x+[1]",min_time,max_time);
    fLv2_t0[i]->SetLineColor(2);    
    fLv2_t0[i]->SetNpx(2000);    
    fRu1_t0[i]=new TF1(Form("fRu1_t0_%d",i),"[0]*x+[1]",min_time,max_time);
    fRu1_t0[i]->SetLineColor(2);
    fRu1_t0[i]->SetNpx(2000);    
    fRu2_t0[i]=new TF1(Form("fRu2_t0_%d",i),"[0]*x+[1]",min_time,max_time);
    fRu2_t0[i]->SetLineColor(2);
    fRu2_t0[i]->SetNpx(2000);        
    fRv1_t0[i]=new TF1(Form("fRv1_t0_%d",i),"[0]*x+[1]",min_time,max_time);
    fRv1_t0[i]->SetLineColor(2);
    fRv1_t0[i]->SetNpx(2000);        
    fRv2_t0[i]=new TF1(Form("fRv2_t0_%d",i),"[0]*x+[1]",min_time,max_time);
    fRv2_t0[i]->SetLineColor(2);
    fRv2_t0[i]->SetNpx(2000);
    
    fLu1_rt0[i]=new TF1(Form("fLu1_rt0_%d",i),"[0]*x+[1]",min_rtime,max_rtime);
    fLu1_rt0[i]->SetLineColor(2);
    fLu1_rt0[i]->SetNpx(2000);    
    fLu2_rt0[i]=new TF1(Form("fLu2_rt0_%d",i),"[0]*x+[1]",min_rtime,max_rtime);
    fLu2_rt0[i]->SetLineColor(2);
    fLu2_rt0[i]->SetNpx(2000);        
    fLv1_rt0[i]=new TF1(Form("fLv1_rt0_%d",i),"[0]*x+[1]",min_rtime,max_rtime);
    fLv1_rt0[i]->SetLineColor(2);
    fLv1_rt0[i]->SetNpx(2000);        
    fLv2_rt0[i]=new TF1(Form("fLv2_rt0_%d",i),"[0]*x+[1]",min_rtime,max_rtime);
    fLv2_rt0[i]->SetLineColor(2);
    fLv2_rt0[i]->SetNpx(2000);            
    fRu1_rt0[i]=new TF1(Form("fRu1_rt0_%d",i),"[0]*x+[1]",min_rtime,max_rtime);
    fRu1_rt0[i]->SetLineColor(2);
    fRu1_rt0[i]->SetNpx(2000);         
    fRu2_rt0[i]=new TF1(Form("fRu2_rt0_%d",i),"[0]*x+[1]",min_rtime,max_rtime);
    fRu2_rt0[i]->SetLineColor(2);
    fRu2_rt0[i]->SetNpx(2000);             
    fRv1_rt0[i]=new TF1(Form("fRv1_rt0_%d",i),"[0]*x+[1]",min_rtime,max_rtime);
    fRv1_rt0[i]->SetLineColor(2);
    fRv1_rt0[i]->SetNpx(2000);             
    fRv2_rt0[i]=new TF1(Form("fRv2_rt0_%d",i),"[0]*x+[1]",min_rtime,max_rtime);            
    fRv2_rt0[i]->SetLineColor(2);
    fRv2_rt0[i]->SetNpx(2000);         
  }

  
  
}

///////////////////////////////////////////////////

void VDCt0::Fill(){

  int counts=0;
  bool T1,T4,T5;


  int Ru1_nwire,Ru2_nwire,Rv1_nwire,Rv2_nwire;
  int Lu1_nwire,Lu2_nwire,Lv1_nwire,Lv2_nwire;

  
  ///===== Fill =======//
  //  int test=10000;
  //   ENum=test;

  for(int k=0;k<ENum;k++){
    for(int i=0;i<nmax;i++){
      Ru1_wire[i]=0.0;
      Ru2_wire[i]=0.0;
      Rv1_wire[i]=0.0;
      Rv2_wire[i]=0.0;
      Lu1_wire[i]=0.0;
      Lu2_wire[i]=0.0;
      Lv1_wire[i]=0.0;
      Lv2_wire[i]=0.0;
      Ru1_rtime[i]=0.0;
      Ru2_rtime[i]=0.0;
      Rv1_rtime[i]=0.0;
      Rv2_rtime[i]=0.0;
      Lu1_rtime[i]=0.0;
      Lu2_rtime[i]=0.0;
      Lv1_rtime[i]=0.0;
      Lv2_rtime[i]=0.0;
      Ru1_time[i]=0.0;
      Ru2_time[i]=0.0;
      Rv1_time[i]=0.0;
      Rv2_time[i]=0.0;
      Lu1_time[i]=0.0;
      Lu2_time[i]=0.0;
      Lv1_time[i]=0.0;
      Lv2_time[i]=0.0;            
    }

    for(int i=0;i<500;i++){
      Ru1_rt_p[i]=-1000.;
      Ru2_rt_p[i]=-1000.;
      Rv1_rt_p[i]=-1000.;
      Rv2_rt_p[i]=-1000.;
      Lu1_rt_p[i]=-1000.;
      Lu2_rt_p[i]=-1000.;
      Lv1_rt_p[i]=-1000.;
      Lv2_rt_p[i]=-1000.;      
    }
    
    T1=false;
    T4=false;
    T5=false;
    T->GetEntry(k);
    
    if(evtype==2.0)T1=true;
    if(evtype==16.0)T4=true;
    if(evtype==32)T5=true;

    Ru1_nwire=(int)Ru1_wire[0];
    Ru2_nwire=(int)Ru2_wire[0];
    Rv1_nwire=(int)Rv1_wire[0];
    Rv2_nwire=(int)Rv2_wire[0];
    Lu1_nwire=(int)Lu1_wire[0];
    Lu2_nwire=(int)Lu2_wire[0];
    Lv1_nwire=(int)Lv1_wire[0];
    Lv2_nwire=(int)Lv2_wire[0];

    

    if(T1){
    if(NLu1_wire>0 && 100>NLu1_wire)for(int i=0;i<NLu1_wire;i++){
	hLu1->Fill(Lu1_wire[i],Lu1_rtime[i]);
	Lu1_rt_p[(int)Lu1_wire[i]]=Lu1_rtime[i]-off[(int)Lu1_wire[i] ][4];    }
    if(NLu2_wire>0 && 100>NLu2_wire)for(int i=0;i<NLu2_wire;i++){
	hLu2->Fill(Lu2_wire[i],Lu2_rtime[i]);
	Lu2_rt_p[(int)Lu2_wire[i]]=Lu2_rtime[i]-off[(int)Lu2_wire[i] ][5];    }
    if(NLv1_wire>0 && 100>NLv1_wire)for(int i=0;i<NLv1_wire;i++){
	hLv1->Fill(Lv1_wire[i],Lv1_rtime[i]);
	Lv1_rt_p[(int)Lv1_wire[i]]=Lv1_rtime[i]-off[(int)Lv1_wire[i] ][6];    }
    if(NLv2_wire>0 && 100>NLv2_wire)for(int i=0;i<NLv2_wire;i++){
	hLv2->Fill(Lv2_wire[i],Lv2_rtime[i]);
      	Lv2_rt_p[(int)Lv2_wire[i]]=Lv2_rtime[i]-off[(int)Lv2_wire[i] ][7];    }

    if(NLu1_wire>0 && 100>NLu1_wire)for(int i=0;i<NLu1_wire;i++){hLu1_c->Fill(Lu1_wire[i],Lu1_time[i]*1.0e9);}
    if(NLu2_wire>0 && 100>NLu2_wire)for(int i=0;i<NLu2_wire;i++){hLu2_c->Fill(Lu2_wire[i],Lu2_time[i]*1.0e9);}
    if(NLv1_wire>0 && 100>NLv1_wire)for(int i=0;i<NLv1_wire;i++){hLv1_c->Fill(Lv1_wire[i],Lv1_time[i]*1.0e9);}
    if(NLv2_wire>0 && 100>NLv2_wire)for(int i=0;i<NLv2_wire;i++){hLv2_c->Fill(Lv2_wire[i],Lv2_time[i]*1.0e9);}
    

    }


  
    if(T4){

    if(NRu1_wire>0 && 100>NRu1_wire)for(int i=0;i<NRu1_wire;i++){
	hRu1->Fill(Ru1_wire[i],Ru1_rtime[i]);
	Ru1_rt_p[(int)Ru1_wire[i]]=Ru1_rtime[i]-off[(int)Ru1_wire[i] ][0];          }
    if(NRu2_wire>0 && 100>NRu2_wire)for(int i=0;i<NRu2_wire;i++){
	hRu2->Fill(Ru2_wire[i],Ru2_rtime[i]);
	Ru2_rt_p[(int)Ru2_wire[i]]=Ru2_rtime[i]-off[(int)Ru2_wire[i] ][1];          }
    if(NRv1_wire>0 && 100>NRv1_wire)for(int i=0;i<NRv1_wire;i++){
	hRv1->Fill(Rv1_wire[i],Rv1_rtime[i]);
      	Rv1_rt_p[(int)Ru1_wire[i]]=Rv1_rtime[i]-off[(int)Rv1_wire[i] ][2];          }
    if(NRv2_wire>0 && 100>NRv2_wire)for(int i=0;i<NRv2_wire;i++){
	hRv2->Fill(Rv2_wire[i],Rv2_rtime[i]);
      	Rv2_rt_p[(int)Rv2_wire[i]]=Rv2_rtime[i]-off[(int)Rv2_wire[i] ][3];          }

    
    if(NRu1_wire>0 && 100>NRu1_wire)for(int i=0;i<NRu1_wire;i++){hRu1_c->Fill(Ru1_wire[i],Ru1_time[i]*1.0e9);}
    if(NRu2_wire>0 && 100>NRu2_wire)for(int i=0;i<NRu2_wire;i++){hRu2_c->Fill(Ru2_wire[i],Ru2_time[i]*1.0e9);}
    if(NRv1_wire>0 && 100>NRv1_wire)for(int i=0;i<NRv1_wire;i++){hRv1_c->Fill(Rv1_wire[i],Rv1_time[i]*1.0e9);}
    if(NRv2_wire>0 && 100>NRv2_wire)for(int i=0;i<NRv2_wire;i++){hRv2_c->Fill(Rv2_wire[i],Rv2_time[i]*1.0e9);}


    
    }


    tnew->Fill();
    
    if(k%(ENum/10)==0){
      int fill=k/(ENum/10)*10;
      cout<<"Filled : "<<fill<<" % ("<<k<<" / "<<ENum<<")"<<endl;}



    
  }//end Fill

 ///===== Project =======//

  for(int i=0;i<nwire;i++){
    //========= raw time hist======================================//
    hRu1_rtime[i]=hRu1->ProjectionY(Form("hRu1_rtime_%d",i),i+1,i+1);
    set->SetTH1(hRu1_rtime[i],Form("RVDC U1 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");    
    hRu2_rtime[i]=hRu2->ProjectionY(Form("hRu2_rtime_%d",i),i+1,i+1);
    set->SetTH1(hRu2_rtime[i],Form("RVDC U2 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");        
    hRv1_rtime[i]=hRv1->ProjectionY(Form("hRv1_rtime_%d",i),i+1,i+1);
    set->SetTH1(hRv1_rtime[i],Form("RVDC V1 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");        
    hRv2_rtime[i]=hRv2->ProjectionY(Form("hRv2_rtime_%d",i),i+1,i+1);
    set->SetTH1(hRv2_rtime[i],Form("RVDC V2 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");        
    hLu1_rtime[i]=hLu1->ProjectionY(Form("hLu1_rtime_%d",i),i+1,i+1);
    set->SetTH1(hLu1_rtime[i],Form("LVDC U1 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");        
    hLu2_rtime[i]=hLu2->ProjectionY(Form("hLu2_rtime_%d",i),i+1,i+1);
    set->SetTH1(hLu2_rtime[i],Form("LVDC U2 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");        
    hLv1_rtime[i]=hLv1->ProjectionY(Form("hLv1_rtime_%d",i),i+1,i+1);
    set->SetTH1(hLv1_rtime[i],Form("LVDC V1 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");
    hLv2_rtime[i]=hLv2->ProjectionY(Form("hLv2_rtime_%d",i),i+1,i+1);
    set->SetTH1(hLv2_rtime[i],Form("LVDC V2 wire-%d raw tdc hist ",i),"rawtime [ch]","Counts");


    //=========  time hist======================================//
    hRu1_time[i]=hRu1_c->ProjectionY(Form("hRu1_time_%d",i),i+1,i+1);
    set->SetTH1(hRu1_time[i],Form("RVDC U1 wire-%d  tdc hist ",i),"time [ns]","Counts");    
    hRu2_time[i]=hRu2_c->ProjectionY(Form("hRu2_time_%d",i),i+1,i+1);
    set->SetTH1(hRu2_time[i],Form("RVDC U2 wire-%d  tdc hist ",i),"time [ns]","Counts");        
    hRv1_time[i]=hRv1_c->ProjectionY(Form("hRv1_time_%d",i),i+1,i+1);
    set->SetTH1(hRv1_time[i],Form("RVDC V1 wire-%d  tdc hist ",i),"time [ns]","Counts");        
    hRv2_time[i]=hRv2_c->ProjectionY(Form("hRv2_time_%d",i),i+1,i+1);
    set->SetTH1(hRv2_time[i],Form("RVDC V2 wire-%d  tdc hist ",i),"time [ns]","Counts");        
    hLu1_time[i]=hLu1_c->ProjectionY(Form("hLu1_time_%d",i),i+1,i+1);
    set->SetTH1(hLu1_time[i],Form("LVDC U1 wire-%d  tdc hist ",i),"time [ns]","Counts");        
    hLu2_time[i]=hLu2_c->ProjectionY(Form("hLu2_time_%d",i),i+1,i+1);
    set->SetTH1(hLu2_time[i],Form("LVDC U2 wire-%d  tdc hist ",i),"time [ns]","Counts");        
    hLv1_time[i]=hLv1_c->ProjectionY(Form("hLv1_time_%d",i),i+1,i+1);
    set->SetTH1(hLv1_time[i],Form("LVDC V1 wire-%d  tdc hist ",i),"time [ns]","Counts");
    hLv2_time[i]=hLv2_c->ProjectionY(Form("hLv2_time_%d",i),i+1,i+1);
    set->SetTH1(hLv2_time[i],Form("LVDC V2 wire-%d  tdc hist ",i),"time [ns]","Counts");            


   
  }

  
  cout<<"Filled !! "<<endl;
  
  }


/////////////////////////////////////////////////////

void VDCt0::Deft0(string ifname){

  ifstream ifparam(ifname.c_str());
  if (ifparam.fail()){ cerr << "failed open files " <<ifname<<endl;}
  cout<<"def file : "<<ifname<<endl;

  
  int plane=-1;
  int i=0;
  string buf;
  
  while( getline(ifparam,buf) ){
 
    if( buf[0]=='#' ){ i=0; plane++;}
    
    if( ifparam.eof() || i>nwire-1)continue;
    ifparam >> T0_def[i][plane] >> T0_def[i+1][plane] >> T0_def[i+2][plane] >> T0_def[i+3][plane] >> T0_def[i+4][plane] >> T0_def[i+5][plane] >> T0_def[i+6][plane] >> T0_def[i+7][plane];
    //    cout<<"T0: "<<"i " <<i<<" "<< T0_def[i][plane] <<" "<< T0_def[i+1][plane] <<" "<< T0_def[i+2][plane] <<" "<< T0_def[i+3][plane] <<" "<< T0_def[i+4][plane] <<" "<< T0_def[i+5][plane] <<" "<< T0_def[i+6][plane] <<" "<< T0_def[i+7][plane]<<endl;    
    i=i+8;


    
  }//end while dat file

  

}

////////////////////////////////////////////////////

  double VDCt0::Findt0(bool rarm,char* plane,int wire){
      double dy,y,yb,ya;

    t0_min=0.0; t0_max=0.0;
    t0_1=0.0; t0_2=0.0;
    slope_1=0.0; slope_2=0.0;
    dt0=0.0;

    TH1D* hnew;
  if(rarm && plane=="U1")hnew=(TH1D*)hRu1_rtime[wire]->Clone();
  else if(rarm && plane=="U2")hnew=(TH1D*)hRu2_rtime[wire]->Clone();
  else if(rarm && plane=="V1")hnew=(TH1D*)hRv1_rtime[wire]->Clone();
  else if(rarm && plane=="V2")hnew=(TH1D*)hRv2_rtime[wire]->Clone();
  else if(rarm==0 && plane=="U1")hnew=(TH1D*)hLu1_rtime[wire]->Clone();
  else if(rarm==0 && plane=="U2")hnew=(TH1D*)hLu2_rtime[wire]->Clone();
  else if(rarm==0 && plane=="V1")hnew=(TH1D*)hLv1_rtime[wire]->Clone();
  else if(rarm==0 && plane=="V2")hnew=(TH1D*)hLv2_rtime[wire]->Clone();
  else {cout<<"Failed to call hist "<<endl;}
  
  hnew->SetName("hnew");
  
    if (!hnew) {
        Error("Findt0","Empty time spectrum to calibrate???");
        return 0;
    }
    double x;
    double y0=0.0;
    Double_t dx = hnew->GetBinWidth(1);
    Double_t slope_min = 0.0, slope = 0.0;
    Double_t sdbin = -1; // Bin number with steepest descent
    double slope_0=0.0;
    Int_t nbins = hnew->GetNbinsX();
    Int_t maxbin = hnew->GetMaximumBin();
    Double_t maxcont = hnew->GetBinContent(maxbin);
    int bin_rmin;
    if(rarm==0)bin_rmin=hnew->GetXaxis()->FindBin(2900);
    if(rarm)bin_rmin=hnew->GetXaxis()->FindBin(2600);
    int kmax=5;
    int jmax=3;
    int nelement=kmax*jmax;    
    double slope_el[kmax];
    double slope_av;
    double slopes[nbins];
    double n[nwire];

    //    for (Int_t i=maxbin+5; i<=nbins; i++) {
        for (Int_t i=bin_rmin; i<=nbins; i++) {

	  if (hnew->GetBinContent(i)>0.5*maxcont) continue;
        if (hnew->GetBinContent(i-1)>0.5*maxcont) continue;
	
        if (i>bin_rmin+1 && i<nbins) {
            slope = (hnew->GetBinContent(i+1)-hnew->GetBinContent(i-1))/(2.*dx);

	    
	    //========= Average of slope ====================//
	    if(i>kmax && i<nbins-kmax){
	      slope_av=0.0;
	      for(int k=0;k<kmax;k++){

	      slope_el[k]=(hnew->GetBinContent(i+1)-hnew->GetBinContent(i-k+1))/(2.*dx);
	      slope_av += slope_el[k]/kmax;


		    }
	    }
	    //================================================//
	    
        } else if (i==1) {
            slope = (hnew->GetBinContent(i+1)-hnew->GetBinContent(i))/dx;
        } else if (i==nbins) {
            slope = (hnew->GetBinContent(i)-hnew->GetBinContent(i-1))/dx;
        }

        if (TMath::Abs(slope)<1.e-4) slope = 0.;

	slopes[i]=slope;

          if (slope<slope_min) {
	//      if (slope_av<slope_min) {
	    
	  if(i>3 && i<nbins-3){
            slope_min = slope;
            sdbin = i;
	  }
	}
    }

    Double_t t0 = 0.;
    if (slope_min<-1 && hnew->GetEntries()>1000.) {


	
      t0 = hnew->GetBinCenter(sdbin) -
	  hnew->GetBinContent(sdbin)/slope_min;

	
	x =  hnew->GetBinCenter(sdbin);
	y  = hnew->GetBinContent(sdbin);
	yb = hnew->GetBinContent(sdbin-1);
	ya = hnew->GetBinContent(sdbin+1);
	slope = (ya-yb)/(2.*dx);
	y0= y - slope*x;

	
	//       	if(ya==0 || yb==0 || slope_min==0 || y==0 || dt0<-1.0)dt0=100.0;
	//	else{dt0=(t0-x)*(sqrt(y)/y+(sqrt(yb)+sqrt(ya))/fabs(ya-yb));	}
	double dy;
	if(y==0)dy=0.0;
	else dy=1./sqrt(y);
	dt0=y/fabs(slope)*sqrt(  dy +pow( (sqrt(ya) + sqrt(yb))/fabs(slope) ,2 )  );
	
	//       	cout<<"dt0: "<<dt0<<" yn "<<y<<" xn "<<x<<" a "<<slope<<" yn+1 "<<ya<<" yn-1 "<<yb<<" calc dt0 "<<DT0<<endl;
	  //y/fabs(slope)*sqrt(  dy +pow( (sqrt(ya) + sqrt(yb))/fabs(slope) ,2 )  )<<endl;
    }


  if(rarm && plane=="U1")fRu1_rt0[wire]->SetParameters(slope_min,y0);
  else if(rarm && plane=="U2")fRu2_rt0[wire]->SetParameters(slope_min,y0);
  else if(rarm && plane=="V1")fRv1_rt0[wire]->SetParameters(slope_min,y0);
  else if(rarm && plane=="V2")fRv2_rt0[wire]->SetParameters(slope_min,y0);
  else if(rarm==0 && plane=="U1")fLu1_rt0[wire]->SetParameters(slope_min,y0);
  else if(rarm==0 && plane=="U2")fLu2_rt0[wire]->SetParameters(slope_min,y0);
  else if(rarm==0 && plane=="V1")fLv1_rt0[wire]->SetParameters(slope_min,y0);
  else if(rarm==0 && plane=="V2")fLv2_rt0[wire]->SetParameters(slope_min,y0);

    
    
    return t0;
}

//////////////////////////////////////////////////////


  double VDCt0::Findt0_time(bool rarm,char* plane,int wire){
      double dy,y,yb,ya;


    t0_min=0.0; t0_max=0.0;
    t0_1=0.0; t0_2=0.0;
    slope_1=0.0; slope_2=0.0;
    dt0=0.0;

    TH1D* hnew;
  if(rarm && plane=="U1")hnew=(TH1D*)hRu1_time[wire]->Clone();
  else if(rarm && plane=="U2")hnew=(TH1D*)hRu2_time[wire]->Clone();
  else if(rarm && plane=="V1")hnew=(TH1D*)hRv1_time[wire]->Clone();
  else if(rarm && plane=="V2")hnew=(TH1D*)hRv2_time[wire]->Clone();
  else if(rarm==0 && plane=="U1")hnew=(TH1D*)hLu1_time[wire]->Clone();
  else if(rarm==0 && plane=="U2")hnew=(TH1D*)hLu2_time[wire]->Clone();
  else if(rarm==0 && plane=="V1")hnew=(TH1D*)hLv1_time[wire]->Clone();
  else if(rarm==0 && plane=="V2")hnew=(TH1D*)hLv2_time[wire]->Clone();
  else {cout<<"Failed to call hist "<<endl;}
  
  hnew->SetName("hnew");
  
    if (!hnew) {
        Error("Findt0","Empty time spectrum to calibrate???");
        return 0;
    }
    double x=0.0;
    double y0=0.0;
    Double_t dx =(double)hnew->GetBinWidth(1);
    Double_t slope_max = 0.0;
    Double_t  slope = 0.0;
    Double_t sdbin = -1; // Bin number with steepest descent
    double slope_0=0.0;
    Int_t nbins = hnew->GetNbinsX();
    Int_t maxbin = hnew->GetMaximumBin();
    Double_t maxcont = hnew->GetBinContent(maxbin);
    int bin_tmax  = hnew-> GetXaxis()->FindBin(50.);
    int kmax=5;
    int jmax=3;
    int nelement=kmax*jmax;    
    double slope_el[kmax];
    double slope_av;
    double slopes[nbins];
    double n[nwire];

    //    for (Int_t i=0; i<=maxbin; i++) {
    for (Int_t i=0; i<=bin_tmax; i++) {
        if (hnew->GetBinContent(i)>0.5*maxcont) continue;
        if (hnew->GetBinContent(i+1)>0.5*maxcont) continue;
	
	//        if (i>1 && i<maxbin) {
	        if (i>1 && i<bin_tmax) {
            slope = (hnew->GetBinContent(i+1)-hnew->GetBinContent(i-1))/(2.*dx);
	    
        } else if (i==1) {
            slope = (hnew->GetBinContent(i+1)-hnew->GetBinContent(i))/dx;
	    //        } else if (i==maxbin) {
		} else if (i==bin_tmax) {
            slope = (hnew->GetBinContent(i)-hnew->GetBinContent(i-1))/dx;
        }

        if (TMath::Abs(slope)<1.e-4) slope = 0.;

          if (slope>slope_max) {

            slope_max = slope;
            sdbin = i;

	  }
    }

    
    Double_t t0 = 0.;
    if (slope_max>+1.0 && hnew->GetEntries()>1000.) {

      t0 = hnew->GetBinCenter(sdbin) -
	hnew->GetBinContent(sdbin)/slope_max;

	
      x =  (double)hnew->GetBinCenter(sdbin);
      y  = (double)hnew->GetBinContent(sdbin);
      yb = (double)hnew->GetBinContent(sdbin-1);
      ya = (double)hnew->GetBinContent(sdbin+1);
      slope = (ya-yb)/(2.*dx);
      
	y0= y - slope_max*x;

	//	dt0=fabs((t0-x)*(sqrt(y)/y+(sqrt(yb)+sqrt(ya))/fabs(ya-yb)));

	double dy;
	if(y==0)dy=0.0;
	else dy=1./sqrt(y);
	dt0=y/fabs(slope)*sqrt(  dy +pow( (sqrt(ya) + sqrt(yb))/fabs(slope) ,2 )  );
	
    }

    
    if(rarm && plane=="U1"){
      Ru1_p0[wire]=slope_max;
      Ru1_p1[wire]=y0;
      fRu1_t0[wire]->SetParameters(slope_max,y0);
    } else if(rarm && plane=="U2"){
      Ru2_p0[wire]=slope_max;
      Ru2_p1[wire]=y0;
      fRu2_t0[wire]->SetParameters(slope_max,y0);
    } else if(rarm && plane=="V1"){
      Rv1_p0[wire]=slope_max;
      Rv1_p1[wire]=y0;
      fRv1_t0[wire]->SetParameters(slope_max,y0);
    }else if(rarm && plane=="V2"){
      Rv2_p0[wire]=slope_max;
      Rv2_p1[wire]=y0;
      fRv2_t0[wire]->SetParameters(slope_max,y0);
    }else if(rarm==0 && plane=="U1"){
      Lu1_p0[wire]=slope_max;
      Lu1_p1[wire]=y0;
      fLu1_t0[wire]->SetParameters(slope_max,y0);
    }else if(rarm==0 && plane=="U2"){
      Lu2_p0[wire]=slope_max;
      Lu2_p1[wire]=y0;
      fLu2_t0[wire]->SetParameters(slope_max,y0);
    }else if(rarm==0 && plane=="V1"){
      Lv1_p0[wire]=slope_max;
      Lv1_p1[wire]=y0;
      fLv1_t0[wire]->SetParameters(slope_max,y0);
    }else if(rarm==0 && plane=="V2"){
      Lv2_p0[wire]=slope_max;
      Lv2_p1[wire]=y0;      
      fLv2_t0[wire]->SetParameters(slope_max,y0);}
  
   t0offset=t0;
  
      return t0;
  }


/////////////////////////////////////////////////////

void VDCt0::SetRawt0(){

  //==== Fill ====//
  
    for(int i=0;i<nwire;i++){
     Ru1t0[i]=Findt0(true,"U1",i);
     if(Ru1t0[i]==0)Ru1t0[i]=T0_def[i][0];
     Ru1t0_err[i]=dt0;
     gRu1->SetPoint(i,i,Ru1t0[i]);
     gRu1->SetPointError(i,0,Ru1t0_err[i]);     

     Ru2t0[i]=Findt0(true,"U2",i);     
     if(Ru2t0[i]==0)Ru2t0[i]=T0_def[i][1];
     Ru2t0_err[i]=dt0;     
     gRu2->SetPoint(i,i,Ru2t0[i]);
     gRu2->SetPointError(i,0,Ru2t0_err[i]);     

     Rv1t0[i]=Findt0(true,"V1",i);
     if(Rv1t0[i]==0)Rv1t0[i]=T0_def[i][2];
     Rv1t0_err[i]=dt0;
     gRv1->SetPoint(i,i,Rv1t0[i]);
     gRv1->SetPointError(i,0,Rv1t0_err[i]);     
     
     Rv2t0[i]=Findt0(true,"V2",i);
     if(Rv2t0[i]==0)Rv2t0[i]=T0_def[i][3];
     Rv2t0_err[i]=dt0;
     gRv2->SetPoint(i,i,Rv2t0[i]);
     gRv2->SetPointError(i,0,Rv2t0_err[i]);     
     
     Lu1t0[i]=Findt0(false,"U1",i);
     if(Lu1t0[i]==0)Lu1t0[i]=T0_def[i][4];
     Lu1t0_err[i]=dt0;
     gLu1->SetPoint(i,i,Lu1t0[i]);
     gLu1->SetPointError(i,0,Lu1t0_err[i]);     
     
     Lu2t0[i]=Findt0(false,"U2",i);
     if(Lu2t0[i]==0)Lu2t0[i]=T0_def[i][5];     
     Lu2t0_err[i]=dt0;
     gLu2->SetPoint(i,i,Lu2t0[i]);
     gLu2->SetPointError(i,0,Lu2t0_err[i]);     
     
     Lv1t0[i]=Findt0(false,"V1",i);
     if(Lv1t0[i]==0)Lv1t0[i]=T0_def[i][6];
     Lv1t0_err[i]=dt0;          
     gLv1->SetPoint(i,i,Lv1t0[i]);
     gLv1->SetPointError(i,0,Lv1t0_err[i]);     

     Lv2t0[i]=Findt0(false,"V2",i);
     if(Lv2t0[i]==0)Lv2t0[i]=T0_def[i][7];     
     Lv2t0_err[i]=dt0;     
     gLv2->SetPoint(i,i,Lv2t0[i]);
     gLv2->SetPointError(i,0,Lv2t0_err[i]);     

    }



}

/////////////////////////////////////////////////////

void VDCt0::Sett0(){

  //==== Fill ====//

    for(int i=0;i<nwire;i++){
      Ru1t0_c[i]=0.0;
      Ru2t0_c[i]=0.0;      
      Rv1t0_c[i]=0.0;
      Rv2t0_c[i]=0.0;
      Lu1t0_c[i]=0.0;
      Lu2t0_c[i]=0.0;      
      Lv1t0_c[i]=0.0;
      Lv2t0_c[i]=0.0;      
    }
  
    for(int i=0;i<nwire;i++){
     Ru1t0_c[i]=Findt0_time(true,"U1",i);
     if(Ru1t0_c[i]==0)Ru1t0_c[i]=100.;
     Ru1t0_err_c[i]=dt0;
     gRu1_c->SetPoint(i,i,Ru1t0_c[i]);
     gRu1_c->SetPointError(i,0,Ru1t0_err_c[i]);     

     Ru2t0_c[i]=Findt0_time(true,"U2",i);     
     if(Ru2t0_c[i]==0)Ru2t0_c[i]=100.;
     Ru2t0_err_c[i]=dt0;     
     gRu2_c->SetPoint(i,i,Ru2t0_c[i]);
     gRu2_c->SetPointError(i,0,Ru2t0_err_c[i]);     

     Rv1t0_c[i]=Findt0_time(true,"V1",i);
     if(Rv1t0_c[i]==0)Rv1t0_c[i]=100.;
     Rv1t0_err_c[i]=dt0;
     gRv1_c->SetPoint(i,i,Rv1t0_c[i]);
     gRv1_c->SetPointError(i,0,Rv1t0_err_c[i]);     
     
     Rv2t0_c[i]=Findt0_time(true,"V2",i);
     if(Rv2t0_c[i]==0)Rv2t0_c[i]=100.;
     Rv2t0_err_c[i]=dt0;
     gRv2_c->SetPoint(i,i,Rv2t0_c[i]);
     gRv2_c->SetPointError(i,0,Rv2t0_err_c[i]);     
     
     Lu1t0_c[i]=Findt0_time(false,"U1",i);
     if(Lu1t0_c[i]==0)Lu1t0_c[i]=100.;
     Lu1t0_err_c[i]=dt0;
     gLu1_c->SetPoint(i,i,Lu1t0_c[i]);
     gLu1_c->SetPointError(i,0,Lu1t0_err_c[i]);     
     
     Lu2t0_c[i]=Findt0_time(false,"U2",i);
     if(Lu2t0_c[i]==0)Lu2t0_c[i]=100.;     
     Lu2t0_err_c[i]=dt0;
     gLu2_c->SetPoint(i,i,Lu2t0_c[i]);
     gLu2_c->SetPointError(i,0,Lu2t0_err_c[i]);     
     
     Lv1t0_c[i]=Findt0_time(false,"V1",i);
     if(Lv1t0_c[i]==0)Lv1t0_c[i]=100.;
     Lv1t0_err_c[i]=dt0;          
     gLv1_c->SetPoint(i,i,Lv1t0_c[i]);
     gLv1_c->SetPointError(i,0,Lv1t0_err_c[i]);     

     Lv2t0_c[i]=Findt0_time(false,"V2",i);
     if(Lv2t0_c[i]==0)Lv2t0_c[i]=100.;     
     Lv2t0_err_c[i]=dt0;     
     gLv2_c->SetPoint(i,i,Lv2t0_c[i]);
     gLv2_c->SetPointError(i,0,Lv2t0_err_c[i]);     

    }

}

//////////////////////////////////////////////////////

void VDCt0::Write(string ofname){



  //  string ofname_main = ofname.substr(0,21);
  //  

  string ofname_main=ofname;
  ofname= ofname_main + ".dat"; 
  string ofname_err =ofname_main + "_err.dat";
  string ofname_def="./param/def_t0.dat";
  ofstream ofs(ofname.c_str());
  ofstream ofs_err(ofname_err.c_str());

  cout<<"Pram file : "<<ofname<<endl;
  cout<<"Param Error : "<<ofname_err<<endl;
  
  double T0,T0_min,T0_max;
  double dT0;
  string vdc_name;
  
  
  
  for(int j=0;j<8;j++){
    if(j==0)vdc_name="# RVDC U1 t0 parameters ";
    else if(j==1)vdc_name="# RVDC U2 t0 parameters ";
    else if(j==2)vdc_name="# RVDC V1 t0 parameters ";
    else if(j==3)vdc_name="# RVDC V2 t0 parameters ";
    else if(j==4)vdc_name="# LVDC U1 t0 parameters ";
    else if(j==5)vdc_name="# LVDC U2 t0 parameters ";
    else if(j==6)vdc_name="# LVDC V1 t0 parameters ";
    else if(j==7)vdc_name="# LVDC V2 t0 parameters ";
    else {cout<<"faled to write"<<endl; break;}
         ofs << vdc_name  <<endl;
	 ofs_err << vdc_name <<"T0 Error "<<endl;    
    for(int i=0;i<nwire;i++){
      if(j==0){T0= Ru1t0[i]; dT0=Ru1t0_err[i];}
      else if(j==1){T0= Ru2t0[i]; dT0=Ru2t0_err[i];}
      else if(j==2){T0= Rv1t0[i]; dT0=Rv1t0_err[i];}
      else if(j==3){T0= Rv2t0[i]; dT0=Rv2t0_err[i];}
      else if(j==4){T0= Lu1t0[i]; dT0=Lu1t0_err[i];}
      else if(j==5){T0= Lu2t0[i]; dT0=Lu2t0_err[i];}
      else if(j==6){T0= Lv1t0[i]; dT0=Lv1t0_err[i];}
      else if(j==7){T0= Lv2t0[i]; dT0=Lv2t0_err[i];}
      else {cout<<"faled to write"<<endl; break;}      
      
      if(fabs(T0)<10000) ofs << T0 <<" ";
      else ofs << 0.0 <<" ";
      if(fabs(T0_min)<1000)  ofs_err << dT0 <<" ";    
      else ofs_err << 0.0 <<" ";          
   
    if((i+1)%8==0){
      ofs << endl;
      ofs_err << endl;}
    
    }
  }
  

  
  ofs.close();
  ofs_err.close();
  
}


//////////////////////////////////////////////////////

void VDCt0::Draw(){

  cout<<"==================================="<<endl;
  cout<<"========= Draw ===================="<<endl;
  cout<<"==================================="<<endl;  


  for(int j=0;j<11;j++){
    c0[j] =new TCanvas(Form("c0_%d",j),Form("LVDC-U1_%d",j));
    c0[j]->Divide(4,8);
    c1[j] =new TCanvas(Form("c1_%d",j),Form("LVDC-U2_%d",j));
    c1[j]->Divide(4,8);
    c2[j] =new TCanvas(Form("c2_%d",j),Form("LVDC-V1_%d",j));
    c2[j]->Divide(4,8);
    c3[j] =new TCanvas(Form("c3_%d",j),Form("LVDC-V2_%d",j));
    c3[j]->Divide(4,8);

    
    for(int i=0;i<32;i++){
    c0[j]->cd(i+1);
    hLu1_rtime[32*j+i]->Draw();
    line_rtime[32*j+i]= new TLine(Lu1t0[32*j+i],0,Lu1t0[32*j+i],hLu1_rtime[32*j+i]->GetMaximum());
    line_rtime[32*j+i]->SetLineColor(kRed);
    line_rtime[32*j+i]->SetLineWidth(2);
    line_rtime[32*j+i]->Draw("same");        
    c1[j]->cd(i+1);
    hLu2_rtime[32*j+i]->Draw();
    line_rtime[32*j+i]= new TLine(Lu2t0[32*j+i],0,Lu2t0[32*j+i],hLu2_rtime[32*j+i]->GetMaximum());
    line_rtime[32*j+i]->SetLineColor(kRed);
    line_rtime[32*j+i]->SetLineWidth(2);
    line_rtime[32*j+i]->Draw("same");        
    c2[j]->cd(i+1);
    hLv1_rtime[32*j+i]->Draw();
    line_rtime[32*j+i]= new TLine(Lv1t0[32*j+i],0,Lv1t0[32*j+i],hLv1_rtime[32*j+i]->GetMaximum());
    line_rtime[32*j+i]->SetLineColor(kRed);
    line_rtime[32*j+i]->SetLineWidth(2);
    line_rtime[32*j+i]->Draw("same");    
    c3[j]->cd(i+1);
    hLv2_rtime[32*j+i]->Draw();
    line_rtime[32*j+i]= new TLine(Lv2t0[32*j+i],0,Lv2t0[32*j+i],hLv2_rtime[32*j+i]->GetMaximum());
    line_rtime[32*j+i]->SetLineColor(kRed);
    line_rtime[32*j+i]->SetLineWidth(2);
    line_rtime[32*j+i]->Draw("same");        
    

    }
  }
  

  //================ time hist Draw ==============================//

  for(int j=0;j<11;j++){
    c4[j] =new TCanvas(Form("c4_%d",j),Form("LVDC-U1_%d",j));
    c4[j]->Divide(4,8);
    c5[j] =new TCanvas(Form("c5_%d",j),Form("LVDC-U2_%d",j));
    c5[j]->Divide(4,8);
    c6[j] =new TCanvas(Form("c6_%d",j),Form("LVDC-V1_%d",j));
    c6[j]->Divide(4,8);
    c7[j] =new TCanvas(Form("c7_%d",j),Form("LVDC-V2_%d",j));
    c7[j]->Divide(4,8);

    
    for(int i=0;i<32;i++){
    c4[j]->cd(i+1);
    hLu1_time[32*j+i]->Draw();
    line_time[32*j+i]= new TLine(Lu1t0_c[32*j+i],0,Lu1t0_c[32*j+i],hLu1_time[32*j+i]->GetMaximum());
    line_time[32*j+i]->SetLineColor(kRed);
    line_time[32*j+i]->SetLineWidth(2);
    line_time[32*j+i]->Draw("same");        
    c5[j]->cd(i+1);
    hLu2_time[32*j+i]->Draw();
    line_time[32*j+i]= new TLine(Lu2t0_c[32*j+i],0,Lu2t0_c[32*j+i],hLu2_time[32*j+i]->GetMaximum());
    line_time[32*j+i]->SetLineColor(kRed);
    line_time[32*j+i]->SetLineWidth(2);
    line_time[32*j+i]->Draw("same");        
    c6[j]->cd(i+1);
    hLv1_time[32*j+i]->Draw();
    line_time[32*j+i]= new TLine(Lv1t0_c[32*j+i],0,Lv1t0_c[32*j+i],hLv1_time[32*j+i]->GetMaximum());
    line_time[32*j+i]->SetLineColor(kRed);
    line_time[32*j+i]->SetLineWidth(2);
    line_time[32*j+i]->Draw("same");    
    c7[j]->cd(i+1);
    hLv2_time[32*j+i]->Draw();
    line_time[32*j+i]= new TLine(Lv2t0_c[32*j+i],0,Lv2t0_c[32*j+i],hLv2_time[32*j+i]->GetMaximum());
    line_time[32*j+i]->SetLineColor(kRed);
    line_time[32*j+i]->SetLineWidth(2);
    line_time[32*j+i]->Draw("same");        
    

    }
  }
  




  
  c10=new TCanvas("c10","c10");  
  c10->Divide(2,2);
  c10->cd(1);
  gLu1->Draw("AP");
  c10->cd(2);
  gLu2->Draw("AP");  
  c10->cd(3);
  gLv1->Draw("AP");
  c10->cd(4);
  gLv2->Draw("AP");
  
  cout<<" Drawn Pictures "<<endl;
  
}

/////////////////////////////////////////////////////

void VDCt0::Print(string ofname){

  cout<<"Print is starting "<<endl;
  cout<<"pdf name : "<<ofname<<endl;

  for(int j=0;j<4;j++){
  for(int i=0;i<11;i++){
    if(i==0 && j==0)c0[i]->Print(Form("%s[",ofname.c_str()));
    if(j==0)c0[i]->Print(Form("%s",ofname.c_str()));
    if(j==1)c1[i]->Print(Form("%s",ofname.c_str()));
    if(j==2)c2[i]->Print(Form("%s",ofname.c_str()));
    if(j==3)c3[i]->Print(Form("%s",ofname.c_str()));    
    //    if(i==10 && j==3)c3[i]->Print(Form("%s]",ofname.c_str()));
  }
  }
  c10->Print(Form("%s",ofname.c_str()));      
  c10->Print(Form("%s]",ofname.c_str()));
  cout<<"Print is done !"<<endl;

}

/////////////////////////////////////////////////////

void VDCt0::Print_c(string ofname){

  cout<<"Print is starting "<<endl;
  cout<<"pdf name : "<<ofname<<endl;

  for(int j=0;j<4;j++){
  for(int i=0;i<11;i++){
    if(i==0 && j==0)c4[i]->Print(Form("%s[",ofname.c_str()));
    if(j==0)c4[i]->Print(Form("%s",ofname.c_str()));
    if(j==1)c5[i]->Print(Form("%s",ofname.c_str()));
    if(j==2)c6[i]->Print(Form("%s",ofname.c_str()));
    if(j==3)c7[i]->Print(Form("%s",ofname.c_str()));    
   if(i==10 && j==3)c7[i]->Print(Form("%s]",ofname.c_str()));
  }
  }
  //  c10->Print(Form("%s",ofname.c_str()));      
  //  c10->Print(Form("%s]",ofname.c_str()));
  cout<<"Print is done !"<<endl;

}

/////////////////////////////////////////////////////


void VDCt0::MakeRoot(string ofname){

  
  //======= Write ========//


  for(int i=0;i<nwire;i++){

    hRu1_rtime[i]->Write();
    hRu2_rtime[i]->Write();
    hRv1_rtime[i]->Write();
    hRv2_rtime[i]->Write();
    hLu1_rtime[i]->Write();
    hLu2_rtime[i]->Write();
    hLv1_rtime[i]->Write();
    hLv2_rtime[i]->Write();
    fLu1_rt0[i]->Write();
    fLu2_rt0[i]->Write();
    fLv1_rt0[i]->Write();
    fLv2_rt0[i]->Write();    
    fRu1_rt0[i]->Write();
    fRu2_rt0[i]->Write();
    fRv1_rt0[i]->Write();
    fRv2_rt0[i]->Write();    
  }
  
  hLu1->Write();
  hLu2->Write();  
  hLv1->Write();  
  hLv2->Write();
  hRu1->Write();
  hRu2->Write();  
  hRv1->Write();  
  hRv2->Write();
    gRu1->SetName("gRu1");
    gRu1->Write();
    gRu2->SetName("gRu2");
    gRu2->Write();
    gRv1->SetName("gRv1");
    gRv1->Write();
    gRv2->SetName("gRv2");
    gRv2->Write();
    gLu1->SetName("gLu1");
    gLu1->Write();
    gLu2->SetName("gLu2");
    gLu2->Write();
    gLv1->SetName("gLv1");
    gLv1->Write();
    gLv2->SetName("gLv2");
    gLv2->Write();  

    tnew->Write();
    
   fnew->Close();
}

/////////////////////////////////////////////////////


void VDCt0::MakeRoot_c(string ofname){

  fnew = new TFile(Form("%s",ofname.c_str()),"recreate");
  tnew =new TTree("T",ofname.c_str());

  //======= Write ========//



  for(int i=0;i<nwire;i++){

    hRu1_time[i]->Write();
    hRu2_time[i]->Write();
    hRv1_time[i]->Write();
    hRv2_time[i]->Write();
    hLu1_time[i]->Write();
    hLu2_time[i]->Write();
    hLv1_time[i]->Write();
    hLv2_time[i]->Write();
    fLu1_t0[i]->SetTitle(Form("%lf * x + %lf",Lu1_p0[i],Lu1_p1[i]));
    fLu1_t0[i]->Write();
    fLu2_t0[i]->SetTitle(Form("%lf * x + %lf",Lu2_p0[i],Lu2_p1[i]));    
    fLu2_t0[i]->Write();
    fLv1_t0[i]->SetTitle(Form("%lf * x + %lf",Lv1_p0[i],Lv1_p1[i]));
    fLv1_t0[i]->Write();
    fLv2_t0[i]->SetTitle(Form("%lf * x + %lf",Lv2_p0[i],Lv2_p1[i]));    
    fLv2_t0[i]->Write();    
    fRu1_t0[i]->Write();
    fRu2_t0[i]->Write();
    fRv1_t0[i]->Write();
    fRv2_t0[i]->Write();
    
  }
  
  hLu1_c->Write();
  hLu2_c->Write();  
  hLv1_c->Write();  
  hLv2_c->Write();
  hRu1_c->Write();
  hRu2_c->Write();  
  hRv1_c->Write();  
  hRv2_c->Write();
    gRu1_c->SetName("gRu1");
    gRu1_c->Write();
    gRu2_c->SetName("gRu2");
    gRu2_c->Write();
    gRv1_c->SetName("gRv1");
    gRv1_c->Write();
    gRv2_c->SetName("gRv2");
    gRv2_c->Write();
    gLu1_c->SetName("gLu1");
    gLu1_c->Write();
    gLu2_c->SetName("gLu2");
    gLu2_c->Write();
    gLv1_c->SetName("gLv1");
    gLv1_c->Write();
    gLv2_c->SetName("gLv2");
    gLv2_c->Write();  

   fnew->Close();
}


/////////////////////////////////////////////////////

#endif
