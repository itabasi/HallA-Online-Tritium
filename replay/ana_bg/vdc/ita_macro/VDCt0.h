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


class VDCt0{

 public:
  VDCt0();
  ~VDCt0();

  
  
 public:
  void SetRunList(string ifname);
  void SetRun(int runnum);
  void SetBranch();
  void MakeHist();
  void Fill();
  double Findt0(bool rarm,char* plane, int wire);
  void Write(string ofname);
  void MakeRoot(string ofname);
  void Draw();
  void Print(string ofname);

  Setting* set;
  //== SetRunList ====//
    int ENum;
    TChain* T;

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

  //===== MakeHist =====//
  TH2F* hRu1;
  TH2F* hRu2;
  TH2F* hRv1;
  TH2F* hRv2;
  TH2F* hLu1;
  TH2F* hLu2;
  TH2F* hLv1;
  TH2F* hLv2;
  TH1D* hRu1_rtime[nwire];
  TH1D* hRu2_rtime[nwire];
  TH1D* hRv1_rtime[nwire];
  TH1D* hRv2_rtime[nwire];
  TH1D* hLu1_rtime[nwire];
  TH1D* hLu2_rtime[nwire];
  TH1D* hLv1_rtime[nwire];
  TH1D* hLv2_rtime[nwire];  

  
  double min_rtime,max_rtime;
  int bin_rtime;


  //==== FIndT0 ======//
  //    double dy,y,yb,ya;
    double a_min,a_max,slope_1,slope_2;
    double t0_min,t0_max,t0_1,t0_2;  
    double dt0;  
  //==== Draw =====//

  TLine* line; 
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
  
  //  TCanvas* c12= new TCanvas("c12","c12");
  //  TCanvas* c13= new TCanvas("c13","c13");
  //  TCanvas* c14= new TCanvas("c14","c14");
  //  TCanvas* c15= new TCanvas("c15","c15");  

  

  //==== Write ====//
  double Ru1t0[nwire],Ru2t0[nwire],Rv1t0[nwire],Rv2t0[nwire];
  double Lu1t0[nwire],Lu2t0[nwire],Lv1t0[nwire],Lv2t0[nwire];
  double Ru1t0_err[nwire],Ru2t0_err[nwire],Rv1t0_err[nwire],Rv2t0_err[nwire];
  double Lu1t0_err[nwire],Lu2t0_err[nwire],Lv1t0_err[nwire],Lv2t0_err[nwire];  

  double Ru1t0_min[nwire],Ru1t0_max[nwire],Ru2t0_min[nwire],Ru2t0_max[nwire],Rv1t0_min[nwire],Rv1t0_max[nwire],Rv2t0_min[nwire],Rv2t0_max[nwire];
  double Lu1t0_min[nwire],Lu1t0_max[nwire],Lu2t0_min[nwire],Lu2t0_max[nwire],Lv1t0_min[nwire],Lv1t0_max[nwire],Lv2t0_min[nwire],Lv2t0_max[nwire];
  

  //==== MakeRoot =====//
  TFile* fnew;
  TTree* tnew;
};


VDCt0::VDCt0(){
  cout<<"start VDC T0 tuning "<<endl;
  set=new Setting();
  set->Initialize();
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
    runnum = runname.substr(35,6);
    num =atoi(runnum.c_str());
    runnum_group = runname.substr(35,5);    
    run_group=atoi(runnum_group.c_str());
    
    T->Add(runname.c_str());
  }
  ENum=T->GetEntries();
  cout<<"Events: "<<ENum<<endl; 
  
  
}

//////////////////////////////////////////////////////////

void VDCt0::SetRun(int runnum){
  T=new TChain("T");
  int sum_run=10;
  cout<<"TChain run number : "<<runnum<<" - "<<runnum+sum_run-1<<endl;

  const string ROOTfilePath="/data/opt_small/VDC_tuning/";
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
	//        cout << "ROOT file " << rootfile << " added to TChain." << endl;
        rootfile = basename + "_" + sub + ".root";
	sub++;
   }
  }//end for

    ENum=T->GetEntries();
  cout<<"Events: "<<ENum<<endl; 

  
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

void VDCt0::MakeHist(){

  
  min_rtime=1500.;
  max_rtime=3000.;
  bin_rtime=1500;
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
  
  for(int i=0;i<nwire;i++){
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
  }

}

///////////////////////////////////////////////////

void VDCt0::Fill(){

  int counts=0;
  bool T1,T4,T5;


  int Ru1_nwire,Ru2_nwire,Rv1_nwire,Rv2_nwire;
  int Lu1_nwire,Lu2_nwire,Lv1_nwire,Lv2_nwire;


  
  ///===== FIll =======//
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

    T1=false;
    T4=false;
    T5=false;
    T->GetEntry(k);
    //    cout<<"evtype: "<<evtype<<endl;
    //    if(evtype==16)for(int i=0;i<nmax;i++)cout<<Form("R.vdc.u1.rawtime[%d]: ",i)<<Ru1_rtime[i]<<endl;
    //    cout<<"R.vdc.u1.rawtime[0]: "<<Ru1_rtime[0]<<endl;
    //    cout<<"R.vdc.u1.rawtime[1]: "<<Ru1_rtime[1]<<endl;
    //    cout<<"R.vdc.u1.rawtime[2]: "<<Ru1_rtime[2]<<endl;
    //    cout<<"R.vdc.u1.rawtime[3]: "<<Ru1_rtime[3]<<endl;
    
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

    

    if(T5){
      /*
    cout<<"Light VDC : "<<k<<"  T1: "<<T1<<endl;
    cout<<"U1 wire: "<<Lu1_wire[0]<<"  U1 Nhits :"<<NLu1_wire<<endl;
    cout<<"U2 wire: "<<Lu2_wire[0]<<"  U2 Nhits :"<<NLu2_wire<<endl;
    cout<<"V1 wire: "<<Lv1_wire[0]<<"  V1 Nhits :"<<NLv1_wire<<endl;
    cout<<"V2 wire: "<<Lv2_wire[0]<<"  V2 Nhits :"<<NLv2_wire<<endl;    
      */
      
    if(NLu1_wire>0 && 100>NLu1_wire)for(int i=0;i<NLu1_wire;i++){hLu1->Fill(Lu1_wire[i],Lu1_rtime[i]);}
    if(NLu2_wire>0 && 100>NLu2_wire)for(int i=0;i<NLu2_wire;i++){hLu2->Fill(Lu2_wire[i],Lu2_rtime[i]);}
    if(NLv1_wire>0 && 100>NLv1_wire)for(int i=0;i<NLv1_wire;i++){hLv1->Fill(Lv1_wire[i],Lv1_rtime[i]);}
    if(NLv2_wire>0 && 100>NLv2_wire)for(int i=0;i<NLv2_wire;i++){hLv2->Fill(Lv2_wire[i],Lv2_rtime[i]);}

    }


  
    if(T4){


      /*
    cout<<"Right VDC : "<<k<<endl;
    cout<<"U1 wire: "<<Ru1_wire[0]<<"  U1 Nhits :"<<NRu1_wire<<" TDC : "<<Ru1_rtime[0]<< endl;
    cout<<"U2 wire: "<<Ru2_wire[0]<<"  U2 Nhits :"<<NRu2_wire<<" TDC : "<<Ru2_rtime[0]<< endl;
    cout<<"V1 wire: "<<Rv1_wire[0]<<"  V1 Nhits :"<<NRv1_wire<<" TDC : "<<Rv1_rtime[0]<< endl;
    cout<<"V2 wire: "<<Rv2_wire[0]<<"  V2 Nhits :"<<NRv2_wire<<" TDC : "<<Rv2_rtime[0]<< endl;    
      */

    if(NRu1_wire>0 && 100>NRu1_wire)for(int i=0;i<NRu1_wire;i++){hRu1->Fill(Ru1_wire[i],Ru1_rtime[i]);}
    if(NRu2_wire>0 && 100>NRu2_wire)for(int i=0;i<NRu2_wire;i++){hRu2->Fill(Ru2_wire[i],Ru2_rtime[i]);}
    if(NRv1_wire>0 && 100>NRv1_wire)for(int i=0;i<NRv1_wire;i++){hRv1->Fill(Rv1_wire[i],Rv1_rtime[i]);}
    if(NRv2_wire>0 && 100>NRv2_wire)for(int i=0;i<NRv2_wire;i++){hRv2->Fill(Rv2_wire[i],Rv2_rtime[i]);}
    
    }



    if(k%(ENum/10)==0){
      int fill=k/(ENum/10)*10;
      cout<<"Filled : "<<fill<<" % ("<<k<<" / "<<ENum<<")"<<endl;}
    
  }//end Fill

 ///===== Project =======//

  for(int i=0;i<nwire;i++){
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
  }

  
  cout<<"Have Done !! "<<endl;
  
  }


/////////////////////////////////////////////////////


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
    Double_t dx = hnew->GetBinWidth(1);
    Double_t slope_min = 0.0, slope = 0.0;
    Double_t sdbin = -1; // Bin number with steepest descent
    double slope_0=0.0;
    Int_t nbins = hnew->GetNbinsX();
    Int_t maxbin = hnew->GetMaximumBin();
    Double_t maxcont = hnew->GetBinContent(maxbin);

    for (Int_t i=maxbin+1; i<=nbins; i++) {
        if (hnew->GetBinContent(i)>0.5*maxcont) continue;
	
        if (i>1 && i<nbins) {
            slope = (hnew->GetBinContent(i+1)-hnew->GetBinContent(i-1))/(2.*dx);

        } else if (i==1) {
            slope = (hnew->GetBinContent(i+1)-hnew->GetBinContent(i))/dx;
        } else if (i==nbins) {
            slope = (hnew->GetBinContent(i)-hnew->GetBinContent(i-1))/dx;
        }

        if (TMath::Abs(slope)<1.e-4) slope = 0.;


	
        if (slope<slope_min) {
            slope_min = slope;
            sdbin = i;
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

	if(ya==0 || yb==0 || slope_min==0 || y==0 || dt0<-1.0)dt0=100.0;
	else{dt0=(t0-x)*(sqrt(y)/y+(sqrt(yb)+sqrt(ya))/fabs(ya-yb));	}
	//       	cout<<"dt0 : "<<dt0 <<" t0 "<<t0 << " y "<<y<<" yb "<<yb <<" ya "<<ya<<" slope_min "<<slope_min<<" slope "<<slope<<endl;


	/*
	slope_0 = (ya-yb)/(2.*dx);
	slope_1 = ((ya-sqrt(ya))-(yb+sqrt(yb)))/2.*dx; 
	slope_2 = ((ya+sqrt(ya))-(yb-sqrt(yb)))/2.*dx;
	t0_1 = x-sqrt(x) -y/(slope_1); // minimum 
	t0_2 = x+sqrt(x) -y/(slope_2);  // maximum

	cout<<endl;
	cout<<"ya: "<<ya <<" y "<<y <<" yb "<<yb<<endl;
	cout<<"slope_1 "<<slope_1<<" slope "<<slope_0<<" slope_2 "<<slope_2<<endl;
	cout<<"t0_1: "<<t0_1<<" t0: "<<t0<<" t0_2 "<<t0_2<<endl;
	*/

	//	if(t0_1>t0_2){t0_max=t0_1; t0_min=t0_2;}
	//	else {t0_min=t0_1; t0_max=t0_2;}

    }

    
    return t0;
}


//////////////////////////////////////////////////////

void VDCt0::Write(string ofname){



  string ofname_main = ofname.substr(0,21);
  string ofname_err =ofname_main + "_err.dat";
  //  string ofname_min  = ofname_main + "_min.dat";
  //  string ofname_max  = ofname_main + "_max.dat";

  ofstream ofs(ofname.c_str());
  ofstream ofs_err(ofname_err.c_str());
  
  //  ofstream ofs_min(ofname_min.c_str());
  //  ofstream ofs_max(ofname_max.c_str());  
  
  double T0,T0_min,T0_max;
  double dT0;
  string vdc_name;

  
  //==== Fill ====//
  
    for(int i=0;i<nwire;i++){
     Ru1t0[i]=Findt0(true,"U1",i);
     //     Ru1t0_min[i]=t0_min; Ru1t0_max[i]=t0_max;
     Ru1t0_err[i]=dt0;
     Ru2t0[i]=Findt0(true,"U2",i);
     Ru2t0_err[i]=dt0;     
     //     Ru2t0_min[i]=t0_min; Ru2t0_max[i]=t0_max;     
     Rv1t0[i]=Findt0(true,"V1",i);
     Rv1t0_err[i]=dt0;     
     //     Rv1t0_min[i]=t0_min; Rv1t0_max[i]=t0_max;     
     Rv2t0[i]=Findt0(true,"V2",i);
     Rv2t0_err[i]=dt0;     
     //     Rv2t0_min[i]=t0_min; Rv2t0_max[i]=t0_max;     
     Lu1t0[i]=Findt0(false,"U1",i);
     Lu1t0_err[i]=dt0;     
     //     Lu1t0_min[i]=t0_min; Lu1t0_max[i]=t0_max;
     Lu2t0[i]=Findt0(false,"U2",i);
     Lu2t0_err[i]=dt0;     
     //     Lu2t0_min[i]=t0_min; Lu2t0_max[i]=t0_max;     
     Lv1t0[i]=Findt0(false,"V1",i);
     Lv1t0_err[i]=dt0;          
     //     Lv1t0_min[i]=t0_min; Lv1t0_max[i]=t0_max;          
     Lv2t0[i]=Findt0(false,"V2",i);
     Lv2t0_err[i]=dt0;     
     //     Lv2t0_min[i]=t0_min; Lv2t0_max[i]=t0_max;

    }
  
  
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
    //    ofs_min << vdc_name <<"T0 min error "<<endl;
    //    ofs_max << vdc_name <<"T0 max error "<<endl;
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
      /*
      if(j==0){T0= Ru1t0[i]; T0_min=Ru1t0_min[i];  T0_max=Ru1t0_max[i];}
      else if(j==1){T0= Ru2t0[i]; T0_min=Ru2t0_min[i];  T0_max=Ru2t0_max[i];}
      else if(j==2){T0= Rv1t0[i]; T0_min=Rv1t0_min[i];  T0_max=Rv1t0_max[i];}
      else if(j==3){T0= Rv2t0[i]; T0_min=Rv2t0_min[i];  T0_max=Rv2t0_max[i];}
      else if(j==4){T0= Lu1t0[i]; T0_min=Lu1t0_min[i];  T0_max=Lu1t0_max[i];}
      else if(j==5){T0= Lu2t0[i]; T0_min=Lu2t0_min[i];  T0_max=Lu2t0_max[i];}
      else if(j==6){T0= Lv1t0[i]; T0_min=Lv1t0_min[i];  T0_max=Lv1t0_max[i];}
      else if(j==7){T0= Lv2t0[i]; T0_min=Lv2t0_min[i];  T0_max=Lv2t0_max[i];}
      else {cout<<"faled to write"<<endl; break;}      
      */
      
      if(fabs(T0)<10000) ofs << T0 <<" ";
      else ofs << 0.0 <<" ";
      if(fabs(T0_min)<1000)  ofs_err << dT0 <<" ";    
      else ofs_err << 0.0 <<" ";          
      //      if(fabs(T0_min)<10000)  ofs_min << T0_min <<" ";    
      //      else ofs_min << 0.0 <<" ";    
      //      if(fabs(T0_max)<10000)   ofs_max << T0_max <<" ";
      //      else  ofs_max << 0.0 <<" ";    
   
    if((i+1)%8==0){
      ofs << endl;
      ofs_err << endl;}
      //      ofs_min << endl;
      //      ofs_max << endl;}
    
    }
  }
  

  
  ofs.close();
  ofs_err.close();
  //  ofs_min.close();
  //  ofs_max.close();  
  
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
    line= new TLine(Lu1t0[32*j+i],0,Lu1t0[32*j+i],hLu1_rtime[32*j+i]->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw("same");        
    c1[j]->cd(i+1);
    hLu2_rtime[32*j+i]->Draw();
    line= new TLine(Lu2t0[32*j+i],0,Lu2t0[32*j+i],hLu2_rtime[32*j+i]->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw("same");        
    c2[j]->cd(i+1);
    hLv1_rtime[32*j+i]->Draw();
    line= new TLine(Lv1t0[32*j+i],0,Lv1t0[32*j+i],hLv1_rtime[32*j+i]->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw("same");    
    c3[j]->cd(i+1);
    hLv2_rtime[32*j+i]->Draw();
    line= new TLine(Lv2t0[32*j+i],0,Lv2t0[32*j+i],hLv2_rtime[32*j+i]->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw("same");        
    

    }
  }
  




  
 
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
    if(i==10 && j==3)c3[i]->Print(Form("%s]",ofname.c_str()));
  }
  }
  cout<<"Print is done !"<<endl;

}

/////////////////////////////////////////////////////

void VDCt0::MakeRoot(string ofname){

  fnew = new TFile(Form("%s",ofname.c_str()),"recreate");
  tnew =new TTree("T",ofname.c_str());

  //======= Write ========//

  hLu1->Write();
  hLu2->Write();  
  hLv1->Write();  
  hLv2->Write();
  hRu1->Write();
  hRu2->Write();  
  hRv1->Write();  
  hRv2->Write();  

  for(int i=0;i<nwire;i++){

    hRu1_rtime[i]->Write();
    hRu2_rtime[i]->Write();
    hRv1_rtime[i]->Write();
    hRv2_rtime[i]->Write();
    hLu1_rtime[i]->Write();
    hLu2_rtime[i]->Write();
    hLv1_rtime[i]->Write();
    hLv2_rtime[i]->Write();
    
  }
   fnew->Close();
}


/////////////////////////////////////////////////////

#endif
