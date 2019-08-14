#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
using namespace std;

const int nmax=100;
void vdc_vel(){

  string filename="/data/nnL_smallroot/tritium_111180.root";
  TChain* T=new TChain("T");
  T->Add(filename.c_str());

 T->SetBranchStatus("*",0);

  //==========================//
  //=========== VDC ==========//
  //==========================//

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
  double Rtr_th[100];
  double Ru1_dist[nmax]; 
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

  T->SetBranchStatus("R.tr.th"          ,1);         T->SetBranchAddress("R.tr.th"          ,Rtr_th);         


  T->SetBranchStatus("Ndata.R.vdc.u1.time"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.u2.time"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.v1.time"          ,1);
  T->SetBranchStatus("Ndata.R.vdc.v2.time"          ,1);  
  T->SetBranchStatus("R.vdc.u1.dist"          ,1);  T->SetBranchAddress("R.vdc.u1.dist"          ,Ru1_dist);         
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

 

  //  TF1* fRu1_vel=new TF1("[0] + [1]/x + [2]/pow(x,2) + [3]/pow(x,3)",1.0,1.5);
  TH2F* hRu1_vel_tan=new TH2F("hRu1_vel_tan","",1000,0.0,2.0,1000,-0.2,0.2);
  TH2F* hRu1_vel_tan_c=new TH2F("hRu1_vel_tan_c","",1000,0.0,2.0,1000,-10,10);

  const  double a1=2.12e-3;
  double a2=0.0;
  double a2_p[4]={-4.2e-4, 1.3e-3, 1.06e-4, 0.0};
  double val=0.0;
  double vel =49280.; //velocity [m/s]
  double val_cor=0.0;
  int ENum=T->GetEntries();
  cout<<"Entries : "<<ENum<<endl;
  for(int i=0;i<ENum;i++){

    T->GetEntry(i);
    
    double x=1./tan(Rtr_th[0]);
    double xx=tan(Rtr_th[0]);
    a2=a2_p[0] + a2_p[1]*x + a2_p[2]*x*x + a2_p[3]*x*x*x;
    val_cor=(1.+a2/a1);


    for(int j=0;j<NRu1_wire;j++){
      if(xx!=0. && Ru1_time[j]>0){
	val=Ru1_dist[j]/Ru1_time[j]/vel;
	hRu1_vel_tan->Fill(val,xx);
	hRu1_vel_tan_c->Fill(val*val_cor,xx);}
    }
  }
  

  TCanvas* c0=new TCanvas("c0","c0");
  c0->Divide(2,1);
  c0->cd(1);
  hRu1_vel_tan->Draw("colz");
  c0->cd(2);
  hRu1_vel_tan_c->Draw("colz");    
  
}
