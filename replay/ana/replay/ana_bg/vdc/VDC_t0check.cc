#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

void VDC_t0check(){

  
  int nmax=15;

  string buf,runname;
  char* test;
  int wmax=400;
  int nplane=4;// U1, V1, U2, V2
  int NRUN=nmax;
  double t0[wmax][nplane][NRUN];
  int RUN=0;
  int run_init=111150;
  char* arm;
   arm="R";
  //  arm="L";
   //   bool T5=true;
  
  for(int i=0;i<nmax;i++){
  int run;
  
   run=run_init+i*50;
   RUN=i;
   //   RUN=run-111000;// last three digits
   /*
   //====coin tirgger in RHRS ====//
   ifstream ifs(Form("DB/db_%s_T5_.vdc.%d.dat",arm,run));
   runname=Form("DB/db_%s_T5_.vdc.%d.dat",arm,run);
   */

   ifs(Form("DB/db_%s.vdc.%d.dat",arm,run));
   runname=Form("DB/db_%s.vdc.%d.dat",arm,run);}

  if(ifs.fail()){ 
    cout<<Form("Can not open file :%s",runname.c_str())<<endl; continue;}
  else cout<<Form("Open file : %s",runname.c_str())<<endl;


  int k=-1; // Wire plane U1, U2, V1 ,V2
  int nwire=0; // # of wire
  string T0;


  while(!ifs.eof()){

    ifs>>T0;
    if(T0[0]=='R' || T0[0]=='L'){k=k+1;nwire=0; continue;}
    t0[nwire][k][RUN]=stod(T0);
    t0[nwire][k][RUN]=t0[nwire][k][RUN]*0.5;//convert bin to ns 0.5 [ns/bin]
    //    cout<<Form("t0[%d][%d][%d] : ",nwire,k,RUN)<<t0[nwire][k][RUN]<<endl;
    nwire++;
  }
 
  }//end for



  TGraph* gu1=new TGraph();
  TGraph* gu2=new TGraph();  
  TGraph* gv1=new TGraph();
  TGraph* gv2=new TGraph();

  int wire_num=300;

  for(int i=0;i<nmax;i++){

    int run=150;
   run=run+i*50;
   gu1->SetPoint(i,run,t0[wire_num][0][i]);
   gv1->SetPoint(i,run,t0[wire_num][1][i]);   
   gu2->SetPoint(i,run,t0[wire_num][2][i]);
   gv2->SetPoint(i,run,t0[wire_num][3][i]);   
  }


  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  gStyle->SetTitleYOffset(1.5);
  gu1->SetMarkerStyle(21);
  gu1->SetMarkerColor(1);
  gu1->SetMarkerSize(1);
  gv1->SetMarkerStyle(21);
  gv1->SetMarkerColor(kGreen);
  gv1->SetMarkerSize(1);
  gu2->SetMarkerStyle(21);
  gu2->SetMarkerColor(kBlue);
  gu2->SetMarkerSize(1);
  gv2->SetMarkerStyle(21);
  gv2->SetMarkerColor(kRed);
  gv2->SetMarkerSize(1);
  if(arm=="R")gu1->GetYaxis()->SetRangeUser(2640*0.5,2660*0.5);
  else gu1->GetYaxis()->SetRangeUser(5880/4.,5920/4.);
  if(arm=="R" && T5)gu1->GetYaxis()->SetRangeUser(1300,1400);
  gu1->SetTitle(Form("%s-HRS VDC Offset run dependence;Run number ;TDC [ns]",arm));
  
  gu1->Draw("AP");
  gv1->Draw("P");
  gu2->Draw("P");
  gv2->Draw("P");
  
  
}
