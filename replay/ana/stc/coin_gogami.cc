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

void coin_gogami(){

  TFile *f1=new TFile("../rootfiles/coin_H2.root","read");
  TTree *T =(TTree*)f1->Get("tree");
 //============= Set Branch Status ==================//

  int max=100; 
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
  double DRevtype;
  double Rbeta[max],Lbeta[max];
  double rs2pathl[max],rs0pathl[max],rtrpathl[max];
  double ls2pathl[max],ls0pathl[max],ltrpathl[max];
  double trigger[100];
  double hallap;
  double Rs2tpads[100],Ls2tpads[100];
  double Rs2trpad[100],Ls2trpad[100];
  double Ls2rt_c[100],Ls2lt_c[100];
  //---- Gogami root ---------//
  double ctime[1000];
  double DRT5;
  double Rs2ra[100];
  double Rs2la[100];
  double Ls2ra[100];
  double Ls2la[100];
  double Ls2la_c[100];
  double Ls2ra_c[100];


 T->SetBranchAddress("ctime",ctime);
 T->SetBranchAddress("LTDC.F1FirstHit",LF1); 
 T->SetBranchAddress("L.s2.t_pads",Ls2tpads);
 T->SetBranchAddress("L.s2.trpad",Ls2trpad);
 
 int nth=10;
 double ac1_adc[nth],ac2_adc[nth];
 double min_coin,max_coin,min_coin_c,max_coin_c;
 double min_ac1,max_ac1,min_ac2,max_ac2,min_adc,max_adc;
 double th1_max,th2_max;
 double ac1_kcut,ac2_kcut_min,ac2_kcut_max;


 min_coin=-20;
 max_coin=20.0;
 min_coin_c=-20;
 max_coin_c=20.0;
 min_ac1=0.0;
 max_ac1=30.;
 min_ac2=0.0;
 max_ac2=50.;
 min_adc=-5.0;
 max_adc=20.;
//==== AC Threshold variable (ACTH)===//
 ac1_adc[0]=max_ac1;
 ac1_adc[1]=1.3;
 ac1_adc[2]=1.0;
 ac2_adc[0]=min_ac2;
 ac2_adc[1]=5.;
 ac2_adc[2]=7;
 

 int nev=T->GetEntries();
 cout<<"Get Entry: "<<nev<<endl;


 
 //===== FIll =======//
 TH2F* hcoin_event=new TH2F("hcoin_event","hist",1000,0,nev,1000,-5,5);
 TH2F* hLs2r_event[16];
 TH2F* hLs2l_event[16]; 
 double tdc_time=58.0e-3;
 double min_tdc=-750.0;
   //-13000*tdc_time;
 double max_tdc=-700.0;
   //-12000*tdc_time;
 for(int i=0;i<16;i++){
   hLs2r_event[i]=new TH2F(Form("hLs2r_event[%d]",i),Form("HRS-L S2 %d Seg  R-PMT Events",i),1000,0,nev,1000,min_tdc,max_tdc);
   hLs2l_event[i]=new TH2F(Form("hLs2l_event[%d]",i),Form("HRS-L S2 %d Seg  L-PMT Events",i),1000,0,nev,1000,min_tdc,max_tdc);
 }
 
 double ct;
 int Ls2pads;
 double f1_s2r;
 double f1_s2l;
 for(int k=0;k<nev;k++){
   T->GetEntry(k);
   ct=ctime[0];

   hcoin_event->Fill(k,ct);
   for(int j=0;j<16;j++){
     f1_s2r=(+LF1[j]-LF1[30])*tdc_time;
     f1_s2l=(+LF1[j+48]-LF1[37])*tdc_time;
     hLs2r_event[j]->Fill(k,f1_s2r);
     hLs2l_event[j]->Fill(k,f1_s2l);
   }

 }

 
 TCanvas* c0=new TCanvas("c0","c0");
 c0->cd();
 hcoin_event->Draw("colz");

 TCanvas* c1=new TCanvas("c1","c1");
 c1->Divide(2,2);
 c1->cd(1);
 hLs2r_event[7]->Draw("colz");
 c1->cd(2);
 hLs2r_event[12]->Draw("colz");
 c1->cd(3);
 hLs2r_event[13]->Draw("colz");
 c1->cd(4);
 hLs2r_event[14]->Draw("colz"); 
}
