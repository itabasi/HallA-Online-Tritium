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
  cout<<"Pion Flag:"<<Pion_flag<<endl;
 
 

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



 //==================================================//  
 // Time scale [ns] //
 // Energy scale [GeV] //

 double ac1_adc[nth],ac2_adc[nth];
 double min_coin,max_coin,min_coin_c,max_coin_c;
 double min_ac1,max_ac1,min_ac2,max_ac2,min_adc,max_adc;
 double th1_max,th2_max;
 double ac1_kcut,ac2_kcut_min,ac2_kcut_max;
 if(mode=="H" || mode=="T"){
 min_coin=-10;
 max_coin=20.0;
 min_coin_c=-10;
 max_coin_c=20.0;
 //testing
 min_coin_c=0.0;
 max_coin_c=10.0;
 min_ac1=0.0;
 max_ac1=5000.;
 min_ac2=0.0;
 max_ac2=20000.;
 min_adc=-500.0;
 max_adc=20000.;
//=== AC Threshold variable ===//
 th1_max=500.;
 ac1_adc[0]=th1_max;
 ac1_adc[1]=300.;
 ac1_adc[2]=100.;

 ac2_adc[0]=min_ac2;
 ac2_adc[1]=1000.;
 ac2_adc[2]=3000.;
 th2_max=6000.;

 //---Kaon Cut ----//
 ac1_kcut=100.;
 ac2_kcut_min=1000.;
 ac2_kcut_max=5000.;
 //----------------//
 }else if(mode=="G"){
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

}


 double bin_coin=(max_coin-min_coin)/tdc_time;
        bin_coin=(int)bin_coin;
 double bin_coin_c=(max_coin_c-min_coin_c)/tdc_time;
        bin_coin_c=(int)bin_coin_c;
 int bin_beta=6000;
 int bin_adc=max_adc-min_adc;
 int bin_ac1=(max_ac1-min_ac1)*3; 
 int bin_ac2=(max_ac2-min_ac2)*3; 

 TH1F* hcoin[1000];
 TF1* fpi[1000];
 TGraph* fcoin=new TGraph();
 int j=0; int l=0;
 double mean[1000];



 
 int bin_Rs2adc,min_Rs2adc,max_Rs2adc,bin_Ls2adc,min_Ls2adc,max_Ls2adc; 
 min_Rs2adc=5000;
 max_Rs2adc=6000;
 min_Ls2adc=5000;
 max_Ls2adc=6000;
 bin_Rs2adc=(max_Rs2adc-min_Rs2adc);
 bin_Ls2adc=(max_Ls2adc-min_Ls2adc);
 
 
 double mtr;
 mtr=938.27e-3;// proton mass [GeV/c^2]
 // double mh;
 double m2; 
 double Ee,Ee_,Ek,Epi;
 double pe,pe_,pk,ppi;
 double coin_t,coin_tc;
 double rtof[16];
 double rbeta,rbeta_k,lbeta;
 double Rs2_off,Ls2_off; 
 double Rs2_tcorr,Ls2_tcorr;
 int Ls2pads,Rs2pads;
 bool cut_ac1,cut_ac2,cut_beta;
 int nac1,nac2,nac3,n;
 double tof_r,tof_l; 
 double rpathl,lpathl;
 double corr_R,corr_L;
 double rpath_corr,lpath_corr;


 //--- Coin Offset -----//
 double pathl_off,s2_offset,coin_offset;
 pathl_off=0.0;

 pathl_off=0.0;
 if(mode=="H"&&kine==2){
 coin_offset=498.;
 s2_offset=-489.0;
 pathl_off=-485.5; 

 }else if(mode=="H" && kine==1){
 pathl_off=-498.+30.-3.0+0.5;
 s2_offset=-500.0+25.;
 coin_offset=-41.35+498.;
 if(Pion_flag==1){coin_offset=coin_offset;cout<<"Pion flag"<<Pion_flag<<endl;}else{coin_offset=coin_offset;cout<<"Pion flag"<<Pion_flag<<endl;}

}


 // double mm;
 //----- Cut Parameters ----------//
 double coin_cutmin=-248;
 double coin_cutmax=-244; 
 double rpathl_cutmin=28.7;
 double rpathl_cutmax=29.4;
 double lpathl_cutmin=28.6;
 double lpathl_cutmax=29.2;
 double rbeta_cutmin=0.0;
 double rbeta_cutmax=1.0;
 double lbeta_cutmin=0.9;
 double lbeta_cutmax=1.0;
 double Rvz_cutmin=-0.1;
 double Rvz_cutmax= 0.1;
 double Lvz_cutmin=-0.1;
 double Lvz_cutmax= 0.1;
 double Rx_cutmin= -0.4;
 double Rx_cutmax= 0.4;
 //-------------------------------//
 bool cut_Rs2,cut_Ls2,cut_rpathl,cut_lpathl,cut_coin,cut_rbeta,cut_lbeta,cut_vz,cut_Rx,cut_trig,coin_trig,right_trig,cut_track,cut_s0;


 double top_x[1000],top_y[1000];
 char* r="r";

    TCanvas* can[100];
    int mv=0;
  //===========================================================//
  //============== RUN STUDY ==================================//
  //===========================================================//

 int runnum,sub_num,i;
  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname, run_num,sub,run_i;

  while(1){

     //TChain //
  TChain* T;
  if(mode=="G"|| mode=="T"){T=new TChain("tree"); }
  else {T=new TChain("T");}

  
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> runname;
    run_num=buf.substr(28,6);
    runnum=atoi(run_num.c_str());
    sub=buf.substr(35,1);
    runnum=atoi(run_num.c_str());
    run_i=buf.substr(31,3);
    i=atoi(run_i.c_str());
    //cout<<"run num: "<<runnum<<endl;
    //cout<<"sub num: "<<sub<<endl;

    if(sub!=r){T->Add(runname.c_str());}
    else if(sub==r){
           T->Add(runname.c_str());
	   cout<<"run num: "<<runnum<<endl;
    int evnt=T->GetEntries();
    if(test_flag){evnt=100000;}
    cout<<"Get Entries: "<<evnt<<endl;

  T->SetBranchStatus("*",0);  
 //------ Right Arm -------------//

 T->SetBranchStatus("RTDC.F1FirstHit",1);
 T->SetBranchAddress("RTDC.F1FirstHit",RF1); 
 T->SetBranchStatus("R.s2.t_pads",1);
 T->SetBranchAddress("R.s2.t_pads",Rs2tpads);
 T->SetBranchStatus("R.s2.trpad",1);
 T->SetBranchAddress("R.s2.trpad",Rs2trpad);
 T->SetBranchStatus("R.s2.ra",1);
 T->SetBranchAddress("R.s2.ra",Rs2ra);
 T->SetBranchStatus("R.s2.la",1);
 T->SetBranchAddress("R.s2.la",Rs2la);
 T->SetBranchStatus("R.a1.asum_c",1);
 T->SetBranchAddress("R.a1.asum_c",&Ra1sum);
 T->SetBranchStatus("R.a2.asum_c",1);
 T->SetBranchAddress("R.a2.asum_c",&Ra2sum);

  // path length//
 T->SetBranchStatus("R.s2.trpath",1); 
 T->SetBranchAddress("R.s2.trpath",rs2pathl); 
 T->SetBranchStatus("R.tr.pathl",1);  
 T->SetBranchAddress("R.tr.pathl",rtrpathl);
 // Target positon information //
 T->SetBranchStatus("R.tr.p",1);
 T->SetBranchAddress("R.tr.p",Rp);
 T->SetBranchStatus("R.tr.vz",1);    
 T->SetBranchAddress("R.tr.vz",Rvz); 

 //------ Left Arm ---------------//
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

 if(mode=="G"){
 T->SetBranchStatus("ctime",1);    
 T->SetBranchAddress("ctime",ctime);
 T->SetBranchStatus("DR.T5",1);    
 T->SetBranchAddress("DR.T5",&DRT5);

}



 //========= Definition of Hist ========================// 
  hcoin[i]=new TH1F(Form("hcoin[%d]",i),"",bin_coin_c,min_coin_c,max_coin_c);
  set->SetTH1(hcoin[i],"Coincidence time w/ Pion cut","coin time [ns]","Counts");


 double Ltof;
 int ev=0;

 
 for(int k=0;k<evnt;k++){
   T->GetEntry(k);

   /*
   if(k==ev*100000){
     cout<<"Fill Event: "<<k<<"/"<<evnt<<endl; 
     ev=ev+1;

     }*/


 pe_=Lp[0];//*sqrt(1+pow(Lth[0],2)+pow(Lph[0],2));
 pk=Rp[0];//*sqrt(1+pow(Rth[0],2)+pow(Rph[0],2));
 ppi=Rp[0];//*sqrt(1+pow(Rth[0],2)+pow(Rph[0],2));
 pe=4.313; // GeV   //hallap*1.0e-3;
 Ee=sqrt(pow(pe,2)+pow(me,2));
 Ee_=sqrt(pow(pe_,2)+pow(me,2));
 Epi=sqrt(pow(ppi,2)+pow(mpi,2));
 Ek=sqrt(pow(pk,2)+pow(mk,2));
 Ls2pads=Ls2tpads[0];
 Rs2pads=Rs2tpads[0];
 rpathl=rtrpathl[0]+rs2pathl[0]; // R-HRS path length S2 -RF
 lpathl=ltrpathl[0]+ls2pathl[0]; // L-HRS path length S2 -RF
 rpath_corr=rpathl/rbeta/c;

 lbeta=1.0;//pe_/Ee_; 
 lpath_corr=lpathl/lbeta/c;
 Rs2_off=s2f1_off(Rs2pads,"R",mode,kine);
 Ls2_off=s2f1_off(Ls2pads,"L",mode,kine);
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 Ltof=tof_l-(-LF1[27]+LF1[28])/2.0*tdc_time;
 if(mode=="G"){
   coin_t=ctime[0];
   coin_tc=ctime[0];
 }else{

  coin_t=tof_r-tof_l-coin_offset-s2_offset; //coin time
  coin_tc=tof_r-tof_l+rpath_corr-lpath_corr-pathl_off-coin_offset; //  coin Path & Offset  correction 

}
//====== Cut condition ========================// 
   cut_rpathl=false;
   cut_lpathl=false;
   cut_Rs2=false;
   cut_Ls2=false;
   cut_vz=false;
   cut_track=false;
   cut_s0=false;
   coin_trig=false;
   if(Rvz_cutmin<Rvz[0] && Rvz[0]<Rvz_cutmax && Lvz_cutmin<Lvz[0] && Lvz[0]<Lvz_cutmax)cut_vz=true;
  
   if(mode=="G"){
     if(ctime[0]>-100)coin_trig=true;
     cut_track=true;
     cut_rpathl=true;
     cut_lpathl=true;
   }else{
 
  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads)cut_track=true;
  if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax)cut_rpathl=true;
  if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax)cut_lpathl=true;
  if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0 && -LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0)cut_s0=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;
  if(cut_Rs2 && cut_Ls2 && cut_s0)coin_trig=true;  
}
 //=======================================//
  


   //==========================================//
   //========= Fill Hist =====================//
   //========================================//


   //--------- Pion Cut Hist ---------------//

   //   if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track && Ra1sum>ac1_kcut 
   // && Ra2sum>ac2_kcut_max){hcoin[i]->Fill(coin_tc);}
   if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track){hcoin[i]->Fill(coin_tc);}

 

 }//end fill

 
 fpi[i]=new TF1(Form("fpi[%d]",i),"gaus(0)",0.0,10.0);
 top_x[i]=hcoin[i]->GetBinContent(hcoin[i]->GetMaximumBin());
 top_y[i]=hcoin[i]->GetMaximum();
 fpi[i]->SetParameter(0,top_y[i]);
 fpi[i]->SetParameter(1,top_x[i]);	    
 hcoin[i]->Fit(Form("fpi[%d]",i),"","",top_x[i]-2.0,top_x[i]+2.0);
 mean[i]=fpi[i]->GetParameter(1);
 cout<<"pion position [ns]: "<<mean[i]<<endl;
 cout<<"l: "<<l<<endl;

 

 cout<<"j: "<<j<<endl;
 if(l==12*j){
   can[j]=new TCanvas(Form("can[%d]",j),Form("can[%d]",j));
   can[j]->Divide(4,3);
   can[j]->cd(1);
   mv=2;
   hcoin[i]->Draw();}
 else{
   cout<<"mv: "<<mv<<endl;
   can[j]->cd(mv);
   hcoin[i]->Draw();
   mv=mv+1;}
 //if()
 
 l=l+1; 
    }//end TChain
    
    delete T;
  }// end while


 TString name;
 if(output_flag){
 name.Form(ofname.c_str());
 can[0]->Print(name+"[","pdf");
 can[0]->Print(name,"pdf");
 for(int k=1;k<j;k++){
 can[k]->Print(name,"pdf");}
 can[j]->Print(name,"pdf");
 can[j]->Print(name +"]","pdf");
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
