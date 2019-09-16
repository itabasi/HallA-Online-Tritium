// Author K. Itabashi Aug. 30th
// HRS nnL experiment missing mass analysis
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
//==================================================================//

int main(int argc, char** argv){

  gStyle->SetOptFit(111111111);
  int ch; char* mode;
  int kine=1;// 1: hydrogen kinematics 2:tritium kinematics
  double tdc_time=58.0e-3;//[ns]
  string ifname = "/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/scripts/ita_scripts/run_list/Lambda_test.list";
  string ofname = "/pdf/hydro1_AC_eff_test.pdf";
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  bool test_flag=false;
  bool pra_flag=false;
  bool ac2_min=true;
  //bool ac2_min=false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"h:f:w:n:bcop:GHTtPu12"))!=-1){
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

    case 'P':
      pra_flag=true;    
	break;

    case 't':
      test_flag=true;    
	break;

    case 'u':
      ac2_min=false;    
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


 //TChain //


  TChain* T;
  if(mode=="G"||mode=="T" ){T=new TChain("tree"); }
  else {T=new TChain("T");}


  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> runname;
    T->Add(runname.c_str());
    cout<<buf<<endl;
  }

  cout<<"mode :"<<mode<<endl;
  cout<<"tdc_time[ns]: "<<tdc_time<<endl;
  int evnt=T->GetEntries();
  if(test_flag){evnt=10000;}
  cout<<"Get Entries: "<<evnt<<endl;
  
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
  double Rth[max],Rph[max],Rx[max],Rvz[max],Lth[max],Lph[max],Lx[max],Lvz[100];
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
  //---- Toyama ana_Lambda ----//
  double Rz,Lz,Rpz,Lpz; 
  double tcoin_t;
  //------------------------//
 double mm; 
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
 double ct_acc; 
 double acc;


 T->SetBranchStatus("*",0);  
 
//------ Right Arm -------------//

 T->SetBranchStatus("RTDC.F1FirstHit",1);
 T->SetBranchAddress("RTDC.F1FirstHit",RF1); 
 T->SetBranchStatus("R.s2.t_pads",1);
 T->SetBranchAddress("R.s2.t_pads",Rs2tpads);
 T->SetBranchStatus("R.s2.trpad",1);
 T->SetBranchAddress("R.s2.trpad",Rs2trpad);
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
  // path length//
 T->SetBranchStatus("L.s2.trpath",1); 
 T->SetBranchAddress("L.s2.trpath",ls2pathl); 
 T->SetBranchStatus("L.tr.pathl",1);   
 T->SetBranchAddress("L.tr.pathl",ltrpathl);
 T->SetBranchStatus("L.tr.p",1);
 T->SetBranchAddress("L.tr.p",Lp);  
 T->SetBranchStatus("L.tr.vz",1);    
 T->SetBranchAddress("L.tr.vz",Lvz);

 if(mode=="G"){
 T->SetBranchStatus("ctime",1);    
 // T->SetBranchAddress("ctime",ctime);
 T->SetBranchAddress("ctime",ctime);
 // T->SetBranchAddress("ctime",&coin_t);
 T->SetBranchStatus("DR.T5",1);    
 T->SetBranchAddress("DR.T5",&DRT5);

 }else if(mode=="T"){
 T->SetBranchStatus("ct",1);    
 T->SetBranchAddress("ct",&tcoin_t);
 // T->SetBranchAddress("ct",&coin_t);
 T->SetBranchStatus("ac1_sum",1);
 T->SetBranchAddress("ac1_sum",&Ra1sum);
 T->SetBranchStatus("ac2_sum",1);
 T->SetBranchAddress("ac2_sum",&Ra2sum);
 T->SetBranchStatus("mm",1);
 T->SetBranchAddress("mm",&mm);
 T->SetBranchStatus("Lz",1);    
 T->SetBranchAddress("Lz",&Lz);
 T->SetBranchStatus("Rz",1);    
 T->SetBranchAddress("Rz",&Rz);
 T->SetBranchStatus("Rp",1);
 T->SetBranchAddress("Rp",&Rpz);
 T->SetBranchStatus("Lp",1);
 T->SetBranchAddress("Lp",&Lpz);
 T->SetBranchStatus("ct_acc",1);
 T->SetBranchAddress("ct_acc",&ct_acc);
}

 //==================================================//  
 // Time scale [ns] //
 // Energy scale [GeV] //

 double ac1_adc[nth],ac2_adc[nth];
 double min_coin,max_coin,min_coin_c,max_coin_c;
 double min_ac1,max_ac1,min_ac2,max_ac2,min_adc,max_adc;
 double th1_max,th2_max,th2_min;
 double ac1_kcut,ac2_kcut_min,ac2_kcut_max;
 // if(mode=="H" || mode=="T"){

 //Paremater Setting //
 // AC set ACP//
   if(mode=="H"){
   min_coin=-10;
   max_coin=20.0;
   min_coin_c=-10;
   max_coin_c=20.0;

 //  min_coin_c=-100;
 //  max_coin_c=1000.0;
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




   }else if(mode=="T"){
   min_coin=-20;
   max_coin=10.0;
   min_coin_c=-20;
   max_coin_c=10.0;

 //==== NPE =======//
 min_ac1=0.0;
 max_ac1=10.0;
 //min_ac2=0.0;
 min_ac2=0.0;
 max_ac2=60.;
 //th1_max=5.0;
 th1_max=max_ac1;
 ac1_adc[0]=0.1;
 ac1_adc[1]=th1_max;
 ac1_adc[2]=th1_max;
 //ac1_adc[1]=0.05;
 //ac1_adc[2]=0.01;


 if(ac2_min){
 ac2_adc[0]=5.0;
 ac2_adc[1]=6.0;
 th2_max=60.;
 ac2_adc[0]=60;
 ac2_adc[1]=60;
 ac2_adc[2]=60;
 th2_min=3;
 // ac2_adc[2]=8.0;
 
 }else{

 th2_max=10.;
 max_ac2=th2_max;
 min_ac2=3.6; 
 ac2_adc[0]=10.0;
 ac2_adc[1]=8.0;
 ac2_adc[2]=7.0;
 th2_min=3.6;
 min_ac2=3.0;

}


 //---Kaon Cut ----//
 ac1_kcut=1.0;
 ac2_kcut_min=3.0;
 ac2_kcut_max=7.0;

 /*
 ac2_adc[0]=max_ac2;
 ac2_adc[1]=6000.;
 ac2_adc[2]=5000.;
 */
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


 double bin_vdc,min_vdc,max_vdc;
 min_vdc=-0.2e-6;
 max_vdc= 1.2e-6;
 bin_vdc=(max_vdc-min_vdc)/tdc_time;
 bin_vdc=(int)bin_vdc;
 double min_s0,max_s0;
 min_s0=-10;
 max_s0=10000;
 int bin_s0=int(max_s0-min_s0);
 
 double min_s2,max_s2;
 min_s2=-10;
 max_s2=5000;
 int bin_s2=max_s2-min_s2;


 double bin_coin=(max_coin-min_coin)/tdc_time;
        bin_coin=(int)bin_coin;
 double bin_coin_c=(max_coin_c-min_coin_c)/tdc_time;
        bin_coin_c=(int)bin_coin_c;
 int bin_beta=6000;
 int bin_adc=max_adc-min_adc;
  int bin_ac1=(max_ac1-min_ac1)*3; 
  int bin_ac2=(max_ac2-min_ac2)*3; 

  if(mode=="T"){
  bin_ac1=1000.; 
  bin_ac2=(max_ac2-min_ac2)*100; 
}

 TH2F* hcoin_ac1[nth];
 TH2F* hcoin_ac2[nth];
 TH2F* hcoin_ac1_acc[nth];
 TH2F* hcoin_ac2_acc[nth];
 TH2F* hvdc1_ac1[nth];
 TH2F* hvdc2_ac1[nth];
 TH2F* hvdc1_ac2[nth];
 TH2F* hvdc2_ac2[nth];
 TH2F* hs0_ac1[nth];
 TH2F* hs2_ac1[nth];
 TH2F* hs0_ac2[nth];
 TH2F* hs2_ac2[nth];
 TH1F* hcoin_t1[nth];
 TH1F* hcoin_t2[nth];
 TH1F* hcoin_ac1_max[nth];
 TH1F* hcoin_ac2_max[nth];
 TH1F* hcoin_t3[nth][nth];

 TH1F* hcoin_t=new TH1F("hcoin_t","Coincidence time S2R-S2L[ns] ",(int)bin_coin,min_coin,max_coin);
 TH1F* hcoin_tc=new TH1F("hcoin_tc","Coincidence time w/ Path Length Correction  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_acc_ac1[nth];
 TH1F* hcoin_acc_ac2[nth];

 TH2F* ha1_a2=new TH2F("ha1_a2","ac1 vs ac2 ADC sum hist",bin_ac1,min_ac1,max_ac1,bin_ac2,min_ac2,max_ac2);
 set->SetTH2(ha1_a2,"","","");

 TH1F* hcoin_k=new TH1F("hcoin_k","Coincidence time w/ Correction Kaon Cut  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_pi=new TH1F("hcoin_pi","Coincidence time w/ Correction Pion  Cut S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_p=new TH1F("hcoin_p","Coincidence time w/ Correction Proton  Cut S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 TH2F* hcoin_ac1_all=new TH2F("hcoin_ac1_all","Coinc time vs AC1 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac1,min_ac1,max_ac1);
 set->SetTH2(hcoin_ac1_all,"","","");
 TH2F* hcoin_ac2_all=new TH2F("hcoin_ac2_all","Coinc time vs AC2 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac2,min_ac2,max_ac2);
 set->SetTH2(hcoin_ac2_all,"","","");
 for(int i=0;i<nth;i++){
 hcoin_ac1[i]=new TH2F(Form("hcoin_ac1[%d]",i),"Coinc time vs AC1 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac1,min_ac1,max_ac1);
 set->SetTH2(hcoin_ac1[i],"","","");
 hcoin_ac2[i]=new TH2F(Form("hcoin_ac2[%d]",i),"Coinc time vs AC2 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac2,min_ac2,max_ac2);
 set->SetTH2(hcoin_ac2[i],"","","");
 hcoin_ac1_acc[i]=new TH2F(Form("hcoin_ac1_acc[%d]",i),"Coinc time (ACC) vs AC1 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac1,min_ac1,max_ac1);
 set->SetTH2(hcoin_ac1_acc[i],"","","");
 hcoin_ac2_acc[i]=new TH2F(Form("hcoin_ac2_acc[%d]",i),"Coinc time (ACC) vs AC2 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac2,min_ac2,max_ac2);
 set->SetTH2(hcoin_ac2_acc[i],"","","");
 hvdc1_ac1[i]=new TH2F(Form("hvdc1_ac1[%d]",i),"VDC1 Efficincy vs AC1 threshold",bin_vdc,min_vdc,max_vdc,bin_ac1,min_ac1,max_ac1);
 hvdc1_ac2[i]=new TH2F(Form("hvdc1_ac2[%d]",i),"VDC1 Efficincy vs AC2 threshold",bin_vdc,min_vdc,max_vdc,bin_ac2,min_ac2,max_ac2);
 hvdc2_ac1[i]=new TH2F(Form("hvdc2_ac1[%d]",i),"VDC2 Efficincy vs AC1 threshold",bin_vdc,min_vdc,max_vdc,bin_ac1,min_ac1,max_ac1);
 hvdc2_ac2[i]=new TH2F(Form("hvdc2_ac2[%d]",i),"VDC2 Efficincy vs AC2 threshold",bin_vdc,min_vdc,max_vdc,bin_ac2,min_ac2,max_ac2);
 hs0_ac1[i]=new TH2F(Form("hs0_ac1[%d]",i),"S0 Efficiency vs AC1 threshold",bin_s0,min_s0,max_s0,bin_ac1,min_ac1,max_ac1);
 hs0_ac2[i]=new TH2F(Form("hs0_ac2[%d]",i),"S0 Efficiency vs AC2 threshold",bin_s0,min_s0,max_s0,bin_ac2,min_ac2,max_ac2);
 hs2_ac1[i]=new TH2F(Form("hs2_ac1[%d]",i),"S2 Efficiency vs AC1 threshold",bin_s2,min_s2,max_s2,bin_ac1,min_ac1,max_ac1);
 hs2_ac2[i]=new TH2F(Form("hs2_ac2[%d]",i),"S2 Efficiency vs AC2 threshold",bin_s2,min_s2,max_s2,bin_ac2,min_ac2,max_ac2);

 hcoin_t1[i]=new TH1F(Form("hcoin_t1[%d]",i), Form("Coincidence 0<AC1<%lf cut",ac1_adc[i]),bin_coin_c,min_coin_c,max_coin_c);
 hcoin_t2[i]=new TH1F(Form("hcoin_t2[%d]",i), Form("Coincidence %lf<AC2<%lf  cut",ac1_adc[i],th2_max),bin_coin_c,min_coin_c,max_coin_c);
 hcoin_ac1_max[i]=new TH1F(Form("hcoin_ac1_max[%d]",i), Form("Coincidence time %lf<AC2<%lf cut",ac2_adc[i],th2_max),bin_coin_c,min_coin_c,max_coin_c);
 set->SetTH1(hcoin_ac1_max[i],"","","");
 hcoin_ac2_max[i]=new TH1F(Form("hcoin_ac2_max[%d]",i), Form("Coincidence time AC1<%lf  cut",ac1_adc[i]),bin_coin_c,min_coin_c,max_coin_c);
 set->SetTH1(hcoin_ac2_max[i],"","","");
 hcoin_acc_ac1[i]=new TH1F(Form("hcoin_acc_ac1[%d]",i),"Coincidence time ACC BG ",bin_coin_c,min_coin_c,max_coin_c);
 set->SetTH1(hcoin_acc_ac1[i],"","","");
 hcoin_acc_ac2[i]=new TH1F(Form("hcoin_acc_ac2[%d]",i),"Coincidence time ACC BG ",bin_coin_c,min_coin_c,max_coin_c);
 set->SetTH1(hcoin_acc_ac2[i],"","","");
 for(int j=0;j<nth;j++){
  hcoin_t3[i][j]=new TH1F(Form("hcoin_t3[%d][%d]",i,j),Form("Coincidence time S2R-S2L[ns] ac1_adc< %lf, ac2_adc<%lf; Coin-Time [ns];Counts ",ac1_adc[i],ac2_adc[j]),bin_coin_c,min_coin_c,max_coin_c);
  set->SetTH1(hcoin_t3[i][j],"","","");}

 }


 
 double mtr;
 mtr=938.27e-3;// proton mass [GeV/c^2]
 // double mh;

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



 int ev=0;
 double ct;
 for(int k=0;k<evnt;k++){
   T->GetEntry(k);

   if(k==ev*100000){
     cout<<"Fill Event: "<<k<<"/"<<evnt<<endl; 
     ev=ev+1;
   }


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
 rbeta=pk/Ek; 
 rpath_corr=rpathl/rbeta/c;
 lbeta=1.0;//pe_/Ee_; 
 lpath_corr=lpathl/lbeta/c;
 
 

 
 if(mode=="H"){
 Rs2_off=s2f1_off(Rs2pads,"R",mode,kine);
 Ls2_off=s2f1_off(Ls2pads,"L",mode,kine);
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
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

   if(mode=="H" || mode=="G"){
     Rz=Rvz[0];
     Lz=Lvz[0];
    }
   // if(Rvz_cutmin<Rvz[0] && Rvz[0]<Rvz_cutmax && Lvz_cutmin<Lvz[0] && Lvz[0]<Lvz_cutmax)cut_vz=true;
  
 
   if(Rvz_cutmin<Rz && Rz<Rvz_cutmax && Lvz_cutmin<Lz && Lz<Lvz_cutmax)cut_vz=true;

   if(mode=="G"){
     coin_t=ctime[0];
     coin_tc=ctime[0];
     if(coin_t>-100)coin_trig=true;
     cut_track=true;
     cut_rpathl=true;
     cut_lpathl=true;
  }else if(mode=="T"){
     coin_tc=tcoin_t;
     coin_t=tcoin_t;
  if(-80.<coin_tc && coin_tc<80.)coin_trig=true;
     cut_track=true;
     cut_rpathl=true;
     cut_lpathl=true; 
  }else if(mode=="H"){
 
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


 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track){
 hcoin_t->Fill(coin_t);
 hcoin_tc->Fill(coin_tc);
 ha1_a2->Fill(Ra1sum,Ra2sum);// AC1 vs AC2
 hcoin_ac1_all->Fill(coin_tc,Ra1sum);
 hcoin_ac2_all->Fill(coin_tc,Ra2sum);
 //--------- with AC1 Cut ---------------// 
  for(int j=0;j<nth;j++){
   cut_ac1=false;
   if(Ra1sum<ac1_adc[j])cut_ac1=true;
   //if(Ra1sum<0.04)cut_ac1=true;

   if(ac2_min && cut_ac1){
    hcoin_t2[j]->Fill(coin_tc); //AC1 cut   
    hcoin_ac2[j]->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
    hcoin_ac2_max[j]->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
 //------- ACC BG -------------------------------------------------------------//
     
  if(mode=="T" &&((-55.<coin_tc && coin_tc <-15.) || (5.<coin_tc && coin_tc<65.)) ){
       ct=coin_tc;       
        while(1){
	  if(-15.0<ct && ct<5.0){
		 hcoin_acc_ac2[j]->Fill(ct);
                 hcoin_acc_ac2[j]->Fill(ct-20);
		 hcoin_ac2_acc[j]->Fill(ct,Ra2sum);
		 break;}
	       else if(ct<-15.0){ct=ct+20.0;}
	       else if(5.0<ct){ct=ct-20.;}
	 }
  }
    //----------------------------------------------------------------------------//	  

   }else if(ac2_min==0 && cut_ac1){
	    // && Ra2sum>th2_min){
    hcoin_t2[j]->Fill(coin_tc); //AC1 cut   
    hcoin_ac2[j]->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
    hcoin_ac2_max[j]->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
 //------- ACC BG -------------------------------------------------------------//
     
  if(mode=="T" &&((-55.<coin_tc && coin_tc <-15.) || (5.<coin_tc && coin_tc<65.)) ){
       ct=coin_tc;       
        while(1){
	  if(-15.0<ct && ct<5.0){
		 hcoin_acc_ac2[j]->Fill(ct);
                 hcoin_acc_ac2[j]->Fill(ct-20);
		 hcoin_ac2_acc[j]->Fill(ct,Ra2sum);
		 break;}
	       else if(ct<-15.0){ct=ct+20.0;}
	       else if(5.0<ct){ct=ct-20.;}
	 }
     }
    //----------------------------------------------------------------------------//	
   }  
  }

 
 
  //-------with AC2 Cut --------------------//
   for(int i=0;i<nth;i++){
       cut_ac2=false;
       // if(Ra2sum>ac2_adc[i])cut_ac2=true;
    //    cut_ac2=true;

       //    if(ac2_min &&Ra2sum>ac2_adc[i]&&Ra2sum< ac2_kcut_max)cut_ac2=true;
    if(ac2_min &&Ra2sum>3.0 && Ra2sum<7.0)cut_ac2=true;
    if(ac2_min==0 && Ra2sum<ac2_adc[i] &&Ra2sum>th2_min)cut_ac2=true;
    // if(6.0>Ra2sum && Ra2sum>3.0)cut_ac2=true;
    if(cut_ac2){
    hcoin_t1[i]->Fill(coin_tc); // AC2 cut
    hcoin_ac1[i]->Fill(coin_tc,Ra1sum);//AC2 cut && (AC1 variable cut)  B
    hcoin_ac1_max[i]->Fill(coin_tc,Ra1sum);//AC2 cut && (AC1 variable cut)  


 //------- ACC BG -------------------------------------------------------------//
 
   if(mode=="T" &&((-55.<coin_tc && coin_tc <-15.) || (5.<coin_tc && coin_tc<65.)) ){
        ct=coin_tc;       
        while(1){
	  if(-15.0<ct && ct<5.0){
		 hcoin_acc_ac1[i]->Fill(ct);
                 hcoin_acc_ac1[i]->Fill(ct-20);
		 hcoin_ac1_acc[i]->Fill(ct,Ra1sum);
		 break;}
	       else if(ct<-15.0){ct=ct+20.0;}
	       else if(5.0<ct){ct=ct-20.0;}
	 }
    }
 
    //----------------------------------------------------------------------------//

    }
   }
 }



 //-------with AC1 && AC2 Cut --------------------//
  cut_ac1=false;
  cut_ac2=false; 
   for(int k=0;k<nth;k++){
   for(int j=0;j<nth;j++){ 
  if(Ra1sum<ac1_adc[j])cut_ac1=true;
  if(Ra2sum<ac2_adc[k])cut_ac2=true;
  if(cut_ac1 && cut_ac2)hcoin_t3[j][k]->Fill(coin_tc); //AC1 & AC2 cut
  } 
   }
   //--------- Kaon Cut Hist ---------------//
 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track && Ra1sum<ac1_kcut 
    && ac2_kcut_min<Ra2sum && Ra2sum<ac2_kcut_max)hcoin_k->Fill(coin_tc);
   //--------- Pion Cut Hist ---------------//
 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track && Ra1sum>ac1_kcut 
    && Ra2sum>ac2_kcut_max)hcoin_pi->Fill(coin_tc);
   //-------- Proton Cut Hist --------------//
 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track && Ra1sum<ac1_kcut 
    && Ra2sum<ac2_kcut_min)hcoin_p->Fill(coin_tc);

 }

 for(int i=0;i<nth;i++){
 hcoin_acc_ac1[i]->Scale(20./100.);
 hcoin_acc_ac2[i]->Scale(20./100.);
 // hcoin_ac1_acc[i]->Scale(20./100.);
 // hcoin_ac2_acc[i]->Scale(20./100.);
 }


 cout<<"Filled Hist "<<endl;

 //==========================================================//
 //============== Survival Rate (SR) analysis =======================//
 //==========================================================//
 int iter_ac1=50; int iter_ac2=50; int iter_max=iter_ac1+iter_ac2;
 if(test_flag){iter_ac1=1;}
 if(pra_flag){iter_ac1=5;}
 TF1* facc[iter_max][nth][2];
 TF1* fpi[iter_max][nth][2];
 TF1* fk[iter_max][nth][2];
 TF1* fcoin[iter_max][nth][2];
 TF1* fp[iter_max][nth][2];
 TH2F* hfom_ac[nth][nth];
 TH1D* hcoin_ac1_p[iter_max][nth];
 TH1D* hcoin_ac2_p[iter_max][nth]; 
 TH1D* hcoin_ac1_all_p[iter_max][nth];
 TH1D* hcoin_ac2_all_p[iter_max][nth]; 
 TH1D* hcoin_ac1_acc_p[iter_max][nth];
 TH1D* hcoin_ac2_acc_p[iter_max][nth]; 
 TGraphErrors* gsum_pi_ac1[nth][nth];
 TGraphErrors* gsum_p_ac1[nth][nth];
 TGraphErrors* gsum_k_ac1[nth][nth];
 TGraphErrors* grate_k_ac1[nth][nth];
 TGraphErrors* grate_p_ac1[nth][nth];
 TGraphErrors* grate_pi_ac1[nth][nth];
 TGraphErrors* gsum_pi_ac2[nth][nth];
 TGraphErrors* gsum_p_ac2[nth][nth];
 TGraphErrors* gsum_k_ac2[nth][nth];
 TGraphErrors* grate_k_ac2[nth][nth];
 TGraphErrors* grate_pi_ac2[nth][nth];
 TGraphErrors* grate_p_ac2[nth][nth];
 TGraphErrors* gSN_k_ac1[nth][nth];
 TGraphErrors* gSN_k_ac2[nth][nth];
 TGraphErrors* gfom_ac1[nth];
 TGraphErrors* gfom_ac2[nth];
 TGraphErrors* gfom=new TGraphErrors();
 set->SetGr(gfom,"","","");
 for(int k=0;k<nth;k++){
   gfom_ac1[k]=new TGraphErrors();
   gfom_ac2[k]=new TGraphErrors();
   set->SetGr(gfom_ac1[k],"","","");
   set->SetGr(gfom_ac2[k],"","","");
   gfom_ac1[k]->SetTitle(Form("FOM Estimation Graph %lf<AC2<%lf;AC1 ADC ;FOM ",ac2_adc[k],th2_max));
   gfom_ac2[k]->SetTitle(Form("FOM Estimation Graph AC1<%lf;AC2 ADC;FOM",ac1_adc[k]));
 gfom_ac1[k]->SetMarkerStyle(21);
 gfom_ac1[k]->SetMarkerColor(kRed);
 gfom_ac1[k]->SetMarkerSize(0.5);
 gfom_ac2[k]->SetMarkerStyle(21);
 gfom_ac2[k]->SetMarkerColor(kRed);
 gfom_ac2[k]->SetMarkerSize(0.5);
 for(int j=0;j<nth;j++){
 gSN_k_ac1[j][k]=new TGraphErrors();
 set->SetGr(gSN_k_ac1[j][k],"","","");
 gSN_k_ac2[j][k]=new TGraphErrors();
 set->SetGr(gSN_k_ac2[j][k],"","","");
 gsum_pi_ac1[j][k]=new TGraphErrors();
 gsum_p_ac1[j][k]=new TGraphErrors();
 gsum_k_ac1[j][k]=new TGraphErrors();
 set->SetGr(gsum_k_ac1[j][k],"","","");
 grate_k_ac1[j][k]=new TGraphErrors();
 set->SetGr(grate_k_ac1[j][k],"","","");
 grate_p_ac1[j][k]=new TGraphErrors();
 set->SetGr(grate_p_ac1[j][k],"","","");
 grate_pi_ac1[j][k]=new TGraphErrors();
 set->SetGr(grate_pi_ac1[j][k],"","","");
 gsum_pi_ac2[j][k]=new TGraphErrors();
 set->SetGr(gsum_pi_ac2[j][k],"","","");
 gsum_p_ac2[j][k]=new TGraphErrors();
 set->SetGr(gsum_pi_ac2[j][k],"","","");
 gsum_k_ac2[j][k]=new TGraphErrors();
 set->SetGr(gsum_k_ac2[j][k],"","","");
 grate_k_ac2[j][k]=new TGraphErrors();
 set->SetGr(grate_k_ac2[j][k],"","","");
 grate_pi_ac2[j][k]=new TGraphErrors();
 set->SetGr(grate_pi_ac2[j][k],"","","");
 grate_p_ac2[j][k]=new TGraphErrors();
 set->SetGr(grate_p_ac2[j][k],"","","");
 //--- TGraphErrors Setting ----//
 gsum_pi_ac1[j][k]->SetTitle(Form("Sum of Pion (AC2<%lf);AC1 ADC-th [ch];Pion Events [Counts]",ac2_adc[k]));
 gsum_p_ac1[j][k]->SetTitle(Form("SUM of Proton (AC2<%lf);AC1 ADC-th [ch];Proton Events [Counts]",ac2_adc[k]));
 gsum_k_ac1[j][k]->SetTitle(Form("SUM of Kaon (AC2<%lf);AC1 ADC-th [ch];Kaon Events [Counts]",ac2_adc[k]));
 grate_k_ac1[j][k]->SetTitle(Form("Kaoin Survival rate ((AC2<%lf));AC1 ADC-th [ch]; Survival rate ",ac2_adc[k]));
 grate_pi_ac1[j][k]->SetTitle(Form("Pion Survuval rate (AC<%lf);AC1 ADC-th [ch]; Surival rate ",ac2_adc[k]));
 grate_p_ac1[j][k]->SetTitle("Proton Survival rate  vs AC1 Threshold ;AC1 ADC-th [ch]; Survival rate "); 
 gSN_k_ac1[j][k]->SetTitle(Form("Kaon S/N rate AC1 [%d][%d]",j,k));
 gsum_pi_ac1[j][k]->SetMarkerStyle(21);
 gsum_pi_ac1[j][k]->SetMarkerColor(kRed);
 gsum_pi_ac1[j][k]->SetMarkerSize(0.5);
 gsum_p_ac1[j][k]->SetMarkerStyle(21);
 gsum_p_ac1[j][k]->SetMarkerColor(kRed);
 gsum_p_ac1[j][k]->SetMarkerSize(0.5);
 gsum_k_ac1[j][k]->SetMarkerStyle(21);
 gsum_k_ac1[j][k]->SetMarkerColor(kRed);
 gsum_k_ac1[j][k]->SetMarkerSize(0.5);
 grate_k_ac1[j][k]->SetMarkerStyle(21);
 grate_k_ac1[j][k]->SetMarkerColor(kBlue);
 grate_k_ac1[j][k]->SetMarkerSize(0.5);
 grate_p_ac1[j][k]->SetMarkerStyle(21);
 grate_p_ac1[j][k]->SetMarkerColor(kBlue);
 grate_p_ac1[j][k]->SetMarkerSize(0.5);
 grate_pi_ac1[j][k]->SetMarkerStyle(21);
 grate_pi_ac1[j][k]->SetMarkerColor(kBlue);
 grate_pi_ac1[j][k]->SetMarkerSize(0.5);
 gSN_k_ac1[j][k]->SetMarkerStyle(21);
 gSN_k_ac1[j][k]->SetMarkerColor(kRed);
 gSN_k_ac1[j][k]->SetMarkerSize(0.5);
 gsum_pi_ac2[j][k]->SetTitle("SUM of Pion vs AC2 Threshold;AC2 ADC-th [ch];Pion Events [Counts]");
 gsum_p_ac2[j][k]->SetTitle("SUM of Proton vs AC2 Threshold;AC2 ADC-th [ch];Proton Events [Counts]");
 gsum_k_ac2[j][k]->SetTitle("SUM of Kaon vs AC2 Threshold;AC2 ADC-th [ch];Kaon Events [Counts]");
 grate_k_ac2[j][k]->SetTitle("Kaoin Survival rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 grate_pi_ac2[j][k]->SetTitle("Pion Survival rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 grate_p_ac2[j][k]->SetTitle("Proton Survival rate  vs AC2 Threshold ;AC2 ADC-th [ch]; Survival rate ");
 gSN_k_ac2[j][k]->SetTitle(Form("Kaon S/N rate AC2 [%d][%d]",j,k));
 gsum_pi_ac2[j][k]->SetMarkerStyle(21);
 gsum_pi_ac2[j][k]->SetMarkerColor(kRed);
 gsum_pi_ac2[j][k]->SetMarkerSize(0.5);
 gsum_p_ac2[j][k]->SetMarkerStyle(21);
 gsum_p_ac2[j][k]->SetMarkerColor(kRed);
 gsum_p_ac2[j][k]->SetMarkerSize(0.5);
 gsum_k_ac2[j][k]->SetMarkerStyle(21);
 gsum_k_ac2[j][k]->SetMarkerColor(kRed);
 gsum_k_ac2[j][k]->SetMarkerSize(0.5);
 grate_k_ac2[j][k]->SetMarkerStyle(21);
 grate_k_ac2[j][k]->SetMarkerColor(kBlue);
 grate_k_ac2[j][k]->SetMarkerSize(0.5);
 grate_p_ac2[j][k]->SetMarkerStyle(21);
 grate_p_ac2[j][k]->SetMarkerColor(kBlue);
 grate_p_ac2[j][k]->SetMarkerSize(0.5);
 grate_pi_ac2[j][k]->SetMarkerStyle(21);
 grate_pi_ac2[j][k]->SetMarkerColor(kBlue);
 grate_pi_ac2[j][k]->SetMarkerSize(0.5);
 gSN_k_ac2[j][k]->SetMarkerStyle(21);
 gSN_k_ac2[j][k]->SetMarkerColor(kRed);
 gSN_k_ac2[j][k]->SetMarkerSize(0.5);
 }
 }



//--- Parameters -----//

 double kmin[iter_max][nth][2],kmax[iter_max][nth][2];
 double inte_ktot[iter_max][nth][2], inte_ksig[iter_max][nth][2];
 double p0_acc[iter_max][nth][2], p1_acc[iter_max][nth][2];
 double n_p[iter_max][nth][2],sig_p[iter_max][nth][2],mean_p[iter_max][nth][2];
 double n_pi[iter_max][nth][2],sig_pi[iter_max][nth][2],mean_pi[iter_max][nth][2];
 double n_k[iter_max][nth][2],sig_k[iter_max][nth][2],mean_k[iter_max][nth][2];
 int bin_ac1_adc[nth][nth],bin_min_ac1,bin_max_ac1,bin_ac2_adc[nth][nth],bin_max_ac2,bin_min_ac2;
 double sum_k[iter_max][nth][2],sum_p[iter_max][nth][2],sum_pi[iter_max][nth][2]; 
  double sum_k_err[iter_max][nth][2],sum_p_err[iter_max][nth][2],sum_pi_err[iter_max][nth][2]; 
double inte_acc[iter_max][nth][2];
 double th_ac1[iter_max],th_ac2[iter_max];
 int bin_th_ac1[iter_max][nth],bin_th_ac2[iter_max][nth]; 
 double nk[iter_max][nth][iter_max][nth][2],npi[iter_max][nth][iter_max][nth][2],np[iter_max][nth][iter_max][nth][2];
 double max_nk[nth][nth][2],max_npi[nth][nth][2],max_np[nth][nth][2];
 double n_p_err[iter_max][nth][2],n_pi_err[iter_max][nth][2],n_k_err[iter_max][nth][2];
 double FOM_ac1[iter_max][nth],FOM_ac2[iter_max][nth];
 double max_fom_ac1,max_fom_ac2;
 int fom_th1,fom_th2;

 max_fom_ac1=0.0;
 max_fom_ac2=0.0;

//---- Defolt parameters -----------//

 TF1* facc_t1def[nth][nth];
 TF1* fpi_t1def[nth][nth];
 TF1* fk_t1def[nth][nth];
 TF1* fcoin_t1def[nth][nth];
 TF1* fp_t1def[nth][nth];
 TF1* facc_t2def[nth][nth];
 TF1* fpi_t2def[nth][nth];
 TF1* fk_t2def[nth][nth];
 TF1* fcoin_t2def[nth][nth];
 TF1* fp_t2def[nth][nth];
 TF1* facc_t3def[nth][nth];
 TF1* fpi_t3def[nth][nth];
 TF1* fk_t3def[nth][nth];
 TF1* fcoin_t3def[nth][nth];
 TF1* fp_t3def[nth][nth];
 TF1* fcoin_t1[nth][nth];
 TF1* fcoin_t2[nth][nth];  
 TF1* fcoin_t3[nth][nth]; 


 double def_sig_p,def_mean_p,def_sig_pi,def_mean_pi,def_sig_k,def_mean_k,def_acc;
 double def_num_k,def_num_p,def_num_pi,def_acc_k,def_acc_pi,def_acc_p;

 if(mode=="H" && kine==2){
 def_sig_p=0.852; def_mean_p=0.0;
 def_sig_pi=0.443; def_mean_pi=11;
 def_sig_k=0.644; def_mean_k=8.0;
 def_acc=27.7;
 }else if(mode=="H" && kine==1){
 def_sig_p=0.852; def_mean_p=0.0;
 def_sig_pi=0.443; def_mean_pi=11;
 def_sig_k=0.644; def_mean_k=8.;
 def_acc=27.7;
 }else if(mode=="G"){
 def_sig_p=0.852; def_mean_p=11.06;
 def_sig_pi=0.4; def_mean_pi=0.0;
 def_sig_k=0.644; def_mean_k=3.16;
 def_acc=22.7;
 }else if(mode=="T"){
 def_sig_p=0.841; def_mean_p=-8.21;
 def_sig_pi=0.275; def_mean_pi=3.049;
 def_sig_k=0.454; def_mean_k=-0.095;
 def_acc=22.7;
}



 double def_t1_k[nth][nth],def_t1_pi[nth][nth],def_t1_p[nth][nth],def_t1_acc[nth][nth];
 double def_t1_k_err[nth][nth],def_t1_pi_err[nth][nth],def_t1_p_err[nth][nth],def_t1_acc_err[nth][nth];
 double t1sig_k[nth][nth],t1sig_p[nth][nth],t1sig_pi[nth][nth],t1mean_p[nth][nth],t1mean_k[nth][nth],t1mean_pi[nth][nth];
 double t1sum_k[nth],t1sum_pi[nth],t1sum_p[nth];
 double t1sum_k_err[nth],t1sum_pi_err[nth],t1sum_p_err[nth];
 double def_t2_k[nth][nth],def_t2_pi[nth][nth],def_t2_p[nth][nth],def_t2_acc[nth][nth];
 double t2sig_k[nth][nth],t2sig_p[nth][nth],t2sig_pi[nth][nth],t2mean_p[nth][nth],t2mean_k[nth][nth],t2mean_pi[nth][nth];
 double def_t2_k_err[nth][nth],def_t2_pi_err[nth][nth],def_t2_p_err[nth][nth],def_t2_acc_err[nth][nth];
 double t2sum_k[nth],t2sum_pi[nth],t2sum_p[nth];
 double t2sum_k_err[nth],t2sum_pi_err[nth],t2sum_p_err[nth];
 double def_t3_k[nth][nth],def_t3_pi[nth][nth],def_t3_p[nth][nth],def_t3_acc[nth][nth];
 double t3sig_k[nth][nth],t3sig_p[nth][nth],t3sig_pi[nth][nth],t3mean_p[nth][nth],t3mean_k[nth][nth],t3mean_pi[nth][nth];
 double t3sum_k[nth][nth],t3sum_pi[nth][nth],t3sum_p[nth][nth];
 double emp[iter_max];

 double rate_k[iter_max][nth][2],rate_p[iter_max][nth][2],rate_pi[iter_max][nth][2];
 double rate_k_err[iter_max][nth][2],rate_p_err[iter_max][nth][2],rate_pi_err[iter_max][nth][2];
 double sum_acc[iter_max][nth][2];
 double max_SN_ac1[nth],max_SN_ac2[nth];
 int SN_ac1[nth],SN_ac2[nth];

 double sum_k_max=1250.;
 double fom_max=0.0;

 for(int i=0;i<nth;i++){
 max_SN_ac1[i]=0.0;
 max_SN_ac2[i]=0.0;
 SN_ac1[i]=0;
 SN_ac2[i]=0;
 }


 for(int i=0;i<iter_max;i++)emp[i]=0.0;
 //bool point_err=true;


cout<<"defined parameters in SR analysis"<<endl;
 
         //----------------------------------//

 TF1* facc_kc=new TF1("facc_kc","[0]",min_coin_c,max_coin_c);
 facc_kc->SetNpx(2000);
 TF1* fk_kc=new TF1("fk_kc","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fk_kc->SetNpx(2000);
 TF1* fpi_pic=new TF1("fpi_pic","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fpi_pic->SetNpx(2000);
 TF1* fp_pc=new TF1("fp_pc","gausn(0)+pol0(3)",min_coin_c,max_coin_c);
 fp_pc->SetNpx(2000);


 fpi_pic->SetParameter(1,def_mean_pi);
 fpi_pic->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fpi_pic->SetParameter(2,def_sig_pi);
 fpi_pic->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fpi_pic->SetParameter(3,def_acc);

 fp_pc->SetParameter(1,def_mean_p);
 fp_pc->SetParLimits(1,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fp_pc->SetParameter(2,def_sig_p);
 fp_pc->SetParLimits(2,0.8*def_sig_p,1.2*def_sig_p);
 fp_pc->SetParameter(3,def_acc);

 hcoin_k->Fit("facc_kc","Rq","",min_coin_c,min_coin_c+3.0);
 def_acc_k=facc_kc->GetParameter(0);

 fk_kc->SetParameter(1,def_mean_k);
 fk_kc->SetParLimits(1,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fk_kc->SetParameter(2,def_sig_k);
 fk_kc->SetParLimits(2,0.8*def_sig_k,1.2*def_sig_k);
 fk_kc->FixParameter(3,def_acc_k);

 hcoin_k->Fit("fk_kc","Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_pi->Fit("fpi_pic","Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 hcoin_p->Fit("fp_pc","Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);

 def_num_k=fk_kc->GetParameter(0);
 def_mean_k=fk_kc->GetParameter(1);
 def_sig_k=fk_kc->GetParameter(2);

 def_num_p=fp_pc->GetParameter(0);
 def_mean_p=fp_pc->GetParameter(1);
 def_sig_p=fp_pc->GetParameter(2);

 def_num_pi=fpi_pic->GetParameter(0);
 def_mean_pi=fpi_pic->GetParameter(1);
 def_sig_pi=fpi_pic->GetParameter(2);



 //============= SRA ======================//

     for(int th1=0;th1<nth;th1++){
     for(int th2=0;th2<nth;th2++){
     for(int i=0;i<iter_ac1;i++){

    cout<<"th1: "<<th1<<"/"<<nth<<":th2: "<<th2<<"/"<<nth<<": "
    <<"i "<<i<<"/"<<iter_ac1<<endl;

   
    th_ac1[i]=th1_max-(th1_max-min_ac1)/iter_ac1*i;
bin_th_ac1[i][th2]=hcoin_ac1[th2]->GetYaxis()->FindBin(th_ac1[i]);//hcoin_ac1[th2] is ac2 fixed threshold
 bin_min_ac1=hcoin_ac1[th2]->GetYaxis()->FindBin(min_ac1);

 hcoin_ac1_all_p[i][th2]=hcoin_ac1[th2]->ProjectionX(Form("hcoin_all_ac1_p[%d][%d]",i,th2),bin_min_ac1,bin_th_ac1[i][th2]);
 hcoin_ac1_all_p[i][th2]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",th_ac1[i],ac2_adc[th2],th2_max));
 set->SetTH1(hcoin_ac1_all_p[i][th2],"","","");

 hcoin_ac1_acc_p[i][th2]=hcoin_ac1_acc[th2]->ProjectionX(Form("hcoin_ac1_acc_p[%d][%d]",i,th2),bin_min_ac1,bin_th_ac1[i][th2]);
 hcoin_ac1_acc_p[i][th2]->SetTitle(Form("Coin-time (ACC) AC1<%lf & %lf<AC2<%lf",th_ac1[i],ac2_adc[th2],th2_max));
 set->SetTH1(hcoin_ac1_acc_p[i][th2],"","","");
 hcoin_ac1_acc_p[i][th2]->Scale(20./100.);

 hcoin_ac1_p[i][th2]=hcoin_ac1[th2]->ProjectionX(Form("hcoin_ac1_p[%d][%d]",i,th2),bin_min_ac1,bin_th_ac1[i][th2]);
 hcoin_ac1_p[i][th2]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",th_ac1[i],ac2_adc[th2],th2_max));
 set->SetTH1(hcoin_ac1_p[i][th2],"","","");
 hcoin_ac1_p[i][th2]->Add(hcoin_ac1_all_p[i][th2], hcoin_ac1_acc_p[i][th2],1.0,-1.0);
  

 //============ AC2===================//

  th_ac2[i]=min_ac2+(th2_max-min_ac2)/iter_ac1*i;

  bin_th_ac2[i][th1]=hcoin_ac2[th1]->GetYaxis()->FindBin(th_ac2[i]);
  bin_min_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(min_ac2);
  bin_max_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(max_ac2);
  hcoin_ac2_all_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_all_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2);
  set->SetTH1(hcoin_ac2_all_p[i][th1],"","","");
  hcoin_ac2_all_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));

  hcoin_ac2_acc_p[i][th1]=hcoin_ac2_acc[th1]->ProjectionX(Form("hcoin_ac2_acc_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2); 
  hcoin_ac2_acc_p[i][th1]->SetTitle(Form("Coin-time (ACC) AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));
  set->SetTH1(hcoin_ac2_acc_p[i][th1],"","","");
  hcoin_ac2_acc_p[i][th1]->Scale(20./100.);

  hcoin_ac2_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2); 
  hcoin_ac2_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));
  hcoin_ac2_p[i][th1]->Add(hcoin_ac2_all_p[i][th1],hcoin_ac2_acc_p[i][th1],1.0,-1.0);
  set->SetTH1(hcoin_ac2_p[i][th1],"","","");

 
 if(ac2_min){
  th_ac2[i]=min_ac2+(th2_max-min_ac2)/iter_ac1*i;

  bin_th_ac2[i][th1]=hcoin_ac2[th1]->GetYaxis()->FindBin(th_ac2[i]);
  bin_min_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(min_ac2);
  bin_max_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(max_ac2);
  hcoin_ac2_all_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_all_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2);
  set->SetTH1(hcoin_ac2_all_p[i][th1],"","","");
  hcoin_ac2_all_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));

  hcoin_ac2_acc_p[i][th1]=hcoin_ac2_acc[th1]->ProjectionX(Form("hcoin_ac2_acc_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2); 
  hcoin_ac2_acc_p[i][th1]->SetTitle(Form("Coin-time (ACC) AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));
  set->SetTH1(hcoin_ac2_acc_p[i][th1],"","","");
  hcoin_ac2_acc_p[i][th1]->Scale(20./100.);

  hcoin_ac2_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_p[%d][%d]",i,th1),bin_th_ac2[i][th1],bin_max_ac2); 
  hcoin_ac2_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));
  hcoin_ac2_p[i][th1]->Add(hcoin_ac2_all_p[i][th1],hcoin_ac2_acc_p[i][th1],1.0,-1.0);
  set->SetTH1(hcoin_ac2_p[i][th1],"","","");
 }else{

   th_ac2[i]= th2_max-(th2_max-min_ac2)/iter_ac1*i;
  bin_th_ac2[i][th1]=hcoin_ac2[th1]->GetYaxis()->FindBin(th_ac2[i]);
  bin_min_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(min_ac2);
  bin_max_ac2=hcoin_ac2[th1]->GetYaxis()->FindBin(max_ac2);
  hcoin_ac2_all_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_all_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]);  
  hcoin_ac2_all_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],th_ac2[i],th2_max));
  set->SetTH1(hcoin_ac2_all_p[i][th1],"","","");
  hcoin_ac2_acc_p[i][th1]=hcoin_ac2_acc[th1]->ProjectionX(Form("hcoin_ac2_acc_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]); 
  hcoin_ac2_acc_p[i][th1]->SetTitle(Form("Coin-time (ACC) AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],bin_min_ac2,th_ac2[i]));
  hcoin_ac2_acc_p[i][th1]->Scale(20./100.);
  set->SetTH1(hcoin_ac2_acc_p[i][th1],"","","");
  hcoin_ac2_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]); 
  hcoin_ac2_p[i][th1]->SetTitle(Form("Coin-time AC1<%lf & %lf<AC2<%lf",ac1_adc[th1],bin_min_ac2,bin_th_ac2[i][th1])); 
  hcoin_ac2_p[i][th1]->Add(hcoin_ac2_all_p[i][th1],hcoin_ac2_acc_p[i][th1],1.0,-1.0);
  set->SetTH1(hcoin_ac2_p[i][th1],"","","");
}
 


 // hcoin_ac2_p[i][th1]=hcoin_ac2[th1]->ProjectionX(Form("hcoin_ac2_p[%d][%d]",i,th1),bin_min_ac2,bin_th_ac2[i][th1]); //hcoin_ac2[th1] is ac1 fixed threshold

 //--- Initial Parameters -----------//

  facc[i][th2][0]=new TF1(Form("facc[%d][%d][0]",i,th2),"[0]",min_coin_c,max_coin_c);
  facc[i][th2][0]->SetNpx(2000);
  hcoin_ac1_p[i][th2]->Fit(Form("facc[%d][%d][0]",i,th2),"Rq","",min_coin_c,min_coin_c+3.);
  p0_acc[i][th2][0]=facc[i][th2][0]->GetParameter(0);
 
  facc[i][th1][1]=new TF1(Form("facc[%d][%d][1]",i,th1),"[0]",min_coin_c,max_coin_c);
  facc[i][th1][1]->SetNpx(2000);  
  hcoin_ac2_p[i][th1]->Fit(Form("facc[%d][%d][1]",i,th1),"Rq","",min_coin_c,min_coin_c+3.);
  p0_acc[i][th1][1]=facc[i][th1][1]->GetParameter(0);

 //------- AC1 --------------// 

 fp[i][th2][0]=new TF1(Form("fp[%d][%d][0]",i,th2),"gausn(0)",min_coin_c,max_coin_c);
 fp[i][th2][0]->SetNpx(2000);
 fpi[i][th2][0] =new TF1(Form("fpi[%d][%d][0]",i,th2),"gausn(0)",min_coin_c,max_coin_c);
 fpi[i][th2][0]->SetNpx(2000);
 fk[i][th2][0]=new TF1(Form("fk[%d][%d][0]",i,th2),"gausn(0)",min_coin_c,max_coin_c);
 fk[i][th2][0]->SetNpx(2000);


 fp[i][th2][0]->FixParameter(1,def_mean_p); 
 fp[i][th2][0]->FixParameter(2,def_sig_p); 
 
 // fpi[i][th2][0]->FixParameter(1,def_mean_pi); 
 //fpi[i][th2][0]->FixParameter(2,def_sig_pi); 
 //fk[i][th2][0]->FixParameter(1,def_mean_k); 
 //fk[i][th2][0]->FixParameter(2,def_sig_k);

 fpi[i][th2][0]->FixParameter(1,def_mean_pi); 
 fk[i][th2][0]->FixParameter(1,def_mean_k); 
 fpi[i][th2][0]->SetParLimits(2,0.95*def_sig_pi,1.05*def_sig_pi); 
 fk[i][th2][0]->SetParLimits(2,0.75*def_sig_k,1.25*def_sig_k);



//----- AC1 Fitting -----------// 
 hcoin_ac1_p[i][th2]->Fit(Form("fp[%d][%d][0]",i,th2),"Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 n_p[i][th2][0]=fp[i][th2][0]->GetParameter(0);
 mean_p[i][th2][0]=fp[i][th2][0]->GetParameter(1);
 sig_p[i][th2][0]=fp[i][th2][0]->GetParameter(2);
 hcoin_ac1_p[i][th2]->Fit(Form("fpi[%d][%d][0]",i,th2),"Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 n_pi[i][th2][0]=fpi[i][th2][0]->GetParameter(0);
 mean_pi[i][th2][0]=fpi[i][th2][0]->GetParameter(1);
 sig_pi[i][th2][0]=fpi[i][th2][0]->GetParameter(2);
 hcoin_ac1_p[i][th2]->Fit(Form("fk[%d][%d][0]",i,th2),"Rq","",def_mean_p-3*def_sig_k,def_mean_k+3*def_sig_k);
 n_k[i][th2][0]=fk[i][th2][0]->GetParameter(0);
 mean_k[i][th2][0]=fk[i][th2][0]->GetParameter(1);
 sig_k[i][th2][0]=fk[i][th2][0]->GetParameter(2);

 //------- Get Error Paramters ---//
 n_pi_err[i][th2][0]=fpi[i][th2][0]->GetParError(0); 
 n_k_err[i][th2][0]=fk[i][th2][0]->GetParError(0); 
 n_p_err[i][th2][0]=fp[i][th2][0]->GetParError(0); 

 //----- AC1 Coint Fitting ---------//

 fcoin[i][th2][0] =new TF1(Form("fcoin[%d][%d][0]",i,th2),"gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin[i][th2][0]->SetNpx(2000);
 fcoin[i][th2][0]->SetTitle(Form("Coin-Time w AC cut  (AC1<%d ch && AC2>%d ch);Coin time [ns];Counts [1/56 ns]",th_ac1[i],ac2_adc));
 fcoin[i][th2][0]->SetParameters(n_pi[i][th2][0],mean_pi[i][th2][0],sig_pi[i][th2][0],n_k[i][th2][0],mean_k[i][th2][0],sig_k[i][th2][0],n_p[i][th2][0],mean_p[i][th2][0],sig_p[i][th2][0]);





 //------- AC2 -------------//

 fp[i][th1][1]=new TF1(Form("fp[%d][%d][1]",i,th1),"gausn(0)",min_coin_c,max_coin_c);
 fp[i][th1][1]->SetNpx(2000);
 fpi[i][th1][1] =new TF1(Form("fpi[%d][%d][1]",i,th1),"gausn(0)",min_coin_c,max_coin_c);
 fpi[i][th1][1]->SetNpx(2000);
 fk[i][th1][1]=new TF1(Form("fk[%d][%d][1]",i,th1),"gausn(0)",min_coin_c,max_coin_c);
 fk[i][th1][1]->SetNpx(2000);
 
 fpi[i][th1][1]->FixParameter(1,def_mean_pi); 
 fpi[i][th1][1]->FixParameter(2,def_sig_pi);
 
 fk[i][th1][1]->FixParameter(1,def_mean_k); 
 // fk[i][th1][1]->FixParameter(2,def_sig_k); 
 fk[i][th1][1]->SetParLimits(2,0.75*def_sig_k,1.25*def_sig_k); 
 fp[i][th1][1]->FixParameter(1,def_mean_p); 
 fp[i][th1][1]->FixParameter(2,def_sig_p); 


  //----- AC2 Fitting -----------// 
 hcoin_ac2_p[i][th1]->Fit(Form("fp[%d][%d][1]",i,th1),"Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);
 n_p[i][th1][1]=fp[i][th1][1]->GetParameter(0);
 mean_p[i][th1][1]=fp[i][th1][1]->GetParameter(1);
 sig_p[i][th1][1]=fp[i][th1][1]->GetParameter(2);
 hcoin_ac2_p[i][th1]->Fit(Form("fpi[%d][%d][1]",i,th1),"Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 n_pi[i][th1][1]=fpi[i][th1][1]->GetParameter(0);
 mean_pi[i][th1][1]=fpi[i][th1][1]->GetParameter(1);
 sig_pi[i][th1][1]=fpi[i][th1][1]->GetParameter(2);
 hcoin_ac2_p[i][th1]->Fit(Form("fk[%d][%d][1]",i,th1),"Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 n_k[i][th1][1]=fk[i][th1][1]->GetParameter(0);
 mean_k[i][th1][1]=fk[i][th1][1]->GetParameter(1);
 sig_k[i][th1][1]=fk[i][th1][1]->GetParameter(2);

 // n_k[i][th1][1]=hcoin_ac2_p[i][th1]->Integral(hcoin_ac2_p[i][th1]->FindBin(-3*sig_k[i][th1][1]+mean_k[i][th1][1])
 //		     ,hcoin_ac2_p[i][th1]->FindBin(+3*sig_k[i][th1][1]+mean_k[i][th1][1]));
 //------- Get Error Paramters ---//
 n_pi_err[i][th2][1]=fpi[i][th1][1]->GetParError(0); 
 n_k_err[i][th2][1]=fk[i][th1][1]->GetParError(0); 
 n_p_err[i][th2][1]=fp[i][th1][1]->GetParError(0); 

 //----- AC2 Coint Fitting ---------//

 fcoin[i][th1][1] =new TF1(Form("fcoin[%d][%d][1]",i,th1),"gausn(0)+gausn(3)+gausn(6)",min_coin_c,max_coin_c);
 fcoin[i][th1][1]->SetNpx(2000);
 fcoin[i][th1][1]->SetTitle(Form("Coin-Time w AC cut (AC1 Cut<%lf ch && AC2 Cut>%lf ch);Coin time [ns];Counts [1/56 ns]",ac1_adc,th_ac2[i]));
 fcoin[i][th1][1]->SetParameters(n_pi[i][th1][1],mean_pi[i][th1][1],sig_pi[i][th1][1],n_k[i][th1][1],mean_k[i][th1][1],sig_k[i][th1][1],n_p[i][th1][1],mean_p[i][th1][1],sig_p[i][th1][1]);


 sum_k_max=1200.;
  //---- AC1 ----//
 sum_k[i][th2][0]=n_k[i][th2][0]/tdc_time;
 //sum_k[0][th2][0]=1200; 
 
 sum_pi[i][th2][0]=n_pi[i][th2][0]/tdc_time;
 sum_p[i][th2][0]=n_p[i][th2][0]/tdc_time;
 //sum_acc[i][th2][0]=p0_acc[i][th2][0]*(6*sig_k[i][th2][0]/tdc_time);
 sum_acc[i][th2][0]=hcoin_ac1_acc_p[i][th2]->Integral(hcoin_ac1_acc_p[i][th2]->FindBin(-3*sig_k[i][th2][0]+mean_k[i][th2][0])
					              ,hcoin_ac1_acc_p[i][th2]->FindBin(+3*sig_k[i][th2][0]+mean_k[i][th2][0]));
 if( sum_acc[i][th2][0]<0.0)sum_acc[i][th2][0]=0.0;
 sum_k_err[i][th2][0]=n_k_err[i][th2][0]/tdc_time; 
 sum_pi_err[i][th2][0]=n_pi_err[i][th2][0]/tdc_time;
 sum_p_err[i][th2][0]=n_p_err[i][th2][0]/tdc_time;
 //---- AC2 ----//
 sum_k[i][th1][1]=n_k[i][th1][1]/tdc_time; 
 //sum_k[0][th1][1]=1200.;
 sum_pi[i][th1][1]=n_pi[i][th1][1]/tdc_time;
 sum_p[i][th1][1]=n_p[i][th1][1]/tdc_time;
 //sum_acc[i][th1][1]=p0_acc[i][th1][1]*(6*sig_k[i][th1][1]/tdc_time);
 sum_acc[i][th1][1]=hcoin_ac2_acc_p[i][th1]->Integral(hcoin_ac2_acc_p[i][th1]->FindBin(-3*sig_k[i][th1][1]+mean_k[i][th1][1])
			     ,hcoin_ac2_acc_p[i][th1]->FindBin(+3*sig_k[i][th1][1]+mean_k[i][th1][1]));
 if(sum_acc[i][th1][1]<0.0)sum_acc[i][th1][1]=0.0;
 sum_k_err[i][th1][1]=n_k_err[i][th1][1]/tdc_time; 
 sum_pi_err[i][th1][1]=n_pi_err[i][th1][1]/tdc_time;
 sum_p_err[i][th1][1]=n_p_err[i][th1][1]/tdc_time;
 

 //---- AC1 ----//
 if(sum_pi[i][th2][0]<1.05*sum_pi[0][th2][0]){
 gsum_pi_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_pi[i][th2][0]);
 gsum_pi_ac1[th1][th2]->SetPointError(i,emp[i],sum_pi_err[i][th2][0]);}
 if(sum_p[i][th2][0]<1.05*sum_p[0][th2][0]){
 gsum_p_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_p[i][th2][0]);
 gsum_p_ac1[th1][th2]->SetPointError(i,emp[i],sum_p_err[i][th2][0]);}
 if(sum_k[i][th2][0]<1.05*sum_k[0][th2][0]){
 gsum_k_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_k[i][th2][0]);
 gsum_k_ac1[th1][th2]->SetPointError(i,emp[i],sum_k_err[i][th2][0]);
 gSN_k_ac1[th1][th2]->SetPoint(i,th_ac1[i],sum_k[i][th2][0]/sum_acc[i][th2][0]);
 gSN_k_ac1[th1][th2]->SetPointError(i,emp[i],emp[i]);}
 //---- AC2 ----//
 if(sum_pi[i][th1][1]<1.05*sum_pi[0][th1][1]){
 gsum_pi_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_pi[i][th1][1]);
 gsum_pi_ac2[th1][th2]->SetPointError(i,emp[i],sum_pi_err[i][th2][1]);}
 if(sum_p[i][th1][1]<1.05*sum_p[0][th1][1]){
 gsum_p_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_p[i][th1][1]);
 gsum_p_ac2[th1][th2]->SetPointError(i,emp[i],sum_p_err[i][th1][1]);}
 if(sum_k[i][th1][1]<1.05*sum_k[0][th1][1]){
 gsum_k_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_k[i][th1][1]);
 gsum_k_ac2[th1][th2]->SetPointError(i,emp[i],sum_k_err[i][th1][1]);
 gSN_k_ac2[th1][th2]->SetPoint(i,th_ac2[i],sum_k[i][th1][1]/sum_acc[i][th1][1]);
 gSN_k_ac2[th1][th2]->SetPointError(i,emp[i],emp[i]);}
 


 //----------- Rate TGraph ---------------------//
  if(rate_k[i][th2][0]<1.1){
    //rate_k[i][th2][0]=sum_k[i][th2][0]/sum_k[0][th2][0];
   rate_k[i][th2][0]=sum_k[i][th2][0]/sum_k_max;
   rate_k_err[i][th2][0]=sqrt(1./sum_k[i][th2][0]*rate_k[i][th2][0]*(1.-rate_k[i][th2][0]));}
  else{rate_k[i][th2][0]=0.0;
       rate_k_err[i][th2][0]=0.0;}


 if(rate_k[i][th1][1]<1.1){
   //  rate_k[i][th1][1]=sum_k[i][th1][1]/sum_k[0][th1][1];
 rate_k[i][th1][1]=sum_k[i][th1][1]/sum_k_max;
 rate_k_err[i][th1][1]=sqrt(1./sum_k[i][th1][1]*rate_k[i][th1][1]*(1.-rate_k[i][th1][1])); 
 }else{
 rate_k[i][th1][1]=0.0;
 rate_k_err[i][th1][1]=0.0;}

 rate_pi[i][th2][0]=sum_pi[i][th2][0]/sum_pi[0][th2][0];
 rate_p[i][th2][0]=sum_p[i][th2][0]/sum_p[0][th2][0];

 rate_pi[i][th1][1]=sum_pi[i][th1][1]/sum_pi[0][th1][1];
 rate_p[i][th1][1]=sum_p[i][th1][1]/sum_p[0][th1][1];

 rate_p_err[i][th2][0]=sqrt(1./sum_p[i][th2][0]*rate_p[i][th2][0]*(1.-rate_p[i][th2][0]));
 rate_pi_err[i][th2][0]=sqrt(1./sum_pi[i][th2][0]*rate_pi[i][th2][0]*(1.-rate_pi[i][th2][0]));

 rate_p_err[i][th1][1]=sqrt(1./sum_p[i][th1][1]*rate_p[i][th1][1]*(1.-rate_p[i][th1][1]));
 rate_pi_err[i][th1][1]=sqrt(1./sum_pi[i][th1][1]*rate_pi[i][th1][1]*(1.-rate_pi[i][th1][1]));
 //---- AC1 ----//
 if(0.0<rate_k[i][th2][0] && rate_k[i][th2][0]<1.05){
 grate_k_ac1[th1][th2]->SetPoint(i,th_ac1[i],rate_k[i][th2][0]);
 grate_k_ac1[th1][th2]->SetPointError(i,0,rate_k_err[i][th2][0]);
 }else{
 grate_k_ac1[th1][th2]->SetPoint(i,th_ac1[i],0.0);
 grate_k_ac1[th1][th2]->SetPointError(i,0,0.0);
}

 if(rate_p[i][th2][0]<1.05){
 grate_p_ac1[th1][th2]->SetPoint(i,th_ac1[i],rate_p[i][th2][0]);
 grate_p_ac1[th1][th2]->SetPointError(i,0,rate_p_err[i][th2][0]);}
 if(rate_pi[i][th2][0]<1.05){
 grate_pi_ac1[th1][th2]->SetPoint(i,th_ac1[i],rate_pi[i][th2][0]);
 grate_pi_ac1[th1][th2]->SetPointError(i,0,rate_pi_err[i][th2][0]);}



 //---- AC2 ----//

 if(0.0<rate_k[i][th1][1]&& rate_k[i][th1][1]<1.05){
  grate_k_ac2[th1][th2]->SetPoint(i,th_ac2[i],rate_k[i][th1][1]);
  grate_k_ac2[th1][th2]->SetPointError(i,0,rate_k_err[i][th1][1]);
 }else{
  grate_k_ac2[th1][th2]->SetPoint(i,th_ac2[i],0);
  grate_k_ac2[th1][th2]->SetPointError(i,0,0);}

 if(rate_p[i][th1][1]<1.05){
  grate_p_ac2[th1][th2]->SetPoint(i,th_ac2[i],rate_p[i][th1][1]);
  grate_p_ac2[th1][th2]->SetPointError(i,0,rate_p_err[i][th1][1]);}
 if(rate_pi[i][th1][1]<1.05){  
 grate_pi_ac2[th1][th2]->SetPoint(i,th_ac2[i],rate_pi[i][th1][1]);
 grate_pi_ac2[th1][th2]->SetPointError(i,0,rate_pi_err[i][th1][1]);}





//-----------------------------//

 if(n_k[i][th2][0]>max_nk[th1][th2][0])max_nk[th1][th2][0]=n_k[i][th2][0];   if(n_k[i][th1][1]>max_nk[th1][th2][1])max_nk[th1][th2][1]=n_k[i][th1][1];   
 if(n_p[i][th2][0]>max_np[th1][th2][0])max_np[th1][th2][0]=n_p[i][th2][0];   if(n_p[i][th1][1]>max_np[th1][th2][1])max_np[th1][th2][1]=n_p[i][th1][1];
 if(n_pi[i][th2][0]>max_npi[th1][th2][0])max_npi[th1][th2][0]=n_pi[i][th2][0];  if(n_pi[i][th1][1]>max_npi[th1][th2][1])max_npi[th1][th2][1]=n_pi[i][th1][1];
          

 if(max_SN_ac1[th2]<sum_k[i][th2][0]/sum_acc[i][th2][0] && sum_k[i][th2][0]/sum_acc[i][th2][0]<6.0){
   max_SN_ac1[th2]=sum_k[i][th2][0]/sum_acc[i][th2][0];
   SN_ac1[th2]=i;}

 if(max_SN_ac2[th1]<sum_k[i][th1][1]/sum_acc[i][th1][1] && sum_k[i][th1][1]/sum_acc[i][th1][1]<6.0 ){
   max_SN_ac2[th1]=sum_k[i][th1][1]/sum_acc[i][th1][1];
   SN_ac2[th1]=i;}


 // if(FOM_ac1[i][th2]<1.0)
  FOM_ac1[i][th2]=sqrt(pow(sum_k[i][th2][0],2)/sum_acc[i][th2][0]);
// if(FOM_ac2[i][th1]<1.0)
  FOM_ac2[i][th1]=sqrt(pow(sum_k[i][th1][1],2)/sum_acc[i][th1][1]);
  //    th_ac1[i]=th1_max-(th1_max-min_ac1)/iter_ac1*i;
 gfom_ac1[th2]->SetPoint(i,th_ac1[i],FOM_ac1[i][th2]);
 gfom_ac2[th1]->SetPoint(i,th_ac2[i],FOM_ac2[i][th1]);
 gfom->SetPoint(i,th_ac1[i],sqrt(pow(FOM_ac1[i][th2],2)+pow(FOM_ac2[i][th1],2)));
 hfom_ac[th1][th2]=new TH2F(Form("hform_ac[%d][%d]",th1,th2),"",iter_ac1,min_ac1,th1_max,iter_ac1,th2_min,th2_max);
 hfom_ac[th1][th2]->Fill(th_ac1[i],th_ac2[i],sqrt(pow(FOM_ac1[i][th2],2)+pow(FOM_ac2[i][th1],2)));
 if(pra_flag){
   cout<<"FOM_ac1: "<<FOM_ac1[i][th2]<<endl;
   cout<<"FOM_ac1 Max"<<max_fom_ac1<<endl;
   cout<<"FOM_ac2: "<<FOM_ac2[i][th1]<<endl;
   cout<<"FOM_ac2 Max"<<max_fom_ac2<<endl;
}

 if(th1==2 && th2==2){
 if(max_fom_ac1<FOM_ac1[i][th2]){
   max_fom_ac1=FOM_ac1[i][th2];
   fom_th1=i;}
 if(max_fom_ac2<FOM_ac2[i][th1]){
   max_fom_ac2=FOM_ac2[i][th1];
   fom_th2=i;}
 

 }

if(fom_max<sqrt(pow(FOM_ac1[i][th2],2)+pow(FOM_ac2[i][th1],2))){
   fom_max=sqrt(pow(FOM_ac1[i][th2],2)+pow(FOM_ac2[i][th1],2));}

 // cout<<"Max FOM K Events: "<<


     }

	}
     }

        cout<<" Fitting is done !"<<endl;


	//=========== FOMANA ========================//

	TH1F* hcoin_fom=new TH1F("hcoin_fom","hcoin_fom",bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_fom,"","","");
	TH1F* hcoin_acc=new TH1F("hcoin_acc","hcoin_acc",bin_coin_c,min_coin_c,max_coin_c);
	set->SetTH1(hcoin_acc,"","","");
	bool ac2_up,ac2_down,ac2_flag;
	ac2_up=false;
        ac2_down=false;
   if(ac2_min)ac2_up=true;
   if(ac2_min==0)ac2_down=true;
	ev=0;



	for(int k=0;k<evnt;k++){
	  T->GetEntry(k);
   if(k==ev*100000){
     cout<<"Fill Event: "<<k<<"/"<<evnt<<endl; 
     ev=ev+1;
   }
   ac2_flag=false;
   if(ac2_up &&Rvz_cutmin<Rz && Rz<Rvz_cutmax && Lvz_cutmin<Lz && Lz<Lvz_cutmax 
      && Ra1sum<th_ac1[fom_th1] &&th_ac2[fom_th2]<Ra2sum)ac2_flag=true;

   if(ac2_down &&Rvz_cutmin<Rz && Rz<Rvz_cutmax && Lvz_cutmin<Lz && Lz<Lvz_cutmax 
       && Ra1sum<th_ac1[fom_th1] &&th_ac2[fom_th2]>Ra2sum && Ra2sum>th2_min)ac2_flag=true;

   //	  if(Rvz_cutmin<Rz && Rz<Rvz_cutmax && Lvz_cutmin<Lz && Lz<Lvz_cutmax 
   //     && Ra1sum<th_ac1[fom_th1] &&th_ac2[fom_th2]<Ra2sum){

   if(ac2_flag){
	    hcoin_fom->Fill(tcoin_t);
   //------- ACC BG -------------------------------------------------------------//     
  if(mode=="T" &&((-55.<tcoin_t && tcoin_t <-15.) || (5.<tcoin_t && tcoin_t<65.)) ){
       ct=tcoin_t;       
        while(1){
	  if(-15.0<ct && ct<5.0){
		 hcoin_acc->Fill(ct);
                 //hcoin_acc->Fill(ct-20);
		 break;}
	       else if(ct<-15.0){ct=ct+20.0;}
	       else if(5.0<ct){ct=ct-20.;}
	 }
     }
    //----------------------------------------------------------------------------//	
	  }
	}



	TH1F* hcoin_fom_p=new TH1F("hcoin_fom_p","hcoin_fom_p",bin_coin_c,min_coin_c,max_coin_c);
	hcoin_acc->Scale(20./100.);
	hcoin_fom_p->Add(hcoin_fom,hcoin_acc,1.0,-1.0);
	set->SetTH1(hcoin_fom_p,"","","");
        TF1* fk_fom=new TF1("fk_fom","gausn(0)",min_coin_c,max_coin_c);
        fk_fom->SetNpx(2000);
        fk_fom->FixParameter(1,def_mean_k);
        fk_fom->FixParameter(2,def_sig_k);
        //fk_fom->SetParLimits(1,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
	hcoin_fom_p->Fit("fk_fom","","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
	double sumk_fom=fk_fom->GetParameter(0);
	double meank_fom=fk_fom->GetParameter(1);
        double sigk_fom=fk_fom->GetParameter(2);
	sumk_fom=sumk_fom/tdc_time;
	//double sk_fom=sumk_fom;
        double sk_fom=hcoin_fom_p->Integral(hcoin_fom_p->FindBin(-3*sigk_fom+meank_fom)
					    ,hcoin_fom_p->FindBin(3*sigk_fom+meank_fom));

	double sk_fom_ct=hcoin_fom_p->Integral(hcoin_fom_p->FindBin(-1.0)
					    ,hcoin_fom_p->FindBin(1.0));


	  //hcoin_fom_p->Integral(hcoin_fom_p->FindBin(-3*sigk_fom+meank_fom)
	  //	    ,hcoin_fom_p->FindBin(3*sigk_fom+meank_fom));

        double nk_fom=hcoin_acc->Integral(hcoin_acc->FindBin(-3*sigk_fom+meank_fom)
					    ,hcoin_acc->FindBin(3*sigk_fom+meank_fom));
	//nk_fom=nk_fom/tdc_time;
	double snk_fom=sk_fom/nk_fom; 
	double fom=sqrt(snk_fom*sumk_fom);


//======== Draw TCanvas ==============//

TLine* lac=new TLine(min_coin_c,ac1_adc[0],max_coin_c,ac1_adc[0]);
 lac->SetLineWidth(2);
 lac->SetLineColor(2);

 TCanvas* c13=new TCanvas("c13","FOM");
 c13->cd(1);
 set->SetTH2(hfom_ac[2][2],"","","");
 //c13->SetLogz(1);
 hfom_ac[2][2]->Draw("colz");

 TCanvas*c12 =new TCanvas("c12","c12");
 c12->Divide(3,1);
 c12->cd(1);
      hcoin_ac1_p[fom_th1][2]->Draw();
 c12->cd(2);
      hcoin_ac2_p[fom_th2][2]->Draw();
 c12->cd(3);
      hcoin_fom_p->Draw();

 TCanvas*c14 =new TCanvas("c14","c14");
 c14->cd();
      gfom_ac1[2]->SetFillColor(2);
      gfom_ac1[2]->SetMarkerColor(2);
      gfom_ac1[2]->SetFillStyle(3005);
      //gfom_ac1[2]->SetMarkerSize(1.0); 
      gfom_ac1[2]->Draw("AP");
    
 TCanvas*c15 =new TCanvas("c15","c15");
 c15->cd();
      gfom_ac2[2]->SetFillColor(2);
      gfom_ac2[2]->SetMarkerColor(2);
      gfom_ac2[2]->SetFillStyle(3005);
      //gfom_ac2[2]->SetMarkerSize(1.0); 
      gfom_ac2[2]->Draw("AP");
    

      TCanvas* c11=new TCanvas("c11","c11");
      c11->cd();
      hcoin_fom->SetTitle("Coin time AC1 && AC2 Cut; coin time [ns];Counts/56 ns");
      hcoin_fom->SetTitleSize(0.06,"x");
      hcoin_fom->SetTitleSize(0.06,"y");
      hcoin_fom->SetLabelSize(0.04,"x");
      hcoin_fom->SetLabelSize(0.04,"y");
      hcoin_fom->GetXaxis()->SetTitleOffset(0.7);
      hcoin_fom->GetYaxis()->SetTitleOffset(0.7);
      hcoin_fom->GetXaxis()->CenterTitle();
      hcoin_fom->GetYaxis()->CenterTitle();
      hcoin_fom->SetLineColor(1);
      //hcoin_fom->SetStats(111111111);
      hcoin_fom->Draw();
      fk_fom->SetLineColor(4);
      fk_fom->SetFillColor(4);
      fk_fom->SetFillStyle(3001);
      //lac->DrawLine(-3*sigk_fom+meank_fom,0,-3*sigk_fom+meank_fom,1000.);
      //lac->DrawLine(3*sigk_fom+meank_fom,0,3*sigk_fom+meank_fom,1000.);
      fk_fom->Draw("same");


      TCanvas* c9=new TCanvas("c9","Kaon AC1 Efficiency study"); 
   
      c9->Divide(1,3);
      c9->cd(1);
      grate_k_ac1[2][2]->SetMarkerColor(2);
      grate_k_ac1[2][2]->SetFillColor(2);
      grate_k_ac1[2][2]->SetFillStyle(3005);
      grate_k_ac1[2][2]->Draw("Ap");

      c9->cd(2);       
      gSN_k_ac1[2][2]->SetMarkerColor(2);
      gSN_k_ac1[2][2]->SetFillColor(2);
      gSN_k_ac1[2][2]->SetFillStyle(3005);
      gSN_k_ac1[2][2]->SetMinimum(0.0);
      gSN_k_ac1[2][2]->SetMaximum(10.0);
      gSN_k_ac1[2][2]->Draw("AP");
      c9->cd(3);
      //      gfom_ac1[2]->SetLineColor(2);
   
      gfom_ac1[2]->SetFillColor(2);
      gfom_ac1[2]->SetMarkerColor(2);
      gfom_ac1[2]->SetFillStyle(3005);
      //gfom_ac1[2]->SetMarkerSize(1.0); 
      gfom_ac1[2]->Draw("AP");
     




      TCanvas* c10=new TCanvas("c10","Kaon AC2 Efficiency study"); 
   
      c10->Divide(1,3);
      c10->cd(1);
      grate_k_ac2[2][2]->SetMarkerColor(2);
      grate_k_ac2[2][2]->SetFillColor(2);
      grate_k_ac2[2][2]->SetFillStyle(3005);
      grate_k_ac2[2][2]->Draw("Ap");

      c10->cd(2);       
      gSN_k_ac2[2][2]->SetMarkerColor(2);
      gSN_k_ac2[2][2]->SetFillColor(2);
      gSN_k_ac2[2][2]->SetFillStyle(3005);
      gSN_k_ac2[2][2]->SetMinimum(0.0);
      gSN_k_ac2[2][2]->SetMaximum(10.0);
      gSN_k_ac2[2][2]->Draw("AP");
      c10->cd(3);
      //      gfom_ac2[2]->SetLineColor(2);
      gfom_ac2[2]->SetFillColor(2);
      gfom_ac2[2]->SetMarkerColor(2);
      gfom_ac2[2]->SetFillStyle(3005);
      //gfom_ac2[2]->SetMarkerSize(1.0); 
      gfom_ac2[2]->Draw("AP");
     

      TCanvas* c3=new TCanvas("c3","FOM Study Graph");
      c3->Divide(3,2);
      for(int i=0;i<nth;i++){
	//	if(test_flag || pra_flag){break;}
      c3->cd(i+1);
      //gfom_ac1[i]->SetLineColor(2);
      gfom_ac1[i]->SetFillColor(2);
      gfom_ac1[i]->SetFillStyle(3005);
      //gfom_ac1[i]->SetMarkerSize(1.0);
      gfom_ac1[i]->SetMarkerColor(2);
      gfom_ac1[i]->Draw("AP");
      c3->cd(i+4);
      gfom_ac2[i]->SetLineColor(2);
      gfom_ac2[i]->SetFillColor(2);
      gfom_ac2[i]->SetMarkerColor(2);
      gfom_ac2[i]->SetFillStyle(3005);
      //gfom_ac2[i]->SetMarkerSize(1.0); 
      gfom_ac2[i]->Draw("AP");
      }

      

      TCanvas* c2=new TCanvas("c2","S/N Max Coin Hist");
      c2->Divide(3,2);
      for(int i=0;i<nth;i++){
      c2->cd(i+1);
      hcoin_ac1_p[SN_ac1[i]][i]->Draw();
      c2->cd(i+4);
      hcoin_ac2_p[SN_ac2[i]][i]->Draw();
      }


      TCanvas* ck1[nth];
      	for(int th2=0;th2<nth;th2++){
	  //	  if(th2==1){
 ck1[th2]=new TCanvas(Form("ck1[%d]",th2),Form("ck1[%d]",th2)); 
 ck1[th2]->cd();
     fk[0][th2][0]->SetLineColor(4);
     fk[0][th2][0]->SetFillColor(4);
     fk[0][th2][0]->SetFillStyle(3001);
     fpi[0][th2][0]->SetFillStyle(3001);
     fpi[0][th2][0]->SetLineColor(3);
     fpi[0][th2][0]->SetFillColor(3);
     fp[0][th2][0]->SetFillStyle(3001);
     fp[0][th2][0]->SetLineColor(5);
     fp[0][th2][0]->SetFillColor(5);
     hcoin_ac1_p[0][th2]->Draw();
     fk[0][th2][0]->Draw("same");
     //fp[0][th2][0]->Draw("same");
     //fpi[0][th2][0]->Draw("same");
     fcoin[0][th2][0]->SetLineColor(1);
     fcoin[0][th2][0]->Draw("same");
     //	}
      }



     
      TCanvas* ck2[nth];
       for(int th1=0;th1<nth;th1++){
	 //	  if(th1==1){
     ck2[th1]=new TCanvas(Form("ck2[%d]",th1),Form("ck2[%d]",th1)); 
     ck2[th1]->cd();
     fk[0][th1][1]->SetLineColor(4);
     fk[0][th1][1]->SetFillColor(4);
     fk[0][th1][1]->SetFillStyle(3001);
     fpi[0][th1][1]->SetFillStyle(3001);
     fpi[0][th1][1]->SetLineColor(3);
     fpi[0][th1][1]->SetFillColor(3);
     fp[0][th1][1]->SetFillStyle(3001);
     fp[0][th1][1]->SetLineColor(5);
     fp[0][th1][1]->SetFillColor(5);
     hcoin_ac2_p[0][th1]->Draw();
     fk[0][th1][1]->Draw("same");
     //fp[0][th1][1]->Draw("same");
     //fpi[0][th1][1]->Draw("same");
     fcoin[0][th1][1]->SetLineColor(1);
     fcoin[0][th1][1]->Draw("same");
     //	}
      }



       TCanvas* cSN=new TCanvas("cSN","Kaon S/N rate"); 
       cSN->Divide(2,1);
       cSN->cd(1);
       gSN_k_ac1[2][1]->SetMarkerColor(2);
       gSN_k_ac1[2][1]->SetFillColor(2);
       gSN_k_ac1[2][1]->SetFillStyle(3005);
       gSN_k_ac1[2][0]->SetMarkerColor(1);
       gSN_k_ac1[2][0]->SetFillColor(1);
       gSN_k_ac1[2][0]->SetFillStyle(3005);
       gSN_k_ac1[2][2]->SetMarkerColor(3);
       gSN_k_ac1[2][2]->SetFillColor(3);
       gSN_k_ac1[2][2]->SetFillStyle(3005);
      
       gSN_k_ac1[2][2]->SetMinimum(0.0);
       gSN_k_ac1[2][2]->SetMaximum(5.0);
       gSN_k_ac1[2][0]->SetMinimum(0.0);
       gSN_k_ac1[2][0]->SetMaximum(5.0);
       gSN_k_ac1[2][1]->SetMinimum(0.0);
       gSN_k_ac1[2][1]->SetMaximum(5.0);
       gSN_k_ac1[2][2]->Draw("AP");
       //gSN_k_ac1[2][1]->Draw("P");
       //gSN_k_ac1[2][0]->Draw("P");
       cSN->cd(2);
       gSN_k_ac2[1][2]->SetMarkerColor(2);
       gSN_k_ac2[1][2]->SetFillColor(2);
       gSN_k_ac2[1][2]->SetFillStyle(3005);
       gSN_k_ac2[0][2]->SetMarkerColor(1);
       gSN_k_ac2[0][2]->SetFillColor(1);
       gSN_k_ac2[0][2]->SetFillStyle(3005);
       gSN_k_ac2[2][2]->SetMarkerColor(3);
       gSN_k_ac2[2][2]->SetFillColor(3);
       gSN_k_ac2[2][2]->SetFillStyle(3005);
 
       gSN_k_ac2[0][2]->SetMinimum(0.0);
       gSN_k_ac2[0][2]->SetMaximum(10.0);
       gSN_k_ac2[1][2]->SetMinimum(0.0);
       gSN_k_ac2[1][2]->SetMaximum(10.0);
       gSN_k_ac2[2][2]->SetMinimum(0.0);
       gSN_k_ac2[2][2]->SetMaximum(10.0);
       gSN_k_ac2[2][2]->Draw("AP");
       // gSN_k_ac2[1][2]->Draw("P");
       // gSN_k_ac2[0][2]->Draw("P");



 
 TCanvas* ccoin_ac=new TCanvas("ccoin_ac","ccoin_ac");
 ccoin_ac->Divide(1,2);
 ccoin_ac->SetLogz(1);
 ccoin_ac->cd(1); 
 ccoin_ac->cd(1)->SetLogz(1);
 // hcoin_ac1[0]->SetLogz();
 //hcoin_ac1[0]->Draw("colz");
 hcoin_ac1_all->Draw("colz");
 ccoin_ac->cd(2);
ccoin_ac->cd(2)->SetLogz(1);
 //ccoin_ac->SetLogz();
 // hcoin_ac2[0]->SetLogz();
 //hcoin_ac2[0]->Draw("colz");
 hcoin_ac2_all->Draw("colz");

 for(int j=0;j<nth;j++){
   ccoin_ac->cd(1);
 lac->SetLineColor(j+1);
 lac->DrawLine(min_coin_c,ac1_adc[j],max_coin_c,ac1_adc[j]);
 ccoin_ac->cd(2);
 lac->DrawLine(min_coin_c,ac2_adc[j],max_coin_c,ac2_adc[j]); 
 
}

 TCanvas* cAC=new TCanvas("cAC","cAC");
 cAC->cd();
 ha1_a2->Draw("colz");
 for(int j=0;j<nth;j++){
 lac->SetLineColor(j+1);
 lac->DrawLine(min_ac1,ac2_adc[j],max_ac1,ac2_adc[j]);
 lac->DrawLine(ac1_adc[j],min_ac2,ac1_adc[j],max_ac2);
 }



 lac->SetLineColor(2);
 

  TCanvas* ccoin_t1=new TCanvas("ccoin_t1","ccoin_t1");
  ccoin_t1->Divide(3,2);    
  int mv1=1;

    for(int th2=0;th2<nth;th2++){
     ccoin_t1->cd(mv1);
      hcoin_ac1_max[th2]->Draw();
      mv1=mv1+1;
  }
    for(int th1=0;th1<nth;th1++){
      ccoin_t1->cd(mv1);
      hcoin_ac2_max[th1]->Draw();
      mv1=mv1+1;
    }

 

      
     TCanvas* cAC1_cut=new TCanvas("cAC1_cut","cAC1_cut");  
    cAC1_cut->Divide(4,3);
    int l;
    for(int i=0;i<12;i++){
     for(int j=0;j<nth;j++){
       if(test_flag || pra_flag){break;}
       else{l=4*i;}
       if(j==1){
	 cAC1_cut->cd(i+1);
     fk[l][j][0]->SetLineColor(4);
     fk[l][j][0]->SetFillColor(4);
     fk[l][j][0]->SetFillStyle(3001);
     fpi[l][j][0]->SetFillStyle(3001);
     fpi[l][j][0]->SetLineColor(46);
     fpi[l][j][0]->SetFillColor(46);
     fp[l][j][0]->SetFillStyle(3001);
     fp[l][j][0]->SetLineColor(8);
     fp[l][j][0]->SetFillColor(8);
     hcoin_ac1_p[l][j]->Draw();
     fk[l][j][0]->Draw("same");
     fp[l][j][0]->Draw("same");
     fpi[l][j][0]->Draw("same");
     fcoin[l][j][0]->Draw("same");
      
     }
    }
     }
     





 

  TCanvas* cAC2_cut=new TCanvas("cAC2_cut","cAC2_cut"); 
     cAC2_cut->Divide(2,3);
     
  for(int i=0;i<6;i++){
     for(int j=0;j<nth;j++){
       if(j==1){
       if(test_flag || pra_flag){break;}
       //else{l=4*i;}
       else{l=i;}
     cAC2_cut->cd(i+1);
     hcoin_ac2_p[l][j]->Draw();
     fk[l][j][1]->SetLineColor(4);
     fk[l][j][1]->SetFillColor(4);
     fk[l][j][1]->SetFillStyle(3001);
     fpi[l][j][1]->SetFillStyle(3001);
     fpi[l][j][1]->SetLineColor(46);
     fpi[l][j][1]->SetFillColor(46);
     fp[l][j][1]->SetFillStyle(3001);
     fp[l][j][1]->SetLineColor(8);
     fp[l][j][1]->SetFillColor(8);
     fk[l][j][1]->Draw("same");
     fp[l][j][1]->Draw("same");
     fpi[l][j][1]->Draw("same");
     fcoin[l][j][1]->Draw("same");
     

       }
     }
     
  }
     


  TCanvas* c8=new TCanvas("c8","AC2 Cut Hist");
  c8->Divide(2,2);
  for(int l=0;l<4;l++){
   if(test_flag || pra_flag){break;}
    c8->cd(l+1);
     hcoin_ac2_p[l][1]->Draw();
     fk[l][1][1]->SetLineColor(4);
     fk[l][1][1]->SetFillColor(4);
     fk[l][1][1]->SetFillStyle(3001);
     fpi[l][1][1]->SetFillStyle(3001);
     fpi[l][1][1]->SetLineColor(46);
     fpi[l][1][1]->SetFillColor(46);
     fp[l][1][1]->SetFillStyle(3001);
     fp[l][1][1]->SetLineColor(8);
     fp[l][1][1]->SetFillColor(8);
     fk[l][1][1]->Draw("same");
     fp[l][1][1]->Draw("same");
     fpi[l][1][1]->Draw("same");
     fcoin[l][1][1]->Draw("same");

    }
    
  
     TCanvas* cAC1_num=new TCanvas("cAC1_num","cAC1_num");
      cAC1_num->Divide(1,3);
     TCanvas* cAC1_rate=new TCanvas("cAC1_rate","cAC1_rate");
          cAC1_rate->Divide(1,3);
      	  TCanvas* c6=new TCanvas("c6","AC K rate");
	  c6->Divide(2,1);
       for(int k=0;k<nth;k++){
    for(int j=0;j<nth;j++){
      if(j==1){
       if(k==0){
   cAC1_num->cd(1);
   gsum_k_ac1[j][k]->SetMarkerColor(1);
   gsum_k_ac1[j][k]->SetFillColor(1);
   gsum_k_ac1[j][k]->SetFillStyle(3005);
   gsum_k_ac1[j][k]->Draw("Ap");
   cAC1_num->cd(2);
   gsum_pi_ac1[j][k]->SetFillColor(1);
   gsum_pi_ac1[j][k]->SetFillStyle(3005);
   gsum_pi_ac1[j][k]->SetMarkerColor(1);
   gsum_pi_ac1[j][k]->Draw("Ap");
   cAC1_num->cd(3);
   gsum_p_ac1[j][k]->SetFillColor(1);
   gsum_p_ac1[j][k]->SetFillStyle(3005);
   gsum_p_ac1[j][k]->SetMarkerColor(1);
   gsum_p_ac1[j][k]->Draw("Ap");
   
    cAC1_rate->cd(1);
    //grate_k_ac1[j][k]->SetFillColor(1);
   // grate_k_ac1[j][k]->SetFillStyle(3005);
    grate_k_ac1[j][k]->GetYaxis()->SetRangeUser(0.0,0.5);
   grate_k_ac1[j][k]->SetMarkerColor(1);
   grate_k_ac1[j][k]->Draw("Ap");
   c6->cd(1);
   grate_k_ac1[j][k]->Draw("AP");
   cAC1_rate->cd(2);
   // grate_pi_ac1[j][k]->SetFillColor(1);
   // grate_pi_ac1[j][k]->SetFillStyle(3005);
   grate_pi_ac1[j][k]->SetMarkerColor(1); 
   grate_pi_ac1[j][k]->Draw("Ap");
   cAC1_rate->cd(3);
   // grate_k_ac1[j][k]->SetFillColor(1);
   // grate_k_ac1[j][k]->SetFillStyle(3005);
   grate_p_ac1[j][k]->SetMarkerColor(1);
   grate_p_ac1[j][k]->Draw("Ap");
     
   
      }else{
   cAC1_num->cd(1);
   gsum_k_ac1[j][k]->SetMarkerColor(k+1);
   gsum_k_ac1[j][k]->Draw("P");   
   cAC1_num->cd(2);   
   gsum_pi_ac1[j][k]->SetMarkerColor(k+1);
   gsum_pi_ac1[j][k]->Draw("P");
   cAC1_num->cd(3);   
   gsum_p_ac1[j][k]->SetMarkerColor(k+1);
   gsum_p_ac1[j][k]->Draw("P");
   
   cAC1_rate->cd(1);
   grate_k_ac1[j][k]->SetMarkerColor(k+1);
   grate_k_ac1[j][k]->Draw("P");
   c6->cd(1);
   grate_k_ac1[j][k]->Draw("P");
   cAC1_rate->cd(2);   
   grate_pi_ac1[j][k]->SetMarkerColor(k+1);
   grate_pi_ac1[j][k]->Draw("P");
   cAC1_rate->cd(3);   
   grate_p_ac1[j][k]->SetMarkerColor(k+1);
   grate_p_ac1[j][k]->Draw("P");
    }
       }
    }
       }




       TCanvas* cAC_rate=new TCanvas("cAC_rate","cAC_rate");
      
       int x=2; int y=2;

       TGraphErrors* grate_k_ac1_cp=(TGraphErrors*)grate_k_ac1[x][y]->Clone();
       TGraphErrors* grate_pi_ac1_cp=(TGraphErrors*)grate_pi_ac1[x][y]->Clone();
       TGraphErrors* grate_p_ac1_cp=(TGraphErrors*)grate_p_ac1[x][y]->Clone();
       TGraphErrors* grate_k_ac2_cp=(TGraphErrors*)grate_k_ac2[x][y]->Clone();
       TGraphErrors* grate_pi_ac2_cp=(TGraphErrors*)grate_pi_ac2[x][y]->Clone();
       TGraphErrors* grate_p_ac2_cp=(TGraphErrors*)grate_p_ac2[x][y]->Clone();

       cAC_rate->Divide(2,1);
       cAC_rate->cd(1);
       grate_k_ac1_cp->SetMinimum(0.0);
       grate_k_ac1_cp->SetMaximum(1.05);
       grate_pi_ac1_cp->SetMinimum(0.0);
       grate_pi_ac1_cp->SetMaximum(1.05);
       grate_p_ac1_cp->SetMinimum(0.0);
       grate_p_ac1_cp->SetMaximum(1.05);
       
       grate_k_ac1_cp->SetMarkerColor(2);
       grate_pi_ac1_cp->SetMarkerColor(3);
       grate_p_ac1_cp->SetMarkerColor(4);
       grate_pi_ac1_cp->Draw("AP");
       grate_k_ac1_cp->Draw("P");   
       grate_p_ac1_cp->Draw("P");
       cAC_rate->cd(2);
       grate_k_ac2_cp->SetMinimum(0.0);
       grate_k_ac2_cp->SetMaximum(1.05);
       grate_pi_ac2_cp->SetMinimum(0.0);
       grate_pi_ac2_cp->SetMaximum(1.05);
       grate_p_ac2_cp->SetMinimum(0.0);
       grate_p_ac2_cp->SetMaximum(1.05);
       grate_k_ac2_cp->SetMarkerColor(2);
       grate_pi_ac2_cp->SetMarkerColor(3);
       grate_p_ac2_cp->SetMarkerColor(4);
       grate_pi_ac2_cp->Draw("AP");
       grate_k_ac2_cp->Draw("P");   
       grate_p_ac2_cp->Draw("P");
       /*
       cAC_rate->Divide(1,2);
       cAC_rate->cd(1);
       grate_k_ac1[1][1]->SetMinimum(0.0);
       grate_k_ac1[1][1]->SetMaximum(1.2);
       grate_pi_ac1[1][1]->SetMinimum(0.0);
       grate_pi_ac1[1][1]->SetMaximum(1.2);
       grate_k_ac1[1][1]->SetMarkerColor(2);
       grate_k_ac1[1][1]->Draw("AP");   
       grate_pi_ac1[1][1]->SetMarkerColor(3);
       grate_pi_ac1[1][1]->Draw("P");
       cAC_rate->cd(2);
       grate_k_ac2[1][1]->SetMinimum(0.0);
       grate_k_ac2[1][1]->SetMaximum(1.1);
       grate_pi_ac2[1][1]->SetMinimum(0.0);
       grate_pi_ac2[1][1]->SetMaximum(1.1);
       grate_k_ac2[1][1]->SetMarkerColor(2);
       grate_k_ac2[1][1]->Draw("AP");   
       grate_pi_ac2[1][1]->SetMarkerColor(3);
       grate_pi_ac2[1][1]->Draw("P");	
       */



  TCanvas* c7=new TCanvas("c7","CAC2 Kaon Survival Rrate");
   TGraphErrors* grate_k_ac2_clone[nth][nth];

 TCanvas* cAC2_rate=new TCanvas("cAC2_rate","cAC2_rate");
          cAC2_rate->Divide(1,3);
 TCanvas* cAC2_num=new TCanvas("cAC2_num","cAC2_num");
          cAC2_num->Divide(1,3);
	         for(int j=0;j<nth;j++){
       for(int k=0;k<nth;k++){
	 if(k==1){
	  if(j==0){
   cAC2_num->cd(1);
   gsum_k_ac2[j][k]->SetMarkerColor(1);
   gsum_k_ac2[j][k]->Draw("AP");
   cAC2_num->cd(2);
   gsum_pi_ac2[j][k]->SetMarkerColor(1);
   gsum_pi_ac2[j][k]->Draw("AP");
   cAC2_num->cd(3);
   gsum_p_ac2[j][k]->SetMarkerColor(1);
   gsum_p_ac2[j][k]->Draw("AP");

   cAC2_rate->cd(1);
   grate_k_ac2[j][k]->SetMarkerColor(1);
   grate_k_ac2[j][k]->Draw("AP");
   c6->cd(2); 
   grate_k_ac2[j][k]->Draw("AP");
   c7->cd();
   grate_k_ac2_clone[j][k]=(TGraphErrors*)grate_k_ac2[j][k]->Clone();
   grate_k_ac2_clone[j][k]->SetMarkerColor(1);
   grate_k_ac2_clone[j][k]->GetXaxis()->SetRangeUser(0.0,5.0);
   grate_k_ac2_clone[j][k]->GetYaxis()->SetRangeUser(0.5,1.2);
   grate_k_ac2_clone[j][k]->SetMarkerSize(1.5);
   grate_k_ac2_clone[j][k]->Draw("AP");
   cAC2_rate->cd(2);
   grate_pi_ac2[j][k]->SetMarkerColor(1);
   grate_pi_ac2[j][k]->Draw("AP");
   cAC2_rate->cd(3);
   grate_p_ac2[j][k]->SetMarkerColor(1);
   grate_p_ac2[j][k]->Draw("AP");
	  }else{
   cAC2_num->cd(1);
   gsum_k_ac2[j][k]->SetMarkerColor(j+1);
   gsum_k_ac2[j][k]->Draw("P");
   cAC2_num->cd(2);
   gsum_pi_ac2[j][k]->SetMarkerColor(j+1);
   gsum_pi_ac2[j][k]->Draw("P");
   cAC2_num->cd(3);
   gsum_p_ac2[j][k]->SetMarkerColor(j+1);
   gsum_p_ac2[j][k]->Draw("P");

   cAC2_rate->cd(1);
   grate_k_ac2[j][k]->SetMarkerColor(j+1);
   grate_k_ac2[j][k]->Draw("P");
   c6->cd(2); 
   grate_k_ac2[j][k]->Draw("P");
   c7->cd();
   grate_k_ac2_clone[j][k]=(TGraphErrors*)grate_k_ac2[j][k]->Clone();
   grate_k_ac2_clone[j][k]->SetMarkerColor(j+1);
   grate_k_ac2_clone[j][k]->SetMarkerSize(1.5);
   grate_k_ac2_clone[j][k]->Draw("P");
   cAC2_rate->cd(2);
   grate_pi_ac2[j][k]->SetMarkerColor(j+1);
   grate_pi_ac2[j][k]->Draw("P");
   cAC2_rate->cd(3);
   grate_p_ac2[j][k]->SetMarkerColor(j+1);
   grate_p_ac2[j][k]->Draw("P");

	  }
	 }
       }
        }



   TCanvas* ck_ac=new TCanvas("ck_ac","ck_ac");
   ck_ac->Divide(1,2);
   ck_ac->cd(1);
   gsum_k_ac1[1][1]->Draw("AP");
   gsum_pi_ac1[1][1]->Draw("P");
		 ck_ac->cd(2);
   gsum_k_ac2[1][1]->Draw("AP");
   gsum_pi_ac2[1][1]->Draw("P");

   
   TCanvas* cx=new TCanvas("cx","Coin time");
   cx->Divide(1,2); 
   cx->cd(1);
   hcoin_t->Draw();
   cx->cd(2);
   hcoin_tc->Draw();

   TCanvas* ck_def=new TCanvas("ck_def","Kaon Cut Hist");
   ck_def->Divide(1,3);
   ck_def->cd(1);
     fk_kc->SetLineColor(4);
     fk_kc->SetFillColor(4);
     fk_kc->SetFillStyle(3001);
     /*
     fpi_kc->SetFillStyle(3001);
     fpi_kc->SetLineColor(46);
     fpi_kc->SetFillColor(46);
     fp_kc->SetFillStyle(3001);
     fp_kc->SetLineColor(8);
     fp_kc->SetFillColor(8);
     */
   hcoin_k->Draw();
   fk_kc->Draw("same");
   // fp_kc->Draw("same");
   // fpi_kc->Draw("same");
    facc_kc->Draw("same");
   // fcoin_kc->Draw("same");
  
  ck_def->cd(2);
  hcoin_pi->Draw();
  fpi_pic->SetFillStyle(3001);
  fpi_pic->SetLineColor(46);
  fpi_pic->SetFillColor(46);
  fpi_pic->Draw("same");

  ck_def->cd(3);
  hcoin_p->Draw();
  fp_pc->SetFillStyle(3001);
  fp_pc->SetLineColor(46);
  fp_pc->SetFillColor(46);
  fp_pc->Draw("same");


  TCanvas *c4=new TCanvas("c4","Coincidence ACC Hist");
  c4->Divide(2,3); 
  for(int k=0;k<nth;k++){
    c4->cd(k+1);
    hcoin_t1[k]->GetXaxis()->SetRangeUser(-5.0,5.0);
    hcoin_t1[k]->Draw();
    hcoin_acc_ac1[k]->SetLineColor(2);
    hcoin_acc_ac1[k]->SetFillColor(2);
    hcoin_acc_ac1[k]->SetFillStyle(3001);
    hcoin_acc_ac1[k]->Draw("same");
    c4->cd(k+4);
    hcoin_t2[k]->GetXaxis()->SetRangeUser(-5.0,5.0);
    hcoin_t2[k]->Draw();
    hcoin_acc_ac2[k]->SetLineColor(2);
    hcoin_acc_ac2[k]->SetFillColor(2);
    hcoin_acc_ac2[k]->SetFillStyle(3001);
    hcoin_acc_ac2[k]->Draw("same");
  }

  TCanvas *c5=new TCanvas("c5","Coin Hist");
  c5->Divide(2,3);
  for(int i=0;i<nth;i++){
 c5->cd(i+1);
 hcoin_ac1_all_p[0][i]->SetLineColor(3);
 hcoin_ac1_all_p[0][i] ->SetFillColor(3);
 hcoin_ac1_all_p[0][i] ->SetFillStyle(3001);
 hcoin_ac1_acc_p[0][i]->SetLineColor(4);
 hcoin_ac1_acc_p[0][i] ->SetFillColor(4);
 hcoin_ac1_acc_p[0][i] ->SetFillStyle(3001);
 hcoin_ac1_p[0][i]->SetLineColor(2);
 hcoin_ac1_p[0][i] ->SetFillColor(2);
 hcoin_ac1_p[0][i] ->SetFillStyle(3001);

  hcoin_ac1_all_p[0][i]->Draw();
  hcoin_ac1_acc_p[0][i]->Draw("same");
  hcoin_ac1_p[0][i]->Draw("same");
  c5->cd(i+4);
 hcoin_ac2_all_p[0][i]->SetLineColor(3);
 hcoin_ac2_all_p[0][i] ->SetFillColor(3);
 hcoin_ac2_all_p[0][i] ->SetFillStyle(3001);
 hcoin_ac2_acc_p[0][i]->SetLineColor(4);
 hcoin_ac2_acc_p[0][i] ->SetFillColor(4);
 hcoin_ac2_acc_p[0][i] ->SetFillStyle(3001);
 hcoin_ac2_p[0][i]->SetLineColor(2);
 hcoin_ac2_p[0][i] ->SetFillColor(2);
 hcoin_ac2_p[0][i] ->SetFillStyle(3001);
  hcoin_ac2_all_p[0][i]->Draw();
  hcoin_ac2_acc_p[0][i]->Draw("same");
  hcoin_ac2_p[0][i]->Draw("same");
  }



      cout<<"drawing is done !!"<<endl;       

 //================ Print Canvas =================================//

        
 TString name;
 if(output_flag){
 
 name.Form(ofname.c_str());
 ccoin_ac->Print(name+"[","pdf");

ccoin_ac->Print(name,"pdf");
// cAC_rate->->Print(name,"pdf"); 
cAC->Print(name,"pdf"); 
 ck_def->Print(name,"pdf");
 ck_ac->Print(name,"pdf"); 
 ccoin_t1->Print(name,"pdf"); 
 


 for(int j=0;j<nth;j++){
ck1[j]->Print(name,"pdf");
ck2[j]->Print(name,"pdf");
}

 cSN->Print(name,"pdf");
 c2->Print(name,"pdf");
 c3->Print(name,"pdf");
 c4->Print(name,"pdf");
 c5->Print(name,"pdf");
 c6->Print(name,"pdf");
 c7->Print(name,"pdf");
 c8->Print(name,"pdf");
 c9->Print(name,"pdf");
 c10->Print(name,"pdf");
 c11->Print(name,"pdf");
 c12->Print(name,"pdf");
 c14->Print(name,"pdf");
 c15->Print(name,"pdf");
 // c13->Print(name,"pdf");
 cAC1_cut->Print(name,"pdf");
 cAC2_cut->Print(name,"pdf");
 cAC_rate->Print(name,"pdf");
 cAC1_num->Print(name,"pdf");
 cAC2_num->Print(name,"pdf");
 cAC1_rate->Print(name,"pdf");
 cAC2_rate->Print(name,"pdf");
 cAC2_rate->Print(name+"]","pdf");

 }
    
 cout<<"Print is done "<<endl;
   


 //======= COMMENT OUT ==================//

 for(int j=0;j<nth;j++){
   for(int k=0;k<nth;k++){

 cout<<"========= j:"<<j<<" k:"<<k<<"========================"<<endl; 

 cout<<Form("ac1_adc[%d]:",j)<<ac1_adc[j]<<Form(" : ac2_adc[%d]:",k)<<ac2_adc[k]<<endl;
 cout<<"max k[ac1]: "<<max_nk[j][k][0]/tdc_time<<" : max k[ac2]: "<<max_nk[j][k][1]/tdc_time<<endl;

 cout<<"========================================================="<<endl;
   }
 }

 cout<<"Kaon Coincidence Get Parameters"<<endl;
 cout<<"def_num_k "<<def_num_k<<endl;
 cout<<"def_sig_k "<<def_sig_k<<endl;
 cout<<"def_mean_k "<<def_mean_k<<endl;
 cout<<"Pion Coincidence Get Parameters"<<endl;
 cout<<"def_num_pi "<<def_num_pi<<endl;
 cout<<"def_sig_pi "<<def_sig_pi<<endl;
 cout<<"def_mean_pi "<<def_mean_pi<<endl;
 cout<<"Proton Coincidence Get Parameters"<<endl;
 cout<<"def_num_p "<<def_num_p<<endl;
 cout<<"def_sig_p "<<def_sig_p<<endl;
 cout<<"def_mean_p "<<def_mean_p<<endl;





	cout<<"======   Maximum FOM   ===="<<endl;
	cout<<"===== AC1 ====="<<endl;
	cout<<"AC1 < "<<th_ac1[fom_th1]<<endl;
	cout<< ac2_adc[2]<<" < AC2 < "<<th2_max<<endl;
	cout<<"AC1 Kaon : "<<sum_k[fom_th1][2][0]<<endl;
	cout<<"AC1 ACC : "<<sum_acc[fom_th1][2][0]<<endl;
	cout<<"AC1 Survival Rate: "<<rate_k[fom_th1][2][0]<<endl;
	cout<<"AC1 S/N ratio: "<<sum_k[fom_th1][2][0]/sum_acc[fom_th1][2][0]<<endl;
	cout<<"AC1 Max FOM: "<<max_fom_ac1<<endl;
	cout<<"===== AC2 ====="<<endl;
	cout<<"AC1 <"<< ac1_adc[2]<<endl;
	cout<<th_ac2[fom_th2]<<" < AC2 < "<<th2_max<<endl;
        cout<<"AC2 Kaon : "<<sum_k[fom_th2][2][1]<<endl;
	cout<<"AC2 ACC : "<<sum_acc[fom_th2][2][1]<<endl;
	cout<<"AC2 Survival Rate: "<<rate_k[fom_th2][2][1]<<endl;
	cout<<"AC2 S/N ratio: "<<sum_k[fom_th2][2][1]/sum_acc[fom_th2][2][1]<<endl;
        cout<<"AC2 Max FOM: "<<max_fom_ac2<<endl;
	cout<<"===== AC1 & AC2 ====="<<endl;
	//cout<<"AC1 < "<<th_ac1[fom_th1]<<endl;
	//cout<<th_ac2[fom_th2]<<" < AC2 < "<<th2_max<<endl;
	cout<<"Kaon events: "<<sumk_fom<<endl;
	cout<<"Kaon Integral:"<<sk_fom<<endl;
	cout<<"Kaon ACC: "<<nk_fom<<endl;
	cout<<"Kaon survival rate: "<<sumk_fom/sum_k_max<<endl;
	cout<<"S/N ratio : "<<snk_fom<<endl;
	cout<<"FOM : "<<fom<<endl;
	cout<<"FOM MAX: "<<fom_max<<endl;
	cout<<"Kaon (-1.0<ct<1.0): "<<sk_fom_ct<<endl;
	// double nk_ct=h
 

 if(draw_flag==0)gSystem->Exit(1);
 theApp->Run();
 return 0;

}


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

double  RS2_off_H1[16]={-16828.7,-16863,-16894,-16893.3,-16870.9,-16867.2,-16900.3,-16876.8,17554.1,-16861.6,-16895,-16890.7,-16854.6,-16852.9,-16850.5,-16861.9};
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
