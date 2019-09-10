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
  bool root_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  bool test_flag=false;
  bool Pion_flag=false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"h:f:w:n:r:bcop:GHTt12P"))!=-1){
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

    case 'r':
      root_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
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
  // if(draw_flag==0)
   gROOT->SetBatch(1);


 Setting *set = new Setting();
 set->Initialize();


 //TChain //
  TChain* T;
  if(mode=="G"|| mode=="T"){T=new TChain("tree"); }
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
    //    cout<<buf<<endl;
  }


  cout<<"mode :"<<mode<<endl;
  cout<<"tdc_time[ns]: "<<tdc_time<<endl;
  cout<<"Pion Flag:"<<Pion_flag<<endl;
  int evnt=T->GetEntries();
  if(test_flag){evnt=1000;}
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
  double Ls2lt[100], Ls2rt[100];
  double Ls2rtc_fadc[100],Ls2ltc_fadc[100];
  double Rs2lt[100], Rs2rt[100];
  double Rs2rtc_fadc[100],Rs2ltc_fadc[100];
  //---- Gogami root ---------//
  double ctime[1000];
  double DRT5;
  double Rs2ra[100];
  double Rs2la[100];
  double Ls2ra[100];
  double Ls2la[100];
  double Ls2la_c[100];
  double Ls2ra_c[100];
  T->SetBranchStatus("*",0);  
 //------ Right Arm -------------//

 T->SetBranchStatus("DR.evtype",1);
 T->SetBranchAddress("DR.evtype",&DRevtype); 
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
 T->SetBranchStatus("R.s2.lt_c",1);
 T->SetBranchAddress("R.s2.lt_c",Rs2l_tc);
 T->SetBranchStatus("R.s2.rt_c",1);
 T->SetBranchAddress("R.s2.rt_c",Rs2r_tc);
 //--- FBUS TDC --------//
 T->SetBranchStatus("R.s2.lt",1);
 T->SetBranchAddress("R.s2.lt",Rs2lt);
 T->SetBranchStatus("R.s2.rt",1);
 T->SetBranchAddress("R.s2.rt",Rs2rt);
 
 T->SetBranchStatus("R.s2.lt_c",1);
 T->SetBranchAddress("R.s2.lt_c",Rs2l_tc); 
 T->SetBranchStatus("R.s2.rt_c",1);
 T->SetBranchAddress("R.s2.rt_c",Rs2r_tc);
 
 //---FADC TDC --------//
 T->SetBranchStatus("R.s2.rtc_fadc",1);
 T->SetBranchAddress("R.s2.rtc_fadc",Rs2rtc_fadc);
 T->SetBranchStatus("R.s2.ltc_fadc",1);
 T->SetBranchAddress("R.s2.ltc_fadc",Rs2ltc_fadc);



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
 T->SetBranchStatus("L.s2.ra_c",1);
 T->SetBranchAddress("L.s2.ra_c",Ls2ra_c);
 T->SetBranchStatus("L.s2.la_c",1);
 T->SetBranchAddress("L.s2.la_c",Ls2la_c);

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
 //--- FBUS TDC --------//
 T->SetBranchStatus("L.s2.lt",1);
 T->SetBranchAddress("L.s2.lt",Ls2lt);
 T->SetBranchStatus("L.s2.rt",1);
 T->SetBranchAddress("L.s2.rt",Ls2rt);
 
 T->SetBranchStatus("L.s2.lt_c",1);
 T->SetBranchAddress("L.s2.lt_c",Ls2l_tc); 
 T->SetBranchStatus("L.s2.rt_c",1);
 T->SetBranchAddress("L.s2.rt_c",Ls2r_tc);
 
 //---FADC TDC --------//
 T->SetBranchStatus("L.s2.rtc_fadc",1);
 T->SetBranchAddress("L.s2.rtc_fadc",Ls2rtc_fadc);
 T->SetBranchStatus("L.s2.ltc_fadc",1);
 T->SetBranchAddress("L.s2.ltc_fadc",Ls2ltc_fadc);
 
 if(mode=="G"){
 T->SetBranchStatus("ctime",1);    
 T->SetBranchAddress("ctime",ctime);
 T->SetBranchStatus("DR.T5",1);    
 T->SetBranchAddress("DR.T5",&DRT5);

}



 //========== New Tree ==============================//
 TFile* fnew =new TFile("%s",ofname.c_str(),"recreate");
 TTree* tnew=new TTree("T","Coin time events shift hist");
 double f1_s2r;//in real this is left PMT
 double f1_s2l;
 double f1_Rs2r[16],f1_Rs2l[16],f1_Ls2r[16],f1_Ls2l[16];
 double fbus_Rs2r[16],fbus_Rs2l[16],fbus_Ls2r[16],fbus_Ls2l[16];
 double fadc_Rs2r[16],fadc_Rs2l[16],fadc_Ls2r[16],fadc_Ls2l[16];
 double coin_f1,coin_fadc,coin_fbus;
 double coin_f1_c,coin_fadc_c,coin_fbus_c;
 int event,trig;
 int Ls2pads,Rs2pads;

 tnew->Branch("nev",&event);
   tnew->Branch("trig",&trig);
   tnew->Branch("Rs2_pads",&Rs2pads);
   tnew->Branch("Ls2_pads",&Ls2pads);
   
 //-------- F1TDC -----------//
  tnew->Branch("Ls2_f1_rt", f1_Ls2r,"Ls2_f1_rt[16]/D");
  tnew->Branch("Ls2_f1_lt", f1_Ls2l,"Ls2_f1_lt[16]/D");
  tnew->Branch("Rs2_f1_rt", f1_Rs2r,"Rs2_f1_rt[16]/D");
  tnew->Branch("Rs2_f1_lt", f1_Rs2l,"Rs2_f1_lt[16]/D");
  tnew->Branch("coin_f1", &coin_f1,"coin_f1/D");
  tnew->Branch("coin_f1_c", &coin_f1_c,"coin_f1_c/D");  
  //------ FBUS -------------//
  tnew->Branch("Ls2_fbus_rt", fbus_Ls2r,"Ls2_fbus_rt[16]/D");
  tnew->Branch("Ls2_fbus_lt", fbus_Ls2l,"Ls2_fbus_lt[16]/D");
  tnew->Branch("Rs2_fbus_rt", fbus_Rs2r,"Rs2_fbus_rt[16]/D");
  tnew->Branch("Rs2_fbus_lt", fbus_Rs2l,"Rs2_fbus_lt[16]/D");
  tnew->Branch("coin_fbus", &coin_fbus,"coin_fbus/D");
  tnew->Branch("coin_fbus_c", &coin_fbus_c,"coin_fbus_c/D");   
  //------- FADC -------------//
  tnew->Branch("Ls2_fadc_rt", fadc_Ls2r,"Ls2_fadc_rt[16]/D");
  tnew->Branch("Ls2_fadc_lt", fadc_Ls2l,"Ls2_fadc_lt[16]/D");
  tnew->Branch("Rs2_fadc_rt", fadc_Rs2r,"Rs2_fadc_rt[16]/D");
  tnew->Branch("Rs2_fadc_lt", fadc_Rs2l,"Rs2_fadc_lt[16]/D");
  tnew->Branch("coin_fadc", &coin_fadc,"coin_fadc/D");
  tnew->Branch("coin_fadc_c", &coin_fadc_c,"coin_fadc_c/D");   
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
 int bin_ac1=(max_ac1-min_ac1); 
 int bin_ac2=(max_ac2-min_ac2); 

 TH2F* hcoin_ac1[nth];
 TH2F* hcoin_ac2[nth];
 TH1F* hcoin_t1[nth];
 TH1F* hcoin_t2[nth];
 TH1F* hcoin_ac1_max[nth];
 TH1F* hcoin_ac2_max[nth];
 TH1F* hcoin_t3[nth][nth];
 TH1F* hcoin_t=new TH1F("hcoin_t","Coincidence time S2R-S2L[ns] ",bin_coin,min_coin,max_coin);
 TH1F* hcoin_tc=new TH1F("hcoin_tc","Coincidence time w/ Path Length Correction  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);


 TH2F* ha1_a2=new TH2F("ha1_a2","ac1 vs ac2 ADC sum hist",bin_ac1,min_ac1,max_ac1,bin_ac2,min_ac2,max_ac2);

 TH1F* hcoin_k=new TH1F("hcoin_k","Coincidence time w/ Correction Kaon Cut  S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);

 TH1F* hcoin_pi=new TH1F("hcoin_pi","Coincidence time w/ Correction Pion  Cut S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_pi_cut=new TH1F("hcoin_pi_cut","Coincidence time w/ Shift cut",bin_coin_c,min_coin_c,max_coin_c);
 TH1F* hcoin_pi_cut_all= new TH1F("hcoin_pi_cut_all","Coincidence time w/ Shift cut all segments",bin_coin_c,min_coin_c,max_coin_c);
  
 TH1F* hcoin_p=new TH1F("hcoin_p","Coincidence time w/ Correction Proton  Cut S2R-S2L[ns] ",bin_coin_c,min_coin_c,max_coin_c);


 
 TH2F* hcoin_Rs2r[16];
 TH2F* hcoin_Rs2l[16];
 TH2F* hcoin_Ls2r[16];
 TH2F* hcoin_Ls2l[16];
 TH2F* hcoin_Rs2[16];
 TH2F* hLs2r_trig[16]; 
 //TH1F* hcoin_R=new TH1F("hcoin_R","",1000,-750,-700);
 TH1F* hcoin_R=new TH1F("hcoin_R","",1000,-100,0); 
 TH1F* hcoin_2=new TH1F("hcoin_2","",1000,-20,20);
 TH2F* hLs2_y[16];
 TH2F* hLs2[16];
 TH2F* hLs2_x[16];
 TH2F* hLs2_th[16];
 TH2F* hLs2_ph[16];
 TH2F* hLs2_path[16];
 TH2F* hLs2_tof[16];
 TH1F* hLtof[16];
 TH2F* hLtof_adc[16];
 TH1F* hRtof[16];
 TH2F* hcoin_event[16];
 TH2F* hLs2_f1_event[16];
 TH2F* hLs2r_y[16];
 int bin_Rs2adc,min_Rs2adc,max_Rs2adc,bin_Ls2adc,min_Ls2adc,max_Ls2adc; 
 min_Rs2adc=5000;
 max_Rs2adc=6000;
 min_Ls2adc=5000;
 max_Ls2adc=6000;
 bin_Rs2adc=(max_Rs2adc-min_Rs2adc)/5;
 bin_Ls2adc=(max_Ls2adc-min_Ls2adc)/5;
 for(int i=0;i<16;i++){
   hcoin_Rs2r[i]=new TH2F(Form("hcoin_Rs2r[%d]",i),Form("Coin vs R-HRS S2 CH%d R ADC",i),bin_Rs2adc,min_Rs2adc,max_Rs2adc,bin_coin_c,min_coin_c,max_coin_c);
   hcoin_Rs2l[i]=new TH2F(Form("hcoin_Rs2l[%d]",i),Form("Coin vs R-HRS S2 CH%d L ADC",i),bin_Rs2adc,min_Rs2adc,max_Rs2adc,bin_coin_c,min_coin_c,max_coin_c);
   hcoin_Ls2r[i]=new TH2F(Form("hcoin_Ls2r[%d]",i),Form("Coin vs L-HRS S2 CH%d R ADC",i),bin_Ls2adc,min_Ls2adc,max_Ls2adc,bin_coin_c,min_coin_c,max_coin_c);
   hcoin_Ls2l[i]=new TH2F(Form("hcoin_Ls2l[%d]",i),Form("Coin vs L-HRS S2 CH%d L ADC",i),bin_Ls2adc,min_Ls2adc,max_Ls2adc,bin_coin_c,min_coin_c,max_coin_c);
   //   hcoin_Rs2[i]=new TH2F(Form("hcoin_Rs2_%d",i),Form("R-HRS S2 seg%d  Mean time ",i),200,0,20);

   hLs2[i]=new TH2F(Form("hLs2[%d]",i),Form("L-HRS S2 CH%d ADC R&L PMT Hist",i),bin_Ls2adc,min_Ls2adc,max_Ls2adc,bin_Ls2adc,min_Ls2adc,max_Ls2adc);
   hLs2_x[i]=new TH2F(Form("hLs2_x[%d]",i),Form("L-HRS S2[%d] Coincidence vs FPx ",i),bin_coin_c,min_coin_c,max_coin_c,500,-1.,1.);
   set->SetTH2(hLs2_x[i],"","Coincidence [ns]","FP x [m]");
   hLs2_y[i]=new TH2F(Form("hLs2_y[%d]",i),Form("L-HRS S2[%d] Coincidence vs FPy ",i),bin_coin_c,min_coin_c,max_coin_c,500,-0.2,0.2);
   set->SetTH2(hLs2_y[i],"","Coincidence [ns]","FP y [m]");
    hLs2_th[i]=new TH2F(Form("hLs2_th[%d]",i),Form("L-HRS S2[%d] Coincidence vs FPth ",i),bin_coin_c,min_coin_c,max_coin_c,500,-0.2,0.2);
   set->SetTH2(hLs2_th[i],"","Coincidence [ns]","FP th [m]");
    hLs2_ph[i]=new TH2F(Form("hLs2_ph[%d]",i),Form("L-HRS S2[%d] Coincidence vs FPph ",i),bin_coin_c,min_coin_c,max_coin_c,500,-0.1,0.1);
   set->SetTH2(hLs2_ph[i],"","Coincidence [ns]","FP ph [m]");
   hLs2_path[i]=new TH2F(Form("hLs2_path[%d]",i),"",bin_coin_c,min_coin_c,max_coin_c,500,25,30);
   set->SetTH2(hLs2_path[i],"","Coincidence [ns]","Path Length [m]");
   hLs2_tof[i]=new TH2F(Form("hLs2_tof[%d]",i),"",bin_coin_c,min_coin_c,max_coin_c,1000,-100,100);
   hLtof[i]=new TH1F(Form("hLtof[%d]",i),"",200,-20,-5.);
   set->SetTH1(hLtof[i],Form("S2[%d]-S0 ",i),"TOF [ns]","Counts");
   hRtof[i]=new TH1F(Form("hRtof[%d]",i),"",200,-20,-5.);
   set->SetTH1(hRtof[i],Form("S2[%d]-S0 ",i),"TOF [ns]","Counts");
   hLtof_adc[i]=new TH2F(Form("hLtof_adc[%d]",i),"",500,-20,-5.,bin_Ls2adc,min_Ls2adc,max_Ls2adc);
   set->SetTH2(hLtof_adc[i],Form("S2[%d]-S0 vs S2r[%d] ADC",i,i),"TOF [ns]","ADC [ch]");


   
}



 for(int i=0;i<nth;i++){
 hcoin_ac1[i]=new TH2F(Form("hcoin_ac1[%d]",i),"Coinc time vs AC1 ADC Hist w/ correction",min_coin_c,max_coin_c,bin_ac1,min_ac1,max_ac1,bin_coin_c);
 hcoin_ac2[i]=new TH2F(Form("hcoin_ac2[%d]",i),"Coinc time vs AC2 ADC Hist w/ correction",bin_coin_c,min_coin_c,max_coin_c,bin_ac2,min_ac2,max_ac2);
 hcoin_t1[i]=new TH1F(Form("hcoin_t1[%d]",i), Form("Coincidence 0<AC1<%lf cut",ac1_adc[i]),bin_coin_c,min_coin_c,max_coin_c);
 hcoin_t2[i]=new TH1F(Form("hcoin_t2[%d]",i), Form("Coincidence %lf<AC2<%lf  cut",ac1_adc[i],th2_max),bin_coin_c,min_coin_c,max_coin_c);

 hcoin_ac1_max[i]=new TH1F(Form("hcoin_ac1_max[%d]",i), Form("Coincidence time %lf<AC2<%lf cut",ac2_adc[i],th2_max),bin_coin_c,min_coin_c,max_coin_c);
 hcoin_ac2_max[i]=new TH1F(Form("hcoin_ac2_max[%d]",i), Form("Coincidence time AC1<%lf  cut",ac1_adc[i]),bin_coin_c,min_coin_c,max_coin_c);

 for(int j=0;j<nth;j++)hcoin_t3[i][j]=new TH1F(Form("hcoin_t3[%d][%d]",i,j),Form("Coincidence time S2R-S2L[ns] ac1_adc< %lf, ac2_adc<%lf; Coin-Time [ns];Counts ",ac1_adc[i],ac2_adc[j]),bin_coin_c,min_coin_c,max_coin_c);


 }


 int bin_seg=16;
 double min_seg=0.0; double max_seg=15.;
 TH2F*hcoin_rseg=new TH2F("hcoin_rseg","Coinc Time vs R-S2seg Hist",bin_seg,min_seg,max_seg,bin_coin_c,min_coin_c,max_coin_c);
 TH2F*hcoin_lseg=new TH2F("hcoin_lseg","Coinc Time vs L-S2seg Hist",bin_seg,min_seg,max_seg,bin_coin_c,min_coin_c,max_coin_c);
 
 TH1F* hcoin_rs2[16]; 
 TH1F* hcoin_ls2[16];
 TH1F* hcoin_rs2_def[16]; 
 TH1F* hcoin_ls2_def[16];
 TGraph* gcoin_event=new TGraph();
 TGraph* gLS2_event[16]; 
 set->SetGr(gcoin_event,"Coin time event correlation plots ","events","coin time [ns]" );
 

 
 TH2F* hLs2r_event[16];
 TH1F* hLs2r_tdc[16];
 TH1F* hLs2l_tdc[16];
 
 TH2F* hLs2l_event[16];
 TH2F* hLs2r_f1_event[16];
 TH2F* hLs2l_f1_event[16];
 TH2F* hLs2r_f1_tw[16];
 TH2F* hLs2l_f1_tw[16];
 TH2F* hLs2_event[16];
 TH1F* hLs2r_f1_tdc[16];
 TH1F* hLs2l_f1_tdc[16];
 TH1F* hLs2r_fbus[16];
 TH1F* hLs2l_fbus[16];
 TH2F* hf1_fbus_r[16];
 TH2F* hf1_fbus[16];
 TH1F* hLs2r_fadc[16];
 TH1F* hLs2l_fadc[16];
 TH2F* hf1_fadc_r[16];
 TH2F* hf1_fadc_l[16];
 TH2F* hfadc_r_event[16];
 TH2F* hfbus_r_event[16];

 //----- FBUS -----------------//
 double fbus_time=0.5;
 double min_fbus=2700*fbus_time;
 double max_fbus=2900*fbus_time;
 int bin_fbus=(int)(max_fbus - min_fbus)/fbus_time;


 //------ FADC ---------------//

 double fadc_time=1.0;//branch was already changed to ns 
   //0.0625; // [ns]
 double min_fadc=140.0;
 double max_fadc=150.0;
 int bin_fadc=(int)(max_fadc-min_fadc)/fadc_time;
 
 


 
 
 double min_ls2t=640;
 double max_ls2t=670;
 double min_rs2t=590;
 double max_rs2t=620;
 int bin_s2t=60;
   //(int)(max_ls2t-min_ls2t)/fbus_time;
 //  bin_s2t=180.;
 double min_Ls2a_c=10.;
 double max_Ls2a_c=3000.;
 int bin_Ls2a_c=1000;

 //----- F1TDC --------------//
 double min_f1_s2t=0.0*tdc_time;
 double max_f1_s2t=54000*tdc_time;
 min_f1_s2t=12000*tdc_time;
 max_f1_s2t=13500*tdc_time;
 int bin_f1_s2t=200;//(int)max_f1_s2t-min_f1_s2t;
 double min_s2,max_s2;
 int bin_s2;
 min_s2=50.0;
 max_s2=70.;
 bin_s2=40;


 
 int bin_ev=10000;
 TH2F* hy_event=new TH2F("hy_event","",5500,0,evnt,100,-0.2,0.2);
 set->SetTH2(hy_event,"Y events shift","Events","y [m]");
 for(int i=0;i<16;i++){


   //=========== FBUS =============//

 hLs2l_fbus[i]=new TH1F(Form("hLs2l_fbus_%d",i),"",bin_fbus,min_fbus,max_fbus);
 set->SetTH1(hLs2l_fbus[i],Form("LHRS S2 %d Seg L-PMT FBUS TDC infomaiton",i),"[ns]","counts");
 
 hf1_fbus[i]=new TH2F(Form("hf1_fbus_%d",i),"",bin_fbus,min_fbus,max_fbus,bin_f1_s2t,min_f1_s2t,max_f1_s2t);
 set->SetTH2(hf1_fbus[i],Form("LHRS S2 %d Seg  FBUS vs F1TDC L-PMT TDC hist ",i),"FBUS [ns]","F1 [ns]");

  hLs2r_fbus[i]=new TH1F(Form("hLs2r_fbus_%d",i),"",bin_fbus,min_fbus,max_fbus);
 set->SetTH1(hLs2r_fbus[i],Form("LHRS S2 %d Seg R-PMT FBUS TDC infomaiton",i),"[ns]","counts");
 
 hf1_fbus_r[i]=new TH2F(Form("hf1_fbus_r_%d",i),"",bin_fbus,min_fbus,max_fbus,bin_f1_s2t,min_f1_s2t,max_f1_s2t);
 set->SetTH2(hf1_fbus_r[i],Form("LHRS S2 %d Seg R-PMT FBUS vs F1TDC TDC hist ",i),"FBUS [ns]","F1 [ns]");

 //============== FADC =============//

 hLs2l_fadc[i]=new TH1F(Form("hLs2l_fadc_%d",i),"",bin_fadc,min_fadc,max_fadc);
 set->SetTH1(hLs2l_fadc[i],Form("LHRS S2 %d Seg L-PMT FADC TDC infomaiton",i),"[ns]","counts");
 
 hf1_fadc_l[i]=new TH2F(Form("hf1_fadc_l_%d",i),"",bin_fadc,min_fadc,max_fadc,bin_f1_s2t,min_f1_s2t,max_f1_s2t);
 set->SetTH2(hf1_fadc_l[i],Form("LHRS S2 %d Seg  FADC vs F1TDC L-PMT TDC hist ",i),"FADC [ns]","F1 [ns]");

  hLs2r_fadc[i]=new TH1F(Form("hLs2r_fadc_%d",i),"",bin_fadc,min_fadc,max_fadc);
 set->SetTH1(hLs2r_fadc[i],Form("LHRS S2 %d Seg R-PMT FADC TDC infomaiton",i),"[ns]","counts");
 
 hf1_fadc_r[i]=new TH2F(Form("hf1_fadc_r_%d",i),"",bin_fadc,min_fadc,max_fadc,bin_f1_s2t,min_f1_s2t,max_f1_s2t);
 set->SetTH2(hf1_fadc_r[i],Form("LHRS S2 %d Seg R-PMT FADC vs F1TDC TDC hist ",i),"FADC [ns]","F1 [ns]");



 
 ///////////////////////////////////////////////////////////
 
   hcoin_rs2[i]=new TH1F(Form("hcoin_rs2[%d]",i),Form("Coincidence Time with R-S2 seg%d",i),bin_coin_c,min_coin_c,max_coin_c);
   hcoin_ls2[i]=new TH1F(Form("hcoin_ls2[%d]",i),Form("Coincidence Time with L-S2 seg%d",i),bin_coin_c,min_coin_c,max_coin_c);
   hcoin_rs2_def[i]=new TH1F(Form("hcoin_rs2_def[%d]",i),Form("Coincidence Time with R-S2 seg%d",i),bin_coin_c,min_coin_c,max_coin_c);
   hcoin_ls2_def[i]=new TH1F(Form("hcoin_ls2_def[%d]",i),Form("Coincidence Time with L-S2 seg%d",i),bin_coin_c,min_coin_c,max_coin_c);
   gLS2_event[i]=new TGraph();
   set->SetGr(gLS2_event[i],Form("HRS-L S2[%d] event correlation plots ",i),"events","coin time [ns]");


   hcoin_event[i]=new TH2F(Form("hcoin_event[%d]",i),"",bin_ev,0,evnt,bin_coin_c,min_coin_c,max_coin_c);
   set->SetTH2(hcoin_event[i],"Coincidence event correlation","Event number","coin time [ns]");
   

   hLs2r_tdc[i]=new TH1F(Form("hLs2r_tdc[%d]",i),"",bin_s2t,min_rs2t,max_rs2t);
   hLs2l_tdc[i]=new TH1F(Form("hLs2l_tdc[%d]",i),"",bin_s2t,min_ls2t,max_ls2t);
   set->SetTH1(hLs2r_tdc[i],Form("Ls2[%d] R-PMT TDC",i),"TDC [ns]","Counts");
   set->SetTH1(hLs2l_tdc[i],Form("Ls2[%d] L-PMT TDC",i),"TDC [ns]","Counts");
   
   hLs2r_event[i]=new TH2F(Form("hLs2r_event[%d]",i),"",bin_ev,0,evnt,bin_s2t,min_rs2t,max_rs2t);
   set->SetTH2(hLs2r_event[i],Form("Ls2[%d]-R PMT event correlation",i),"events","R-PMT TDC [sec]");
   hLs2l_event[i]=new TH2F(Form("hLs2l_event[%d]",i),"",bin_ev,0,evnt,bin_s2t,min_ls2t,max_ls2t);
   set->SetTH2(hLs2l_event[i],Form("LS2[%d]-L PMT event correlation",i),"events","L-PMT TDC [sec]");


   
   hLs2r_y[i]=new TH2F(Form("hLs2r_y[%d]",i),"",500,-0.2,0.2,bin_f1_s2t,min_f1_s2t,max_f1_s2t);
   set->SetTH2(hLs2r_y[i],Form("L-HRS S2 R-PMT [%d] TDC vs FPy ",i),"FP y [m]","TDC [ns]");   
   hLs2r_f1_event[i]=new TH2F(Form("hLs2r_f1_event[%d]",i),"",bin_ev,0,evnt,bin_f1_s2t,min_f1_s2t,max_f1_s2t);
   set->SetTH2(hLs2r_f1_event[i],Form("F1 Ls2[%d]-R PMT event correlation",i),"events","R-PMT TDC [ns]");
   hLs2l_f1_event[i]=new TH2F(Form("hLs2l_f1_event[%d]",i),"",bin_ev,0,evnt,bin_f1_s2t,min_f1_s2t,max_f1_s2t);
   set->SetTH2(hLs2l_f1_event[i],Form("F1 LS2[%d]-L PMT event correlation",i),"events","L-PMT TDC [ns]");
   hLs2_event[i]=new TH2F(Form("hLs2_event[%d]",i),"",bin_ev,0,evnt,bin_s2,min_s2,max_s2);
     set->SetTH2(hLs2_event[i],Form("LS2[%d]( L-R TDC) event correlation",i),"events","L-R TDC [ns]");
   //   hLs2_f1_event[i]=new TH2F(Form("hLs2_f1_event[%d]",i),"",bin_ev,0,evnt,30000,50000,80000);
      hLs2_f1_event[i]=new TH2F(Form("hLs2_f1_event[%d]",i),"",bin_ev,0,evnt,300,50000,80000);
   set->SetTH2(hLs2_f1_event[i],Form("LS2[%d] F1(  L-R TDC) event correlation",i),"events","L-R TDC [ns]");

   hLs2r_f1_tdc[i]=new TH1F(Form("hLs2r_f1_tdc[%d]",i),"",bin_f1_s2t,min_f1_s2t,max_f1_s2t);
   set->SetTH1(hLs2r_f1_tdc[i],Form("L-HRS seg %d R-PMT TDC ",i),"TDC [ns]","Counts");
   
      hLs2l_f1_tdc[i]=new TH1F(Form("hLs2l_f1_tdc[%d]",i),"",bin_f1_s2t,min_f1_s2t,max_f1_s2t);
   set->SetTH1(hLs2l_f1_tdc[i],Form("L-HRS seg %d L-PMT TDC ",i),"TDC [ns]","Counts");

   hLs2r_f1_tw[i]=new TH2F(Form("hLs2r_f1_tw[%d]",i),"",bin_Ls2a_c,min_Ls2a_c,max_Ls2a_c,bin_f1_s2t,min_f1_s2t,max_f1_s2t);
   set->SetTH2(hLs2r_f1_tw[i],"L-HRS S2 R-PMT Time walk","ADC [ch]","TDC [ns]");
   hLs2l_f1_tw[i]=new TH2F(Form("hLs2l_f1_tw[%d]",i),"",bin_Ls2a_c,min_Ls2a_c,max_Ls2a_c,bin_f1_s2t,min_f1_s2t,max_f1_s2t);
   set->SetTH2(hLs2l_f1_tw[i],"L-HRS S2 L-PMT Time walk","ADC [ch]","TDC [ns]");

   hLs2r_trig[i]=new TH2F(Form("hLs2r_trig[%d]",i),"",10,0,10,bin_f1_s2t,min_f1_s2t,max_f1_s2t);   
      set->SetTH2(hLs2r_trig[i],"L-HRS S2 R-PMT Triger ","Trigger number","TDC [ns]");
 }



  TH2F* hLs2r_f1_seg=new TH2F("hLs2r_f1_seg","",bin_ev,0,evnt,bin_seg,min_seg,max_seg);
 set->SetTH1(hLs2r_f1_seg,"HRS-L R-PMT Segment dependence of event","event","Segment");
  TH2F* hLs2l_f1_seg=new TH2F("hLs2l_f1_seg","",bin_ev,0,evnt,bin_seg,min_seg,max_seg);
 set->SetTH1(hLs2l_f1_seg,"HRS-L L-PMT Segment dependence of event","event","Segment");


 cout<<"define hist "<<endl;
 
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
 bool cut_ac1,cut_ac2,cut_beta;
 int nac1,nac2,nac3,n;
 double tofs_r,tof_l; 
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
 bool cut_Rs2,cut_Ls2,cut_rpathl,cut_lpathl,cut_coin,cut_rbeta,cut_lbeta,cut_vz,cut_Rx,cut_trig,coin_trig,right_trig,cut_track,cut_s0,tri_T5;
 double Rtof;
 double Ltof;
 double Rtof_cut;
 double tof_r;
 int ev=0;
 int e=0;
 int ev_k=0;




 
 cout<<"start Fill "<<endl;
 
 for(int k=0;k<evnt;k++){

   
   if(k==ev*100000){
     cout<<"Fill Event: "<<k<<"/"<<evnt<<endl; 
     ev=ev+1;}


   /*
   //=== Initialize
   for(int l=0;l<16;l++){

     f1_Rs2r[l]=-2222.0;
     f1_Rs2l[l]=-2222.0;          
     f1_Ls2r[l]=-2222.0;
     f1_Ls2l[l]=-2222.0;

     fbus_Rs2r[l]=-2222.0;
     fbus_Rs2l[l]=-2222.0;          
     fbus_Ls2r[l]=-2222.0;
     fbus_Ls2l[l]=-2222.0;          

     fadc_Rs2r[l]=-2222.0;
     fadc_Rs2l[l]=-2222.0;          
     fadc_Ls2r[l]=-2222.0;
     fadc_Ls2l[l]=-2222.0;          
   }

   */
   coin_f1=-2222.0;
   coin_fbus=-2222.0;
   coin_fadc=-2222.0;
   coin_f1_c=-2222.0;
   coin_fbus_c=-2222.0;
   coin_fadc_c=-2222.0;   

   


   //======= Get Entry ====//
   
   T->GetEntry(k); 
   event=k;
   trig=(int)DRevtype;



   
 pe_=Lp[0];//*sqrt(1+pow(Lth[0],2)+pow(Lph[0],2));
 pk=Rp[0];//*sqrt(1+pow(Rth[0],2)+pow(Rph[0],2));
 ppi=Rp[0];//*sqrt(1+pow(Rth[0],2)+pow(Rph[0],2));
 pe=4.313; // GeV   //hallap*1.0e-3;
 Ee=sqrt(pow(pe,2)+pow(me,2));
 Ee_=sqrt(pow(pe_,2)+pow(me,2));
 Epi=sqrt(pow(ppi,2)+pow(mpi,2));
 Ek=sqrt(pow(pk,2)+pow(mk,2));
 Ls2pads=(int)Ls2tpads[0];
 Rs2pads=(int)Rs2tpads[0];


 
 rpathl=rtrpathl[0]+rs2pathl[0]; // R-HRS path length S2 -RF
 lpathl=ltrpathl[0]+ls2pathl[0]; // L-HRS path length S2 -RF

 if(Pion_flag==1){rbeta=ppi/Epi;
   if(k==1){cout<<"Pion mass beta"<<endl;}} //Pion
 else{rbeta=pk/Ek; 
   if(k==1){cout<<"Kaon mass beta"<<endl;}}//Kaon}
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
  //  Rtof=tof_r-(-LF1[Ls2pads+48]+LF1[37])*tdc_time-rpath_corr+lpath_corr;
  //  Rtof=tof_r;//-rpath_corr+lpath_corr;
  //  Rtof_cut=(-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9])/2.0*tdc_time;  
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
   tri_T5=false;
   if(Rvz_cutmin<Rvz[0] && Rvz[0]<Rvz_cutmax && Lvz_cutmin<Lvz[0] && Lvz[0]<Lvz_cutmax)cut_vz=true;
  
   if(mode=="G"){
     if(ctime[0]>-100)coin_trig=true;
     cut_track=true;
     cut_rpathl=true;
     cut_lpathl=true;
     tri_T5=true;
   }else{
 
  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads)cut_track=true;
  if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax)cut_rpathl=true;
  if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax)cut_lpathl=true;
  if(RF1[43]>0 && RF1[44]>0 && LF1[27]>0 && LF1[28]>0)cut_s0=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;
  if(cut_Rs2 && cut_Ls2 && cut_s0)coin_trig=true; 
  if(DRevtype==5)tri_T5=true; 
}
 //=======================================//
  


   //==========================================//
   //========= Fill Hist =====================//
   //========================================//


 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track){
 hcoin_t->Fill(coin_t);
 hcoin_tc->Fill(coin_tc);
 ha1_a2->Fill(Ra1sum,Ra2sum);// AC1 vs AC2

 //--------- with AC1 Cut ---------------// 
  for(int j=0;j<nth;j++){
   cut_ac1=false;
   if(Ra1sum<ac1_adc[j])cut_ac1=true;
   if(cut_ac1 && Ra2sum<th2_max){
    hcoin_t2[j]->Fill(coin_tc); //AC1 cut   
    hcoin_ac2[j]->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
    hcoin_ac2_max[j]->Fill(coin_tc,Ra2sum); //AC1 cut && (AC2 variable cut)
     }
    }
  //-------with AC2 Cut --------------------//
   for(int k=0;k<nth;k++){
    cut_ac2=false;
    if(Ra2sum>ac2_adc[k])cut_ac2=true;
    if(cut_ac2 && Ra1sum<th1_max){
    hcoin_t1[k]->Fill(coin_tc); // AC2 cut
    hcoin_ac1[k]->Fill(coin_tc,Ra1sum);//AC2 cut && (AC1 variable cut)  B
    hcoin_ac1_max[k]->Fill(coin_tc,Ra1sum);//AC2 cut && (AC1 variable cut)  
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


 //  Pion Cut //
 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track && Ra1sum>ac1_kcut 
    && Ra2sum>ac2_kcut_max && tri_T5){
  
    hcoin_rs2[Rs2pads]->Fill(coin_tc); 
    hcoin_ls2[Ls2pads]->Fill(coin_tc);
    hcoin_Rs2r[Rs2pads]->Fill(Rs2ra[Rs2pads],coin_tc);
    hcoin_Rs2l[Rs2pads]->Fill(Rs2la[Rs2pads],coin_tc);
    hcoin_Ls2r[Ls2pads]->Fill(Ls2ra[Ls2pads],coin_tc);
    hcoin_Ls2l[Ls2pads]->Fill(Ls2la[Ls2pads],coin_tc);
    hLs2[Ls2pads]->Fill(Ls2la[Ls2pads],Ls2ra[Ls2pads]);
    hLs2_x[Ls2pads]->Fill(coin_tc,Lx[0]);
    hLs2_y[Ls2pads]->Fill(coin_tc,Ly[0]);
    hLs2_th[Ls2pads]->Fill(coin_tc,Lth[0]);
    hLs2_ph[Ls2pads]->Fill(coin_tc,Lph[0]);    
    hLs2_path[Ls2pads]->Fill(coin_tc,lpathl);
    hLtof[Ls2pads]->Fill(Ltof);
    hLtof_adc[Ls2pads]->Fill(Ltof,Ls2ra[Ls2pads]);

 }

         //------ Coincidence time event correlation -----//
 //    if(k==e){
 //   if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track && Ra1sum>ac1_kcut   && Ra2sum>ac2_kcut_max && tri_T5){
     //    if(cut_lpathl){


       f1_s2r=(-LF1[Ls2pads]+LF1[30])*tdc_time; //in real this is left PMT
       f1_s2l=(-LF1[Ls2pads+48]+LF1[37])*tdc_time;


       //========================//
       //======== R-HRS =========//
       //========================//

       //--------- F1TDC --------//
       f1_Rs2r[Rs2pads]=(-RF1[Rs2pads+48]+RF1[46])*tdc_time;
       f1_Rs2l[Rs2pads]=(-RF1[16+Rs2pads]+RF1[9])*tdc_time;

       
       //--------- FBUS --------//
       fbus_Rs2r[Rs2pads]=Rs2rt[Rs2pads]*fbus_time;
       fbus_Rs2l[Rs2pads]=Rs2lt[Rs2pads]*fbus_time;
   
       //--------- FADC --------//
       fadc_Rs2r[Rs2pads]=Rs2rtc_fadc[Rs2pads];
       fadc_Rs2l[Rs2pads]=Rs2ltc_fadc[Rs2pads];

       

       
       //========================//
       //======== L-HRS =========//
       //========================//
       
       //--------- F1TDC --------//
       f1_Ls2r[Ls2pads]=(-LF1[Ls2pads+48]+LF1[37])*tdc_time;
       f1_Ls2l[Ls2pads]=(-LF1[Ls2pads]+LF1[30])*tdc_time;
       
       //--------- FBUS --------//
       fbus_Ls2r[Ls2pads]=Ls2rt[Ls2pads]*fbus_time;
       fbus_Ls2l[Ls2pads]=Ls2lt[Ls2pads]*fbus_time;

       //--------- FADC --------//
       fadc_Ls2r[Ls2pads]=Ls2rtc_fadc[Ls2pads];
       fadc_Ls2l[Ls2pads]=Ls2ltc_fadc[Ls2pads];
       

       //------- Coin time ---------//

       coin_f1=(f1_Rs2r[Rs2pads]+f1_Rs2l[Rs2pads])/2.0
	 -(f1_Ls2r[Ls2pads]+f1_Ls2l[Ls2pads])/2.0;
       coin_f1_c=coin_tc; // offset & path correction
       
       coin_fbus=(fbus_Rs2r[Rs2pads]+fbus_Rs2l[Rs2pads])/2.0
	 -(fbus_Ls2r[Ls2pads]+fbus_Ls2l[Ls2pads])/2.0;
       coin_fbus_c=coin_fbus+rpath_corr-lpath_corr;
       
       coin_fadc=(fadc_Rs2r[Rs2pads]+fadc_Rs2l[Rs2pads])/2.0
		    -(fadc_Ls2r[Ls2pads]+fadc_Ls2l[Ls2pads])/2.0;

       coin_fadc_c=coin_fadc+rpath_corr-lpath_corr;
       
  if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track && tri_T5){
      hcoin_R->Fill(Rtof);
      hcoin_2->Fill(coin_tc);
      hcoin_event[Ls2pads]->Fill(k,coin_tc);
      gcoin_event->SetPoint(k,k,coin_tc); 
      hy_event->Fill(k,Ly[0]);     
      gLS2_event[Ls2pads]->SetPoint(k,k,coin_tc);
      hLs2r_tdc[Ls2pads]->Fill(Ls2r_tc[Ls2pads]*1.0e9);
      hLs2l_tdc[Ls2pads]->Fill(Ls2l_tc[Ls2pads]*1.0e9);
      hLs2r_event[Ls2pads]->Fill(k,Ls2r_tc[Ls2pads]*1.0e9);
      hLs2l_event[Ls2pads]->Fill(k,Ls2l_tc[Ls2pads]*1.0e9);
      hLs2_event[Ls2pads]->Fill(k,(Ls2l_tc[Ls2pads]-Ls2r_tc[Ls2pads])*1.0e9);
      hLs2r_y[Ls2pads]->Fill(Ly[0],f1_s2r);      
      hLs2r_f1_event[Ls2pads]->Fill(k,f1_s2r);
      hLs2l_f1_event[Ls2pads]->Fill(k,f1_s2l);
      hLs2_f1_event[Ls2pads]->Fill(k,f1_s2r+f1_s2l);
      hLs2r_f1_seg->Fill(k,Ls2pads,f1_s2r);
      hLs2l_f1_seg->Fill(k,Ls2pads,f1_s2l);
      hLs2r_f1_tdc[Ls2pads]->Fill(f1_s2r);
      hLs2l_f1_tdc[Ls2pads]->Fill(f1_s2l);
      hLs2r_f1_tw[Ls2pads]->Fill(Ls2ra_c[Ls2pads],f1_s2r);
      hLs2l_f1_tw[Ls2pads]->Fill(Ls2la_c[Ls2pads],f1_s2l);

      //----------- FBUS -------------//
      hLs2l_fbus[Ls2pads]->Fill(Ls2lt[Ls2pads]*fbus_time);
      hf1_fbus[Ls2pads]->Fill(Ls2lt[Ls2pads]*fbus_time,f1_s2r);
      hLs2r_fbus[Ls2pads]->Fill(Ls2rt[Ls2pads]*fbus_time);
      hf1_fbus_r[Ls2pads]->Fill(Ls2rt[Ls2pads]*fbus_time,f1_s2l);
 
      //---------- FADC -----------//

      hLs2l_fadc[Ls2pads]->Fill(Ls2ltc_fadc[Ls2pads]);
      hf1_fadc_l[Ls2pads]->Fill(Ls2ltc_fadc[Ls2pads],f1_s2r);
      hLs2r_fadc[Ls2pads]->Fill(Ls2rtc_fadc[Ls2pads]);
      hf1_fadc_r[Ls2pads]->Fill(Ls2rtc_fadc[Ls2pads],f1_s2l);

      
    }
  hLs2r_trig[Ls2pads]->Fill(DRevtype,f1_s2r);  

    
 //-------- Proton Cut Hist --------------//
 if(coin_trig && cut_vz && cut_lpathl && cut_rpathl && cut_track && Ra1sum<ac1_kcut 
    && Ra2sum<ac2_kcut_min)hcoin_p->Fill(coin_tc);

 tnew->Fill();

 }
 
 cout<<"Filled Hist "<<endl;



 //============ COINCIDENCE SEGMENT CORRECTION ================//

 double min_pi,max_pi,rn_pi[16],rsig_pi[16],rmean_pi[16],ln_pi[16],lsig_pi[16],lmean_pi[16],rmean_pi_err[16],lmean_pi_err[16];
 double diff_r[16],diff_l[16],r_pmax[16],l_pmax[16],nl_pmax[16],nr_pmax[16],bin_rmax[16],bin_lmax[16],bin_rxmax[16],bin_rymax[16],bin_lxmax[16],bin_lymax[16];
 double rn_pi_def[16],rsig_pi_def[16],rmean_pi_def[16],ln_pi_def[16],lsig_pi_def[16],lmean_pi_def[16],rmean_pi_err_def[16],lmean_pi_err_def[16];
 double diff_r_def[16],diff_l_def[16],r_pmax_def[16],l_pmax_def[16],nl_pmax_def[16],nr_pmax_def[16],bin_rmax_def[16],bin_lmax_def[16],bin_rxmax_def[16],bin_rymax_def[16],bin_lxmax_def[16],bin_lymax_def[16];
 min_pi=0., max_pi=25;
 TF1* ffit_r[16];
 TF1* ffit_l[16];


 TGraphErrors* gs2r_mean=new TGraphErrors();
 gs2r_mean->SetName("gs2r_mean");
 gs2r_mean->SetTitle("coincidence tiem of r-s2 segmments; segment num ;coin time [ns]");
 gs2r_mean->GetXaxis()->SetTitleSize(0.05);
 gs2r_mean->GetYaxis()->SetTitleSize(0.05);
 gs2r_mean->GetXaxis()->SetLabelSize(0.05);
 gs2r_mean->GetYaxis()->SetLabelSize(0.04);
 TGraphErrors* gs2l_mean=new TGraphErrors();
 gs2l_mean->SetName("gs2l_mean");
 gs2l_mean->SetTitle("coincidence tiem of l-s2 segments dependence ; segment num ;coin time [ns] ");
 gs2l_mean->GetXaxis()->SetTitleSize(0.05);
 gs2l_mean->GetYaxis()->SetTitleSize(0.05);
 gs2l_mean->GetXaxis()->SetLabelSize(0.05);
 gs2l_mean->GetYaxis()->SetLabelSize(0.04);

 gs2r_mean->GetYaxis()->SetRangeUser(-2.0,2.0);
 gs2l_mean->GetYaxis()->SetRangeUser(-2.0,2.0);

 TGraphErrors* gs2r_nevent=new TGraphErrors();
 set->SetGr(gs2r_nevent,"","","");
 gs2r_nevent->SetName("gs2r_nevent");
 gs2r_nevent->SetTitle("coincidence tiem of r-s2 segmments; segment num ;# of Events");
 gs2r_nevent->GetXaxis()->SetTitleSize(0.05);
 gs2r_nevent->GetYaxis()->SetTitleSize(0.05);
 gs2r_nevent->GetXaxis()->SetLabelSize(0.05);
 gs2r_nevent->GetYaxis()->SetLabelSize(0.04);

 TGraphErrors* gs2r_nevent2=new TGraphErrors();
 set->SetGr(gs2r_nevent2,"","","");
 gs2r_nevent2->SetName("gs2r_nevent2");
 gs2r_nevent2->SetTitle("coincidence tiem of r-s2 segmments; segment num ;# of Events");
 gs2r_nevent2->GetXaxis()->SetTitleSize(0.05);
 gs2r_nevent2->GetYaxis()->SetTitleSize(0.05);
 gs2r_nevent2->GetXaxis()->SetLabelSize(0.05);
 gs2r_nevent2->GetYaxis()->SetLabelSize(0.04);

 //=======================================================//



ffit_r[7] =new TF1(Form("ffit_r[%d]",7),"gaus",min_coin_c,max_coin_c);
ffit_l[7] =new TF1(Form("ffit_l[%d]",7),"gaus",min_coin_c,max_coin_c);

 hcoin_rs2[7]->Fit(Form("ffit_r[%d]",7),"0","",min_pi,max_pi);
 rmean_pi[7]=ffit_r[7]->GetParameter(1);
 rmean_pi_err[7]=ffit_r[7]->GetParError(1);
 hcoin_ls2[7]->Fit(Form("ffit_l[%d]",7),"0","",min_pi,max_pi);
 lmean_pi[7]=ffit_l[7]->GetParameter(1);
 lmean_pi_err[7]=ffit_l[7]->GetParError(1);


 //=======================================================//


   double rs2_off[16],ls2_off[16];




   for(int i=0;i<16;i++){   
   // for(int j=0;j<2;j++){
   ffit_r[i] =new TF1(Form("ffit_r[%d]",i),"gaus",min_coin_c,max_coin_c);
   ffit_r[i]->SetNpx(2000);
   ffit_l[i] =new TF1(Form("ffit_l[%d]",i),"gaus",min_coin_c,max_coin_c);
   ffit_l[i]->SetNpx(2000);

   //----------- l-hrs -------------------------//

   bin_rmax[i]= hcoin_rs2[i]->GetMaximumBin();
   r_pmax[i]= hcoin_rs2[i]->GetXaxis()->GetBinCenter(bin_rmax[i]);
   nr_pmax[i]= hcoin_rs2[i]->GetBinContent(bin_rmax[i]);
   rsig_pi[i]=3.5;   

   ffit_r[i]->SetParameter(0,nr_pmax[i]);
   ffit_r[i]->SetParameter(1,r_pmax[i]);
   ffit_r[i]->FixParameter(2,rsig_pi[i]);
   ffit_r[i]->SetParLimits(0,0.8*nr_pmax[i],1.2*nr_pmax[i]);
   ffit_r[i]->SetParLimits(1,r_pmax[i]-1.0,r_pmax[i]+1.0);
   ffit_r[i]->SetParLimits(2,0.5,4.0);
    
   hcoin_rs2[i]->Fit(Form("ffit_r[%d]",i),"Rb","",min_pi,max_pi);
   rmean_pi[i]=ffit_r[i]->GetParameter(1);
   rmean_pi_err[i]=ffit_r[i]->GetParError(1);
   nr_pmax[i]=ffit_r[i]->GetParameter(0);
   rsig_pi[i]=ffit_r[i]->GetParameter(2);
   nl_pmax[i]=hcoin_ls2[i]->GetYaxis()->GetBinCenter(hcoin_ls2[i]->GetMaximumBin());
   l_pmax[i]= hcoin_ls2[i]->GetXaxis()->GetBinCenter(hcoin_ls2[i]->GetMaximumBin());

   //----------- l-hrs -------------------------//

   bin_lmax[i]= hcoin_ls2[i]->GetMaximumBin();
   l_pmax[i]= hcoin_ls2[i]->GetXaxis()->GetBinCenter(bin_lmax[i]);
   nl_pmax[i]= hcoin_ls2[i]->GetBinContent(bin_lmax[i]);
   lsig_pi[i]=3.5;   




   ffit_l[i]->SetParameter(0,nl_pmax[i]);
   ffit_l[i]->SetParameter(1,l_pmax[i]);
   ffit_l[i]->FixParameter(2,lsig_pi[i]);
   ffit_l[i]->SetParLimits(0,0.8*nl_pmax[i],1.2*nl_pmax[i]);
   ffit_l[i]->SetParLimits(1,l_pmax[i]-1.0,l_pmax[i]+1.0);
   ffit_l[i]->SetParLimits(2,0.5,4.0);


   hcoin_ls2[i]->Fit(Form("ffit_l[%d]",i),"Rb","",min_pi,max_pi);
   lmean_pi[i]=ffit_l[i]->GetParameter(1);
   lmean_pi_err[i]=ffit_l[i]->GetParError(1);
   nl_pmax[i]=ffit_l[i]->GetParameter(0);
   lsig_pi[i]=ffit_l[i]->GetParameter(2);
   nl_pmax[i]=hcoin_ls2[i]->GetYaxis()->GetBinCenter(hcoin_ls2[i]->GetMaximumBin());
   l_pmax[i]= hcoin_ls2[i]->GetXaxis()->GetBinCenter(hcoin_ls2[i]->GetMaximumBin());



 if(0<i){
 gs2r_mean->SetPoint(i,i,(rmean_pi[i]-rmean_pi[7]));
 gs2r_mean->SetPointError(i,0.0,rmean_pi_err[i]);

 }else{
 gs2r_mean->SetPoint(i,i,0.0);
 gs2r_mean->SetPointError(i,0.0,0.0);

}
// gs2r_mean->setpointerror(i,0,rmean_pi_err[i]);

 
 if(i<15){
 gs2l_mean->SetPoint(i,i,lmean_pi[i]-lmean_pi[7]);
 gs2l_mean->SetPointError(i,0.0,lmean_pi_err[i]);

 }else{
   gs2l_mean->SetPoint(i,i,0.0);
}

 // gs2l_mean->setpointerror(i,0,lmean_pi_err[i]);

 diff_r[i]=-(rmean_pi[i]-rmean_pi[7])/tdc_time;
 diff_l[i]=-(lmean_pi[i]-lmean_pi[7])/tdc_time;
 rs2_off[i]=-(rmean_pi[i]-rmean_pi[7])/tdc_time*2+s2f1_off(i,"R",mode,kine);
 ls2_off[i]=(lmean_pi[i]-lmean_pi[7])/tdc_time*2+s2f1_off(i,"L",mode,kine);


}


   //===========L-HRS time shift analysis =========//

   TF1* fLs2r_f1_tdc[16];
   TGraphErrors* gs2r_shift=new TGraphErrors();
   TGraphErrors* gs2r_sig=new TGraphErrors();

   set->SetGr(gs2r_shift,"L-HRS TDC Shift ratio","L-HRS S2 R-PMT Segments","ratio");
   set->SetGr(gs2r_sig,"L-HRS TDC Sigma ","Segments","sigma [ns]");
   double peak[16];
   double mean_s2r[16];
   double sig_s2r[16];
   double sig_err_s2r[16];
   double inte_s2r[16];
   double sum_s2r[16];
   double shift_rate[16];
   for(int i=0;i<16;i++){
   fLs2r_f1_tdc[i]=new TF1(Form("fLs2r_f1_tdc[%d]",i),"gausn(0)",min_f1_s2t,max_f1_s2t);


   peak[i]=hLs2r_f1_tdc[i]->GetXaxis()->GetBinCenter(hLs2r_f1_tdc[i]->GetMaximumBin());
     //GetXaxis()->GetBinCenter(bin_rmax[i]);
     //GetXaxis()->GetBinCenter(hLs2r_f1_tdc[i]->GetMaximum());
   //   bin_rmax[i]= hcoin_rs2[i]->GetMaximumBin();
   //   r_pmax[i]= hcoin_rs2[i]->GetXaxis()->GetBinCenter(bin_rmax[i]);
   //   nr_pmax[i]= hcoin_rs2[i]->GetBinContent(bin_rmax[i]);
   hLs2r_f1_tdc[i]->Fit(Form("fLs2r_f1_tdc[%d]",i),"","Rq",peak[i]-5.,peak[i]+10);
   mean_s2r[i]=fLs2r_f1_tdc[i]->GetParameter(1);
   sig_s2r[i]=fLs2r_f1_tdc[i]->GetParameter(2);
   sig_err_s2r[i]=fLs2r_f1_tdc[i]->GetParError(2);



   //inte_s2r[i]=hLs2r_f1_tdc[i]->Integral(hLs2r_f1_tdc[i]->GetXaxis()->FindBin(mean_s2r[i]-3*sig_s2r[i]),hLs2r_f1_tdc[i]->GetXaxis()->FindBin(mean_s2r[i]+3*sig_s2r[i]));
   // inte_s2r[i]=inte_s2r[i]*(max_f1_s2t-min_f1_s2t)/bin_f1_s2t;//3sig Events

   //   inte_s2r[i]=fLs2r_f1_tdc[i]->Integral(mean_s2r[i]-3*sig_s2r[i],mean_s2r[i]+3*sig_s2r[i]);
   //Fit func Integral

   inte_s2r[i]=hLs2r_f1_tdc[i]->Integral(hLs2r_f1_tdc[i]->GetXaxis()->FindBin(mean_s2r[i]-2.0),hLs2r_f1_tdc[i]->GetXaxis()->FindBin(mean_s2r[i]+2.0));
  inte_s2r[i]=inte_s2r[i]*(max_f1_s2t-min_f1_s2t)/bin_f1_s2t;//+-2 ns Events
   
   sum_s2r[i]=hLs2r_f1_tdc[i]->Integral(hLs2r_f1_tdc[i]->GetXaxis()->FindBin(min_f1_s2t)
					,hLs2r_f1_tdc[i]->GetXaxis()->FindBin(max_f1_s2t)); 
   sum_s2r[i]=sum_s2r[i]*(max_f1_s2t-min_f1_s2t)/bin_f1_s2t; //Bin Integral
 
 shift_rate[i]=inte_s2r[i]/sum_s2r[i];
 gs2r_shift->SetPoint(i,i,shift_rate[i]);
 gs2r_sig->SetPoint(i,i,sig_s2r[i]);
 gs2r_sig->SetPointError(i,0,sig_err_s2r[i]);
 gs2r_nevent->SetPoint(i,i,sum_s2r[i]);
 gs2r_nevent2->SetPoint(i,i,inte_s2r[i]);
   }
 




   //=============================================//
   //============= ReFill =======================//
   //============================================//



   
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
 Ls2pads=(int)Ls2tpads[0];
 Rs2pads=(int)Rs2tpads[0];
 rpathl=rtrpathl[0]+rs2pathl[0]; // R-HRS path length S2 -RF
 lpathl=ltrpathl[0]+ls2pathl[0]; // L-HRS path length S2 -RF

 if(Pion_flag==1){rbeta=ppi/Epi;
   if(k==1){cout<<"Pion mass beta"<<endl;}} //Pion
 else{rbeta=pk/Ek; 
   if(k==1){cout<<"Kaon mass beta"<<endl;}}//Kaon}
 rpath_corr=rpathl/rbeta/c;

 lbeta=1.0;//pe_/Ee_; 
 lpath_corr=lpathl/lbeta/c;
 Rs2_off=s2f1_off(Rs2pads,"R",mode,kine);
 Ls2_off=s2f1_off(Ls2pads,"L",mode,kine);
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 Ltof=tof_l-(-LF1[27]+LF1[28])/2.0*tdc_time;

  coin_t=tof_r-tof_l-coin_offset-s2_offset; //coin time
  coin_tc=tof_r-tof_l+rpath_corr-lpath_corr-pathl_off-coin_offset; //  coin Path & Offset  correction 
//====== Cut condition ========================// 
   cut_rpathl=false;
   cut_lpathl=false;
   cut_Rs2=false;
   cut_Ls2=false;
   cut_vz=false;
   cut_track=false;
   cut_s0=false;
   coin_trig=false;
   tri_T5=false;
   if(Rvz_cutmin<Rvz[0] && Rvz[0]<Rvz_cutmax && Lvz_cutmin<Lvz[0] && Lvz[0]<Lvz_cutmax)cut_vz=true;
   
  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads)cut_track=true;
  if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax)cut_rpathl=true;
  if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax)cut_lpathl=true;
  if(RF1[43]>0 && RF1[44]>0 && LF1[27]>0 && LF1[28]>0)cut_s0=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)cut_Rs2=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)cut_Ls2=true;
  if(cut_Rs2 && cut_Ls2 && cut_s0)coin_trig=true; 
  if(DRevtype==5)tri_T5=true; 
  

   double f1_s2r=(-LF1[Ls2pads]+LF1[30])*tdc_time;
   double f1_s2l=(-LF1[Ls2pads+48]+LF1[37])*tdc_time;


     if(coin_trig && cut_lpathl && cut_rpathl && tri_T5){
        hcoin_pi->Fill(coin_tc);
        if(Ls2pads<12 || (f1_s2r>mean_s2r[Ls2pads]-2.0 && f1_s2r<mean_s2r[Ls2pads]+2.0))
	   hcoin_pi_cut->Fill(coin_tc);
	if(f1_s2r>mean_s2r[Ls2pads]-2.0 && f1_s2r<mean_s2r[Ls2pads]+2.0)
	  hcoin_pi_cut_all->Fill(coin_tc);
	   }
     
}




   
   
   //===== Draw ========//

   
TLine* ltime=new TLine(min_coin_c,ac1_adc[0],max_coin_c,ac1_adc[0]);
 ltime->SetLineWidth(2);
 ltime->SetLineColor(2);
 //ltime->DrawLine();

 TCanvas* c4=new TCanvas("c4","hcoin piton cut");
 c4->cd();
 hcoin_pi->Draw();
 hcoin_pi_cut->SetLineColor(2);
 hcoin_pi_cut_all->SetLineColor(6);
 hcoin_pi_cut->Draw("same");
 hcoin_pi_cut_all->Draw("same");
 TCanvas* c5=new TCanvas("c5","c5");
 c5->Divide(4,4);
 TCanvas* c6=new TCanvas("c6","c6");
 c6->Divide(4,4);
 for(int i=0;i<16;i++){
 c5->cd(i+1);
 hcoin_rs2[i]->SetLineColor(kBlue);
 hcoin_rs2[i]->Draw();
 c6->cd(i+1);
 hcoin_ls2[i]->SetLineColor(kBlue);
 hcoin_ls2[i]->Draw();
 
}


 TCanvas* c7=new TCanvas("c7","c7");
 c7->Divide(2,1);
 c7->cd(1);
 gs2r_mean->Draw("AP");
 gs2r_mean->GetYaxis()->SetRangeUser(-0.2,0.2);
 c7->cd(2);
 gs2l_mean->Draw("AP");
 gs2l_mean->GetYaxis()->SetRangeUser(-0.2,0.2);

 TCanvas*c8=new TCanvas("c8","c8");
 c8->Divide(2,2);
 c8->cd(1);
 hLs2[12]->Draw("colz");
 c8->cd(2);
 hLs2[13]->Draw("colz");
 c8->cd(3);
 hLs2[14]->Draw("colz");
 c8->cd(4);
 hLs2[15]->Draw("colz");

 TCanvas*c9=new TCanvas("c9","c9");
 c9->Divide(2,2);
 c9->cd(1);
 hLs2_x[12]->Draw("colz");
 c9->cd(2);
 hLs2_y[12]->Draw("colz");
 c9->cd(3);
 hLs2_th[12]->Draw("colz");
 c9->cd(4);
 hLs2_ph[12]->Draw("colz");
 
 TCanvas*c10=new TCanvas("c10","c10");
 c10->cd();
 hLs2_path[12]->Draw("colz");

 TCanvas*c11=new TCanvas("c11","c11");
 c11->Divide(2,2);
 c11->cd(1);
 hLtof[12]->Draw();
 c11->cd(2);
 hLtof[13]->Draw();
 c11->cd(3);
 hLtof[14]->Draw();
 c11->cd(4);
 hLtof[15]->Draw();

 TCanvas*c12=new TCanvas("c12","c12");
 c12->Divide(2,2);
 c12->cd(1);
 hLtof_adc[12]->Draw("colz");
 c12->cd(2);
 hLtof_adc[13]->Draw("colz");
 c12->cd(3);
 hLtof_adc[14]->Draw("colz");
 c12->cd(4);
 hLtof_adc[15]->Draw("colz");

 TCanvas* c13=new TCanvas("c13","c13");
 c13->Divide(2,2);
 c13->cd(1);
 c13->cd(1)->SetLogz(1);
 hcoin_event[12]->GetYaxis()->SetRangeUser(0,10.);
 hcoin_event[12]->Draw("colz");
 c13->cd(2);
 c13->cd(2)->SetLogz(1);
 hcoin_event[13]->GetYaxis()->SetRangeUser(0,10.);
 hcoin_event[13]->Draw("colz");
 c13->cd(3);
 c13->cd(3)->SetLogz(1);
 hcoin_event[14]->GetYaxis()->SetRangeUser(0,10.);
 hcoin_event[14]->Draw("colz");
 c13->cd(4);
 c13->cd(4)->SetLogz(1);
 hcoin_event[7]->GetYaxis()->SetRangeUser(0,10.);
 hcoin_event[7]->Draw("colz");

 /*
 TCanvas* c14=new TCanvas("c14","c14");
 c14->Divide(2,2);
 c14->cd(1);
 gLS2_event[12]->GetYaxis()->SetRangeUser(0,10.);
 // gLS2_event[12]->GetXaxis()->SetRangeUser(0,100000.);
 gLS2_event[12]->SetFillColor(4);
 gLS2_event[12]->SetMarkerSize(0.5);
 gLS2_event[12]->SetMarkerColor(2);
 gLS2_event[12]->SetFillStyle(3005);
 
 gLS2_event[12]->Draw("AP");
 c14->cd(2);
 gLS2_event[13]->SetFillColor(6);
 gLS2_event[13]->SetMarkerSize(0.5);
 gLS2_event[13]->SetMarkerColor(6);
 gLS2_event[13]->SetFillStyle(3005);
 gLS2_event[13]->GetYaxis()->SetRangeUser(0,10.);
 gLS2_event[13]->Draw("AP");
 c14->cd(3);
 gLS2_event[14]->SetFillColor(8);
 gLS2_event[14]->SetMarkerSize(0.5);
 gLS2_event[14]->SetMarkerColor(8);
 gLS2_event[14]->SetFillStyle(3005);
 gLS2_event[14]->GetYaxis()->SetRangeUser(0,10.);
 gLS2_event[14]->Draw("AP");
  c14->cd(4);
 gLS2_event[7]->GetYaxis()->SetRangeUser(0,10.);
 gLS2_event[7]->Draw("AP");
 
 TCanvas* c15=new TCanvas("c15","c15");
 c15->cd();
 gLS2_event[7]->SetFillColor(2);
 gLS2_event[7]->SetMarkerSize(0.5);
 gLS2_event[7]->SetMarkerColor(2);
 gLS2_event[7]->SetFillStyle(3005);
 gLS2_event[7]->Draw("AP");
 gLS2_event[12]->Draw("P");
 gLS2_event[13]->Draw("P");
 gLS2_event[14]->Draw("P");
 cout<<"Drawing is done !"<<endl;
 */

 TCanvas* c16=new TCanvas("c16","c16");
 c16->Divide(2,2);
 c16->cd(1);
 hLs2r_event[12]->Draw("colz");
 c16->cd(2);
 hLs2r_event[13]->Draw("colz");
 c16->cd(3);
 hLs2r_event[14]->Draw("colz");
 c16->cd(4);
 hLs2r_event[7]->Draw("colz");

 TCanvas* c17=new TCanvas("c17","c17");
 c17->Divide(2,2);
 c17->cd(1);
 hLs2l_event[12]->Draw("colz");
 c17->cd(2);
 hLs2l_event[13]->Draw("colz");
 c17->cd(3);
 hLs2l_event[14]->Draw("colz");
 c17->cd(4);
 hLs2l_event[7]->Draw("colz");

 TCanvas* c18=new TCanvas("c18","c18");
 c18->Divide(2,2);
 c18->cd(1);
 hLs2r_f1_event[12]->Draw("colz");
 c18->cd(2);
 hLs2r_f1_event[13]->Draw("colz");
 c18->cd(3);
 hLs2r_f1_event[14]->Draw("colz");
 c18->cd(4);
 hLs2r_f1_event[7]->Draw("colz");
 TCanvas* c19=new TCanvas("c19","c19");
 c19->Divide(2,2);
 c19->cd(1);
 hLs2l_f1_event[12]->Draw("colz");
 c19->cd(2);
 hLs2l_f1_event[13]->Draw("colz");
 c19->cd(3);
 hLs2l_f1_event[14]->Draw("colz");
 c19->cd(4);
 hLs2l_f1_event[7]->Draw("colz");

 TCanvas* c20=new TCanvas("c20","c20");
 c20->Divide(2,2);
 c20->cd(1);
 hLs2_event[12]->Draw("colz");
 c20->cd(2);
 hLs2_event[13]->Draw("colz");
 c20->cd(3);
 hLs2_event[14]->Draw("colz");
 c20->cd(4);
 hLs2_event[7]->Draw("colz");

 TCanvas* c21=new TCanvas("c21","c21");
 c21->Divide(1,2);
 c21->cd(1);
 c21->SetLogz(1);
 hLs2r_f1_seg->Draw("colz");
 c21->cd(2);
 c21->SetLogz(1);
 hLs2l_f1_seg->Draw("colz");

TCanvas* c22=new TCanvas("c22","c22");
TCanvas* c23=new TCanvas("c23","c23"); 
TCanvas* c24=new TCanvas("c24","c24");
TCanvas* c25=new TCanvas("c25","c25");
TCanvas* c26=new TCanvas("c26","c26");
TCanvas* c27=new TCanvas("c27","c27"); 
 c22->Divide(4,4);
 c23->Divide(4,4);
 c24->Divide(4,4);
 c25->Divide(4,4);
 c26->Divide(4,4);
 c27->Divide(4,4);

 
 for(int i=0;i<16;i++){
   c22->cd(1+i);
   ltime->DrawLine(-2.+peak[i],0,-2.+peak[i],10000.);
   ltime->DrawLine(+2.+peak[i],0.,+2.+peak[i],10000.);
   hLs2r_f1_event[i]->Draw("colz");
   ///   line->Draw("same");
   c23->cd(i+1);
   hLs2l_f1_event[i]->Draw("colz");
   c24->cd(i+1);
   hLs2r_f1_tdc[i]->Draw();
   c25->cd(i+1);
   hLs2l_f1_tdc[i]->Draw();
   c26->cd(i+1);
   hLs2l_f1_tw[i]->Draw("colz");
   c27->cd(i+1);
   hLs2r_f1_tw[i]->Draw("colz");


 }
 
TCanvas* c28=new TCanvas("c28","c28");
 c28->cd();
 gs2r_shift->Draw("AP");

TCanvas* c29=new TCanvas("c29","c29");
 c29->cd();
 gs2r_sig->Draw("AP");

 TCanvas* c30=new TCanvas("c30","c30");
 c30->cd();
 gs2r_nevent->SetMarkerColor(2);
 gs2r_nevent2->SetMarkerColor(3);
 gs2r_nevent->Draw("AP");
 gs2r_nevent2->Draw("P");

 TCanvas* c31=new TCanvas("c31","c31"); 
 c31->cd();
 hcoin_R->Draw();
 TCanvas* c32=new TCanvas("c32","c32"); 
 c32->cd();
 hcoin_2->Draw();
	
  TCanvas* c33=new TCanvas("c33","c33"); 
  c33->Divide(2,2);
  c33->cd(1);
  hLs2r_y[7]->Draw("colz");
  c33->cd(2);
  hLs2r_y[12]->Draw("colz");	
  c33->cd(3);
  hLs2r_y[13]->Draw("colz");
  c33->cd(4);
  hLs2r_y[14]->Draw("colz");		  

 TCanvas* c34=new TCanvas("c34","c34"); 
 c34->cd();
 hy_event->Draw("colz");    

 TCanvas* c35=new TCanvas("c35","c35"); 
 c35->Divide(2,1);
 c35->cd(1);
 hLs2r_trig[7]->Draw("colz");
 c35->cd(2); 
 hLs2r_trig[12]->Draw("colz");
 TString name;
 if(output_flag){
 name.Form(ofname.c_str());
 // c0->Print(name+"[","pdf");//
 // c0->Print(name,"pdf");
 // c1->Print(name,"pdf");
 // c2->Print(name,"pdf");
 // c3->Print(name,"pdf");
 c4->Print(name+"[","pdf");
 c4->Print(name,"pdf");
 // c5->Print(name,"pdf");
 // c6->Print(name,"pdf");
 //c7->Print(name,"pdf");
 //c8->Print(name,"pdf");
 // c9->Print(name,"pdf");
 // c10->Print(name,"pdf");
 // c11->Print(name,"pdf");
 // c12->Print(name,"pdf");
 // c14->Print(name,"pdf");
 // c15->Print(name,"pdf"); //
 //c16->Print(name,"pdf");
  // c17->Print(name,"pdf");
 c18->Print(name,"pdf");
 c19->Print(name,"pdf");
 // c20->Print(name,"pdf");
 // c21->Print(name,"pdf");
 c22->Print(name,"pdf");
 c23->Print(name,"pdf");
 c24->Print(name,"pdf");
 c25->Print(name,"pdf"); 
 c26->Print(name,"pdf"); 
 c27->Print(name,"pdf");
 c28->Print(name,"pdf");  
 c29->Print(name,"pdf");
 c30->Print(name,"pdf");
 c31->Print(name,"pdf");
 c32->Print(name,"pdf");
 cout<<"c32 is drawn"<<endl;
 // c33->Print(name,"pdf");
 cout<<"c33 is drawn"<<endl; 
 c34->Print(name,"pdf");    
 cout<<"c34 is drawn"<<endl;
 c35->Print(name,"pdf");    
 c13->Print(name,"pdf");
 c13->Print(name +"]","pdf");
 cout<<"Print is done !"<<endl;
}


 // make root //
 if(root_flag){
   TFile* fout =new TFile(ofname.c_str(),"recreate");
     tnew->Write();
   for(int i=0;i<16;i++){
     
     hLs2l_fadc[i]->Write();
     hf1_fadc_l[i]->Write();
     hLs2r_fadc[i]->Write();
     hf1_fadc_r[i]->Write();

     
     hLs2l_fbus[i]->Write();
     hf1_fbus[i]->Write();
     hLs2r_fbus[i]->Write();
     hf1_fbus_r[i]->Write();

     hcoin_event[i]->Write();
     hLs2l_event[i]->Write();
     hLs2r_event[i]->Write(); 
     hLs2l_tdc[i]->Write();
     hLs2r_tdc[i]->Write(); 
     hLs2l_f1_event[i]->Write();
     hLs2r_f1_event[i]->Write();
     hLs2l_f1_tdc[i]->Write();
     hLs2r_f1_tdc[i]->Write();
     hLs2_f1_event[i]->Write();
     hLs2_event[i]->Write();
     hcoin_rs2[i]->Write();
     hcoin_ls2[i]->Write();
     hLs2l_f1_tw[i]->Write(); 
     hLs2r_f1_tw[i]->Write();
     hLs2r_y[i]->Write();

   }

     hLs2r_f1_seg->Write(); 
     hLs2l_f1_seg->Write();
     hcoin_pi->Write();
     hcoin_pi_cut->Write();
     hcoin_pi_cut_all->Write();
     gs2r_sig->SetName("gs2r_sig");
     gs2r_mean->Write();
     gs2r_sig->Write();
     gs2r_shift->SetName("gs2r_shift");
     gs2r_shift->Write();
     gs2r_nevent->Write();
     gs2r_nevent2->Write();
     hy_event->Write();    
   fout->Close();
 }




 
 for(int i=0;i<16;i++){
 cout<<Form("rmin_pi[%d]",i)<<rmean_pi[i]<<endl;
 cout<<Form("lmin_pi[%d]",i)<<lmean_pi[i]<<endl;
 cout<<Form("rmin_pi_diff[%d]-[7]",i)<<rmean_pi[i]-rmean_pi[7]<<endl;
 cout<<Form("lmin_pi_diff[%d]-[7]",i)<<lmean_pi[i]-lmean_pi[7]<<endl;
 cout<<Form("rs2[%d] t0 diff [ch",i)<<diff_r[i]<<endl;
 cout<<Form("ls2[%d] t0 diff [ch]: ",i)<<diff_l[i]<<endl;
 cout<<Form("rs2[%d] [ch]: ",i)<<rs2_off[i]<<endl;
 cout<<Form("ls2[%d] [ch]: ",i)<<ls2_off[i]<<endl;
 cout<<Form("r_pmax[%d]",i)<<r_pmax[i]<<endl; 
 cout<<Form("l_pmax[%d]",i)<<l_pmax[i]<<endl;
 cout<<Form("nr_pmax[%d]",i)<<nr_pmax[i]<<endl;
 cout<<Form("rsig_pi[%d]",i)<< rsig_pi[i]<<endl;

 }

 cout<<"====== parameters set ================="<<endl;
  cout<<"===== Rs2_off [1]======"<<endl;
 for(int i=0;i<16;i++){
  
   cout<<rs2_off[i]<<",";
}
 cout<<endl;
  cout<<"===== Ls2_off [1] ======"<<endl;
for(int i=0;i<16;i++){
  
   cout<<ls2_off[i]<<",";
}
 cout<<endl; 


 double sum_nev,sum_nev2,sum_cut;
 sum_nev=0.0;
 sum_nev2=0.0;
 sum_cut=0.0;
 cout<<"========== TDC shift study ========"<<endl;
 for(int i=0;i<16;i++){
   cout<<"i:"<<i<<endl;
   cout<<Form("Maximum point",i)<<peak[i]<<endl;
   cout<<"mean: "<<mean_s2r[i]<<endl;
   cout<<"sigma: "<<sig_s2r[i]<<endl;
   cout<<"Integral range "<<mean_s2r[i]-3*sig_s2r[i]<<"< tdc [ns] <"<<mean_s2r[i]+3*sig_s2r[i]<<endl;
   cout<<"shift rate [%]: "<<shift_rate[i]*100<<endl;
   cout<<"Fitting counts: "<<sum_s2r[i]<<endl;
   cout<<"Integral counts: "<<inte_s2r[i]<<endl;
   sum_nev=sum_nev+sum_s2r[i];
   sum_nev2=sum_nev2+inte_s2r[i];
   sum_cut=sum_cut+(sum_s2r[i]-inte_s2r[i]);
   }



 cout<<"Sum of Total Event (within 3 sigma): "<<sum_nev<<endl;
 cout<<"Sum of Total Event (within 2 ns)   : "<<sum_nev2<<endl;
 cout<<"Cut ecent (N_3sig-N_2ns)           : "<<sum_cut<<endl;
 

 gSystem->Exit(1);
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
