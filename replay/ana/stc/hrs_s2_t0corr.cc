using namespace std;
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
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
#include <TVector3.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>
#include <TText.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TColor.h>
#include <TPaveText.h>
#include <TRandom.h>
#include <TMinuit.h>
#include "Setting.h"
#include "Tree.h"


//=================================================//
//=============       Main     ====================//
//=================================================//


int main(int argc, char** argv){


  string ifname = "/data1/rootfiles/tritium_111721.root";
  string ofname = "ang_rhrs.root";
  string matrix_name="./matrix /momcalib_matrix.dat";
  string ofMTPname="./matrix/momcalib_test.dat";
  string opt_file="./scale_offset_20190210.dat";  
  string iteration;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = false;
  bool root_flag=false;
  bool matrix_flag=false;
  // bool RHRS_flag=true;
  bool RHRS_flag=false; 
  bool tuning_flag=false;
  string pngname;
  extern char *optarg;
  char* Target="H";
  string F1tdc="1";
  int f1tdc=1;


  
  
  while((ch=getopt(argc,argv,"h:s:w:t:p:f:x:r:y:m:o:O:i:RLbcop"))!=-1){
    switch(ch){
      
      
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;

    case 's':
      root_flag = true;
      draw_flag = false;
      single=true;
      ifname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;      

      
    case 'r':
      root_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;


    case 't':
      Target  = optarg;
      if(Target=="H")cout<<"Target : Hydrogen "<<endl;
      if(Target=="T")cout<<"Target : Tritium "<<endl;      
      break;

      
    case 'm':
      mode  = optarg;
      if(mode=="C")cout<<"Both arm momentum tuning "<<endl;
      if(Target=="L")cout<<"LHRS momentum tuning "<<endl;
      if(Target=="R")cout<<"RHRS momentum tuning "<<endl;
      if(mode=="CA")cout<<"Both arm momentum & RHRS angle tuning "<<endl;
      if(mode=="A")cout<<"RHRS angle tuning"<<endl;
      if(mode=="I"){
	cout<<"Both arm momentum tuning && Initial matrix"<<endl;
	Initial=true;
      }
      
      break;
      

    case 'O':
      F1tdc= optarg;
      f1tdc=atoi(F1tdc.c_str());
      //if F1tdc=1 tdc resolution 0.056
      //if F1tdc=2 tdc resolution 0.058
      //if F1tdc=3 tdc resolution 0.058 & Lp scale ON
      cout<<"F1 resolution mode : "<<f1tdc<<endl;      
      break;


      
    case 'w':
      ofMTPname  = optarg;
      cout<<"output new parameter filename : "<<ofMTPname<<endl;      
      matrix_flag=true;
      break;

    case 'i':
      iteration =optarg;
      tuning_flag=true;
      nite=atoi(iteration.c_str());
      cout<<"#iteration : "<<nite<<endl;
      break;
      
    case 'o':
      root_flag = true;
      draw_flag = false;
      matrix_flag=true;
      ofname = optarg;
      ofname =root_init+ ofname+ root_end;
      ofMTPname=optarg;
      ofMTPname = dat_init + ofMTPname;//+ dat_end;

      cout<<"output root filename : "<<ofname<<endl;      
      cout<<"output new parameter filename : "<<ofMTPname + ".dat"<<endl;      
      break;
      
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;

    case 'p':
      matrix_name =optarg;
      break;
      


    case 'h':
      cout<<"-f : input root  filename"<<endl;
      cout<<"-p : input matrix filename"<<endl;      
      cout<<"-w : output matrix filename"<<endl;
      cout<<"-r : output root filename"<<endl;
      cout<<"-t : target mode H:hydrogen T:trititum"<<endl;
      cout<<"-m : tuning mode C:R&L-HRS, L:L-HRS, R:R-HRS "<<endl;      
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




  
  TApplication *theApp =new TApplication("App",&argc,argv);
  gSystem->Load("libMinuit");

  momcalib* Mom=new momcalib();


  
   gSystem->Exit(1);
   theApp->Run();
 
  return 0;
  
}//end main

    //=======================================================//
    //===================    Class   ========================//



s2t0::s2t0(){
 set->Initialize();
  if( root_out ) ofp = new TFile(Form("%s",ofroot.c_str()),"recreate"); 
};

s2t0::~s2t0(){};

//////////////////////////////////////////////////////////////

void s2t0::SetRoot(string ifname){
  add_tree(ifname);
  pack_tree();
  readtreeHRSR();
  readtreeHRSL();
}

///////////////////////////////////////////////////////////////

void s2t0::SetRunList(string ifname){

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    add_tree(runname);
  }

  pack_tree();
  readtreeHRSR();
  readtreeHRSL();
}
//////////////////////////////////////////////////////////////

void MakeHist(){

  
  for(int i=0;i<s2seg;i++){
    hRs2_rt[i]=new TH1D(Form("hRs2_rt_%d",i),Form("RHRS S2 Seg%d R-PMT ; TDC [ns] ;Counts",i),bin_s2,min_s2,max_s2);
    hLs2_rt[i]=new TH1D(Form("hLs2_rt_%d",i),Form("LHRS S2 Seg%d R-PMT ; TDC [ns] ;Counts",i),bin_s2,min_s2,max_s2);
    hcoin_rs2_rt[i]=new TH2D(Form("hcoin_rs2_rt_%d",i),Form("RHRS S2 Seg%d R-PMT vs coin time; TDC [ns] ; coin time [ns] ",i),bin_);

    hcoin_ls2_rt[i]=new TH2D(Form("hcoin_ls2_rt_%d",i),Form("LHRS S2 Seg%d R-PMT vs coin time; TDC [ns] ; coin time [ns] ",i),bin_coin,min_coin,max_coin,bin_s2,min_s2,max_s2);

    hRs2_lt[i]=new TH1D(Form("hRs2_lt_%d",i),Form("RHRS S2 Seg%d L-PMT ; TDC [ns] ;Counts",i),bin_s2,min_s2,max_s2);
    hLs2_lt[i]=new TH1D(Form("hLs2_lt_%d",i),Form("LHRS S2 Seg%d L-PMT ; TDC [ns] ;Counts",i),bin_s2,min_s2,max_s2);
    hcoin_rs2_lt[i]=new TH2D(Form("hcoin_rs2_lt_%d",i),Form("RHRS S2 Seg%d L-PMT vs coin time; TDC [ns] ; coin time [ns] ",i),bin_);

    hcoin_ls2_lt[i]=new TH2D(Form("hcoin_ls2_lt_%d",i),Form("LHRS S2 Seg%d L-PMT vs coin time; TDC [ns] ; coin time [ns] ",i),bin_coin,min_coin,max_coin,bin_s2,min_s2,max_s2);    
    
  }

  hcoin=new TH2D("hcoin","Coin time ; coin time [ns];Counts",bin_coin,min_coin,max_coin);
 
}



void Fill(){

  ENum=T->GetEntries();
  cout<<"Events : "<<ENum<<endl;


  for(int k=0;k<ENum;k++){

    //==== Inititalization ===//
    rs2_seg=-1;
    ls2_Seg=-1;
    for(int i=0;i<s2seg;i++){
      rs2[i]=-100;
      ls2[i]=-100;
    }
    
    T->GetEntry(k);

    rs2_seg=(int)R_s2_t_pads[0];
    ls2_seg=(int)L_s2_t_pads[0];
    

    //    rs2_rt[rs2_seg]=-RTDC_F1_FirstHit[48+rs2_seg]+RTDC_F1_FirstHit[46];

    rs2_rt[rs2_seg]=-RTDC_F1_FirstHit[16+rs2_seg]+RTDC_F1_FirstHit[9];
    rs2_lt[rs2_seg]=-RTDC_F1_FirstHit[48+rs2_seg]+RTDC_F1_FirstHit[46];
    ls2_rt[rs2_seg]=-RTDC_F1_FirstHit[48+ls2_seg]+RTDC_F1_FirstHit[37];
    ls2_lt[rs2_seg]=-RTDC_F1_FirstHit[ls2_seg]+RTDC_F1_FirstHit[30];
    
  }





}
