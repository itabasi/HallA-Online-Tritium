#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
using namespace std;
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
#include "momcalib.h"
#include "Setting.h"
#include "define.h"

//===================================================================//
//============================= Main ================================//
//==================================================================//

int main(int argc, char** argv){


  gStyle->SetOptFit(111111111);
  int ch; char* mode="C";
  //  double tdc_time=58.0e-3;//[ns]
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
  //  char *root_init="/w/halla-scifs17exp/triton/itabashi/rootfiles/calib_root/";//ifarm
  string root_init="../rootfiles/";
  string root_end=".root";
  string dat_init="../matrix/";
  string dat_end=".dat";
  string matrix="matrix/zt_RHRS_2.dat";
  
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
  Mom->mode(mode,Target,f1tdc);
  Mom->MTParam(matrix_name);
  if(single)Mom->SingleRoot(ifname);
  else   Mom->SetRoot(ifname);
  if(root_flag)Mom->NewRoot(ofname);
  Mom->MakeHist();
  if(tuning_flag)Mom->EventSelection();
  if(tuning_flag && nite>0)Mom->MomTuning(ofMTPname);
  if(root_flag && ( (tuning_flag && nite>0) || tuning_flag==0)) Mom->Fill();  
  if(root_flag) Mom->Close();  


  cout<<endl;
  cout<<"============== Input files ==============="<<endl;
  cout<<"run files   : "<<ifname<<endl;
  cout<<"matrix param : "<<matrix_name<<endl;
  cout<<"analysis mode : "<<"Coincidence mode "<<endl;
  cout<<endl;
  cout<<"============== Output files ==============="<<endl;
  cout<<"root_flag "<<root_flag<<endl;
  if(root_flag)cout<<"Root file   : "<<ofname<<endl;
  if(matrix_flag)cout<<"New matrix file : "<<ofMTPname<<endl;
  cout<<"==========================================="<<endl;
  cout<<endl;


  
   gSystem->Exit(1);
   theApp->Run();

 
  return 0;
  
}//end main

