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
#include "Setting.h"
#include "Param.h"
#include "rascalib.h"


string ofname("output.pdf");
string ofroot("output.root");

//===================================================================//
//============================= Main ================================//
//==================================================================//

int main(int argc, char** argv){


  gStyle->SetOptFit(111111111);
  int ch; char* mode;
  double tdc_time=58.0e-3;//[ns]
  string ifname = "/data1/rootfiles/tritium_111721.root";
  string ofname = "ang_rhrs.root";
  string matrix_name="./matrix /zt_RHRS_opt.dat";
  string ofMTPname="./matrix/zt_RHRS_opt_new.dat";
  string matrix_x="sample_matrix/newpar_xpt_2.dat";
  string matrix_y="sample_matrix/newpar_ypt_2.dat";
  string opt_file="./scale_offset_20190210.dat";
  string ofparam;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = false;
  bool root_flag=false;
  bool calib_x=false;
  bool calib_y=false;
  // bool RHRS_flag=true;
  bool RHRS_flag=false; 
  string pngname;
  extern char *optarg;
  //  char *root_init="/w/halla-scifs17exp/triton/itabashi/rootfiles/calib_root/";//ifarm
  string root_init="../rootfiles/";
  
  string root_end=".root";
  string dat_init="../matrix/";
  string dat_end=".dat";
  string matrix="matrix/zt_RHRS_2.dat";
  string itel;
  
  while((ch=getopt(argc,argv,"h:f:x:y:o:r:i:RLbcop"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
      
    case 't':
      ofMTPname = optarg;
      cout<<"output new parameter filename : "<<ofMTPname<<endl;      
      break;

    case 'i':

      itel= optarg;
      nite= atoi(itel.c_str());
      //  cout<<"nite: "<<nite<<endl;
      break;


    case 'r':
      root_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;
      
    case 'o':
      root_flag = true;
      //      cout<<"root flag: "<<root_flag<<endl;
      draw_flag = false;
      ofname = optarg;
      ofname =root_init+ ofname+ root_end;
      ofMTPname=optarg;
      ofMTPname = dat_init + ofMTPname;//+ dat_end;


      cout<<"output root filename : "<<ofname<<endl;      
      cout<<"output new parameter filename : "<<ofMTPname<<endl;      
      break;

    case 'x':
    matrix_x=optarg;
    calib_x=true;
    cout<<"input X matrix filename: "<<matrix_x<<endl;
    break;

    case 'y':
    matrix_y=optarg;
    calib_y=true;
    cout<<"input Y matrix filename: "<<matrix_y<<endl;
    break;

      
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
  

    case 'R':
      RHRS_flag = true;
      //     matrix_name = optarg;   
     cout<<"R-HRS analysis"<<endl;
      break;

    case 'L':
      RHRS_flag =false;
      //      matrix_name = optarg;
      cout<<"L-HRS analysis"<<endl;
      break;


    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
      cout<<"-RL: input HRS R or L Matrix Parameter name"<<endl;
      cout<<"-r : output root filename"<<endl;
      cout<<"-t : output tuning matrix parameter name"<<endl;
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

  rascalib * Ras=new rascalib();

  Ras->SetRoot(ifname,RHRS_flag);
  Ras->MakeHist();
  Ras->NewRoot(ofname);
  if(calib_x)Ras->MTParam_x(matrix_x,RHRS_flag);
  if(calib_y)Ras->MTParam_y(matrix_y,RHRS_flag);
  Ras-> EventSelection(RHRS_flag, calib_x);
  Ras->Tuning(ofMTPname);
  Ras->Fill(RHRS_flag, calib_x);
  if(root_flag)Ras->Close();

   gSystem->Exit(1);
   theApp->Run();


  
  return 0;

}//end main



