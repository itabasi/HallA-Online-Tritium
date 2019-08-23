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
#include "angcalib.h"



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
  string matrix_xp="sample_matrix/newpar_xpt_2.dat";
  string matrix_yp="sample_matrix/newpar_ypt_2.dat";
  string opt_file="./scale_offset_20190210.dat";  
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = false;
  bool root_flag=false;
  // bool RHRS_flag=true;
  bool RHRS_flag=false; 
  bool fill_flag=false;
  string pngname;
  extern char *optarg;
  string itel;  
  //  char *root_init="/w/halla-scifs17exp/triton/itabashi/rootfiles/calib_root/";//ifarm
  string root_init="../rootfiles/angcalib/";
  string root_end=".root";
  string dat_init="../matrix/";
  string dat_end=".dat";
  string matrix="matrix/zt_RHRS_2.dat";
  
  while((ch=getopt(argc,argv,"h:f:i:r:x:y:t:o:s:RLdbcop"))!=-1){
    switch(ch){

    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;

    case 's':
      ifname = optarg;
      cout<<"input rootname : "<<ifname<<endl;
      break;

      
    case 't':
      ofMTPname = optarg;
      cout<<"output new parameter filename : "<<ofMTPname<<endl;      
      break;

    case 'i':
      fill_flag=true;
      itel= optarg;
      nite= atoi(itel.c_str());
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
    matrix_xp=optarg;
    cout<<"input Xp matrix filename: "<<matrix_xp<<endl;
    break;

    case 'y':
    matrix_yp=optarg;
    cout<<"input Yp matrix filename: "<<matrix_yp<<endl;
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

    case 'd':
      draw_flag =true;
      BreakTrue=false; 
      cout<<"Draw analysis"<<endl;
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
  angcalib * Ang=new angcalib();
  cout<<"root flag "<<root_flag<<endl;
   Ang->SetBranch(ifname,RHRS_flag);
   if(root_flag && draw_flag==0)Ang->NewBranch(ofname,RHRS_flag);
   Ang->HolePosi(RHRS_flag);
   Ang->MakeHist();
   Ang->Scale_corr(opt_file);
   Ang->Mxpt(matrix_xp);
   Ang->Mypt(matrix_yp);
   if(nite != 0 || draw_flag)Ang->EventSelect(RHRS_flag);
   //Ang->EventSelect(RHRS_flag);
   if(nite>0 && draw_flag==0)Ang->Tuning(ofMTPname);
   if(nite >=0 && draw_flag==0 && fill_flag) Ang->Fill(RHRS_flag);
   //Ang->Fill(RHRS_flag);
   if(draw_flag)Ang->Draw();
   //Ang->Draw();
   if (root_flag && draw_flag==0)Ang->Write();
   if(draw_flag==0)Ang->Close_tree();
   cout<<"input root: "<<ifname<<endl;
   cout<<"output root: "<<ofname<<endl;
   if(draw_flag==0)gSystem->Exit(1);
   theApp->Run();


  
  


  return 0;

}//end main




