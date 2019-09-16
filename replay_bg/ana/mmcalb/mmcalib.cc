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

#include "mmcalib.h"
#include "Setting.h"
#include "define.h"

//===================================================================//
//============================= Main ================================//
//==================================================================//

int main(int argc, char** argv){


  gStyle->SetOptFit(111111111);
  int ch; char* mode;
  string ifname = "/data1/rootfiles/tritium_111721.root";
  string ofname = "ang_rhrs.root";
  string matrix_name = "./matrix/matrix_param.list";
  bool output_flag = false;
  bool output_tree_flag = false;
  bool root_flag=false;
  bool RHRS_flag=false;
  string pngname;
  extern char *optarg;
  

  while((ch=getopt(argc,argv,"hb:f:x:y:r:m:RLbcop"))!=-1){
    switch(ch){

    case 'R':
      RHRS= true;
      break;
    case 'L':
      LHRS= true;
      break;

    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
      
    case 'r':
      root_flag = true;
      ofname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;

    case 'm':
      matrix_name =optarg;
      break;

      
    case 'b':
      cout<<"BACH MODE!"<<endl;
      break;
  

      


    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
      cout<<"-RL: input HRS R or L"<<endl;
      cout<<"-r : output root filename"<<endl;
      cout<<"-m : input matrix parameters"<<endl;
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


   mmcalib* mmc=new mmcalib();
   mmc->SetRunList(ifname);
   mmc->NewRoot(ofname,ifname);
   mmc->MakeHist();
   mmc->MTParam(matrix_name);
   mmc->Fill();
   mmc->Close_tree();


   cout<<"========================================="<<endl;
   cout<<"=========== OutPut information =========="<<endl;
   cout<<"========================================="<<endl;
   cout<<"Set Run List: "<<ifname<<endl;
   cout<<"Input Matrix: "<<matrix_name<<endl;
   cout<<"New Root File: "<<ofname<<endl;

   
    gSystem->Exit(1);
   theApp->Run();

  return 0;

}//end main

