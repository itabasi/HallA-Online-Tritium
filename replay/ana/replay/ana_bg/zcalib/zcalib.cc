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
#include "zcalib.h"

string ofname("output.pdf");
string ofroot("output.root");

//===================================================================//
//============================= Main ================================//
//==================================================================//

int main(int argc, char** argv){


  gStyle->SetOptFit(111111111);
  int ch; char* mode;
  double tdc_time=58.0e-3;//[ns]
  string ifname = "/data1/rootfiles/tritium_111132.root";
  string ofname = "zt_rhrs.root";
  string matrix_name="./matrix/zt_RHRS_opt.dat";
  string ofMTPname="./matrix/zt_RHRS_opt_new.dat";
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = false;
  bool root_flag=false;
  // bool RHRS_flag=true;
  bool RHRS_flag=true; 
  string pngname;
  extern char *optarg;
  char *root_init="/w/halla-scifs17exp/triton/itabashi/rootfiles/calib_root/";//ifarm
  root_init="../rootfiles/";
  
  string root_end=".root";
  string dat_init="../matrix/";
  string dat_end=".dat";
  string matrix="matrices/zt_RHRS_2.dat";
  string itel;
  bool single_flag=false;
  while((ch=getopt(argc,argv,"h:s:f:w:n:r:m:R:L:t:i:o:bcop"))!=-1){
    switch(ch){
      
      
    case 's':
      ifname = optarg;
      cout<<"input Root : "<<ifname<<endl;
      single_flag=true;
      break;

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

    case 'r':
      root_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;

    case 'm':
      matrix_name =optarg;
      cout<<"input matrx parameter filename: "<<matrix_name<<endl;
      
    case 't':
      ofMTPname = optarg;
      cout<<"output new parameter filename : "<<ofMTPname<<endl;      
      break;

    case 'o':
      root_flag = true;
      cout<<"root flag: "<<root_flag<<endl;
      draw_flag = false;
      ofname = optarg;
      ofname =root_init+ ofname+ root_end;
      ofMTPname=optarg;
      ofMTPname = dat_init + ofMTPname;

    case 'i':

      itel= optarg;
      nite= atoi(itel.c_str());
      //  cout<<"nite: "<<nite<<endl;
      break;


    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
  

    case 'R':
      RHRS_flag = true;
      matrix_name = optarg;   
      cout<<"R-HRS analysis"<<endl;
      break;

    case 'L':
      RHRS_flag =false;
      matrix_name = optarg;
      cout<<"L-HRS analysis"<<endl;
      break;


    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
      cout<<"-m : input matrix filename"<<endl;
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


  zcalib * Zcalib=new zcalib();

   TApplication *theApp =new TApplication("App",&argc,argv);
 gSystem->Load("libMinuit"); 

 if(single_flag) Zcalib->SetRoot(ifname);
   else Zcalib->ChainRoot(ifname);

  if(root_flag)Zcalib->NewBranch(ofname, RHRS_flag);
  Zcalib->Mzt(matrix_name, RHRS_flag);
  Zcalib->Mzt_L(matrix_name, RHRS_flag);
  Zcalib->MakeHist();
  Zcalib->EventSelection(RHRS_flag);
  Zcalib->MTtuning(ofMTPname,matrix_name);
  Zcalib->Fill(RHRS_flag);
  if(draw_flag)Zcalib->Draw();
  if(root_flag)Zcalib->MakeTree();
  Zcalib->Close_tree();
  gSystem->Exit(1);
  theApp->Run();



  return 0;

}//end main




