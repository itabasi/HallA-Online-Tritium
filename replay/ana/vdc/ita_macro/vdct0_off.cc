#include <iostream>
#include <fstream>
using namespace std;
#include "VDCt0.h"
#include "TApplication.h"
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
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
#include "VDCt0.h"


int main(int argc, char** argv){

  int ch;
  int runnum;
  string ifname ="";
  string ofname = "./param/test.dat";
  bool root_flag=false;
  bool print_flag=false;
  bool draw_flag=false;
  bool batch_flag=true;
  string root_name;
  string print_name;
   extern char *optarg;
   //  extern string optarg;   
  while((ch=getopt(argc,argv,"h:f:w:s:n:i:r:p:o:bcop"))!=-1){
    switch(ch){
            


    case 'i':
      runnum=atoi(optarg);
      cout<<"Run num : "<<runnum<<endl;
      break;      

    case 'h':
      cout<<"-i :Run number"<<endl;
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

  
  string param_init;
  string param_end;
  string root_init;
  string root_end;
  

  param_init="./param/initial_2run/";
  param_end=".dat";
  root_init="../../rootfiles/VDC/initial_2run/";
  root_end=".root";
  ostringstream run;
  run<<runnum<<"-"<<runnum+1;
  string def_param = "./param/def_t0.dat";

  string paraname = param_init + run.str() + param_end;
  root_name = root_init  + run.str() + "_t0tuned" + root_end;
  cout<<"input param file : "<<paraname<<endl;
  cout<<"output root file : "<<root_name<<endl;



  
  TApplication *theApp =new TApplication("App",&argc,argv);  
  VDCt0* vdct0=new VDCt0();


   gROOT->SetBatch(1);
   vdct0->Deft0(def_param);   
   vdct0->SetRun(runnum);
   vdct0->GetOffset(paraname);
   vdct0->NewRoot(root_name);
   vdct0->SetBranch();
   vdct0->MakeHist();
   vdct0->Fill();
   vdct0->SetRawt0();   
   vdct0->MakeRoot(root_name);

   cout<<endl;
   cout<<"========= output file =============="<<endl;
   cout<<"rootfile: "<<root_name<<endl;



  
  
  gSystem->Exit(1);
  theApp->Run();

  return 0;
  
  


  }
