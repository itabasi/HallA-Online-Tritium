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


//=======================================================//
//================     Main       =======================//
//=======================================================//


int main(int argc, char** argv){

  gStyle->SetOptFit(111111111);
  int ch;
  string ifname="";
  string ofname = "./param/test.dat";
  string   print_name, root_name;
  bool sigle_flag=false;
  bool runlist_flag=false;
  //  istringstream runnum;
   int runnum; 
  
   extern char *optarg;
   //  extern string optarg;   
  while((ch=getopt(argc,argv,"h:f:w:s:n:i:r:o:m:bcop"))!=-1){
    switch(ch){

    case 'i':
      runnum=atoi(optarg);
      cout<<"Run num : "<<runnum<<endl;
      break;      
      
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      runlist_flag=true;
      break;
    case 'o':
      ofname = optarg;
      root_name="./../../rootfiles/VDC/t0tuned/" + ofname + ".root";
      print_name="./../../pdf/VDC/ita_mac/t0tuned/" +ofname + ".pdf";
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
  string print_init;
  string param_end;
  string print_end;
  string root_init;
  string root_end;
  param_init="./param/t0tuned/";
  print_init="../../pdf/VDC/ita_mac/t0tuned";
  root_init="../../rootfiles/VDC/t0tuned";
  print_end=".pdf";
  param_end=".dat";
  root_end=".root";
  ostringstream run;
  run<<runnum<<"-"<<runnum+9;
  string def_param = "./param/def_t0.dat";
  string param_name = param_init + run.str() + param_end;

  if(runlist_flag==0){
  print_name = print_init + run.str() + print_end;
  root_name = root_init + run.str() + root_end;

  cout<<"print file: "<<print_name<<endl;
  cout<<"param file: "<<param_name<<endl;
  }
 TApplication *theApp =new TApplication("App",&argc,argv);  
  
   VDCt0* vdct0= new VDCt0();
   vdct0->Deft0(def_param);   
   if(runlist_flag) vdct0->SetRunList(ifname);
   else   vdct0->SetRun(runnum);
   vdct0->SetBranch();
   vdct0->MakeHist();
   vdct0->Fill();
   vdct0->Write(param_name);
   vdct0->Draw();
   vdct0->Print(print_name);
   vdct0->MakeRoot(root_name);
   gROOT->SetBatch(1);
   gSystem->Exit(1);
   theApp->Run();

 
 return 0;

}


