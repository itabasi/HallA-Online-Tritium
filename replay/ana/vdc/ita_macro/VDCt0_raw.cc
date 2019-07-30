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
  string param_name;
  bool sigle_flag=false;
  bool runlist_flag=false;
  bool outname=false;
  //  istringstream runnum;
   int runnum; 
  
   extern char *optarg;
   //  extern string optarg;   
  while((ch=getopt(argc,argv,"h:f:w:s:n:i:r:o:m:bcop"))!=-1){
    switch(ch){

    case 'i':
      runnum=atoi(optarg);
      cout<<"Run num : "<<runnum<<endl;
      runlist_flag=false;
      outname=true;
      break;      
      
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      runlist_flag=true;
      break;
    case 'o':
      ofname = optarg;

      root_name="./../../rootfiles/VDC/t0ttuned_1ns/" + ofname + ".root";
      print_name="./../../pdf/VDC/ita_mac/t0tuned_1ns/" +ofname + ".pdf";
      param_name="./param/t0tuned_1ns/" + ofname;
      
      //      root_name="./../../rootfiles/VDC/initial/" + ofname + ".root";
      //      print_name="./../../pdf/VDC/ita_mac/initial/" +ofname + ".pdf";
      //      param_name="./param/initial/" + ofname;
      outname=false;
      break;
      
      
    case 'h':
      cout<<"-i :Run number"<<endl;
      cout<<"-o :Set File names (root & pdf & param)"<<endl;
      cout<<"Example :  ./bin/VDCt0_raw -i #run -o filename"<<endl;
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
  
  //  param_init="./param/t0tuned/";
  //  print_init="../../pdf/VDC/ita_mac/t0tuned";
  //  root_init="../../rootfiles/VDC/t0tuned";
  //  param_init="./param/";
  //  print_init="../../pdf/VDC/ita_mac/initial/";
  //  root_init="../../rootfiles/VDC/initial/";
  param_init="./param/t0tuned_1ns/";
  print_init="../../pdf/VDC/ita_mac/t0tuned_1ns/";
  root_init="../../rootfiles/VDC/t0tuned_1ns/";    
  print_end=".pdf";
  param_end=".dat";
  root_end=".root";
  ostringstream run;
  run<<runnum<<"-"<<runnum+9;
  string def_param = "./param/def_t0.dat";


 
  if(outname==true){
  print_name = print_init + run.str() + print_end;
  root_name = root_init + run.str() + root_end;
  param_name = param_init + run.str();// + param_end;

  }

  cout<<"output root file: "<<root_name<<endl;
  cout<<"output print file: "<<print_name<<endl;
  cout<<"output param file: "<<param_name+".dat"<<endl;

  
 TApplication *theApp =new TApplication("App",&argc,argv);  
 VDCt0* vdct0= new VDCt0();
   vdct0->Deft0(def_param);   
   if(runlist_flag) vdct0->SetRunList(ifname);
   else   vdct0->SetRun(runnum);
   gROOT->SetBatch(1);
   vdct0->SetBranch();
   vdct0->MakeHist();
   vdct0->Fill();
   vdct0->SetRawt0();   
   vdct0->Write(param_name);
   vdct0->Draw();
   vdct0->Print(print_name);
   vdct0->MakeRoot(root_name);
   cout<<endl;
   cout<<"========= output file =============="<<endl;
   cout<<"rootfile: "<<root_name<<endl;
   cout<<"pdffile: "<<print_name<<endl;
   cout<<"paramfile : "<<param_name<<endl;
   gSystem->Exit(1);
   theApp->Run();
 
 return 0;

}


