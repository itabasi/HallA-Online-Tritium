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
  string param_name ;
  bool sigle_flag=false;
  bool runlist_flag=false;
  bool root_flag=false;
  bool print_flag=false;
  bool param_flag=false;
  //  istringstream runnum;
   int runnum; 
   int Nrun=0;
   bool test_flag=false;
   extern char *optarg;
   //  extern string optarg;   
  while((ch=getopt(argc,argv,"h:f:w:s:n:i:r:o:m:Tbcop"))!=-1){
    switch(ch){

    case 'i':
      runnum=atoi(optarg);
      cout<<"Run num : "<<runnum<<endl;
      root_flag=true;
      print_flag=true;
      runlist_flag=false;
      break;      

    case 'n':
      Nrun=atoi(optarg);
      cout<<"Branch num "<<Nrun<<endl;
      break;
      
    case 'T':
      cout<<"Test mode !"<<endl;
      test_flag=true;
      break;
      
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      runlist_flag=true;
      break;

    case 'r':
      root_name = optarg;
      root_flag=true;
      print_flag=false;
      cout<<"output rootname : "<<root_name<<endl;
      break;

    case 'w':
      print_name = optarg;
      print_flag=true;
      cout<<"output pdfname : "<<print_name<<endl;
      break;      
      
    case 'o':
      ofname = optarg;
      print_flag=true;
      root_flag=true;
      param_flag=true;
      root_name="./../../rootfiles/VDC/t0tuned_all/" + ofname + ".root";
      print_name="./../../pdf/VDC/ita_mac/t0tuned_all/" +ofname + ".pdf"; 
      param_name="./param/t0tuned_all/"+ofname;//+ ".dat";     
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
  
  //  param_init="./param/t0tuned/";
  //  print_init="../../pdf/VDC/ita_mac/t0tuned";
  //  root_init="../../rootfiles/VDC/t0tuned";
  
  //  param_init="./param/initial/";
  //  print_init="../../pdf/VDC/ita_mac/initial/";
  //  root_init="../../rootfiles/VDC/initial/";

  //  param_init="./param/t0tuned_1ns/";
  //  print_init="../../pdf/VDC/ita_mac/t0tuned_1ns/";
  //  root_init="../../rootfiles/VDC/t0tuned_1ns/";

  param_init=Form("./param/t0tuned_%dns/",Nrun  );
  print_init=Form("../../pdf/VDC/ita_mac/t0tuned_%dns/",Nrun);
  root_init=Form("../../rootfiles/VDC/t0tuned_%dns/",Nrun);    

  //======== T0 tuned (all) ========//  
  param_init=Form("./param/t0tuned_all/",Nrun  );
  print_init=Form("../../pdf/VDC/ita_mac/t0tuned_all/",Nrun);
  root_init=Form("../../rootfiles/VDC/t0tuned_all/",Nrun);    

  
  if(test_flag){
    param_init="./param/test/";
    print_init="../../pdf/VDC/ita_mac/test/";
    root_init="../../rootfiles/VDC/test/";
  }
  
  print_end=".pdf";
  param_end=".dat";
  root_end=".root";
  
  ostringstream run;
  run<<runnum<<"-"<<runnum+Nrun-1;
  string def_param = "./param/def_t0.dat";
  //  if(param_flag)  param_name+=param_end;


  
  if(runlist_flag==0){
  print_name = print_init + run.str() + print_end;
  root_name = root_init + run.str() + root_end;
  param_name = param_init + run.str() + param_end;
  cout<<"root file: "<<root_name<<endl;
  cout<<"print file: "<<print_name<<endl;
  //  cout<<"param file: "<<param_name<<endl;
  
  }

  
 TApplication *theApp =new TApplication("App",&argc,argv);  

 gStyle->SetOptStat(000001111);
   VDCt0* vdct0= new VDCt0();
    gROOT->SetBatch(1);
   vdct0->Deft0(def_param);   
   if(runlist_flag) vdct0->SetRunList(ifname);
   else   vdct0->SetRun(runnum,Nrun);
   vdct0->SetBranch();
   vdct0->MakeHist();
   vdct0->Fill();
   vdct0->Sett0();   
   if(param_flag)   vdct0->Write_time(param_name);
   vdct0->Draw();
   if(print_flag)vdct0->Print_c(print_name);
   if(root_flag)vdct0->MakeRoot_c(root_name);

   cout<<"======= Output File ========"<<endl;
   if(root_flag)  cout<<"root files : "<<root_name<<endl;
   if(print_flag) cout<<"pdf name : "<<print_name<<endl;
   if(param_flag) cout<<"param name : "<<param_name +".dat"<<endl;

   gSystem->Exit(1);
   theApp->Run();

   
 
 return 0;

}


