const double c=299792458e-9;// [m/ns]
const double mk=493.7e-3;// Kaon mass [GeV/c^2]
const double me=0.511e-3;// electron mass [GeV/c^2] 
const double ml=1115.7e-3;//Lambda mass [GeV/c^2]
const double mn=939.6e-3; // neutron mass [GeV/c^2]
const double mpi=139.6e-3;// pion mass [GeV/c^2]

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
using namespace std;
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
#include "TPaveText.h"
#include "TRandom.h"




//=====================================================================//
//============================= Main =================================//
//===================================================================//

int main(int argc, char** argv){
  //------ Initial Parameters Setting ---------------//
  int ch; char* mode="H";
  int kine=1;// 1: hydrogen kinematics 2:tritium kinematics
  double tdc_time=58.0e-3;//[ns]
  string ifname = "/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/scripts/ita_scripts/run_list/Lambda_test.list";
  string ofname = "/pdf/hydro1_AC_eff_test.pdf";


  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"h:f:w:n:bcop:GHT12"))!=-1){
    switch(ch){
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
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
  
    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
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


 if(draw_flag==0)gROOT->SetBatch(1);


 //TChain //
  TChain* T=new TChain("T");

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> runname;
    T->Add(runname.c_str());
    //     cout<<buf<<endl;
  }

  int MAX=100;
  int R_Nu1_wire,R_Nu2_wire,R_Nv1_wire,R_Nv2_wire;
  int L_Nu1_wire,L_Nu2_wire,L_Nv1_wire,L_Nv2_wire;  
  double R_u1_wire[MAX],R_u2_wire[MAX],R_v1_wire[MAX],R_v2_wire[MAX];
  double L_u1_wire[MAX],L_u2_wire[MAX],L_v1_wire[MAX],L_v2_wire[MAX];
  int R_Nu1_wire,R_Nu2_wire,R_Nv1_wire,R_Nv2_wire;
  int L_Nu1_time,L_Nu2_time,L_Nv1_time,L_Nv2_time;  
  double R_u1_time[MAX],R_u2_time[MAX],R_v1_time[MAX],R_v2_time[MAX];
  double L_u1_time[MAX],L_u2_time[MAX],L_v1_time[MAX],L_v2_time[MAX];

  
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("Ndata.R.vdc.u1.wire"          ,1);
  T->Branch("Ndata.R.vdc.u1.wire"          ,&R_Nu1_wire,"R_Nu1_wire/I");  
  T->SetBranchStatus("Ndata.R.vdc.u2.wire"          ,1);
  T->Branch("Ndata.R.vdc.u2.wire"          ,&R_Nu2_wire,"R_Nu2_wire/I");    
  T->SetBranchStatus("Ndata.R.vdc.v1.wire"          ,1);
  T->Branch("Ndata.R.vdc.v1.wire"          ,&R_Nv1_wire,"R_Nv1_wire/I");  
  T->SetBranchStatus("Ndata.R.vdc.v2.wire"          ,1);
  T->Branch("Ndata.R.vdc.v2.wire"          ,&R_Nv2_wire,"R_Nv2_wire/I");  
  T->SetBranchStatus("R.vdc.u1.wire"          ,1);
  T->Branch("R.vdc.u1.wire"                   ,R_vdc_u1_wire     ,"R_vdc_u1_wire[Ndata.R.vdc.u1.wire]/D"  );  
  T->SetBranchStatus("R.vdc.u2.wire"          ,1);
  T->Branch("R.vdc.u2.wire"                   ,R_vdc_u2_wire     ,"R_vdc_u2_wire[Ndata.R.vdc.u2.wire]/D"  );    
  T->SetBranchStatus("R.vdc.v1.wire"          ,1);
  T->Branch("R.vdc.v1.wire"                   ,R_vdc_v1_wire     ,"R_vdc_v1_wire[Ndata.R.vdc.v1.wire]/D"  );    
  T->SetBranchStatus("R.vdc.v2.wire"          ,1);
  T->Branch("R.vdc.v2.wire"                   ,R_vdc_v2_wire     ,"R_vdc_v2_wire[Ndata.R.vdc.v2.wire]/D"  );  

  T->SetBranchStatus("Ndata.R.vdc.u1.time"          ,1);
  T->Branch("Ndata.R.vdc.u1.time"          ,&R_Nu1_time,"R_Nu1_time/I");  
  T->SetBranchStatus("Ndata.R.vdc.u2.time"          ,1);
  T->Branch("Ndata.R.vdc.u2.time"          ,&R_Nu2_time,"R_Nu2_time/I");    
  T->SetBranchStatus("Ndata.R.vdc.v1.time"          ,1);
  T->Branch("Ndata.R.vdc.v1.time"          ,&R_Nv1_time,"R_Nv1_time/I");  
  T->SetBranchStatus("Ndata.R.vdc.v2.time"          ,1);
  T->Branch("Ndata.R.vdc.v2.time"          ,&R_Nv2_time,"R_Nv2_time/I");  
  T->SetBranchStatus("R.vdc.u1.time"          ,1);
  T->Branch("R.vdc.u1.time"                   ,R_vdc_u1_time     ,"R_vdc_u1_time[Ndata.R.vdc.u1.time]/D"  );  
  T->SetBranchStatus("R.vdc.u2.time"          ,1);
  T->Branch("R.vdc.u2.time"                   ,R_vdc_u2_time     ,"R_vdc_u2_time[Ndata.R.vdc.u2.time]/D"  );    
  T->SetBranchStatus("R.vdc.v1.time"          ,1);
  T->Branch("R.vdc.v1.time"                   ,R_vdc_v1_time     ,"R_vdc_v1_time[Ndata.R.vdc.v1.time]/D"  );    
  T->SetBranchStatus("R.vdc.v2.time"          ,1);
  T->Branch("R.vdc.v2.time"                   ,R_vdc_v2_time     ,"R_vdc_v2_time[Ndata.R.vdc.v2.time]/D"  );  



  

  //======== LHRS ===========//

  T->SetBranchStatus("Ndata.L.vdc.u1.wire"          ,1);
  T->Branch("Ndata.L.vdc.u1.wire"          ,&L_Nu1_wire,"L_Nu1_wire/I");  
  T->SetBranchStatus("Ndata.L.vdc.u2.wire"          ,1);
  T->Branch("Ndata.L.vdc.u2.wire"          ,&L_Nu2_wire,"L_Nu2_wire/I");    
  T->SetBranchStatus("Ndata.L.vdc.v1.wire"          ,1);
  T->Branch("Ndata.L.vdc.v1.wire"          ,&L_Nv1_wire,"L_Nv1_wire/I");  
  T->SetBranchStatus("Ndata.L.vdc.v2.wire"          ,1);
  T->Branch("Ndata.L.vdc.v2.wire"          ,&L_Nv2_wire,"L_Nv2_wire/I");  
  T->SetBranchStatus("L.vdc.u1.wire"          ,1);
  T->Branch("L.vdc.u1.wire"                   ,L_vdc_u1_wire     ,"L_vdc_u1_wire[Ndata.L.vdc.u1.wire]/D"  );  
  T->SetBranchStatus("L.vdc.u2.wire"          ,1);
  T->Branch("L.vdc.u2.wire"                   ,L_vdc_u2_wire     ,"L_vdc_u2_wire[Ndata.L.vdc.u2.wire]/D"  );    
  T->SetBranchStatus("L.vdc.v1.wire"          ,1);
  T->Branch("L.vdc.v1.wire"                   ,L_vdc_v1_wire     ,"L_vdc_v1_wire[Ndata.L.vdc.v1.wire]/D"  );    
  T->SetBranchStatus("L.vdc.v2.wire"          ,1);
  T->Branch("L.vdc.v2.wire"                   ,L_vdc_v2_wire     ,"L_vdc_v2_wire[Ndata.L.vdc.v2.wire]/D"  );  

  T->SetBranchStatus("Ndata.L.vdc.u1.time"          ,1);
  T->Branch("Ndata.L.vdc.u1.time"          ,&L_Nu1_time,"L_Nu1_time/I");  
  T->SetBranchStatus("Ndata.L.vdc.u2.time"          ,1);
  T->Branch("Ndata.L.vdc.u2.time"          ,&L_Nu2_time,"L_Nu2_time/I");    
  T->SetBranchStatus("Ndata.L.vdc.v1.time"          ,1);
  T->Branch("Ndata.L.vdc.v1.time"          ,&L_Nv1_time,"L_Nv1_time/I");  
  T->SetBranchStatus("Ndata.L.vdc.v2.time"          ,1);
  T->Branch("Ndata.L.vdc.v2.time"          ,&L_Nv2_time,"L_Nv2_time/I");  
  T->SetBranchStatus("L.vdc.u1.time"          ,1);
  T->Branch("L.vdc.u1.time"                   ,L_vdc_u1_time     ,"L_vdc_u1_time[Ndata.L.vdc.u1.time]/D"  );  
  T->SetBranchStatus("L.vdc.u2.time"          ,1);
  T->Branch("L.vdc.u2.time"                   ,L_vdc_u2_time     ,"L_vdc_u2_time[Ndata.L.vdc.u2.time]/D"  );    
  T->SetBranchStatus("L.vdc.v1.time"          ,1);
  T->Branch("L.vdc.v1.time"                   ,L_vdc_v1_time     ,"L_vdc_v1_time[Ndata.L.vdc.v1.time]/D"  );    
  T->SetBranchStatus("L.vdc.v2.time"          ,1);
  T->Branch("L.vdc.v2.time"                   ,L_vdc_v2_time     ,"L_vdc_v2_time[Ndata.L.vdc.v2.time]/D"  );  

  


  int evnt=T->GetEntries();

  double min_time = -0.1e-6;
  double max_time =  0.5e-6;
  int bin_time=300;
  int min_wire = 0;
  int max_wire = 400;
  bin_wire=max_wire-min_wire;
  
  TH2F* hLu1_time_wire=new TH2F("hLu1_time_wire","",bin_time,min_time,max_time,bin_wire,min_mire,max_wire);
  double time;
  int wire;
  
  for(int i=0;i<evnt;i++){

    wire=-100;
    time=-100;
    T->GetEntry(i);

    wire=
    
    

  }









}
