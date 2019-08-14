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
#include "VDCt0_plot.h"

//=======================================================//
//================     Main       =======================//
//=======================================================//


int main(int argc, char** argv){

  //  gStyle->SetOptStat(stat);
  
  int ch;
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
  while((ch=getopt(argc,argv,"h:f:w:s:n:i:r:o:bcop"))!=-1){
    switch(ch){
            
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;

    case 'o':
      root_flag=true;
      draw_flag=false;
      print_flag=true;
      ofname = optarg;
      root_name="./../../rootfiles/VDC/" + ofname + ".root";
      print_name="./../../pdf/VDC/ita_mac/" +ofname + ".pdf";
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

  

 TApplication *theApp =new TApplication("App",&argc,argv);  
  
   VDCt0_plot* plot= new VDCt0_plot();
   //   plot->SetPoint(ifname);
   if(draw_flag)    gROOT->SetBatch(1);
   plot->SetPointError(ifname);   
   plot->Draw();
   if(print_flag)   plot->Print(print_name);
   if(root_flag)    plot->MakeRoot(root_name);

   cout<<"pdf file "<<print_name<<endl;
   cout<<"root file "<<root_name<<endl;   
   gSystem->Exit(1);
   theApp->Run();

 
 return 0;

}


