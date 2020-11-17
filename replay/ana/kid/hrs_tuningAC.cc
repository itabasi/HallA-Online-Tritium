#include <iostream>
#include <fstream>
using namespace std;
#include "TApplication.h"
#include "hrs_tuningAC.h"


//=======================================================//
//================     Main       =======================//
//=======================================================//


int main(int argc, char** argv){

  gStyle->SetOptFit(111111111);
  int ch;
  string ifname = "/adaqfs/home/a-onl/tritium_work/itabashi/nnL/HallA-Online-Tritium/replay/scripts/ita_scripts/run_list/Lambda_test.list";
  string ofname = "/pdf/hydro1_AC_eff_test.pdf";
  string root_name;
  string print_name;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag  = true;
  bool coin_flag  = false;
  bool print_flag = false;
  bool root_flag  = false;
  //  bool ac2_min=true;
  string itel;  
  string pngname;
  //  bool single_flag=false;

  
  extern char *optarg;
  while((ch=getopt(argc,argv,"h:f:w:s:n:r:i:o:bcop:GHT12"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 's':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      single_flag=true;
      root_flag = true;
      draw_flag = false;            
      break;

      
    case 'w':
      print_flag = true;
       draw_flag = false;
      print_name = optarg;
      cout<<"output PDF filename : "<<print_name<<endl;
      break;


    case 'r':
      root_flag = true;
      draw_flag = false;      
      root_name = optarg;
      cout<<"output root filename : "<<root_name<<endl;      
      break;

    case 'o':
      root_flag=true;
      draw_flag=false;
      print_flag=true;
      ofname = optarg;
      root_name="./../rootfiles/ACtuning/" + ofname + ".root";
      print_name="./../pdf/ACtuning/" +ofname + ".pdf";
      break;
      
    case 'i':
      itel= optarg;
      nth= atoi(itel.c_str());
      break;

      
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
  
    case 'G':
    mode="G";
      break;
  
    case 'H':
    mode="H";
      break;

    case 'T':
      mode="T";    
	break;

    case 'U':
      ac2_min=false;    
	break;

  case '1':
    tdc_time=56.23e-3;//[ns]
    kine=1;
      break;

  case '2':
    tdc_time=58e-3;//[ns]
    kine=2;
      break;


    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
      cout<<"-r : output root filename"<<endl;      
      cout<<"-o : output pdf & root  filename"<<endl;
      cout<<"-1 : F1TDC resolution 56 ns"<<endl;
      cout<<"-2 : F1TDC resolution 58 ns"<<endl;      
      cout<<"-T or H or G : Mode of root"<<endl;      
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


  
  tuningAC* AC=new tuningAC();
  if(single_flag)AC->SetRun(ifname);
  else AC->SetRunList(ifname);
  if(root_flag)AC->SetRoot(root_name);
  AC->SetBranch();
  AC->SetParam();
  AC->MakeHist();
  AC->Fill();
  AC->Fitting();
  AC->Tuning();
  AC->Draw();
  if(print_flag)AC->Print(print_name);
  if(root_flag)AC->Write();
  AC->Comment();
  
  cout<<"=========== Output files ============="<<endl;
  cout<<"output rootfile: "<<root_name<<endl;
  
  
	
  TApplication *theApp =new TApplication("App",&argc,argv);
 if(draw_flag==0)gROOT->SetBatch(1);
 if(draw_flag==0)gSystem->Exit(1);
 theApp->Run();
 return 0;

}
