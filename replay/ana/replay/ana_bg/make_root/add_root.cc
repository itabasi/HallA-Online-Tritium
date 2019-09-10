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
//#include "Setting.h"

using namespace std;
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  int ch;
  extern char *optarg;
  string ifname("input.root");
  string ofname("output.root");

  while((ch=getopt(argc,argv,"hf:wr:"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      break;
    case 'w':
      ofname = optarg;
      break;
    case 'r':
      ofname = optarg;
      break;
    case 'h':
      std::cout<<"-f (inputfile): input ROOT file name"<<std::endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;

    }
  }


  // TApplication *theApp=new Tapplication("App", &argc,argv);
   TApplication theApp("App", &argc, argv);

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf,runname;
  TChain *oldtree = new TChain("T");
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> runname;


 //oldtree->Add(Form("%s",ifname.c_str()));
 //  oldtree->Add(Form("/adaqfs/home/a-onl/tritium/replay/t2root/Rootfiles/%s",runname.c_str()));
  oldtree->Add(Form(runname.c_str()));
  cout<<buf<<endl;
  }

  cout<<"start ADD Tree "<<endl;
  
  TFile *fout=new TFile(ofname.c_str(),"recreate");
  TTree* newtree=oldtree->CloneTree();
  cout<<"clone is done "<<endl;

 newtree->Write();
 cout<<"witten new tree  "<<endl;
 fout->Close();
 gSystem->Exit(1);
 theApp.Run();
 
 return 0;

}

