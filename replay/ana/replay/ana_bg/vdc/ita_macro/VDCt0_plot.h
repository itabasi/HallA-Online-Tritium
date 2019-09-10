#ifndef VDCt0_plot_h
#define VDCt0_plot_h 1
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
using namespace std;
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
#include "Setting.h"
#include "VDCt0.h"


const  int run_max=100;

class VDCt0_plot : public VDCt0{

 public:
  VDCt0_plot();
  ~VDCt0_plot();

 public:
  void SetPoint(string ifname);
  void SetPointError(string ifname);  
  void MakeRoot(string ofname);
  void Draw();
  void Print(string ofname);


  //=== GetPlot ====//
  TGraph* gLu1[nwire];
  TGraph* gLu2[nwire];
  TGraph* gLv1[nwire];
  TGraph* gLv2[nwire];
  TGraph* gRu1[nwire];
  TGraph* gRu2[nwire];
  TGraph* gRv1[nwire];
  TGraph* gRv2[nwire];

  TGraphAsymmErrors* gLu1_ac[nwire];
  TGraphAsymmErrors* gLu2_ac[nwire];
  TGraphAsymmErrors* gLv1_ac[nwire];
  TGraphAsymmErrors* gLv2_ac[nwire];
  TGraphAsymmErrors* gRu1_ac[nwire];
  TGraphAsymmErrors* gRu2_ac[nwire];
  TGraphAsymmErrors* gRv1_ac[nwire];
  TGraphAsymmErrors* gRv2_ac[nwire];  
  TH1F* hLu1[nwire];
  
  double T0[nwire][8][run_max];
  double t0_min[nwire][8][run_max];
  double t0_max[nwire][8][run_max];
  double dT0[nwire][8][run_max];
  int runnum[run_max];  
  double dT0_l[nwire][8][run_max],dT0_h[nwire][8][run_max];
  
};

VDCt0_plot::VDCt0_plot(){
  set=new Setting();
  set->Initialize();

}

void VDCt0_plot::SetPoint(string ifname){


  string line;
  int j=0;
  int jmax=0;
  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  
  string buf, runname, run;
  string buf_min,buf_max;
  while(    getline(ifp,buf)  ){
    //    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    //    cout<<"param file : "<<runname<<endl;

    //===== Param file ======//


  ifstream ifparam(runname.c_str());  
    
  if (ifparam.fail()){ cerr << "failed open files" <<runname<<endl; break;}

  run=runname.substr(82,6);
  runnum[j]=atoi(run.c_str());

  //  cout<<"run number: "<<runnum[j]<<endl;
  int plane=-1;
  int i=0;
  while( getline(ifparam,buf) ){
    //     getline(ifparam,buf);

    if( buf[0]=='#' ){
      i=0;
      plane++;}
    if( ifparam.eof() || i+1 > nwire){continue;}

    ifparam >> T0[i][plane][j] >> T0[i+1][plane][j] >> T0[i+2][plane][j] >> T0[i+3][plane][j] >> T0[i+4][plane][j] >> T0[i+5][plane][j] >> T0[i+6][plane][j] >> T0[i+7][plane][j];


    
    i=i+8;
  }//end while dat file

  j++;
  jmax=j;
  }//end while open file




  //====== Set Point ====//

    //      cout<<"k : "<<k<<endl;
      //      cout<<"runnum : "<<runnum[k]<<endl;      

  for(int i=0;i<nwire;i++){

      gRu1[i]=new TGraph();
      gRu1[i]->SetMarkerStyle(21);
      gRu1[i]->SetMarkerColor(kRed);
      gRu1[i]->SetMarkerSize(0.5);
      gRu2[i]=new TGraph();
      gRu2[i]->SetMarkerStyle(21);
      gRu2[i]->SetMarkerColor(kRed);
      gRu2[i]->SetMarkerSize(0.5);      
      gRv1[i]=new TGraph();
      gRv1[i]->SetMarkerStyle(21);
      gRv1[i]->SetMarkerColor(kRed);
      gRv1[i]->SetMarkerSize(0.5);
      gRv2[i]=new TGraph();
      gRv2[i]->SetMarkerStyle(21);
      gRv2[i]->SetMarkerColor(kRed);
      gRv2[i]->SetMarkerSize(0.5);      
      gLu1[i]=new TGraph();
      gLu1[i]->SetMarkerStyle(21);
      gLu1[i]->SetMarkerColor(kRed);
      gLu1[i]->SetMarkerSize(0.5);
      gLu2[i]=new TGraph(); 
      gLu2[i]->SetMarkerStyle(21);
      gLu2[i]->SetMarkerColor(kRed);
      gLu2[i]->SetMarkerSize(0.5);     
      gLv1[i]=new TGraph();
      gLv1[i]->SetMarkerStyle(21);
      gLv1[i]->SetMarkerColor(kRed);
      gLv1[i]->SetMarkerSize(0.5);      
      gLv2[i]=new TGraph();
      gLv2[i]->SetMarkerStyle(21);
      gLv2[i]->SetMarkerColor(kRed);
      gLv2[i]->SetMarkerSize(0.5);

      
  for(int k=0;k<jmax;k++){      

      //      cout<<"i: "<<i<<" k : "<<k<<" run : "<<runnum[k]<<" T0 : "<<T0[i][0][k]<<endl;
      gRu1[i]->SetPoint(k,runnum[k]-111000,T0[i][0][k]);
      set->SetGr(gRu1[i],Form("RVDC U1 Wire %d T0 graph ",i),"# Run","Offset [ch (0.5 ns)]");
      gRu2[i]->SetPoint(k,runnum[k]-111000,T0[i][1][k]);
      set->SetGr(gRu2[i],Form("RVDC U2 Wire %d T0 graph ",i),"# Run","Offset [ch (0.5 ns)]");      
      gRv1[i]->SetPoint(k,runnum[k]-111000,T0[i][2][k]);
      set->SetGr(gRv1[i],Form("RVDC V1 Wire %d T0 graph ",i),"# Run","Offset [ch (0.5 ns)]");
      gRv2[i]->SetPoint(k,runnum[k]-111000,T0[i][3][k]);
      set->SetGr(gRv2[i],Form("RVDC V2 Wire %d T0 graph ",i),"# Run","Offset [ch (0.5 ns)]");      
      gLu1[i]->SetPoint(k,runnum[k]-111000,T0[i][4][k]);
      set->SetGr(gLu1[i],Form("LVDC U1 Wire %d T0 graph ",i),"# Run","Offset [ch (0.5 ns)]");
      gLu2[i]->SetPoint(k,runnum[k]-111000,T0[i][5][k]);
      set->SetGr(gLu2[i],Form("LVDC U2 Wire %d T0 graph ",i),"# Run","Offset [ch (0.5 ns)]");      
      gLv1[i]->SetPoint(k,runnum[k]-111000,T0[i][6][k]);
      set->SetGr(gLv1[i],Form("LVDC V1 Wire %d T0 graph ",i),"# Run","Offset [ch (0.5 ns)]");
      gLv2[i]->SetPoint(k,runnum[k]-111000,T0[i][7][k]);      
      set->SetGr(gLv2[i],Form("LVDC V2 Wire %d T0 graph ",i),"# Run","Offset [ch (0.5 ns)]");
      
    }
  }

 
}

////////////////////////////////////////////////////////////


void VDCt0_plot::SetPointError(string ifname){

  cout<<"======================================"<<endl;
  cout<<"========== SetPoint Error ============"<<endl;  
  cout<<"======================================"<<endl;
  string line;
  int j=0;
  int jmax=0;
  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  
  string buf, runname, run;
  string buf_min,buf_max;
  string buf_err;
  while(    getline(ifp,buf)  ){
    //    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    //    cout<<"param file : "<<runname<<endl;

    //===== Param file ======//

    string main = runname.substr(0,95);
    string ifname_err = main + "_err.dat";
    //  string ifname_min  = main + "_min.dat";
    //  string ifname_max  = main + "_max.dat";
  ifstream ifparam(runname.c_str());
  ifstream ifparam_err(Form("%s",ifname_err.c_str()),ios::in);  

  //  ifstream ifparam_min(Form("%s",ifname_min.c_str()),ios::in);  
  //  ifstream ifparam_max(Form("%s",ifname_max.c_str()),ios::in);
  //  cout<<"error file: "<<ifname_err<<endl;
  //  cout<<"min file: "<<ifname_min<<endl;
  //  cout<<"max file: "<<ifname_max<<endl;
  
  if (ifparam.fail()){ cerr << "failed open files " <<runname<<endl; break;}
  if (ifparam_err.fail()){ cerr << "failed open dT0 files " <<ifname_err.c_str()<<endl; break;}  
  //  if (ifparam_min.fail()){ cerr << "failed open T0min files " <<ifname_min.c_str()<<endl; break;}
  //  if (ifparam_max.fail()){ cerr << "failed open T0max files " <<ifname_max.c_str()<<endl; break;}  
  run=runname.substr(82,6);
  runnum[j]=atoi(run.c_str());

  int plane=-1;
  int i=0;


  
  while( getline(ifparam,buf) ){
    if( buf[0]=='#' ){
      i=0;
      plane++;}
    if( ifparam.eof() || i+1 > nwire){continue;}
    ifparam >> T0[i][plane][j] >> T0[i+1][plane][j] >> T0[i+2][plane][j] >> T0[i+3][plane][j] >> T0[i+4][plane][j] >> T0[i+5][plane][j] >> T0[i+6][plane][j] >> T0[i+7][plane][j]; 
    i=i+8;
  }//end while dat file

  
  //=== dT0 ====//
  plane=-1;
  while( getline(ifparam_err,buf_err) ){
    if( buf_err[0]=='#' ){
      i=0;
      plane++;}
      //      cout<<"plane : "<<plane<<" j: "<<j<<" : ifparam: "<<buf_max<<endl;}
    if( ifparam_err.eof() || i+1 > nwire){continue;}
    ifparam_err >> dT0[i][plane][j] >> dT0[i+1][plane][j] >> dT0[i+2][plane][j] >> dT0[i+3][plane][j] >> dT0[i+4][plane][j] >> dT0[i+5][plane][j] >> dT0[i+6][plane][j] >> dT0[i+7][plane][j];

    //    cout<<"i: "<<i<<" "<<dT0[i][plane][j] <<" "<< dT0[i+1][plane][j] <<" "<< dT0[i+2][plane][j] <<" "<< dT0[i+3][plane][j] <<" "<< dT0[i+4][plane][j] <<" "<< dT0[i+5][plane][j] <<" "<< dT0[i+6][plane][j] <<" "<< dT0[i+7][plane][j]<<endl; 

    
    i=i+8;


    
  }//end while dat file



 
  /*
  //==== Max Error ====//
    plane=-1;
  while( getline(ifparam_max,buf_max) ){
    if( buf_max[0]=='#' ){i=0;   plane++;}
      //     cout<<"plane : "<<plane<<" j: "<<j<<" : ifparam: "<<buf_max<<endl;}
    if( ifparam_max.eof() || i+1 > nwire){ continue;}
    ifparam_max >> t0_max[i][plane][j] >> t0_max[i+1][plane][j] >> t0_max[i+2][plane][j] >> t0_max[i+3][plane][j] >> t0_max[i+4][plane][j] >> t0_max[i+5][plane][j] >> t0_max[i+6][plane][j] >> t0_max[i+7][plane][j];


    dT0_h[i][plane][j]=+t0_max[i][plane][j]-T0[i][plane][j];
    dT0_h[i+1][plane][j]=+t0_max[i+1][plane][j]-T0[i+1][plane][j];      
    dT0_h[i+2][plane][j]=+t0_max[i+2][plane][j]-T0[i+2][plane][j];      
    dT0_h[i+3][plane][j]=+t0_max[i+3][plane][j]-T0[i+3][plane][j];      
    dT0_h[i+4][plane][j]=+t0_max[i+4][plane][j]-T0[i+4][plane][j];      
    dT0_h[i+5][plane][j]=+t0_max[i+5][plane][j]-T0[i+5][plane][j];      
    dT0_h[i+6][plane][j]=+t0_max[i+6][plane][j]-T0[i+6][plane][j];
    dT0_h[i+7][plane][j]=+t0_max[i+7][plane][j]-T0[i+7][plane][j];    
    
    i=i+8;
  }//end while dat file




  //===== Min Error =====//
  plane=-1;
  while( getline(ifparam_min,buf_min) ){
    if( buf_min[0]=='#' ){ i=0;   plane++;}//  cout<<"plane : "<<plane<<" j "<<j <<" : buf_min: "<<buf_min<<endl;}
    if( ifparam_min.eof() || i+1 > nwire){  continue;}
    ifparam_min >> t0_min[i][plane][j] >> t0_min[i+1][plane][j] >> t0_min[i+2][plane][j] >> t0_min[i+3][plane][j] >> t0_min[i+4][plane][j] >> t0_min[i+5][plane][j] >> t0_min[i+6][plane][j] >> t0_min[i+7][plane][j];


    dT0_l[i][plane][j]=fabs(-t0_min[i][plane][j]+T0[i][plane][j]);
    dT0_l[i+1][plane][j]=fabs(-t0_min[i+1][plane][j]+T0[i+1][plane][j]);      
    dT0_l[i+2][plane][j]=fabs(-t0_min[i+2][plane][j]+T0[i+2][plane][j]);      
    dT0_l[i+3][plane][j]=fabs(-t0_min[i+3][plane][j]+T0[i+3][plane][j]);      
    dT0_l[i+4][plane][j]=fabs(-t0_min[i+4][plane][j]+T0[i+4][plane][j]);      
    dT0_l[i+5][plane][j]=fabs(-t0_min[i+5][plane][j]+T0[i+5][plane][j]);      
    dT0_l[i+6][plane][j]=fabs(-t0_min[i+6][plane][j]+T0[i+6][plane][j]);
    dT0_l[i+7][plane][j]=fabs(-t0_min[i+7][plane][j]+T0[i+7][plane][j]);    
 
    i=i+8;
  }//end while dat file
  */
  
  
   
  j++;
  jmax=j;
  }//end while open file



  //====== Set Point ====//


  for(int i=0;i<nwire;i++){

      gRu1_ac[i]=new TGraphAsymmErrors();
      gRu1_ac[i]->SetMarkerStyle(21);
      gRu1_ac[i]->SetMarkerColor(kRed);
      gRu1_ac[i]->SetMarkerSize(0.5);
      gRu2_ac[i]=new TGraphAsymmErrors();
      gRu2_ac[i]->SetMarkerStyle(21);
      gRu2_ac[i]->SetMarkerColor(kRed);
      gRu2_ac[i]->SetMarkerSize(0.5);      
      gRv1_ac[i]=new TGraphAsymmErrors();
      gRv1_ac[i]->SetMarkerStyle(21);
      gRv1_ac[i]->SetMarkerColor(kRed);
      gRv1_ac[i]->SetMarkerSize(0.5);
      gRv2_ac[i]=new TGraphAsymmErrors();
      gRv2_ac[i]->SetMarkerStyle(21);
      gRv2_ac[i]->SetMarkerColor(kRed);
      gRv2_ac[i]->SetMarkerSize(0.5);      
      gLu1_ac[i]=new TGraphAsymmErrors();
      gLu1_ac[i]->SetMarkerStyle(21);
      gLu1_ac[i]->SetMarkerColor(kRed);
      gLu1_ac[i]->SetMarkerSize(0.5);
      gLu2_ac[i]=new TGraphAsymmErrors(); 
      gLu2_ac[i]->SetMarkerStyle(21);
      gLu2_ac[i]->SetMarkerColor(kRed);
      gLu2_ac[i]->SetMarkerSize(0.5);     
      gLv1_ac[i]=new TGraphAsymmErrors();
      gLv1_ac[i]->SetMarkerStyle(21);
      gLv1_ac[i]->SetMarkerColor(kRed);
      gLv1_ac[i]->SetMarkerSize(0.5);      
      gLv2_ac[i]=new TGraphAsymmErrors();
      gLv2_ac[i]->SetMarkerStyle(21);
      gLv2_ac[i]->SetMarkerColor(kRed);
      gLv2_ac[i]->SetMarkerSize(0.5);

      
  for(int k=0;k<jmax;k++){      

    
    gRu1_ac[i]->SetPoint(k,runnum[k]-111000,T0[i][0][k]);
    gRu1_ac[i]->SetPointError(k,0,0,dT0[i][0][k],dT0[i][0][k]);            
    //    gRu1_ac[i]->SetPointError(k,0,0,dT0_l[i][0][k],dT0[i][0][k]);            
    set->SetGr(gRu1_ac[i],Form("RVDC U1 Wire %d t0 graph ",i),"# Run","Offset [ch (0.5 ns)]");
    gRu2_ac[i]->SetPoint(k,runnum[k]-111000,T0[i][1][k]);
    //    gRu2_ac[i]->SetPointError(k,0,0,dT0_l[i][1][k],dT0[i][1][k]);
    gRu2_ac[i]->SetPointError(k,0,0,dT0[i][1][k],dT0[i][1][k]);
    set->SetGr(gRu2_ac[i],Form("RVDC U2 Wire %d t0 graph ",i),"# Run","Offset [ch (0.5 ns)]");   
    gRv1_ac[i]->SetPoint(k,runnum[k]-111000,T0[i][2][k]);
    gRv1_ac[i]->SetPointError(k,0,0,dT0[i][2][k],dT0[i][2][k]);
    //    gRv1_ac[i]->SetPointError(k,0,0,dT0_l[i][2][k],dT0[i][2][k]);                
    set->SetGr(gRv1_ac[i],Form("RVDC V1 Wire %d t0 graph ",i),"# Run","Offset [ch (0.5 ns)]");
    gRv2_ac[i]->SetPoint(k,runnum[k]-111000,T0[i][3][k]);
    gRv2_ac[i]->SetPointError(k,0,0,dT0[i][3][k],dT0[i][3][k]);
    //    gRv2_ac[i]->SetPointError(k,0,0,dT0_l[i][3][k],dT0[i][3][k]);                      
    set->SetGr(gRv2_ac[i],Form("RVDC V2 Wire %d t0 graph ",i),"# Run","Offset [ch (0.5 ns)]");   
    gLu1_ac[i]->SetPoint(k,runnum[k]-111000,T0[i][4][k]);
    gLu1_ac[i]->SetPointError(k,0,0,dT0[i][4][k],dT0[i][4][k]);
    //    gLu2_ac[i]->SetPointError(k,0,0,dT0_l[i][5][k],dT0[i][5][k]);                  
    set->SetGr(gLu1_ac[i],Form("LVDC U1 Wire %d t0 graph ",i),"# Run","Offset [ch (0.5 ns)]");
    gLu2_ac[i]->SetPoint(k,runnum[k]-111000,T0[i][5][k]);
    gLu2_ac[i]->SetPointError(k,0,0,dT0[i][5][k],dT0[i][5][k]);
    //    gLu2_ac[i]->SetPointError(k,0,0,dT0_l[i][5][k],dT0[i][5][k]);                      
    set->SetGr(gLu2_ac[i],Form("LVDC U2 Wire %d t0 graph ",i),"# Run","Offset [ch (0.5 ns)]");   
    gLv1_ac[i]->SetPoint(k,runnum[k]-111000,T0[i][6][k]);
    gLv1_ac[i]->SetPointError(k,0,0,dT0[i][6][k],dT0[i][6][k]);
    //    gLu2_ac[i]->SetPointError(k,0,0,dT0_l[i][5][k],dT0[i][5][k]);                  
    set->SetGr(gLv1_ac[i],Form("LVDC V1 Wire %d t0 graph ",i),"# Run","Offset [ch (0.5 ns)]");
    gLv2_ac[i]->SetPoint(k,runnum[k]-111000,T0[i][7][k]);
    gLv2_ac[i]->SetPointError(k,0,0,dT0[i][7][k],dT0[i][7][k]);
    //    gLv2_ac[i]->SetPointError(k,0,0,dT0_l[i][7][k],dT0[i][7][k]);                
    set->SetGr(gLv2_ac[i],Form("LVDC V2 Wire %d t0 graph ",i),"# Run","Offset [ch (0.5 ns)]");
      
    }
  }

 
}


////////////////////////////////////////////////////////////

void VDCt0_plot::MakeRoot(string ofname){

  fnew =new  TFile(ofname.c_str(),"recreate");

  for(int i=0;i<nwire;i++){
    /*
    gRu1[i]->SetName(Form("gRu1_%d",i));
    gRu1[i]->Write();
    gRu2[i]->SetName(Form("gRu2_%d",i));
    gRu2[i]->Write();
    gRv1[i]->SetName(Form("gRv1_%d",i));
    gRv1[i]->Write();
    gRv2[i]->SetName(Form("gRv2_%d",i));
    gRv2[i]->Write();
    gLu1[i]->SetName(Form("gLu1_%d",i));
    gLu1[i]->Write();
    gLu2[i]->SetName(Form("gLu2_%d",i));
    gLu2[i]->Write();
    gLv1[i]->SetName(Form("gLv1_%d",i));
    gLv1[i]->Write();
    gLv2[i]->SetName(Form("gLv2_%d",i));
    gLv2[i]->Write();
    */
    
    gRu1_ac[i]->SetName(Form("gRu1_ac_%d",i));
    gRu1_ac[i]->Write();
    gRu2_ac[i]->SetName(Form("gRu2_ac_%d",i));
    gRu2_ac[i]->Write();
    gRv1_ac[i]->SetName(Form("gRv1_ac_%d",i));
    gRv1_ac[i]->Write();
    gRv2_ac[i]->SetName(Form("gRv2_ac_%d",i));
    gRv2_ac[i]->Write();
    gLu1_ac[i]->SetName(Form("gLu1_ac_%d",i));
    gLu1_ac[i]->Write();
    gLu2_ac[i]->SetName(Form("gLu2_ac_%d",i));
    gLu2_ac[i]->Write();
    gLv1_ac[i]->SetName(Form("gLv1_ac_%d",i));
    gLv1_ac[i]->Write();
    gLv2_ac[i]->SetName(Form("gLv2_ac_%d",i));
    gLv2_ac[i]->Write();        
    
    
  }

  fnew->Close();
  
}

///////////////////////////////////////////////////////////

void VDCt0_plot::Draw(){

  /*
  for(int j=0;j<11;j++){
 
    c0[j] =new TCanvas(Form("c0_%d",j),Form("LVDC-U1_%d",j));
    c0[j]->Divide(4,8);
    c1[j] =new TCanvas(Form("c1_%d",j),Form("LVDC-U2_%d",j));
    c1[j]->Divide(4,8);
    c2[j] =new TCanvas(Form("c2_%d",j),Form("LVDC-V1_%d",j));
    c2[j]->Divide(4,8);
    c3[j] =new TCanvas(Form("c3_%d",j),Form("LVDC-V2_%d",j));
    c3[j]->Divide(4,8);

   for(int i=0;i<32;i++){
    c0[j]->cd(i+1);
      gLu1[32*j+i]->SetMarkerStyle(21);
      gLu1[32*j+i]->SetMarkerColor(2);
      gLu1[32*j+i]->SetMarkerSize(0.1);
    gLu1[32*j+i]->Draw("ALP");
    gLu1[32*j+i]->GetXaxis()->SetRangeUser(160,500);
    gLu1[32*j+i]->GetYaxis()->SetRangeUser(2800,3000);    
    c1[j]->cd(i+1);
      gLu2[32*j+i]->SetMarkerStyle(21);
      gLu2[32*j+i]->SetMarkerColor(2);
      gLu2[32*j+i]->SetMarkerSize(0.1);    
    gLu2[32*j+i]->Draw("ALP");
    gLu2[32*j+i]->GetXaxis()->SetRangeUser(160,500);    
    gLu2[32*j+i]->GetYaxis()->SetRangeUser(2800,3000);    
    c2[j]->cd(i+1);
      gLv1[32*j+i]->SetMarkerStyle(21);
      gLv1[32*j+i]->SetMarkerColor(2);
      gLv1[32*j+i]->SetMarkerSize(0.1);
    gLv1[32*j+i]->GetXaxis()->SetRangeUser(160,500);    
    gLv1[32*j+i]->GetYaxis()->SetRangeUser(2800,3000);    
    gLv1[32*j+i]->Draw("ALP");
    c3[j]->cd(i+1);
      gLv2[32*j+i]->SetMarkerStyle(21);
      gLv2[32*j+i]->SetMarkerColor(2);
      gLv2[32*j+i]->SetMarkerSize(0.1);    
    gLv2[32*j+i]->GetXaxis()->SetRangeUser(160,500);    
    gLv2[32*j+i]->GetYaxis()->SetRangeUser(2800,3000);    
    gLv2[32*j+i]->Draw("ALP"); 
    }

}

  */
  
  //======TGraphError ========//


  for(int j=0;j<11;j++){
 
    c4[j] =new TCanvas(Form("c4_%d",j),Form("LVDC-U1_%d",j));
    c4[j]->Divide(4,8);
    c5[j] =new TCanvas(Form("c5_%d",j),Form("LVDC-U2_%d",j));
    c5[j]->Divide(4,8);
    c6[j] =new TCanvas(Form("c6_%d",j),Form("LVDC-V1_%d",j));
    c6[j]->Divide(4,8);
    c7[j] =new TCanvas(Form("c7_%d",j),Form("LVDC-V2_%d",j));
    c7[j]->Divide(4,8);

   for(int i=0;i<32;i++){
    c4[j]->cd(i+1);
      gLu1_ac[32*j+i]->SetMarkerStyle(21);
      gLu1_ac[32*j+i]->SetMarkerColor(2);
      gLu1_ac[32*j+i]->SetMarkerSize(0.5);
    gLu1_ac[32*j+i]->Draw("AP");
    gLu1_ac[32*j+i]->GetXaxis()->SetRangeUser(160,500);
    gLu1_ac[32*j+i]->GetYaxis()->SetRangeUser(2800,3000);    
    c5[j]->cd(i+1);
      gLu2_ac[32*j+i]->SetMarkerStyle(21);
      gLu2_ac[32*j+i]->SetMarkerColor(2);
      gLu2_ac[32*j+i]->SetMarkerSize(0.5);    
    gLu2_ac[32*j+i]->Draw("AP");
    gLu2_ac[32*j+i]->GetXaxis()->SetRangeUser(160,500);    
    gLu2_ac[32*j+i]->GetYaxis()->SetRangeUser(2800,3000);    
    c6[j]->cd(i+1);
      gLv1_ac[32*j+i]->SetMarkerStyle(21);
      gLv1_ac[32*j+i]->SetMarkerColor(2);
      gLv1_ac[32*j+i]->SetMarkerSize(0.5);
    gLv1_ac[32*j+i]->GetXaxis()->SetRangeUser(160,500);    
    gLv1_ac[32*j+i]->GetYaxis()->SetRangeUser(2800,3000);    
    gLv1_ac[32*j+i]->Draw("AP");
    c7[j]->cd(i+1);
      gLv2_ac[32*j+i]->SetMarkerStyle(21);
      gLv2_ac[32*j+i]->SetMarkerColor(2);
      gLv2_ac[32*j+i]->SetMarkerSize(0.5);    
    gLv2_ac[32*j+i]->GetXaxis()->SetRangeUser(160,500);    
    gLv2_ac[32*j+i]->GetYaxis()->SetRangeUser(2800,3000);    
    gLv2_ac[32*j+i]->Draw("AP"); 
    }

}
  


  

}

///////////////////////////////////////////////////////////

void VDCt0_plot::Print(string ofname){


  cout<<"Print is starting "<<endl;
  cout<<"pdf name : "<<ofname<<endl;

  /*
  for(int j=0;j<4;j++){
  for(int i=0;i<11;i++){
    if(i==0 && j==0)c0[i]->Print(Form("%s[",ofname.c_str()));
    if(j==0)c0[i]->Print(Form("%s",ofname.c_str()));
    if(j==1)c1[i]->Print(Form("%s",ofname.c_str()));
    if(j==2)c2[i]->Print(Form("%s",ofname.c_str()));
    if(j==3)c3[i]->Print(Form("%s",ofname.c_str()));    
    if(i==10 && j==3)c3[i]->Print(Form("%s]",ofname.c_str()));
  }
  }
  */

  for(int j=0;j<4;j++){
  for(int i=0;i<11;i++){
    if(i==0 && j==0)c4[i]->Print(Form("%s[",ofname.c_str()));
    if(j==0)c4[i]->Print(Form("%s",ofname.c_str()));
    if(j==1)c5[i]->Print(Form("%s",ofname.c_str()));
    if(j==2)c6[i]->Print(Form("%s",ofname.c_str()));
    if(j==3)c7[i]->Print(Form("%s",ofname.c_str()));    
    if(i==10 && j==3)c7[i]->Print(Form("%s]",ofname.c_str()));
  }
  }
  
  cout<<"Print is done !"<<endl;


}

////////////////////////////////////////////////////////////

#endif
