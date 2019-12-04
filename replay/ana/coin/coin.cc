/*
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
*/

void coin(){

  //  string ifname="";

  /*
  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> runname;
    T->Add(runname.c_str());

  }
  */

  string ifname="../rootfiles/mmass/ana_Lambda/Lambda_small_optH_1119.root";  
  TChain* T=new TChain("T");
  T->Add(ifname.c_str());
  int ENum =T->GetEntries();
  cout<<"Get Entries: "<<ENum<<endl;

  TCanvas* c0=new TCanvas("c0","c0");

  double ct,Rx_fp,ac1,ac2,trig;
  int z_cut,runnum,pid_cut;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("ct_b",1);
  T->SetBranchAddress("ct_b",&ct);
  T->SetBranchStatus("runnum",1);
  T->SetBranchAddress("runnum",&runnum);  
  T->SetBranchStatus("Rx_fp",1);
  T->SetBranchAddress("Rx_fp",&Rx_fp);  
  T->SetBranchStatus("Rx_fp",1);
  T->SetBranchAddress("Rx_fp",&Rx_fp);
  T->SetBranchStatus("trig",1);
  T->SetBranchAddress("trig",&trig);
  T->SetBranchStatus("z_cut",1);
  T->SetBranchAddress("z_cut",&z_cut);
  T->SetBranchStatus("pid_cut",1);
  T->SetBranchAddress("pid_cut",&pid_cut);      
  T->SetBranchStatus("ac1_npe_sum",1);
  T->SetBranchAddress("ac1_npe_sum",&ac1);
  T->SetBranchStatus("ac2_npe_sum",1);
  T->SetBranchAddress("ac2_npe_sum",&ac2);


  int bin_run,min_run,max_run;

  double min_ct=-5.0;
  double max_ct=5.0;

  int bin_ct=(int)(max_ct-min_ct)*50;
  min_run=111160;
  max_run=111220;
  bin_run=max_run-min_run;
  
  TH2D*hcoin_run=new TH2D("hcoin_run","Coin time Run dependece ; Runnum ; coin time [ns]",bin_run,min_run,max_run,bin_ct,min_ct,max_ct);
  TH1D* hcoin=new TH1D("hcoin","Coin time ; coin time [ns] ;Counts",bin_ct,min_ct,max_ct);
  TH1D* hcoin_c=new TH1D("hcoin_c","Coin time ; coin time [ns] ;Counts",bin_ct,min_ct,max_ct);  

  bool pi_cut=false;
  double off[bin_run];
  int nrun[bin_run];
  string pname="param/test.param";
  ifstream ifp(pname.c_str());
  int s=0;
  string buf;
  string p;
  string nruns,offs;
  

  while(1){
    if( ifp.eof() ) break;
    ifp >> nrun[s] >> off[s] ;
    //    cout<<"s "<<s<<" run "<<nrun[s]<<" off "<<off[s]<<endl;
    s++;
  }  

  //test
  ENum=10000000;

  
  for(int k=0; k<ENum;k++){

    ct=-1000.;
    Rx_fp=-10.;
    z_cut=-1;
    pid_cut=-1;
    trig=0;
    pi_cut=false;
    T->GetEntry(k);
    if(ac1>1.0 && ac2>5.0)pi_cut=true;
    if(pi_cut && z_cut>0 && trig==5)hcoin_run->Fill(runnum,ct);
    
    if(pid_cut && z_cut>0 && trig==5){

      //      hcoin->Fill(ct);
      int nn=0;

      while(nn<bin_run){
	if(runnum==nrun[nn]){
	  hcoin->Fill(ct);
	  hcoin_c->Fill(ct-off[nn]);
	  break;
	}else{ nn++;}
	  //if(nn<bin_run)nn++;
      } // while
      
    }
  }//end fill
  

  ///////////////////////
  // Draw
  ///////////////////////

  hcoin->SetLineColor(4);
  hcoin->SetFillColor(4);
  hcoin->SetFillStyle(3002);
  hcoin_c->SetLineColor(2);
  hcoin_c->SetFillColor(2);
  hcoin_c->SetFillStyle(3002);
  
  c0->cd();
  hcoin_c->Draw();
  hcoin->Draw("same");
}

