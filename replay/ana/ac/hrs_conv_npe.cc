const double c=299792458e-9;// [m/ns]
const double mk=493.7e-3;// Kaon mass [GeV/c^2]
const double me=0.511e-3;// electron mass [GeV/c^2] 
const double ml=1115.7e-3;//Lambda mass [GeV/c^2]
const double mn=939.6e-3; // neutron mass [GeV/c^2]
const double mpi=139.6e-3;// pion mass [GeV/c^2]
const  int nth=3; //th num
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




//===================================================================//
//============================= Main ================================//
//===================================================================//

int main(int argc, char** argv){
  gStyle->SetOptLogy(1);
  int ch; char* mode="H";
  string ifname ="/data/opt_small/pre_run/tritium_93403.root";
  string ofname = "..//pdf/ac/hydro1_AC_eff_test.pdf";
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  bool param_flag =false;
  bool single_flag=true;
  bool root_flag=false;
  string pngname;
  string paramname="./param/offset_ac.dat";
  string ofroot;
  extern char *optarg;
  while((ch=getopt(argc,argv,"h:f:r:w:n:p:s:bcop:GHT"))!=-1){
    switch(ch){

    case 'f':
      ifname = optarg;
      single_flag=false;
      cout<<"input filename : "<<ifname<<endl;
      break;

    case 's':
      ifname = optarg;
      single_flag=true;
      cout<<"input filename : "<<ifname<<endl;
      break;

    case 'r':
      ofroot = optarg;
      root_flag=true;
      cout<<"output root filename : "<<ofroot<<endl;
      break;      

    case 'w':
      output_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;

    case 'p':
     paramname = optarg;
      cout<<"Param filename : "<<paramname<<endl;
      param_flag =true;      
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



  //=============== ROOT File Mode ================//
  /*if(ifname.c_str()=="../run_list/coin_H2_1.root")mode="G";
  else if(ifname.c_str()=="../run_list/Lambda_small.list" ||ifname.c_str()=="../run_list/Lambda_test.list")mode="H_1";
  else{cout<<"false to read mode types ";};
  */


 double tdc_time=56.23e-3;//[ns]
 // double tdc_time=58.e-3;//[ns] 

/*
 if(mode=="G"){tdc_time=56.23e-3;}
 else if(mode=="H"){tdc_time=56.23e-3;}
 else if(mode=="T"){tdc_time=58e-3;}
 */


 //TChain //
  TChain* T;
  if(mode=="G"){T=new TChain("tree"); }
  else {T=new TChain("T");}
  if(single_flag==0){
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
    } 
  }else{
    T->Add(ifname.c_str()); //single root
  }
  
  double pe1_a1[24],pe1_a2[26],n_a1[24],n_a2[26],ped_a1[24],ped_a2[26],conv_a1[24],conv_a2[26];
  int dec=0;
  int seg =-1;
  double pe1=0.0, ped=0.0;
  if(param_flag){
  ifstream ifpa(Form("%s",paramname.c_str()),ios::in);

  if(!ifpa){ cout<<"no input file "<<paramname<<endl; exit(1); }
  string buf2;
    while(1){
      getline(ifpa,buf2);
      if( buf2[0]=='#' ){ continue; }
      if( ifpa.eof() ) break;
      dec=0;
      seg=-1;
      pe1=0.0;
      ped=0.0;
      stringstream ssbuf(buf2);
      ssbuf >> dec >> seg >> ped >> pe1;
      //      cout<<"AC "<<dec<<" "<<seg<<" "<<" "<<pe1-ped<<endl;
      
      if(dec==1 && seg>=0){
	pe1_a1[seg]=pe1;
	ped_a1[seg]=ped;
	conv_a1[seg]=(pe1_a1[seg]-ped_a1[seg]);
      }else if(dec==2 && seg>=0){
	pe1_a2[seg]=pe1;
	ped_a2[seg]=ped;
	conv_a2[seg]=(pe1_a2[seg]-ped_a2[seg]);	
      }else if(seg<0){}
    }
  }

  



  
  cout<<"mode :"<<mode<<endl;
  cout<<"tdc_time[ns]: "<<tdc_time<<endl;
  int evnt=T->GetEntries();
  cout<<"Get Entries: "<<evnt<<endl;


  double Ra1a_p[100],Ra2a_p[100],Ra1sum,Ra2sum,Ra1a[100],Ra2a[100],trig,Ra1sum_p,Ra2sum_p;
  double Ra1a_c[100],Ra2a_c[100];
 T->SetBranchStatus("*",0);
 T->SetBranchStatus("R.a1.a",1);
 T->SetBranchAddress("R.a1.a",Ra1a);
 T->SetBranchStatus("R.a2.a",1);
 T->SetBranchAddress("R.a2.a",Ra2a); 
 T->SetBranchStatus("R.a1.a_p",1);
 T->SetBranchAddress("R.a1.a_p",Ra1a_p);
 T->SetBranchStatus("R.a2.a_p",1);
 T->SetBranchAddress("R.a2.a_p",Ra2a_p);
 T->SetBranchStatus("R.a1.a_c",1);
 T->SetBranchAddress("R.a1.a_c",Ra1a_c);
 T->SetBranchStatus("R.a2.a_c",1);
 T->SetBranchAddress("R.a2.a_c",Ra2a_c); 
 
 T->SetBranchStatus("R.a1.asum_c",1);
 T->SetBranchAddress("R.a1.asum_c",&Ra1sum);
 T->SetBranchStatus("R.a2.asum_c",1);
 T->SetBranchAddress("R.a2.asum_c",&Ra2sum); 
 T->SetBranchStatus("R.a2.asum_p",1);
 T->SetBranchAddress("R.a2.asum_p",&Ra2sum_p); 
 T->SetBranchStatus("R.a2.asum_p",1);
 T->SetBranchAddress("R.a2.asum_p",&Ra2sum_p);
 // T->SetBranchStatus("DR.evtypebits",1);
 // T->SetBranchAddress("DR.evtypebits",&trig);
 
 TH1F *ha1_adc[24];
 TH1F *ha2_adc[26];
 TH1F *ha1_npe[24];
 TH1F *ha2_npe[26];
 TH1F *ha1_adc_scaled[24];
 TH1F *ha2_adc_scaled[26];
 TF1 *fconv_a1[24];
 TF1 *fconv_a2[26];
 
 double bin_ac1,min_ac1,max_ac1,bin_ac2,min_ac2,max_ac2;
 min_ac1=-100.0;
 max_ac1=10000.;
 min_ac2=-100;
 max_ac2=10000.;
 bin_ac1=(max_ac1-min_ac1);
 bin_ac1=(int)bin_ac1;
 bin_ac2=(max_ac2-min_ac2);
 bin_ac2=(int)bin_ac2;

 double bin_npe1,min_npe1,max_npe1,bin_npe2,min_npe2,max_npe2;
 min_npe1=-1.0;
 max_npe1=30.0;
 bin_npe1=(max_npe1-min_npe1)*200;
 bin_npe1=(int)bin_npe1;
 min_npe2=-1.0;
 max_npe2=30.0;
 bin_npe2=(max_npe2-min_npe2)*200;
 bin_npe2=(int)bin_npe2;

 TH1F *ha1_adc_sum=new TH1F("ha1_adc_sum","AC1 ADC SUM HIST",bin_ac1,min_ac1,max_ac1);
 TH1F *ha2_adc_sum=new TH1F("ha2_adc_sum","AC2 ADC SUM HIST",bin_ac2,min_ac2,max_ac2);
 TH1F *ha1_npe_sum=new TH1F("ha1_npe_sum","AC1 NPE SUM HIST",bin_npe1,min_npe1,max_npe1);
 TH1F *ha2_npe_sum=new TH1F("ha2_npe_sum","AC2 NPE SUM HIST",bin_npe2,min_npe2,max_npe2);
 
 for(int i=0;i<24;i++){
   ha1_adc[i]=new TH1F(Form("ha1_adc[%d]",i),"AC1 ADC HIST",bin_ac1,min_ac1,max_ac1);
   ha1_adc_scaled[i]=new TH1F(Form("ha1_adc_scaled[%d]",i),"AC1 ADC HIST",bin_npe1,min_npe1,max_npe1);   
   ha1_npe[i]=new TH1F(Form("ha1_npe[%d]",i),"AC1 NPE HIST",bin_npe1,min_npe1,max_npe1);

   fconv_a1[i]=new TF1(Form("fconv_a1[%d]",i),"gaus",min_ac1,max_ac1);
 }

for(int i=0;i<26;i++){
   ha2_adc[i]=new TH1F(Form("ha2_adc[%d]",i),"AC2 ADC HIST",bin_ac2,min_ac2,max_ac2);
   ha2_adc_scaled[i]=new TH1F(Form("ha2_adc_scaled[%d]",i),"AC2 ADC SCALED HIST",bin_npe2,min_npe2,max_npe2);   
   ha2_npe[i]=new TH1F(Form("ha2_npe[%d]",i),"AC2 NPE HIST",bin_npe2,min_npe2,max_npe2);   
   fconv_a2[i]=new TF1(Form("fconv_a2[%d]",i),"gaus",min_ac2,max_ac2);
 }
 


  //======== New Root ==========//

 TFile* fnew=new TFile(ofroot.c_str(),"recreate");
 TTree* tnew=new TTree("T",ofroot.c_str());
  tnew = T->CloneTree(0);
  double ac1_p[24],ac1_c[24],ac2_p[26],ac2_c[26],ac1_sum_c,ac2_sum_c,ac2_sum_p,ac1_sum_p;
  tnew->Branch("ac1_p",ac1_p,"ac1_p[24]/D");
  tnew->Branch("ac1_c",ac1_c,"ac1_c[24]/D");   
  tnew->Branch("ac2_p",ac2_p,"ac2_p[26]/D");
  tnew->Branch("ac2_c",ac2_c,"ac2_c[26]/D");
  tnew->Branch("ac1_sum_p",&ac1_sum_p,"ac1_sum_p/D");   
  tnew->Branch("ac2_sum_p",&ac2_sum_p,"ac2_sum_p/D");     
  tnew->Branch("ac1_sum_c",&ac1_sum_c,"ac1_sum_c/D");   
  tnew->Branch("ac2_sum_c",&ac2_sum_c,"ac2_sum_c/D");   
 
 
 for(int k=0;k<evnt;k++){
   T->GetEntry(k);
   ac1_sum_c=0.0,ac2_sum_c=0.0;
   ac1_sum_p=0.0,ac2_sum_p=0.0;     

   for(int i=0;i<24;i++){
     ac1_p[i]=Ra1a_p[i];
     ac1_c[i]=Ra1a_p[i]/conv_a1[i];     
     ha1_adc[i]->Fill(Ra1a_p[i]);
     ha1_adc_scaled[i]->Fill(Ra1a_c[i]/400);     
     ha1_npe[i]->Fill(Ra1a_p[i]/conv_a1[i]);
     ha1_npe_sum->Fill(Ra1a_p[i]/conv_a1[i]);
     ac1_sum_p+=ac1_p[i];     
     ac1_sum_c+=ac1_c[i];
   }
   for(int i=0;i<26;i++){
     ac2_p[i]=Ra2a_p[i];
     ac2_c[i]=Ra2a_p[i]/conv_a2[i];          
     ha2_adc[i]->Fill(Ra2a_p[i]);
     ha2_adc_scaled[i]->Fill(Ra2a_c[i]/400.);     
     ha2_npe[i]->Fill(Ra2a_p[i]/conv_a2[i]);
     ha2_npe_sum->Fill(Ra2a_p[i]/conv_a2[i]);
     ac2_sum_p+=ac2_p[i];     
     ac2_sum_c+=ac2_c[i];}
     //     ha2_npe[i]->Fill(Ra2a_p[i]/pe1_a2[i]);
     //     ha2_npe_sum->Fill(Ra2a_p[i]/pe1_a2[i]);}

   ha1_adc_sum->Fill(Ra1sum);
   ha2_adc_sum->Fill(Ra2sum);
   tnew->Fill();
   if(k%100000==0)cout<<"Filled "<<k<<" / "<<evnt<<endl;
}



 /*
 
 for(int i=0;i<24;i++){
   // pe1_a1[i]=ha1_adc[i]->GetBinContent(ha1_adc[i]->GetMaximumBin());
 //n_a1[i]=ha1_adc[i]->GetYaxis()->GetBinContent(ha1_adc[i]->GetMaximumBin());
 // fconv_a1->SetParameter(0,n_a1[i]);
 //fconv_a1[i]->SetParameter(1,pe1_a1[i]);
   fconv_a1[i]->SetParameter(1,400.);
   ha1_adc[i]->Fit(Form("fconv_a1[%d]",i),"Rq","",200,500.);
   pe1_a1[i]=fconv_a1[i]->GetParameter(1);
 }


 for(int i=0;i<26;i++){
   // pe1_a2[i]=ha2_adc[i]->GetBinContent(ha2_adc[i]->GetMaximumBin());
 //n_a2[i]=ha2_adc[i]->GetYaxis()->GetBinContent(ha2_adc[i]->GetMaximumBin());
 //fconv_a2->SetParameter(0,n_a2[i]);
   // fconv_a2[i]->SetParameter(1,pe1_a2[i]);
 fconv_a2[i]->SetParameter(1,400.);
 ha2_adc[i]->Fit(Form("fconv_a2[%d]",i),"Rq","",200.,500.);
 pe1_a2[i]=fconv_a2[i]->GetParameter(1);
 }

 */
 
 /*
 for(int k=0;k<evnt;k++){
   T->GetEntry(k);
   for(int i=0;i<24;i++){
     ha1_npe[i]->Fill(Ra1a_p[i]/pe1_a1[i]);
     ha1_npe_sum->Fill(Ra1a_p[i]/pe1_a1[i]);}
   //  ha1_npe[i]->Fill(Ra1a_p[i]/400.);
   //     ha1_npe_sum->Fill(Ra1a_p[i]/400.);}
   for(int i=0;i<26;i++){
     ha2_npe[i]->Fill(Ra2a_p[i]/pe1_a2[i]);
     ha2_npe_sum->Fill(Ra2a_p[i]/pe1_a2[i]);}
   //     ha2_npe[i]->Fill(Ra2a_p[i]/400.);
   //     ha2_npe_sum->Fill(Ra2a_p[i]/400.);}
 }

 */
 

 TCanvas* c0=new TCanvas("c0","AC1 ADC HIST");
 c0->Divide(4,6);
 TCanvas* c1=new TCanvas("c1","AC1 NPE HIST");
 c1->Divide(4,6);
 TCanvas* c2=new TCanvas("c2","AC2 ADC HIST");
 c2->Divide(4,7);
 TCanvas* c3=new TCanvas("c3","AC2 NPE HIST");
 c3->Divide(4,7);
 TCanvas* c4=new TCanvas("c4","AC SUM HIST");
 c4->Divide(2,2);

 
 for(int i=0;i<24;i++){
   c0->cd(i+1);
   ha1_adc[i]->Draw();
   //   fconv_a1[i]->Draw("same");
   c1->cd(i+1);
   ha1_npe[i]->Draw();
 }
 
 for(int i=0;i<26;i++){
   c2->cd(i+1);
   ha2_adc[i]->Draw();
   //   fconv_a2[i]->Draw("same");
   c3->cd(i+1);
   ha2_npe[i]->Draw();
 }

   c4->cd(1);
   ha1_adc_sum->Draw();
   c4->cd(2);
   ha2_adc_sum->Draw();
   c4->cd(3);
   ha1_npe_sum->Draw();
   c4->cd(4);
   ha2_npe_sum->Draw();



   
TString name;
 if(output_flag){
 name.Form(ofname.c_str());
 c0->Print(name+"[","pdf");
 c0->Print(name,"pdf");
 c1->Print(name,"pdf");
 c2->Print(name,"pdf");
 c3->Print(name,"pdf");
 c4->Print(name,"pdf");
 c4->Print(name+"]","pdf");

 cout<<"Print is done "<<endl;
 }     
    


  //======== New Root ==========//
 if(root_flag){

  for(int i=0;i<24;i++){ 
    ha1_adc[i]->SetName(Form("ha1_adc_%d",i));                   
    ha1_adc[i]->Write();
    ha1_adc_scaled[i]->SetName(Form("ha1_scaled_%d",i));        
    ha1_adc_scaled[i]->Write();
    ha1_npe[i]->SetName(Form("ha1_npe_%d",i));                    
    ha1_npe[i]->Write();}
  for(int i=0;i<26;i++){
    ha2_adc[i]->SetName(Form("ha2_adc_%d",i));                
    ha2_adc[i]->Write();
    ha2_adc_scaled[i]->SetName(Form("ha2_scaled_%d",i));            
    ha2_adc_scaled[i]->Write();
    ha2_npe[i]->SetName(Form("ha2_npe_%d",i));                        
    ha2_npe[i]->Write();
  }

   ha1_adc_sum->Write();
   ha2_adc_sum->Write();
   ha1_npe_sum->Write();
   ha2_npe_sum->Write();
   tnew->Write();
 }

   fnew->Clone();
   
 if(draw_flag==0 || root_flag)gSystem->Exit(1);
 theApp->Run();
 return 0;

		    
}
