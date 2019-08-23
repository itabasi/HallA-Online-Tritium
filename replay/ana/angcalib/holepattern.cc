#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <TChain.h>
#include <TMinuit.h>
#include <iostream>
#include <fstream>
//#include "Param.h"




void holepattern(){
  string main = "ang_RHRS_4th_0824_0_0";
  string end_root = ".root";
  string end_dat = ".dat";
  string ifname="./rootfile/" + main + ".root";
  string ofname="./param/" + main + "dat";

  TFile* f = new TFile(ifname.c_str() );

  //  Tree* t = (TTree*)f->Get("T");

  TH2D* h2[10];
  for(int i=0;i<10;i++)h2[i] =(TH2D*)f->Get(Form("h2_%d",i));


					    

  const double step = 0.492 * 2.54;
  const int nrow = 11; // the number of row in SS pattern
  const int ncol = 7;  // the number of column in SS pattern
  const int nsshole = nrow*ncol; // the number of holes to consider 
  const int nfoil =10;



  
  TMarker* mark[nsshole];
  TMarker* mark_real[nfoil][nsshole];  
  double w[nfoil][nsshole];
  double ssx_off[nfoil][nsshole];
  double ssy_off[nfoil][nsshole];
  bool flag[nfoil][nsshole];
  double refx[nsshole],refy[nsshole],refx_real[nfoil][nsshole],refy_real[nfoil][nsshole];
  double ssy_cent_real[nsshole],ssx_cent_real[nsshole];
  
  int nhole=0;
   for(int i=0; i<ncol ; i++){
    for(int j=0; j<nrow; j++){
      ssy_cent_real[i] = -3.0*step + step*i;
      if(j%2==0)ssy_cent_real[i] = ssy_cent_real[i] - step/2.0;
      ssx_cent_real[j] = 5.0*step - step*j;
      refx[nhole] = ssx_cent_real[j];
      refy[nhole] = ssy_cent_real[i];
      nhole++;
    }
   }


   

  for(int i=0;i<nfoil;i++){
    for(int j=0;j<nsshole;j++){
      
      ssx_off[i][j]=0.0;
      ssy_off[i][j]=0.0;
      w[i][j]=1.0;
      flag[i][j]=false;
      
      
      
      if((1 <= j && j <= 76) && (j % 11 != 0)
	   && (3 <= i && i <= 7)){flag[i][j]=true;}
      
      
      
    }
    
  }


  
    for(int j=0;j<nsshole;j++){

      mark[j] = new TMarker(refy[j],refx[j],28);
      mark[j]->SetMarkerColor(1);

      for(int i=0;i<nfoil;i++){

	refx_real[i][j] =refx[j]+ssx_off[i][j];
	refy_real[i][j] =refy[j]+ssy_off[i][j];      

	mark_real[i][j] = new TMarker(refy_real[i][j],refx_real[i][j],20);
	mark_real[i][j] -> SetMarkerColor(1);

      }

    }


    	mark[23]->SetMarkerColor(2);
	mark[38]->SetMarkerColor(2);

      for(int i=0;i<nfoil;i++){
      	mark_real[i][23] -> SetMarkerColor(2);
	mark_real[i][38] -> SetMarkerColor(2);  }

      

    ofstream * ofs=new ofstream(ofname.c_str());

  *ofs << "# "<<ifname.c_str()<<" sieve slit offset parameters "<<endl;
  *ofs << "### "<<"( #foil, #hole, flag, w, ssx_off, ssy_off ) ### "<<endl;
  for(int i=0;i<nfoil;i++){  
    for(int j=0;j<nsshole;j++){

      *ofs <<  i << " " << j << " " << flag[i][j] << " " <<
	w[i][j] <<  " " << ssx_off[i][j] << " " <<  ssy_off[i][j] <<endl;  
    }
  }
  
  ofs->close();

  TCanvas* c1=new TCanvas("c1","c1");  
  c1->Divide(4,3);
  for(int i=0;i<nfoil;i++){
    c1->cd(i+1);
    h2[i]->Draw("colz");
  }


  
  TCanvas* c0=new TCanvas("c0","c0");  
  c0->Divide(2,2);  
  c0->cd(1);
  h2[4]->Draw("colz");
  c0->cd(2);
  h2[5]->Draw("colz");
  c0->cd(3);
  h2[6]->Draw("colz");
  c0->cd(4);
  h2[7]->Draw("colz");

  TCanvas* c2=new TCanvas("c2","c2");  
  c2->Divide(2,2);  
  c2->cd(1);
  h2[2]->Draw("colz");
  c2->cd(2);
  h2[3]->Draw("colz");
  c2->cd(3);
  h2[8]->Draw("colz");
  c2->cd(4);
  h2[9]->Draw("colz");
  

  
  
    for(int j=0;j<nsshole;j++){
      for(int i=0;i<nfoil;i++){
	if(flag[i][j]){
	  c0->cd(1);
	  mark[j]->Draw("same");
 	  mark_real[4][j]->Draw("same");
	  c0->cd(2);
	  mark[j]->Draw("same");
	  mark_real[5][j]->Draw("same");
	  c0->cd(3);
	  mark[j]->Draw("same");
	  mark_real[6][j]->Draw("same");
	  c0->cd(4);
	  mark[j]->Draw("same");
	  mark_real[7][j]->Draw("same");
	  c2->cd(1);
	  mark[j]->Draw("same");
	  mark_real[2][j]->Draw("same");
	  c2->cd(2);
	  mark[j]->Draw("same");
	  mark_real[3][j]->Draw("same");
	  c2->cd(3);
	  mark[j]->Draw("same");
	  mark_real[8][j]->Draw("same");
	  c2->cd(4);
	  mark[j]->Draw("same");	  
	  mark_real[9][j]->Draw("same");
	}
	c1->cd(i+1);
	mark[j]->Draw("same");
      }
    }

  

  
}
