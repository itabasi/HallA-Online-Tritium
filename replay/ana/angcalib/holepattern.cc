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

      //====== offset paramters======= //

      ssx_off[3][29]= 0.8;
      ssx_off[3][31]= 0.8;
      ssx_off[3][40]= 0.2;

      ssx_off[3][56]= -0.75;
      ssx_off[3][68]= -0.6;
      ssx_off[3][69]= -0.8;
      ssx_off[3][70]= -0.7;

      ssy_off[3][34]= -0.2;
      ssy_off[3][36]= -0.2;
      ssy_off[3][38]= -0.2;
      ssy_off[3][40]= -0.2;
      ssy_off[3][42]= -0.2;
      ssy_off[3][45]= -0.2;
      ssy_off[3][46]= -0.2;
      ssy_off[3][47]= -0.2;
      ssy_off[3][48]= -0.2;
      ssy_off[3][49]= -0.2;
      ssy_off[3][50]= -0.2;
      ssy_off[3][51]= -0.2;
      ssy_off[3][52]= -0.2;
      ssy_off[3][53]= -0.2;


      ssy_off[3][56]= -0.45;
      ssy_off[3][57]= -0.25;
      ssy_off[3][58]= -0.5;
      ssy_off[3][59]= -0.25;
      ssy_off[3][60]= -0.5;
      ssy_off[3][61]= -0.25;
      ssy_off[3][62]= -0.5;
      ssy_off[3][63]= -0.25;
      ssy_off[3][64]= -0.5;


      ssy_off[3][68]= -0.65;
      ssy_off[3][69]= -1.25;
      ssy_off[3][70]= -0.5;
      ssy_off[3][71]= -1.0;
      ssy_off[3][72]= -0.65;
      ssy_off[3][73]= -1.25;
      ssy_off[3][74]= -0.85;

      ssy_off[4][64]= -0.4;
      ssy_off[4][68]= -0.35;
      ssy_off[4][69]= -0.7;
      ssy_off[4][70]= -0.3;
      ssy_off[4][71]= -0.7;
      ssy_off[4][72]= -0.35;
      ssy_off[4][73]= -0.7;
      ssy_off[4][74]= -0.4;
      
      ssy_off[8][3]= 0.75;
      ssy_off[8][5]= 0.75;
      ssy_off[8][7]= 1.25;
      ssy_off[8][12]= 0.8;
      ssy_off[8][13]= 0.8;
      ssy_off[8][14]= 0.8;
      ssy_off[8][15]= 0.8;
      ssy_off[8][16]= 0.8;
      ssy_off[8][17]= 0.8;
      ssy_off[8][18]= 0.8;
      ssy_off[8][19]= 0.8;
      ssy_off[8][20]= 0.8;

      ssy_off[8][23]= 0.7;
      ssy_off[8][24]= 0.6;
      ssy_off[8][25]= 0.7;
      ssy_off[8][26]= 1.0;
      ssy_off[8][27]= 0.85;
      ssy_off[8][28]= 1.0;
      ssy_off[8][29]= 0.9;
      ssy_off[8][30]= 1.0;
      ssy_off[8][35]= 0.65;
      ssy_off[8][36]= 1.0;
      ssy_off[8][37]= 0.8;
      ssy_off[8][38]= 0.9;
      ssy_off[8][39]= 1.0;
      ssy_off[8][40]= 1.0;
      ssy_off[8][41]= 1.0;
      ssy_off[8][42]= 1.0;

      ssy_off[8][48]= 0.7;
      ssy_off[8][50]= 0.9;
      ssy_off[8][52]= 1.0;
      
      if((1 <= j && j <= 76) && (j % 11 != 0)
	   && (3 <= i && i <= 8)){flag[i][j]=true;}


      //===== foil 3 =======//

      for(int k=0;k<=21;k++)flag[3][k]=false;

      flag[3][22] =false;
      flag[3][24] =false;
      flag[3][26] =false;
      flag[3][28] =false;
      flag[3][30] =false;
      flag[3][32] =false;
      flag[3][43] =false;
      flag[3][54] =false;
      flag[3][65] =false;
      flag[3][67] =false;
      //      flag[3][69] =false;
      //      flag[3][71] =false;
      //      flag[3][73] =false;
      flag[3][75] =false;
      flag[3][76] =false;

      //======== foil 4 =========//
      for(int k=0;k<=21;k++)flag[4][k]=false;
      flag[4][12] =false;
      flag[4][20] =false;
      flag[4][21] =false;
      flag[4][32] =false;
      flag[4][43] =false;
      flag[4][21] =false;
      flag[4][65] =false;
      flag[4][67] =false;
      flag[4][75] =false;
      flag[4][76] =false;

       flag[4][33] =true;
       flag[4][44] =true;
       flag[4][55] =true;

       //========= foil 5 =======//
      for(int k=0;k<=11;k++)flag[5][k]=false;
      for(int k=66;k<=nsshole;k++)flag[5][k]=false;
      flag[5][11] =false;
      flag[5][12] =false;
      flag[5][13] =false;
      flag[5][15] =false;
      flag[5][17] =false;
      flag[5][19] =false;
      flag[5][20] =false;
      flag[5][21] =false;
      flag[5][32] =false;
      flag[5][43] =false;
      flag[5][54] =false;
      flag[5][21] =false;
      flag[5][65] =false;
      //      flag[5][67] =false;
      //      flag[5][75] =false;
      //      flag[5][76] =false;

      flag[5][33] =true;
      flag[5][44] =true;

      
      //========== foil 6 ========//
      for(int k=0;k<=11;k++)flag[6][k]=false;
      for(int k=66;k<=nsshole;k++)flag[6][k]=false;
      flag[6][11] =false;
      //      flag[6][12] =false;
      //      flag[6][13] =false;
      //      flag[6][15] =false;
      //      flag[6][17] =false;
      flag[6][19] =false;
      flag[6][20] =false;
      flag[6][21] =false;
      flag[6][32] =false;
      flag[6][43] =false;
      flag[6][54] =false;
      flag[6][64] =false;
      flag[6][65] =false;
      flag[6][67] =false;
      flag[6][75] =false;
      flag[6][76] =false;

      flag[6][33] =true;
      flag[6][44] =true;


      //========== foil 7 =========//

      for(int k=0;k<=11;k++)flag[7][k]=false;
      for(int k=66;k<=nsshole;k++)flag[7][k]=false;
      flag[7][11] =false;
      //      flag[7][12] =false;
      //      flag[7][13] =false;
      //      flag[7][15] =false;
      //      flag[7][17] =false;
      flag[7][19] =false;
      //      flag[7][20] =false;
      flag[7][21] =false;
      flag[7][32] =false;
      flag[7][43] =false;
      flag[7][54] =false;
      flag[7][56] =false;
      flag[7][58] =false;
      flag[7][60] =false;
      flag[7][62] =false;
      flag[7][64] =false;
      flag[7][65] =false;
      flag[7][67] =false;
      flag[7][75] =false;
      flag[7][76] =false;

      flag[7][33] =true;

      
       
       //========= foil 8 =======//
      
      flag[8][0]  =false;
      flag[8][1]  =false;			   
      flag[8][2]  =false;
      flag[8][4]  =false;
      flag[8][6]  =false;
      flag[8][8]  =false;
      flag[8][9]  =false;
      flag[8][10] =false;
      flag[8][20] =false;
      flag[8][21] =false;
      flag[8][31] =false;
      flag[8][32] =false;
      flag[8][34] =false;
      flag[8][42] =false;
      flag[8][43] =false;
      flag[8][45] =false;
      flag[8][46] =false;
      flag[8][47] =false;
      flag[8][49] =false;
      flag[8][51] =false;
      flag[8][53] =false;
      flag[8][54] =false;
      
      for(int k=55;k<nsshole;k++){flag[8][k]=false;}


      
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

      if(5<=i && i<=7) w[i][j]=1.0;
      else if(flag[i][j])w[i][j]=3.0;

	
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
	  if(flag[4][j])mark_real[4][j]->Draw("same");
	  c0->cd(2);
	  mark[j]->Draw("same");
	  if(flag[5][j])mark_real[5][j]->Draw("same");
	  c0->cd(3);
	  mark[j]->Draw("same");
	  if(flag[6][j])mark_real[6][j]->Draw("same");
	  c0->cd(4);
	  mark[j]->Draw("same");
	  if(flag[7][j])mark_real[7][j]->Draw("same");
	  c2->cd(1);
	  mark[j]->Draw("same");
	  if(flag[2][j])mark_real[2][j]->Draw("same");
	  c2->cd(2);
	  mark[j]->Draw("same");
	  if(flag[3][j])mark_real[3][j]->Draw("same");
	  c2->cd(3);
	  mark[j]->Draw("same");
	  if(flag[8][j])mark_real[8][j]->Draw("same");
	  c2->cd(4);
	  mark[j]->Draw("same");	  
	  if(flag[9][j])mark_real[9][j]->Draw("same");
	}
	c1->cd(i+1);
	mark[j]->Draw("same");
      }
    }

  

  
}
