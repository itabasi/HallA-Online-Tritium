#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <TChain.h>
#include <TMinuit.h>
#include <iostream>
#include <fstream>
//#include "Param.h"




void holepattern(){


  string main = "angcalib_4th_0915_4";
  //    string main = "test_R";
  //  string main = "ang_RHRS_sieve_init";
  string end_root = ".root";
  string end_dat = ".dat";
  string ifname="../rootfiles/angcalib/" + main + ".root";
  string ofname="./param/" + main + ".dat";

  TFile* f = new TFile(ifname.c_str() );

  //  Tree* t = (TTree*)f->Get("T");

  TH2D* h2[10];
  //    for(int i=0;i<10;i++)h2[i] =(TH2D*)f->Get(Form("hss_%d",i));
  //       for(int i=0;i<10;i++)h2[i] =(TH2D*)f->Get(Form("h2_%d",i));
        for(int i=0;i<10;i++)h2[i] =(TH2D*)f->Get(Form("h2_new_%d",i));

					    

  const double step = 0.492 * 2.54;
  const int nrow = 11; // the number of row in SS pattern
  const int ncol = 8;  // the number of column in SS pattern
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
      /*
      ssx_off[0][45]=-1.1;
      ssx_off[0][47]=-0.7;
      ssx_off[0][53]= 0.6;
      ssx_off[0][54]= 0.85;
      ssx_off[0][55]=-1.3;
      ssx_off[0][56]=-1.05;
      ssx_off[0][57]=-0.7;
      ssx_off[0][58]=-0.6;
      ssx_off[0][64]= 0.85;
      ssx_off[0][65]= 0.85;
      ssx_off[0][66]=-1.3;
      ssx_off[0][67]=-1.05;
      ssx_off[0][74]= 0.6;
      ssx_off[0][76]= 1.0;

      ssx_off[1][44]= -1.0;
      ssx_off[1][45]= -0.7;
      ssx_off[1][46]= -0.8;
      ssx_off[1][55]= -1.0;
      ssx_off[1][56]= -0.7;
      ssx_off[1][66]= -1.0;
      ssx_off[1][68]= -0.8;
      */
      //      ssx_off[0][45]=-1.0;
      //      ssx_off[0][56]=-1.0;



      
      /*
      // RHRS initial matrix //      
      ssy_off[1][j]=-2.5;
      ssx_off[1][j]=-0.75; // x

      ssy_off[2][j]=-1.95;
      ssx_off[2][j]=-0.7; // x

      ssy_off[3][j]=-1.4;
      ssx_off[3][j]=-0.7; // x            

      ssy_off[4][j]=-0.9;
      ssx_off[4][j]=-0.8; // x      

      ssy_off[5][j]=-0.275;
      ssx_off[5][j]=-0.7; // x      

      ssy_off[6][j]= 0.2;
      ssx_off[6][j]=-0.9; // x

      ssy_off[7][j]= 0.8;
      ssx_off[7][j]=-0.7; // x      

      ssy_off[8][j]= 1.8;
      ssx_off[8][j]=-0.8; // x      

      ssy_off[9][j]= 2.35;
      ssx_off[9][j]=-0.7; // x      
      */

      /*
      ssy_off[0][45]= -0.25;
      ssy_off[0][47]= -0.25;
      ssy_off[0][49]= -0.25;
      ssy_off[0][51]= -0.25;
      ssy_off[0][53]= -0.25;
      ssy_off[0][55]= -0.20;
      ssy_off[0][56]= -0.25;
      ssy_off[0][57]= -0.20;
      ssy_off[0][58]= -0.25;
      ssy_off[0][59]= -0.20;
      ssy_off[0][60]= -0.25;
      ssy_off[0][61]= -0.20;
      ssy_off[0][62]= -0.25;
      ssy_off[0][63]= -0.20;
      ssy_off[0][64]= -0.25;

      ssy_off[0][66]= -0.25;
      ssy_off[0][67]= -0.30;
      ssy_off[0][68]= -0.25;
      ssy_off[0][69]= -0.30;
      ssy_off[0][70]= -0.25;
      ssy_off[0][71]= -0.30;
      ssy_off[0][72]= -0.25;
      ssy_off[0][73]= -0.30;
      ssy_off[0][74]= -0.25;
      ssy_off[0][75]= -0.30;
      ssy_off[0][76]= -0.25;

      ssy_off[0][77]= -0.4;
      ssy_off[0][78]= -0.5; 
      ssy_off[0][79]= -0.4;
      ssy_off[0][80]= -0.5; 
      ssy_off[0][81]= -0.4; 
      ssy_off[0][82]= -0.4; 
      ssy_off[0][83]= -0.4;
      ssy_off[0][84]= -0.2; 
      ssy_off[0][85]= -0.2;

      //==== foil0 x ===//

      ssx_off[0][55]= -0.9;
      ssx_off[0][66]= -0.9;
      
      ssx_off[0][45]= -0.8;
      ssx_off[0][56]= -0.8;
      ssx_off[0][67]= -0.8;
      ssx_off[0][78]= -0.7;

      ssx_off[0][46]= -0.5;
      ssx_off[0][57]= -0.5;
      ssx_off[0][68]= -0.5;
      ssx_off[0][79]= -0.5;

      ssx_off[0][42]=  0.55;
      ssx_off[0][53]=  0.55;
      ssx_off[0][54]=  0.7;
      ssx_off[0][64]=  0.55;
      ssx_off[0][65]=  0.7;
      ssx_off[0][75]=  0.55;
      ssx_off[0][76]=  0.7;
      
      //      ssy_off[0][86]= -0.5; 

      //====== foil 1 ========//
      ssx_off[1][45]=  -0.6;
      ssy_off[1][45]=  -0.1;
      ssx_off[1][56]=  -0.6;
      ssx_off[1][67]=  -0.6;
      ssy_off[1][67]=  -0.25;
      ssy_off[1][69]=  -0.3;
      ssx_off[1][55]=  -0.2;
      ssy_off[1][55]=  -0.2;
      ssy_off[1][66]=  -0.25;
      ssx_off[1][66]=  -0.8;
      ssy_off[1][68]=  -0.25;
      ssy_off[1][78]=  -0.3;
      ssy_off[1][79]=  -0.3;
      ssy_off[1][80]=  -0.3;
      ssy_off[1][81]=  -0.3;
      ssy_off[1][82]=  -0.3;
      ssy_off[1][83]=  -0.2;
      ssy_off[1][84]=  -0.3;
      ssx_off[1][54]=   0.6;
      //      ssy_off[1][85]=  -0.3;
      //      ssy_off[1][86]=  -0.3;
      */

      
      if((1 <= j && j <= 88) //&& (j % 11 != 0)
	 //	 && (1 <= i && i <= 8)){flag[i][j]=true;}
	 && (0 <= i && i <= 9)){flag[i][j]=true;}

      //	 && (4 <= i && i <= 7)){flag[i][j]=true;}


      //=====foil 0 ========//
      for(int k=0;k<=33;k++)flag[0][k]=false;
      
      flag[0][34] =false;
      flag[0][35] =false;
      flag[0][37] =false;
      flag[0][39] =false;
      flag[0][43] =false;
      flag[0][44] =false;
      
      //      flag[0][45] =false;
      //      flag[0][54] =false;
      //      flag[0][56] =false;
      //      flag[0][65] =false;
      //      flag[0][67] =false;
      flag[0][77] =false;
      flag[0][78] =false;
      flag[0][87] =false;
      flag[0][86] =false;

      //=====foil 1 ========//
      //      flag[1][38]=true;
      for(int k=0;k<=32;k++)flag[1][k]=false;
      //      for(int k=77;k<=88;k++)flag[1][k]=false;
      flag[1][22] =false;
      flag[1][23] =false;
      flag[1][24] =false;
      flag[1][25] =false;
      flag[1][26] =false;
      flag[1][28] =false;
      flag[1][30] =false;
      flag[1][32] =false;
      flag[1][33] =false;
      flag[1][35] =false;

      //      flag[1][42] =false;
      flag[1][43] =false;
      //      flag[1][54] =false;
      flag[1][65] =false;
      //      flag[1][67] =false;
      //      flag[1][75] =false;
      flag[1][76] =false;
      flag[1][77] =false;
      flag[1][78] =false;
      //      flag[1][79] =true;
      //      flag[1][81] =true;
      //      flag[1][83] =true;
      //      flag[1][85] =true;      
      flag[1][86] =false;
      flag[1][87] =false;
      //      flag[1][33] =true;
      //      flag[1][44] =true;
      
      //===== foil 2 =======//
      for(int k=0;k<=22;k++)flag[2][k]=false;
      for(int k=77;k<=88;k++)flag[2][k]=false;
      flag[2][22] =false;
      flag[2][23] =false;
      flag[2][24] =false;
      flag[2][26] =false;
      flag[2][28] =false;
      flag[2][30] =false;
      flag[2][32] =false;
      flag[2][43] =false;
      flag[2][54] =false;
      //      flag[2][63] =false;
      //      flag[2][65] =false;
      //      flag[2][67] =false;
      //      flag[2][75] =false;
      flag[2][76] =false;

      flag[2][44] =true;
      flag[2][55] =true;
      flag[2][78] =true;
      flag[2][79] =true;
      flag[2][81] =true;
      flag[2][83] =true;
      flag[2][85] =true;
      //===== foil 3 =======//

      for(int k=0;k<=21;k++)flag[3][k]=false;
      for(int k=77;k<=88;k++)flag[3][k]=false;
      //      flag[3][22] =false;
      //      flag[3][24] =false;
      //      flag[3][26] =false;
      //      flag[3][28] =false;
      //      flag[3][30] =false;
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

      //      flag[3][22] =true;
      flag[3][33] =true;
      flag[3][44] =true;
      flag[3][55] =true;
      flag[3][81] =true;
      flag[3][83] =true;
      //======== foil 4 =========//
      for(int k=0;k<=21;k++)flag[4][k]=false;
      for(int k=77;k<=88;k++)flag[4][k]=false;
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

      flag[4][16] =true;
      //       flag[4][44] =true;
      //       flag[4][55] =true;

       //========= foil 5 =======//
      for(int k=0;k<=11;k++)flag[5][k]=false;
      for(int k=77;k<=88;k++)flag[5][k]=false;
      //      for(int k=64;k<=nsshole;k++)flag[5][k]=false;
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
      flag[5][67] =false;
      flag[5][69] =false;
      flag[5][71] =false;
      flag[5][73] =false;
      flag[5][75] =false;
      flag[5][76] =false;

      //      flag[5][33] =true;
      //      flag[5][44] =true;

      
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

      //      flag[6][33] =true;
      //      flag[6][44] =true;


      //========== foil 7 =========//

      //      for(int k=0;k<=11;k++)flag[7][k]=false;
      for(int k=66;k<=nsshole;k++)flag[7][k]=false;
      flag[7][0] =false;
      flag[7][1] =false;
      flag[7][2] =false;
      flag[7][4] =false;
      flag[7][6] =false;
      //      flag[7][8] =false;
      flag[7][9] =false;
      flag[7][10] =false;
      //      flag[7][12] =false;

      //      flag[7][12] =false;
      //      flag[7][13] =false;
      //      flag[7][15] =false;
      //      flag[7][17] =false;
      //      flag[7][19] =false;
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

      //      flag[7][33] =true;

      
       
       //========= foil 8 =======//

      flag[8][11]  =false;
      flag[8][22]  =false;
      flag[8][33]  =false;
      flag[8][44]  =false;
      flag[8][0]  =false;
      flag[8][1]  =false;			   
      flag[8][2]  =false;
      //      flag[8][4]  =false;
      //      flag[8][6]  =false;
      flag[8][8]  =false;
      flag[8][9]  =false;
      flag[8][10] =false;
      //      flag[8][20] =false;
      flag[8][21] =false;
      //      flag[8][31] =false;
      flag[8][32] =false;
      //      flag[8][34] =false;
      //      flag[8][42] =false;
      flag[8][43] =false;
      flag[8][45] =false;
      //      flag[8][46] =false;
      flag[8][47] =false;
      flag[8][49] =false;
      flag[8][51] =false;
      flag[8][53] =false;
      flag[8][54] =false;

      //      flag[8][22] =true;
      //      flag[8][33] =true;
      //      flag[8][44] =true;

      for(int k=55;k<nsshole;k++){flag[8][k]=false;}

      //====== foil 9 =========//
      
      for(int k=43;k<nsshole;k++){flag[9][k]=false;}
      //      flag[9][1]=false;
      flag[9][2]=false;
      //      flag[9][3]=false;
      //      flag[9][8]=false;
      flag[9][9] =false;
      flag[9][10]=false;
      //      flag[9][13]=false;
      //      flag[9][19]=false;
      flag[9][21]=false;
      flag[9][32]=false;
      flag[9][42]=false;
      flag[9][43]=false;

      
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

	//=== weight ======//
	//      if(0<=i && i<=9)
	w[i][j]=1.0;
	if( 7<=i ) w[i][j]=5.0; 
       	if( i<=1 ) w[i][j]=20.0; 
	//if(i==1 || i==8)w[i][j]=10.0; 
      //      if(i<=1 && 67<=j ) w[i][j]=10.0;
      //      if(i<=8 && j<=10 ) w[i][j]=10.0;

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


  TCanvas* c3=new TCanvas("c3","c3");  
  c3->Divide(1,2);  
  c3->cd(1);
  h2[0]->Draw("colz");
  c3->cd(2);
  h2[1]->Draw("colz");

  
  
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

	  c3->cd(1);
	  mark[j]->Draw("same");	  
	  if(flag[0][j])mark_real[0][j]->Draw("same");
	  c3->cd(2);
	  mark[j]->Draw("same");	  
	  if(flag[1][j])mark_real[1][j]->Draw("same");
	  
	  
	}
	c1->cd(i+1);
	mark[j]->Draw("same");
      }
    }

  

  
}
