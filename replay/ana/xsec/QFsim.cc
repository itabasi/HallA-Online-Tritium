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
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TColor.h>
#include <TPaveText.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;
using namespace Eigen;
#include "QFsim.h"
#include "define.h"
#include "Param.h"
#include "Math/ProbFunc.h"
#include "Math/DistFunc.h"

////////////////////////////////////////////////////////////////////



void QFsim::SetRootList(string ifrname){

  cout<<endl;
  cout<<"==============================="<<endl;
  cout<<"=====< Set Root File >=========="<<endl;
  cout<<"==============================="<<endl;


  string buf;
  int s=0;
  ifstream ifp(Form("%s",ifrname.c_str()),ios::in);
  if (ifp.fail()){ cerr << "Failed read file" <<ifrname.c_str()<<endl; exit(1);}
  string rname[100];
  string names;
  int id;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    //    sbuf >>rname[s];
    sbuf >>id >>names;
    rname[id] = names;
    cout<<" id "<<id <<rname[id]<<endl;
    s++;
  }
  


  fexp   = new TFile(rname[0].c_str()); // Exp. Root File
  fH_sim = new TFile(rname[1].c_str()); // SIMC Lambda Root File
  fT_sim = new TFile(rname[2].c_str()); // SIMC nnL    Root File
  
  fH_sim_10keV = new TFile(rname[3].c_str()); // SIMC Lam(10 keV) Root File  
  fT_sim_10keV = new TFile(rname[4].c_str()); // SIMC nnL(10 keV) Root File
  fexp_10keV = new TFile(rname[5].c_str()); // Exp. Root File (10 keV)
  
  ofr    = new TFile(rname[6].c_str(),"recreate"); // OutPut Root File
  
  Tnew = new TTree("T","QF Fitting ");
  cout<<"Out Put Root File "<<rname[6]<<endl;
}


//////////////////////////////////////////////////////////////////////////


void QFsim::SetHist(){

  // Experimental Data
  hexp_peak_L = (TH1D*)fexp -> Get("h_peak_L");
  hexp_peak_L ->SetName("hexp_peak_L");
  hexp_peak = (TH1D*)fexp -> Get("h_peak_nnL");
  hexp_peak ->SetName("hexp_peak");  
  hexp        = (TH1D*)fexp -> Get("h_peak_nnL");
  hexp ->SetName("hexp");
  hexp->SetTitle("Exp data w/o Accidental B.G.; Missing mass [GeV]; Counts/2MeV");
  hexp_c        = (TH1D*)hexp -> Clone();
  hexp_c ->SetName("hexp_c");
  hexp_c->SetTitle("Exp data w/ Accidental B.G.; Missing mass [GeV]; Counts/2MeV");
  hexp_acc        = (TH1D*)fexp -> Get("h_acc_nnL");
  hexp_acc ->SetName("hexp_acc");
  hexp_acc->SetTitle("Exp data Accidental B.G.; Missing mass [GeV]; Counts/2MeV");  
  hexp_nnL   = (TH1D*)fexp -> Get("h_mm_nnL");
  hexp_nnL ->SetName("hexp_nnL");
  hexp_nnL->SetTitle("Exp data w/ Accidental B.G.(total); Missing mass [GeV]; Counts/2MeV");  
  hexp_L     = (TH1D*)fexp -> Get("h_peak_L");
  hexp_L ->SetName("hexp_L");
  //  hexp->SetTitle("Exp data (Hydrogen data) w/o Accidental B.G.; Missing mass [GeV]; Counts/2MeV");

  hmm       =(TH1D*)fT_sim->Get("hmm");
  hmm       ->SetName("hmm");
  hmm->SetTitle("Simulation w/o B.G.; Missing mass [GeV]; Counts/2MeV");
  hmm_fsi1  =(TH1D*)fT_sim->Get("hmm_fsi1");
  hmm_fsi1->SetName("hmm_fsi1");
  hmm_fsi1->SetTitle("Simulation w/o B.G. with model 1 FSI; Missing mass [GeV]; Counts/2MeV");  
  hmm_fsi2  =(TH1D*)fT_sim->Get("hmm_fsi2");
  hmm_fsi2->SetName("hmm_fsi2");
  hmm_fsi2->SetTitle("Simulation w/o B.G. with model 2 FSI; Missing mass [GeV]; Counts/2MeV");  
  hmm_fsi3  =(TH1D*)fT_sim->Get("hmm_fsi3");
  hmm_fsi3->SetName("hmm_fsi3");
  hmm_fsi3->SetTitle("Simulation w/o B.G. with model 3 FSI; Missing mass [GeV]; Counts/2MeV");  

  hmm_s       = (TH1D*)hmm      -> Clone();
  hmm_s       ->SetName("hmm_s");
  hmm_s->SetTitle("Simulation w/o B.G. with scaling ; Missing mass [GeV]; Counts/2MeV");    
  hmm_fsi1_s  = (TH1D*)hmm_fsi1 -> Clone();
  hmm_fsi1_s ->SetName("hmm_fsi1_s");
  hmm_fsi1_s->SetTitle("Simulation w/o B.G. with scaling (model 1 FSI); Missing mass [GeV]; Counts/2MeV");    
  hmm_fsi2_s  = (TH1D*)hmm_fsi2 -> Clone();
  hmm_fsi2_s ->SetName("hmm_fsi2_s");
  hmm_fsi2_s->SetTitle("Simulation w/o B.G. with scaling (model 2 FSI); Missing mass [GeV]; Counts/2MeV");      
  hmm_fsi3_s  = (TH1D*)hmm_fsi3-> Clone();
  hmm_fsi3_s ->SetName("hmm_fsi3_s");
  hmm_fsi3_s->SetTitle("Simulation w/o B.G. with scaling (model 3 FSI); Missing mass [GeV]; Counts/2MeV");    



  //  hexp_10keV =(TH1D*)fexp_10keV->Get("hmm");


  
  hexp_10keV =(TH1D*)fexp_10keV->Get("hmm_peak");  // w/o ACC. B.G.
  hexp_10keV -> SetName("hexp_10keV");
  hexp_10keV_s =(TH1D*)hexp_10keV->Clone();  // w/o ACC. B.G.
  hexp_10keV_s -> SetName("hexp_10keV_s");
  hexp_10keV_c =(TH1D*)hexp_10keV->Clone();  // w/o ACC. B.G.
  hexp_10keV_c -> SetName("hexp_10keV_c");  
  
  hexp_acc_10keV =(TH1D*)fexp_10keV->Get("hmm_acc");
  hexp_acc_10keV -> SetName("hexp_acc_10keV");

  
  hmm_10keV =(TH1D*)fT_sim_10keV->Get("hmm");
  hmm_10keV -> SetName("hmm_10keV");
  hmm_10keV_s =(TH1D*)hmm_10keV->Clone();
  hmm_10keV_s -> SetName("hmm_10keV_s");  
  hmm_10keV_c =(TH1D*)hmm_10keV->Clone();
  hmm_10keV_c -> SetName("hmm_10keV_c");
  
  hH_nnLsim_10keV  = (TH1D*)fH_sim_10keV->Get("hmm_nnL");
  hH_nnLsim_10keV -> SetName("hH_nnLsim_10keV");
  hH_nnLsim_10keV_s =(TH1D*)hH_nnLsim_10keV->Clone();
  hH_nnLsim_10keV_s -> SetName("hH_nnLsim_10keV_s");  
  
  hmm_fsi1_10keV =(TH1D*)fT_sim_10keV->Get("hmm_fsi1");
  hmm_fsi1_10keV -> SetName("hmm_fsi1_10keV");
  hmm_fsi2_10keV =(TH1D*)fT_sim_10keV->Get("hmm_fsi2");
  hmm_fsi2_10keV -> SetName("hmm_fsi2_10keV");
  hmm_fsi3_10keV =(TH1D*)fT_sim_10keV->Get("hmm_fsi3");
  hmm_fsi3_10keV -> SetName("hmm_fsi3_10keV");


  
  hexp_10keV_s =(TH1D*)hexp_10keV->Clone();
  hexp_10keV_s -> SetName("hexp_10keV_s");
  hmm_10keV_s =(TH1D*)hmm_10keV->Clone();
  hmm_10keV_s -> SetName("hmm_10keV_s");
  hmm_fsi1_10keV_s =(TH1D*)hmm_fsi1_10keV->Clone();
  hmm_fsi1_10keV_s -> SetName("hmm_fsi1_10keV_s");
  hmm_fsi2_10keV_s =(TH1D*)hmm_fsi2_10keV->Clone();
  hmm_fsi2_10keV_s -> SetName("hmm_fsi2_10keV_s");
  hmm_fsi3_10keV_s =(TH1D*)hmm_fsi3_10keV->Clone();
  hmm_fsi3_10keV_s -> SetName("hmm_fsi3_10keV_s");
  
  hexp_10keV_c =(TH1D*)hexp_10keV->Clone();
  hexp_10keV_c -> SetName("hexp_10keV_c");
  hmm_10keV_c =(TH1D*)hmm_10keV->Clone();
  hmm_10keV_c -> SetName("hmm_10keV_c");
  hmm_fsi1_10keV_c =(TH1D*)hmm_fsi1_10keV->Clone();
  hmm_fsi1_10keV_c -> SetName("hmm_fsi1_10keV_c");
  hmm_fsi2_10keV_c =(TH1D*)hmm_fsi2_10keV->Clone();
  hmm_fsi2_10keV_c -> SetName("hmm_fsi2_10keV_c");
  hmm_fsi3_10keV_c =(TH1D*)hmm_fsi3_10keV->Clone();
  hmm_fsi3_10keV_c -> SetName("hmm_fsi3_10keV_c");
  

  hH_Lsim   = (TH1D*)fH_sim ->Get("hmm_L");
  hH_Lsim -> SetName("hH_Lsim");
  hH_Lsim_s  =(TH1D*)hH_Lsim->Clone();
  hH_Lsim_s -> SetName("hH_Lsim_s");
  hH_nnLsim   = (TH1D*)fH_sim ->Get("hmm_nnL");
  hH_nnLsim  -> SetName("hH_nnLsim");
  hH_nnLsim_s =(TH1D*)hH_nnLsim->Clone();
  hH_nnLsim_s -> SetName("hH_nnLsim_s"); 
  
  hT_Lsim   = (TH1D*)fT_sim ->Get("hmm_L");
  hT_Lsim -> SetName("hT_Lsim");
  hT_Lsim_s =(TH1D*)hT_Lsim->Clone();
  hT_Lsim_s -> SetName("hT_Lsim_s");
  hT_nnLsim   = (TH1D*)fT_sim ->Get("hmm_nnL");
  hT_nnLsim  -> SetName("hT_nnLsim");
  hT_nnLsim_s =(TH1D*)hT_nnLsim->Clone();
  hT_nnLsim_s -> SetName("hT_nnLsim_s"); 

  hnnL_sim   =(TH1D*)fT_sim->Get("hmm");
  hnnL_sim ->SetName("hnnL_sim");
  hnnL_sim_s =(TH1D*)hnnL_sim->Clone();
  hnnL_sim_s ->SetName("hnnL_sim_s");
  
  gchi_L = new TGraphErrors();
  gchi_L->SetName("gchi_L");
  gchi_L->SetMarkerStyle(7);
  
  gchi_nnL = new TGraphErrors();
  gchi_nnL->SetName("gchi_nnL");  
  gchi_nnL->SetMarkerStyle(7);

  gdiff_L = new TGraphErrors();
  gdiff_L->SetName("gdiff_L");
  gdiff_L->SetMarkerStyle(7);
  
  for(int i=0;i<nvp;i++){
  gchi_fsi[i] = new TGraphErrors();
  gchi_fsi[i]->SetName(Form("gchi_fsi_%d",i));
  gdiff_fsi[i] = new TGraphErrors();
  gdiff_fsi[i]->SetName(Form("gdiff_fsi_%d",i));


  g_s[i] = new TGraph();
  g_s[i]->SetName(Form("g_s_%d",i));
  g_n[i] = new TGraph();
  g_n[i]->SetName(Form("g_n_%d",i));
  g_val[i] = new TGraph();
  g_val[i]->SetName(Form("g_val_%d",i));
  g_ps[i] = new TGraph();
  g_ps[i]->SetName(Form("g_ps_%d",i));

  g_ps_sig[i] = new TGraph();
  g_ps_sig[i]->SetName(Form("g_ps_sig_%d",i));  
  g_val_sig[i] = new TGraph();
  g_val_sig[i]->SetName(Form("g_val_sig_%d",i));
  g_mm_sig[i] = new TGraph();
  g_mm_sig[i]->SetName(Form("g_mm_sig_%d",i));
  
  g_ps[i]->SetMarkerStyle(7);
  g_ps[i]->SetMarkerColor(i+1);
  g_val[i]->SetMarkerStyle(7);
  g_val[i]->SetMarkerColor(i+1);

  g_ps_sig[i]->SetMarkerStyle(7);
  g_ps_sig[i]->SetMarkerColor(i+1);
  g_val_sig[i]->SetMarkerStyle(7);
  g_val_sig[i]->SetMarkerColor(i+1);  
  g_mm_sig[i]->SetMarkerStyle(7);
  g_mm_sig[i]->SetMarkerColor(i+1);  
  
  
  }



  //=== CalcPS =====//
  for(int i=0;i<nvp;i++){

  g_val[i]->SetMinimum(1E-10);
  g_val[i]->GetXaxis()->SetRangeUser(-100,150);
  g_ps[i] ->GetXaxis()->SetRangeUser(-100,150);
  g_s[i]  ->GetXaxis()->SetRangeUser(-100,150);
  g_n[i]  ->GetXaxis()->SetRangeUser(-100,150);

  }

  
}


////////////////////////////////////////////////////////////


void QFsim::FitQFL(){


  int nbins = hexp_peak_L->GetXaxis()->GetNbins();
  int nmin  = hexp_peak_L->GetXaxis()->GetXmin();
  int nmax  = hexp_peak_L->GetXaxis()->GetXmax();

  double ymax0 = hexp_peak_L -> GetBinContent(hexp_peak_L->GetMaximumBin());
  double ymax1 = hT_Lsim  -> GetBinContent(hT_Lsim->GetMaximumBin());

  double y0[nbins],y1[nbins];
  //  int bin_x0=0;
  //  int bin_x1=nbins;;
  double x0;
  bool max= true;
  int bin_x0 = hexp_peak_L->GetXaxis()->FindBin(-50.0);
  int bin_x1 = hexp_peak_L->GetXaxis()->FindBin( 80.0);
  int bin_x2 = hexp_peak_L->GetXaxis()->FindBin(-10.0);
  int bin_x3 = hexp_peak_L->GetXaxis()->FindBin( 10.0);
  for(int ibin = 0;ibin<nbins;ibin++){
    x0        =  hexp_peak_L  -> GetBinCenter(ibin);
    y0[ibin]  =  hexp_peak_L  -> GetBinContent(ibin);
    y1[ibin]  =  hT_Lsim      -> GetBinContent(ibin);
    //    cout<<"ibin "<<ibin<<" x0 "<<x0<<" y0 "<<y0[ibin]<<" y1 "<<y1[ibin]<<endl;
  }


  double chi2,chi2_min,w,wmin;
  int wmax=1000;
  double diff[nbins];
  double width = 0.01;
  chi2_min=1e20;
  for(int wi=0; wi<wmax;wi++){
    chi2=0.0;
    w = ymax0/ymax1*(1.0 +(double)(-wmax/2. + wi)*(double)width);
    for(int ibin = bin_x0;ibin<bin_x1;ibin++){
      if(bin_x2< ibin && ibin < bin_x3)continue; // Lambda peak
      //      if(y0[ibin]!=0)chi2 +=  pow((y0[ibin] -w*y1[ibin])/y0[ibin],2.0)/double(bin_x1-bin_x0+1-(bin_x3 - bin_x2+1));
      if(y0[ibin]!=0)chi2 +=  pow((y0[ibin] -w*y1[ibin])/sqrt(fabs(y0[ibin])),2.0)/double(bin_x1-bin_x0+1-(bin_x3 - bin_x2+1));
      //      cout<<"w "<<w<<" ibin "<<ibin<<" y0 "<<y0[ibin]<<" y1*w "<<w*y1[ibin]<<" chi2 "<<chi2<<endl;
    }

    //    cout<<"bin_x0 "<<bin_x0<<" bin_x1 "<<bin_x1<<" bin_x2 "<<bin_x2<<" bin_x3 "<<bin_x3<<endl;

    if(chi2_min > chi2){
      chi2_min = chi2;
      wmin = w;
      for(int i=0;i<nbins;i++)
	diff[i]=(y0[i] -wmin*y1[i])/sqrt(fabs(y0[i]));
    }
    gchi_L->SetPoint(wi,w,chi2);
  } // for w
  

  for(int i=0;i<nbins;i++)gdiff_L->SetPoint(i,hT_Lsim_s->GetBinCenter(i),diff[i]);
  
  wmin_L = wmin;
  hT_Lsim_s ->Scale(wmin);

  //  hexp_L->Scale(wmin);
  hexp_L->Add(hT_Lsim_s,-1);  // Lambda B.G Estimation
  // Counts # of Lambda B.G.

  NL_bg=0.0;
  NL_bg =hexp_L->Integral(hexp_L->FindBin(-5.0),hexp_L->FindBin(5.0));

  cout<<"wmin_L "<<wmin_L<<" # of Labmda B.G "<<NL_bg<<endl;

}

////////////////////////////////////////////////////////////////////////////
/*
void QFsim::FitQFnnL(){




  int nbins = hnnL_sim->GetXaxis()->GetNbins();
  int nmin  = hnnL_sim->GetXaxis()->GetXmin();
  int nmax  = hnnL_sim->GetXaxis()->GetXmax();

  // Substract Lambda B.G.

  double ratio = NL_bg/(double)hH_nnLsim->GetEntries();
  hH_nnLsim_s->Scale(ratio);
  hexp_c->Add(hH_nnLsim_s,-1.0);
  
  // hexp_c->Add(hexp_L,-1);  
  //  hexp_c->Add(hT_Lsim_s,-1);



  
  double ymax0 = hexp_c -> GetBinContent(hexp_c->GetMaximumBin());
  double ymax1 = hnnL_sim  -> GetBinContent(hnnL_sim->GetMaximumBin());



  
  double y0[nbins],y1[nbins];
  int bin_x0=0;
  int bin_x1=nbins;
  double x0;
  bool max= true;
  for(int ibin = 0;ibin<nbins;ibin++){
    x0 = hnnL_sim -> GetXaxis()->GetBinCenter(ibin);
    if( x0 < 40.0 )bin_x0  = ibin; // Start Chi calc
    if( x0 < 120. )bin_x1  = ibin; // end   Chi calc
    y0[ibin]  =  hexp_c     -> GetBinContent(ibin);
    y1[ibin]  =  hnnL_sim   -> GetBinContent(ibin);
    
  }


  
  double chi2,chi2_min,w,wmin;
  int wmax=100;
  double width = 0.01;
  chi2_min=1e20;
  for(int wi=0; wi<wmax;wi++){
    chi2=0.0;
    w = ymax0/ymax1*(1.0 +(double)(-wmax/2. + wi)*(double)width);
    for(int ibin = bin_x0;ibin<bin_x1;ibin++){
      if(y0[ibin]!=0)chi2 +=  pow((y0[ibin] -w*y1[ibin])/sqrt(fabs(y0[ibin])),2.0)/double(nbins-bin_x0);

    }
    if(chi2_min > chi2){
      chi2_min = chi2;
      wmin = w;
    }
    gchi_nnL->SetPoint(wi,w,chi2);
  } // for w



  hnnL_sim_s->Scale(wmin);
  hmm_s->Scale(wmin);



  double scale_fsi1 = hnnL_sim_s ->GetBinContent(hnnL_sim_s->GetMaximumBin())/hmm_fsi1->GetBinContent(hmm_fsi1->GetMaximumBin());
  hmm_fsi1_s->Scale(scale_fsi1);
  
  double scale_fsi2 = hnnL_sim_s ->GetBinContent(hnnL_sim_s->GetMaximumBin())/hmm_fsi2->GetBinContent(hmm_fsi2->GetMaximumBin());
  hmm_fsi2_s->Scale(scale_fsi2);
  
  double scale_fsi3 = hnnL_sim_s ->GetBinContent(hnnL_sim_s->GetMaximumBin())/hmm_fsi3->GetBinContent(hmm_fsi3->GetMaximumBin());
  hmm_fsi3_s->Scale(scale_fsi3);


  
  hexp_c->Add(hH_nnLsim_s, 1.0);  
  hmm_s->Add(hH_nnLsim_s,1.0);
  hmm_fsi1_s->Add(hH_nnLsim_s,1.0);
  hmm_fsi2_s->Add(hH_nnLsim_s,1.0);
  hmm_fsi3_s->Add(hH_nnLsim_s,1.0);


  hmm_c  = (TH1D*)hmm_s -> Clone();
  hmm_c       ->SetName("hmm_c");
  hmm_c->SetTitle("Simulation w/ B.G. with scaling ; Missing mass [GeV]; Counts/2MeV");    
  hmm_fsi1_c  = (TH1D*)hmm_fsi1_s -> Clone();
  hmm_fsi1_c ->SetName("hmm_fsi1_c");
  hmm_fsi1_c->SetTitle("Simulation w/ B.G. with scaling (model 1 FSI); Missing mass [GeV]; Counts/2MeV");    
  hmm_fsi2_c  = (TH1D*)hmm_fsi2_s -> Clone();
  hmm_fsi2_c ->SetName("hmm_fsi2_c");
  hmm_fsi2_c->SetTitle("Simulation w/ B.G. with scaling (model 2 FSI); Missing mass [GeV]; Counts/2MeV");    
  hmm_fsi3_c  = (TH1D*)hmm_fsi3_s-> Clone();
  hmm_fsi3_c ->SetName("hmm_fsi3_c");
  hmm_fsi3_c->SetTitle("Simulation w/ B.G. with scaling (model 3 FSI); Missing mass [GeV]; Counts/2MeV");
  
  hexp_c     ->Add(hexp_acc,1.0);
  hmm_c      ->Add(hexp_acc,1.0);
  hmm_fsi1_c ->Add(hexp_acc,1.0);
  hmm_fsi2_c ->Add(hexp_acc,1.0);
  hmm_fsi3_c ->Add(hexp_acc,1.0);

  
}

*/


/////////////////////////////////////////////////////////////////////////////

void QFsim::FitQFnnL_new(){




  int nbins = hnnL_sim->GetXaxis()->GetNbins();
  int nmin  = hnnL_sim->GetXaxis()->GetXmin();
  int nmax  = hnnL_sim->GetXaxis()->GetXmax();
  int nevnt = hnnL_sim->GetEntries();

  //------ 10 keV nnL (H run ) hist ------//

  double nbins_Lam_10keV =(double) hH_nnLsim_10keV->GetXaxis()->GetNbins();
  double nbins_Lam       =(double) hH_nnLsim      ->GetXaxis()->GetNbins();



  
  
  // Substract Lambda B.G.

  double ratio = NL_bg/(double)hH_nnLsim->GetEntries();
  hH_nnLsim_s->Scale(ratio);
  hexp_c->Add(hH_nnLsim_s,-1.0);

  // Lambda B.G in 10 keV hist
  double ratio_10keV = NL_bg/(double)hH_nnLsim_10keV->GetEntries();
  hH_nnLsim_10keV_s->Scale(ratio_10keV);
  hexp_10keV_c->Add(hH_nnLsim_10keV_s,-1.0); // w/o Lam. B.G.
  
  
  
  // Scale FSI MM simulation
  double ENum[nvp], scale[nvp],ENum_exp;
  double ENum_10keV[nvp], scale_10keV[nvp],ENum_exp_10keV;



  
  ENum_exp = hexp_c->Integral(hexp_c->GetXaxis()->FindBin(-1000.),hexp_c->GetXaxis()->FindBin(1000.));  // w/o Lambda B.G.
  
  //  ENum[0] =hnnL_sim_s->Integral(hnnL_sim_s->GetXaxis()->FindBin(-1000.),hnnL_sim_s->GetXaxis()->FindBin(1000.));
  
    ENum[0] =hmm->Integral(hmm->GetXaxis()->FindBin(-1000.),hmm->GetXaxis()->FindBin(1000.));

  ENum[1] = hmm_fsi1->Integral(hmm_fsi1->GetXaxis()->FindBin(-1000.),hmm_fsi1->GetXaxis()->FindBin(1000.));

  ENum[2] = hmm_fsi2->Integral(hmm_fsi2->GetXaxis()->FindBin(-1000.),hmm_fsi2->GetXaxis()->FindBin(1000.));

  ENum[3] = hmm_fsi3->Integral(hmm_fsi3->GetXaxis()->FindBin(-1000.),hmm_fsi3->GetXaxis()->FindBin(1000.));  
  


  
  //   ENum_exp_10keV = hexp_10keV_s->Integral(hexp_10keV_s->GetXaxis()->FindBin(-1000.),hexp_10keV_s->GetXaxis()->FindBin(1000.));


  
  ENum_exp_10keV = (double)hexp_10keV_c->GetEntries(); // w/o Lambda B.G.

  
  ENum_10keV[0] =hmm_10keV->Integral(hmm_10keV->GetXaxis()->FindBin(-1000.),hmm_10keV->GetXaxis()->FindBin(1000.));

  ENum_10keV[1] = hmm_fsi1_10keV->Integral(hmm_fsi1_10keV->GetXaxis()->FindBin(-1000.),hmm_fsi1_10keV->GetXaxis()->FindBin(1000.));

  ENum_10keV[2] = hmm_fsi2_10keV->Integral(hmm_fsi2_10keV->GetXaxis()->FindBin(-1000.),hmm_fsi2_10keV->GetXaxis()->FindBin(1000.));

  ENum_10keV[3] = hmm_fsi3_10keV->Integral(hmm_fsi3_10keV->GetXaxis()->FindBin(-1000.),hmm_fsi3_10keV->GetXaxis()->FindBin(1000.));  


  //  cout<<"exp 10 keV "<<ENum_exp_10keV<<" exp "<<ENum_exp<<endl;


  
  scale[0] = ENum_exp/ENum[0];
  scale[1] = ENum_exp/ENum[1];
  scale[2] = ENum_exp/ENum[2];
  scale[3] = ENum_exp/ENum[3];

  hnnL_sim_s->Scale(scale[0]);
  hmm_s     ->Scale(scale[0]);
  hmm_fsi1_s->Scale(scale[1]);
  hmm_fsi2_s->Scale(scale[2]);
  hmm_fsi3_s->Scale(scale[3]);
  
  scale_10keV[0] = ENum_exp_10keV/ENum_10keV[0];
  scale_10keV[1] = ENum_exp_10keV/ENum_10keV[1];
  scale_10keV[2] = ENum_exp_10keV/ENum_10keV[2];
  scale_10keV[3] = ENum_exp_10keV/ENum_10keV[3];


  //  cout<<"scale "<<scale[0]<<" scale (10keV ) "<<scale_10keV[0]<<endl;
  //  hnnL_sim_10keV_s->Scale(scale_10keV[0]);

  
  hmm_10keV_s     ->Scale(scale_10keV[0]);
  hmm_fsi1_10keV_s->Scale(scale_10keV[1]);
  hmm_fsi2_10keV_s->Scale(scale_10keV[2]);
  hmm_fsi3_10keV_s->Scale(scale_10keV[3]);
  

  //  cout<<"scale 10 keV  "<<scale_10keV[0]<<endl;


  // Get Counts each Hist 
  
  
  double ymax0 = hexp_c    -> GetBinContent(hexp_c->GetMaximumBin());
  double ymax1 = hnnL_sim  -> GetBinContent(hnnL_sim->GetMaximumBin());
  //  double ymax1 = hnnL_sim_s  -> GetBinContent(hnnL_sim->GetMaximumBin());



  double yexp[nbins],yfsi[nbins][nvp];

  //  for(int ibin = 0;ibin<nbins;ibin++){
  for(int ibin = 0;ibin<nbins;ibin++){

      yexp[ibin]     =  hexp_c     -> GetBinContent(ibin);    
      yfsi[ibin][0]  =  hnnL_sim_s -> GetBinContent(ibin);
      yfsi[ibin][1]  =  hmm_fsi1_s -> GetBinContent(ibin);
      yfsi[ibin][2]  =  hmm_fsi2_s -> GetBinContent(ibin);
      yfsi[ibin][3]  =  hmm_fsi3_s -> GetBinContent(ibin);
    
  }




  //   Calculation of minimalization of chi2 (Fitting)
  //   wmin : weighting factor
  
  double chi2,chi2_min[nvp],w,wmin[nvp];
  int wmax=100;
  double width = 0.02;
  double diff[nbins][nvp];
  for(int i=0;i<nvp;i++)chi2_min[i]=1e20;
  
  int bin_st  = hexp_c->GetXaxis()->FindBin(0.0);  // Start Chi calc
  int bin_end = hexp_c->GetXaxis()->FindBin(40.0); // End   Chi calc

  
  // w/o FSI hist Fit not enhance reigon MX>40 MeV
  int bin_st_wofsi = hexp_c->GetXaxis()->FindBin(  40.0); // Start Chi calc
  int bin_end_wofsi = hexp_c->GetXaxis()->FindBin(160.0); // End   Chi calc
  
  for(int i=0;i<nvp;i++){ // V-potential
    for(int wi=0; wi<wmax;wi++){ // Set weight Factor
      chi2=0.0;
      w = (1.0 +(double)(-wmax/2. + wi)*(double)width);
      for(int ibin = bin_st;ibin<bin_end;ibin++){

	if(i==0 && (ibin<bin_st_wofsi || bin_end_wofsi<ibin)
	   && yfsi[ibin][i]!=0)chi2 +=  pow((yexp[ibin] -w*yfsi[ibin][i])/sqrt(fabs(yexp[ibin])),2.0)/double(bin_end_wofsi - bin_st_wofsi +1);
	
	if(yfsi[ibin][i]!=0)chi2 +=  pow((yexp[ibin] -w*yfsi[ibin][i])/sqrt(fabs(yexp[ibin])),2.0)/double(bin_end - bin_st +1);
      }
      
      if(chi2_min[i] > chi2){
	chi2_min[i] = chi2;
	wmin[i] = w;
	for(int j=0;j<nbins;j++)
	  diff[j][i] = (yexp[j] -wmin[i]*yfsi[j][i])/sqrt(fabs(yexp[j]));
    }
      gchi_fsi[i]->SetPoint(wi,w,chi2);
    } // end for w
  } // end for i (nvp)
    


    for(int ibin=0;ibin<nbins;ibin++){ 
      gdiff_fsi[1] ->SetPoint(ibin,hmm_fsi1->GetBinCenter(ibin), diff[ibin][1]);
      gdiff_fsi[2] ->SetPoint(ibin,hmm_fsi2->GetBinCenter(ibin), diff[ibin][2]);
      gdiff_fsi[3] ->SetPoint(ibin,hmm_fsi3->GetBinCenter(ibin), diff[ibin][3]);
    }
  
  hnnL_sim_s->Scale(wmin[0]);
  hmm_s->Scale(wmin[0]);
  hmm_10keV_s ->Scale(wmin[0]);

  
  hmm_fsi1_s->Scale(wmin[1]);
  hmm_fsi2_s->Scale(wmin[2]);
  hmm_fsi3_s->Scale(wmin[3]);

  hmm_fsi1_10keV_s->Scale(wmin[1]);
  hmm_fsi2_10keV_s->Scale(wmin[2]);
  hmm_fsi3_10keV_s->Scale(wmin[3]);    

  
  hmm_c  = (TH1D*)hmm_s -> Clone();
  hmm_c       ->SetName("hmm_c");
  hmm_fsi1_c  = (TH1D*)hmm_fsi1_s -> Clone();
  hmm_fsi1_c ->SetName("hmm_fsi1_c");
  hmm_fsi2_c  = (TH1D*)hmm_fsi2_s -> Clone();
  hmm_fsi2_c ->SetName("hmm_fsi2_c");
  hmm_fsi3_c  = (TH1D*)hmm_fsi3_s-> Clone();
  hmm_fsi3_c ->SetName("hmm_fsi3_c");

  
  hmm_10keV_c  = (TH1D*)hmm_10keV_s -> Clone();
  hmm_10keV_c       ->SetName("hmm_10keV_c");
  hmm_fsi1_10keV_c  = (TH1D*)hmm_fsi1_10keV_s -> Clone();
  hmm_fsi1_10keV_c ->SetName("hmm_fsi1_10keV_c");
  hmm_fsi2_10keV_c  = (TH1D*)hmm_fsi2_10keV_s -> Clone();
  hmm_fsi2_10keV_c ->SetName("hmm_fsi2_10keV_c");
  hmm_fsi3_10keV_c  = (TH1D*)hmm_fsi3_10keV_s-> Clone();
  hmm_fsi3_10keV_c ->SetName("hmm_fsi3_10keV_c");
  

  // Add Lambda B.G. 
  hexp_c->Add(hH_nnLsim_s, 1.0);
  hmm_c ->Add(hH_nnLsim_s, 1.0);
  hmm_fsi1_c->Add(hH_nnLsim_s,1.0);
  hmm_fsi2_c->Add(hH_nnLsim_s,1.0);
  hmm_fsi3_c->Add(hH_nnLsim_s,1.0);

  // Add Acc. B.G.
  hexp_c     ->Add(hexp_acc,1.0);
  hmm_c      ->Add(hexp_acc,1.0);
  hmm_fsi1_c ->Add(hexp_acc,1.0);
  hmm_fsi2_c ->Add(hexp_acc,1.0);
  hmm_fsi3_c ->Add(hexp_acc,1.0);


  // Add Lambda B.G. 
  hexp_10keV_c->Add(hH_nnLsim_10keV_s, 1.0);
  hmm_10keV_c ->Add(hH_nnLsim_10keV_s, 1.0);
  hmm_fsi1_10keV_c->Add(hH_nnLsim_10keV_s,1.0);
  hmm_fsi2_10keV_c->Add(hH_nnLsim_10keV_s,1.0);
  hmm_fsi3_10keV_c->Add(hH_nnLsim_10keV_s,1.0);

  // Add Acc. B.G.
  hexp_10keV_c     ->Add(hexp_acc_10keV,1.0);
  hmm_10keV_c      ->Add(hexp_acc_10keV,1.0);
  hmm_fsi1_10keV_c ->Add(hexp_acc_10keV,1.0);
  hmm_fsi2_10keV_c ->Add(hexp_acc_10keV,1.0);
  hmm_fsi3_10keV_c ->Add(hexp_acc_10keV,1.0);  


  //===== Calc Peak Signifiance ==========//



  double sim_obin = (hmm_c->GetBinCenter(2) - hmm_c->GetBinCenter(1));
  double sim_10keV_obin =(hmm_10keV->GetBinCenter(2) - hmm_10keV->GetBinCenter(1));
  
  CalcPS_c(hexp_10keV_c,hmm_10keV_c, 0);  // w/o FSI 
  CalcPS_c(hexp_10keV_c,hmm_fsi1_10keV_c, 1);  // w/ FSI Verma Potential 
  CalcPS_c(hexp_10keV_c,hmm_fsi2_10keV_c, 2);  // w/ FSI Verma Potential
  CalcPS_c(hexp_10keV_c,hmm_fsi3_10keV_c, 3);  // w/ FSI Verma Potential 

 
  
}


////////////////////////////////////////////////////////////////////////////



void QFsim::CalcPS(TH1D* hmm){


  double min_mm = hmm->GetXaxis()->GetXmin();
  double max_mm = hmm->GetXaxis()->GetXmax();
  const int Nbin = (int)(max_mm - max_mm)*30;

  double S[Nbin],N[Nbin],Pval[Nbin],PS[Nbin], MM[Nbin];
  double bin_width = (hmm->GetBinCenter(2) - hmm->GetBinCenter(1));
  
  double width = 1.0;
  double window = width/bin_width;
  double bgr = 2.0; // BackGround Range with sigma
  cout<<"Bin size "<<bin_width<<" MeV"<<endl;
  cout<<"Peak Resolution: "<<window<<" MeV"<<endl;
  
  TF1* f = new TF1("f_pois","TMath::Poisson(x,[0])",0,1000);

  
  for(int i=0;i<Nbin;i++){

    Pval[i] =  PS[i] = 0.0;
    Pval[i] =  PS[i] = 0.0;
    double bg1 = hmm->Integral(i-int((2+bgr)*window),i-int(2*window));
    double bg2 = hmm->Integral(i+int(2*window),i+int((2+bgr)*window));
    double bg = (bg1 + bg2) / (bgr);
    double pk = hmm->Integral(i - int(window), i + int(window));
    double mm = hmm->GetBinCenter(i);

    
    f ->SetParameter(0,bg);
    double pval = f->Integral(pk,1000.,1.0e-12);
    if( pval == 0. ){ pval = 1.; }
    if( i-int(6*window) <= 0    ){ pval = 1.; }
    if( i+int(6*window) >= Nbin ){ pval = 1.; }
    double ps = ROOT::Math::gaussian_quantile_c(pval,1.);
    if( pval>=1 || pval<=0 ){ ps = 0; }
    if( ps < 0 ){ ps = 0; }

    S[i]    = pk;
    N[i]    = bg;
    Pval[i] = pval;
    PS[i]   = ps;
    MM[i]   = mm;
    f->Clear();

    //    g_s   ->SetPoint(i,MM[i],S[i]);
    //    g_n   ->SetPoint(i,MM[i],N[i]);
    //    g_val ->SetPoint(i,MM[i],Pval[i]);
    //    g_ps  ->SetPoint(i,MM[i],PS[i]);
    
  }
    

  
} 



//////////////////////////////////////////////////////////////////////////////

void QFsim::CalcPS_c(TH1D* hexp, TH1D* hsim, int model){


  
  double min_mm = hexp->GetXaxis()->GetXmin();
  double max_mm = hexp->GetXaxis()->GetXmax();
  //  const int Nbin = (int)(max_mm - min_mm)*30;
  double bin_width = (hexp->GetBinCenter(2) - hexp->GetBinCenter(1));
  const int Nbin = (int)( (max_mm - min_mm)/bin_width );
  double S[Nbin],N[Nbin],Pval[Nbin],PS[Nbin], MM[Nbin];
  TF1* f = new TF1("f_pois","TMath::Poisson(x,[0])",0,1000);

  
  //  double window = width/bin_width;

  double max_ps,min_val,mean;
  double peak_sigma = 1.0; //MeV
  
  int jmax =10;
  
  for(int j=0;j<jmax;j++){


    //    double width = 2.0*peak_sigma; //MeV

    min_val = 100.;
    max_ps  = 0.0 ;
    mean    = -100.;
    
    double width = ((double)j*0.1 +1.0) *peak_sigma; //MeV
    double window = width/bin_width;


    for(int i=0;i<Nbin;i++){
      
      Pval[i] =  PS[i] = 0.0;
      Pval[i] =  PS[i] = 0.0;
      
      double bg = hsim -> Integral(i - int(window), i + int(window));
      double pk = hexp -> Integral(i - int(window), i + int(window));
      double mm = hexp -> GetBinCenter(i);

      
      f ->SetParameter(0,bg);
      double pval = f->Integral(pk,1000.,1.0e-12);
      if( pval == 0. ){ pval = 1.; }
      if( i-int(6*window) <= 0    ){ pval = 1.; }
      if( i+int(6*window) >= Nbin ){ pval = 1.; }
      //    double ps = ROOT::Math::gaussian_quantile_c(pval,1.);
      
      double ps = ROOT::Math::gaussian_quantile_c(pval,1.0);
      if( pval>=1 || pval<=0 ){ ps = 0; }
      if( ps < 0 ){ ps = 0; }
      
      S[i]    = pk;
      N[i]    = bg;
      Pval[i] = pval;
      PS[i]   = ps;
      MM[i]   = mm;
      f->Clear();

      if(j==0 ){
      g_s[model]   ->SetPoint(i,MM[i],S[i]);
      g_n[model]   ->SetPoint(i,MM[i],N[i]);
      g_val[model] ->SetPoint(i,MM[i],Pval[i]);
      g_ps[model]  ->SetPoint(i,MM[i],PS[i]);

      }
      
      if(fabs(MM[i])<2.0 && max_ps < PS[i]){
	max_ps = PS[i];
	min_val = Pval[i];
	mean    = MM[i];
      }
    }// END ibin

    g_ps_sig[model]-> SetPoint(j,width,max_ps);
    g_val_sig[model]-> SetPoint(j,width,min_val); 
    g_mm_sig[model]-> SetPoint(j,width,mean); 
  }// END j sigma
  
}



/////////////////////////////////////////////////////////////////////////////


void QFsim::Write(){

  // Histgram


  hexp->Write();
  hexp_c->Write();
  hexp_L->Write();
  hexp_nnL->Write();
  hexp_peak_L->Write();
  hexp_peak->Write();
  hexp_acc->Write();
  hH_Lsim->Write();
  hH_Lsim_s->Write();
  hT_Lsim->Write();
  hT_Lsim_s->Write();

  hH_nnLsim->Write();
  hH_nnLsim_s->Write();
  hT_nnLsim->Write();
  hT_nnLsim_s->Write();  


  hmm     ->Write();
  hmm_fsi1->Write();
  hmm_fsi2->Write();
  hmm_fsi3->Write();
  hmm_s     ->Write();  
  hmm_fsi1_s->Write();
  hmm_fsi2_s->Write();
  hmm_fsi3_s->Write();
  hmm_c     ->Write();
  hmm_fsi1_c->Write();
  hmm_fsi2_c->Write();
  hmm_fsi3_c->Write();  

  // 10 keV 
  hexp_10keV    ->Write();
  hmm_10keV     ->Write();
  hmm_fsi1_10keV->Write();
  hmm_fsi2_10keV->Write();
  hmm_fsi3_10keV->Write();
  hexp_10keV_s    ->Write();
  hmm_10keV_s     ->Write();  
  hmm_fsi1_10keV_s->Write();
  hmm_fsi2_10keV_s->Write();
  hmm_fsi3_10keV_s->Write();
  hexp_10keV_c    ->Write();
  hmm_10keV_c     ->Write();
  hmm_fsi1_10keV_c->Write();
  hmm_fsi2_10keV_c->Write();
  hmm_fsi3_10keV_c->Write();  
  
  // TGraph
  gchi_L->Write();
  gchi_nnL->Write();
  gdiff_L->Write();

  for(int i=0;i<nvp;i++){
    gchi_fsi[i]->Write();
    gdiff_fsi[i]->Write();
    g_s[i]   ->Write();
    g_n[i]   ->Write();
    g_val[i] ->Write();
    g_ps[i]  ->Write();
    g_val_sig[i]->Write();
    g_ps_sig[i]->Write();
    g_mm_sig[i]->Write();
  }
  
  ofr->Write();
  ofr->Close();
}




////////////////////////////////////////////////////////
///////////////////// Nain /////////////////////////////
////////////////////////////////////////////////////////




int main(int argc, char** argv){


  int ch;
  //  string ifname = "../rootfiles/simc/nnL_simc.root";
    string ifname = "./param/test.in";

  extern char *optarg;

  while((ch=getopt(argc,argv,"h:p:s:r:w:f:s:l:TDHeISbcop"))!=-1){
    switch(ch){
            
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;

    case 'h':
      cout<<"-f : input root  filename"<<endl;

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
  QFsim* QF =new QFsim();

  QF -> SetRootList(ifname);
  QF -> SetHist();
  QF -> FitQFL();
  //QF -> FitQFnnL();
  QF -> FitQFnnL_new();
  QF -> Write();
  gSystem->Exit(1);
  theApp->Run(); 
  return 0;

}//end main
