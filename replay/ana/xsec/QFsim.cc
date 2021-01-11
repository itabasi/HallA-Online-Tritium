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
  ofr    = new TFile(rname[3].c_str(),"recreate"); // OutPut Root File
  Tnew = new TTree("T","QF Fitting ");
  cout<<"Out Put Root File "<<rname[3]<<endl;
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
  hexp_c        = (TH1D*)hexp -> Clone();
  hexp_c ->SetName("hexp_c");
  hexp_acc        = (TH1D*)fexp -> Get("h_acc_nnL");
  hexp_acc ->SetName("hexp_acc");
  hexp_nnL   = (TH1D*)fexp -> Get("h_mm_nnL");
  hexp_nnL ->SetName("hexp_nnL");
  hexp_L     = (TH1D*)fexp -> Get("h_peak_L");
  hexp_L ->SetName("hexp_L");
  

  hmm       =(TH1F*)fT_sim->Get("hmm");
  hmm       ->SetName("hmm");
  hmm_fsi1  =(TH1F*)fT_sim->Get("hmm_fsi1");
  hmm_fsi1->SetName("hmm_fsi1");
  hmm_fsi2  =(TH1F*)fT_sim->Get("hmm_fsi2");
  hmm_fsi2->SetName("hmm_fsi2");
  hmm_fsi3  =(TH1F*)fT_sim->Get("hmm_fsi3");
  hmm_fsi3->SetName("hmm_fsi3");


  hmm_s       = (TH1F*)hmm      -> Clone();
  hmm_s       ->SetName("hmm_s");
  hmm_fsi1_s  = (TH1F*)hmm_fsi1 -> Clone();
  hmm_fsi1_s ->SetName("hmm_fsi1_s");
  hmm_fsi2_s  = (TH1F*)hmm_fsi2 -> Clone();
  hmm_fsi2_s ->SetName("hmm_fsi2_s");
  hmm_fsi3_s  = (TH1F*)hmm_fsi3-> Clone();
  hmm_fsi3_s ->SetName("hmm_fsi3_s");



  
  hH_Lsim   = (TH1F*)fH_sim ->Get("hmm_L");
  hH_Lsim -> SetName("hH_Lsim");
  hH_Lsim_s  =(TH1F*)hH_Lsim->Clone();
  hH_Lsim_s -> SetName("hH_Lsim_s");
  hH_nnLsim   = (TH1F*)fH_sim ->Get("hmm_nnL");
  hH_nnLsim  -> SetName("hH_nnLsim");
  hH_nnLsim_s =(TH1F*)hH_nnLsim->Clone();
  hH_nnLsim_s -> SetName("hH_nnLsim_s"); 
  
  hT_Lsim   = (TH1F*)fT_sim ->Get("hmm_L");
  hT_Lsim -> SetName("hT_Lsim");
  hT_Lsim_s =(TH1F*)hT_Lsim->Clone();
  hT_Lsim_s -> SetName("hT_Lsim_s");
  hT_nnLsim   = (TH1F*)fT_sim ->Get("hmm_nnL");
  hT_nnLsim  -> SetName("hT_nnLsim");
  hT_nnLsim_s =(TH1F*)hT_nnLsim->Clone();
  hT_nnLsim_s -> SetName("hT_nnLsim_s"); 

  hnnL_sim   =(TH1F*)fT_sim->Get("hmm");
  hnnL_sim ->SetName("hnnL_sim");
  hnnL_sim_s =(TH1F*)hnnL_sim->Clone();
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
    if( x0 < 0.0 )bin_x0 = ibin;
    if( x0 < 100)bin_x1 = ibin;
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


  hmm_c  = (TH1F*)hmm_s -> Clone();
  hmm_c       ->SetName("hmm_c");
  hmm_fsi1_c  = (TH1F*)hmm_fsi1_s -> Clone();
  hmm_fsi1_c ->SetName("hmm_fsi1_c");
  hmm_fsi2_c  = (TH1F*)hmm_fsi2_s -> Clone();
  hmm_fsi2_c ->SetName("hmm_fsi2_c");
  hmm_fsi3_c  = (TH1F*)hmm_fsi3_s-> Clone();
  hmm_fsi3_c ->SetName("hmm_fsi3_c");
  
  hexp_c     ->Add(hexp_acc,1.0);
  hmm_c      ->Add(hexp_acc,1.0);
  hmm_fsi1_c ->Add(hexp_acc,1.0);
  hmm_fsi2_c ->Add(hexp_acc,1.0);
  hmm_fsi3_c ->Add(hexp_acc,1.0);

  
}


/////////////////////////////////////////////////////////////////////////////

void QFsim::FitQFnnL_new(){




  int nbins = hnnL_sim->GetXaxis()->GetNbins();
  int nmin  = hnnL_sim->GetXaxis()->GetXmin();
  int nmax  = hnnL_sim->GetXaxis()->GetXmax();

  // Substract Lambda B.G.

  double ratio = NL_bg/(double)hH_nnLsim->GetEntries();
  hH_nnLsim_s->Scale(ratio);
  hexp_c->Add(hH_nnLsim_s,-1.0);


  

  // Scale FSI MM simulation
  double ENum[nvp], scale[nvp],ENum_exp;

  ENum_exp = hexp_c->Integral(hexp_c->GetXaxis()->FindBin(-1000.),hexp_c->GetXaxis()->FindBin(1000.));
  
  ENum[0] =hnnL_sim_s->Integral(hnnL_sim_s->GetXaxis()->FindBin(-1000.),hnnL_sim_s->GetXaxis()->FindBin(1000.));

  ENum[1] = hmm_fsi1->Integral(hmm_fsi1->GetXaxis()->FindBin(-1000.),hmm_fsi1->GetXaxis()->FindBin(1000.));

  ENum[2] = hmm_fsi2->Integral(hmm_fsi2->GetXaxis()->FindBin(-1000.),hmm_fsi2->GetXaxis()->FindBin(1000.));

  ENum[3] = hmm_fsi3->Integral(hmm_fsi3->GetXaxis()->FindBin(-1000.),hmm_fsi3->GetXaxis()->FindBin(1000.));  
  

  scale[0] = ENum_exp/ENum[0];
  scale[1] = ENum_exp/ENum[1];
  scale[2] = ENum_exp/ENum[2];
  scale[3] = ENum_exp/ENum[3];

  hnnL_sim_s->Scale(scale[0]);
  hmm_s     ->Scale(scale[0]);
  hmm_fsi1_s->Scale(scale[1]);
  hmm_fsi2_s->Scale(scale[2]);
  hmm_fsi3_s->Scale(scale[3]);
  
  
  

  
  double ymax0 = hexp_c -> GetBinContent(hexp_c->GetMaximumBin());
  double ymax1 = hnnL_sim  -> GetBinContent(hnnL_sim->GetMaximumBin());




  double yexp[nbins],yfsi[nbins][nvp];

  //  for(int ibin = 0;ibin<nbins;ibin++){
  for(int ibin = 0;ibin<nbins;ibin++){

      yexp[ibin]     =  hexp_c     -> GetBinContent(ibin);    
      yfsi[ibin][0]  =  hnnL_sim_s -> GetBinContent(ibin);
      yfsi[ibin][1]  =  hmm_fsi1_s -> GetBinContent(ibin);
      yfsi[ibin][2]  =  hmm_fsi2_s -> GetBinContent(ibin);
      yfsi[ibin][3]  =  hmm_fsi3_s -> GetBinContent(ibin);
    
  }


  
  double chi2,chi2_min[nvp],w,wmin[nvp];
  int wmax=100;
  double width = 0.02;
  double diff[nbins][nvp];
  for(int i=0;i<nvp;i++)chi2_min[i]=1e20;
  
  int bin_st  = hexp_c->GetXaxis()->FindBin(0.0);
  int bin_end = hexp_c->GetXaxis()->FindBin(160.0);
  // w/o FSI hist Fit not enhance reigon MX>40 MeV
  int bin_st_wofsi = hexp_c->GetXaxis()->FindBin(40.0);
  int bin_end_wofsi = hexp_c->GetXaxis()->FindBin(160.0);
  
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
      gdiff_fsi[2] ->SetPoint(ibin,hmm_fsi1->GetBinCenter(ibin), diff[ibin][2]);
      gdiff_fsi[3] ->SetPoint(ibin,hmm_fsi1->GetBinCenter(ibin), diff[ibin][3]);
    }
  
  hnnL_sim_s->Scale(wmin[0]);
  hmm_s->Scale(wmin[0]);

  //  hmm_fsi1_s->Scale(wmin[1]);
  //  hmm_fsi2_s->Scale(wmin[2]);
  //  hmm_fsi3_s->Scale(wmin[3]);  

  hmm_c  = (TH1F*)hmm_s -> Clone();
  hmm_c       ->SetName("hmm_c");
  hmm_fsi1_c  = (TH1F*)hmm_fsi1_s -> Clone();
  hmm_fsi1_c ->SetName("hmm_fsi1_c");
  hmm_fsi2_c  = (TH1F*)hmm_fsi2_s -> Clone();
  hmm_fsi2_c ->SetName("hmm_fsi2_c");
  hmm_fsi3_c  = (TH1F*)hmm_fsi3_s-> Clone();
  hmm_fsi3_c ->SetName("hmm_fsi3_c");


  hexp_c->Add(hH_nnLsim_s, 1.0);
  hmm_c ->Add(hH_nnLsim_s, 1.0);
  hmm_fsi1_c->Add(hH_nnLsim_s,1.0);
  hmm_fsi2_c->Add(hH_nnLsim_s,1.0);
  hmm_fsi3_c->Add(hH_nnLsim_s,1.0);


  
  hexp_c     ->Add(hexp_acc,1.0);
  hmm_c      ->Add(hexp_acc,1.0);
  hmm_fsi1_c ->Add(hexp_acc,1.0);
  hmm_fsi2_c ->Add(hexp_acc,1.0);
  hmm_fsi3_c ->Add(hexp_acc,1.0);

  
}


////////////////////////////////////////////////////////////////////////////

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
  
  // TGraph
  gchi_L->Write();
  gchi_nnL->Write();
  gdiff_L->Write();

  for(int i=0;i<nvp;i++){
    gchi_fsi[i]->Write();
    gdiff_fsi[i]->Write();
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
  //  QF -> FitQFnnL();
  QF -> FitQFnnL_new();
  QF -> Write();
  gSystem->Exit(1);
  theApp->Run(); 
  return 0;
  
}//end main
