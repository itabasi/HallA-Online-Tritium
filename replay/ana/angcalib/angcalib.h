#ifndef angcalib_h
#define angcalib_h 1
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <TChain.h>
#include <TMinuit.h>
#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"

int nite=0;

bool RHRSTrue=false;
bool BreakTrue=true;

extern double calcf2t_4th_2(double* P,
			    double xf, double xpf, 
			    double yf, double fpf,
			    double z);
extern double calcf2t_zt(double* P, 
			 double xf, double xpf,
			 double yf, double ypf);

extern double tune(double* pa, int j, int angflag);
extern void fcn1(int &nPar, double* /*grad*/, 
		 double &fval, double* param, int /*iflag*/);
extern void fcn2(int &nPar, double* /*grad*/, 
		 double &fval, double* param, int /*iflag*/);




//====================================================================//
//=========================== Angcalib Class =========================//
//===================================================================//
class angcalib
{
 public:
  angcalib();
  ~angcalib();

 public:
  void HolePosi(bool rarm);
  void Mxpt(string matrix_xp);  
  void Mypt(string matrix_yp);  
  void Scale_corr(string offset_file);  
  void EventSelect(bool rarm);
  void Tuning(string ofMTPname);
  void SetRoot(string ifname);
  //  void SetRunList(string ifname);
  void SetBranch(string ifname, bool rarm);
  void NewBranch(string ofname, bool rarm);
  void MakeHist();
  //  void MTtuning(string ofMTPname, string matrix_name);
  void Fill(bool rarm);
  void Write();
  void Fill_tuned(bool rarm);
  void Close_tree();
  void Draw();

  bool Xp_flag;
  bool Yp_flag;
  //-------- SetBranch ------------//
  TFile* f2;
  TTree* t2;

  double runnum;
  double hallap;
  Double_t trig5;
  Double_t trig4;
  Double_t trig1; 
  int Nrvz,Nlvz;
  int NRvz,NLvz;
  double a1, a2;
  double r_s2_t_pads[max];
  double l_s2_t_pads[max];
  double r_s2_nthit;
  double l_s2_nthit;
  double r_th_fp[max];
  double l_th_fp[max];
  double r_ph_fp[max];
  double l_ph_fp[max];
  double l_x_fp[max];
  double r_x_fp[max];
  double l_y_fp[max];
  double r_y_fp[max];
  double ztR_opt[max];
  double lx[max];
  double ly[max];
  double rx[max];
  double ry[max];
  double ps_asum,sh_asum;
  double gs_asum;
  double Zt[max];
  double Rvz[max];
  double Lvz[max];
  double Rth[max], Rph[max];
  double Lth[max], Lph[max];

  //-------- NewBranch ------------// 
  TFile* fnew;
  TTree* tnew;
  
  double Xpt[max], Ypt[max];
  double Xpt_tuned[max], Ypt_tuned[max];
  double Xpt_init,Ypt_init;
  double Xt, Yt;
  double Zt_tuned[max];
  double XFP[max],  YFP[max];
  double XpFP[max], YpFP[max];
  double ss_x,ss_y;
  int tflag;
  int zfoil;
  int tuning_num;
  //------ Make Hist -----------//
  TH2F* h1;
  TH2F* h2[nfoil];
  TH2F* h2_new[nfoil];
  TH2F* h2_;
  TH2F* h3[nfoil];
  TH2F* h3_;
  TH2F* h3_a;
  TH2F* h3_b;
  TH2F* h3_c;
  TH1F* hz[nfoil];
  TH1F* hz_all;
  TH1F* hph_cut; 
  TH1F* hssy_cut; 
  TH1F* hth;
  TH1F* hth_c;  
  char tempc[500];
  char tempc2[500];
  const double l0 = 100.3;
  double dth[nfoil];
  TGraphErrors* gchi_xp=new TGraphErrors();
  TGraphErrors* gchi_yp=new TGraphErrors();
  //-------- Hole Posi --------//
  int nhole;
  TMarker* mark[nsshole];
  TMarker* mark_real[nsshole][nfoil];  
  TMarker* mark_ang[nsshole];
 //-------- Mxpt && Mypt------------//
  //  double Pxpt[nParamT];  
  //  double Pypt[nParamT];


  //-------- Scale Correction ----//
  double offs_xp[nfoil], offs_yp[nfoil];
  double scal_xp[nfoil], scal_yp[nfoil];
  double temp;
  
 
  //------ Event Selection ---------//
  int ent;
  div_t d;
  int holeg_temp;
  int foilg_temp;
  bool holethrough;
  bool filled;
  bool rtrig;
  bool ltrig;
  bool rarm;
  bool gs_cut;
  double ssx, ssy;
  //------- Tuning -------------//
  bool offset_flag[nfoil];
  //  double x[nite];
  double chi_sq1[100], chi_sq2[100]; // A chi square for each tuning

  //-------- Fill -------------//
  int evshift = 30000;

 //------- Draw ---------------//
  TCanvas* c1;
  TCanvas* c0;
  
};

//=============================================//
//=========== anglecalib ======================//
//=============================================//
angcalib::angcalib(){
  gROOT->SetStyle("Plain");
  //  gStyle->SetOptStat(0);
  Xp_flag=true;
  Yp_flag=true;
  cout<<"=========== Start Angle Calibration ==========="<<endl;
  cout<<" Xp tuning : "<<Xp_flag<<endl;
  cout<<" Yp tuning : "<<Yp_flag<<endl;
  cout<<" Matrix order : "<<nn<<" z "<<nnz<<" # Parameters "<<nParamT<<endl;
  cout<<endl;
};

  angcalib::~angcalib(){};



//========= SetBranch ========================//
void angcalib::SetBranch(string ifname, bool rarm){

  
  f2 = new TFile(ifname.c_str());
  t2 = (TTree*)f2->Get("T");
  
  t2->SetBranchAddress("fEvtHdr.fRun", &runnum);
  t2->SetBranchAddress("HALLA_p", &hallap);
  t2->SetBranchAddress("DR.T1", &trig1);
  t2->SetBranchAddress("DR.T4", &trig4);
  t2->SetBranchAddress("R.ps.asum_c", &ps_asum);
  t2->SetBranchAddress("R.sh.asum_c", &sh_asum);  
  t2->SetBranchAddress("DR.T5", &trig5);
  t2->SetBranchAddress("R.a1.asum_c", &a1);
  t2->SetBranchAddress("R.a2.asum_c", &a2);  
  t2->SetBranchAddress("Ndata.R.tr.vz", &NRvz);
  t2->SetBranchAddress("Ndata.L.tr.vz", &NLvz);
  t2->SetBranchAddress("L.cer.asum_c",&gs_asum);  
  t2->SetBranchAddress("R.tr.x",   &r_x_fp);
  t2->SetBranchAddress("L.tr.x",   &l_x_fp);
  t2->SetBranchAddress("R.tr.y",   &r_y_fp);
  t2->SetBranchAddress("L.tr.y",   &l_y_fp);
  t2->SetBranchAddress("R.tr.th",  &r_th_fp);
  t2->SetBranchAddress("L.tr.th",  &l_th_fp);
  t2->SetBranchAddress("R.tr.ph",  &r_ph_fp);
  t2->SetBranchAddress("L.tr.ph",  &l_ph_fp);

  //------ At Target ------------//
  t2->SetBranchAddress("R.tr.vx",   &rx);
  t2->SetBranchAddress("L.tr.vx",   &lx);
  t2->SetBranchAddress("R.tr.vy",   &ry);
  t2->SetBranchAddress("L.tr.vy",   &ly);
  t2->SetBranchAddress("R.tr.vz", &Rvz);
  t2->SetBranchAddress("L.tr.vz", &Lvz);
  t2->SetBranchAddress("R.tr.tg_th", &Rth);
  t2->SetBranchAddress("R.tr.tg_ph", &Rph);
  t2->SetBranchAddress("L.tr.tg_th", &Lth);
  t2->SetBranchAddress("L.tr.tg_ph", &Lph);
  
  if(rarm==true){
   t2->SetBranchAddress("R.tr.vz_opt",ztR_opt);
   t2->SetBranchAddress("R.tr.vz_tuned",Zt);
   //  t2->SetBranchAddress("R.tr.vz",Zt); //not tuned 
  } else{
   t2->SetBranchAddress("L.tr.vz_opt",ztR_opt);
   //    t2->SetBranchAddress("L.tr.vz",Zt);
   t2->SetBranchAddress("L.tr.vz_tuned",Zt);
  }

  
  

};

//=========== SetNewBranch ==================//
void angcalib::NewBranch(string ofname, bool rarm){


  fnew = new TFile(ofname.c_str(),"recreate");

  if(rarm==true) tnew = new TTree("T","For angle calibration (RHRS)");
    else tnew = new TTree("T","For angle calibration (LHRS)");
 
  tnew->Branch("fEvtHdr.fRun", &runnum,"fEvtHdr.fRun/D"    );
  tnew->Branch("HALLA_p", &hallap,"HALLA_p/D" );
  tnew->Branch("DR.T1", &trig1, "DR.T1/D"    );
  tnew->Branch("DR.T4", &trig4, "DR.T4/D"   );
  tnew->Branch("DR.T5", &trig5,  "DR.T5/D"  );
  tnew->Branch("R.ps.asum_c", &ps_asum,"R.ps.asum_c/D");
  tnew->Branch("R.sh.asum_c", &sh_asum,"R.sh,asum_c/D");
  tnew->Branch("R.a1.asum_c", &a1,"R.a1.asum_c/D");
  tnew->Branch("R.a2.asum_c", &a2,"R.a2.asum_c/D");
  tnew->Branch("L.cer.asum_c",&gs_asum,"L.cer.asum_c/D");
  tnew->Branch("R.tr.x",   r_x_fp, "R.tr.x[100]/D"  );
  tnew->Branch("L.tr.x",   l_x_fp, "L.tr.x[100]/D"  );
  tnew->Branch("R.tr.y",   r_y_fp, "R.tr.y[100]/D"  );
  tnew->Branch("L.tr.y",   l_y_fp, "L.tr.y[100]/D"  );
  tnew->Branch("R.tr.th",  r_th_fp,"R.tr.th[100]/D" );
  tnew->Branch("L.tr.th",  l_th_fp,"L.tr.th[100]/D" );
  tnew->Branch("R.tr.ph",  r_ph_fp,"R.tr.ph[100]/D" );
  tnew->Branch("L.tr.ph",  l_ph_fp,"L.tr.ph[100]/D" );  

  //----- At Target ---------------//
  tnew->Branch("R.tr.vx",   &rx,"R.tr.vx[100]/D");
  tnew->Branch("L.tr.vx",   &lx,"L.tr.vx[100]/D");
  tnew->Branch("R.tr.vy",   &ry,"R.tr.vy[100]/D");
  tnew->Branch("L.tr.vy",   &ly,"L.tr.vy[100]/D");
  tnew->Branch("Ndata.R.tr.vz", &NRvz, "Ndata.R.tr.vz/I");
  tnew->Branch("Ndata.L.tr.vz", &NLvz, "Ndata.L.tr.vz/I");
  tnew->Branch("R.tr.vz", Rvz, "R.tr.vz[100]/D");
  tnew->Branch("L.tr.vz", Lvz, "L.tr.vz[100]/D");
  tnew->Branch("R.tr.tg_th", Rth,"R.tr.tg_th[100]/D");
  tnew->Branch("R.tr.tg_ph", Rph,"R.tr.th_ph[100]/D");
  tnew->Branch("L.tr.tg_th", Lth,"L.tr.tg_th[100]/D");
  tnew->Branch("L.tr.tg_ph", Lph,"L.tr.tg_ph[100]/D");
  tnew->Branch("ss_x",&ss_x,"ss_x/D " );
  tnew->Branch("ss_y",&ss_y,"ss_y/D " );
  tnew->Branch("tflag",&tflag,"tflag/I " );
  tnew->Branch("zfoil",&zfoil,"zfoil/I " );    
  if(rarm==true){
    tnew->Branch("R.tr.vz_opt",ztR_opt, "R.tr.vz_opt[100]/D" );
    tnew->Branch("R.tr.tg_th_opt",Xpt,"R.tr.tg_th_opt[100]/D " );
    tnew->Branch("R.tr.tg_ph_opt",Ypt,"R.tr.tg_ph_opt[100]/D" );
    tnew->Branch("R.tr.vz_tuned",Zt,"R.tr.vz_tuned[100]/D" );
    tnew->Branch("R.tr.tg_th_tuned",Xpt_tuned,"R.tr.tg_th_tuned[100]/D " );
    tnew->Branch("R.tr.tg_ph_tuned",Ypt_tuned,"R.tr.tg_ph_tuned[100]/D " );
    

  }else{
    tnew->Branch("L.tr.vz_opt",ztR_opt, "L.tr.vz_opt[100]/D" );
    tnew->Branch("L.tr.th_tg_opt",Xpt,"L.tr.tg_th_opt[100]/D " );
    tnew->Branch("L.tr.tg_ph_opt",Ypt,"L.tr.tg_ph_opt[100]/D" );
    tnew->Branch("L.tr.vz_tuned",Zt,"L.tr.vz_tuned[100]/D" );
    tnew->Branch("L.tr.tg_th_tuned",Xpt_tuned,"L.tr.tg_th_tuned[100]/D " );
    tnew->Branch("L.tr.tg_ph_tuned",Ypt_tuned,"L.tr.tg_ph_tuned[100]/D " );    
  }



  
};

//========= Make Hist =======================//
void angcalib::MakeHist(){
  
  h1= new TH2F("h1","",100,-5.0,5.0,100,-8.0,8.0);
  //    h1= new TH2F("h1","",100,-5.0,5.0,100,-8.0,8.0); 
  h1->GetXaxis()->SetTitle("y (rad)");
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle("x");
  
  for(int i=0;i<nfoil;i++){
    sprintf(tempc,"h2_%d",i);
    h2[i] = new TH2F(tempc, tempc,
		     h1->GetXaxis()->GetNbins(),
		     h1->GetXaxis()->GetXmin(),
		     h1->GetXaxis()->GetXmax(),
		     h1->GetYaxis()->GetNbins(),
		     h1->GetYaxis()->GetXmin(),
		     h1->GetYaxis()->GetXmax());
    h2[i]->GetXaxis()->SetTitle("y_SS (cm)");
    h2[i]->GetYaxis()->SetTitle("x_SS (cm)");
    l[i] = 0;
    l[i] = sqrt(pow(l0,2.0) + pow(fcent_real[i]*100.,2.0) -2.0*l0*fcent_real[i]*100.*cos(hrs_ang));
    dth[i] = asin(l0/l[i]*sin(hrs_ang)) - hrs_ang;
    projectf[i] = cos( dth[i] );

    
    sprintf(tempc,"h2_new_%d",i);
    h2_new[i] = (TH2F*)h2[i]->Clone(tempc);
    sprintf(tempc,"h3_%d",i);
    h3[i] = (TH2F*)h2[i]->Clone(tempc);
    h3[i]->SetMarkerColor(i+1);

    hz[i]=new TH1F(Form("hz_%d",i),"",400,-0.2,0.2);
    hz[i]->SetFillColor(i+2);
    hz[i]->SetFillStyle(3001);
  }
    hz_all=new TH1F("hz_all","",400,-0.2,0.2);  
  h2_ = (TH2F*)h2[0]->Clone("h2_");
  h2_->SetTitle("SS pattern before event selection");
  h3_ = (TH2F*)h3[0]->Clone("h3_");
  h3_->SetTitle("SS pattern after event selection");
  h3_a = (TH2F*)h3[0]->Clone("h3_a");
  h3_a->SetTitle("SS pattern w/o   matrix tuning");  
  h3_b = (TH2F*)h3[0]->Clone("h3_b");
  h3_b->SetTitle("SS pattern input matrix tuning");  
  h3_c = (TH2F*)h3[0]->Clone("h3_c");
  h3_c->SetTitle("SS pattern after matrix tuning");  
  gchi_xp->SetFillColor(2);
  gchi_xp->SetMarkerStyle(20);
  gchi_xp->SetMarkerColor(2);
  gchi_xp->SetFillStyle(3005);
  gchi_xp->SetMarkerSize(1.0);

  gchi_yp->SetFillColor(2);
  gchi_yp->SetMarkerStyle(20);
  gchi_yp->SetMarkerColor(2);
  gchi_yp->SetFillStyle(3005);
  gchi_yp->SetMarkerSize(1.0);

  hth=new TH1F("hth","",1000,-0.1,0.1);
  hth_c=new TH1F("hth_c","",1000,-0.1,0.1);  

  hph_cut = new TH1F("hph_cut","",1000,-0.1,0.1);
  hssy_cut = new TH1F("hssy_cut","",1000,-5.0,5.0);
  
};



//========== HolePosi =======================//
void angcalib::HolePosi(bool rarm){

  if(rarm)RHRSTrue=true;
  else RHRSTrue=false;

  cout<<"RHRS_Flag "<<RHRSTrue<<endl;
  
  nhole = 0;
   double ssy_cent_real[nrow];
   double ssx_cent_real[ncol];

   //===== Positon offset  =======//

   for(int k=0;k<nsshole;k++){
     for(int i=0;i<nfoil;i++){
     ssx_off[k][i]=0.0;
     ssy_off[k][i]=0.0;


     /*
     ssy_off[13][4]=0.7;     
     ssy_off[15][4]=0.6;     
     ssy_off[17][4]=0.6;     
     ssy_off[19][4]=0.6;     
     ssy_off[21][4]=0.7;
     ssy_off[24][4]=0.2;     
     ssy_off[26][4]=0.2;     
     ssy_off[28][4]=0.2;     
     ssy_off[30][4]=0.2;     
     ssy_off[32][4]=0.2;
     */

     ssy_off[23][3]= -0.25;
     ssy_off[25][3]= -0.25;
     ssy_off[27][3]= -0.25;
     ssy_off[29][3]= -0.25;
     ssy_off[31][3]= -0.25;
     ssy_off[34][3]= -0.5;
     ssy_off[35][3]= -0.35;
     ssy_off[36][3]= -0.5;
     ssy_off[37][3]= -0.35;
     ssy_off[38][3]= -0.5;
     ssy_off[39][3]= -0.35;
     ssy_off[40][3]= -0.5;
     ssy_off[41][3]= -0.35;
     ssy_off[42][3]= -0.5;
     ssy_off[45][3]= -0.75;
     ssy_off[46][3]= -0.4;
     ssy_off[47][3]= -0.75;
     ssy_off[48][3]= -0.4;
     ssy_off[49][3]= -0.75;
     ssy_off[50][3]= -0.4;
     ssy_off[51][3]= -0.75;
     ssy_off[52][3]= -0.4;
     ssy_off[53][3]= -0.75;


     
     ssx_off[56][4]= -0.75;
   
     ssy_off[36][4]= -0.2;     
     ssy_off[38][4]= -0.2;     
     ssy_off[46][4]= -0.2;
     ssy_off[47][4]= -0.2;     
     ssy_off[48][4]= -0.2;     
     ssy_off[49][4]= -0.2;     
     ssy_off[50][4]= -0.2;
     ssy_off[51][4]= -0.2;     
     ssy_off[52][4]= -0.2;
     ssy_off[53][4]= -0.2;
     ssy_off[56][4]= -0.6;
     ssy_off[57][4]= -0.4;     
     ssy_off[58][4]= -0.65;
     ssy_off[59][4]= -0.4;
     ssy_off[60][4]= -0.65; 
     ssy_off[61][4]= -0.4;
     ssy_off[62][4]= -0.65;
     ssy_off[63][4]= -0.4;
     ssy_off[70][4]= -0.9;
     ssy_off[72][4]= -0.9;
	  
     ssy_off[57][5]= -0.1;     
     ssy_off[59][5]= -0.1;     
     ssy_off[61][5]= -0.1;     
     ssy_off[63][5]= -0.1;     
     ssy_off[65][5]= -0.1;

    
     ssy_off[56][5]= -0.4;     
     ssy_off[58][5]= -0.4;     
     ssy_off[60][5]= -0.4;     
     ssy_off[62][5]= -0.4;     
     ssy_off[64][5]= -0.4;
     ssy_off[70][5]= -0.65;     
     ssy_off[72][5]= -0.65;

     
     ssy_off[56][6]= -0.4;
     ssy_off[57][6]= -0.2;
     ssy_off[58][6]= -0.4;
     ssy_off[59][6]= -0.2;
     ssy_off[60][6]= -0.4;
     ssy_off[61][6]= -0.2;
     ssy_off[62][6]= -0.4;
     ssy_off[63][6]= -0.2;


     ssy_off[1][7]=  0.2;
     ssy_off[3][7]=  0.2;
     ssy_off[5][7]=  0.2;
     ssy_off[7][7]=  0.2;
     ssy_off[9][7]=  0.2;
     ssy_off[57][7]= -0.35;
     ssy_off[59][7]= -0.35;
     ssy_off[61][7]= -0.35;
     ssy_off[63][7]= -0.35;

     ssx_off[20][7]=  0.5;
     ssx_off[31][7]=  0.5;
     ssx_off[42][7]=  0.5;
     ssx_off[53][7]=  0.5;
     ssx_off[64][7]=  0.5;

     //     ssy_off[k][4]  = 0.0;
     
     /*
     ssx_off[k][5]  = -0.75;
     ssy_off[k][5]  =  0.3;
     ssx_off[45][5] = -1.5;
     ssy_off[45][5] =  0.1;
     ssx_off[57][5] = -1.2;
     ssy_off[59][5] =  0.2;
     ssy_off[61][5] =  0.2;
     ssy_off[63][5] =  0.2;


     
     ssx_off[k][6]  = -0.75;
     ssy_off[k][6]  = -0.25;     
     ssx_off[45][6] = -1.25;
     ssy_off[45][6] = -0.5;
     ssy_off[47][6] = -0.4;
     ssy_off[49][6] = -0.2;
     ssy_off[51][6] = -0.2;
     ssy_off[53][6] = -0.2;
     ssx_off[53][6] = -0.3;
     
     ssy_off[57][6] = -0.4;
     ssy_off[58][6] = -0.5;     
     ssy_off[59][6] = -0.4;
     ssy_off[60][6] = -0.5;
     ssy_off[61][6] = -0.4;
     ssy_off[62][6] = -0.5;
     ssy_off[63][6] = -0.4;
     ssy_off[68][6] = -0.5;
     ssy_off[70][6] = -0.5;
     ssy_off[72][6] = -0.5;



     
     ssx_off[k][7] = -0.75;
     ssy_off[k][7] = -0.75;

     ssx_off[45][7] = -1.0;
     ssy_off[45][7] = -1.0;          

     
     ssy_off[47][7] = -1.0;
     ssy_off[49][7] = -1.0;          
     ssy_off[56][7] = -1.0;          
     ssy_off[57][7] = -1.0;     
     ssy_off[58][7] = -1.0;          
     ssy_off[59][7] = -1.0;     
     ssy_off[60][7] = -1.0;          
     ssy_off[61][7] = -1.0;     
     ssy_off[62][7] = -1.0;          
     ssy_off[63][7] = -1.0;     
     ssy_off[70][7] = -1.0;          
     ssy_off[72][7] = -1.0;          
     */
     


     //     ssy_off[k][5]  = + 0.3;
     //     ssy_off[45][5] =   0.0;
     //     ssy_off[47][5] =   0.0;
     //     ssy_off[49][5] =   0.0;     
     //     ssy_off[35][5] =   0.5;
     //     ssy_off[37][5] =   0.5;          
     //     ssy_off[39][5] =   0.5;
     //     ssy_off[41][5] =   0.5;


     
     //     ssy_off[11][4] =  0.6;
     //     ssy_off[13][4] =  0.6;
     //     ssy_off[15][4] =  0.6;
     //     ssy_off[17][4] =  0.6;
     //     ssy_off[19][4] =  0.6;

     /*
     
     ssy_off[11][4] =  0.6;
     ssy_off[13][4] =  0.6;
     ssy_off[15][4] =  0.6;
     ssy_off[17][4] =  0.6;
     ssy_off[19][4] =  0.6;
     ssy_off[60][4] = -0.3;


     ssy_off[11][5] =  0.6;
     ssy_off[13][5] =  0.6;
     ssy_off[15][5] =  0.6;
     ssy_off[17][5] =  0.6;
     ssy_off[19][5] =  0.6;

     */

     
     //     ssy_off[12][5] =  0.5;
     //     ssy_off[14][5] =  0.5;
     //     ssy_off[16][5] =  0.5;
     //     ssy_off[18][5] =  0.5;
     //     ssy_off[20][5] =  0.5;

     


     //     ssy_off[64][6]= -0.5;
     //     ssy_off[68][6]= -0.5;
     //     ssy_off[70][6]= -0.15;
     //     ssy_off[72][6]= -0.15;
     //     ssy_off[74][6]= -0.2;			
       

     //     ssy_off[k][7]=-0.3;  
     

 
     //       ssy_off[58][4] = -0.35;	
     //       ssy_off[60][4] = -0.35;
     //       ssy_off[12][4] = 0.4;
     //       ssy_off[14][4] = 0.4;			
     //       ssy_off[16][4] = 0.4;
     //       ssy_off[18][4] = 0.4;
     //       ssy_off[20][4] = 0.4;

     /*       
       ssy_off[12][5]= 0.5;
       ssy_off[14][5]= 0.5;
       ssy_off[16][5]= 0.5;
       ssy_off[18][5]= 0.5;
       ssy_off[20][5]= 0.5;	
       ssy_off[58][5]= -0.2;
       ssy_off[60][5]= -0.2;
       
       
       ssy_off[58][6]= -0.25;
       ssy_off[60][6]= 0.0;
       ssy_off[64][6]= -0.5;
       ssy_off[68][6]= -0.5;
       ssy_off[70][6]= -0.5;
       ssy_off[72][6]= -0.5;
       ssy_off[74][6]= -0.7;			
       //     ssy_off[k][6]=-0.5;
       */
     
     
     //	ssy_off[12][5] =  0.5;
     //	ssy_off[14][5]   =  0.5;
     //	ssy_off[16][5] =  0.5;
     //	ssy_off[18][5] =  0.5;
     //	ssy_off[20][5] =  0.5;
     
     /*
       ssy_off[34][5] =  -0.15;
       ssy_off[36][5] =  -0.15;		
       ssy_off[38][5] =  -0.15;
       ssy_off[40][5] =  -0.15;	
       ssy_off[42][5] =  -0.15;	
       
       ssy_off[35][5] =  -0.25;
       ssy_off[37][5] =  -0.25;		
       ssy_off[39][5] =  -0.25;
       ssy_off[41][5] =  -0.25;	
       
       ssy_off[23][5] =  -0.25;
       ssy_off[25][5] =  -0.25;		
       ssy_off[27][5] =  -0.25;
       ssy_off[29][5] =  -0.25;	
       ssy_off[31][5] =  -0.25;
       
       ssy_off[24][5] =  -0.2;
       ssy_off[26][5] =  -0.2;		
       ssy_off[28][5] =  -0.2;
       ssy_off[30][5] =  -0.2;	
       ssy_off[32][5] =  -0.2;	
     */
     
     
     
     TFlag[k][i]=false;
     
     //     if((11 <= k && k <= 76) && (k % 11 != 0 &&(k-10) % 11 != 0)
     //     	&& (4 <= i && i <= 7)){TFlag[k][i]=true;}
     //	   && (i==5) ){TFlag[k][i]=true;}

         if((1 <= k && k <= 76) && (k % 11 != 0)
	    //	    && (i==5) ){TFlag[k][i]=true;}
	    && (3 <= i && i <= 7)){TFlag[k][i]=true;}

	 //     TFlag[11][i]=false;
	 //     TFlag[13][i]=false;
	 //     TFlag[15][i]=false;     
	 //     TFlag[17][i]=false;
	 //     TFlag[19][i]=false;

	 TFlag[0][i]=false;
	 TFlag[2][i]=false;
	 TFlag[4][i]=false;
	 TFlag[6][i]=false;
	 TFlag[8][i]=false;
	 TFlag[10][i]=false;
	 TFlag[32][i]=false;
	 TFlag[43][i]=false;
	 TFlag[54][i]=false;
	 TFlag[65][i]=false;
	 TFlag[74][i]=false;
	 TFlag[76][i]=false;
	 TFlag[64][i]=false;          
	 TFlag[67][i]=false;
	 TFlag[69][i]=false;
	 TFlag[71][i]=false;     
	 TFlag[73][i]=false;
	 TFlag[75][i]=false;

	 //     TFlag[11][4]=true;
	 //     TFlag[13][4]=true;
	 //     TFlag[15][4]=true;     
	 //     TFlag[17][4]=true;
	 //     TFlag[19][4]=true;
	 

    for(int j=0;j<=21;j++)TFlag[j][3]=false;	 
    for(int j=57;j<=nhole;j++)TFlag[j][3]=false;	 
    TFlag[56][3]=false;
    TFlag[57][3]=false;
    TFlag[58][3]=false;
    TFlag[59][3]=false;
    TFlag[60][3]=false;
    TFlag[61][3]=false;	 
    TFlag[62][3]=false;
    TFlag[63][3]=false;
    TFlag[64][3]=false;
    TFlag[65][3]=false;
    TFlag[66][3]=false;
    TFlag[67][3]=false;
    TFlag[68][3]=false;
    TFlag[69][3]=false;	 
    TFlag[70][3]=false;
    TFlag[71][3]=false;
    TFlag[72][3]=false;
     
    
     
     TFlag[1][4]=false;
     TFlag[3][4]=false;
     TFlag[5][4]=false;
     TFlag[7][4]=false;
     TFlag[9][4]=false;     
     TFlag[11][4]=false;
     TFlag[12][4]=false;
     TFlag[13][4]=false;
     TFlag[14][4]=false;
     TFlag[15][4]=false;     
     TFlag[16][4]=false;
     TFlag[17][4]=false;
     TFlag[18][4]=false;
     TFlag[19][4]=false;
     TFlag[20][4]=false;
     TFlag[32][4]=false;
     TFlag[43][4]=false;
     TFlag[54][4]=false;
     TFlag[65][4]=false;
     //     TFlag[56][4]=false;
     //     TFlag[58][4]=false;     
     //     TFlag[62][4]=false;          
     TFlag[68][4]=false;
     //     TFlag[70][4]=false;
     //     TFlag[72][4]=false;
     TFlag[74][4]=false;     


	 
     TFlag[1][5]=false;
     TFlag[3][5]=false;
     TFlag[5][5]=false;
     TFlag[7][5]=false;
     TFlag[9][5]=false;     
     TFlag[11][5]=false;
     TFlag[13][5]=false;
     TFlag[15][5]=false;     
     TFlag[17][5]=false;
     TFlag[19][5]=false;
     TFlag[32][5]=false;
     TFlag[43][5]=false;
     TFlag[54][5]=false;
     TFlag[65][5]=false;
     //     TFlag[56][5]=false;     
     //     TFlag[58][5]=false;
     //     TFlag[60][5]=false;     
     //     TFlag[62][5]=false;     
     TFlag[68][5]=false;
     //     TFlag[70][5]=false;
     //     TFlag[72][5]=false;
     TFlag[74][5]=false;     


     
     
     TFlag[1][6]=false;
     TFlag[3][6]=false;
     TFlag[5][6]=false;
     TFlag[7][6]=false;
     TFlag[9][6]=false;     
     //     TFlag[11][6]=false;
     //     TFlag[13][6]=false;
     //     TFlag[15][6]=false;     
     //     TFlag[17][6]=false;
     //     TFlag[19][6]=false;
     TFlag[68][6]=false;
     TFlag[70][6]=false;
     TFlag[72][6]=false;
     TFlag[74][6]=false;     
     TFlag[76][6]=false;     

     
     
     //     TFlag[11][7]=false;
     //     TFlag[13][7]=false;
     //     TFlag[15][7]=false;     
     //     TFlag[17][7]=false;
     //     TFlag[19][7]=false;
     TFlag[56][7]=false;
     TFlag[58][7]=false;
     TFlag[60][7]=false;
     TFlag[62][7]=false;     
     TFlag[64][7]=false;
     TFlag[68][7]=false;
     TFlag[70][7]=false;
     TFlag[72][7]=false;
     TFlag[74][7]=false;     
     TFlag[76][7]=false;     

     

     




     
         
     /*
     TFlag[12][6]=false;
     TFlag[14][6]=false;
     TFlag[16][6]=false;
     TFlag[18][6]=false;     
     TFlag[20][6]=false;

     
     TFlag[12][7]=false;
     TFlag[14][7]=false;
     TFlag[16][7]=false;
     TFlag[18][7]=false;     
     TFlag[20][7]=false;
     */
     
     /*
     TFlag[24][7]=false;
     TFlag[26][7]=false;
     TFlag[28][7]=false;
     TFlag[30][7]=false;     
     TFlag[32][7]=false;     
     TFlag[68][7]=false;
     */
     
   
     }
   }

   
   //========= Sieve Slit offset tuning ========//
   //   refy_real[38]=0.5;

   
   
   for(int i=0; i<ncol ; i++){
    for(int j=0; j<nrow; j++){
      ssy_cent_real[i] = -3.0*step + step*i;
      if(j%2==0)ssy_cent_real[i] = ssy_cent_real[i] - step/2.0;
      ssx_cent_real[j] = 5.0*step - step*j;
      refx[nhole] = ssx_cent_real[j];
      refy[nhole] = ssy_cent_real[i];
      mark[nhole] = new TMarker(refy[nhole],refx[nhole],28);
      mark[nhole]->SetMarkerColor(1);

      for(int k=0;k<nfoil;k++){

      refx_real[nhole][k] =refx[nhole]+ssx_off[nhole][k];
      refy_real[nhole][k] =refy[nhole]+ssy_off[nhole][k];      
      mark_real[nhole][k] =new TMarker(refy_real[nhole][k],refx_real[nhole][k],20);
      mark_real[nhole][k] -> SetMarkerColor(1);

      }
      // k is nfoil
      //      int k=5;
      //mark_ang[nhole] = new TMarker(ssy_cent_real[i]/l[k]/projectf[k],ssx_cent_real[i]/l[k]/projectf[k]);


      
      //      cout<<"nhole "<<nhole<<" i "<<i<<" j "<<j <<" refx "<<ssx_cent_real[j]<<" refy "<<ssy_cent_real[i]<<endl;
      //		if(pow(ssx-(refx[nhole]-0.8),2.0)/pow(selec_widthx,2.0)
      //	   + pow(ssy-(refy[nhole]-0.3),2.0)/pow(selec_widthy,2.0)<0.64){

      nhole++;

      
    }
  }
      for(int k=0;k<nfoil;k++){

	mark_real[38][k]->SetMarkerColor(2);
      	mark_real[23][k]->SetMarkerColor(2);
       

	
      }
      mark[38]->SetMarkerColor(2);
      mark[23]->SetMarkerColor(2);
};

//========== Mxpt ============================//
void angcalib::Mxpt(string matrix_xp){

  ifstream Mxpt(matrix_xp);
  

  for (int i=0;i<nParamT;i++){
    double par=0.0;
    int p=0;
    Mxpt >> par >> p >> p >> p >> p >> p;
    Pxpt[i]  = par;
    OptPar1[i] = par;
    //    cout<<Form("OptPar1[%d] ",i)<<OptPar1[i]<<endl;
  }
  Mxpt.close();  
};

//========== Mypt ============================//
void angcalib::Mypt(string matrix_yp){

  ifstream Mypt(matrix_yp);

  for (int i=0;i<nParamT;i++){
    double par=0.0;
    int p=0;
    Mypt >> par >> p >> p >> p >> p >> p;
    Pypt[i]  = par;
    OptPar2[i] = par;
  }
  Mypt.close();  
};

//========== Scaler Correction ==================//

void angcalib::Scale_corr(string offset_file){

  ifstream ifs(offset_file);
  for (int i=0 ; i<nfoil ; i++){
    ifs >> temp
	 >> offs_xp[i] >> scal_xp[i]
	 >> offs_yp[i] >> scal_yp[i];
    if(scal_yp[i]>0.5) offset_flag[i] = true;
    else offset_flag[i] = false;
  }

};



//============== Event Selection =========================//
void angcalib::EventSelect(bool rarm){

  bool ac_cut=false;
  bool PID_flag=false;
  bool Xpt_cut=false;
  bool Ypt_cut=false;
  tuning_num=0;
  ent=t2->GetEntries();
  cout<<"=============================================="<<endl;
  cout<<"======== Event Selection is starting! ========"<<endl;
  cout<<"=============================================="<<endl;
  cout<<"Event: "<<ent<<endl;
  cout<<"Tuning Events: "<<nmax<<endl;
  cout<<"rarm "<<rarm<<endl;
  //   ent=500000;
  d=div(ent,10000);

    // ----- Initialization ------- //
  for (int i=0 ; i< ent ; i++){

    if(i<nmax){
    x[i]    = -2222.0;
    y[i]    = -2222.0;
    xp[i]   = -2222.0;
    yp[i]   = -2222.0;
    foil_flag[i] = -1;}
    
    for(int j=0 ; j<max ; j++){
      r_x_fp[j]  = -2222.0;
      r_th_fp[j] = -2222.0;
      r_y_fp[j]  = -2222.0;
      r_ph_fp[j] = -2222.0;
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0;
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;
      Xpt[j]     = -2222.0;
      Ypt[j]     = -2222.0;
      Zt[j]      = -2222.0;
      Rth[j]     = -2222.0;
      Rph[j]     = -2222.0;
      Lth[j]     = -2222.0;
      Lph[j]     = -2222.0;             
    }




    ssx=-2222.0;
    ssy=-2222.0;
    Xt = -2222.0;
    Yt = -2222.0;
    holeg_temp = -1;
    foilg_temp = -1;
    holethrough = false;
    filled = false;
    trig1 = 0.0;
    trig4 = 0.0;
    trig5 = 0.0;
    rtrig = false;
    ltrig = false;
    gs_cut=false; //L-HRS gas cut
    ac_cut=false;
    PID_flag=false;
    Xpt_cut=false;
    Ypt_cut=false;
    
    t2->GetEntry(i);
    
    //if(i+evshift<ent) t2->GetEntry(i+evshift); 
    //else t2->GetEntry(i-ent+evshift);

    // ------------------------------------------- //
    // ------- Event selection for tuning -------- //
    // ------------------------------------------- //

    if(trig4>1.0) rtrig = true;
    else rtrig = false;
    if(trig1>1.0) ltrig = true;
    else ltrig = false;
    if(gs_asum>1000 || rarm==true) gs_cut=true;
    if(a1>100 && a2>2000)ac_cut=true;
    if((gs_cut && rarm == 0) || (rarm && ac_cut) )PID_flag=true;

    
    if(rarm==true){
    XFP[0]   = r_x_fp[0];
    XpFP[0]  = r_th_fp[0];
    YFP[0]   = r_y_fp[0];
    YpFP[0]  = r_ph_fp[0];
    Xpt[0]   = Rth[0];
    Ypt[0]   = Rph[0];    
    }
    else{
    XFP[0]   = l_x_fp[0];
    XpFP[0]  = l_th_fp[0];
    YFP[0]   = l_y_fp[0];
    YpFP[0]  = l_ph_fp[0];
    Xpt[0]   = Lth[0];
    Ypt[0]   = Lph[0];
    }


    
    XFP[0]  = (XFP[0]-XFPm)/XFPr;
    XpFP[0] = (XpFP[0]-XpFPm)/XpFPr;
    YFP[0]  = (YFP[0]-YFPm)/YFPr;
    YpFP[0] = (YpFP[0]-YpFPm)/YpFPr;

    Zt[0]   = (Zt[0]-Ztm)/Ztr;
    Xpt[0]  = (Xpt[0]-Xptm)/Xptr;
    Ypt[0] = (Ypt[0]-Yptm)/Yptr;





      Xpt[0]  = calcf2t_4th_2(Pxpt,
			   XFP[0], XpFP[0],
			   YFP[0], YpFP[0],
			   Zt[0]);

      Ypt[0] = calcf2t_4th_2(Pypt,
			   XFP[0], XpFP[0],
			   YFP[0], YpFP[0],
			   Zt[0]);


    
      


    XFP[0]  = XFP[0]  * XFPr + XFPm;
    XpFP[0] = XpFP[0] * XpFPr + XpFPm;
    YFP[0]  = YFP[0]  * YFPr + YFPm;
    YpFP[0] = YpFP[0] * YpFPr + YpFPm;    
    
    Ypt[0]  = Ypt[0]*Yptr +Yptm;
    Xpt[0]  = Xpt[0]*Xptr +Xptm;
    Zt[0]   = Zt[0]*Ztr +Ztm;
    

    h1->Fill(Ypt[0],Xpt[0]);

    
        if(fabs(Xpt[0]) < 0.08 
           && fabs(Ypt[0]) < 0.06 && PID_flag){
	  
	  for(int j=0 ; j<nfoil ; j++){

	    // int j=5;

 
	if(fcent[j]-selection_width<Zt[0]
	   && Zt[0]<fcent[j]+selection_width){

	  hz[j]->Fill(Zt[0]);
	  hz_all->Fill(Zt[0]);
	  //if(offset_flag[j]==true){ 
	  if(offset_flag[j]==true || offset_flag[j]==false){
	    // (in case you don't need scale+offset for event selection)
	    //ssx = (-Xpt*l[j]*projectf[j] + offs_xp[j])*scal_xp[j]; // for initial parameters (xpt_LHRS_4.dat)
	    //ssy = (-Ypt*l[j]*projectf[j] + offs_yp[j])*scal_yp[j]; // for initial parameters (ypt_LHRS_4.dat)

	    
	    ssx = -Xpt[0]*l[j]*projectf[j];
	    if(RHRSTrue==0)ssy = -Ypt[0]*l[j]*projectf[j];
	    if(RHRSTrue)ssy = -Ypt[0]*l[j]*projectf[j];
	    //if(j==8) ssy = ssy * 1.08;  // for second parameters
	    //if(j==9) ssy = ssy * 1.126; // for second parameters

	    hph_cut->Fill(Ypt[0]);
	    hssy_cut->Fill(ssy);
	    
	    h2_new[j]->Fill(ssy,ssx);
     	    h2_      ->Fill(ssy,ssx);
	    nhole=0;
	  
	    //      refx[nhole] = ssx_cent_real[j]-0.8;
	    //      refy[nhole] = ssy_cent_real[i]-0.3;
	    
	    for(int col=0 ; col<ncol ; col++){
	      for(int row=0 ; row<nrow ; row++){

		//		if(pow(ssx-(refx[nhole]-0.8),2.0)/pow(selec_widthx,2.0)
		//		   + pow(ssy-(refy[nhole]-0.3),2.0)/pow(selec_widthy,2.0)<0.64){

		//		if(pow(ssx-(refx[nhole]),2.0)/pow(selec_widthx,2.0)
		//		   + pow(ssy-(refy[nhole]),2.0)/pow(selec_widthy,2.0)<0.64){


		 // RHRS initial tuning setting //
		if(pow(ssx-(refx_real[nhole][j]),2.0)/pow(selec_widthx,2.0)
		   + pow(ssy-(refy_real[nhole][j]),2.0)/pow(selec_widthy,2.0)<0.64){

		  holeg_temp = nhole;
		  foilg_temp = j;
		  if(TFlag[holeg_temp][foilg_temp])holethrough = true;
		  
		  //	    holethrough = true;		  
		  //	  if(refx[nhole]==0.0)holethrough = true;

		}
		nhole++;
	      }
	    }
	  }


	  
	  if(ntune_event<nmax && holethrough==true
	     && filled==false ){
	    //	    foil_flag[ntune_event] = j;
	    //	    foil_flag[ntune_event] = 5;

	    foil_flag[ntune_event] = foilg_temp;
	    holegroup[ntune_event] = holeg_temp;
	    x[ntune_event]  = XFP[0];  // scaled 
	    y[ntune_event]  = YFP[0];  // scaled
	    xp[ntune_event] = XpFP[0]; // scaled
	    yp[ntune_event] = YpFP[0]; // scaled
	    z_recon[ntune_event] = Zt[0]; // not scaled
	    th[ntune_event] = Xpt[0];
	    ph[ntune_event] = Ypt[0];	    
	    ntune_event++;
	    filled=true;
	    //	    cout<<"ssx "<<ssx<<" ssy "<<ssy<<endl;
	    h3[j]->Fill(ssy,ssx);
	    h3_  ->Fill(ssy,ssx);
	    tuning_num=i;
	  }
	}
	  } //for (j)
	}
    //        if(i % 100000 == 0)cout<<i<<" / "<<ent<<endl;
   if(i % (d.quot*1000) == 0)cout<<i<<" / "<<ent<<endl;
   if(ntune_event >= nmax && BreakTrue)break;

  }

	    
  cout << " The number of events selected to be used for tuning: "
       << ntune_event << endl;
  cout << " The nnumver of total tuning events : "<<tuning_num<<endl;


};


//============ Tuning =========================//
void angcalib::Tuning(string ofMTPname){

  if (nite>0){
    cout << "======================================================" <<endl;
    cout << "=================  Tuning started =================== " << endl;
    cout << "======================================================" <<endl;}

  cout<<endl;
  cout<<" Matrix order : nn "<<nn<<" nnz "<<nnz<<" nParamT "<<nParamT<<endl;
  cout<<endl;
  
    const  char* new_tempc=ofMTPname.c_str();
    //    cout<<"new marix file: "<<new_tempc<<endl;

    
  for(int i=0 ; i<nite ; i++){

    cout<<"tuning i: "<<i+1<<" /"<<nite<<endl;

    // --------------------------- //
    // ---- Parameter tuning ----- //
    // --------------------------- //
    x[i] = i+1;

    ofstream * ofs1;
    ofstream * ofs2;

    if(Xp_flag){
      cout<<"------- Xp tuning -----"<<endl;
      chi_sq1[i] = tune(OptPar1,i,1);
      gchi_xp->SetPoint(i,i,chi_sq1[i]);    
    cout << " Xp Tuning# = " << i+1 << ": chisq = "
	 << chi_sq1[i] <<endl;
    sprintf(tempc,  "%s_xpt_%d.dat",new_tempc,i);
    cout<<"new matrix xp: "<<tempc<<endl;    
    ofs1 = new ofstream(tempc);
    }// xpt
    if(Yp_flag){cout<<"------- Yp tuning -----"<<endl;    
      chi_sq2[i] = tune(OptPar2,i,2);
        gchi_yp->SetPoint(i,i,chi_sq2[i]);
    cout << "YXp Tuning# = " << i+1 << ": chisq = "
	 << chi_sq2[i] <<endl;
    sprintf(tempc2, "%s_ypt_%d.dat",new_tempc,i);         
    cout<<"new matrix yp: "<<tempc2<<endl;
    ofs2 = new ofstream(tempc2);
    } // ypt

    
  
    int nppp = 0;

    for(int i=0 ; i<nn+1 ; i++){
      for(int e=0 ; e<nn+1 ; e++){
	for(int d=0 ; d<nn+1 ; d++){
	  for(int c=0 ; c<nn+1 ; c++){
	    for(int b=0 ; b<nn+1 ; b++){
	      for(int a=0 ; a<nn+1 ; a++){  
		if(a+b+c+d+e==i){
		  if(Xp_flag){
		  *ofs1 << OptPar1[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  /*
		  cout <<nppp
		       << " "<<OptPar1[nppp] 
		       << " " << a 
		       << " " << b
		       << " " << c
		       << " " << d
		       << " " << e << endl;*/
		  
		  }
		  if(Yp_flag){
		  *ofs2 << OptPar2[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  cout <<nppp
		       << " "<<OptPar2[nppp] 
		       << " " << a 
		       << " " << b
		       << " " << c
		       << " " << d
		       << " " << e << endl;
		  }
		  nppp++;

		}
	      }
	    }
	  }
	}
      }
    }
    if(Xp_flag){
    ofs1->close();
    ofs1->clear();
    }
    if(Yp_flag){
    ofs2->close();
    ofs2->clear();
    }


    
  }

  cout<<"========== Tuning is done ============="<<endl;
  
};


//============ Fill ======================//
void angcalib::Fill(bool rarm){

  cout<<"==========================================="<<endl;
  cout<<"========== Fill Tree & Hist =============== "<<endl;
  cout<<"=========================================== "<<endl;

  ent=t2->GetEntries();
  d=div(ent,10000);
  cout<<"Event: "<<ent<<endl;

  
  // ----- Initialization ------- //
  for (int i=0 ; i< ent ; i++){
    
    //    if(i<nmax)foil_flag[i] = -1; 
    for(int j=0 ; j<max ; j++){
      Xpt[j]=-2222.0;
      Ypt[j]=-2222.0;
      Xpt_tuned[j]=-2222.0;
      Ypt_tuned[j]=-2222.0;
      Rth[j]     = -2222.0;
      Rph[j]     = -2222.0;      
      r_x_fp[j]  = -2222.0;
      r_th_fp[j] = -2222.0;
      r_y_fp[j]  = -2222.0;
      r_ph_fp[j] = -2222.0;
      Lth[j]     = -2222.0;
      Lph[j]     = -2222.0;            
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0;
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;
      Zt[j]=-2222.0;
      
    }

    ssx=-2222.0;
    ssy=-2222.0;
    Xt = -2222.0;
    Yt = -2222.0;
    holeg_temp = -1;
    foilg_temp = -1;
    trig1 = 0.0;
    trig4 = 0.0;
    trig5 = 0.0;


    holethrough = false;
    filled = false;
    rtrig = false;
    ltrig = false;

    
      t2->GetEntry(i);

      // if(i+evshift<ent) t2->GetEntry(i+evshift); 
      //      else t2->GetEntry(i-ent+evshift);


      
      if(i <= tuning_num)tflag=1;
      else tflag=0;
      
    if(trig4>1.0) rtrig = true;
    else rtrig = false;
    if(trig1>1.0) ltrig = true;
    else ltrig = false;

    
    if(rarm==true){
    XFP[0]   = r_x_fp[0];
    XpFP[0]  = r_th_fp[0];
    YFP[0]   = r_y_fp[0];
    YpFP[0]  = r_ph_fp[0];
    Xpt[0]   =  Rth[0];
    Ypt[0]   =  Rph[0];    
    Xpt_init   = Rth[0];
    Ypt_init   = Rph[0];
    } else {
    XFP[0]   = l_x_fp[0];
    XpFP[0]  = l_th_fp[0];
    YFP[0]   = l_y_fp[0];
    YpFP[0]  = l_ph_fp[0];
    //    Xpt[0]   =  Lth[0];
    //    Ypt[0]   =  Lph[0];        
    Xpt_init   = Lth[0];
    Ypt_init   = Lph[0];
    }

    
    h1->Fill(Ypt[0],Xpt[0]);
    hth->Fill(Xpt[0]);

    
    XFP[0]  = (XFP[0]-XFPm)/XFPr;
    XpFP[0] = (XpFP[0]-XpFPm)/XpFPr;
    YFP[0]  = (YFP[0]-YFPm)/YFPr;
    YpFP[0] = (YpFP[0]-YpFPm)/YpFPr;
    Zt[0]   = (Zt[0]-Ztm)/Ztr;


    
    Xpt[0]  = (Xpt[0]-Xptm)/Xptr;
    Xpt_tuned[0] = (Xpt_tuned[0]-Xptm)/Xptr;
    Ypt[0]  = (Ypt[0]-Yptm)/Yptr;
    Ypt_tuned[0] = (Ypt_tuned[0]-Yptm)/Yptr;


    Xpt[0]  = calcf2t_4th_2(Pxpt,
			   XFP[0], XpFP[0],
			   YFP[0], YpFP[0],
			   Zt[0]);    

    Ypt[0] = calcf2t_4th_2(Pypt,
			   XFP[0], XpFP[0],
			   YFP[0], YpFP[0],
			   Zt[0]);
    

    Xpt_tuned[0] = calcf2t_4th_2(OptPar1,
			   XFP[0], XpFP[0],
			   YFP[0], YpFP[0],
			   Zt[0]);
    

    Ypt_tuned[0] = calcf2t_4th_2(OptPar2,
			   XFP[0], XpFP[0],
			   YFP[0], YpFP[0],
			   Zt[0]);


    
    Xpt[0]  = Xpt[0]*Xptr +Xptm; // scaled
    Ypt[0]  = Ypt[0]*Yptr +Yptm; // scaled    
    Xpt_tuned[0]  = Xpt_tuned[0]*Xptr +Xptm; // scaled   
    Ypt_tuned[0]  = Ypt_tuned[0]*Yptr +Yptm; // scaled
    Zt[0]   = Zt[0]*Ztr +Ztm; //scaled


    
    
    hth_c->Fill(Xpt[0]);
    

    for(int j=0 ; j<nfoil ; j++){
	if(fcent[j]-selection_width<Zt[0]
	   && Zt[0]<fcent[j]+selection_width){

	  zfoil=j;

	  /*
	  //===== w/o matrix tuing =====// 
    	    ssx = -Xpt_init*l[j]*projectf[j];
	    ssy = -Ypt_init*l[j]*projectf[j];
	    h3_a  ->Fill(ssy,ssx);
	    ssx=-2222.0;ssy=-2222.0;
	  //=== Input matrix tuing =====// 
    	    ssx = -Xpt[0]*l[j]*projectf[j];
	    ssy = -Ypt[0]*l[j]*projectf[j];
	    h3_b  ->Fill(ssy,ssx);
	    ssx=-2222.0;ssy=-2222.0;	    
	  */
	  
	    //===== w/  matrix tuing =====//
	  
    	    ssx = -Xpt_tuned[0]*l[j]*projectf[j];

	    if(RHRSTrue==0)ssy = -Ypt_tuned[0]*l[j]*projectf[j];
	    if(RHRSTrue)ssy = -Ypt_tuned[0]*l[j]*projectf[j];
	    h3_c  ->Fill(ssy,ssx);
	    ss_x=ssx;
	    ss_y=ssy;
	    if(RHRSTrue==0)h2[j]->Fill(ssy,ssx);
	    if(RHRSTrue)h2[j]->Fill(ssy,ssx);
	}
      }

    tnew->Fill();

   if(i % (d.quot*1000) == 0)cout<<i<<" / "<<ent<<endl;    
  }


  cout<<"FIll is Done !"<<endl;
};

//=========== Write Tree ==================//
void angcalib::Write(){

 cout<<"===================================="<<endl;
 cout<<"========= Write Tree & Hist ========"<<endl;
 cout<<"===================================="<<endl;


 

  tnew->Write();

  //----- Write Hist ----------//
  gchi_xp->SetName("gchi_xp");
  gchi_xp->Write();
  gchi_yp->SetName("gchi_yp");
  gchi_yp->Write();

  hph_cut->Write();
  hssy_cut->Write();
  h1  ->Write();
  h2_ ->Write();
  h3_->Write();
  h3_a->Write();        
  h3_b->Write();        
  h3_c->Write();    

  for(int i=0;i<nfoil;i++){
    h2[i]->Write();
    h2_new[i]->Write();
    h3[i]->Write();
    hz[i]->Write();
  }
  
  hz_all->Write();
  hth->Write();
  hth_c->Write();
  
};


//============ Close_tree  ================//

void angcalib::Close_tree(){

  f2->Close();
  fnew->Close();
  cout<<"===== closed root ======"<<endl;
}

//============ Draw ===================//
void angcalib::Draw(){
  c1=new TCanvas("c1","c1");
  c1->Divide(2,2);
  c1->cd(1);
  h2_new[4]->Draw("colz");
  c1->cd(2);
  h2_new[5]->Draw("colz");  
  c1->cd(3);
  h2_new[6]->Draw("colz");
  c1->cd(4);
  h2_new[7]->Draw("colz");

  TCanvas*  c10=new TCanvas("c10","c10");
  c10->Divide(2,2);
  c10->cd(1);
  h2_new[0]->Draw("colz");
  c10->cd(2);
  h2_new[1]->Draw("colz");  
  c10->cd(3);
  h2_new[2]->Draw("colz");
  c10->cd(4);
  h2_new[3]->Draw("colz");
  //  hz_all->Draw();
  //  hz[4]->Draw("same");
  //  hz[5]->Draw("same");
  //  hz[6]->Draw("same");  

  TCanvas* cx=new TCanvas("cx","cx");
  cx->cd();
  h2_new[5]->Draw();
  
  c0=new TCanvas("c0","c0");
  c0->Divide(2,2);
  c0->cd(1);
  h3[4]->Draw("colz");
  c0->cd(2);
  h3[5]->Draw("colz");  
  c0->cd(3);
  h3[6]->Draw("colz");
  c0->cd(4);
  h3_->Draw("colz");

  
  
  int nhole=0;

  
  for(int i=0; i<ncol ; i++){
    for(int j=0; j<nrow; j++){
      for(int k=1;k<5;k++){
	if(TFlag[nhole][k+3]){
	  c0->cd(k);
	  mark[nhole]->Draw("same");
	  mark_real[nhole][k+3]->Draw("same");
	  c1->cd(k);
	  mark[nhole]->Draw("same");
	  mark_real[nhole][k+3]->Draw("same");
	  cx->cd();
	  mark[nhole]->Draw("same");
	}
	if(TFlag[nhole][k-1]){
	  c10->cd(k);
	  mark[nhole]->Draw("same");
	  mark_real[nhole][k-1]->Draw("same");
	}

	
      }
      nhole++;
    }
  }

  
  TCanvas* c3=new TCanvas("c3","c3");
  c3->Divide(4,3);
  c3->cd(11);
  hz_all->Draw();  
  
  for(int i=0;i<nfoil;i++){
    c3->cd(i+1);
    h2_new[i]->Draw("colz");
    for(int j=0; j<nhole; j++){ mark[j]->Draw("same");}
    c3->cd(11);
    hz[i]->Draw("same");
  }


  
  cout<<"================================="<<endl;
  cout<<"========== Drawn ================"<<endl;
  cout<<"================================="<<endl;
};


//////////////////////////////////////////////////
double calcf2t_4th_2(double* P, double xf, double xpf, 
		     double yf, double ypf, double zt)
//////////////////////////////////////////////////
{
  // ------------------------------------------------ //
  // ----- 4rd order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //


  const int nMatT=nn;   
  const int nXf=nn;
  const int nXpf=nn;
  const int nYf=nn;
  const int nYpf=nn;
  const int nZt=nnz;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		  //cout<<"n "<<n<<" i "<<npar<<" x "<<x<<" a "<<a<<" b "<<b<<" c "<<c<<" d "<<d<<" e "<<e<<endl;
		}
		else{
		  x = 0.;
		}
		Y += x*P[npar];
		if(P[0]==OptPar2[0]){
		  //		  cout<<"npar "<<npar<<" Y "<<Y<<" x "<<x<<" Param "<<P[npar]<<endl;
		  //		  cout<<"xf "<<xf<<" xpf "<<xpf<<" yf "<<yf<<" ypf "<<ypf<<endl; 
		}
		npar++;
	      }
	      
	    }
	  }
	}
      }    
    }
  }


  return Y; 
  
}


// ######################################################
double calcf2t_zt(double* P, double xf, double xpf, 
                 double yf, double ypf){
// ######################################################
  const int nMatT=nnz; 
  const int nXf=nnz;
  const int nXpf=nnz;
  const int nYf=nnz;
  const int nYpf=nnz;

  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0;
  
  for (int n=0;n<nMatT+1;n++){
  	for (d=0;d<n+1;d++){
	  for (c=0;c<n+1;c++){
	    for (b=0;b<n+1;b++){
	      for (a=0;a<n+1;a++){
		
		if (a+b+c+d==n){
		  if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf){
		    x = pow(xf,double(a))*pow(xpf,double(b))*
		      pow(yf,double(c))*pow(ypf,double(d));
		  }
		  else{
		    x = 0.;
		  }
		  Y += x*P[npar];
		  npar++;
		}
		
	      }
	    }
	  }
  	}
  }
  
  return Y;
}



// #############################################################
double tune(double* pa, int j, int angflag) 
// #############################################################
{
  
  double chi2 = 0.0;
  double arglist[10]; 
  int ierflg = 0;
  int allparam = nParamT;
  TMinuit* minuit = new TMinuit(allparam);
  if(angflag==1){
    minuit->SetFCN(fcn1);
  }else{ minuit->SetFCN(fcn2);}

 
  double start[allparam];
  double step[allparam];
  const int nMatT=nn;  
  const int nXf=nn;
  const int nXpf=nn;
  const int nYf=nn;
  const int nYpf=nn;
  const int nZt=nnz; // The number of order is reduced for test (4-->2)
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  start[npar] = pa[npar];
		  step[npar] = 1.0e-3;
		}
		else{
		  start[npar] = 0.0;
		  step[npar]  = 0.0;
		}
		npar++;
	      }
	    }
	  }
	}
      }    
    }
  }
  
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1);
  
  double LLim[allparam];// Lower limit for all of the parameter
  double ULim[allparam];// Upper limit for all of the parameter
  char pname[500];
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);

    //    start[i] = pa[i]; 
    //LLim[i] = pa[i] - pa[i]*0.8;
    //ULim[i] = pa[i] + pa[i]*0.8;

    //    if(i<126)step[i]=0.0;
    //    else step[i]=1.0e-3;
    //    else step[i]=0.0;
    
    LLim[i] = pa[i] - 5.0; // temp
    ULim[i] = pa[i] + 5.0; // temp     
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
    
  }

  
  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0; // original
  //arglist[0] = 1.0; // test
  //arglist[0] = 0.0;   // test
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  
  // ~~~~ Migrad + Simplex  ~~~~ 
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg); // Chi-square minimization
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double er;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);

  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){
    //minuit -> GetParameter(i,par[i],e);
    if(angflag==1){minuit -> GetParameter(i,OptPar1[i],er);}
    else {minuit -> GetParameter(i,OptPar2[i],er);}
    
  }
  
  return chi2;
}

// #############################################################
void fcn1(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{

  
  const double sigma = 1.0;
  double ztR      = 0.0;
  double refpos   = 0.0;
  double residual = 0.0;
  double ang      = 0.0;
  double sspos    = 0.0;
  double total_chi2 = 0.0;
  
  double nev[nfoil][nsshole];
  double chi2[nfoil][nsshole];
  double w[nfoil][nsshole];

  
  
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      nev[i][j]  = 0.0;
      chi2[i][j] = 0.0;
      w[i][j]    = 1.0;
    }
  }

  
  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;
    ang    = 0.0;
    sspos  = 0.0;
    refpos = 0.0;  refpos = refx[holegroup[i]];
    ztR    = 0.0;  ztR    = z_recon[i];


    

    
    //    cout<<"ztR "<<ztR<<" foil "<<foil_flag[i]<<" hole "<<holegroup[i]<<endl;
    
    //if(foil_flag[i]==i) nev[i]++;

    x[i]  = (x[i]-XFPm)/XFPr;
    xp[i] = (xp[i]-XpFPm)/XpFPr;
    y[i]  = (y[i]-YFPm)/YFPr;
    yp[i] = (yp[i]-YpFPm)/YpFPr;
    ztR= (ztR-Ztm)/Ztr; // only zt was not scaled, so apply scaling here
    ang  = (ang-Xptm)/Xptr;

    
    ang = calcf2t_4th_2(param,
			x[i], xp[i],
			y[i], yp[i],
			ztR);

    
    x[i]  = x[i]  * XFPr + XFPm;    
    xp[i] = xp[i] * XpFPr + XpFPm;
    y[i]  = y[i]  * YFPr + YFPm;
    yp[i] = yp[i] * YpFPr + YpFPm;    
    ztR = ztR*Ztr +Ztm;
    ang = ang*Xptr +Xptm;

    
    sspos = -ang*l[foil_flag[i]]*projectf[foil_flag[i]]; // in centimeter

  

    // ------------------- //
    // --- Residual ------ //
    // ------------------- //
    
    residual = sspos - refpos;

    chi2[foil_flag[i]][holegroup[i]]
      = chi2[foil_flag[i]][holegroup[i]] + pow(residual,2.0);

    

    nev[foil_flag[i]][holegroup[i]]++;
  }

  
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      //      int i=5;   
      //if(nev[i][j]>0){
      if(nev[i][j]>10){ 
      //if(nev[i][j]>50){ // using only holes with more than 50 events
	chi2[i][j] = chi2[i][j]/nev[i][j]/pow(sigma,2.0);
	//      cout<<"i "<<i<<" j "<<j <<" chi2 "<<chi2[i][j]<<" nev "<<nev[i][j]<<endl;	
      }
      else chi2[i][j] = 0.0;
      total_chi2 = total_chi2 + chi2[i][j]*w[i][j];

						   
    }
    
  }
  
  fval = total_chi2;//(double)nfoil/(double)nsshole;
  
}


// #############################################################
void fcn2(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  
  const double sigma = 1.0;
  double ztR      = 0.0;
  double refpos   = 0.0;
  double residual = 0.0;
  double ang      = 0.0;
  double sspos    = 0.0;
  double total_chi2 = 0.0;
  
  double nev[nfoil][nsshole];
  double chi2[nfoil][nsshole];
  double w[nfoil][nsshole];
  
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      nev[i][j]  = 0.0;
      chi2[i][j] = 0.0;
      w[i][j]    = 1.0;
            
    }
  }
  
  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;
    ang    = 0.0;
    sspos  = 0.0;
    refpos = 0.0;  refpos = refy[holegroup[i]];
    ztR    = 0.0;  ztR    = z_recon[i];
    

    
    x[i]  = (x[i]-XFPm)/XFPr;
    xp[i] = (xp[i]-XpFPm)/XpFPr;
    y[i]  = (y[i]-YFPm)/YFPr;
    yp[i] = (yp[i]-YpFPm)/YpFPr;
    ztR= (ztR-Ztm)/Ztr; // only zt was not scaled, so apply scaling here
    ang  = (ang-Xptm)/Xptr;


    ang = calcf2t_4th_2(param,
			x[i], xp[i],
			y[i], yp[i],
			ztR);


    x[i]  = x[i]  * XFPr + XFPm;    
    xp[i] = xp[i] * XpFPr + XpFPm;
    y[i]  = y[i]  * YFPr + YFPm;
    yp[i] = yp[i] * YpFPr + YpFPm;    
    ztR = ztR*Ztr +Ztm;
    ang = ang*Yptr +Yptm;

    
    if(RHRSTrue==0)sspos = -ang*l[foil_flag[i]]*projectf[foil_flag[i]]; // in centimeter
    if(RHRSTrue)sspos = -ang*l[foil_flag[i]]*projectf[foil_flag[i]]; // in centimeter

    // ------------------- //
    // --- Residual ------ //
    // ------------------- //

    residual = sspos - refpos;
    chi2[foil_flag[i]][holegroup[i]]
      = chi2[foil_flag[i]][holegroup[i]] + pow(residual,2.0);
    
    nev[foil_flag[i]][holegroup[i]]++;
  }

  
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){


      if(nev[i][j]>0){ 
      //if(nev[i][j]>50){ // using only holes with more than 50 events
	chi2[i][j] = chi2[i][j]/nev[i][j]/pow(sigma,2.0);

      }
      else chi2[i][j] = 0.0;      
      total_chi2 = total_chi2 + chi2[i][j]*w[i][j];
    }
     }
  
  fval = total_chi2;


}




#endif

