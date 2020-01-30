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
  //  void HolePosi(bool rarm);
  //  void SSHole(string paraname, bool rarm);
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
  void SSHole(string paramname, bool rarm);
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
  double rasx,rasy;
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
  TH2F* hss[nfoil];
  TH2F* hang[nfoil];
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

  
  if(rarm)RHRSTrue=true;
  cout<<" HRS ARM ";
  if(RHRSTrue)cout<<" RHRS "<<endl;
  if(RHRSTrue==false)cout<<" LHRS "<<endl;
  
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
   t2->SetBranchAddress("Rrb.Raster2.rawcur.x",  &rasx);
   t2->SetBranchAddress("Rrb.Raster2.rawcur.y",  &rasy);
   //  t2->SetBranchAddress("R.tr.vz",Zt); //not tuned 
  } else{
   t2->SetBranchAddress("L.tr.vz_opt",ztR_opt);
   //    t2->SetBranchAddress("L.tr.vz",Zt);
   t2->SetBranchAddress("Lrb.Raster2.rawcur.x",  &rasx);
   t2->SetBranchAddress("Lrb.Raster2.rawcur.y",  &rasy);
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
    tnew->Branch("Rrb.Raster2.rawcur.x",&rasx,"rasx/D " );
    tnew->Branch("Rrb.Raster2.rawcur.y",&rasy,"rasy/D " );
  }else{
    tnew->Branch("L.tr.vz_opt",ztR_opt, "L.tr.vz_opt[100]/D" );
    tnew->Branch("L.tr.th_tg_opt",Xpt,"L.tr.tg_th_opt[100]/D " );
    tnew->Branch("L.tr.tg_ph_opt",Ypt,"L.tr.tg_ph_opt[100]/D" );
    tnew->Branch("L.tr.vz_tuned",Zt,"L.tr.vz_tuned[100]/D" );
    tnew->Branch("L.tr.tg_th_tuned",Xpt_tuned,"L.tr.tg_th_tuned[100]/D " );
    tnew->Branch("L.tr.tg_ph_tuned",Ypt_tuned,"L.tr.tg_ph_tuned[100]/D " );    
    tnew->Branch("Lrb.Raster2.rawcur.x",&rasx,"rasx/D " );
    tnew->Branch("Lrb.Raster2.rawcur.y",&rasy,"rasy/D " );
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

    hss[i] = new TH2F(Form("hss_%d",i), Form("hss_%d",i),
		     h1->GetXaxis()->GetNbins(),
		     h1->GetXaxis()->GetXmin(),
		     h1->GetXaxis()->GetXmax(),
		     h1->GetYaxis()->GetNbins(),
		     h1->GetYaxis()->GetXmin(),
		     h1->GetYaxis()->GetXmax());
    hss[i]->GetXaxis()->SetTitle("y_SS (cm)");
    hss[i]->GetYaxis()->SetTitle("x_SS (cm)");


    hang[i] = new TH2F(Form("hang_%d",i), Form("hang_%d",i),
		       200,-0.06,0.06,200,-0.06,0.06);

    hang[i]->GetXaxis()->SetTitle("#phi [rad]");
    hang[i]->GetYaxis()->SetTitle("#theta [rad]");

    

    
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


//========== SSHole =========================//

void angcalib::SSHole(string paramname, bool rarm){

  ifstream ifp(paramname.c_str());
  string buf;
  int hole, foil;
   double ssy_cent_real[nrow];
   double ssx_cent_real[ncol];


   for(int i=0;i<nfoil;i++){
     for(int j=0;j<nsshole;j++){
       w[i][j]=0.0;
       TFlag[j][i]=false;
     }
    l[i] = 0;
    l[i]=(l0-fcent_real[i]/cos(hrs_ang))*100.;    
    dth[i] = asin(l0*sin(hrs_ang)/(l0*cos(hrs_ang) -fcent_real[i]));
    //    cout<<"i "<<i<<" dth "<<dth[i]*180./3.14<<endl;

   }




   
   while(getline(ifp, buf)){
     if( buf[0]=='#' ){ continue; }
     if( ifp.eof() ) break;	
     hole=0; foil=0; 
     stringstream sbuf;
     sbuf <<buf;
     int flag=0;
     double a=0.,b=0.,c=0.;

     sbuf >> foil >> hole >> flag >> a >> b >> c;

     w[foil][hole] = a;
     ssx_off[hole][foil] = b;
     ssy_off[hole][foil] = c;

    if(flag==1)TFlag[hole][foil]=true;
    else if(flag==0)TFlag[hole][foil]=false;
    

    
  }//while
  

  
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

      refx_real[nhole][k] =refx[nhole] + ssx_off[nhole][k];
      refy_real[nhole][k] =refy[nhole] + ssy_off[nhole][k];      
      mark_real[nhole][k] =new TMarker(refy_real[nhole][k],refx_real[nhole][k],20);
      mark_real[nhole][k] -> SetMarkerColor(1);

      }

      nhole++;

      
    }
  }
      for(int k=0;k<nfoil;k++){

	mark_real[38][k]->SetMarkerColor(2);
	if(rarm)mark_real[23][k]->SetMarkerColor(2);
	else  mark_real[45][k]->SetMarkerColor(2);
       
      }
      mark[38]->SetMarkerColor(2);
      if(rarm) mark[23]->SetMarkerColor(2);
      else mark[45]->SetMarkerColor(2);
  


}




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
      if(fcent[j]-selection_width<Zt[0]
	 && Zt[0]<fcent[j]+selection_width){
	
	hz[j]->Fill(Zt[0]);
	hz_all->Fill(Zt[0]);
	//if(offset_flag[j]==true){ 
	if(offset_flag[j]==true || offset_flag[j]==false){
	 
	  if(RHRSTrue==0)ssy=l[j]*sin(atan(-Ypt[0]))/cos(dth[j]-atan(-Ypt[0]));
	  if(RHRSTrue)ssy=l[j]*sin(atan(-Ypt[0]))/cos(dth[j]+atan(-Ypt[0]));
	  
      double lx;
      if(ssy>0)lx=sqrt(pow(l[j],2.0) + pow(ssy,2.0) + 2.0*l[j]*ssy*sin(dth[j]));
      else     lx=sqrt(pow(l[j],2.0) + pow(ssy,2.0) - 2.0*l[j]*ssy*sin(dth[j]));
      ssx =-Xpt[0]*lx;
      
      
      hss[j]->Fill(ssy,ssx);
      hph_cut  ->Fill(Ypt[0]);
      hssy_cut ->Fill(ssy);
      
	  h2[j]->Fill(ssy,ssx);
	  h2_      ->Fill(ssy,ssx);
	  nhole=0;
	  
	  for(int col=0 ; col<ncol ; col++){
	    for(int row=0 ; row<nrow ; row++){
	      
	      
	      
	      // RHRS initial tuning setting //

	      if(pow(ssx-(refx_real[nhole][j]),2.0)/pow(selec_widthx,2.0)
		 + pow(ssy-(refy_real[nhole][j]),2.0)/pow(selec_widthy,2.0)<0.64){
		
		  holeg_temp = nhole;
		  foilg_temp = j;

		  if(TFlag[holeg_temp][foilg_temp])holethrough = true;
		  
		}
	      nhole++;
	    }
	  }
	}
	
	
	  
	  if(ntune_event<nmax && holethrough==true
	     && filled==false ){

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
	    h3[j]->Fill(ssy,ssx);
	    h3_  ->Fill(ssy,ssx);
	    tuning_num=i;
	  }
	}
	  } //for (j)
	}
	
   if(i % (d.quot*1000) == 0)cout<<i<<" / "<<ent<<endl;
   if(ntune_event >= nmax && BreakTrue)break;
   
  }
  
  
  cout << " The number of events selected to be used for tuning: "<< ntune_event << endl;
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

  
  for(int i=0 ; i<nite ; i++){

    cout<<"tuning i: "<<i+1<<" /"<<nite<<endl;

    // --------------------------- //
    // ---- Parameter tuning ----- //
    // --------------------------- //
    x[i] = i+1;

    ofstream * ofs1;
    ofstream * ofs2;
    int nppp = 0;
    
    if(Xp_flag){
      cout<<"------- Xp tuning -----"<<endl;
      chi_sq1[i] = tune(OptPar1,i,1);
      gchi_xp->SetPoint(i,i,chi_sq1[i]);    
      cout << " Xp Tuning# = " << i+1 << ": chisq = "
	 << chi_sq1[i] <<endl;
    }
      
    sprintf(tempc,  "%s_xpt_%d.dat",new_tempc,i);
    cout<<"new matrix xp: "<<tempc<<endl;    
    ofs1 = new ofstream(tempc);
    nppp = 0;

    for(int i=0 ; i<nn+1 ; i++){
      for(int e=0 ; e<nn+1 ; e++){
	for(int d=0 ; d<nn+1 ; d++){
	  for(int c=0 ; c<nn+1 ; c++){
	    for(int b=0 ; b<nn+1 ; b++){
	      for(int a=0 ; a<nn+1 ; a++){  
		if(a+b+c+d+e==i){

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
		  

		  nppp++;

		}
	      }
	    }
	  }
	}
      }
    }


    ofs1->close();
    ofs1->clear();
    


    
    if(Yp_flag){cout<<"------- Yp tuning -----"<<endl;    
      chi_sq2[i] = tune(OptPar2,i,2);
      gchi_yp->SetPoint(i,i,chi_sq2[i]);
      cout << "YXp Tuning# = " << i+1 << ": chisq = "
	   << chi_sq2[i] <<endl;
    }

    
    sprintf(tempc2, "%s_ypt_%d.dat",new_tempc,i);         
    cout<<"new matrix yp: "<<tempc2<<endl;
    ofs2 = new ofstream(tempc2);


    
    nppp = 0;
    for(int i=0 ; i<nn+1 ; i++){
      for(int e=0 ; e<nn+1 ; e++){
	for(int d=0 ; d<nn+1 ; d++){
	  for(int c=0 ; c<nn+1 ; c++){
	    for(int b=0 ; b<nn+1 ; b++){
	      for(int a=0 ; a<nn+1 ; a++){  
		if(a+b+c+d+e==i){

		  *ofs2 << OptPar2[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  /*
		  cout <<nppp
		       << " "<<OptPar2[nppp] 
		       << " " << a 
		       << " " << b
		       << " " << c
		       << " " << d
		       << " " << e << endl;
		  */
		  

		  nppp++;

		}
	      }
	    }
	  }
	}
      }
    }

    ofs2->close();
    ofs2->clear();
   // ypt

    
  

    
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
  cout<<"RHRS mode "<<RHRSTrue<<endl;
  
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
    Xpt[0]   =  Lth[0];
    Ypt[0]   =  Lph[0];        
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
	  
	    //===== w/  matrix tuing =====//
	  hang[j]->Fill(-Ypt_tuned[0],-Xpt_tuned[0]);
	  if(RHRSTrue==0) ssy=l[j]*sin(atan(-Ypt_tuned[0]))/cos(dth[j]-atan(-Ypt_tuned[0]));
	  if(RHRSTrue)    ssy=l[j]*sin(atan(-Ypt_tuned[0]))/cos(dth[j]+atan(-Ypt_tuned[0]));
	  
	  double lx;
	  if(ssy>0)lx=sqrt(pow(l[j],2.0) + pow(ssy,2.0) + 2.0*l[j]*ssy*sin(dth[j]));
	  else     lx=sqrt(pow(l[j],2.0) + pow(ssy,2.0) - 2.0*l[j]*ssy*sin(dth[j]));
	  ssx =-Xpt_tuned[0]*lx;


	    

	    //	    	  cout<<"l "<<l[j]<<" Ypt "<<Ypt[0]*180./3.14<<" ssy "<<ssy<<" ssx "<<ssx<<endl;	  
	    //    	    ssx = -Xpt_tuned[0]*l[j]*projectf[j];
	    //	    if(RHRSTrue==0)ssy = -Ypt_tuned[0]*l[j]*projectf[j];
	    //	    if(RHRSTrue)ssy = -Ypt_tuned[0]*l[j]*projectf[j];


	    h2_new[j]->Fill(ssy,ssx);
	    h3_c  ->Fill(ssy,ssx);
	    ss_x=ssx;
	    ss_y=ssy;
	    //	    if(RHRSTrue==0)h2[j]->Fill(ssy,ssx);
	    //	    if(RHRSTrue)h2[j]->Fill(ssy,ssx);
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
    hss[i]->Write();
    hang[i]->Write();
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
  h2[4]->Draw("colz");
  c1->cd(2);
  h2[5]->Draw("colz");  
  c1->cd(3);
  h2[6]->Draw("colz");
  c1->cd(4);
  h2[7]->Draw("colz");

  TCanvas*  c10=new TCanvas("c10","c10");
  c10->Divide(2,2);
  c10->cd(1);
  h2[0]->Draw("colz");
  c10->cd(2);
  h2[1]->Draw("colz");  
  c10->cd(3);
  h2[2]->Draw("colz");
  c10->cd(4);
  h2[3]->Draw("colz");
  //  hz_all->Draw();
  //  hz[4]->Draw("same");
  //  hz[5]->Draw("same");
  //  hz[6]->Draw("same");  

  TCanvas* cx=new TCanvas("cx","cx");
  cx->cd();
  h2_->Draw("colz");
  for(int j=0; j<nhole; j++){ mark[j]->Draw("same");}
  
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
  c3->cd(12);
  h2_->Draw("colz");
  for(int j=0; j<nhole; j++){ mark[j]->Draw("same");}
  
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
  const int nZt=nn;
  
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

    //    if(i<56)step[i]=0.0;
    
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
  double ypt=0.0;
  double nev[nfoil][nsshole];
  double chi2[nfoil][nsshole];
  double lx=0.0;
  double ssx,ssy;
  
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      nev[i][j]  = 0.0;
      chi2[i][j] = 0.0;
    }
  }


  
      for(int i=0 ; i<ntune_event ; i++){
	residual = 0.0;
	ang    = 0.0;
	sspos  = 0.0;
	refpos = 0.0;  refpos = refx[holegroup[i]];
	ztR    = 0.0;  ztR    = z_recon[i];
	ypt=0.0;  
	lx=0.0;
	ssx=0.0;
	ssy=0.0;


	  
    //if(foil_flag[i]==i) nev[i]++;

  
    x[i]  = (x[i]-XFPm)/XFPr;
    xp[i] = (xp[i]-XpFPm)/XpFPr;
    y[i]  = (y[i]-YFPm)/YFPr;
    yp[i] = (yp[i]-YpFPm)/YpFPr;
    ztR   = (ztR-Ztm)/Ztr; // only zt was not scaled, so apply scaling here
    ang  =  (ang-Xptm)/Xptr;

   
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

    ypt=ph[i];    

    if(RHRSTrue==0) ssy=l[foil_flag[i]]*sin(atan(-ypt))/cos(dth[foil_flag[i]] - atan(-ypt));
    if(RHRSTrue)    ssy=l[foil_flag[i]]*sin(atan(-ypt))/cos(dth[foil_flag[i]] + atan(-ypt));
    
    if(ssy>0)lx=sqrt(pow(l[foil_flag[i]],2.0) + pow(ssy,2.0) + 2.0*l[foil_flag[i]]*ssy*sin(dth[foil_flag[i]]));
    else     lx=sqrt(pow(l[foil_flag[i]],2.0) + pow(ssy,2.0) - 2.0*l[foil_flag[i]]*ssy*sin(dth[foil_flag[i]]));

    ssx =-ang*lx;

    // ------------------- //
    // --- Residual ------ //
    // ------------------- //
        
    residual = ssx - refpos;

    chi2[foil_flag[i]][holegroup[i]]
      = chi2[foil_flag[i]][holegroup[i]] + pow(residual,2.0);

    nev[foil_flag[i]][holegroup[i]]++;

      }

      
      for(int i=0 ; i<nfoil ; i++){
	for(int j=0 ; j<nsshole ; j++){
	  if(nev[i][j]>10){ 
	    chi2[i][j] = chi2[i][j]/nev[i][j]/pow(sigma,2.0);
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
  double ssy=0.0;
  double nev[nfoil][nsshole];
  double chi2[nfoil][nsshole];
  //  double w[nfoil][nsshole];


  
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      nev[i][j]  = 0.0;
      chi2[i][j] = 0.0;
    }
  }
  

  for(int i=0 ; i<ntune_event ; i++){

    residual = 0.0;
    ang    = 0.0;
    sspos  = 0.0;
    refpos = 0.0;  refpos = refy[holegroup[i]];
    ztR    = 0.0;  ztR    = z_recon[i];
    ssy=0.0;

    
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
    
    //    if(RHRSTrue==0)sspos = -ang*l[foil_flag[i]]*projectf[foil_flag[i]]; // in centimeter
    //    if(RHRSTxsrue)sspos = -ang*l[foil_flag[i]]*projectf[foil_flag[i]]; // in centimeter

    if(RHRSTrue==0) ssy=l[foil_flag[i]]*sin(atan(-ang))/cos(dth[foil_flag[i]]-atan(-ang));
    if(RHRSTrue)    ssy=l[foil_flag[i]]*sin(atan(-ang))/cos(dth[foil_flag[i]]+atan(-ang));
    
    
    // ------------------- //
    // --- Residual ------ //
    // ------------------- //


    
    //    residual = sspos - refpos;
    
    residual = ssy - refpos;
    chi2[foil_flag[i]][holegroup[i]]
      = chi2[foil_flag[i]][holegroup[i]] + pow(residual,2.0);
    
    nev[foil_flag[i]][holegroup[i]]++;


  }
  

  
  
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      if(nev[i][j]>10){
      //if(nev[i][j]>50){ // using only holes with more than 50 events
	chi2[i][j] = chi2[i][j]/nev[i][j]/pow(sigma,2.0);
      }else chi2[i][j] = 0.0;

      total_chi2 = total_chi2 + chi2[i][j]*w[i][j];
      //      cout<<"i "<<i<<" j "<<j<<"chi2 "<<chi2[i][j]<<" total_chi2 "<<total_chi2<<" w "<<w[i][j]<<endl;
    }
  }
  

  
  fval = total_chi2;


}




#endif

