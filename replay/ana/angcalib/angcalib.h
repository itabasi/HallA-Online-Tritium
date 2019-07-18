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
  void HolePosi();
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

 //-------- Mxpt && Mypt------------//
  double Pxpt[nParamT];  
  double Pypt[nParamT];


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

  
};

//=============================================//
//=========== anglecalib ======================//
//=============================================//
angcalib::angcalib(){
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);};

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
  tnew->Branch("R.tr.tg_th", &Rth,"R.tr.tg_th[100]/D");
  tnew->Branch("R.tr.tg_ph", &Rph,"R.tr.th_ph[100]/D");
  tnew->Branch("L.tr.tg_th", &Lth,"L.tr.tg_th[100]/D");
  tnew->Branch("L.tr.tg_ph", &Lph,"L.tr.tg_ph[100]/D");
  tnew->Branch("ss_x",&ss_x,"ss_x/D " );
  tnew->Branch("ss_y",&ss_y,"ss_y/D " );
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
  
  h1= new TH2F("h1","",200,-5.0,5.0,200,-8.0,8.0); 
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
    //    cout << i << " " << l[i] << " " << dth[i] << " " << projectf[i] << endl;
    sprintf(tempc,"h2_new_%d",i);
    h2_new[i] = (TH2F*)h2[i]->Clone(tempc);
    sprintf(tempc,"h3_%d",i);
    h3[i] = (TH2F*)h2[i]->Clone(tempc);
    h3[i]->SetMarkerColor(i+1);
  }
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

};



//========== HolePosi =======================//
void angcalib::HolePosi(){
  nhole = 0;
   double ssy_cent_real[nrow];
   double ssx_cent_real[ncol];
  for(int i=0; i<ncol ; i++){
    for(int j=0; j<nrow; j++){
      ssy_cent_real[i] = -3.0*step + step*i;
      if(j%2==0)ssy_cent_real[i] = ssy_cent_real[i] - step/2.0;
      ssx_cent_real[j] = 5.0*step - step*j;
      refx[nhole] = ssx_cent_real[j];
      refy[nhole] = ssy_cent_real[i];
      mark[nhole] = new TMarker(refy[nhole],refx[nhole],28);
         nhole++;

    }
  }

};

//========== Mxpt ============================//
void angcalib::Mxpt(string matrix_xp){

  ifstream Mxpt(matrix_xp);
  

  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt >> par >> p >> p >> p >> p >> p;
    Pxpt[i]  = par;
    OptPar1[i] = par;

  }
  Mxpt.close();  
};

//========== Mypt ============================//
void angcalib::Mypt(string matrix_yp){

  ifstream Mypt(matrix_yp);

  for (int i=0;i<nParamT;i++){
    double par=0.;
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
//ifs->std::close();
};



//============== Event Selection =========================//
void angcalib::EventSelect(bool rarm){

  
  ent=t2->GetEntries();
  cout<<"=============================================="<<endl;
  cout<<"======== Event Selection is starting! ========"<<endl;
  cout<<"=============================================="<<endl;
  cout<<"Event: "<<ent<<endl;
  cout<<"Tuning Events: "<<nmax<<endl;
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
    if(gs_asum>1000 || rarm==true)gs_cut=true;

    
    if(rarm==true){
    XFP[0]   = r_x_fp[0];
    XpFP[0]  = r_th_fp[0];
    YFP[0]   = r_y_fp[0];
    YpFP[0]  = r_ph_fp[0];

    }
    else{
    XFP[0]   = l_x_fp[0];
    XpFP[0]  = l_th_fp[0];
    YFP[0]   = l_y_fp[0];
    YpFP[0]  = l_ph_fp[0];

    }


    XFP[0]  = (XFP[0]-XFPm)/XFPr;
    XpFP[0] = (XpFP[0]-XpFPm)/XpFPr;
    YFP[0]  = (YFP[0]-YFPm)/YFPr;
    YpFP[0] = (YpFP[0]-YpFPm)/YpFPr;
    Zt[0]   = (Zt[0]-Ztm)/Ztr;

    
    Xpt[0]  = calcf2t_4th_2(Pxpt,
			   XFP[0], XpFP[0],
			   YFP[0], YpFP[0],
			   Zt[0]);
    Ypt[0] = calcf2t_4th_2(Pypt,
			   XFP[0], XpFP[0],
			   YFP[0], YpFP[0],
			   Zt[0]);

    Ypt[0]  = Ypt[0]*Yptr +Yptm;
    Xpt[0]  = Xpt[0]*Xptr +Xptm;
    Zt[0]   = Zt[0]*Ztr +Ztm;


    
    if(fabs(Xpt[0]) < 0.08 
       && fabs(Ypt[0]) < 0.06 && gs_cut){
      

      for(int j=0 ; j<nfoil ; j++){
	
	if(fcent[j]-selection_width<Zt[0]
	   && Zt[0]<fcent[j]+selection_width){
	  
	  h2[j]->Fill(-Ypt[0]*l[j]*projectf[j],-Xpt[0]*l[j]*projectf[j]);



	  
	  //if(offset_flag[j]==true){ 
	  if(offset_flag[j]==true || offset_flag[j]==false){
	    // (in case you don't need scale+offset for event selection)
	    //ssx = (-Xpt*l[j]*projectf[j] + offs_xp[j])*scal_xp[j]; // for initial parameters (xpt_LHRS_4.dat)
	    //ssy = (-Ypt*l[j]*projectf[j] + offs_yp[j])*scal_yp[j]; // for initial parameters (ypt_LHRS_4.dat)

	    
	    ssx = -Xpt[0]*l[j]*projectf[j];
	    ssy = -Ypt[0]*l[j]*projectf[j];

	    //if(j==8) ssy = ssy * 1.08;  // for second parameters
	    //if(j==9) ssy = ssy * 1.126; // for second parameters
	    
	    h2_new[j]->Fill(ssy,ssx);
     	    h2_      ->Fill(ssy,ssx);
	    nhole=0;
	  

	    
	    for(int col=0 ; col<ncol ; col++){
	      for(int row=0 ; row<nrow ; row++){
		if(pow(ssx-refx[nhole],2.0)/pow(selec_widthx,2.0)
		   + pow(ssy-refy[nhole],2.0)/pow(selec_widthy,2.0)<1.0){
		  holeg_temp = nhole;
		  foilg_temp = j;
		  holethrough = true;
		}
		nhole++;
	      }
	    }
	  }

	  if(ntune_event<nmax && holethrough==true
	     && filled==false ){
	    //foil_flag[ntune_event] = j;
	    foil_flag[ntune_event] = foilg_temp;
	    holegroup[ntune_event] = holeg_temp;
	    x[ntune_event]  = XFP[0];  // scaled 
	    y[ntune_event]  = YFP[0];  // scaled
	    xp[ntune_event] = XpFP[0]; // scaled
	    yp[ntune_event] = YpFP[0]; // scaled
	    z_recon[ntune_event] = Zt[0]; // not scaled
	    ntune_event++;

	    filled=true;
	    h3[j]->Fill(ssy,ssx);
	    h3_  ->Fill(ssy,ssx);
	  }
	}
      }
    }
    //        if(i % 100000 == 0)cout<<i<<" / "<<ent<<endl;
   if(i % (d.quot*1000) == 0)cout<<i<<" / "<<ent<<endl;
  }

	    
  cout << " The number of events selected to be used for tuning: "
       << ntune_event << endl;



};


//============ Tuning =========================//
void angcalib::Tuning(string ofMTPname){

  if (nite>0){
    cout << "======================================================" <<endl;
    cout << "=================  Tuning started =================== " << endl;
    cout << "======================================================" <<endl;}

    const  char* new_tempc=ofMTPname.c_str();
    cout<<"new marix file: "<<new_tempc<<endl;

    
  for(int i=0 ; i<nite ; i++){

    cout<<"tuning i: "<<i+1<<" /"<<nite<<endl;
    // --------------------------- //
    // ---- Parameter tuning ----- //
    // --------------------------- //
    x[i] = i+1;

    cout<<"------- Xp tuning -----"<<endl;
    chi_sq1[i] = tune(OptPar1,i,1);   // xpt
    cout<<"------- Yp tuning -----"<<endl;    
    chi_sq2[i] = tune(OptPar2,i,2);  // ypt

    
    cout << " Tuning# = " << i << ": chisq = "
	 << chi_sq1[i] << " "
	 << chi_sq2[i] << endl;
    cout << endl;
    gchi_xp->SetPoint(i,i,chi_sq1[i]);    
    gchi_yp->SetPoint(i,i,chi_sq2[i]);


    //    sprintf(tempc,  "../matrix/newpar_xpt_%d.dat",i); 
    //    sprintf(tempc2, "../matrix/newpar_ypt_%d.dat",i);
    sprintf(tempc,  "%s_xpt_%d.dat",new_tempc,i); 
    sprintf(tempc2, "%s_ypt_%d.dat",new_tempc,i);     

    cout<<"new matrix xp: "<<tempc<<endl;
    cout<<"new matrix yp: "<<tempc<<endl;
    ofstream * ofs1 = new ofstream(tempc);
    ofstream * ofs2 = new ofstream(tempc2);
    int nppp = 0;

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
		  *ofs2 << OptPar2[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  nppp++;
		  //	      cout << Plen_opt[nppp] 
		  //		   << " " << a 
		  //		   << " " << b
		  //		   << " " << c
		  //		   << " " << d << endl;
		}
	      }
	    }
	  }
	}
      }
    }
    ofs1->close();
    ofs1->clear();
    ofs2->close();
    ofs2->clear();



    
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

      

    if(trig4>1.0) rtrig = true;
    else rtrig = false;
    if(trig1>1.0) ltrig = true;
    else ltrig = false;

    
    if(rarm==true){
    XFP[0]   = r_x_fp[0];
    XpFP[0]  = r_th_fp[0];
    YFP[0]   = r_y_fp[0];
    YpFP[0]  = r_ph_fp[0];
    Xpt_init   = Rth[0];
    Ypt_init   = Rph[0];
    } else {
    XFP[0]   = l_x_fp[0];
    XpFP[0]  = l_th_fp[0];
    YFP[0]   = l_y_fp[0];
    YpFP[0]  = l_ph_fp[0];
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

    hth_c->Fill(Xpt[0]);

     Zt[0]   = Zt[0]*Ztr +Ztm; //scaled
    
      for(int j=0 ; j<nfoil ; j++){
	if(fcent[j]-selection_width<Zt[0]
	   && Zt[0]<fcent[j]+selection_width){


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
	  //===== w/  matrix tuing =====//	  
    	    ssx = -Xpt_tuned[0]*l[j]*projectf[j];
	    ssy = -Ypt_tuned[0]*l[j]*projectf[j];
	    h3_c  ->Fill(ssy,ssx);
	    ss_x=ssx;
	    ss_y=ssy;
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
  } 

  hth->Write();
  hth_c->Write();
  
};


//============ Close_tree  ================//

void angcalib::Close_tree(){

  f2->Close();
  fnew->Close();
}

//============ Draw ===================//
void angcalib::Draw(){
  c1=new TCanvas("c1","c1");
  c1->cd();
  h1->Draw();
  for(int i=0 ; i<nfoil ; i++){
    h2[i]->Draw("same");
  }

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
  //cout << allparam << endl;
  TMinuit* minuit = new TMinuit(allparam);
  if(angflag==1){
    minuit->SetFCN(fcn1);
  }
  else minuit->SetFCN(fcn2);


  double start[allparam];
  double step[allparam];
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=2; // The number of order is reduced for test (4-->2)
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
		  step[npar] = 0.0;
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
    //start[i] = pa[i]; 
    //step[i] = 1.0e-3;
    
    //LLim[i] = pa[i] - pa[i]*0.8;
    //ULim[i] = pa[i] + pa[i]*0.8;
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
    if(angflag==1){
      minuit -> GetParameter(i,OptPar1[i],er);
    }
    else minuit -> GetParameter(i,OptPar2[i],er);
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
    
    //if(foil_flag[i]==i) nev[i]++;

    ztR= (ztR-Ztm)/Ztr; // only zt was not scaled, so apply scaling here
    ang = calcf2t_4th_2(param,
			x[i], xp[i],
			y[i], yp[i],
			ztR);
    ang = ang*Xptr +Xptm;
    sspos = -ang*l[foil_flag[i]]*projectf[foil_flag[i]]; // in centimeter
    
    // ------------------- //
    // --- Residual ------ //
    // ------------------- //
    residual = sspos - refpos;
    //cout << i << ": " << sspos << "-(" << refpos << ")=" << residual << endl;
    
    chi2[foil_flag[i]][holegroup[i]]
      = chi2[foil_flag[i]][holegroup[i]] + pow(residual,2.0);
    
    nev[foil_flag[i]][holegroup[i]]++;
  }
  
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      
      //if(nev[i][j]>0){
      if(nev[i][j]>10){ 
      //if(nev[i][j]>50){ // using only holes with more than 50 events
	chi2[i][j] = chi2[i][j]/nev[i][j]/pow(sigma,2.0);
      }
      else chi2[i][j] = 0.0;
      
      total_chi2 = total_chi2 + chi2[i][j]*w[i][j];
    }
  }
  
  fval = total_chi2/(double)nfoil/(double)nsshole;
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
    
    //if(foil_flag[i]==i) nev[i]++;

    ztR= (ztR-Ztm)/Ztr; // only zt was not scaled, so apply scaling here
    ang = calcf2t_4th_2(param,
			x[i], xp[i],
			y[i], yp[i],
			ztR);
    ang = ang*Yptr +Yptm;
    sspos = -ang*l[foil_flag[i]]*projectf[foil_flag[i]]; // in centimeter
    
    // ------------------- //
    // --- Residual ------ //
    // ------------------- //
    residual = sspos - refpos;
    //cout << i << ": " << sspos << "-(" << refpos << ")=" << residual << endl;
    
    chi2[foil_flag[i]][holegroup[i]]
      = chi2[foil_flag[i]][holegroup[i]] + pow(residual,2.0);
    
    nev[foil_flag[i]][holegroup[i]]++;
  }
  
  for(int i=0 ; i<nfoil ; i++){
    for(int j=0 ; j<nsshole ; j++){
      
      //if(nev[i][j]>0){
      if(nev[i][j]>10){ 
      //if(nev[i][j]>50){ // using only holes with more than 50 events
	chi2[i][j] = chi2[i][j]/nev[i][j]/pow(sigma,2.0);
      }
      else chi2[i][j] = 0.0;
      
      total_chi2 = total_chi2 + chi2[i][j]*w[i][j];
    }
  }
  
  fval = total_chi2/(double)nfoil/(double)nsshole;
}




#endif

