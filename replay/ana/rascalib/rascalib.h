#ifndef rascalib_h
#define rascalib_h 1
#include <iostream>
#include <fstream>
//#include "Tree.h"
//#include "Tuning.h"
#include "tree.h"
#include "Param.h"
#include "Setting.h"
#include "define.h"
#include "zcalib.h"

bool single_flag=false;
double Ras_curx[nmax],Ras_cury[nmax];
double Ras_cur[nmax];
double parRaster[nParamTx];


//double foil_flag[nmax],x[nmax],y[nmax],xp[nmax],yp[nmax],z_recon[nmax],Zt_wRC[nmax];
 
//int nite=0;

extern  void func(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
extern  double Calc_ras(double a,double b,double c){return  a *b + c;};
extern  double Calc_z(double* P, double xf, double xpf,double yf, double ypf);

//=================================================//
//============ Raster calib Class =================//
//=================================================//

class rascalib :public tree{

 public:
  rascalib();
  ~rascalib();

  Setting* set;
  zcalib* zcorr;
  
  void Swich_Ras(bool ras_x);
  void MakeHist();
  void SetRoot(string ifname, bool sin);
  void NewRoot(string ofname);
  void MTRead(string ifname, bool rarm);
  void MTParam_x(string matrix_name_x, bool rarm);
  void MTParam_y(string matrix_name_y, bool rarm);
  void RasCorr(bool rarm);
  void zCorr(bool rarm);
  void EventSelection(bool rarm, bool ras_x);
  double Chi(double* pa,int j);
  void Tuning(string ofparam);
  void Fill(bool rarm, bool ras_x);
  void Close();


  
  /////////////////
  //// NewRoot ///
  ////////////////

  bool sieve=true;
  double ssx,ssy;

  
  TFile* ofp;
  TTree* tnew;
  double R_tr_vx_c[MAX],R_tr_vy_c[MAX];
  double L_tr_vx_c[MAX],L_tr_vy_c[MAX];
  double R_tr_vz_ras[MAX],L_tr_vz_ras[MAX];  

  //////////////////
  //// Make Hist ///
  //////////////////

  TH1F* h2[nfoil];
  TH1F* hz;
  TH1F* hz_cut;
  TH1F* hz_cut_r;
  TH1F* hz_cut_c;
  TH1F* hz_cut_rc;  
  TH1F* hz_ac;
  TH1F* hz_r;
  TH1F* hz_rc;
  TH1F* hz_c;
  TGraphErrors* gchi=new TGraphErrors();  
  ///////////////
  //// MTParam //
  //////////////
  string param[10];
  double parX[nParamTx];
  double parY[nParamTy];
  double parX_0,parX_1,parX_2,parX_3;
  double parY_0,parY_1,parY_2,parY_3;
  double Opt_Par_x[nParamTx],Opt_Par_y[nParamTy];
  double par[nParamT_ras];
  bool ras_x;
  /////////////////////////
  ///// Event Selection ///
  /////////////////////////
  double XFP,XpFP,YFP,YpFP;
  bool rtrig,ltrig;
  double Zt,Zt_r,Zt_tuned;
  int evshift=10000;
  double Ras_x,Ras_y;
  div_t d;
  double RasterCor;
  //////////////////////////
  ///// func //////////////
  /////////////////////////
  double Opt_Par[nParamTx];

  ///////////////////////
  //// Tuning //////////
  //////////////////////

  double chi_sq[100];
  
  
};


////////////////////////////////////////////////////

rascalib::rascalib(){
  
  set =new Setting();
  set->Initialize();
  zcorr=new zcalib();

  
}

////////////////////////////////////////////////////

rascalib::~rascalib(){}

////////////////////////////////////////////////////

void rascalib::MTRead(string ifname,bool rarm){

  cout<<endl;
  cout<<"==============================="<<endl;
  cout<<"=== Input Matrix Parameters ==="<<endl;
  cout<<"==============================="<<endl;

  string buf;
  int s=0;
  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >>param[s];
    cout<<param[s]<<endl;
    s++;
  }  



    zcorr->Mzt(param[0],rarm);
    zcorr->Mzt_L(param[0],rarm);    
    MTParam_x(param[1],rarm);
    //    MTParam_y(param[2],rarm);    
    ras_x=true;//raster-x correction 

    
}



///////////////////////////////////////////////////

void rascalib::MTParam_x(string matrix_x_name, bool rarm){


    char name_Mx[500];
    sprintf(name_Mx, matrix_x_name.c_str()); // optimized
    cout<<"HRS X parameters file: "<<name_Mx<<endl;
  ifstream Mx(name_Mx);

  for (int i=0;i<nParamTx;i++){
    //    Mx >> parX[i];
    Mx >> Opt_Par_x[i];
    cout<<"ParX: "<<Opt_Par_x[i]<<endl;    
  // ---- Raster x is tuned by this macro ---- //
    if(i==1 || i==3) {parX[i] = 0.0; Opt_Par_x[i] =0.0;}

  }


  for(int i=0;i<nParamTx;i++){
    par[i]=parX[i];
    Opt_Par[i]=Opt_Par_x[i];
  }
  

  Mx.close();    


}


////////////////////////////////////////////////////////

void rascalib::MTParam_y(string matrix_y_name, bool rarm){

    char name_My[500];
    sprintf(name_My, matrix_y_name.c_str()); // optimized
    cout<<"HRS Y  parameters file: "<<name_My<<endl;
   ifstream My(name_My);
   
  for (int i=0;i<nParamTy;i++){
    //    My >> parY[i];
    My >> Opt_Par_x[i];
  // ---- Raster y is tuned by this macro ---- //
    if(i==1 || i==3) {parY[i] = 0.0;  Opt_Par_x[i] = 0.0; }
  }


    for(int i=0;i<nParamTy;i++){
    par[i]=parY[i];
    Opt_Par[i]=Opt_Par_y[i];}


  My.close();    


}

///////////////////////////////////////////////////

void rascalib::SetRoot(string ifname,bool sin){

  //  SetTree(ifname);

  if(sin){ cout<<"===== Single mode ======="<<endl;SetTree(ifname);}
  else {  Chain_Tree(ifname); }
 
  SetBranch();
  t1->SetBranchAddress("R.tr.vz_tuned", ztR_opt); // raster current
  //    t1->SetBranchAddress("L.tr.vz_tuned", ztR_opt); // raster current
  //    if(sieve){
  //    t1->SetBranchAddress("ss_x",&ssx);
  //    t1->SetBranchAddress("ss_y",&ssy);
  //  }

}

///////////////////////////////////////////////////

void rascalib::NewRoot(string ofname){


  
  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");

  
  
  //  tnew=new TTree("T","Raster matrix tuning");
  tnew =t1->CloneTree(0);


  
   //------ RHRS ---------//
  tnew->Branch("R.tr.vx_c"              , R_tr_vx_c    ,"R.tr.vx_c[10]/D"   );
  tnew->Branch("R.tr.vy_c"              , R_tr_vy_c    ,"R.tr.vy_c[10]/D"   );
  tnew->Branch("R.tr.vz_ras"            , R_tr_vz_ras    ,"R.tr.vz_ras[10]/D"   );
  //  tnew->Branch("R.tr.vz_tuned"          , R_tr_vz_ras    ,"R.tr.vz_tuned[10]/D"   );
  //------ LHRS --------// 
  //  tnew->Branch("L.tr.vz"                , L_tr_vz_ras    ,"L.tr.vz[10]/D"   );
  tnew->Branch("L.tr.vx_c"              , L_tr_vx_c    ,"L.tr.vx_c[10]/D"   );
  tnew->Branch("L.tr.vy_c"              , L_tr_vy_c    ,"L.tr.vy_c[10]/D"   );
  //  tnew->Branch("L.tr.vz_tuned"            , L_tr_vz_ras    ,"L.tr.vz_tuned[10]/D"   );
  

}

/////////////////////////////////////////////////

void rascalib::MakeHist(){


  for(int i=0 ; i<nfoil ; i++){
    char hist_name[500];    
    sprintf(hist_name,"h2_%d",i);
    h2[i] = new TH1F(hist_name, hist_name,1000,-0.2,0.2);

    h2[i]->GetXaxis()->SetTitle("z-LHRS (m)");
    h2[i]->GetXaxis()->SetRangeUser(-0.2,0.2);

  }
  hz=new TH1F("hz","",1000,-0.2,0.2);
  set->SetTH1(hz,"w/o Raster correction z hist","[m]","Counts");
  hz->SetLineColor(1);
  hz_cut=new TH1F("hz_cut","",1000,-0.2,0.2);
  set->SetTH1(hz_cut,"w/o Raster correction z hist","[m]","Counts");
  hz_cut->SetLineColor(1);
  hz_cut_r=new TH1F("hz_cut_r","",1000,-0.2,0.2);
  set->SetTH1(hz_cut_r,"w/o Raster correction z hist","[m]","Counts");
  hz_cut_r->SetLineColor(2);


  hz_cut_c=new TH1F("hz_cut_c","",1000,-0.2,0.2);
  set->SetTH1(hz_cut_c,"w/o Raster correction z hist","[m]","Counts");
  hz_cut_c->SetLineColor(1);
  hz_cut_rc=new TH1F("hz_cut_rc","",1000,-0.2,0.2);
  set->SetTH1(hz_cut_rc,"w/o Raster correction z hist","[m]","Counts");
  hz_cut_rc->SetLineColor(2);    

  
  hz_r=new TH1F("hz_r","",1000,-0.2,0.2);
  set->SetTH1(hz_r,"w/ Raster correction z hist","[m]","Counts");
  hz_r->SetLineColor(3);
  hz_c=new TH1F("hz_c","",1000,-0.2,0.2);
  set->SetTH1(hz_c,"w/ z correction z hist","[m]","Counts");
  hz_c->SetLineColor(4);
  
  hz_rc=new TH1F("hz_rc","",1000,-0.2,0.2);
  set->SetTH1(hz_rc,"Afeter Tuning w/ Raster correction z hist","[m]","Counts");
  hz_rc->SetLineColor(2);

  hz_ac=new TH1F("hz_ac","",1000,-0.2,0.2);
  set->SetTH1(hz_ac,"Afeter Tuning w/ Raster correction z hist","[m]","Counts");
  hz_ac->SetLineColor(3);


  gchi->SetFillColor(2);
  gchi->SetMarkerStyle(20);
  gchi->SetMarkerColor(2);
  gchi->SetFillStyle(3005);
  gchi->SetMarkerSize(1.0);

  
}

/////////////////////////////////////////////////////////////

void rascalib::zCorr(bool rarm){


  
      //========== nomalization =================//
      XFP  = (XFP-XFPm)/XFPr; 
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
      
      if(rarm)Zt = Calc_z(Pzt, XFP, XpFP, YFP, YpFP); // calculated z position in RHRS
      else Zt = Calc_z(Pzt, XFP, XpFP, YFP, YpFP); // calculated z position in LHRS
      
      //============ scaling ====================//
      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;
      //      ztR[0] = ztR[0] * Ztr + Ztm;
      Zt = Zt * Ztr + Ztm; 
      //      cout<<"Zt "<<Zt<<endl;
}

/////////////////////////////////////////////////////////////

void rascalib::RasCorr(bool rarm){


      
      //==== Raster Tunng =====//
      if(ras_x) RasterCor = Calc_ras(Ras_x, Opt_Par[2], Opt_Par[0]);
      else RasterCor = Calc_ras(Ras_y, Opt_Par[2], Opt_Par[0]);
      //      cout<<"Opt_Par[2] "<<Opt_Par[2]<<" Opt_Par[0] "<<Opt_Par[0]<<endl;
      // ----------------------------------------------------------- //
      // ------ Converting raster current to raster position ------- //
      // ----------------------------------------------------------- //
      
      RasterCor = RasterCor/tan(hrs_ang);
      Zt_r = Zt + RasterCor;
      
      //      if(rarm)Rvz[0]=Zt_r;
      //      else {Lvz[0]=Zt_r;}


}


/////////////////////////////////////////////////////////////
void rascalib::EventSelection(bool rarm, bool ras_x){

  cout<<"================================="<<endl;
  cout<<"========= Event Selction ========"<<endl;
  cout<<"================================="<<endl;


  
    if(evshift>ENum)evshift=ENum;
   cout<<"Tuning Event: "<<nmax<<endl;
   if(ENum<10000)d=div(nmax,1000);
   else   d=div(nmax,1000);

   //   if(ENum>300000)ENum=300000;
   
  for(int i=0 ; i<nmax ; i++){
    x[i]    = -2222.0; // x at FP
    y[i]    = -2222.0; // y at FP
    xp[i]   = -2222.0; // x' at FP
    yp[i]   = -2222.0; // y' at FP
    foil_flag[i] = -1; // group number of foil
  }



  //  for (int i=0 ; i<nmax ; i++){
  for (int i=0 ; i<ENum ; i++){


    ///===========================//
    //==== Initialization ========//
    //============================//
    
     for(int j=0 ; j<MAX ; j++){
      r_x_fp[j]  = -2222.0;
      r_th_fp[j] = -2222.0;
      r_y_fp[j]  = -2222.0;
      r_ph_fp[j] = -2222.0;
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0;
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;      
    }

    
     Ras_x=-2222.0;
     Ras_y=-2222.0;
     if(sieve){
       ssx=-100.;
       ssy=-100;
     }
    
    rtrig = false;
    ltrig = false;
    
    t1->GetEntry(i);
    //    if(i+evshift<ENum) t1->GetEntry(i+evshift); 
    //    else t1->GetEntry(i-ENum+evshift);

    if(trig1>1.0) ltrig = true;
    else ltrig = false;
    if(trig4>1.0) rtrig = true;
    else rtrig =false;
    
    if(rarm){
    XFP   = r_x_fp[0];
    XpFP  = r_th_fp[0];
    YFP   = r_y_fp[0];
    YpFP  = r_ph_fp[0];
    //    Ras_x = L_Ras_x;
    Ras_x = R_Ras_x;
    Ras_y = R_Ras_y;
    Zt    = Rvz[0];
    }else{
    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    Ras_x = L_Ras_x;
    Ras_y = L_Ras_y;
    Zt    = Lvz[0];      
     }

    

    zCorr(rarm);
    RasCorr(rarm);
    hz_r->Fill(Zt); //w raster correction
    bool Sieve=false;
    //    if(fabs(Rth[0])<0.01 &&fabs(Rph[0])<0.01)Sieve=true;
    if(fabs(ssx)<1. &&fabs(ssy)<1.)Sieve=true;

    //    cout<<"ssx "<<ssx<<" ss "<<ssy<<" Z "<<Zt<<endl;


    if(((ltrig==true && rarm==false)|| (rarm==true && rtrig ==true))
       && fabs(XFP)  <2.0 
       && fabs(XpFP) <0.1
       && fabs(YFP)  <0.5
       && fabs(YpFP) <0.1
       && (L_gs_asum>1800 || (rarm && fabs(Rth[0])<0.01 &&fabs(Rph[0])<0.01 )) // pi^{-} rejection with gas Cherenkov 
       ){
      //    if(rarm && Sieve){
      //      cout<<"z "<<Zt<<" z_rc "<<Zt_r<<" Rasx "<<Ras_x<<endl;
      // ------------------------------------------- //
      // ------- Event selection for tuning -------- //
      // ------------------------------------------- //

      if(nite>0){
	for(int j=0 ; j<nfoil ; j++){
	  if(j<=4 && j<=6 || j==0 || j==9){
	  	    if(fcent[j]-selection_width<Zt+RasterCor 
	  	       && Zt+RasterCor<fcent[j]+selection_width){
	    //	    if((-0.03<Zt && Zt<-0.01 && j==4) || (-0.003<Zt && Zt<0.01 && j==5) || (0.02<Zt && Zt<0.04 && j==6) || (j==0 && -0.13<Zt && Zt<-0.11) || (j==9 && 0.15<Zt && Zt<0.13)){ // RHRS

	    //	    if()
	    
	    h2[j]->Fill(Zt+RasterCor);
	    if(j+2!=10) h2[j]->SetLineColor(j+2);
	    else h2[j]->SetLineColor(2);
	    h2[j]->SetLineStyle(9);
	    if(ntune_event<nmax){
	      foil_flag[ntune_event] = j;
	      x[ntune_event]  = XFP;
	      y[ntune_event]  = YFP;
	      xp[ntune_event] = XpFP;
	      yp[ntune_event] = YpFP;
	      Zrecon[ntune_event] = Zt; // Reconstructed z position
	      Zt_wRC[ntune_event] = Zt + RasterCor; // z position with raster correction
	      Ras_curx[ntune_event] = Ras_x;
	      Ras_cury[ntune_event] = Ras_y;
	      if(ras_x)Ras_cur[ntune_event] = Ras_x;
	      else Ras_cur[ntune_event] = Ras_y;
	      hz_cut->Fill(Zt); // w/o raster correction
	      hz_cut_r->Fill(Zt_r); // w/ raster correction	     
	      ntune_event++;
	      if(ntune_event %  1000 == 0)cout<<"tune_events "<<ntune_event<<" / "<<nmax<<endl;
	    }
	  }
	}
      } // nite>0
    }
    }
    if(i % (100000) == 0)cout<<i<<" / "<<nmax<<endl;
  }

  cout<<"================"<<endl;
  cout<<"   ntune_event "<<ntune_event<<endl;
  cout<<"================"<<endl;
}

//////////////////////////////////////////////////////////////////////////

 
double rascalib::Chi(double *pa, int j){

    //=================================//
    //==== Raster Calibration =========//
    //================================//


  int allparam=nParamTx;
  double chi2 = 0.0;
  double arglist[10]; 
  int ierflg = 0;
  allparam = nParamTx;
  cout << allparam << endl;
  TMinuit* minuit = new TMinuit(allparam);
  //  minuit->SetFCN(func);
  minuit->SetFCN(func);
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1);
  double start[allparam];
  double step[allparam];
  double LLim[allparam];// Lower limit for all of the parameter
  double ULim[allparam];// Upper limit for all of the parameter
  char pname[500];

  
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i]; // initial parameters
    //step[i] = 1.0e-5;      
    if(i==0){       // offset parameter for raster x
      step[i] = 1.0e-3;
    }
    else if(i==2){  // gradient parameters for raster x
      step[i] = 1.0e-8;
    }
    else step[i] = 0.0; // raster y is not tuned here
    
    //LLim[i] = pa[i] - pa[i]*0.8;
    //ULim[i] = pa[i] + pa[i]*0.8;
    
    LLim[i] = pa[i] -0.5; // temporary 
    ULim[i] = pa[i] +0.5; // temporary 

    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }

    // ~~~~ Strategy ~~~~
  //  arglist[0] = 2.0;
  arglist[0] = 1.0;
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  
  // ~~~~ Migrad + Simplex  ~~~~ 
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg); // Chi-square minimization
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double e;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);

  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){
    //minuit -> GetParameter(i,par[i],e);
    minuit -> GetParameter(i,Opt_Par[i],e);
  }
  
  return chi2;
  
}

//////////////////////////////////////////////////////////////////////////


void rascalib::Tuning(string ofparam){


  for(int i=0 ; i<nite ; i++){
    // --------------------------- //
    // ---- Parameter tuning ----- //
    // --------------------------- //
    x[i] = i+1;
    //    if(i==0){ 
    //  chi_sq[i] = Chi(par,i); 
    // }
    //  else { 
      chi_sq[i] = Chi(Opt_Par,i); 
      // }
    cout << " Tuning# = " << i << ": chisq = " << chi_sq[i] << endl; // 
    cout << endl;
    
    char tempc[500];
    
    sprintf(tempc, "%s_%d.dat",ofparam.c_str(),i); //new matrix will be stored here
    ofstream * ofs = new ofstream(tempc);
    
    *ofs << Opt_Par[0] << " " << Opt_Par[1] << " " << Opt_Par[2] << " "<< Opt_Par[3] << endl;
    cout << Opt_Par[0] << " " << Opt_Par[1] << " " << Opt_Par[2] << " "<< Opt_Par[3] << endl;

    gchi->SetPoint(i,i,chi_sq[i]);
    
    ofs->close();
    ofs->clear();
  }



}

////////////////////////////////////////////////////////////////////////

void rascalib::Fill(bool rarm, bool ras_x){

  cout<<"================================="<<endl;
  cout<<"=========== Fill ================"<<endl;
  cout<<"================================="<<endl;



  ENum=t1->GetEntries();
  if(evshift>ENum)evshift=ENum;
  cout<<"Fill Event: "<<ENum<<endl;
  if(ENum<10000)d=div(ENum,1000);
  else   d=div(ENum,1000);
  //  if(ENum>1000000)ENum=1000000;

  for (int i=0 ; i< ENum ; i++){

    
    ///===========================//
    //==== Initialization ========//
    //============================//

     for(int j=0 ; j<MAX ; j++){
      r_x_fp[j]  = -2222.0;
      r_th_fp[j] = -2222.0;
      r_y_fp[j]  = -2222.0;
      r_ph_fp[j] = -2222.0;
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0;
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;
      R_tr_vz_ras[j]=-2222.0;
      ztR_opt[j]=-100.;
    }

     Ras_x=-2222.0;
     Ras_y=-2222.0;
     Zt=-2222.0;
     Zt_r=-2222.0;
     R_a1_asum=-2222.0;
     R_a2_asum=-2222.0;
   //============================//

     t1->GetEntry(i);    
     
    if(rarm){
    XFP   = r_x_fp[0];
    XpFP  = r_th_fp[0];
    YFP   = r_y_fp[0];
    YpFP  = r_ph_fp[0];
    Ras_x = R_Ras_x;
    Ras_y = R_Ras_y;
    Zt    = Rvz[0];

    }else{
    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    Ras_x = L_Ras_x;
    Ras_y = L_Ras_y;
    Zt    = Lvz[0];      
     }



    hz->Fill(Zt); //w/o raster correction
    zCorr(rarm);
    hz_c->Fill(Zt); //after tuning hist
    RasCorr(rarm); //raster tuning
    hz_rc->Fill(Zt_r); //after tuning hist

    if(R_a1_asum<50 && R_a2_asum>2000)hz_ac->Fill(Zt_r);

    //    if(i % (d.quot * 1000) == 0)cout<<i<<" / "<<ENum<<endl;


    if(
       //       ((ltrig==true && rarm==false)|| (rarm==true && rtrig ==true)) &&
       fabs(XFP)  <2.0 
       && fabs(XpFP) <0.1
       && fabs(YFP)  <0.5
       && fabs(YpFP) <0.1
       && (L_gs_asum>1800 || (rarm && fabs(Rth[0])<0.01 &&fabs(Rph[0])<0.01 ))
       ){
	for(int j=0 ; j<nfoil ; j++){
	  if(fcent[j]-selection_width<Zt+RasterCor 
	     && Zt+RasterCor<fcent[j]+selection_width){

	      hz_cut_c->Fill(Zt); //w/o raster correction
	      hz_cut_rc->Fill(Zt_r); //w/o raster correction	     
	  } // nite>0
	}
    }


    if(rarm){
      R_tr_vz_ras[0]=Zt_r;
      //      Rvz[0]=Zt_r;
      ztR_opt[0]=Zt_r;
    }else {
      L_tr_vz_ras[0]=Zt_r;
      ztR_opt[0]=Zt_r;
      //      Lvz[0]=Zt_r;
    }    
    if(i%100000==0)cout<<i<<" / "<<ENum<<endl;
    tnew->Fill();   
  }

}
////////////////////////////////////////////////////////////////////////
void rascalib::Close(){


  tnew->Write();
  hz->Write();
  hz_cut->Write();
  hz_cut_r->Write();
  hz_cut_c->Write();
  hz_cut_rc->Write();
  hz_r->Write();
  hz_c->Write();
  hz_rc->Write();
  hz_ac->Write();
  gchi->SetName("gchi");
  gchi->Write();
  
  ofp->Close();

}


/////////////////////////////////////////////////////////////////////////
///////////////////////// Function //////////////////////////////////////
////////////////////////////////////////////////////////////////////////


double Calc_z(double* P, double xf, double xpf, double yf, double ypf)
{
  
  // -----3rd order ----- 
  // This is the third order claculation byb  using 35 parameter
  const int nMatT=nnz;  // These are for the RHRS use same for the LHRS
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


/////////////////////////////////////////////////////////////////////////

void func(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/){

  
  const double sigma = 0.00035;
  double ztR  = 0.0;
  double refz = 0.0;
  double nev[nfoil];
  double residual;
  double chi2[nfoil];
  double w[nfoil];
  double total_chi2 = 0.0;
  
  for(int i=0 ; i<nfoil ; i++){
    nev[i]  = 0.0;
    chi2[i] = 0.0;
    w[i]    = 1.0;
  }
  w[0]=5.0;
  w[9]=5.0;


  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;    
    refz = 0.0;  refz = fcent_real[foil_flag[i]];
    ztR  = 0.0;  ztR  = Zrecon[i];


    if(foil_flag[i]==i) nev[i]++;
    
    // ----------------------------- //
    // ----- Raster correction ----- //
    // ----------------------------- //
    double Ras_para= Ras_cur[i];//Ras_curx[i]; // x correction
    double Ras_Cor = Calc_ras(Ras_para, param[2], param[0]); 

    Ras_Cor = Ras_Cor/tan(hrs_ang);
    double ztR_wRC; 
    ztR_wRC = ztR + Ras_Cor; // z with raster correction

    //  cout<<"Par0 "<<param[0]<<" Par2 "<<param[2]<<" Ras_para "<<Ras_para
    //      <<" zt "<<ztR<<" Ras "<<Ras_Cor<<" zt_r "<<ztR_wRC<<" refz "<<refz<<endl;
    
    // ------------------- //
    // --- Residual ------ //
    // ------------------- //
    residual = ztR_wRC - refz; 
    chi2[foil_flag[i]] = chi2[foil_flag[i]] + pow(residual,2.0);
  }
  
  for(int i=0 ; i<nfoil ; i++){
    if(nev[i]!=0){
      chi2[i] = sqrt(chi2[i])/nev[i]/sigma;
    }
    total_chi2 = total_chi2 + chi2[i]*w[i];
  }
  
  fval = total_chi2/(double)nfoil;

  
}


//////////////////////////////////////////////////////////////////////////



#endif
