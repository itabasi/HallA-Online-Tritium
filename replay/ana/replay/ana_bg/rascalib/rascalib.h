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


  double Ras_curx[nmax],Ras_cury[nmax];
  double Ras_cur[nmax];
  double parRaster[nParamTx];

int nite=0;

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
  void Swich_Ras(bool ras_x);
  void MakeHist();
  void SetRoot(string ifname, bool rarm);
  void NewRoot(string ofname);
  void MTParam_x(string matrix_name_x, bool rarm);
  void MTParam_y(string matrix_name_y, bool rarm);
  void EventSelection(bool rarm, bool ras_x);

  //  void func(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
  double Chi(double* pa,int j);
  void Tuning(string ofparam);
  void Fill(bool rarm, bool ras_x);
  void Close();


  
  /////////////////
  //// NewRoot ///
  ////////////////

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
  TH1F* hz_r;
  TH1F* hz_rc;
  TGraphErrors* gchi=new TGraphErrors();  
  ///////////////
  //// MTParam //
  //////////////

  double parX[nParamTx];
  double parY[nParamTy];
  double parX_0,parX_1,parX_2,parX_3;
  double parY_0,parY_1,parY_2,parY_3;
  double Opt_Par_x[nParamTx],Opt_Par_y[nParamTy];
  double par[nParamT_ras];
  
  /////////////////////////
  ///// Event Selection ///
  /////////////////////////
  double XFP,XpFP,YFP,YpFP;
  bool rtrig,ltrig;
  double Zt,Zt_r,Zt_tuned;
  int evshift=10000;
  double Ras_x,Ras_y;
  div_t d;
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
  
}

////////////////////////////////////////////////////

rascalib::~rascalib(){}

///////////////////////////////////////////////////


void rascalib::SetRoot(string ifname, bool rarm){

  SetTree(ifname);
  SetBranch();
  
  
}

///////////////////////////////////////////////////

void rascalib::NewRoot(string ofname){

  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  tnew=new TTree("T","Raster matrix tuning");
  tnew =t1->CloneTree(0);


  
   //------ RHRS ---------//
  tnew->Branch("R.tr.vx_c"              , R_tr_vx_c    ,"R.tr.vx_c[10]/D"   );
  tnew->Branch("R.tr.vy_c"              , R_tr_vy_c    ,"R.tr.vy_c[10]/D"   );
  tnew->Branch("R.tr.vx_ras"            , R_tr_vz_ras    ,"R.tr.vz_ras[10]/D"   );
  //------ LHRS --------// 
  tnew->Branch("L.tr.vx_c"              , L_tr_vx_c    ,"L.tr.vx_c[10]/D"   );
  tnew->Branch("L.tr.vy_c"              , L_tr_vy_c    ,"L.tr.vy_c[10]/D"   );
  tnew->Branch("L.tr.vx_ras"            , L_tr_vz_ras    ,"L.tr.vz_ras[10]/D"   );
  

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
  hz_r=new TH1F("hz_r","",1000,-0.2,0.2);
  set->SetTH1(hz_r,"w/ Raster correction z hist","[m]","Counts");

  hz_rc=new TH1F("hz_rc","",1000,-0.2,0.2);
  set->SetTH1(hz_rc,"Afeter Tuning w/ Raster correction z hist","[m]","Counts");
  hz_rc->SetLineColor(2);


  gchi->SetFillColor(2);
  gchi->SetMarkerStyle(20);
  gchi->SetMarkerColor(2);
  gchi->SetFillStyle(3005);
  gchi->SetMarkerSize(1.0);

  
}
/////////////////////////////////////////////////
void rascalib::MTParam_x(string matrix_x_name, bool rarm){


    char name_Mx[500];
    sprintf(name_Mx, matrix_x_name.c_str()); // optimized
    cout<<"HRS X parameters file: "<<name_Mx<<endl;
  ifstream Mx(name_Mx);

  for (int i=0;i<nParamTx;i++){
    Mx >> parX[i];
    Mx >> Opt_Par_x[i];
  // ---- Raster x is tuned by this macro ---- //
    if(i==1 || i==3) {parX[i] = 0.0; Opt_Par_x[i] =0.0;}
    
  }

  
  for(int i=0;i<nParamTx;i++){
    par[i]=parX[i];
    Opt_Par[i]=Opt_Par_x[i];}
  
  
  Mx.close();    


}


////////////////////////////////////////////////////////

void rascalib::MTParam_y(string matrix_y_name, bool rarm){

    char name_My[500];
    sprintf(name_My, matrix_y_name.c_str()); // optimized
    cout<<"HRS Y  parameters file: "<<name_My<<endl;
   ifstream My(name_My);
   
  for (int i=0;i<nParamTy;i++){
    My >> parY[i];
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

void rascalib::EventSelection(bool rarm, bool ras_x){

  cout<<"================================="<<endl;
  cout<<"========= Event Selction ========"<<endl;
  cout<<"================================="<<endl;
  
    if(evshift>ENum)evshift=ENum;
   cout<<"Get Entries: "<<ENum<<endl;
   if(ENum<10000)d=div(ENum,1000);
   else   d=div(ENum,10000);

  for(int i=0 ; i<nmax ; i++){
    x[i]    = -2222.0; // x at FP
    y[i]    = -2222.0; // y at FP
    xp[i]   = -2222.0; // x' at FP
    yp[i]   = -2222.0; // y' at FP
    foil_flag[i] = -1; // group number of foil
  }

  
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

    
    rtrig = false;
    ltrig = false;
    
    //t1->GetEntry(i);
    if(i+evshift<ENum) t1->GetEntry(i+evshift); 
    else t1->GetEntry(i-ENum+evshift);

    if(trig1>1.0) ltrig = true;
    else ltrig = false;
    if(trig4>1.0) rtrig = true;
    else rtrig =false;
    
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

    hz->Fill(Zt);
    
    
    if(((ltrig==true && rarm==false)|| (rarm==true && rtrig ==true))
       && fabs(XFP)  <2.0 
       && fabs(XpFP) <0.1
       && fabs(YFP)  <0.5
       && fabs(YpFP) <0.1
       && L_gs_asum>1800   // pi^{-} rejection with gas Cherenkov 
       ){

      /* 
      XFP  = (XFP-XFPm)/XFPr; 
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
      ztR[0] = Calc_z(Pzt_L, XFP, XpFP, YFP, YpFP); // calculated z position
      
      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;
      ztR[0] = ztR[0] * Ztr + Ztm; 
      */

      
      // ----------------------------------------------------------- //
      // ------ Converting raster current to raster position ------- //
      // ----------------------------------------------------------- //

      
      double RasterCor;
      if(ras_x) RasterCor= Calc_ras(Ras_x, par[2], par[0]);
      else RasterCor= Calc_ras(Ras_y, par[2], par[0]);
    
      RasterCor = RasterCor/tan(hrs_ang);
      Zt_r = Zt+RasterCor;

      hz_r->Fill(Zt_r);
      
      // ------------------------------------------- //
      // ------- Event selection for tuning -------- //
      // ------------------------------------------- //
      if(nite>0){
	for(int j=0 ; j<nfoil ; j++){
	  if(fcent[j]-selection_width<Zt+RasterCor 
	     && Zt+RasterCor<fcent[j]+selection_width){
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
	      z_recon[ntune_event] = Zt; // Reconstructed z position
	      //	      Zt_wRC[ntune_event] = Zt + RasterCor; // z position with raster correction
	      Ras_curx[ntune_event] = Ras_x;
	      Ras_cury[ntune_event] = Ras_y;
	      if(ras_x)Ras_cur[ntune_event] = Ras_x;
	      else Ras_cur[ntune_event] = Ras_y;
	     
	      ntune_event++;
	    }
	  }
	}
      } // nite>0
    }
        if(i % (d.quot * 1000) == 0)cout<<i<<" / "<<ENum<<endl;
  }
  
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
  //cout << allparam << endl;
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
      step[i] = 1.0e-5;
    }
    else if(i==2){  // gradient parameters for raster x
      step[i] = 1.0e-10;
    }
    else step[i] = 0.0; // raster y is not tuned here
    
    //LLim[i] = pa[i] - pa[i]*0.8;
    //ULim[i] = pa[i] + pa[i]*0.8;
    LLim[i] = pa[i] -1.0; // temporary 
    ULim[i] = pa[i] +1.0; // temporary 
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }

    // ~~~~ Strategy ~~~~
  arglist[0] = 2.0;
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
    }

     Ras_x=-2222.0;
     Ras_y=-2222.0;
     Zt=-2222.0;

   //============================//
     

     
    t1->GetEntry(i);    

    if(rarm){
    Ras_x = R_Ras_x;
    Ras_y = R_Ras_y;
    Zt    = Rvz[0];
    }else{
    Ras_x = L_Ras_x;
    Ras_y = L_Ras_y;
    Zt    = Lvz[0];      
     }

    
      // ----------------------------------------------------------- //
      // ------ Converting raster current to raster position ------- //
      // ----------------------------------------------------------- //

 
    double RasterCor;
    if(ras_x) RasterCor = Calc_ras(Ras_x, Opt_Par[2], Opt_Par[0]);
    else RasterCor = Calc_ras(Ras_y, Opt_Par[2], Opt_Par[0]);

    
      RasterCor = RasterCor/tan(hrs_ang);
      Zt_r = Zt+RasterCor;
      if(rarm)Rvz[0]=Zt_r;
      else Lvz[0]=Zt_r;
      hz_rc->Fill(Zt_r); //after tuning hist

    if(i % (d.quot * 1000) == 0)cout<<i<<" / "<<ENum<<endl;  
    tnew->Fill();
      
    }

}
////////////////////////////////////////////////////////////////////////
void rascalib::Close(){

  tnew->Write();
  hz->Write();
  hz_r->Write();
  hz_rc->Write();
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
  
  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;
    refz = 0.0;  refz = fcent_real[foil_flag[i]];
    ztR  = 0.0;  ztR  = z_recon[i];

    if(foil_flag[i]==i) nev[i]++;
    
    // ----------------------------- //
    // ----- Raster correction ----- //
    // ----------------------------- //
    double Ras_para= Ras_cur[i];//Ras_curx[i]; // x correction
    double Ras_Cor = Calc_ras(Ras_para, param[2], param[0]); 
    Ras_Cor = Ras_Cor/tan(hrs_ang);
    double ztR_wRC; 
    ztR_wRC = ztR + Ras_Cor; // z with raster correction
    
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
