#ifndef momcalib_h
#define momcalib_h 1

#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"
#include "define.h"
#include "tree.h"


bool RHRS;
bool LHRS;
const int MAX=1000;



extern double s2f1_off(int i,char* ARM,char* MODE, int KINE);
extern void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
//extern double tune(double* pa, int j, int angflag); 
extern  double Calc_ras(double a,double b,double c){return  a *b + c;};  
extern double calcf2t_ang(double* P,double xf, double xpf, double yf, double fpf,double z);
extern double calcf2t_zt(double* P, double xf, double xpf, double yf, double ypf);
extern double calcf2t_mom(double* P, double xf, double xpf, double yf, double ypf, double zt);

extern double calcf2t_mom_RL( double* RP, double Rxf, double Rxpf, double Ryf, double Rypf, double Rzt,
			      double* LP, double Lxf, double Lxpf, double Lyf, double Lypf, double Lzt	     );



//=================================================//
//============= MomCalib Class ====================//
//=================================================//

class momcalib : public tree
{

 public:
  momcalib();
  ~momcalib();

 public:
  //  void SetBranch(string ifname, bool rarm);
  //  void NewBranch(string ofname, bool rarm);
  void EventSelection();
  void SetRoot(string ifname);
  void NewRoot(string ofname);
  void MTParam(string mtparam);
  void MTParam_R();
  void MTParam_L();
  void MTP_mom();
  int mode(string arm);
  void Mzt(string ifname);
  void Mxpt(string ifname);
  void Mypt(string ifname);
  void MomTuning(string ofname);
  double Eloss(double xp,double z,  char* arm);
  //  void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
  double tune(double* pa, int j, int angflag);   
  // void Mpt(string ifname);
  void Fill();
  void Close();

  int MODE=5;
  
  //------ MTParam --------//
  string param[10];
  

  //------ NewROOT -----------//
  TFile* ofp;
  TTree* tnew;

  //------ MTtuing -----------//
  //  double Pzt[1000],Pzt_L[1000];
  //  double Pxpt[1000],Pypt[1000],Pxpt_L[1000],Pypt_L[1000];
  //  double Prp[1000],Plp[1000];

  //------ EventSelection ---//
  div_t d;
    double mmL_range,mmS_range;
    int nlam=0;
    int nsig=0;
    
 // Cut Parameters //
 const double coin_cutmin=-248;
 const double coin_cutmax=-244; 
 const double rpathl_cutmin=28.7;
 const double rpathl_cutmax=29.4;
 const double lpathl_cutmin=28.6;
 const double lpathl_cutmax=29.2;
 const double rbeta_cutmin=0.0;
 const double rbeta_cutmax=1.0;
 const double lbeta_cutmin=0.9;
 const double lbeta_cutmax=1.0;
 const double Rvz_cutmin=-0.1;
 const double Rvz_cutmax= 0.1;
 const double Lvz_cutmin=-0.1;
 const double Lvz_cutmax= 0.1;
 const double Rx_cutmin= -0.4;
 const double Rx_cutmax= 0.4;
 const double a1_th=100.0;
 const double a2_th=500.0; 
 // Event Flag //
    bool Rs0_flag;
    bool Ls0_flag;
    bool Rs2_flag;
    bool Ls2_flag;
    bool Scinti_flag;
    bool ac1_flag;
    bool ac2_flag;
    bool gs_flag;
    bool sh_flag;
    bool ps_flag;
    bool trig_flag;
    bool track_flag;
    bool Rpathl_flag;
    bool Lpathl_flag;
    bool coin_flag;
    bool x_flag;
    bool y_flag;
    bool z_flag;
    bool RPID_flag;
    bool LPID_flag;

    
 double mm;
 double mm_L;
 double mm_acc;
 double mm_pi;
 double ct;
 double ctime;
 double Ee,Ee_,Ek,Epi;
 double pe,pe_,pk,ppi;
 double coin_t,coin_tc;
 double rbeta,rbeta_k,lbeta;
 double Rs2_off,Ls2_off; 
 double Rs2_tcorr,Ls2_tcorr;
 int Ls2pads,Rs2pads;
 bool cut_ac1,cut_ac2,cut_beta;
 int nac1,nac2,nac3,n;
 double tof_r,tof_l; 
 double rpathl,lpathl;
 double corr_R,corr_L;
 double rpath_corr,lpath_corr;
 //---- MomTuning ----//
 double chi_sq[1000],chi_sq1[100],chi_sq2[100];
 TGraphErrors* gchi_p;
 TGraphErrors* gchi_Rp;
 TGraphErrors* gchi_Lp;
 //----- ELoss ------//
 

};

//=====================================//
//======== Momcalib ==================//
//=====================================//

momcalib::momcalib(){};
momcalib::~momcalib(){};


int momcalib::mode(string arm){
  if(arm=="R")MODE=-1;
  else if(arm=="L")MODE=1;
  else if(arm=="C")MODE=0;
  else{cout<<"Please Select Tuning Mode !"<<endl; exit(1);}
  return MODE;
}





//-----------------------------//
//--------- SetRoot -----------//
//-----------------------------//

void momcalib::SetRoot(string ifname){

  ChainTree(ifname);
  SetBranch();

};

//-----------------------------//
//--------- SetRoot -----------//
//-----------------------------//

void momcalib::NewRoot(string ofname){

  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  tnew=new TTree("T","Momentum matrix tuning");
  tnew =T->CloneTree(0);


  tnew->Branch("mm", &mm,"mm/D");  
  tnew->Branch("mm_acc", &mm_acc,"mm_acc/D");  
  tnew->Branch("mm_L", &mm_L,"mm_L/D");  
  tnew->Branch("mm_pi", &mm_pi,"mm_pi/D");
  tnew->Branch("coin", &coin_t,"coin_t/D");  
  tnew->Branch("coin_acc", &ctime,"coin_acc/D");  

};

//=====================================================================//

void momcalib::MTP_mom(){

  //====== RHRS Momentum parameters ========//
    char name_Mpt[500];
    sprintf(name_Mpt, param[8].c_str()); // optimized
    ifstream Mpt(name_Mpt);

   for(int i=0;i<nParamTp;i++){
    double par=0.;
    int p=0;
    Mpt >> par >> p >> p >> p >> p;
    Prp[i]=par;
    Opt_par_R[i]=par;
    Opt_par[i]=par;
   }
  Mpt.close();

  //====== LHRS Momentum parameters ========//
    char name_Mpt_L[500];
    sprintf(name_Mpt_L, param[9].c_str()); // optimized
    ifstream Mpt_L(name_Mpt_L);

   for(int i=0;i<nParamTp;i++){
    double par=0.;
    int p=0;
    Mpt_L >> par >> p >> p >> p >> p;
    Plp[i]=par;
    Opt_par_L[i]=par;
    Opt_par[i+nParamTp]=par;  // Both momentum paramters  
   }
  Mpt_L.close();





}
//=====================================================================//

void momcalib::MTParam_R(){

  //=================//
  //==== RHRS =======//
  //=================//

  
  //====== RHRS z parameters ======//

    char name_Mzt[500];
    sprintf(name_Mzt, param[0].c_str()); // optimized
    ifstream Mzt(name_Mzt);

   for(int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt >> par >> p >> p >> p >> p;
    Pzt[i]=par;

   }
  Mzt.close();

  
  //====== RHRS raster paramters =======//
    char name_Mras[500];
    sprintf(name_Mras, param[1].c_str()); // optimized
    cout<<"RHRS Raster parameters file: "<<name_Mras<<endl;
  ifstream Mras(name_Mras);

  for (int i=0;i<nParamT_ras;i++){

    Mras >> Pras[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) Pras[i] = 0.0;
  }


  Mras.close();    

  
  //===== RHRS theta parameters ======// 
    char name_Mxpt[500];
    sprintf(name_Mxpt, param[2].c_str()); // optimized
  ifstream Mxpt(name_Mxpt);

    for(int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt >> par >> p >> p >> p >> p >> p;
    Pxpt[i]  = par;
  }
  Mxpt.close();  
 //===== RHRS phi parameters ======//
  char name_Mypt[500];
    sprintf(name_Mypt, param[3].c_str()); // optimized  
  ifstream Mypt(name_Mypt);

  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mypt >> par >> p >> p >> p >> p >> p;
    Pypt[i]  = par;
  }
  Mypt.close();    



  
};
//////////////////////////////////////////////////////////////

void momcalib::MTParam_L(){

  //=================//
  //===== LHRS ======//
  //=================//

  
  //====== LHRS z parameters ======//  
  char name_Mzt_L[500];
  sprintf(name_Mzt_L,param[4].c_str()); // optimized
  ifstream Mzt_L(name_Mzt_L); 
  for (int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt_L >> par >> p >> p >> p >> p;
    Pzt_L[i]=par;
  }
 Mzt_L.close();

  //====== LHRS raster paramters =======//
    char name_Mras_L[500];
    sprintf(name_Mras_L, param[5].c_str()); // optimized
    cout<<"LHRS Raster parameters file: "<<name_Mras_L<<endl;
  ifstream Mras_L(name_Mras_L);

  for (int i=0;i<nParamT_ras;i++){

    Mras_L >> Pras_L[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) Pras_L[i] = 0.0;
  }
  
  Mras_L.close();    

 
 //===== LHRS theta parameters ======// 
  char name_Mxpt_L[500];
    sprintf(name_Mxpt_L, param[6].c_str()); // optimized
  ifstream Mxpt_L(name_Mxpt_L);
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt_L >> par >> p >> p >> p >> p >> p;
    Pxpt_L[i]  = par;
  }
  Mxpt_L.close();
 //===== LHRS phi parameters ===x==//
  char name_Mypt_L[500];
    sprintf(name_Mypt_L, param[7].c_str()); // optimized
    cout<<"LHRS phi parameters file: "<<name_Mypt_L<<endl;
  ifstream Mypt_L(name_Mypt_L);

  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mypt_L >> par >> p >> p >> p >> p >> p;
    Pypt_L[i]  = par;
  }
  Mypt_L.close();    


}



//////////////////////////////////////////////////////////////////
void momcalib::EventSelection(){


  ENum= T->GetEntries();
  int TNum=300000;
  cout<<"Tuning Events :"<<TNum<<endl;
     if(TNum<10000)d=div(TNum,1000);
   else   d=div(TNum,10000);

  for (int i=0;i<TNum;i++){

    // ===== Initialization ====== //    
    for(int j=0;j<Max;j++){
      Lp[j]=-2222.0;
      Rp[j]=-2222.0;
      Ls2tpads[j]=-2222.0;
      Rs2tpads[j]=-2222.0;
      rtrpathl[j]=-2222.0;
      rs2pathl[j]=-2222.0;
      Rx_fp[j]=-2222.0;
      Ry_fp[j]=-2222.0;
      Rth_fp[j]=-2222.0;
      Rph_fp[j]=-2222.0;
      Rx[j]=-2222.0;
      Ry[j]=-2222.0;
      Rz[j]=-2222.0;
      Rth[j]=-2222.0;
      Rph[j]=-2222.0;
      Lx_fp[j]=-2222.0;
      Ly_fp[j]=-2222.0;
      Lth_fp[j]=-2222.0;
      Lph_fp[j]=-2222.0;
      Lx[j]=-2222.0;
      Ly[j]=-2222.0;
      Lz[j]=-2222.0;
      Lth[j]=-2222.0;
      Lph[j]=-2222.0;
    }
    mm    =-2222.0;
    mm_L  =-2222.0;
    mm_pi =-2222.0;
    Ra1sum=-2222.0;
    Ra2sum=-2222.0;
    Rgssum=-2222.0;
    Lpssum=-2222.0;    
    Lshsum=-2222.0;
    R_Ras_x=-2222.0;
    L_Ras_x=-2222.0;
    //---- Events Flag -----//
    Rs2_flag=false;
    Ls2_flag=false;
    Rs0_flag=false;
    Ls0_flag=false;
    Scinti_flag=false;
    ac1_flag=false;
    ac2_flag=false;
    gs_flag=false;
    sh_flag=false;
    ps_flag=false;
    trig_flag=false;
    track_flag=false;
    Rpathl_flag=false;
    Lpathl_flag=false;
    coin_flag=false;
    x_flag=false;
    y_flag=false;
    z_flag=false;
    RPID_flag=false;
    LPID_flag=false;

    
    T->GetEntry(i);

 pe_=Lp[0]; // Scattered Electron Momentum 
 pk=Rp[0];  // Kaon side Momentum 
 pe=hallap/1000; // GeV


 //------ Set Physics value --------//
 Ee=sqrt(pow(pe,2)+pow(Me,2));
 Ee_=sqrt(pow(pe_,2)+pow(Me,2));
 Ek=sqrt(pow(pk,2)+pow(MK,2));
 Ls2pads=(int)Ls2tpads[0];
 Rs2pads=(int)Rs2tpads[0];
 rpathl=rtrpathl[0]+rs2pathl[0]; // R-HRS path length S2 -RF
 lpathl=ltrpathl[0]+ls2pathl[0]; // L-HRS path length S2 -RF
 rbeta=pk/Ek; 
 rpath_corr=rpathl/rbeta/c;
 lbeta=pe_/Ee_; 
 lpath_corr=lpathl/lbeta/c;


 
 //-------- Tuning value ----------//


      RasterCor= Calc_ras(R_Ras_x, Pras[2], Pras[0]);
      RasterCor = RasterCor/tan(hrs_ang);
      Rz[0] = Rz[0]+RasterCor;

 
 // RHRS //
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Rth[0]    = (Rth[0] - Xptm)/Xptr;
    Rph[0]    = (Rph[0] - Yptm)/Yptr;    
    Rz[0]   = (Rz[0]-Ztm)/Ztr;
    
    
    Rz[0] = calcf2t_zt(Pzt, Rx_fp[0], Rth_fp[0], Ry_fp[0], Rph_fp[0]);
    
    Rth[0]  = calcf2t_ang(Pxpt,
			   Rx_fp[0], Rth_fp[0],
			   Ry_fp[0], Rph_fp[0],
			   Rz[0]);
    Rph[0] = calcf2t_ang(Pypt,
			   Rx_fp[0], Rth_fp[0],
			   Ry_fp[0], Rph_fp[0],
			   Rz[0]);

    Rp[0] = calcf2t_mom(Prp, Rx_fp[i], Ry_fp[i], Rth_fp[i], Rph_fp[i],Rz[i]);

    
    Rph[0]  = Rph[0]*Yptr +Yptm;
    Rth[0]  = Rth[0]*Xptr +Xptm;
    Rz[0]  = Rz[0]*Ztr +Ztm;


    Rx_fp[0]  = Rx_fp[0]  * XFPr + XFPm;
    Rth_fp[0] = Rth_fp[0] * XpFPr + XpFPm;
    Ry_fp[0]  = Ry_fp[0]  * YFPr + YFPm;
    Rph_fp[0] = Rph_fp[0] * YpFPr + YpFPm;


 // LHRS //



      RasterCor_L = Calc_ras(L_Ras_x, Pras_L[2], Pras_L[0]);
      RasterCor_L = RasterCor_L/tan(hrs_ang);
      Lz[0] = Lz[0]+RasterCor_L;




    Lx_fp[0]  = (Lx_fp[0]-XFPm)/XFPr;
    Lth_fp[0] = (Lth_fp[0]-XpFPm)/XpFPr;
    Ly_fp[0]  = (Ly_fp[0]-YFPm)/YFPr;
    Lph_fp[0] = (Lph_fp[0]-YpFPm)/YpFPr;
    Lth[0]    = (Lth[0] - Xptm)/Xptr;
    Lph[0]    = (Lph[0] - Yptm)/Yptr;    
    Lz[0]   = (Lz[0]-Ztm)/Ztr;

    
    Lz[0] = calcf2t_zt(Pzt_L, Lx_fp[0], Lth_fp[0], Ly_fp[0], Lph_fp[0]);
    
    Lth[0]  = calcf2t_ang(Pxpt_L,
			   Lx_fp[0], Lth_fp[0],
			   Ly_fp[0], Lph_fp[0],
			   Lz[0]);
    Lph[0] = calcf2t_ang(Pypt_L,
			   Lx_fp[0], Lth_fp[0],
			   Ly_fp[0], Lph_fp[0],
			   Lz[0]);

    Lp[0] = calcf2t_mom(Prp, Lx_fp[i], Ly_fp[i], Lth_fp[i], Lph_fp[i],Lz[i]);

    
    Lph[0]  = Lph[0]*Yptr +Yptm;
    Lth[0]  = Lth[0]*Xptr +Xptm;
    Lz[0]  = Lz[0]*Ztr +Ztm;


    Lx_fp[0]  = Lx_fp[0]  * XFPr + XFPm;
    Lth_fp[0] = Lth_fp[0] * XpFPr + XpFPm;
    Ly_fp[0]  = Ly_fp[0]  * YFPr + YFPm;
    Lph_fp[0] = Lph_fp[0] * YpFPr + YpFPm;

    
    
 //--------------------------------//




 
 TVector3 L_v, R_v, B_v;
 L_v.SetMagThetaPhi( Lp[0], Lth_fp[0], Lph_fp[0] );
 R_v.SetMagThetaPhi( Rp[0], Rth_fp[0], Rph_fp[0] );
 B_v.SetMagThetaPhi( sqrt(Ee*Ee-Me*Me), 0, 0 );
 L_v.RotateZ( -13.2 / 180. * PI );
 R_v.RotateZ(  13.2 / 180. * PI );

 mm = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
                  - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
 
 //--- Set Coincidence time ---------// 
 Rs2_off=RS2_off_H1[Rs2pads];
 Ls2_off=LS2_off_H1[Ls2pads];
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 coin_t=tof_r-tof_l+rpath_corr-lpath_corr-pathl_off-coin_offset; //  coin Path & Offset  correction   


 
 //---- Event Selection ----------//
  if(coin_t<1.0 && -1.0<coin_t)coin_flag=true;
  if(Rvz_cutmin<Rz[0] && Rz[0]<Rvz_cutmax && Lvz_cutmin<Lz[0] && Lz[0]<Lvz_cutmax) z_flag=true; 
  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads) track_flag=true;
  if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax) Rpathl_flag=true;
  if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax) Lpathl_flag=true;
  if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0) Rs0_flag=true;
  if(-LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0) Ls0_flag=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)   Rs2_flag=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)      Ls2_flag=true;
  if(DRevtype==5)                              trig_flag=true;
  if(Ra1sum>a1_th)      ac1_flag=true;
  if(Ra2sum>a2_th)      ac2_flag=true;
  gs_flag=true;
  ps_flag=true;
  sh_flag=true;
  
  if(ac1_flag==true && ac2_flag==true && gs_flag==true) RPID_flag=true;
  if(ps_flag==true && sh_flag==true)                    LPID_flag=true;
  if(Rs2_flag && Rs0_flag && Ls2_flag && Ls0_flag)      Scinti_flag=true;



    // ------------------------------------------- //
    // ------- Event selection for tuning -------- //
    // ------------------------------------------- //

  if(trig_flag && Scinti_flag && track_flag && Rpathl_flag && Lpathl_flag){
 
  //======== Production Events ============//
    if(coin_flag && RPID_flag && LPID_flag && z_flag){

      mm_L=mm;

      if(ML-mmL_range < mm && mm < ML+mmL_range && i<nmax){
      mass_flag[ntune_event]=0;
      mass_ref[ntune_event] = ML;
      MM[ntune_event]=mm;
      beam_p[ntune_event]=pe;
      rx_fp[ntune_event]  = Rx_fp[0];
      rth_fp[ntune_event] = Rth_fp[0];
      ry_fp[ntune_event]  = Ry_fp[0];
      rph_fp[ntune_event] = Rph_fp[0];
      rx[ntune_event]  = Rx[0];
      rth[ntune_event] = Rth[0];
      ry[ntune_event]  = Ry[0];
      rph[ntune_event] = Rph[0];
      rz[ntune_event] = Rz[0];      
      lx_fp[ntune_event]  = Lx_fp[0];
      lth_fp[ntune_event] = Lth_fp[0];
      ly_fp[ntune_event]  = Ly_fp[0];
      lph_fp[ntune_event] = Lph_fp[0];
      lx[ntune_event]  = Lx[0];
      lth[ntune_event] = Lth[0];
      ly[ntune_event]  = Ly[0];
      lph[ntune_event] = Lph[0];
      lz[ntune_event] = Lz[0];
      Rras_x[ntune_event] =R_Ras_x;
      Lras_x[ntune_event] =L_Ras_x;      
      ntune_event++; 
    }else if(MS0-mmS_range < mm && mm < MS0+mmS_range && i<nmax){
      mass_flag[ntune_event]=1;
      mass_ref[ntune_event] = MS0;
      beam_p[ntune_event]=pe;
      MM[ntune_event]=mm;
      rx_fp[ntune_event]  = Rx_fp[0];
      rth_fp[ntune_event] = Rth_fp[0];
      ry_fp[ntune_event]  = Ry_fp[0];
      rph_fp[ntune_event] = Rph_fp[0];
      rx[ntune_event]  = Rx[0];
      rth[ntune_event] = Rth[0];
      ry[ntune_event]  = Ry[0];
      rph[ntune_event] = Rph[0];
      rz[ntune_event] = Rz[0];      
      lx_fp[ntune_event]  = Lx_fp[0];
      lth_fp[ntune_event] = Lth_fp[0];
      ly_fp[ntune_event]  = Ly_fp[0];
      lph_fp[ntune_event] = Lph_fp[0];
      lx[ntune_event]  = Lx[0];
      lth[ntune_event] = Lth[0];
      ly[ntune_event]  = Ly[0];
      lph[ntune_event] = Lph[0];
      lz[ntune_event] = Lz[0];
      Rras_x[ntune_event] =R_Ras_x;
      Lras_x[ntune_event] =L_Ras_x;            
      ntune_event++;      }



    }
  } // event selection
   


  }//End Fill
  
  }

//========================================================//


void momcalib::MTParam(string mtparam){


  cout<<"==============================="<<endl;
  cout<<"=== Input Matrix Parameters ==="<<endl;
  cout<<"==============================="<<endl;

  string buf;
  int s=0;
  ifstream ifp(Form("%s",mtparam.c_str()),ios::in);
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >>param[s];
    cout<<param[s]<<endl;
    s++;
  }

  if(RHRS) {MTParam_R();cout<<" Input RHRS Matrix parameter "<<endl;}
  if(LHRS) {MTParam_L();cout<<" Input LHRS Matrix parameter "<<endl;}
  MTP_mom();
}


//========================================================//


void momcalib::MomTuning(string ofname){


  if (nite>0){
    cout << "======================================================" <<endl;
    cout << "=================  Tuning started =================== " << endl;
    cout << "======================================================" <<endl;}
    char tempc[500],tempc2[500];
    const  char* new_tempc=ofname.c_str();
    cout<<"new marix file: "<<new_tempc<<endl;
    ofstream * ofs1;
    ofstream * ofs2;       
  for(int i=0 ; i<nite ; i++){

    cout<<"tuning i: "<<i+1<<" /"<<nite<<endl;
    // --------------------------- //
    // ---- Parameter tuning ----- //
    // --------------------------- //
    
    if(MODE==-1){
    cout<<"------- Rp tuning -----"<<endl;
    chi_sq1[i] = tune(Opt_par_R,i,-1);   // Rp
    }else if(MODE==1){
    cout<<"------- Lp tuning -----"<<endl;    
    chi_sq2[i] = tune(Opt_par_L,i,1);  // Lp
    }else if(MODE==0){
    cout<<"---------pk & pe tuning ------"<<endl;
    chi_sq[i] = tune(Opt_par,i,0);
    }
    
    cout << " Tuning# = " << i << ": chisq = ";
    if(MODE==-1) cout  << chi_sq1[i] << endl;
    if(MODE== 1) cout  << chi_sq2[i] << endl;
    if(MODE==0)  cout  << chi_sq[i]  << endl;
    cout << endl;
    gchi_p ->SetPoint(i,i,chi_sq[i]);    
    gchi_Rp->SetPoint(i,i,chi_sq1[i]);    
    gchi_Lp->SetPoint(i,i,chi_sq2[i]);

    if(MODE==-1 || MODE==0){
    sprintf(tempc,  "%s_Rp_%d_MODE%d.dat",new_tempc,i,MODE); 
    cout<<"new matrix Rp: "<<tempc<<endl;
    ofs1 = new ofstream(tempc);}

    if(MODE==1 ||MODE==0){
    sprintf(tempc2, "%s_Lp_%d_MODE%d.dat",new_tempc,i,MODE);
    cout<<"new matrix Lp: "<<tempc2<<endl;    
    ofs2 = new ofstream(tempc2);}

    
    int nppp = 0;
    for(int i=0 ; i<nn+1 ; i++){
      for(int e=0 ; e<nn+1 ; e++){
	for(int d=0 ; d<nn+1 ; d++){
	  for(int c=0 ; c<nn+1 ; c++){
	    for(int b=0 ; b<nn+1 ; b++){
	      for(int a=0 ; a<nn+1 ; a++){  
		if(a+b+c+d+e==i){
		  if(MODE==-1){		  
		  *ofs1 << Opt_par_R[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  }else if(MODE==1){
		  *ofs2 << Opt_par_L[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  }else if(MODE==0){
		  *ofs1 << Opt_par_R[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  *ofs2 << Opt_par_L[nppp+nParamTp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  }
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

    if(MODE==-1 || MODE==0){
      ofs1->close();
      ofs1->clear();
    }else if(MODE==0 || MODE==1){
      ofs2->close();
      ofs2->clear();}


    
  }

  
  cout<<"========== Tuning is done ============="<<endl;
  
  

}

//========================================================//

double momcalib::Eloss(double xp, double z,char* arm){


  
  double x = - tan(hrs_ang-xp);
  double ph[3],pl[2];
  double dEloss;
  double dEloss_h = ph[0]*sin(ph[1]*x)+ph[2];
  double dEloss_l = pl[0]*x +pl[1];
  bool high;
  
  if(z>0.08)high=false;
  else high=true;
  
    //==== thickness 0.400 mm ========//

  if(arm=="R"){
    ph[0] = -1.3175;
    ph[1] = -4.6151;
    ph[2] = 2.0369;
    pl[0] = 3.158e-2;
    pl[1] = 4.058e-1;
  }else if(arm=="L"){
    ph[0] = -1.3576;
    ph[1] = -4.5957;
    ph[2] = 2.0909;
    pl[0] = 6.2341e-3;
    pl[1] = 4.0336e-1;
  }

  if(high)dEloss = dEloss_h;
  else dEloss = dEloss_l;
  
  return dEloss;
  
}

//========================================================//

void momcalib::Fill(){


  ENum= T->GetEntries();
  cout<<"Events :"<<ENum<<endl;
     if(ENum<10000)d=div(ENum,1000);
   else   d=div(ENum,10000);

  for (int i=0;i<ENum;i++){

    // ===== Initialization ====== //    
    for(int j=0;j<Max;j++){
      Lp[j]=-2222.0;
      Rp[j]=-2222.0;
      Ls2tpads[j]=-2222.0;
      Rs2tpads[j]=-2222.0;
      rtrpathl[j]=-2222.0;
      rs2pathl[j]=-2222.0;
      Rx_fp[j]=-2222.0;
      Ry_fp[j]=-2222.0;
      Rth_fp[j]=-2222.0;
      Rph_fp[j]=-2222.0;
      Rx[j]=-2222.0;
      Ry[j]=-2222.0;
      Rth[j]=-2222.0;
      Rph[j]=-2222.0;
      Lx_fp[j]=-2222.0;
      Ly_fp[j]=-2222.0;
      Lth_fp[j]=-2222.0;
      Lph_fp[j]=-2222.0;
      Lx[j]=-2222.0;
      Ly[j]=-2222.0;
      Lth[j]=-2222.0;
      Lph[j]=-2222.0;
    }
    mm    =-2222.0;
    mm_L  =-2222.0;
    mm_pi =-2222.0;
    Ra1sum=-2222.0;
    Ra2sum=-2222.0;
    Rgssum=-2222.0;
    Lpssum=-2222.0;    
    Lshsum=-2222.0;
    
    //---- Events Flag -----//
    Rs2_flag=false;
    Ls2_flag=false;
    Rs0_flag=false;
    Ls0_flag=false;
    Scinti_flag=false;
    ac1_flag=false;
    ac2_flag=false;
    gs_flag=false;
    sh_flag=false;
    ps_flag=false;
    trig_flag=false;
    track_flag=false;
    Rpathl_flag=false;
    Lpathl_flag=false;
    coin_flag=false;
    x_flag=false;
    y_flag=false;
    z_flag=false;
    RPID_flag=false;
    LPID_flag=false;

    
    T->GetEntry(i);

 pe_=Lp[0]; // Scattered Electron Momentum 
 pk=Rp[0];  // Kaon side Momentum 
 pe=hallap/1000; // GeV


 //------ Set Physics value --------//
 Ee=sqrt(pow(pe,2)+pow(Me,2));
 Ee_=sqrt(pow(pe_,2)+pow(Me,2));
 Ek=sqrt(pow(pk,2)+pow(MK,2));
 Ls2pads=(int)Ls2tpads[0];
 Rs2pads=(int)Rs2tpads[0];
 rpathl=rtrpathl[0]+rs2pathl[0]; // R-HRS path length S2 -RF
 lpathl=ltrpathl[0]+ls2pathl[0]; // L-HRS path length S2 -RF
 rbeta=pk/Ek; 
 rpath_corr=rpathl/rbeta/c;
 lbeta=pe_/Ee_; 
 lpath_corr=lpathl/lbeta/c;











 
 TVector3 L_v, R_v, B_v;
 L_v.SetMagThetaPhi( Lp[0], Lth_fp[0], Lph_fp[0] );
 R_v.SetMagThetaPhi( Rp[0], Rth_fp[0], Rph_fp[0] );
 B_v.SetMagThetaPhi( sqrt(Ee*Ee-Me*Me), 0, 0 );
 L_v.RotateZ( -13.2 / 180. * PI );
 R_v.RotateZ(  13.2 / 180. * PI );

 mm = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
                  - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
 
 //--- Set Coincidence time ---------// 
 Rs2_off=RS2_off_H1[Rs2pads];
 Ls2_off=LS2_off_H1[Ls2pads];
 tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
 tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
 coin_t=tof_r-tof_l+rpath_corr-lpath_corr-pathl_off-coin_offset; //  coin Path & Offset  correction   


 
 //---- Event Selection ----------//
  if(coin_t<1.0 && -1.0<coin_t)coin_flag=true;
  if(Rvz_cutmin<Rz[0] && Rz[0]<Rvz_cutmax && Lvz_cutmin<Lz[0] && Lz[0]<Lvz_cutmax) z_flag=true; 
  if(Rs2trpad[0]==Rs2pads && Ls2trpad[0]==Ls2pads) track_flag=true;
  if(rpathl_cutmin<rpathl && rpathl<rpathl_cutmax) Rpathl_flag=true;
  if(lpathl_cutmin<lpathl && lpathl<lpathl_cutmax) Lpathl_flag=true;
  if(-RF1[43]+RF1[46]>0 && -RF1[44]+RF1[46]>0) Rs0_flag=true;
  if(-LF1[27]+LF1[30]>0 && -LF1[28]+LF1[30]>0) Ls0_flag=true;
  if(RF1[48+Rs2pads]>0 && RF1[16+Rs2pads]>0)   Rs2_flag=true;
  if(LF1[Ls2pads]>0 && LF1[Ls2pads+48]>0)      Ls2_flag=true;
  if(DRevtype==5)                              trig_flag=true;
  if(Ra1sum>a1_th)      ac1_flag=true;
  if(Ra2sum>a2_th)      ac2_flag=true;
  gs_flag=true;
  ps_flag=true;
  sh_flag=true;
  
  if(ac1_flag==true && ac2_flag==true && gs_flag==true) RPID_flag=true;
  if(ps_flag==true && sh_flag==true)                    LPID_flag=true;
  if(Rs2_flag && Rs0_flag && Ls2_flag && Ls0_flag)      Scinti_flag=true;



  
  if(trig_flag && Scinti_flag && track_flag && Rpathl_flag && Lpathl_flag){
 
  //======== Production Events ============//
    if(coin_flag && RPID_flag && LPID_flag && z_flag)mm_L=mm;

  //======== ACCidencel B.G. ==============///
    if(RPID_flag && LPID_flag && z_flag &&
       (-63<coin_t && coin_t <-15) || (15<coin_t && coin_t<63)){
      ct=coin_t;

              while(1){
	       if(-3.0<ct && ct<3.0){ctime=ct; break;}
	       else if(ct<-3.0) ct=ct+6;
	       else if(3.0<ct)  ct=ct-6;
	      }
    }
    //=====================================//

    
  }

  
  tnew->Fill();
          if(i % (d.quot * 1000) == 0)cout<<i<<" / "<<ENum<<endl;
  }


  
}

//========================================================//


//=======================================================//

void momcalib::Close(){

   tnew->Write();
   ofp->Close();

}







//=========================================================//
//================= Function ==============================//
//=========================================================//

// ####################################################
double s2f1_off(int i,char* ARM,char* MODE, int KINE){
// ####################################################

  double RS2_offset[16],LS2_offset[16];
  if(MODE=="H" && KINE==2){
 
 double  RS2_off_H2[16]={-16911.4,-16864.3,-16900,-16897,-16873.8,-16868.4,-16901.1,-16876.8,-16895.4,-16860.9,-16893.1,-16884.4,-16847.3,-16842.7,-16836.9,-16882.6};
 double  LS2_off_H2[16]={-25336.9,-25386.6,-25367.5,-25392.3,-25391.1,-25386.2,-25422,-25428.9,-25417.3,-25426.8,-25438.7,-25383.4,-25396,-25418.5,-25436.4,-26082.1};
 
  LS2_offset[i]=LS2_off_H2[i];
  RS2_offset[i]=RS2_off_H2[i];
  }


  if(MODE=="H" && KINE==1){

double  RS2_off_H1[16]={-16828.7,-16863,-16894,-16893.3,-16870.9,-16867.2,-16900.3,-16876.8,17554.1,-16861.6,-16895,-16890.7,-16854.6,-16852.9,-16850.5,-16861.9};
double  LS2_off_H1[16]={-25335,-25385.6,-25367,-25392.1,-25391.7,-25386.4,-25422.1,-25428.9,-25414.9,-25424.7,-25436.9, -25381.2,-25390,-25413.4,-25428.7,-26640.8};
  LS2_offset[i]=LS2_off_H1[i];
  RS2_offset[i]=RS2_off_H1[i];
  }

 double s2f1_offset; 
 if(ARM=="R")s2f1_offset=RS2_offset[i];
 else  if(ARM=="L")s2f1_offset=LS2_offset[i];
 else {cout<<"false read out !!"<<endl;}

  return s2f1_offset;

}


// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  
  const double sigma = 1.0;
  double ztR      = 0.0;
  double refpos   = 0.0;
  double residual = 0.0;
  double ang      = 0.0;
  double sspos    = 0.0;
  double total_chi2 = 0.0;
  double mass;
  double Ee,Ee_,Ek;
  double rbeta,lbeta;
  double rp,lp;
  double chi2;
  double rx_c,ry_c,rth_c,rph_c,rz_c,lx_c,ly_c,lth_c,lph_c;
  for(int i=0 ; i<ntune_event ; i++){

    //====== Initialization ========//
    rp=0.0;
    lp=0.0;
    residual =0.0;


    // RHRS //

      RasterCor= Calc_ras(Rras_x[i], Pras[2], Pras[0]);
      RasterCor = RasterCor/tan(hrs_ang);
      rz[i] = rz[i]+RasterCor;
      
    //======== Scaled paramters ===============//

    rx_fp[i]  = (rx_fp[i]-XFPm)/XFPr;
    rth_fp[i] = (rth_fp[i]-XpFPm)/XpFPr;
    ry_fp[i]  = (ry_fp[i]-YFPm)/YFPr;
    rph_fp[i] = (rph_fp[i]-YpFPm)/YpFPr;
    rth[i]    = (rth[i] - Xptm)/Xptr;
    rph[i]    = (rph[i] - Yptm)/Yptr;    
    rz[i]     = (rz[i]-Ztm)/Ztr;

    //======== Tuned parameters ==============///


    rz[i] = calcf2t_zt(Pzt, rx_fp[i], rth_fp[i], ry_fp[i], rph_fp[i]);
    
    rth[i]  = calcf2t_ang(Pxpt,
			   rx_fp[i], rth_fp[i],
			   ry_fp[i], rph_fp[i],
			   rz[i]);
    rph[i] = calcf2t_ang(Pypt,
			   rx_fp[i], rth_fp[i],
			   ry_fp[i], rph_fp[i],
			   rz[i]);

    //==== momentum tuning ==========//

    rp = calcf2t_mom(Opt_par_R, rx_fp[i], ry_fp[i], rth_fp[i], rph_fp[i],rz[i]);



    //==== Scaled to Nomal =======//
    
    rth[i]  = rth[i]*Xptr +Xptm;
    rph[i]  = rph[i]*Yptr +Yptm;
    rz[i]   = rz[i]*Ztr +Ztm;
    rp      = rp   *Momr +Momm;
    
    rx_fp[i]  = rx_fp[i]  * XFPr + XFPm;
    rth_fp[i] = rth_fp[i] * XpFPr + XpFPm;
    ry_fp[i]  = ry_fp[i]  * YFPr + YFPm;
    rph_fp[i] = rph_fp[i] * YpFPr + YpFPm;
    //=============================//
   

    // LHRS //

      RasterCor_L= Calc_ras(Lras_x[i], Pras_L[2], Pras_L[0]);
      RasterCor_L = RasterCor_L/tan(hrs_ang);
      lz[i] = lz[i]+RasterCor_L;
    
    //======== Scaled paramters ===============//

    lx_fp[i]  = (lx_fp[i]-XFPm)/XFPr;
    lth_fp[i] = (lth_fp[i]-XpFPm)/XpFPr;
    ly_fp[i]  = (ly_fp[i]-YFPm)/YFPr;
    lph_fp[i] = (lph_fp[i]-YpFPm)/YpFPr;
    lth[i]    = (lth[i] - Xptm)/Xptr;
    lph[i]    = (lph[i] - Yptm)/Yptr;    
    lz[i]     = (lz[i]-Ztm)/Ztr;

    //======== Tuned parameters ==============///

    lz[i] = calcf2t_zt(Pzt_L, lx_fp[i], lth_fp[i], ly_fp[i], lph_fp[i]);
    
    lth[i]  = calcf2t_ang(Pxpt_L,
			   lx_fp[i], lth_fp[i],
			   ly_fp[i], lph_fp[i],
			   lz[i]);
    lph[i] = calcf2t_ang(Pypt_L,
			   lx_fp[i], lth_fp[i],
			   ly_fp[i], lph_fp[i],
			   lz[i]);

    //==== momentum tuning ==========//
    
    lp= calcf2t_mom(Opt_par_L, lx_fp[i], ly_fp[i], lth_fp[i], lph_fp[i], lz[i]);


    //==== Scaled to Nomal =======//
    
    lth[i]  = lth[i]*Xptr +Xptm;
    lph[i]  = lph[i]*Yptr +Yptm;
    lz[i]  = lz[i]*Ztr +Ztm;
    lp      = lp   *Momr +Momm;
    
    lx_fp[i]  = lx_fp[i]  * XFPr + XFPm;
    lth_fp[i] = lth_fp[i] * XpFPr + XpFPm;
    ly_fp[i]  = ly_fp[i]  * YFPr + YFPm;
    lph_fp[i] = lph_fp[i] * YpFPr + YpFPm;

    //=============================//
   
    

    


    
 //------ Set Physics value --------//
 Ee=sqrt(pow(beam_p[i],2)+pow(Me,2));
 Ee_=sqrt(pow(lp,2)+pow(Me,2));
 Ek=sqrt(pow(rp,2)+pow(MK,2));


 rbeta=rp/Ek; 
 lbeta=lp/Ee_; 


 TVector3 L_v, R_v, B_v;
 L_v.SetMagThetaPhi( lp, lth_fp[i], lph_fp[i] );
 R_v.SetMagThetaPhi( rp, rth_fp[i], rph_fp[i] );
 B_v.SetMagThetaPhi( sqrt(Ee*Ee-Me*Me), 0, 0 );
 L_v.RotateZ( -13.2 / 180. * PI );
 R_v.RotateZ(  13.2 / 180. * PI );

 mass = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
                  - (B_v - L_v - R_v)*(B_v - L_v - R_v) );



    
    //======================//
    //====== Residual ======//
    //======================//



    
    residual = mass - mass_ref[i];
    chi2=chi2 + pow(residual,2.0);
    

    
  }

 
}//end fcn



// #############################################################
double momcalib::tune(double* pa, int j, int MODE) 
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10]; 
  int ierflg = 0;
  int allparam = nParamTp;
  int arm=1;
  if(MODE==0){allparam=nParamTp*2;
    arm=2;
  }
  TMinuit* minuit = new TMinuit(allparam);
    minuit->SetFCN(fcn); // fcn Chi-square function


  double start[allparam];
  double step[allparam];
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=3; // The number of order is reduced for test (4-->2)
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;

  for(int f=0;f<arm;f++){
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
    if(MODE==0){minuit -> GetParameter(i,Opt_par[i],er);
      //      if(i<allparam/2)minuit -> GetParameter(i,OptPar_R[i],er);// RHRS momentum parameters
      //      else minuit -> GetParameter(i,OptPar_L[i],er); // LHRS momentum parameters
    }else if(MODE==-1){minuit -> GetParameter(i,Opt_par_R[i],er);// RHRS momentum parameters
    }else if(MODE==1)minuit -> GetParameter(i,Opt_par_L[i],er); // LHRS momentum parameters

  }
  
  return chi2;
}


// ###################################################
double calcf2t_mom(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt)
// ###################################################
{
  // ----------------------------------------------------------------- //
  // ------ 4rd order using xf, xpf, yf, ypf, zt, xt, xpt, yt, ytp --- //
  // ----------------------------------------------------------------- //

  const int nMatT=nn;  
  const int nXf=nn;
  const int nXpf=nn;
  const int nYf=nn;
  const int nYpf=nn;
  const int nZt=nn;
  const int nXt=nn;
  const int nXpt=nn;
  const int nYt=nn;
  const int nYpt=nn;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0,f=0,g=0,h=0,i=0;
  for (int n=0;n<nMatT+1;n++){
	    for(e=0;e<n+1;e++){
	      for (d=0;d<n+1;d++){
		for (c=0;c<n+1;c++){ 
		  for (b=0;b<n+1;b++){
		    for (a=0;a<n+1;a++){ 
		      if (a+b+c+d+e+f+g+h+i==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt
		    && f<=nXt && g<=nXpt && h<=nYt && i<=nYpt){
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



// ###################################################
double calcf2t_mom_RL(double* RP, double Rxf, double Rxpf, 
		      double Ryf, double Rypf, double Rzt,
		      double* LP, double Lxf, double Lxpf, 
		      double Lyf, double Lypf, double Lzt)		      
// ###################################################
{
  // ------------------------------------------------------------------------------ //
  // ------ 4rd order using xf, xpf, yf, ypf, zt, xt, xpt, yt, ytp IN Both HRS ---- //
  // ------------------------------------------------------------------------------ //

  const int nMatT=nn;  
  const int nXf=nn;
  const int nXpf=nn;
  const int nYf=nn;
  const int nYpf=nn;
  const int nZt=nn;
  const int nXt=nn;
  const int nXpt=nn;
  const int nYt=nn;
  const int nYpt=nn;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0,f=0,g=0,h=0,i=0,j=0;
  for (int n=0;n<nMatT+1;n++){
    for (j=0;j<n+1;j++){ 
      for (i=0;i<n+1;i++){ 
	for (h=0;h<n+1;h++){
	  for (g=0;g<n+1;g++){ 
	    for (f=0;f<n+1;f++){ 
	      for(e=0;e<n+1;e++){
		for (d=0;d<n+1;d++){
		  for (c=0;c<n+1;c++){ 
		    for (b=0;b<n+1;b++){
		      for (a=0;a<n+1;a++){ 
			if (a+b+c+d+e+f+g+h+i+j==n){
			  if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt
			      && f<=nXt && g<=nXpt && h<=nYt && i<=nYpt){
			    x = pow(Rxf,double(a))*pow(Rxpf,double(b))*
			      pow(Ryf,double(c))*pow(Rypf,double(d))*pow(Rzt,double(e))*
			      pow(Lxf,double(f))+pow(Lxpf,double(g))*pow(Lyf,double(h))*pow(Lypf,double(i))*pow(Lzt,double(j));
			  }
			  else{
			    x = 0.;
			  }
			  Y += x*RP[npar]; 
			  npar++;
			}
		      }
		    }
		  }
		}    
	      }
	    }
	  }
	}
      }
    }
  }
  return Y; 
  
}



// ###################################################
double calcf2t_ang(double* P, double xf, double xpf, 
		     double yf, double ypf, double zt)
// ####################################################
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
      	

// #################################################
double calcf2t_zt(double* P, double xf, double xpf, 
                 double yf, double ypf){
// ###############################################

  int nnz=3;  

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







#endif
