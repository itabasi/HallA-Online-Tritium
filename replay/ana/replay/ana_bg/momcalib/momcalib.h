#ifndef momcalib_h
#define momcalib_h 1

#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"
#include "define.h"
#include "tree.h"



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
  int Mode(string arm);
  void Mzt(string ifname);
  void Mxpt(string ifname);
  void Mypt(string ifname);
  // void Mpt(string ifname);
  void Fill();
  void Close();
  //------ NewROOT -----------//
  TFile* ofp;
  TTree* tnew;


  //------ EventSelection ---//
  div_t d;
  
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
    

};

//=====================================//
//======== Momcalib ==================//
//=====================================//

momcalib::momcalib(){};
momcalib::~momcalib(){};


int momcalib::mode(string arm){
  int MODE=5;
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

 // RHRS //
    Rx_fp[0]  = (Rx_fp[0]-XFPm)/XFPr;
    Rth_fp[0] = (Rth_fp[0]-XpFPm)/XpFPr;
    Ry_fp[0]  = (Ry_fp[0]-YFPm)/YFPr;
    Rph_fp[0] = (Rph_fp[0]-YpFPm)/YpFPr;
    Rz[0]   = (Rz[0]-Ztm)/Ztr;

    Rth[0]  = calcf2t_4th_2(Pxpt,
			   Rx_fp[0], Rth_fp[0],
			   Ry_fp[0], Rph_fp[0],
			   Rz[0]);
    Rph[0] = calcf2t_4th_2(Pypt,
			   Rx_fp[0], Rth_fp[0],
			   Ry_fp[0], Rph_fp[0],
			   Rz[0]);

    Rz[0] = calcf2t_zt(Pzt, Rx_fp[0], Ry_fp[0], YFP, YpFP);
    
    Ypt[0]  = Ypt[0]*Yptr +Yptm;
    Xpt[0]  = Xpt[0]*Xptr +Xptm;
    Zt[0]   = Zt[0]*Ztr +Ztm;


      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;

 // LHRS //
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

    ztR[0] = calcf2t_zt(Pzt, XFP, XpFP, YFP, YpFP);
    
    Ypt[0]  = Ypt[0]*Yptr +Yptm;
    Xpt[0]  = Xpt[0]*Xptr +Xptm;
    Zt[0]   = Zt[0]*Ztr +Ztm;


      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;      
    
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
    if(coin_flag && RPID_flag && LPID_flag && z_flag)mm_L=mm;
  }


  }//End Fill
  


  
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
  

  for(int i=0 ; i<ntune_event ; i++){

   // ===== Initialization ====== //    
      for(int j=0;j<MAX;j++){
	if(RHRS){
	  //	  R_tr_vx[j]  = -2222.0;
	  //	  R_tr_vy[j]  = -2222.0;
	  R_tr_vz[j]  = -2222.0;	  
	  R_tr_tg_th[j] = -2222.0;
	  R_tr_tg_ph[j] = -2222.0;
	}
	if(LHRS){
	  //          L_tr_vx[j]  = -2222.0;
	  //	  L_tr_vy[j]  = -2222.0;
	  L_tr_vz[j]  = -2222.0;
	  L_tr_tg_th[j] = -2222.0;
	  L_tr_tg_ph[j] = -2222.0;
	}
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



   tree->GetEntry(i);

   //======== w/o Tuning Hist =======//
   //	hRvz->Fill(R_tr_vz[0]);
   //	hRth->Fill(R_tr_tg_th[0]);
   //	hRph->Fill(R_tr_tg_ph[0]);	
   //	hLvz->Fill(L_tr_vz[0]);
   //	hLth->Fill(L_tr_tg_th[0]);
   //	hLph->Fill(L_tr_tg_ph[0]);
	

    
      // ----------------------------------------------------------- //
      // --------------------- RHRS ------------------------ ------- //
      // ----------------------------------------------------------- //
    

    if(RHRS){
      R_tr_x[0]  = (R_tr_x[0]-XFPm)/XFPr;
      R_tr_th[0] = (R_tr_th[0]-XpFPm)/XpFPr;
      R_tr_y[0]  = (R_tr_y[0]-YFPm)/YFPr;
      R_tr_ph[0] = (R_tr_ph[0]-YpFPm)/YpFPr;

      R_tr_vz[0] = calcf2t_zt(Pzt, R_tr_x[0], R_tr_th[0], R_tr_y[0], R_tr_ph[0]);
      R_tr_vz[0] = R_tr_vz[0] * Ztr + Ztm;

      //-------w/ Tuning z calib ----------//
      //      hRvz_c->Fill(R_tr_vz[0]);

      // ----------------------------------------------------------- //
      // ------ Converting raster current to raster position ------- //
      // ----------------------------------------------------------- //


      double RasterCor = calcRasterCor(R_Ras_x, parRaster_R_2, parRaster_R_0);
      RasterCor = RasterCor/tan(hrs_ang);
      //      R_tr_vz[0]=R_tr_vz[0]+RasterCor; w/ Raster correction
      //      hRvz_Rc->Fill(R_tr_vz[0]);


      
      double Zt=R_tr_vz[0];
      Zt= (Zt-Ztm)/Ztr;
      
      R_tr_tg_th[0]  =  calcf2t_ang(Pxpt, R_tr_x[0], R_tr_th[0], R_tr_y[0], R_tr_ph[0],Zt);
      R_tr_tg_ph[0]  =  calcf2t_ang(Pypt, R_tr_x[0], R_tr_th[0], R_tr_y[0], R_tr_ph[0],Zt);
      R_tr_tg_th[0]  =  R_tr_tg_th[0]*Xptr +Xptm;
      R_tr_tg_ph[0]  =  R_tr_tg_ph[0]*Yptr +Yptm;

      R_tr_x[0]  = R_tr_x[0]  * XFPr  + XFPm;
      R_tr_y[0]  = R_tr_y[0]  * YFPr  + YFPm;
      R_tr_th[0] = R_tr_th[0] * XpFPr + XpFPm;
      R_tr_ph[0] = R_tr_ph[0] * YpFPr + YpFPm;

      //      	hRth_c->Fill(R_tr_tg_th[0]);
      //	hRph_c->Fill(R_tr_tg_ph[0]);


    }//end RHRS Fill
    


	
      // ----------------------------------------------------------- //
      // --------------------- LHRS ------------------------ ------- //
      // ----------------------------------------------------------- //
	

    if(LHRS){
      L_tr_x[0]  = (L_tr_x[0]-XFPm)/XFPr;
      L_tr_th[0] = (L_tr_th[0]-XpFPm)/XpFPr;
      L_tr_y[0]  = (L_tr_y[0]-YFPm)/YFPr;
      L_tr_ph[0] = (L_tr_ph[0]-YpFPm)/YpFPr;
      L_tr_vz[0] =  calcf2t_zt(Pzt_L, L_tr_x[0], L_tr_th[0], L_tr_y[0], L_tr_ph[0]);
      L_tr_vz[0] = L_tr_vz[0] * Ztr + Ztm;


      //-------w/ Tuning z calib ----------//
      //      hLvz_c->Fill(L_tr_vz[0]);

      // ----------------------------------------------------------- //
      // ------ Converting raster current to raster position ------- //
      // ----------------------------------------------------------- //

      double RasterCor_L = calcRasterCor(L_Ras_x, parRaster_L_2, parRaster_L_0);
      RasterCor_L = RasterCor_L/tan(hrs_ang);
      //      L_tr_vz[0]=L_tr_vz[0]+RasterCor;
      //      hLvz_Rc->Fill(L_tr_vz[0]);

      double Zt_L=L_tr_vz[0];
      Zt_L= (Zt_L-Ztm)/Ztr;
      L_tr_tg_th[0]  =  calcf2t_ang(Pxpt_L, L_tr_x[0], L_tr_th[0], L_tr_y[0], L_tr_ph[0],Zt_L);
      L_tr_tg_ph[0] =   calcf2t_ang(Pypt_L, L_tr_x[0], L_tr_th[0], L_tr_y[0], L_tr_ph[0],Zt_L);
      L_tr_tg_th[0]  = L_tr_tg_th[0]*Xptr +Xptm;
      L_tr_tg_ph[0]  = L_tr_tg_ph[0]*Yptr +Yptm;

      L_tr_x[0]  = L_tr_x[0]  * XFPr  + XFPm;
      L_tr_y[0]  = L_tr_y[0]  * YFPr  + YFPm;
      L_tr_th[0] = L_tr_th[0] * XpFPr + XpFPm;
      L_tr_ph[0] = L_tr_ph[0] * YpFPr + YpFPm;


      //	hLth_c->Fill(L_tr_tg_th[0]);
      //	hLph_c->Fill(L_tr_tg_ph[0]);	

    }//End LHRS Fill










    
    //    tnew->Fill();
    if(i % (d.quot*1000) == 0)cout<<i<<" / "<<ENum<<endl;
  }

}//end fcn




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
    minuit->SetFCN(fcn1); // fcn1 Chi-square function
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







#endif
