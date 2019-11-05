//////////////////////////////////////
//t0 calibration of S2(F1TDC & Fbus)//
// by K. Itabashi Oct. 2019           //
//////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <iomanip>
#include <csignal>
#include <stdlib.h>
#include <climits>
using namespace std;
#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "Tree.h"
#include "Setting.h"
#include "ParamMan.h"
#define Calibration


const int NCanvas = 18;//num of canvas
class coin_t0_calib : public Tree
{
public:
  coin_t0_calib();
  ~coin_t0_calib();
  void makehist();
  void loop();
  void fit();
  void draw(); 
  void savecanvas(string ofname); 
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetRoot(string ifname);
  void SetInputParam(string ifname);
  void SetLR(int lr){LR=lr;}
  void OutParameter(string ifname, string ofname);
  double  coin_time(int Rs2_seg, int Ls2_seg); 
private:
  int GetMaxEvent() { return ENumMax; }
  int ENumMax;
  bool anaL_oneevent();
  bool anaR_oneevent();
  
  
  TH1F *h_s2s0_beta[16], *h_s2s0_beta_FB[16];
  TH1F *h_s2s0_tof[16], *h_s2s0_tdiff[16], *h_s2s0_tof_FB[16], *h_s2s0_tdiff_FB[16];
  bool coin_OK=false;
  TH2F *h2_s2s0_f1beta_fbbeta[16], *h2_s2s0_f1tof_fbtof[16];
  TF1 *ga_tdiffF1[16], *ga_tdiffFB[16];
  TF1* fcoin_R[16];
  TF1* fcoin_L[16];
  TF1* fcoin;
  TH2F *h_frame[4];
  TH1F* hcoin;
  TH1F* hcoin_c;
  TH1F* hcoin_Rs2_t[16];
  TH1F* hcoin_Ls2_t[16];
  TH2F* hcoin_Rs2_ra[16];
  TH2F* hcoin_Rs2_la[16];
  TH2F* hcoin_Ls2_ra[16];
  TH2F* hcoin_Ls2_la[16];
  TH2F* hcoin_Rp;
  TH2F* hcoin_Lp;
  int LR;//L = 0, R = 1
  int run_num;
  double coin_F1_offset;
  double F1reso;
  Setting *set;
  ParamMan *param;
  TCanvas *c[NCanvas];
  TGraphErrors *tg_tdiffFB_pos, *tg_tdiffFB_wid;
  TGraphErrors *tg_tdiffF1_pos, *tg_tdiffF1_wid;
  TGraphErrors *tg_coin_R_pos, *tg_coin_L_pos;
  TGraphErrors *tg_coin_R_wid, *tg_coin_L_wid;
  double tdiffFB_pos[16], tdiffFB_wid[16], etdiffFB_pos[16], etdiffFB_wid[16];
  double tdiffF1_pos[16], tdiffF1_wid[16], etdiffF1_pos[16], etdiffF1_wid[16];
  
  int tr_n;//num. of track
  double beta[MAX],betaF1[MAX],s2_trpad[MAX],paths2s0[MAX];
  double S2T_F1TDC[16],S0T_F1TDC[1];
  double S2B_F1TDC[16],S0B_F1TDC[1];
  double S2_F1time[16],S0_F1time[1];
  double S2_lt[16],S2_rt[16],S0_lt[1],S0_rt[1];
  double S2_time[16],S0_time[1];
  double coin_pos,coin_wid;
  double coin_R[16],coin_L[16];  
  double coin_Rw[16],coin_Lw[16];
  
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
coin_t0_calib::coin_t0_calib()
{

  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->SetBatch(1);

  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetOptDate(0);
//  gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);
  gStyle->SetGridWidth(1);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first
  gStyle->SetOptStat("iMenr");
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);

  for(int i=0;i<NCanvas;i++){
    c[i]= new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1400,800 );
  }

  set = new Setting();
}
////////////////////////////////////////////////////////////////////////////
coin_t0_calib::~coin_t0_calib(){
}
////////////////////////////////////////////////////////////////////////////
void coin_t0_calib::SetRoot(string ifname){
  add_tree(ifname);
  pack_tree();

    readtreeHRSL();
    readtreeTrackL();

      readtreeHRSR();
      readtreeTrackR();
      readtreeA1R();
      readtreeA2R();
}
////////////////////////////////////////////////////////////////////////////
void coin_t0_calib::SetInputParam(string ifname){
  param = new ParamMan(ifname.c_str());
  if(param -> SetVal())cout<<"F1TDC parameter setted : really cool acutually"<<endl;
  coin_F1_offset=param->GetF1CoinOffset();
  F1reso=param->GetF1reso();
}
////////////////////////////////////////////////////////////////////////////
void coin_t0_calib::makehist(){
  string LorR;
  if(LR==0)     LorR="L";
  else if(LR==1)LorR="R";
  for(int i=0;i<16;i++){
    h_s2s0_beta[i]     = new TH1F(Form("h_s2s0_beta%d",i)    , Form("h_s2s0_beta%d",i)     ,200,  -1,1.5);
    h_s2s0_tof[i]      = new TH1F(Form("h_s2s0_tof%d",i)     , Form("h_s2s0_tof%d",i)      ,800,-100,100);
    h_s2s0_tdiff[i]    = new TH1F(Form("h_s2s0_tdiff%d",i)   , Form("h_s2s0_tdiff%d",i)    , 80,-10,  10);
    h_s2s0_beta_FB[i]  = new TH1F(Form("h_s2s0_beta_FB%d",i) , Form("h_s2s0_beta_FB%d",i)  ,200,  -1,1.5);
    h_s2s0_tof_FB[i]   = new TH1F(Form("h_s2s0_tof_FB%d",i)  , Form("h_s2s0_tof_FB%d",i)   ,800,-100,100);
    h_s2s0_tdiff_FB[i] = new TH1F(Form("h_s2s0_tdiff_FB%d",i), Form("h_s2s0_tdiff_FB%d",i) ,50, -10, 10);
    set->SetTH1(h_s2s0_beta[i]      ,Form("#beta S2%s%d - S0(F1)"              ,LorR.c_str(),i),"#beta"    ,"counts");
    set->SetTH1(h_s2s0_tof[i]       ,Form("ToF S2%s%d - S0"                    ,LorR.c_str(),i),"ToF[ns]"  ,"counts");
    set->SetTH1(h_s2s0_tdiff[i]     ,Form("TDiff (S2%s%d - S0) - ToF calc"     ,LorR.c_str(),i),"Tdiff[ns]","counts");
    set->SetTH1(h_s2s0_beta_FB[i]   ,Form("#beta S2%s%d - S0 (FB)"             ,LorR.c_str(),i),"#beta"    ,"counts");
    set->SetTH1(h_s2s0_tof_FB[i]    ,Form("ToF S2%s%d - S0 (FB)"               ,LorR.c_str(),i),"ToF[ns]"  ,"counts");
    set->SetTH1(h_s2s0_tdiff_FB[i]  ,Form("TDiff (S2%s%d - S0) - ToF calc (FB)",LorR.c_str(),i),"Tdiff[ns]","counts");
    h2_s2s0_f1tof_fbtof[i]    = new TH2F(Form("h2_s2s0_f1tof_fbtof%d",i)  ,Form("h2_s2s0_f1tof_fbtof%d",i), 100, -20, 20, 100, -20, 20);
    h2_s2s0_f1beta_fbbeta[i]  = new TH2F(Form("h2_s2s0_f1beta_fbbeta%d",i),Form("h2_s2s0_f1beta_fbbeta%d",i), 200, -1, 1.5, 200, -1, 1.5);
    set->SetTH2(h2_s2s0_f1tof_fbtof[i],Form("ToF F1 vs Fbus S2%s%d",LorR.c_str(),i),"ToF(F1)[ns]","ToF(Fbus)[ns]");
    set->SetTH2(h2_s2s0_f1beta_fbbeta[i],Form("#beta F1 vs Fbus S2%s%d",LorR.c_str(),i),"#beta(F1)","#beta(Fbus)");

    hcoin_Rs2_t[i]=new TH1F(Form("hcoin_Rs2_t_%d",i),"",200,-5,5);
    set->SetTH1(hcoin_Rs2_t[i],Form("Coincidence time with RS2 Seg %d",i),"coin time [ns]","Counts");
    hcoin_Ls2_t[i]=new TH1F(Form("hcoin_Ls2_t_%d",i),"",200,-5,5);
    set->SetTH1(hcoin_Ls2_t[i],Form("Coincidence time with LS2 Seg %d",i),"coin time [ns]","Counts");
    
    hcoin_Rs2_ra[i]=new TH2F(Form("hcoin_Rs2_ra_%d",i),"",200,0.0,1000,200,-2,2);
    set->SetTH2(hcoin_Rs2_ra[i],Form("Coincidence time vs RS2 ADC Seg %d R-PMT",i),"ADC [ch]","coin time [ns]");

    hcoin_Rs2_la[i]=new TH2F(Form("hcoin_Rs2_la_%d",i),"",200,0.0,1000,200,-2,2);
    set->SetTH2(hcoin_Rs2_la[i],Form("Coincidence time vs RS2 ADC Seg %d L-PMT",i),"ADC [ch]","coin time [ns]");

    hcoin_Ls2_ra[i]=new TH2F(Form("hcoin_Ls2_ra_%d",i),"",200,0.0,1000,200,-2,2);
    set->SetTH2(hcoin_Ls2_ra[i],Form("Coincidence time vs LS2 ADC Seg %d R-PMT",i),"ADC [ch]","coin time [ns]");

    hcoin_Ls2_la[i]=new TH2F(Form("hcoin_Ls2_la_%d",i),"",200,0.0,1000,200,-2,2);
    set->SetTH2(hcoin_Ls2_la[i],Form("Coincidence time vs LS2 ADC Seg %d L-PMT",i),"ADC [ch]","coin time [ns]");


    
  }

  //==== coin hist =====//
  int bin_pi=200;
  double pi_pos=3.0;
  double min_pi=pi_pos-1.0;
  double max_pi=pi_pos+1.0;
  
  hcoin=new TH1F("hcoin","Coin time ; coin time [ns];Counts",1000,-20,20);
  hcoin_Rp=new TH2F("hcoin_Rp","",400,1.6,2.0,bin_pi,min_pi,max_pi);
  set->SetTH2(hcoin_Rp,"RHRS momemtum  vs coin time","Rp [GeV]","coin time [ns]");
  hcoin_Lp=new TH2F("hcoin_Lp","",400,1.9,2.3,bin_pi,min_pi,max_pi);
  set->SetTH2(hcoin_Lp,"LHRS momemtum  vs coin time","Lp [GeV]","coin time [ns]");  
  
}
////////////////////////////////////////////////////////////////////////////
void coin_t0_calib::loop(){

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    if(n%10000==0)cout<<n <<" / "<<ENum<<endl;

    coin_OK=false;
    tree->GetEntry(n);

    int rs2_seg=(int)R_s2_t_pads[0];
    int ls2_seg=(int)L_s2_t_pads[0];
    double coint=coin_time(rs2_seg,ls2_seg);
    
    if(coin_OK){
      hcoin->Fill(coint);
      hcoin_Rs2_ra[rs2_seg]->Fill(R_s2_ra_p[rs2_seg],coint);
      hcoin_Rs2_la[rs2_seg]->Fill(R_s2_la_p[rs2_seg],coint);
      hcoin_Ls2_ra[rs2_seg]->Fill(L_s2_ra_p[ls2_seg],coint);
      hcoin_Ls2_la[rs2_seg]->Fill(L_s2_la_p[ls2_seg],coint);
      hcoin_Rs2_t[rs2_seg]->Fill(coint);
      hcoin_Ls2_t[ls2_seg]->Fill(coint);
      
    }
    if(LR==0)anaL_oneevent();
    else if(LR==1)anaR_oneevent();

    for(int i=0;i<16;i++){
      for(int j=0;j<tr_n;j++){
        bool F1Hits = false, FbusHits = false;
	bool pid_flag=false;
	double Beta=0.0;

	if(LR==0)Beta=L_tr_p[j]/sqrt(L_tr_p[j]*L_tr_p[j]+Me*Me);
	else if(LR==1)Beta=R_tr_p[j]/sqrt(R_tr_p[j]*R_tr_p[j]+Mpi*Mpi);

	
	//====== PID =======//
	if(LR==0)pid_flag=true;
	else if(LR==1 && R_a1_asum_p>50. && R_a2_asum_p>2000)pid_flag=true; // pion selection


        if(s2_trpad[j]==i){

	  //anaR	  cout<<"S2T "<<S2T_F1TDC[i]<<" S2B "<<S2B_F1TDC[i]<<" S0T "<<S0T_F1TDC[i]<<" S0B "<<S0B_F1TDC[i]<<endl;

          if(S2T_F1TDC[i]>0 && S2B_F1TDC[i]>0 && S0T_F1TDC[0]>0. && S0B_F1TDC[0]>0.)F1Hits=true;
          if(S2_lt[i]>0. && S2_rt[i]>0. && S0_lt[0]>0 && S0_rt[0]>0)FbusHits=true;
          if(F1Hits && pid_flag){
            h_s2s0_tdiff[i] ->Fill(S2_F1time[i] - S0_F1time[0] - paths2s0[j]/(Beta*LightVelocity));
            h_s2s0_tof[i]   ->Fill(S2_F1time[i] - S0_F1time[0]);
            h_s2s0_beta[i]  ->Fill(betaF1[j]);
          }
          if(FbusHits && pid_flag){
            h_s2s0_tdiff_FB[i] ->Fill(S2_time[i] - S0_time[0] + paths2s0[j]/(Beta*LightVelocity));
            h_s2s0_tof_FB[i]   ->Fill(S2_time[i] - S0_time[0]);
            h_s2s0_beta_FB[i]  ->Fill(beta[j]);
          }
          if(F1Hits&&FbusHits){
            h2_s2s0_f1tof_fbtof[i]   ->Fill(S2_F1time[i] - S0_F1time[0],S2_time[i] - S0_time[0]);
            h2_s2s0_f1beta_fbbeta[i] ->Fill(betaF1[j],beta[j]);
          }
        }
      }
    }
  }

}
////////////////////////////////////////////////////////////////////////////
void coin_t0_calib::fit(){

  ////////
  //Fbus//
  ////////
  h_frame[0] = new TH2F("h_frame0","h_frame0",10, -0.8, 15.8,10,-0.5,0.5);
  h_frame[1] = new TH2F("h_frame1","h_frame1",10, -0.8, 15.8,10, 0.0,1.0);
  set->SetTH2(h_frame[0] , "TDiff peak pos each S2(Fbus)","S2 paddle","TDiff peak[ns]" );
  set->SetTH2(h_frame[1] , "TDiff width each S2(Fbus)"   ,"S2 paddle","TDiff width[ns]");
  tg_tdiffFB_pos = new TGraphErrors(); 
  tg_tdiffFB_wid = new TGraphErrors();
  set->SetGrErr(tg_tdiffFB_pos, "TDiff peak pos each S2","S2 paddle","TDiff peak[ns]" ,1,4,23);
  set->SetGrErr(tg_tdiffFB_wid, "TDiff width each S2"   ,"S2 paddle","TDiff width[ns]",1,4,24);


  for(int i=0;i<16;i++){
    ga_tdiffFB[i] = new TF1(Form("ga_tdiffFB%d",i+1),"gaus",-2,2);
    set->SetTF1(ga_tdiffFB[i],2,1,1);
    double min=-50,max=50;
    min = h_s2s0_tdiff_FB[i]->GetXaxis()->GetBinCenter(h_s2s0_tdiff_FB[i]->GetMaximumBin()) -3.;
    max = h_s2s0_tdiff_FB[i]->GetXaxis()->GetBinCenter(h_s2s0_tdiff_FB[i]->GetMaximumBin()) +3.;
    set->FitGaus(h_s2s0_tdiff_FB[i],min,max,1.5,5);
    h_s2s0_tdiff_FB[i]->Fit(ga_tdiffFB[i],"QR","",min,max);
    tdiffFB_pos[i]  = ga_tdiffFB[i]->GetParameter(1);
    tdiffFB_wid[i]  = ga_tdiffFB[i]->GetParameter(2);
    etdiffFB_pos[i] = ga_tdiffFB[i]->GetParError(1);
    etdiffFB_wid[i] = ga_tdiffFB[i]->GetParError(2);;

    tg_tdiffFB_pos ->SetPoint(i,i,tdiffFB_pos[i]); 
    tg_tdiffFB_wid ->SetPoint(i,i,tdiffFB_wid[i]);
    tg_tdiffFB_pos ->SetPointError(i,0,etdiffFB_pos[i]); 
    tg_tdiffFB_wid ->SetPointError(i,0,etdiffFB_wid[i]);
    
    param->SetTimeTune(CID_FbS2,i,LR,0,tdiffFB_pos[i]);
    param->SetTimeTune(CID_FbS2,i,LR,1,tdiffFB_pos[i]);
  }

  //////
  //F1//
  //////
  tg_tdiffF1_pos = new TGraphErrors(); 
  tg_tdiffF1_wid = new TGraphErrors();
  set->SetGrErr(tg_tdiffF1_pos, "TDiff peak pos each S2(F1)","S2 paddle","TDiff peak[ns]" ,1,4,23);
  set->SetGrErr(tg_tdiffF1_wid, "TDiff width each S2(F1)"   ,"S2 paddle","TDiff width[ns]",1,4,24);
  h_frame[2] = new TH2F("h_frame0","h_frame0",10, -0.8, 15.8,10,-1,1);
  h_frame[3] = new TH2F("h_frame1","h_frame1",10, -0.8, 15.8,10, 0.0,2.5);
  set->SetTH2(h_frame[2] , "TDiff peak pos each S2(F1TDC)","S2 paddle","TDiff peak[ns]" );
  set->SetTH2(h_frame[3] , "TDiff width each S2(F1TDC)"   ,"S2 paddle","TDiff width[ns]");
  tg_coin_R_pos=new TGraphErrors();
  tg_coin_L_pos=new TGraphErrors();
  set->SetGrErr(tg_coin_R_pos, "S2-Seg(i) vs pion posi (Red=R : Blue=L)","S2 paddle","mean of pion [ns]" ,1,4,23);
  set->SetGrErr(tg_coin_L_pos, "S2-Seg(i) vs pion posi (Red=R : Blue=L)","S2 paddle","mean of pion[ns]",1,4,24);
  tg_coin_R_wid=new TGraphErrors();
  tg_coin_L_wid=new TGraphErrors();
  set->SetGrErr(tg_coin_R_wid, "S2-Seg(i) vs pion width (Red=R : Blue=L)","S2 paddle","width of pion [ns]" ,1,4,23);
  set->SetGrErr(tg_coin_L_wid, "S2-Seg(i) vs pion width (Red=R : Blue=L) ","S2 paddle","width of pion[ns]",1,4,24);  
  
  for(int i=0;i<16;i++){
    ga_tdiffF1[i] = new TF1(Form("ga_tdiffF1%d",i+1),"gaus",-2,2);
    set->SetTF1(ga_tdiffF1[i],2,1,1);
    double min=-50,max=50;
    min = h_s2s0_tdiff[i]->GetXaxis()->GetBinCenter(h_s2s0_tdiff[i]->GetMaximumBin()) -3.;
    max = h_s2s0_tdiff[i]->GetXaxis()->GetBinCenter(h_s2s0_tdiff[i]->GetMaximumBin()) +3.;
    set->FitGaus(h_s2s0_tdiff[i],min,max,1.5,5);
    h_s2s0_tdiff[i]->Fit(ga_tdiffF1[i],"QR","",min,max);
    tdiffF1_pos[i]  = ga_tdiffF1[i]->GetParameter(1);
    tdiffF1_wid[i]  = ga_tdiffF1[i]->GetParameter(2);
    etdiffF1_pos[i] = ga_tdiffF1[i]->GetParError(1);
    etdiffF1_wid[i] = ga_tdiffF1[i]->GetParError(2);;

    tg_tdiffF1_pos ->SetPoint(i,i,tdiffF1_pos[i]); 
    tg_tdiffF1_wid ->SetPoint(i,i,tdiffF1_wid[i]);
    tg_tdiffF1_pos ->SetPointError(i,0,etdiffF1_pos[i]); 
    tg_tdiffF1_wid ->SetPointError(i,0,etdiffF1_wid[i]);

    //==== coin fit =====//
    fcoin_R[i]=new TF1(Form("fcoin_R_%d",i),"gausn(0)",-1,1);
    fcoin_L[i]=new TF1(Form("fcoin_L_%d",i),"gausn(0)",-1,1);
    fcoin_R[i]->SetLineColor(2);
    fcoin_L[i]->SetLineColor(2);
    hcoin_Rs2_t[i]->Fit(fcoin_R[i],"RQ","",-0.75,0.75);
    hcoin_Ls2_t[i]->Fit(fcoin_L[i],"RQ","",-0.75,0.75);
    coin_R[i]=fcoin_R[i]->GetParameter(1);
    coin_L[i]=fcoin_L[i]->GetParameter(1);
    if(fabs(fcoin_R[i]->GetParameter(1))<1) tg_coin_R_pos->SetPoint(i,i,fcoin_R[i]->GetParameter(1));
    if(fabs(fcoin_L[i]->GetParameter(1))<1) tg_coin_L_pos->SetPoint(i,i,fcoin_L[i]->GetParameter(1));
    if(fabs(fcoin_R[i]->GetParameter(1))<1) tg_coin_R_pos->SetPointError(i,0,fcoin_R[i]->GetParError(1));
    if(fabs(fcoin_L[i]->GetParameter(1))<1) tg_coin_L_pos->SetPointError(i,0,fcoin_L[i]->GetParError(1));
    coin_Rw[i]=fcoin_R[i]->GetParameter(2);
    coin_Lw[i]=fcoin_L[i]->GetParameter(2);
    if(fabs(coin_R[i])<1) tg_coin_R_wid->SetPoint(i,i,coin_Rw[i]);
    if(fabs(coin_L[i])<1) tg_coin_L_wid->SetPoint(i,i,coin_Lw[i]);
    if(fabs(coin_R[i])<1) tg_coin_R_wid->SetPointError(i,0,fcoin_R[i]->GetParError(2));
    if(fabs(coin_L[i])<1) tg_coin_L_wid->SetPointError(i,0,fcoin_L[i]->GetParError(2));
    
    
    param->SetTimeTune(CID_F1S2,i,LR,0,tdiffF1_pos[i]);
    param->SetTimeTune(CID_F1S2,i,LR,1,tdiffF1_pos[i]);

  }

    fcoin=new TF1("fcoin","gausn(0)",-5,5);
    hcoin->Fit(fcoin,"RQ","",-1,1);
    coin_pos=fcoin->GetParameter(1);
    coin_wid=fcoin->GetParameter(2);

    
}
////////////////////////////////////////////////////////////////////////////
void coin_t0_calib::draw(){

  c[0]->Clear();c[0]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[0]->cd(i+1);gPad->SetLogy(1);h_s2s0_tof[i]->Draw();
  }

  c[1]->Clear();c[1]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[1]->cd(i+1);gPad->SetLogy(1);h_s2s0_tdiff[i]->Draw();
  }

  c[2]->Clear();c[2]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[2]->cd(i+1);gPad->SetLogy(1);h_s2s0_tof_FB[i]->Draw();
  }

  c[3]->Clear();c[3]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[3]->cd(i+1);gPad->SetLogy(1);h_s2s0_tdiff_FB[i]->Draw();
  }

  c[4]->Clear();c[4]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[4]->cd(i+1);gPad->SetLogz(1);h2_s2s0_f1tof_fbtof[i] ->Draw("colz");
  }
  c[5]->Clear();c[5]->Divide(2,2);
  c[5]->cd(1);h_frame[0]->Draw("");tg_tdiffFB_pos ->Draw("sameP");///
  c[5]->cd(2);h_frame[1]->Draw("");tg_tdiffFB_wid ->Draw("sameP");///
  c[5]->cd(3);h_frame[2]->Draw("");tg_tdiffF1_pos ->Draw("sameP");///
  c[5]->cd(4);h_frame[3]->Draw("");tg_tdiffF1_wid ->Draw("sameP");///

  c[6]->Clear();c[6]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[6]->cd(i+1);gPad->SetLogy(1);h_s2s0_beta[i]->Draw();
  }

  c[7]->Clear();c[7]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[7]->cd(i+1);gPad->SetLogy(1);h_s2s0_beta_FB[i]->Draw();
  }

  c[8]->Clear();c[8]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[8]->cd(i+1);gPad->SetLogz(1);h2_s2s0_f1beta_fbbeta[i] ->Draw("colz");
  }

  c[9]->Clear();c[9]->cd();
  hcoin->Draw();

  c[10]->Clear();c[10]->Divide(4,4);
  c[11]->Clear();c[11]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[10]->cd(i+1);hcoin_Rs2_t[i]->Draw();
    c[11]->cd(i+1);hcoin_Ls2_t[i]->Draw();
  }




  c[12]->Clear();c[12]->Divide(4,4);
  c[13]->Clear();c[13]->Divide(4,4);
  c[14]->Clear();c[14]->Divide(4,4);
  c[15]->Clear();c[15]->Divide(4,4);

  for(int i=0;i<16;i++){
    c[12]->cd(i+1);hcoin_Rs2_ra[i]->Draw("colz");
    c[13]->cd(i+1);hcoin_Rs2_la[i]->Draw("colz");
    c[14]->cd(i+1);hcoin_Ls2_ra[i]->Draw("colz");
    c[15]->cd(i+1);hcoin_Ls2_la[i]->Draw("colz");
  }

  c[16]->Clear();c[16]->cd();  
  tg_coin_R_pos->SetMarkerSize(2);
  tg_coin_R_pos->SetMarkerColor(2);
  tg_coin_R_pos->SetMarkerStyle(20);
  tg_coin_L_pos->SetMarkerSize(2);
  tg_coin_L_pos->SetMarkerColor(4);
  tg_coin_L_pos->SetMarkerStyle(20);  
  tg_coin_R_pos->Draw("AP");
  tg_coin_L_pos->Draw("P");
  TLine* lcoin=new TLine(0,coin_pos,16,coin_pos);
  lcoin->SetLineWidth(2);
  lcoin->Draw("same");

  c[17]->Clear();c[17]->cd();  
  tg_coin_R_wid->SetMarkerSize(2);
  tg_coin_R_wid->SetMarkerColor(2);
  tg_coin_R_wid->SetMarkerStyle(20);
  tg_coin_L_wid->SetMarkerSize(2);
  tg_coin_L_wid->SetMarkerColor(4);
  tg_coin_L_wid->SetMarkerStyle(20);  
  tg_coin_R_wid->Draw("AP");
  tg_coin_L_wid->Draw("P");
  TLine* lcoin_w=new TLine(0,coin_wid,16,coin_wid);
  lcoin_w->SetLineWidth(2);
  lcoin_w->Draw("same");
  
  
  param->WriteToFile("param/tmp.param");
  
}
////////////////////////////////////////////////////////////////////////////
void coin_t0_calib::savecanvas(string ofname){
  c[0]->Print(Form("%s[",ofname.c_str()) );
  for(int i=0;i<NCanvas;i++){
    c[i]->Print(Form("%s" ,ofname.c_str()) );
  }
  c[NCanvas-1]->Print(Form("%s]",ofname.c_str()) );
  cout<<ofname<<" saved"<<endl;
}
////////////////////////////////////////////////////////////////////////////
bool coin_t0_calib::anaL_oneevent(){
  convertF1TDCL(param);

  tr_n = (int)L_tr_n;
  //cout<<"tr_n" <<tr_n<<endl;
  if(tr_n>MAX)tr_n=MAX;

  for(int i=0;i<tr_n;i++){
  //cout<<"s2_trpad : "<<R_s2_trpad[i]<<endl;
    s2_trpad[i]=L_s2_trpad[i];
    paths2s0[i]=L_s2_trpath[i] - L_s0_trpath[i];
    beta[i]    =L_tr_beta[i];
    //betaF1[i]  =GetBeta_S0S2wF1TDCL(i);
  }

  for(int i=0;i<LS2;i++){
    //F1TDC
    S2T_F1TDC[i] = LS2T_F1TDC[i];
    S2B_F1TDC[i] = LS2B_F1TDC[i];
    S2_F1time[i] = LS2_F1time[i];
    //Fbus TDC
    S2_time[i]   = 1.e+9*L_s2_time[i];
    S2_lt[i]     = L_s2_lt[i];
    S2_rt[i]     = L_s2_rt[i];
  }
  for(int i=0;i<LS0;i++){
    //F1TDC
    S0T_F1TDC[i] = LS0T_F1TDC[i];
    S0B_F1TDC[i] = LS0B_F1TDC[i];
    S0_F1time[i] = LS0_F1time[i];
    //Fbus TDC
    S0_time[i]   = 1.e+9*L_s0_time[i];
    S0_lt[i]     = L_s0_lt[i];
    S0_rt[i]     = L_s0_rt[i];
  }


  if(DR_T5>0. ) return true;
  else return false;
}
////////////////////////////////////////////////////////////////////////////
bool coin_t0_calib::anaR_oneevent(){
  //cout<<"coin_t0_calib::anaR_oneevent"<<endl;

  convertF1TDCR(param);
  tr_n = (int)R_tr_n;
  
  if(tr_n>MAX)tr_n=MAX;

  for(int i=0;i<tr_n;i++){
  //cout<<"s2_trpad : "<<R_s2_trpad[i]<<endl;
    s2_trpad[i]=R_s2_trpad[i];
    paths2s0[i]=R_s2_trpath[i] - R_s0_trpath[i];
    beta[i]    =R_tr_beta[i];
    betaF1[i]  =GetBeta_S0S2wF1TDCR(i);
  }

  for(int i=0;i<RS2;i++){
    //F1TDC
    S2T_F1TDC[i] = RS2T_F1TDC[i];
    S2B_F1TDC[i] = RS2B_F1TDC[i];
    S2_F1time[i] = RS2_F1time[i];
    //Fbus TDC
    S2_time[i]   = 1.e+9*R_s2_time[i];
    S2_lt[i]     = R_s2_lt[i];
    S2_rt[i]     = R_s2_rt[i];
  }
  for(int i=0;i<RS0;i++){
    //F1TDC
    S0T_F1TDC[i] = RS0T_F1TDC[i];
    S0B_F1TDC[i] = RS0B_F1TDC[i];
    S0_F1time[i] = RS0_F1time[i];
    //Fbus TDC
    S0_time[i]   = 1.e+9*R_s0_time[i];
    S0_lt[i]     = R_s0_lt[i];
    S0_rt[i]     = R_s0_rt[i];
  }


  if(DR_T5>0. ) return true;
  else return false;
}

///////////////////////////////////////////////////////////////////////////

double coin_t0_calib:: coin_time(int Rs2_seg, int Ls2_seg){

  convertF1TDCR(param);
  convertF1TDCL(param);
  double coin_offset=coin_F1_offset;
  double Rpathl=R_tr_pathl[0]+R_s2_trpath[0];
  double Lpathl=L_tr_pathl[0]+L_s2_trpath[0];
  double Beta_L=L_tr_p[0]/sqrt(L_tr_p[0]*L_tr_p[0]+Me*Me);
  double Beta_R=R_tr_p[0]/sqrt(R_tr_p[0]*R_tr_p[0]+Mpi*Mpi);
  double tof_r=RS2_F1time[Rs2_seg] - Rpathl/(Beta_R*LightVelocity);
  double tof_l=LS2_F1time[Ls2_seg] - Lpathl/(Beta_L*LightVelocity);
  
  if(RS2_F1time[Rs2_seg]!=-9999. &&LS2_F1time[Ls2_seg]!=-9999.)coin_OK=true;
  else coin_OK=false;
  return -tof_r + tof_l - coin_offset;

    
  
}


///////////////////////////////////////////////////////////////////////////

void coin_t0_calib:: OutParameter(string ifname,string ofname){

    ifstream ifp(ifname.c_str());
    ofstream ofp(ofname.c_str(),std::ios::out);
    string title;
    int l=0;
    string st;
    string p;
    string F1_s2L,F1_s2R;
    string t0_st[16];
    double t0_s2R[16],t0_s2L[16];
    string buf;
    double coin_t0[16];


    //======== Print Offset Value =========//
    
    if(LR==0){F1_s2R="L.s2.R.off_F1";F1_s2L="L.s2.L.off_F1";
      for(int i=0;i<16;i++)coin_t0[16]=coin_R[i];
    }else if(LR==1){F1_s2R="R.s2.R.off_F1";F1_s2L="R.s2.L.off_F1";
      for(int i=0;i<16;i++)coin_t0[16]=coin_L[i];}

    string F1_off="Coin.off_F1";
    

    while(!ifp.eof()){
      ifp >> title;

      if(title==F1_s2R || title==F1_s2L || title==F1_off){	
	ifp >> p;// >> p; >> p;
	ofp<<title<<" = ";

	//============ Coin offset =======================//

	if(title==F1_off){ // Set Coin offset 
	  string coint0;
	  ifp  >> coint0;
	  cout<<title<<" "<<p<<" "<<coint0<<endl;
	  double coin_t0=atoi(coint0.c_str());
	  ofp <<coin_t0*F1reso+coin_pos<<" ";

	}else{
	  
	//============ F1 S2 offset =====================//
	  
	  for(int k=0;k<16;k++){
	    ifp>> t0_st[k];
	    if(title==F1_s2R){ // Set S2-RPMT Offset
	      t0_s2R[k]=atoi(t0_st[k].c_str());
	      ofp<<coin_t0[k]/F1reso+t0_s2R[k]<<" ";
	    }else if(title==F1_s2L){// Set S2-LPMT Offset
	      t0_s2L[k]=atoi(t0_st[k].c_str());
	      ofp<<coin_t0[k]/F1reso+t0_s2L[k]<<" ";
	    }
	  }//for
	}
	
	ofp<<endl;
	
      }else{
      getline(ifp,buf);
      ofp <<title<<buf<<endl;

      }
      
    }//while


}



////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "Rootfiles/tritium_94003.root";
  string ofname = "/pdf/test.pdf";
  string paramname = "param/94003.param";
  string ofparam;
  int ch;
  int lr=0;
  int MaxNum = 0;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  bool new_param_flag=false;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:P:bcop:LR"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 'w':
      output_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 'c':
      coin_flag = true;
      break;
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
    case 'p':
      draw_flag = false;
      paramname = optarg;
      cout<<"input param name : "<<paramname<<endl;
      break;
    case 'P':
      new_param_flag=true;
      ofparam = optarg;
      cout<<"output param name : "<<ofparam<<endl;
      break;
    case 'L':
      lr = 0;
      break;
    case 'R':
      lr = 1;
      break;
    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-p : input param file"<<endl;
      cout<<"-L : Left HRS"<<endl;
      cout<<"-R : Right HRS"<<endl;
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

  TApplication *theApp = new TApplication("App", &argc, argv);
  coin_t0_calib *calib = new coin_t0_calib();

  calib->SetMaxEvent(MaxNum);
  calib->SetInputParam(paramname);
  calib->SetLR(lr);
  calib->SetRoot(ifname);
  calib->makehist();
  calib->loop();
  calib->fit();
  if(new_param_flag)calib->OutParameter(paramname,ofparam);

  calib->draw();
  if(output_flag)calib->savecanvas(ofname);
  delete calib;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

