#ifndef tree_h
#define tree_h 1

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;
#include "TApplication.h"
#include "Param.h"



class tree{

 public:
  tree();
  ~tree();

 public:
  void SetTree(string ifname);
  void Chain_Tree(string ifname);
  void SetBranch();

  
  
   ////////////////
   //// SetRoot ///
  ////////////////
  
  TChain* t1;
  // double xp[nmax], yp[nmax];
  double z_recon[nmax];
  int ENum;
  div_t d;
  Double_t trig5;
  Double_t trig4;
  Double_t trig1; 

  double rtime_s0[MAX], ltime_s0[MAX];
  double rtime_s2[MAX], ltime_s2[MAX];
  double rtime[MAX], ltime[MAX];
  double rpathl[MAX], lpathl[MAX];
  double rpathl_s2[MAX], lpathl_s2[MAX];
  double Rvz[MAX],Lvz[MAX];
  double a1, a2;
  double mom1[MAX], mom2[MAX];
  double rf1tdc[100];
  double lf1tdc[100];
  
  double rvz[MAX], lvz[MAX];
  int Nrvz,Nlvz;
  int NRvz,NLvz;
  double Rth[MAX], Rph[MAX];
  double Lth[MAX], Lph[MAX];
  double runnum;
  double hallap;
  double r_ras_x[MAX],r_ras_y[MAX];
  double l_ras_x[MAX],l_ras_y[MAX];
  double r_s2_t_pads[MAX];
  double l_s2_t_pads[MAX];
  double r_s2_nthit;
  double l_s2_nthit;
  double r_th_fp[MAX];
  double l_th_fp[MAX];
  double r_ph_fp[MAX];
  double l_ph_fp[MAX];
  double l_x_fp[MAX];
  double r_x_fp[MAX];
  double l_y_fp[MAX];
  double r_y_fp[MAX];
  double lx[MAX];
  double ly[MAX];
  double rx[MAX];
  double ry[MAX];
  double r_s2_la_c[16];
  double r_s2_ra_c[16];
  double l_s2_la_c[16];
  double l_s2_ra_c[16];
  double rbeta[MAX];
  double lbeta[MAX];
  double nhit, nhit_R;
  double R_ps_asum,R_sh_asum;
  double L_gs_asum;
  double R_a1_asum,R_a2_asum;
  double a1_tdc[24];
  double a2_tdc[26];
  double ctime[MAX];
  double ztR[MAX];
  double ztR_opt[MAX];
  double ztR_tuned[MAX];
  double R_Ras_x,R_Ras_y,L_Ras_x,L_Ras_y;

  

};



/////////////////////////////////////////////////////////////////////////
tree::tree(){
  t1=new TChain("T");
  
}
tree::~tree(){}
/////////////////////////////////////////////////////////////////////////
void tree::SetTree(string ifname){


  t1->Add(ifname.c_str());
  cout<<ifname.c_str()<<endl;
  ENum=t1->GetEntries();
  
}
/////////////////////////////////////////////////////////////////////////

void tree::Chain_Tree(string ifname)
{

   ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    t1->Add(runname.c_str());
    cout<<runname.c_str()<<endl;
  }

  ENum=t1->GetEntries();
}
////////////////////////////////////////////////////////////////////////

void tree::SetBranch(){

  cout<<"====================================="<<endl;
  cout<<"======== SetBranch  ================="<<endl;
  cout<<"===================================="<<endl;
  
  
  
  t1->SetBranchStatus("*",0);  

  t1->SetBranchStatus("fEvtHdr.fRun",1);
  t1->SetBranchStatus("HALLA_p",     1);
  t1->SetBranchStatus("DR.T1",       1);
  t1->SetBranchStatus("DR.T4",       1);
  t1->SetBranchStatus("DR.T5",       1);
  t1->SetBranchStatus("Ndata.R.tr.vz", 1);
  t1->SetBranchStatus("Ndata.L.tr.vz", 1);
  t1->SetBranchStatus("R.tr.vz",     1);
  t1->SetBranchStatus("L.tr.vz",     1);
  t1->SetBranchStatus("R.tr.x",      1);
  t1->SetBranchStatus("L.tr.x",      1);
  t1->SetBranchStatus("R.tr.y",      1);
  t1->SetBranchStatus("L.tr.y",      1);
  t1->SetBranchStatus("R.tr.vx",     1);
  t1->SetBranchStatus("L.tr.vx",     1);
  t1->SetBranchStatus("R.tr.vy",     1);
  t1->SetBranchStatus("L.tr.vy",     1);
  t1->SetBranchStatus("R.tr.th",     1);
  t1->SetBranchStatus("L.tr.th",     1);
  t1->SetBranchStatus("R.tr.ph",     1);
  t1->SetBranchStatus("L.tr.ph",     1);
  t1->SetBranchStatus("R.tr.tg_th",  1);
  t1->SetBranchStatus("R.tr.tg_ph",  1);
  t1->SetBranchStatus("L.tr.tg_th",  1);
  t1->SetBranchStatus("L.tr.tg_ph",  1);
  t1->SetBranchStatus("R.ps.asum_c", 1);
  t1->SetBranchStatus("R.sh.asum_c", 1);
  t1->SetBranchStatus("L.cer.asum_c",1);
  t1->SetBranchStatus("R.a1.asum_c", 1);
  t1->SetBranchStatus("R.a2.asum_c", 1);
  t1->SetBranchStatus("Lrb.Raster2.rawcur.x", 1); // raster current
  t1->SetBranchStatus("Lrb.Raster2.rawcur.y", 1); // raster current
  t1->SetBranchStatus("Rrb.Raster2.rawcur.x", 1); // raster current
  t1->SetBranchStatus("Rrb.Raster2.rawcur.y", 1); // raster current

  t1->SetBranchAddress("fEvtHdr.fRun", &runnum);
  t1->SetBranchAddress("HALLA_p", &hallap);
  t1->SetBranchAddress("DR.T1", &trig1);
  t1->SetBranchAddress("DR.T4", &trig4);
  t1->SetBranchAddress("DR.T5", &trig5);
  t1->SetBranchAddress("Ndata.R.tr.vz",&Nrvz);
  t1->SetBranchAddress("Ndata.L.tr.vz",&Nlvz);
  t1->SetBranchAddress("R.tr.vz", rvz);
  t1->SetBranchAddress("L.tr.vz", lvz);
  t1->SetBranchAddress("R.tr.x",   &r_x_fp);
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("R.tr.y",   &r_y_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("R.tr.vx",   &rx);
  t1->SetBranchAddress("L.tr.vx",   &lx);
  t1->SetBranchAddress("R.tr.vy",   &ry);
  t1->SetBranchAddress("L.tr.vy",   &ly);
  t1->SetBranchAddress("R.tr.th",  &r_th_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("R.tr.ph",  &r_ph_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);
  t1->SetBranchAddress("R.tr.tg_th", &Rth);
  t1->SetBranchAddress("R.tr.tg_ph", &Rph);
  t1->SetBranchAddress("L.tr.tg_th", &Lth);
  t1->SetBranchAddress("L.tr.tg_ph", &Lph);
  t1->SetBranchAddress("R.ps.asum_c", &R_ps_asum);
  t1->SetBranchAddress("R.sh.asum_c", &R_sh_asum);
  t1->SetBranchAddress("L.cer.asum_c",&L_gs_asum);
  t1->SetBranchAddress("R.a1.asum_c", &R_a1_asum);
  t1->SetBranchAddress("R.a2.asum_c", &R_a2_asum);
  t1->SetBranchAddress("Lrb.Raster2.rawcur.x", &L_Ras_x); // raster current
  t1->SetBranchAddress("Lrb.Raster2.rawcur.y", &L_Ras_y); // raster current
  t1->SetBranchAddress("Rrb.Raster2.rawcur.x", &R_Ras_x); // raster current
  t1->SetBranchAddress("Rrb.Raster2.rawcur.y", &R_Ras_y); // raster current

}

///////////////////////////////////////////////////////////////////////

#endif
