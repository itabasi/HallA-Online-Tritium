#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#include "tree.h"

//============ tree ===================//
tree::tree()
{
t1=new TChain("T");
}
//============ ~tree ===================//
tree::~tree(){}

//=========== chain_tree ===============//
void tree::chain_tree(string ifname)
{

  tree->Add(Form("%s",ifname.c_str()));
  cout<<ifname<<endl;

}
//============ readtree ================//
void tree::readtree(){

  
  t1->SetBranchAddress("fEvtHdr.fRun", &runnum    );
  t1->SetBranchAddress("HALLA_p", &hallap );
  t1->SetBranchAddress("DR.T1", &trig1    );
  t1->SetBranchAddress("DR.T4", &trig4    );
  t1->SetBranchAddress("DR.T5", &trig5    );
  //t1->SetBranchAddress("R.tr.time", &rtime);
  //t1->SetBranchAddress("L.tr.time", &ltime);
  //t1->SetBranchAddress("R.tr.pathl", &rpathl);
  //t1->SetBranchAddress("L.tr.pathl", &lpathl);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);
  t1->SetBranchAddress("R.tr.p", &mom1);
  t1->SetBranchAddress("L.tr.p", &mom2);
  //t1->SetBranchAddress("RTDC.F1FirstHit", &rf1tdc);
  //t1->SetBranchAddress("LTDC.F1FirstHit", &lf1tdc);
  t1->SetBranchAddress("R.tr.vz", &rvz);
  t1->SetBranchAddress("L.tr.vz", &lvz);
  t1->SetBranchAddress("R.tr.tg_th", &th1);
  t1->SetBranchAddress("R.tr.tg_ph", &ph1);
  t1->SetBranchAddress("L.tr.tg_th", &th2);
  t1->SetBranchAddress("L.tr.tg_ph", &ph2);
  t1->SetBranchAddress("R.s0.time", &rtime_s0);
  t1->SetBranchAddress("L.s0.time", &ltime_s0);
  t1->SetBranchAddress("R.s2.time", &rtime_s2);
  t1->SetBranchAddress("L.s2.time", &ltime_s2);
  //t1->SetBranchAddress("R.s2.t_pads", &r_s2_t_pads);
  //t1->SetBranchAddress("L.s2.t_pads", &l_s2_t_pads);
  //t1->SetBranchAddress("R.s2.nthit",   &r_s2_nthit);
  //t1->SetBranchAddress("L.s2.nthit",   &l_s2_nthit);
  t1->SetBranchAddress("R.tr.x",   &r_x_fp);
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("R.tr.y",   &r_y_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("R.tr.th",  &r_th_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("R.tr.ph",  &r_ph_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);
  t1->SetBranchAddress("R.tr.beta",  &rbeta);
  t1->SetBranchAddress("L.tr.beta",  &lbeta);
  //t1->SetBranchAddress("R.s2.trpath",  &rpathl_s2);
  //t1->SetBranchAddress("L.s2.trpath",  &lpathl_s2);
  //t1->SetBranchAddress("L.s2.nthit",&nhit);
  //t1->SetBranchAddress("R.s2.nthit",&nhit_R);
  t1->SetBranchAddress("R.ps.asum_c", &ps_asum);
  //t1->SetBranchAddress("ctime", &ctime);
  t1->SetBranchAddress("R.a1.t_fadc", &a1_tdc);
  t1->SetBranchAddress("R.a2.t_fadc", &a2_tdc);

}


