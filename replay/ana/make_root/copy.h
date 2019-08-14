#ifndef copy_h
#define copy_h 1
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
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
#include "TRandom.h"
#include "Setting.h"

//=====================================================//
//============ Make Rootfiles =========================//
//=====================================================//

class copy{

  Setting* set;

 public:
  copy();
  ~copy();

 public:
  void CoinBranch();
  void ACBranch();
  void VDCBranch();
  void ScintiBranch();
  void SCBranch();
  void FPBranch();
  void TGBranch();

  TChain *oldtree = new TChain("T");
}


copy::copy(){
  set->Initialize();
  oldtree->SetBranchStatus("*",0);   }
copy::~copy(){};

//============================================================//

void copy::ACBranch(){

  oldtree->SetBranchStatus("R.a1.a"          ,1);
  oldtree->SetBranchStatus("R.a1.a_c"        ,1);
  oldtree->SetBranchStatus("R.a1.a_p"        ,1);
  //  oldtree->SetBranchStatus("R.a1.t"          ,1);
  //  oldtree->SetBranchStatus("R.a1.t_c"        ,1);
  oldtree->SetBranchStatus("R.a1.asum_p"     ,1);
  oldtree->SetBranchStatus("R.a1.asum_c"     ,1);
  //  oldtree->SetBranchStatus("R.a1.trpath"     ,1);
  oldtree->SetBranchStatus("R.a1.trx"        ,1);
  oldtree->SetBranchStatus("R.a1.try"        ,1);
  oldtree->SetBranchStatus("R.a1.nhits"      ,1);
  oldtree->SetBranchStatus("R.a1.peak"       ,1);
  //  oldtree->SetBranchStatus("R.a1.t_fadc"     ,1);
  //  oldtree->SetBranchStatus("R.a1.tc_fadc"    ,1);
  oldtree->SetBranchStatus("R.a1.nahit"      ,1);     
  oldtree->SetBranchStatus("R.a2.a"          ,1);
  oldtree->SetBranchStatus("R.a2.a_c"        ,1);
  oldtree->SetBranchStatus("R.a2.a_p"        ,1);
  //  oldtree->SetBranchStatus("R.a2.t"          ,1);
  //  oldtree->SetBranchStatus("R.a2.t_c"        ,1);
  oldtree->SetBranchStatus("R.a2.asum_p"     ,1);
  oldtree->SetBranchStatus("R.a2.asum_c"     ,1);
  oldtree->SetBranchStatus("R.a2.trpath"     ,1);
  oldtree->SetBranchStatus("R.a2.trx"        ,1);
  oldtree->SetBranchStatus("R.a2.try"        ,1);
  oldtree->SetBranchStatus("R.a2.nhits"      ,1);
  oldtree->SetBranchStatus("R.a2.peak"       ,1);
  //  oldtree->SetBranchStatus("R.a2.t_fadc"     ,1);
  //  oldtree->SetBranchStatus("R.a2.tc_fadc"    ,1);
  oldtree->SetBranchStatus("R.a2.nahit"      ,1); 


}

//============================================================//

void copy::VDCBranch(){

  
  //==========================//
  //====== RHRS VDC ==========//
  //==========================//

  oldtree->SetBranchStatus("R.vdc.u1.nhit"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.nhit"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.nhit"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.R.vdc.u1.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.R.vdc.u2.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.R.vdc.v1.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.R.vdc.v2.nhit"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.wire"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.wire"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.wire"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.wire"          ,1);
  oldtree->SetBranchStatus("R.vdc.u1.rawtime"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.rawtime"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.rawtime"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.rawtime"          ,1);
  oldtree->SetBranchStatus("R.vdc.u1.time"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.time"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.time"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.time"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.time"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.time"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.time"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.time"          ,1);
  oldtree->SetBranchStatus("R.vdc.u1.dist"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.dist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.dist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.dist"          ,1);  
  oldtree->SetBranchStatus("R.vdc.u1.ddist"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.ddist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.ddist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.ddist"          ,1);    
  oldtree->SetBranchStatus("R.vdc.u1.trdist"          ,1);
  oldtree->SetBranchStatus("R.vdc.u2.trdist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v1.trdist"          ,1);
  oldtree->SetBranchStatus("R.vdc.v2.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u1.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.u2.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v1.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.R.vdc.v2.trdist"          ,1);  

  //==========================//
  //====== LHRS VDC ==========//
  //==========================//

  oldtree->SetBranchStatus("L.vdc.u1.nhit"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.nhit"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.nhit"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.L.vdc.u1.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.L.vdc.u2.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.L.vdc.v1.nhit"          ,1);
  //  oldtree->SetBranchStatus("Ndata.L.vdc.v2.nhit"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.wire"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.wire"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.wire"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.wire"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.wire"          ,1);
  oldtree->SetBranchStatus("L.vdc.u1.rawtime"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.rawtime"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.rawtime"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.rawtime"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.rawtime"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.time"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.time"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.time"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.time"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.time"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.time"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.time"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.time"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.dist"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.dist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.dist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.dist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.dist"          ,1);  
  oldtree->SetBranchStatus("L.vdc.u1.ddist"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.ddist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.ddist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.ddist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.ddist"          ,1);    
  oldtree->SetBranchStatus("L.vdc.u1.trdist"          ,1);
  oldtree->SetBranchStatus("L.vdc.u2.trdist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v1.trdist"          ,1);
  oldtree->SetBranchStatus("L.vdc.v2.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u1.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.u2.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v1.trdist"          ,1);
  oldtree->SetBranchStatus("Ndata.L.vdc.v2.trdist"          ,1);      
  
}

//============================================================//

void ScintiBranch(){


  
  oldtree->SetBranchStatus("R.s0.la"         ,1);
  oldtree->SetBranchStatus("R.s0.la_c"       ,1);
  oldtree->SetBranchStatus("R.s0.la_p"       ,1);
  oldtree->SetBranchStatus("R.s0.lt"         ,1);
  oldtree->SetBranchStatus("R.s0.lt_c"       ,1);
  oldtree->SetBranchStatus("R.s0.ra"         ,1);
  oldtree->SetBranchStatus("R.s0.ra_c"       ,1);
  oldtree->SetBranchStatus("R.s0.ra_p"       ,1);
  oldtree->SetBranchStatus("R.s0.rt"         ,1);
  oldtree->SetBranchStatus("R.s0.rt_c"       ,1);
  oldtree->SetBranchStatus("R.s0.t_pads"     ,1);
  oldtree->SetBranchStatus("R.s0.time"       ,1);
  oldtree->SetBranchStatus("R.s0.dedx"       ,1);
  oldtree->SetBranchStatus("R.s0.trdy"       ,1);
  oldtree->SetBranchStatus("R.s0.troff"      ,1);
  oldtree->SetBranchStatus("R.s0.trpad"      ,1);
  oldtree->SetBranchStatus("R.s0.trpath"     ,1);
  oldtree->SetBranchStatus("R.s0.trx"        ,1);
  oldtree->SetBranchStatus("R.s0.try"        ,1);
  oldtree->SetBranchStatus("R.s0.lnhits"     ,1);
  oldtree->SetBranchStatus("R.s0.lpeak"      ,1);
  oldtree->SetBranchStatus("R.s0.lt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s0.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("R.s0.rnhits"     ,1);
  oldtree->SetBranchStatus("R.s0.rpeak"      ,1);
  oldtree->SetBranchStatus("R.s0.rt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s0.rtc_fadc"   ,1);
  oldtree->SetBranchStatus("R.s2.la"         ,1);
  oldtree->SetBranchStatus("R.s2.la_c"       ,1);
  oldtree->SetBranchStatus("R.s2.la_p"       ,1);
  oldtree->SetBranchStatus("R.s2.lt"         ,1);
  oldtree->SetBranchStatus("R.s2.lt_c"       ,1);
  oldtree->SetBranchStatus("R.s2.ra"         ,1);
  oldtree->SetBranchStatus("R.s2.ra_c"       ,1);
  oldtree->SetBranchStatus("R.s2.ra_p"       ,1);
  oldtree->SetBranchStatus("R.s2.rt"         ,1);
  oldtree->SetBranchStatus("R.s2.rt_c"       ,1);
  oldtree->SetBranchStatus("R.s2.t_pads"     ,1);
  oldtree->SetBranchStatus("R.s2.time"       ,1);
  oldtree->SetBranchStatus("R.s2.dedx"       ,1);
  oldtree->SetBranchStatus("R.s2.trdx"       ,1);
  oldtree->SetBranchStatus("R.s2.troff"      ,1);
  oldtree->SetBranchStatus("R.s2.trpad"      ,1);
  oldtree->SetBranchStatus("R.s2.trpath"     ,1);
  oldtree->SetBranchStatus("R.s2.trx"        ,1);
  oldtree->SetBranchStatus("R.s2.try"        ,1);
  oldtree->SetBranchStatus("R.s2.lnhits"     ,1);
  oldtree->SetBranchStatus("R.s2.lpeak"      ,1);
  oldtree->SetBranchStatus("R.s2.lt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s2.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("R.s2.rnhits"     ,1);
  oldtree->SetBranchStatus("R.s2.rpeak"      ,1);
  oldtree->SetBranchStatus("R.s2.rt_fadc"    ,1);
  oldtree->SetBranchStatus("R.s2.rtc_fadc"   ,1);


  oldtree->SetBranchStatus("L.s0.la"         ,1);
  oldtree->SetBranchStatus("L.s0.la_c"       ,1);
  oldtree->SetBranchStatus("L.s0.la_p"       ,1);
  oldtree->SetBranchStatus("L.s0.lt"         ,1);
  oldtree->SetBranchStatus("L.s0.lt_c"       ,1);
  oldtree->SetBranchStatus("L.s0.ra"         ,1);
  oldtree->SetBranchStatus("L.s0.ra_c"       ,1);
  oldtree->SetBranchStatus("L.s0.ra_p"       ,1);
  oldtree->SetBranchStatus("L.s0.rt"         ,1);
  oldtree->SetBranchStatus("L.s0.rt_c"       ,1);
  oldtree->SetBranchStatus("L.s0.t_pads"     ,1);
  oldtree->SetBranchStatus("L.s0.time"       ,1);
  oldtree->SetBranchStatus("L.s0.dedx"       ,1);
  oldtree->SetBranchStatus("L.s0.trdy"       ,1);
  oldtree->SetBranchStatus("L.s0.troff"      ,1);
  oldtree->SetBranchStatus("L.s0.trpad"      ,1);
  oldtree->SetBranchStatus("L.s0.trpath"     ,1);
  oldtree->SetBranchStatus("L.s0.trx"        ,1);
  oldtree->SetBranchStatus("L.s0.try"        ,1);
  oldtree->SetBranchStatus("L.s0.lnhits"     ,1);
  oldtree->SetBranchStatus("L.s0.lpeak"      ,1);
  oldtree->SetBranchStatus("L.s0.lt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s0.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("L.s0.rnhits"     ,1);
  oldtree->SetBranchStatus("L.s0.rpeak"      ,1);
  oldtree->SetBranchStatus("L.s0.rt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s0.rtc_fadc"   ,1);
  oldtree->SetBranchStatus("L.s2.la"         ,1);
  oldtree->SetBranchStatus("L.s2.la_c"       ,1);
  oldtree->SetBranchStatus("L.s2.la_p"       ,1);
  oldtree->SetBranchStatus("L.s2.lt"         ,1);
  oldtree->SetBranchStatus("L.s2.lt_c"       ,1);
  oldtree->SetBranchStatus("L.s2.ra"         ,1);
  oldtree->SetBranchStatus("L.s2.ra_c"       ,1);
  oldtree->SetBranchStatus("L.s2.ra_p"       ,1);
  oldtree->SetBranchStatus("L.s2.rt"         ,1);
  oldtree->SetBranchStatus("L.s2.rt_c"       ,1);
  oldtree->SetBranchStatus("L.s2.t_pads"     ,1);
  oldtree->SetBranchStatus("L.s2.time"       ,1);
  oldtree->SetBranchStatus("L.s2.dedx"       ,1);
  oldtree->SetBranchStatus("L.s2.trdx"       ,1);
  oldtree->SetBranchStatus("L.s2.troff"      ,1);
  oldtree->SetBranchStatus("L.s2.trpad"      ,1);
  oldtree->SetBranchStatus("L.s2.trpath"     ,1);
  oldtree->SetBranchStatus("L.s2.trx"        ,1);
  oldtree->SetBranchStatus("L.s2.try"        ,1);
  oldtree->SetBranchStatus("L.s2.lnhits"     ,1);
  oldtree->SetBranchStatus("L.s2.lpeak"      ,1);
  oldtree->SetBranchStatus("L.s2.lt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s2.ltc_fadc"   ,1);
  oldtree->SetBranchStatus("L.s2.rnhits"     ,1);
  oldtree->SetBranchStatus("L.s2.rpeak"      ,1);
  oldtree->SetBranchStatus("L.s2.rt_fadc"    ,1);
  oldtree->SetBranchStatus("L.s2.rtc_fadc"   ,1);

  
}

//============================================================//

void SCBranch(){
  
  oldtree->SetBranchStatus("L.cer.a"         ,1);
  oldtree->SetBranchStatus("L.cer.a_c"       ,1);
  oldtree->SetBranchStatus("L.cer.a_p"       ,1);
  oldtree->SetBranchStatus("L.cer.t"         ,1);
  oldtree->SetBranchStatus("L.cer.t_c"       ,1);
  oldtree->SetBranchStatus("L.cer.asum_p"    ,1);
  oldtree->SetBranchStatus("L.cer.asum_c"    ,1);
  oldtree->SetBranchStatus("L.cer.trpath"    ,1);
  oldtree->SetBranchStatus("L.cer.trx"       ,1);
  oldtree->SetBranchStatus("L.cer.try"       ,1);
  oldtree->SetBranchStatus("L.cer.nhits"     ,1);
  oldtree->SetBranchStatus("L.cer.peak"      ,1);
  oldtree->SetBranchStatus("L.cer.t_fadc"    ,1);
  oldtree->SetBranchStatus("L.cer.tc_fadc"   ,1);
  oldtree->SetBranchStatus("L.ps.asum_p"     ,1);
  oldtree->SetBranchStatus("L.ps.asum_c"     ,1);
  oldtree->SetBranchStatus("L.sh.asum_p"     ,1);
  oldtree->SetBranchStatus("L.sh.asum_c"     ,1);


}

//============================================================//

void copy::FPBranch(){

  oldtree->SetBranchStatus("R.tr.n"          ,1);
  //oldtree->SetBranchAddress("R.tr.n" ,&R_tr_n );
  oldtree->SetBranchStatus("R.tr.flag"       ,1);
  oldtree->SetBranchStatus("R.tr.ndof"       ,1);
  oldtree->SetBranchStatus("R.tr.chi2"       ,1);
  oldtree->SetBranchStatus("R.tr.beta"       ,1);
  oldtree->SetBranchStatus("R.tr.d_x"        ,1);
  oldtree->SetBranchStatus("R.tr.d_y"        ,1);
  oldtree->SetBranchStatus("R.tr.d_th"       ,1);
  oldtree->SetBranchStatus("R.tr.d_ph"       ,1);
  oldtree->SetBranchStatus("R.tr.r_x"        ,1);
  //oldtree->SetBranchAddress("R.tr.x" ,R_tr_x );
  oldtree->SetBranchStatus("R.tr.r_y"        ,1);
  oldtree->SetBranchStatus("R.tr.r_th"       ,1);
  //oldtree->SetBranchAddress("R.tr.th",R_tr_th);
  oldtree->SetBranchStatus("R.tr.r_ph"       ,1);
  oldtree->SetBranchStatus("R.tr.x"          ,1);
  oldtree->SetBranchStatus("R.tr.y"          ,1);
  //oldtree->SetBranchAddress("R.tr.p" ,R_tr_p );
  oldtree->SetBranchStatus("R.tr.th"         ,1);
  oldtree->SetBranchStatus("R.tr.ph"         ,1);
  oldtree->SetBranchStatus("R.tr.time"       ,1);
  oldtree->SetBranchStatus("R.tr.p"          ,1);
  oldtree->SetBranchStatus("R.tr.pathl"      ,1);
  oldtree->SetBranchStatus("R.tr.px"         ,1);
  oldtree->SetBranchStatus("R.tr.py"         ,1);
  oldtree->SetBranchStatus("R.tr.pz"         ,1);  
  oldtree->SetBranchStatus("R.tr.vx"         ,1);
  oldtree->SetBranchStatus("R.tr.vy"         ,1);
  oldtree->SetBranchStatus("R.tr.vz"         ,1);
  oldtree->SetBranchStatus("L.tr.n"          ,1);
  //oldtree->SetBranchAddress("L.tr.n" ,&L_tr_n );
  oldtree->SetBranchStatus("L.tr.flag"       ,1);
  oldtree->SetBranchStatus("L.tr.ndof"       ,1);
  oldtree->SetBranchStatus("L.tr.chi2"       ,1);
  oldtree->SetBranchStatus("L.tr.beta"       ,1);
  oldtree->SetBranchStatus("L.tr.d_x"        ,1);
  oldtree->SetBranchStatus("L.tr.d_y"        ,1);
  oldtree->SetBranchStatus("L.tr.d_th"       ,1);
  oldtree->SetBranchStatus("L.tr.d_ph"       ,1);
  oldtree->SetBranchStatus("L.tr.r_x"        ,1);
  oldtree->SetBranchStatus("L.tr.r_y"        ,1);
  oldtree->SetBranchStatus("L.tr.r_th"       ,1);
  oldtree->SetBranchStatus("L.tr.r_ph"       ,1);
  oldtree->SetBranchStatus("L.tr.x"          ,1);
  //oldtree->SetBranchAddress("L.tr.x" ,L_tr_x );
  oldtree->SetBranchStatus("L.tr.y"          ,1);
  oldtree->SetBranchStatus("L.tr.th"         ,1);
  //oldtree->SetBranchAddress("L.tr.th",L_tr_th);
  oldtree->SetBranchStatus("L.tr.ph"         ,1);
  oldtree->SetBranchStatus("L.tr.time"       ,1);
  oldtree->SetBranchStatus("L.tr.p"          ,1);
  //oldtree->SetBranchAddress("L.tr.p" ,L_tr_p );
  oldtree->SetBranchStatus("L.tr.pathl"      ,1);
  oldtree->SetBranchStatus("L.tr.px"         ,1);
  oldtree->SetBranchStatus("L.tr.py"         ,1);
  oldtree->SetBranchStatus("L.tr.pz"         ,1);
  oldtree->SetBranchStatus("L.tr.vx"         ,1);
  oldtree->SetBranchStatus("L.tr.vy"         ,1);
  oldtree->SetBranchStatus("L.tr.vz"         ,1);
  
}

//============================================================//

#endif
