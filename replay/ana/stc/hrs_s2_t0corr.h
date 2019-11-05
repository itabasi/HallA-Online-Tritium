#ifndef hrs_s2_t0corr_h
#define hrs_s2_t0corr_h 1
#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"
#include "define.h"
#include "Tree.h"


//=================================================//
//=============   s2 t0 Class  ====================//
//=================================================//

class s2t0 : public tree
{

 public:
  s2t0();
  ~s2t0();

 public:
  
  void SetRoot(string ifname);
  void SetRunList(string ifname);
  void GetT0();
  void Fill();
  void MakeHist();
  void NewRoot();
  Setting *set = new Setting();


  //=============================//

  double tdc_time;
  const  int s2seg=16;

  //===== Fill======//
  int  ENum;
  double rs2_rt[s2seg].rs2_lt[s2seg];
  double ls2_rt[s2seg],ls2_lt[s2seg];
  int rs2_seg,ls2_seg;
  
  //=== MakeHist ===//
  double min_s2,max_s2;
  int bin_s2;
  double min_coin,max_coin;
  int bin_coin;

  TH1D* hcoin;
  TH1D* hRs2_rt[s2seg];
  TH1D* hLs2_rt[s2seg];
  TH2D* hcoin_rs2_rt[s2seg];
  TH2D* hcoin_rs2_rt[s2seg];
  TH1D* hRs2_lt[s2seg];
  TH1D* hLs2_lt[s2seg];
  TH2D* hcoin_rs2_lt[s2seg];
  TH2D* hcoin_rs2_lt[s2seg];


  
}







#endif
