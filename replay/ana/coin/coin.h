#ifndef coin_h
#define coin_h 1
#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"
#include "define.h"
#include "tree.h"


class coin{



 public:
  coin();
  ~coin();

 public:
  double CoinT();
  void DefParam();
 Setting* set=new Setting();

 
  
}


coin(){
  set->Initialize();
  DefParam();}
~coin(){}

//====================================================//
void coin::DefParam(){


}

//===================================================//
void coin::CoinT{

     //--- Set Coincidence time ---------//
     int Ls2pads=(int)Ls2tpads[0];
     int Rs2pads=(int)Rs2tpads[0];
     double Rs2_off=RS2_off_H1[Rs2pads];
     double Ls2_off=LS2_off_H1[Ls2pads];
     double rpathl=rtrpathl[0]+rs2pathl[0];
     double lpathl=ltrpathl[0]+ls2pathl[0];
     double rpath_corr = rpathl/c/(Rp[0]/sqrt(Rp[0]*Rp[0] + MK*MK));
     double lpath_corr=lpathl/c/(Lp[0]/sqrt(Lp[0]*Lp[0] + Me*Me));
     double tof_r=(((-RF1[48+Rs2pads]+RF1[46]-RF1[16+Rs2pads]+RF1[9]+Rs2_off)/2.0))*tdc_time;
     double tof_l=((-LF1[Ls2pads]+LF1[30]-LF1[Ls2pads+48]+LF1[37]+Ls2_off)/2.0)*tdc_time;
     coin_t=tof_r-tof_l+rpath_corr-lpath_corr-pathl_off-coin_offset; //  coin Path & Offset  correction   
     //----------------------------------//     



}

