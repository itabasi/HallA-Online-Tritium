using namespace std;
#include "coincalib.h"
#include "t0calib.h"
#include "Param.h"
#include <TMinuit.h>


t0calib::MakeHist(){

  min_ct=-100.;
  max_ct=100;
  bin_ct=1000;
  
  for(int i=0;i<ns2seg;i++){
    hcoinRs2T[i] = new TH1D(Form("hcoinRs2T_%d",i),bin_ct,min_ct,max_ct);
    hcoinRs2B[i] = new TH1D(Form("hcoinRs2B_%d",i),bin_ct,min_ct,max_ct);
    hcoinLs2T[i] = new TH1D(Form("hcoinLs2T_%d",i),bin_ct,min_ct,max_ct);
    hcoinLs2B[i] = new TH1D(Form("hcoinLs2B_%d",i),bin_ct,min_ct,max_ct);
    hs0Rs2T[i] = new TH1D(Form("hs0Rs2T_%d",i),bin_tof,min_tof,max_tof);
    hs0Rs2B[i] = new TH1D(Form("hs0Rs2B_%d",i),bin_tof,min_tof,max_tof);
    hs0Ls2T[i] = new TH1D(Form("hs0Ls2T_%d",i),bin_tof,min_tof,max_tof);
    hs0Ls2B[i] = new TH1D(Form("hs0Ls2B_%d",i),bin_tof,min_tof,max_tof);
    
    }


  
}


