#include "TChain.h"


void MakeHist(){
  
  string  rname ="../rootfiles/mmass/ana_Lambda/mmcalib_nmr_small/nnL_small_Ole_all.root";

  
  TChain* T =new TChain("T");
  
  T->Add(rname.c_str());
  
  int ENum = T->GetEntries();
  double mm_nnL,mm_acc[10];
  T->SetBranchAddress("mm_nnL",&mm_nnL);
  T->SetBranchAddress("mm_mix",mm_acc);
  double min_mm = -300;
  double max_mm =  200;
  int bin_mm = (int)(max_mm - min_mm)*100;
  TH1D* hmm = new TH1D("hmm","",bin_mm, min_mm,max_mm);
  TH1D* hmm_peak = new TH1D("hmm_peak","",bin_mm, min_mm,max_mm);
  TH1D* hmm_acc = new TH1D("hmm_acc","",bin_mm, min_mm,max_mm);
    
  
  
  for(int i=0; i<ENum;i++){
    
    T->GetEntry(i);
    hmm->Fill(mm_nnL);
    hmm_peak->Fill(mm_nnL);
    for(int i=0;i<10;i++)
      hmm_acc->Fill(mm_acc[i],1./10./20.);
    
  }
  
  
  string ofrname ="../rootfiles/fsi/mmass/nnL_exp_hist.root";
  //    ofstream* ofr =new ofstream(ofrname.c_str());
  TFile* ofr = new TFile(ofrname.c_str(),"recreate");

  hmm_peak->Add(hmm_acc,-1);
  hmm->Write();
  hmm_acc->Write();
  hmm_peak->Write();
    

}
