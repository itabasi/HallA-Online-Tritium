#include "TChain.h"

//bool nnL=false;
bool nnL=true;

void MakeHist_simc(){
  
  //string  rname ="../rootfiles/mmass/ana_Lambda/mmcalib_nmr_small/nnL_small_Ole_all.root";
  //  string  rname ="../rootfiles/mmass/ana_Lambda/test/nnL_small_Ole_all.root";
  //string  rname ="../rootfiles/simc/nnL_simc.root";
  
  string rname;// ="../rootfiles/fsi/mmass/nnL_Jl0.root";

  if(nnL) rname ="../rootfiles/fsi/mmass/nnL_Jl0.root";
  else    rname ="../rootfiles/fsi/mmass/Lam_Jl0.root";
  TChain* T =new TChain("SNT");
  T->Add(rname.c_str());
  
  int ENum = T->GetEntries();
  float momL,momR,mm_nnL,mm_L;
  double Prel,fac1,fac2,fac3;
  
  T->SetBranchAddress("mm_nnL",&mm_nnL);
  T->SetBranchAddress("mm_L",&mm_L);
  T->SetBranchAddress("Lp_rec",&momL);
  T->SetBranchAddress("Rp_rec",&momR);
  //  T->SetBranchAddress("Prel",&Prel);
  T->SetBranchAddress("fac1",&fac1);
  T->SetBranchAddress("fac2",&fac2);
  T->SetBranchAddress("fac3",&fac3);
  double min_mm = -300;
  double max_mm =  200;
  int bin_mm     = (int)(max_mm - min_mm)/2; //2 MeV Bin
  TH1D* hmm      = new TH1D("hmm","",bin_mm, min_mm,max_mm);
  TH1D* hmm_L    = new TH1D("hmm_L","",bin_mm, min_mm,max_mm);
  TH1D* hmm_nnL  = new TH1D("hmm_nnL","",bin_mm, min_mm,max_mm);
  TH1D* hmm_fsi1 = new TH1D("hmm_fsi1","",bin_mm, min_mm,max_mm);
  TH1D* hmm_fsi2 = new TH1D("hmm_fsi2","",bin_mm, min_mm,max_mm);
  TH1D* hmm_fsi3 = new TH1D("hmm_fsi3","",bin_mm, min_mm,max_mm);

  int bin_mm_10keV = (int)(max_mm - min_mm)*100;  // 10 keV Bin
  TH1D* hmm_10keV      = new TH1D("hmm","",bin_mm_10keV, min_mm,max_mm);
  TH1D* hmm_L_10keV    = new TH1D("hmm_L","",bin_mm_10keV, min_mm,max_mm);
  TH1D* hmm_nnL_10keV  = new TH1D("hmm_nnL","",bin_mm_10keV, min_mm,max_mm);
  TH1D* hmm_fsi1_10keV = new TH1D("hmm_fsi1","",bin_mm_10keV, min_mm,max_mm);
  TH1D* hmm_fsi2_10keV = new TH1D("hmm_fsi2","",bin_mm_10keV, min_mm,max_mm);
  TH1D* hmm_fsi3_10keV = new TH1D("hmm_fsi3","",bin_mm_10keV, min_mm,max_mm);
  
  

  // === Momentum cut ======= //
  double Lp_min = 2120;  // MeV
  //  double Lp_max = 2250;  // MeV
  double Lp_max = 2500;  // MeV
  double Rp_min = 1730;  // MeV
  double Rp_max = 1900;  // MeV
  // ======================= //

  //  double mom_cut = 4050.;
  
  for(int i=0; i<ENum;i++){
    
    T->GetEntry(i);
    // momentum cut //

    if((Lp_min < momL && momL < Lp_max)
       && (Rp_min < momR && momR < Rp_max)){
      
      hmm->Fill(mm_nnL);
      hmm_L->Fill(mm_L);
      hmm_nnL->Fill(mm_nnL);
      hmm_fsi1->Fill(mm_nnL,fac1);
      hmm_fsi2->Fill(mm_nnL,fac2);
      hmm_fsi3->Fill(mm_nnL,fac3);

      hmm_10keV->Fill(mm_nnL);
      hmm_L_10keV->Fill(mm_L);
      hmm_nnL_10keV->Fill(mm_nnL);
      hmm_fsi1_10keV->Fill(mm_nnL,fac1);
      hmm_fsi2_10keV->Fill(mm_nnL,fac2);
      hmm_fsi3_10keV->Fill(mm_nnL,fac3);      
      
    }
    
  }

 
  //    string ofrname ="../rootfiles/fsi/mmass/nnL_exp_hist_momL_cut.root";
  //  string ofrname ="../rootfiles/fsi/mmass/test.root";
  //    ofstream* ofr =new ofstream(ofrname.c_str());

  
  string ofrname;
  if(nnL) ofrname ="../rootfiles/fsi/mmass/nnL_Jl0_mom_cut.root";
  else    ofrname ="../rootfiles/fsi/mmass/Lam_Jl0_mom_cut.root";
  
  TFile* ofr = new TFile(ofrname.c_str(),"recreate");


  hmm->Write();
  hmm_L->Write();
  hmm_nnL->Write();
  hmm_fsi1->Write();
  hmm_fsi2->Write();
  hmm_fsi3->Write();
  ofr->Clone();
  
  string ofrname2;
  if(nnL) ofrname2 ="../rootfiles/fsi/mmass/nnL_Jl0_10keV_mom_cut.root";
  else    ofrname2 ="../rootfiles/fsi/mmass/Lam_Jl0_10keV_mom_cut.root";
  TFile* ofr2 = new TFile(ofrname2.c_str(),"recreate");

  hmm_10keV->Write();
  hmm_L_10keV->Write();
  hmm_nnL_10keV->Write();
  hmm_fsi1_10keV->Write();
  hmm_fsi2_10keV->Write();
  hmm_fsi3_10keV->Write();  

  ofr2->Clone();
  

  cout<<"---------- Infomation ------------"<<endl;
  cout<<"input root : "<<rname<<endl;
  cout<<"root file : "<<ofrname<<endl;
  cout<<"root file 10 keV: "<<ofrname2<<endl;
  
}
