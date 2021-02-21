#include "TChain.h"

//bool MOMCUT= false;
bool MOMCUT = true;
bool GCOIN  = false; 
void MakeHist(){
  
  //  string  rname ="../rootfiles/mmass/ana_Lambda/mmcalib_nmr_small/nnL_small_Ole_all.root";
  //  string  rname ="../rootfiles/mmass/ana_Lambda/test/nnL_small_Ole_all.root";
  //string  rname ="../rootfiles/mmass/ana_Lambda2/mmcalib_coin_small/nnL_small_Ole_all.root"; // coin


  
  string  rname;


  if(GCOIN) rname ="../rootfiles/mmass/ana_Lambda2/mmcalib_gcoin_small/nnL_small_Ole_all.root"; // gcoin
  else      rname ="../rootfiles/mmass/ana_Lambda2/mmcalib_coin_small/nnL_small_Ole_all.root"; // coin

  cout<<"=== MODE =========="<<endl;
  cout<<"MOMCUT : "<<MOMCUT<<endl;
  cout<<"GCOIN : "<<GCOIN<<endl;
  


  TChain* T =new TChain("T");
  T->Add(rname.c_str());
  
  int ENum = T->GetEntries();
  double mm_nnL,mm_acc[10],momL,momR,acc_scale;
  T->SetBranchAddress("mm_nnL",&mm_nnL);
  T->SetBranchAddress("mm_mix",mm_acc);
  T->SetBranchAddress("Lp",&momL);
  T->SetBranchAddress("Rp",&momR);
  T->SetBranchAddress("acc_scale",&acc_scale);
  int pid_cut,ct_cut,z_cut,acc_cut;
  T->SetBranchAddress("z_cut",&z_cut);
  T->SetBranchAddress("pid_cut",&pid_cut);
  T->SetBranchAddress("ct_cut",&ct_cut);
  T->SetBranchAddress("acc_cut",&acc_cut);

  double mm_nnL_all,mm_L_all;
  T->SetBranchAddress("mm_nnL_all",&mm_nnL_all);
  T->SetBranchAddress("mm_L_all",&mm_L_all);
  
  double min_mm = -300;
  double max_mm =  200;
  //  int bin_mm = (int)(max_mm - min_mm)*100;
  int bin_mm = (int)(max_mm - min_mm)/2; //2 MeV Bin
  TH1D* hmm = new TH1D("hmm","",bin_mm, min_mm,max_mm);
  TH1D* hmm_test = new TH1D("hmm_test","",bin_mm, min_mm,max_mm);
  TH1D* hmm_peak = new TH1D("hmm_peak","",bin_mm, min_mm,max_mm);
  TH1D* hmm_acc = new TH1D("hmm_acc","",bin_mm, min_mm,max_mm);
  
  TH1D* h_peak_L = new TH1D("h_peak_L","",bin_mm, min_mm,max_mm);
  TH1D* h_peak_nnL = new TH1D("h_peak_nnL","",bin_mm, min_mm,max_mm);
  TH1D* h_acc_L = new TH1D("h_acc_L","",bin_mm, min_mm,max_mm);
  TH1D* h_acc_nnL = new TH1D("h_acc_nnL","",bin_mm, min_mm,max_mm);
  TH1D* h_mm_L = new TH1D("h_mm_L","",bin_mm, min_mm,max_mm);
  TH1D* h_mm_nnL = new TH1D("h_mm_nnL","",bin_mm, min_mm,max_mm);


  // 10 keV Bin Hist //
  int bin_mm_10keV = (int)(max_mm - min_mm)*100; //10 keV Bin
  TH1D* hmm_10keV = new TH1D("hmm","",bin_mm_10keV, min_mm,max_mm);
  TH1D* hmm_peak_10keV = new TH1D("hmm_peak","",bin_mm_10keV, min_mm,max_mm);
  TH1D* hmm_acc_10keV = new TH1D("hmm_acc","",bin_mm_10keV, min_mm,max_mm);  
  

  //==== Momentum Cut ==========//

  double Lp_min = 2.12;
  double Lp_max = 2.5;
  double Rp_min = 1.73;
  double Rp_max = 1.9;

  
  //  double mom_cut = 4.05;
  // TF1* fcut = new TF1("fcut","[0]*x+[1]",Lp_min,Lp_max);
  //  fcut->SetParameters(4.05,-1.0);
  
  
  for(int i=0; i<ENum;i++){
    
    T->GetEntry(i);
    hmm_test->Fill(mm_nnL);
    //    if(momL>0 && momL<Lp_cut){
    if((Lp_min < momL && momL<Lp_max && Rp_min <momR && momR <Rp_max) // MOM CUT
       || !MOMCUT ){ // w/o MOMCUT
      hmm->Fill(mm_nnL);
      hmm_peak->Fill(mm_nnL);
      hmm_10keV->Fill(mm_nnL);
      hmm_peak_10keV->Fill(mm_nnL);
      if(z_cut && pid_cut){	
       
	if(acc_cut){      // Accidental B.G. Events
	  
	  h_acc_L   -> Fill(mm_L_all,acc_scale); // Weight : ct range 
	  h_acc_nnL -> Fill(mm_nnL_all,acc_scale);
	  
	}else if(ct_cut){          // Coin cut  Events
	  
	  h_mm_L     -> Fill(mm_L_all);
	  h_mm_nnL     -> Fill(mm_nnL_all);
	  h_peak_L   -> Fill(mm_L_all);
	  h_peak_nnL -> Fill(mm_nnL_all);
	  
      }
	

	for(int i=0;i<10;i++){
	  hmm_acc->Fill(mm_acc[i],acc_scale/10.);
	  hmm_acc_10keV->Fill(mm_acc[i],acc_scale/10.);
	}
      }
      
    } // Mom Cut

    
    
  }
  
  
  
  string ofrname; 

  if(MOMCUT && GCOIN)
    ofrname ="../rootfiles/fsi/mmass/nnL_exp_gcoin_hist_momL_cut.root";
  else if(!MOMCUT && GCOIN)
    ofrname ="../rootfiles/fsi/mmass/nnL_exp_gcoin_hist.root";
  else if(MOMCUT && !GCOIN)
    ofrname ="../rootfiles/fsi/mmass/nnL_exp_coin_hist_momL_cut.root";
  else if(!MOMCUT && !GCOIN)
    ofrname ="../rootfiles/fsi/mmass/nnL_exp_gcoin_hist.root";


  
  TFile* ofr = new TFile(ofrname.c_str(),"recreate");
  

  hmm_peak->Add(hmm_acc,-1);
  hmm_peak_10keV->Add(hmm_acc_10keV,-1);
  h_peak_L->Add(h_acc_L,-1);
  h_peak_nnL->Add(h_acc_nnL,-1);

  hmm_test->Write();
  hmm->Write();
  hmm_acc->Write();
  hmm_peak->Write();
  
  h_mm_L->Write();
  h_acc_L->Write();
  h_peak_L->Write();

  h_mm_nnL->Write();
  h_acc_nnL->Write();
  h_peak_nnL->Write();
  

  ofr->Close();

  //=== 10 keV Bin Hist ====//

  string ofrname2;

  if(GCOIN && MOMCUT)
    ofrname2 = "../rootfiles/fsi/mmass/nnL_exp_gcoin_10keV_hist_momL_cut.root";
  else if(GCOIN && !MOMCUT)
    ofrname2 = "../rootfiles/fsi/mmass/nnL_exp_gcoin_10keV_hist.root";
  else if(!GCOIN && MOMCUT)
    ofrname2 = "../rootfiles/fsi/mmass/nnL_exp_coin_10keV_hist_momL_cut.root";
  else if(!GCOIN && !MOMCUT)
    ofrname2 = "../rootfiles/fsi/mmass/nnL_exp_coin_10keV_hist.root";
  
  TFile* ofr2 = new TFile(ofrname2.c_str(),"recreate");

  hmm_10keV->Write();
  hmm_acc_10keV->Write();
  hmm_peak_10keV->Write();

  ofr2->Close();


  cout<<"=======OUTPUT ROOT========="<<endl;
  cout<<"new root : "<<ofrname<<endl;
  cout<<"new root 10keV: "<<ofrname2<<endl;
}
