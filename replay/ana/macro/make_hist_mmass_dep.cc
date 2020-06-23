//#include <Setting.h>


void make_hist_mmass_dep(){

  //  Setting* set=new Setting();
  string ifname ="../rootfiles/mmass/ana_Lambda/mmcalib_new/Lambda_small_OleH_all.root";

  TFile*ifp = new TFile(ifname.c_str());
  TChain* T =new TChain("T");
  T->Add(ifname.c_str());
  double mm_L,mm_nnL,Lph;

  T->SetBranchStatus("*",0);
  T->SetBranchStatus("mm_L",1);
  T->SetBranchAddress("mm_L",&mm_L);
  T->SetBranchStatus("mm_nnL",1);
  T->SetBranchAddress("mm_nnL",&mm_nnL);  
  T->SetBranchStatus("Lph",1);
  T->SetBranchAddress("Lph",&Lph);  
  
  int ENum=T->GetEntries();
  cout<<"ENum : "<<ENum<<endl;

  double min_phi,max_phi;
  min_phi = -0.03;
  max_phi = 0.03;
  int bin_phi= 100;
  double min_mmL,max_mmL;
  min_mmL = -50.;
  max_mmL = 150.;
  int bin_mmL = 100;
  double min_mm_nnL,max_mm_nnL;
  min_mm_nnL = -50.;
  max_mm_nnL = 250.;
  int bin_mm_nnL = 150;  
  
  TH2D* hLam_Lphi = new TH2D("hLam_Lphi","Hkine vs Lphi ; Lphi [rad]; M_{X}-M_{#Lambda} [MeV]",bin_phi,min_phi,max_phi, bin_mmL,min_mmL,max_mmL);
  TH2D* hnnL_Lphi = new TH2D("hnnL_Lphi","Tkine vs Lphi ; Lphi [rad]; M_{X}-M_{#Lambda} [MeV]",bin_phi,min_phi,max_phi, bin_mm_nnL,min_mm_nnL,max_mm_nnL);
  

  for(int nev=0;nev<ENum;nev++){
    T->GetEntry(nev);
    hLam_Lphi->Fill(Lph,mm_L);
    hnnL_Lphi->Fill(Lph,mm_nnL);
  }

  string ofname ="../rootfiles/macro/mm_vs_Lphi.root";
  TFile* ofp =new TFile(ofname.c_str(),"recreate");
  
  hLam_Lphi->Write();
  hnnL_Lphi->Write();





  cout<<"new Rootfile : "<<ofname<<endl;
  gSystem->Exit(1);
  
}
