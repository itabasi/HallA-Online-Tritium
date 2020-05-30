

void Al_hist(){


  string fname="../../rootfiles/mmass/ana_Lambda/mmcalib_new/nnL_small_Ole_all.root";

  TFile* ifs=new TFile(fname.c_str());
  TH1F* hMgL=(TH1F*)ifs->Get("h_mm_MgL");
  TH1F* hMgL_acc=(TH1F*)ifs->Get("h_mm_MgL_acc");
  TH1F* hMgL_peak=(TH1F*)ifs->Get("h_peak_MgL");
  TH1F* hMgL_100keV=new TH1F("hMgL_100keV","27MgL hist 100 keV bin",6000,-300,300);
  TChain* T=new TChain("T");
  T->Add(fname.c_str());
  double mm_MgL;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("mm_MgL",1);
  T->SetBranchAddress("mm_MgL",&mm_MgL);

  string ofname ="./Al_hist.root";
  TFile* ofs=new TFile(ofname.c_str(),"recreate");  
  TTree* tnew = new TTree("T","Al events");
  tnew =T->CloneTree(1);  
  int ENum=T->GetEntries();

  for(int i=0;i<ENum;i++){
    T->GetEntry(i);
    hMgL_100keV->Fill(mm_MgL);
    
  }

  hMgL->Write();
  hMgL_acc->Write();
  hMgL_peak->Write();
  hMgL_100keV->Write();
  tnew->Write();
}
