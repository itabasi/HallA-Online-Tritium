

void makehist(){

  //  string ifname="../rootfiles/mmass/ana_Lambda/2020-02-25_3/Lambda_small_OleH_all.root";
  string ifname="../rootfiles/mmass/ana_Lambda/2020-02-25_3/nnL_small_Ole_all.root";
  string ofname="nnLhist.root";


  TFile* f=new TFile(ifname.c_str());
  TTree* T=new TTree("T","");

  //==========   ============== =============//

  TH1D* hLam=(TH1D*)f->Get("h_mm_L");
  TH1D* hLam_acc=(TH1D*)f->Get("h_acc_L");
  TH1D* hnnL=(TH1D*)f->Get("h_mm_nnL");
  TH1D* hnnL_acc=(TH1D*)f->Get("h_acc_nnL");

  //==========   ============== =============//



  //  TH1D* hLam_acc=(TH1D*)h_acc_L->Clone();
  //  TH1D* hnnL=(TH1D*)h_mm_nnL->Clone();
  //  TH1D* hnnL_acc=(TH1D*)h_acc_nnL->Clone();


  TFile* fnew=new TFile(ofname.c_str(),"recreate");


  hLam->Write();
  hLam_acc->Write();
  hnnL->Write();
  hnnL_acc->Write();



}
