

void mmass_vs_phi(){

  string ifname = "../rootfiles/macro/mm_vs_Lphi.root";
  TFile* ifp = new TFile(ifname.c_str());


  //---------- YSlieces in Lambda -------------//
  double min_mmL=-50;
  double max_mmL=150;
  int   bin_mmL=100;
  TH2F* hLam_Lphi =((TH2F*)ifp->Get("hLam_Lphi"));
  TF1* fslice_phi=new TF1("fslice_phi","gausn(0)",-10,10);
  TH1F* hslice_phi;
  double min_phi,max_phi;
  min_phi=-0.05,max_phi=0.05;
  TF1* ffit_phi=new TF1("ffit_phi","[0]*x +[1]",min_phi,max_phi);
  //---------- YSlieces in nnL -------------//  
  double min_mm_nnL=-50;
  double max_mm_nnL=250;
  int   bin_mm_nnL=150;
  TH2F* hnnL_Lphi =((TH2F*)ifp->Get("hnnL_Lphi"));
  TF1* fslice_phi_nnL=new TF1("fslice_phi_nnL","gausn(0)",30,80);
  TH1F* hslice_phi_nnL;
  TF1* ffit_phi_nnL=new TF1("ffit_phi_nnL","[0]*x +[1]",min_phi,max_phi);  


  ffit_phi->SetLineColor(2);
  hLam_Lphi->FitSlicesY(fslice_phi,0,-1,0,"QRG2");

  hslice_phi=(TH1F*)gROOT->FindObject("hLam_Lphi_1");
  hslice_phi->Fit("ffit_phi","","",min_phi,max_phi);


  ffit_phi_nnL->SetLineColor(2);
  hnnL_Lphi->FitSlicesY(fslice_phi_nnL,0,-1,0,"QRG2");

  hslice_phi_nnL=(TH1F*)gROOT->FindObject("hnnL_Lphi_1");
  hslice_phi_nnL->Fit("ffit_phi_nnL","","",min_phi,max_phi);


  TCanvas* c0=new TCanvas("c0","c0");
  c0->Divide(2,1);
  c0->cd(1);
  hLam_Lphi->Draw("colz");
  ffit_phi->Draw("same");

  c0->cd(2);
  hnnL_Lphi->Draw("colz");
  ffit_phi_nnL->Draw("same");  
}


