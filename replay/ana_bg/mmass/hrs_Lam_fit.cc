
void hrs_Lam_fit(){
  
  //  TFile *fin=new TFile("./../../rootfiles/Lambda_small1_mm_ana1222.root");
  TFile *fin=new TFile("./../rootfiles/mmass/ana_Lambda/Lambda_small_H_mm.root");  
  // TTree * T=(TTree*)fin->Get("T");
  TH1D* hmm=(TH1D*)fin->Get("h_mm");
  TH1D* hmm_acc=(TH1D*)fin->Get("h_mm_acc");
  TH1D* hpeak=new TH1D("hpeak","Missing Mass w/o AC Cut w/o ACC",500,0.5,1.5); 
  TH1D* hmm_Al=(TH1D*)fin->Get("h_mm_Al");
  hpeak->Add(hmm,hmm_acc,1.,-1.);

  bool bg_pol2=true;
  bool bg_gausn=false;
  TF1* fbg_L;

  fbg_L=new TF1("fbg_L","pol1(0)",1.02,1.3);  
  fbg_L->FixParameter(0,-1137.13); 
  fbg_L->FixParameter(1,1071.8);

  hpeak->Fit("fbg_L","","",1.09,1.13);
  double pbg[3];
  pbg[0]=fbg_L->GetParameter(0);
  pbg[1]=fbg_L->GetParameter(1);

  
  double min_L=1.105;
  double max_L=1.13;
  double pL[3],pL_err[3];
    
  TF1* fpeak_L=new TF1("fpeak_L","gausn(0)",0.5,1.5);  
  fpeak_L->SetParameter(0,3.07090e+00); 
  fpeak_L->SetParameter(1,1.11546e+00);
  fpeak_L->SetParameter(2,3.6814e-03);
  hpeak->Fit("fpeak_L","","",min_L,max_L);
  pL[0]=fpeak_L->GetParameter(0);
  pL[1]=fpeak_L->GetParameter(1); 
  pL[2]=fpeak_L->GetParameter(2);  
  

  TF1* fLam=new TF1("fLam","pol1(0)+gausn(2)",min_L,max_L);
  fLam->SetParameter(0,fbg_L->GetParameter(0)); 
  fLam->SetParameter(1,fbg_L->GetParameter(1));
  //  fLam->FixParameter(0,fbg_L->GetParameter(0));
  //  fLam->FixParameter(1,fbg_L->GetParameter(1));
  fLam->SetParameter(2,pL[0]); 
  fLam->SetParameter(3,pL[1]);
  fLam->SetParameter(4,pL[2]);

  hpeak->Fit("fLam","","",min_L,max_L);
  pbg[0]=fLam->GetParameter(0);
  pbg[1]=fLam->GetParameter(1);
  pL[0]=fLam->GetParameter(2);
  pL[1]=fLam->GetParameter(3); 
  pL[2]=fLam->GetParameter(4);  
  pL_err[0]=fLam->GetParError(2);
  pL_err[1]=fLam->GetParError(3); 
  pL_err[2]=fLam->GetParError(4);


 
  fpeak_L->SetParameters(pL[0],pL[1],pL[2]);
  fbg_L->SetParameters(pbg[0],pbg[1]);

  TCanvas* c0=new TCanvas("c0","MM Hist");
  c0->cd();
  hpeak->Draw();
  //  fbg_L->SetLineColor(kBlue);
  fLam->SetLineColor(kRed);
  fpeak_L->SetLineColor(kBlue);  
  //  fbg_L->GetXaxis()->SetRangeUser(min_L,max_L);
  fLam ->GetXaxis()->SetRangeUser(min_L,max_L);
  //  fpeak_L ->GetXaxis()->SetRangeUser(min_L,max_L);
  fpeak_L->SetFillStyle(3002);
  fpeak_L->SetFillColor(kBlue);
  fpeak_L->SetNpx(2000);
  fLam->Draw("same");
  //  fbg_L->Draw("same");
  fpeak_L->Draw("same");

  TCanvas* c1=new TCanvas("c1","MM Hist");
  c1->cd();  
  hmm->Draw();
  hmm_acc->Draw("same");
  fpeak_L->Draw("same");
  
  cout<<"Lambda events: "<<pL[0]/0.002<<endl;



  
  
}
