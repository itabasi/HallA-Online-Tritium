
void hrs_peak_fit_Ltuned(){
  
 
 
  TFile *fin=new TFile("./../rootfiles/mmass/ana_Lambda/Lambda_small_H_Ltuned_mm.root");  //0th tuning
  //  TH1D* hmm=(TH1D*)fin->Get("h_mm");
  //  TH1D* hmm_acc=(TH1D*)fin->Get("h_mm_acc");
  TH1D* hmm=(TH1D*)fin->Get("h_mm_L");
  TH1D* hmm_acc=(TH1D*)fin->Get("h_acc_L");
  TH1D* hpeak=new TH1D("hpeak","Missing Mass w/ AC Cut w/o ACC",500,0.5,1.5); 
  TH1D* hmm_Al=(TH1D*)fin->Get("h_mm_Al");
  hpeak->Add(hmm,hmm_acc,1.,-1.);

  bool bg_pol2=true;
  bool bg_gausn=false;


  //======== Lambda B.G. ============//
  TF1*  fbg_L=new TF1("fbg_L","pol1(0)",1.1,1.3);  
  //  fbg_L->FixParameter(0,-1137.13); 
  //  fbg_L->FixParameter(1,1071.8);
  fbg_L->SetParameter(0,-1137.13); 
  fbg_L->SetParameter(1,1071.8);

  hpeak->Fit("fbg_L","","",1.12,1.19);
  double pbg[3];
  pbg[0]=fbg_L->GetParameter(0);
  pbg[1]=fbg_L->GetParameter(1);


  //======== Sigma B.G. ==============//
  TF1*  fbg_S=new TF1("fbg_S","pol1(0)",1.17,1.26);  
  //  fbg_S->FixParameter(0,448.634); 
  //  fbg_S->FixParameter(1,-352.65);
  fbg_S->SetParameter(0,448.634); 
  fbg_S->SetParameter(1,-352.65);
  hpeak->Fit("fbg_S","","",1.20,1.27);
  double pbg_S[3];
  pbg_S[0]=fbg_S->GetParameter(0);
  pbg_S[1]=fbg_S->GetParameter(1);

  


  //======== Lambda Peak ============//
  double min_L=1.125;
  double max_L=1.170;
  double pL[3],pL_err[3];
    
  TF1* fpeak_L=new TF1("fpeak_L","gausn(0)",0.5,1.5);  
  fpeak_L->SetParameter(0,1.89942e+00); 
  fpeak_L->SetParameter(1,1.14900e+00);
  fpeak_L->SetParameter(2,6.74142e-03);
  hpeak->Fit("fpeak_L","","",min_L,max_L);
  pL[0]=fpeak_L->GetParameter(0);
  pL[1]=fpeak_L->GetParameter(1); 
  pL[2]=fpeak_L->GetParameter(2);  

  pL_err[0]=fpeak_L->GetParError(0);
  pL_err[1]=fpeak_L->GetParError(1); 
  pL_err[2]=fpeak_L->GetParError(2);  

  //======== Sigma Peak ============//
  double min_S=1.21;
  double max_S=1.24;
  double pS[3],pS_err[3];
    
  TF1* fpeak_S=new TF1("fpeak_S","gausn(0)",0.5,1.5);  
  fpeak_S->SetParameter(0,6.00984e-01); 
  fpeak_S->SetParameter(1,1.22777e+00);
  fpeak_S->SetParameter(2,5.98084e-03);
  hpeak->Fit("fpeak_S","","",min_S,max_S);

  pS[0]=fpeak_S->GetParameter(0);
  pS[1]=fpeak_S->GetParameter(1); 
  pS[2]=fpeak_S->GetParameter(2);  

  pS_err[0]=fpeak_S->GetParError(0);
  pS_err[1]=fpeak_S->GetParError(1); 
  pS_err[2]=fpeak_S->GetParError(2);  


  //======== Lambda Peak + B.G. =============//

  TF1* fLam=new TF1("fLam","pol1(0)+gausn(2)",min_L,max_L);
  fLam->SetParameter(0,fbg_L->GetParameter(0)); 
  fLam->SetParameter(1,fbg_L->GetParameter(1));
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


  
  //======== Sigma Peak + B.G. =============//
  TF1* fSig=new TF1("fSig","pol1(0)+gausn(2)",min_S,max_S);
  fSig->SetParameter(0,fbg_S->GetParameter(0)); 
  fSig->SetParameter(1,fbg_S->GetParameter(1));
  fSig->SetParameter(2,pS[0]); 
  fSig->SetParameter(3,pS[1]);
  fSig->SetParameter(4,pS[2]);

  hpeak->Fit("fSig","","",min_S,max_S);
  pbg_S[0]=fSig->GetParameter(0);
  pbg_S[1]=fSig->GetParameter(1);
  pS[0]=fSig->GetParameter(2);
  pS[1]=fSig->GetParameter(3); 
  pS[2]=fSig->GetParameter(4);  
  pS_err[0]=fSig->GetParError(2);
  pS_err[1]=fSig->GetParError(3); 
  pS_err[2]=fSig->GetParError(4);


 
  fpeak_S->SetParameters(pS[0],pS[1],pS[2]);
  fbg_S->SetParameters(pbg_S[0],pbg_S[1]);



  double bin_min_mmL,bin_max_mmL;
  bin_min_mmL=hmm->GetXaxis()->FindBin(pL[1]-3*pL[2]);
  bin_max_mmL=hmm->GetXaxis()->FindBin(pL[1]+3*pL[2]);
  double sum_L=hmm->Integral(bin_min_mmL,bin_max_mmL);

  double bin_min_mmS,bin_max_mmS;
  bin_min_mmS=hmm->GetXaxis()->FindBin(pS[1]-3*pS[2]);
  bin_max_mmS=hmm->GetXaxis()->FindBin(pS[1]+3*pS[2]);
  double sum_S=hmm->Integral(bin_min_mmS,bin_max_mmS);


  
  
  TCanvas* c0=new TCanvas("c0","MM Hist");
  c0->cd();
  hpeak->Draw();
  fLam->SetLineColor(kRed);
  fpeak_L->SetLineColor(kBlue);  
  fLam ->GetXaxis()->SetRangeUser(min_L,max_L);
  fpeak_L->SetFillStyle(3002);
  fpeak_L->SetFillColor(kBlue);
  fpeak_L->SetNpx(2000);

  fSig->SetLineColor(kRed);
  fpeak_S->SetLineColor(kBlue);  
  fSig ->GetXaxis()->SetRangeUser(min_S,max_S);
  fpeak_S->SetFillStyle(3002);
  fpeak_S->SetFillColor(kBlue);
  fpeak_S->SetNpx(2000);
  
  fLam->Draw("same");
  fpeak_L->Draw("same");
  fpeak_S->Draw("same");  

  
  TCanvas* c1=new TCanvas("c1","MM Hist");
  c1->cd();
  hmm->SetTitle("Missing mass w/ z and coin cut");
  hmm->GetXaxis()->SetRangeUser(1.0,1.3);
  //  hmm->GetYaxis()->SetRangeUser(0.0,450);    
  hmm->Draw();
  hmm_acc->SetLineColor(2);
  hmm_acc->Draw("same");
  fpeak_L->Draw("same");
  fpeak_S->Draw("same");  


 

  
  cout<<"======= Lambda ======="<<endl;
  cout<<"Lambda events: "<<pL[0]/0.002<<endl;
  cout<<"Lambda events Error : "<<pL_err[0]/0.002<<endl;  
  cout<<"total events: "<<sum_L<<endl;
  cout<<"Lambda mean: " <<pL[1]<<endl;
  cout<<"Lambda mean Error: " <<pL_err[1]<<endl;  
  cout<<"Lambda sigma: "<<pL[2]<<endl;
  cout<<"Lambda sigma Error: "<<pL_err[2]<<endl;  
  cout<<"Peak Signifiance S/sq(S+N): "<<pL[0]/0.002/sqrt(sum_L)<<endl;
  cout<<"======= Sigma ======="<<endl;
  cout<<"Sigma events: "<<pS[0]/0.002<<endl;
  cout<<"Sigma events Error: "<<pS_err[0]/0.002<<endl;  
  cout<<"total events: "<<sum_S<<endl;  
  cout<<"Sigma mean: "<<pS[1]<<endl;
  cout<<"Sigma mean Error: "<<pS_err[1]<<endl;  
  cout<<"Sigma sigma: "<<pS[2]<<endl;
  cout<<"Sigma sigma Error: "<<pS_err[2]<<endl;  
  cout<<"Peak Significance : "<<pS[0]/0.002/sqrt(sum_S)<<endl;
  
  
}
