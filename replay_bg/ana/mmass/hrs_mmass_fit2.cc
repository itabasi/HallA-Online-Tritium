
void hrs_mmass_fit2(){
  
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
  TF1* fbg;
  if(bg_gausn){
  fbg=new TF1("fbg","gausn(0)",0.5,1.5);  
  fbg->FixParameter(0,4.82437e+00); 
  fbg->FixParameter(1,1.12980e+00);
  fbg->FixParameter(2,3.28259e-02);
  }
  if(bg_pol2){
  fbg=new TF1("fbg","pol2(0)",1.02,1.3);  
  fbg->FixParameter(0,-1724.51); 
  fbg->FixParameter(1,2962.51);
  fbg->FixParameter(2,-1266.79);
  }



  hpeak->Fit("fbg","","",1.0,1.4);
  double pbg[3];
  pbg[0]=fbg->GetParameter(0);
  pbg[1]=fbg->GetParameter(1);
  pbg[2]=fbg->GetParameter(2);
  
  TF1* fLam=new TF1("fLam","fbg+gausn(3)",0.5,1.5);
  fLam->SetParameter(0,fbg->GetParameter(0)); 
  fLam->SetParameter(1,fbg->GetParameter(1));
  fLam->SetParameter(2,fbg->GetParameter(2));
 //  fLam->FixParameter(0,fbg->GetParameter(0)); 
  //  fLam->SetParLimits(0,0.5*fbg->GetParameter(0),1.1*fbg->GetParameter(0)); 
  fLam->FixParameter(1,fbg->GetParameter(1));
  fLam->FixParameter(2,fbg->GetParameter(2));
  fLam->SetParameter(3,2.35002e+00); 
  fLam->SetParameter(4,1.11560e+00);
  fLam->SetParameter(5,4.67814e-03);
  hpeak->Fit("fLam","","",1.10,1.13);
  double pL[3],pL_err[3];
  pL[0]=fLam->GetParameter(3);
  pL[1]=fLam->GetParameter(4); 
  pL[2]=fLam->GetParameter(5);  
  pL_err[0]=fLam->GetParError(3);
  pL_err[1]=fLam->GetParError(4); 
  pL_err[2]=fLam->GetParError(5);  
  TF1* fSig=new TF1("fSig","fbg+gausn(3)",0.5,1.5);
  fSig->FixParameter(0,fbg->GetParameter(0)); 
  //  fSig->SetParameter(0,fbg->GetParameter(0)); 
  fSig->FixParameter(1,fbg->GetParameter(1));
  fSig->FixParameter(2,fbg->GetParameter(2));
  fSig->SetParameter(3,2.57414e+00); 
  fSig->SetParameter(4,1.20040e+00);
  fSig->SetParameter(5,1.27321e-02);
  hpeak->Fit("fSig","","",1.19,1.21);
  // hpeak->Fit("fSig","","",1.1,1.15);
  double pS[3],pS_err[3];
  pS[0]=fSig->GetParameter(3);
  pS[1]=fSig->GetParameter(4); 
  pS[2]=fSig->GetParameter(5);
  pS_err[0]=fSig->GetParError(3);
  pS_err[1]=fSig->GetParError(4); 
  pS_err[2]=fSig->GetParError(5);


  TF1* fpeak_L=new TF1("fpeak_L","gausn",0.5,1.5);
  fpeak_L->SetParameters(pL[0],pL[1],pL[2]);
  TF1* fpeak_S=new TF1("fpeak_S","gausn",0.5,1.5);
  fpeak_S->SetParameters(pS[0],pS[1],pS[2]);



  TF1* fmm=new TF1("fmm","gausn(0)+gausn(3)+gausn(6)",0.5,1.5);
  fmm->SetParameters(pL[0],pL[1],pL[2],pS[0],pS[1],pS[2],pbg[0],pbg[1],pbg[2]);
  fmm->FixParameter(0,pL[0]);
  fmm->FixParameter(1,pL[1]);
  fmm->FixParameter(2,pL[2]);
  fmm->FixParameter(3,pS[0]);
  fmm->FixParameter(4,pS[1]);
  fmm->FixParameter(5,pS[2]);
  fmm->FixParameter(7,pbg[1]);
  fmm->FixParameter(8,pbg[2]);



  fLam->FixParameter(0,pL[0]);
  fLam->FixParameter(1,pL[1]);
  fLam->FixParameter(2,pL[2]);
  fSig->FixParameter(0,pS[0]);  
  fSig->FixParameter(1,pS[1]);
  fSig->FixParameter(2,pS[2]);

  //  hpeak->Fit("fmm","","",0.5,1.5);
  pL[0]=fmm->GetParameter(0);
  pL[1]=fmm->GetParameter(1); 
  pL[2]=fmm->GetParameter(2);  
  pS[0]=fmm->GetParameter(3);
  pS[1]=fmm->GetParameter(4); 
  pS[2]=fmm->GetParameter(5);


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
  fbg->SetLineColor(kBlue);
  fLam->SetLineColor(kRed);
  fSig->SetLineColor(kGreen);
  fLam->Draw("same");
  fSig->Draw("same");
  fbg->Draw("same");
  
 
  TCanvas* c1=new TCanvas("c1","MM Fit");
  c1->cd();
  hpeak->Draw();
  fpeak_L->SetLineColor(4);
  fpeak_L->SetFillColor(4);
  fpeak_L->SetFillStyle(3001);
  fpeak_S->SetLineColor(3);
  fpeak_S->SetFillColor(3);
  fpeak_S->SetFillStyle(3001);
  fpeak_L->SetNpx(2000);
  fpeak_S->SetNpx(2000);
  fpeak_L->Draw("same");
  fpeak_S->Draw("same");
  fmm->SetLineColor(2);
  fmm->SetNpx(2000);
  fmm->Draw("same");

  
  TCanvas* c2=new TCanvas("c2","MM ACC");
  c2->cd();
  hmm->SetLineColor(2);
  hmm->SetFillColor(2);
  hmm->SetFillStyle(3001);
  hmm->Draw();
  hmm_acc->SetLineColor(4);
  hmm_acc->SetFillColor(4);
  hmm_acc->SetFillStyle(3001);
  hmm_acc->Draw("same");
 

  //======== COMMENT OUT ==========//

  cout<<"========== <COMMENT OUT> ============= "<<endl;

  cout<<"Lambda Peak "<<endl;
  cout<<"Number of Lambda : "<<sum_L<<endl;
  cout<<"Number of Lambda Peak: "<<pL[0]/0.002<<endl;
  cout<<"Number of Lambda ACC: "<<sum_L-pL[0]/0.002<<endl;
  cout<<"Lambda S/N: "<<pL[0]/0.002/(sum_L-pL[0]/0.002)<<endl;
  cout<<"Mean of Lambda : "<<pL[1]<<"+-"<<pL_err[1]<<" GeV/c^2 "<<endl;
  cout<<"Sigma of Lambda : "<<pL[2]<<"+-"<<pL_err[2]<<" GeV/c^2 "<<endl; 
  cout<<"Sigma Peak "<<endl;
  cout<<"Number of Sigma : "<<sum_S<<endl;
  cout<<"Number of Sigma Peak: "<<pS[0]/0.002<<endl;
  cout<<"Number of Sigma ACC: "<<sum_S-pS[0]/0.002<<endl;
  cout<<"Sigma S/N: "<<pS[0]/0.002/(sum_S-pS[0]/0.002)<<endl;
  cout<<"Mean of Sigma : "<<pS[1]<<"+-"<<pS_err[1]<<" GeV/c^2 "<<endl;
  cout<<"Sigma of Sigma : "<<pS[2]<<"+-"<<pS_err[2]<<" GeV/c^2 "<<endl;
 
}
