
void hrs_mmass_fit(){
  

  //  TFile *fin=new TFile("./../../rootfiles/Lambda_small1_mm_ana1223.root");
  TFile *fin=new TFile("./../rootfiles/mmass/ana_Lambda/Lambda_small_H_mm.root");
  TH1D* hmm=(TH1D*)fin->Get("h_mm");
  TH1D* hmm_acc=(TH1D*)fin->Get("h_mm_acc");
  TH1D* hpeak=(TH1D*)fin->Get("h_peak_mm");
  TH1D* hmm_Al=(TH1D*)fin->Get("h_peak_Al");
    //new TH1D("hpeak","Missing Mass w/o AC Cut w/o ACC",500,0.5,1.5); 
  TH1D* hmm_pi=(TH1D*)fin->Get("h_mm_pi");
  //hpeak->Add(h_mm,h_mm_acc,1.,-1.);
  //hmm_Al->Add(h_mm_Al,h_mm_acc,1,-1);
  //  hmm_pi->Add(h_mm_pi,h_mm_acc,1,-1);
  TF1* fpi=new TF1("fpi","[0]*gausn(1)",0.5,1.5);
  TF1* fAl=new TF1("fAl","[0]*pol4(1)",1.0 ,1.3);
  TF1* facc=new TF1("facc","[0]*pol5(1)",1.1,1.3);
 
  facc->SetParameters(1,-7.44848e+06,3.34387e+07,-5.98543e+07,5.33887e+07,-2.37277e+07,4.20303e+06);
  facc->SetParLimits(0,0.5,1.5);
  facc->FixParameter(0,1.0);
  hmm_acc->Fit("facc","","",1.1,1.25);
  double p_acc[7];
  for(int i=0;i<7;i++){
  p_acc[i]=facc->GetParameter(i);}
  fpi->SetParameters(1.,1.78531e+02,1.17365e+00,4.55952e-02);
  fpi->FixParameter(0,1);
  cout<<"Pion Hist Fitting "<<endl;
  hmm_pi->Fit("fpi","","",1.0,1.3);
  double p_pi[4];
  for(int i=0;i<4;i++){
    p_pi[i]=fpi->GetParameter(i);}

  fAl->SetParameters(1,885860,-3.08018e+06,3.99347e+06,-2.28779e+06,488658);
  fAl->SetParLimits(0,0.5,1.5);
  fAl->FixParameter(0,1.0);
  cout<<"Alminum Background fitting "<<endl;
  hmm_Al->Fit("fAl","","",1.0,1.3);
  double p_Al[6];
  for(int i=0;i<6;i++){
    p_Al[i]=fAl->GetParameter(i);}

  TF1* fpi_cp=hmm_pi  ->GetFunction("fpi");
  TF1* fAl_cp=hmm_Al  ->GetFunction("fAl");
  TF1* facc_cp=hmm_acc->GetFunction("facc");    

 

  TF1* fbg=new TF1("fbg","fAl+fpi",1.0,1.3);
  fbg->SetParLimits(0,0.01*p_Al[0],1.0*p_Al[0]);
  for(int i=1;i<6;i++){fbg->FixParameter(i,p_Al[i]);}
  fbg->SetParLimits(6,0.0,5*p_pi[0]);
  for(int i=7;i<10;i++){fbg->FixParameter(i,p_pi[i-6]);}
  hpeak->Fit("fbg","","",1.05,1.3);
  
 //for(int i=7;i<13;i++){fbg->FixParameter(i,p_acc[i-6]);}
  //fbg->SetParLimits(13,0.0,5*p_pi[0]);
  // fbg->FixParameter(13,0.0);
  fAl->FixParameter(0,fbg->GetParameter(0));
  fpi->FixParameter(0,fbg->GetParameter(6));  
//  facc->FixParameter(0,fbg->GetParameter(6));
  
  /*
  ffbg->SetParameter(0,p_acc[0]);
  ffbg->SetParameter(1,p_acc[1]);
  ffbg->SetParameter(2,p_acc[2]);
  ffbg->SetParameter(3,p_acc[3]);
  ffbg->SetParameter(4,p_acc[4]);
  ffbg->SetParameter(5,p_acc[5]);
  

  fbg->SetParLimits(0,0.01,1.);
  fbg->FixParameter(1,1.17986e+00);
  fbg->FixParameter(2,5.22373e-02);
  fbg->FixParameter(4,1.17483e+00);
  fbg->FixParameter(5,5.02931e-02);
  hpeak->Fit("fbg","","",1.0,1.4);
  */

  /*
  TF1* fLam=new TF1("fLam","gaus(0)+gaus(3)+gaus(6)",0.5,1.5);

  fLam->SetParameter(1,1.115);
  fLam->SetParLimits(1,1.14,1.16);
  fLam->SetParameter(2,4.47765e-03);
  fLam->SetParLimits(2,0.0,5.0e-3);
  fLam->FixParameter(3,fbg->GetParameter(0));
  fLam->FixParameter(4,1.17986e+00);
  fLam->FixParameter(5,5.22373e-02);
  fLam->FixParameter(6,fbg->GetParameter(3));
  fLam->FixParameter(7,1.17483e+00);
  fLam->FixParameter(8,5.02931e-02);

  
  hpeak->Fit("fLam","","",1.113,1.117);
  */


 
  TCanvas* c0=new TCanvas("c0","MM Hist");
  c0->cd();
  hpeak->Draw();
  fpi->SetLineColor(3);
  fAl->SetLineColor(4);
  fbg->SetLineColor(2);
  fbg->SetFillColor(2);
  fbg->SetFillStyle(3001);
  fpi->Draw("same");
  fAl->Draw("same");
  fbg->Draw("same");
  

  
  TCanvas* c1=new TCanvas("c1","MM ACC");
  c1->cd();
  hmm->Draw();
  hmm_acc->SetLineColor(4);
  hmm_acc->SetFillColor(4);
  hmm_acc->SetFillStyle(3001);
  hmm_acc->Draw("same");
  

  TCanvas* c2=new TCanvas("c2","Hist");
  c2->Divide(2,2);
  c2->cd(1);
  hmm_acc->Draw();
  facc_cp->SetLineColor(2);
  facc_cp->Draw("same");
  facc->Draw("same");
  c2->cd(2);
  hmm_Al->Draw();
  fAl_cp->SetLineColor(2);
  fAl_cp->Draw("same");
  fAl->Draw("same");
  c2->cd(3);
  hmm_pi->Draw();
  fpi_cp->SetLineColor(2);
  fpi_cp->Draw("same");
  fpi->Draw("same");
  c2->cd(4);
  hpeak->Draw();
  facc_cp->SetLineColor(2);
  fpi_cp->SetLineColor(3);
  fAl_cp->SetLineColor(4);
  fpi->Draw("same");
  fAl->Draw("same");

}
