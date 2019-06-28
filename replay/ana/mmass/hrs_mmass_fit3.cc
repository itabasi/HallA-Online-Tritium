
void hrs_mmass_fit3(){
  

  TFile *fin=new TFile("./../../rootfiles/Lambda_small1_mm_ana1223.root");

  TH1D* hmm=(TH1D*)fin->Get("h_mm");
  TH1D* hmm_acc=(TH1D*)fin->Get("h_mm_acc");
  TH1D* hpeak=(TH1D*)fin->Get("h_peak_mm");
  TH1D* hmm_Al=(TH1D*)fin->Get("h_mm_Al");
  TH1D* hmm_Al_acc=(TH1D*)fin->Get("h_mm_Al_acc");
  TH1D* hpeak_Al=(TH1D*)fin->Get("h_peak_Al");
  TH1D* hmm_pi=(TH1D*)fin->Get("h_mm_pi");
  hmm_pi->Add(h_mm_pi,h_mm_acc,1,-4.0/40.);


  double tcoin_t,Ra1sum,Ra2sum,mm,Lz,Rz,Rpz,Lpz,ct_acc;

  TTree *T=(TTree*)fin->Get("tree");
 T->SetBranchStatus("ct",1);    
 T->SetBranchAddress("ct",&tcoin_t);
 // T->SetBranchAddress("ct",&coin_t);
 T->SetBranchStatus("ac1_sum",1);
 T->SetBranchAddress("ac1_sum",&Ra1sum);
 T->SetBranchStatus("ac2_sum",1);
 T->SetBranchAddress("ac2_sum",&Ra2sum);
 T->SetBranchStatus("mm",1);
 T->SetBranchAddress("mm",&mm);
 T->SetBranchStatus("Lz",1);    
 T->SetBranchAddress("Lz",&Lz);
 T->SetBranchStatus("Rz",1);    
 T->SetBranchAddress("Rz",&Rz);
 T->SetBranchStatus("Rp",1);
 T->SetBranchAddress("Rp",&Rpz);
 T->SetBranchStatus("Lp",1);
 T->SetBranchAddress("Lp",&Lpz);
 T->SetBranchStatus("ct_acc",1);
 T->SetBranchAddress("ct_acc",&ct_acc);





  TF1* fpi=new TF1("fpi","[0]*gausn(1)",1.0,1.3);
  TF1* fAl=new TF1("fAl","[0]*gausn(1)+[4]*gausn(5)",1.0 ,1.3); 
  
  double p_pi[4];
  p_pi[0]=1.0;
  p_pi[1]=3.92404e+01;
  p_pi[2]=1.18391e+00;
  p_pi[3]=3.25563e-02;
  fpi->SetParameters(p_pi[0],p_pi[1],p_pi[2],p_pi[3]);
 fpi->FixParameter(0,1);
  cout<<"=====fpi Fit====="<<endl;
  hmm_pi->Fit("fpi","","",1.0,1.3);
  for(int i=0;i<4;i++){
    p_pi[i]=fpi->GetParameter(i);}

  double p_Al[8];
  p_Al[0]=1.0;
  p_Al[1]=3.45652e+01;
  p_Al[2]=1.13062e+00;
  p_Al[3]=3.85575e-02;
  p_Al[4]=1.0;
  p_Al[5]=2.34929e+00;
  p_Al[6]=1.23090e+00;
  p_Al[7]=2.08151e-02;
  
  fAl->SetParameters(p_Al[0],p_Al[1],p_Al[2],p_Al[3],p_Al[4],p_Al[5],p_Al[6],p_Al[7]);
  cout<<"=====fAl Fit====="<<endl;
  hpeak_Al->Fit("fAl","","",1.0,1.3);

  for(int i=0;i<8;i++){
    p_Al[i]=fAl->GetParameter(i);}






  
 
  //TF1* fbg=new TF1("fbg","[0]*fAl+fpi",1.0,1.3);
    TF1* fbg=new TF1("fbg","[0]*([1]*gausn(2)+[5]*gausn(6)+[9]*gausn(10))",1.0,1.3);

    //  fbg->SetParLimits(0,0.0,1);
  for(int i=1;i<9;i++){fbg->FixParameter(i,p_Al[i-1]);}
  fbg->SetParLimits(9,0.0,10);
  // fbg->FixParameter(9,0.0);
  for(int i=10;i<13;i++){fbg->FixParameter(i,p_pi[i-10]);}
 
  cout<<"======= fbg Fit ========"<<endl;
  hpeak->Fit("fbg","","",1.0,1.3);
  double p_bg[13]; 
  for(int i=0;i<13;i++){p_bg[i]=fbg->GetParameter(i);}


  TF1* fLam=new TF1("fLam","gausn(0)+fbg",1.0,1.3);
  TF1* fSig=new TF1("fSig","gausn(0)+fbg",1.0,1.3);
  for(int i=0;i<13;i++){
  fLam->FixParameter(i+3,p_bg[i]);
  fSig->FixParameter(i+3,p_bg[i]);}
 
  double p_L[3];
  p_L[0]=1.95096;
  p_L[1]=1.11555;
  p_L[2]=4.21431e-03;
  fLam->SetParameter(0,p_L[0]);
  fLam->SetParameter(1,p_L[1]);
  fLam->SetParameter(2,p_L[2]);
  fLam->SetParLimits(0,0.8*p_L[0],1.5*p_L[0]);
  fLam->SetParLimits(1,0.8*p_L[1],1.2*p_L[1]);
  fLam->SetParLimits(2,0.8*p_L[2],1.2*p_L[2]);
  // hpeak->Fit("fLam","","",1.1,1.13);


  p_L[0]=fLam->GetParameter(0);
  p_L[1]=fLam->GetParameter(1);
  p_L[2]=fLam->GetParameter(2);
  

 double p_S[3];
  p_S[0]=5.66467e-01;
  p_S[1]=1.20025e+00;
  p_S[2]=4.81458e-03;
  fSig->SetParameter(0,p_S[0]);
  fSig->SetParameter(1,p_S[1]);
  fSig->SetParameter(2,p_S[2]);
  //  hpeak->Fit("fSig","","",1.19,1.22);


  p_S[0]=fLam->GetParameter(0);
  p_S[1]=fLam->GetParameter(1);
  p_S[2]=fLam->GetParameter(2);
 


 //for(int i=7;i<13;i++){fbg->FixParameter(i,p_acc[i-6]);}
  //fbg->SetParLimits(13,0.0,5*p_pi[0]);
  // fbg->FixParameter(13,0.0);

  //fAl->FixParameter(0,fbg->GetParameter(0));
  // fpi->FixParameter(0,fbg->GetParameter(7));  

  
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
  THStack *hbg=new THStack("hbg","MM Back Ground");  


  TF1* fpi_cp=hmm_pi  ->GetFunction("fpi");
  TF1* fAl_cp=hpeak_Al  ->GetFunction("fAl");
   hpeak_Al->Scale(p_bg[0]);
  fAl_cp->Scale(p_bg[9]);
  hbg->Add(fpi_cp);
  hbg->Add(fAl_cp); 


  TCanvas* c0=new TCanvas("c0","MM Hist");
  c0->cd();
  hpeak->Draw();
  //fSig->SetLineColor(3);
  // fLam->SetLineColor(4);
  //fAl_cp->SetLineColor(1);
  
  fbg->SetLineColor(4);
  fbg->SetFillColor(4);
  fbg->SetFillStyle(3001);
  fAl->Draw("same");
  fbg->Draw("same");

  cout<<"p_BG[0]"<<p_bg[0]<<endl;
  cout<<"p_BG[9]"<<p_bg[9]<<endl;
  // fLam->Draw("same");
  //  fSig->Draw("same");

   
  TCanvas* c1=new TCanvas("c1","MM ACC");
  c1->cd();
  hmm->Draw();
  hacc->SetLineColor(4);
  hmm_acc->SetFillColor(4);
  hmm_acc->SetFillStyle(3001);
  hmm_acc->Scale(1.8/40.);

  hmm_acc->Draw("same");



  TCanvas* c2=new TCanvas("c2","Hist");
  c2->Divide(1,2);
  c2->cd(1);
 
  hpeak_Al->Draw();
  fAl_cp->SetLineColor(1);
  fAl_cp->Draw("same");
  fAl->SetLineColor(2);
  fAl->Draw("same");
  c2->cd(2);
  hmm_pi->Draw();
  fpi_cp->SetLineColor(1);
  fpi->SetLineColor(3);
  fpi_cp->Draw("same");
  fpi->Draw("same");
  
}
