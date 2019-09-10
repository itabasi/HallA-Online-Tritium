
void hrs_mmass_fit4(){
  

  TFile *fin=new TFile("./../../rootfiles/Lambda_small1_mm_ana1230.root");

  TH1D* hmm=(TH1D*)fin->Get("h_mm");
  TH1D* hmm_acc=(TH1D*)fin->Get("h_mm_acc");
  TH1D* hpeak=(TH1D*)fin->Get("h_peak_mm");
  TH1D* hmm_Al=(TH1D*)fin->Get("h_mm_Al");
  TH1D* hmm_Al_acc=(TH1D*)fin->Get("h_mm_Al_acc");
  TH1D* hpeak_Al=(TH1D*)fin->Get("h_peak_Al");
  TH1D* hmm_pi=(TH1D*)fin->Get("h_mm_pi");
  hmm_pi->Add(h_mm_pi,h_mm_acc,1,-4.0/40.);


  double ct,Ra1sum,Ra2sum,mm,Lz,Rz,Rpz,Lpz,ct_acc;

  TTree *T=(TTree*)fin->Get("tree");
 T->SetBranchStatus("ct",1);    
 T->SetBranchAddress("ct",&ct);
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


 bool vz_cut;
 bool ac_cut;
 bool ct_cut;
 bool acc_cut;
 double ac1_th=1.4;
 double ac2_th_b=3.6;
 double ac2_th_t=10.;

  double min_coin=-20.0;
  double max_coin=10.0;
  double bin_coin=(max_coin-min_coin)/0.056;


 bin_coin=(int)bin_coin;
 double min_mm=0.5;
 double max_mm=1.5;
 //double bin_mm=(max_mm-min_mm)/0.002;
 //bin_mm=(int)bin_mm;
 int bin_mm=500;
 TH1F* hcoin_ac=new TH1F("hcoin_ac","Coin AC1 &AC2 cut",bin_coin,min_coin,max_coin);
 TH1F* hcoin_acc=new TH1F("hcoin_acc","Coin ACC AC1 &AC2 cut",bin_coin,min_coin,max_coin);
 TH1F* hcoin_p=new TH1F("hcoin_p","Coin ACC AC1 &AC2 cut",bin_coin,min_coin,max_coin);
 TH1F* hmm_ac=new TH1F("hmm_ac","MM AC1 &AC2 cut",bin_mm,min_mm,max_mm);
 TH1F* hmm_ac_acc=new TH1F("hmm_ac_acc","MM ACC AC1 &AC2 cut",bin_mm,min_mm,max_mm);
 TH1F* hmm_ac_p=new TH1F("hmm_ac_p","MM AC1 &AC2 cut",bin_mm,min_mm,max_mm);
 double ac1_npe,ac2_npe;
 double rz,lz;
 int evnt=T->GetEntries();
 cout<<"Events is "<<evnt<<endl;
 //  evnt=10000;
 int i=0;

 for(int k=0;k<evnt;k++){
   T->GetEntry(k);

   if(k==1000000*i){
   cout<<"k "<<k<<endl;
   i=i+1;}
   vz_cut=false;
   ac_cut=false;
   ct_cut=false;
   acc_cut=false;
   ac1_npe=Ra1sum;
   ac2_npe=Ra2sum;
   rz=Rz;
   lz=Lz;

      if((rz<0.1 && -0.1<rz) && (-0.1<lz && lz<0.1))vz_cut=true;
      if(ac1_npe<ac1_th && ac2_th_b<ac2_npe && ac2_npe<ac2_th_t)ac_cut=true;
      if(-1.0<ct && ct<1.0)ct_cut=true;
      if((-40<ct && ct<-20)||(20<ct && ct<40))acc_cut=true;
   //======== FIll Hist ==================//
      if(vz_cut && ac_cut){
	hcoin_ac->Fill(ct);
	hcoin_acc->Fill(ct_acc);}
   if(vz_cut&&ac_cut&&ct_cut)hmm_ac->Fill(mm);
   if(vz_cut&&ac_cut&&acc_cut)hmm_ac_acc->Fill(mm);


}

 hcoin_p->Add(hcoin_ac,hcoin_acc,1.,-6.0/96.);
 hmm_ac_p->Add(hmm_ac,hmm_ac_acc,1.,-2.0/40.)
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
  hmm_ac->Fit("fbg","","",1.0,1.3);
  double p_bg[13]; 
  for(int i=0;i<13;i++){p_bg[i]=fbg->GetParameter(i);}


  cout<<"p_bg[9]:"<<p_bg[9]<<endl;



T

  TF1* fLam=new TF1("fLam","gausn(0)+pol2(4)",1.0,1.3);
  TF1* fSig=new TF1("fSig","gausn(0)+pol2(4)",1.0,1.3);
  
  /*for(int i=0;i<13;i++){
  fLam->FixParameter(i+3,p_bg[i]);
  fSig->FixParameter(i+3,p_bg[i]);}
  */


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
  hmm_ac_p->Fit("fLam","","",1.1,1.13);


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
 


 
   
  THStack *hbg=new THStack("hbg","MM Back Ground");  

  TF1* fpi_cp=hmm_pi->GetFunction("fpi");
  TF1* fAl_cp=hpeak_Al->GetFunction("fAl");
  //hpeak_Al->Scale(p_bg[0]);
  //hmm_pi->Scale(p_bg[9]);
  //hbg->Add(hpeak_Al);
  //  hbg->Add(hmm_pi); 
  


  /*
  TCanvas* c0=new TCanvas("c0","MM Hist");
  c0->cd();
  hpeak->Draw();
  
  fbg->SetLineColor(4);
  fbg->SetFillColor(4);
  fbg->SetFillStyle(3001);
  fAl->Draw("same");
  fbg->Draw("same");

  cout<<"p_BG[0]"<<p_bg[0]<<endl;
  cout<<"p_BG[9]"<<p_bg[9]<<endl;

  
   
  TCanvas* c1=new TCanvas("c1","MM ACC");
  c1->cd();
  hmm->Draw();
  hmm_acc->SetLineColor(4);
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
  */


  TCanvas* c3=new TCanvas("c3","AC cut Hist");
  c3->Divide(2,1);
  c3->cd(1);
  /*
  hcoin_ac->Draw();
  hcoin_acc->SetLineColor(4);
  hcoin_acc->SetFillColor(4);
  hcoin_acc->SetFillStyle(3001);
  hcoin_acc->Scale(6./96.);
  hcoin_acc->Draw("same");
  */
  hmm->Draw();
  hmm_acc->SetLineColor(4);
  hmm_acc->SetFillColor(4);
  hmm_acc->SetFillStyle(3001);
  hmm_acc->Scale(1.8/40.);
  hmm_acc->Draw("same");
  c3->cd(2);
  hmm_ac->Draw();
  hmm_ac_acc->Scale(2./40.);
  hmm_ac_acc->SetLineColor(4);
  hmm_ac_acc->SetFillColor(4);
  hmm_ac_acc->SetFillStyle(3001);
  hmm_ac_acc->Draw("same");

  TCanvas* c4=new TCanvas("c4","AC cut Hist");
  c4->Divide(2,1);
  c4->cd(1);
  //   hcoin_p->Draw();
  hmm_ac->Draw();
  c4->cd(2);
  hmm_ac_p->Draw();
  fLam->Draw("same");

}
