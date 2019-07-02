
void hrs_mmass_fit6(){
  

  TFile *fin=new TFile("./../../rootfiles/Lambda_small1_mm_ana1230.root");
  //  TFile *fin=new TFile("./../../rootfiles/Lambda_small_mm_ana1230.root");
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
 T->SetBranchStatus("ct_acc",1);
 T->SetBranchAddress("ct_acc",&ct_acc);
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


 bool vz_cut;
 bool ac_cut;
 bool ct_cut;
 bool acc_cut;

 double ac1_th=1.4;
 double ac2_th_b=3.6;
 double ac2_th_t=10.;

 /*
 ac1_th1=1.0;
 ac2_th_t=10.;
 ac2_th_b=3.5;
 */


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
  TH1F* hmm_ac2=new TH1F("hmm_ac2","MM AC1 &AC2 cut",bin_mm,min_mm,max_mm);
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
      ac_cut=true;
      if(-1.0<ct && ct<1.0)ct_cut=true;
      if((-55<ct && ct<-15)||(5<ct && ct<65))acc_cut=true;
      if(acc_cut && ac_cut && vz_cut){
	double ct2=ct;
      while(1){                                                                               
	if(-15.0<ct2 && ct2<5.0){                                                             
	  hmm_ac_acc->Fill(mm);                                               
	  
	  break;}                                                                        
	else if(ct2<-15.0){ct2=ct2+20.0;}                                                   
	else if(5.0<ct2){ct2=ct2-20.;}                                                      
      }                             
      }
                    //======== FIll Hist ==================//
      if(vz_cut && ac_cut){
	hcoin_ac->Fill(ct);
	hcoin_acc->Fill(ct_acc);}
      if(vz_cut&&ac_cut&&ct_cut){
	hmm_ac->Fill(mm);
	hmm_ac2->Fill(mm);}



 }


 hmm_ac_acc->Scale(2.0/100.);
 
 // hcoin_p->Add(hcoin_ac,hcoin_acc,1.,-6.0/96.);
 hmm_ac_p->Add(hmm_ac,hmm_ac_acc,1.,-1.0);
 //hmm_ac->Add(hmm_ac_acc,-2.0/40.);
   double def_sig,def_mean_p,def_sig_pi,def_mean_pi;
   double def_sig_k,def_mean_k;
   double def_num_k,def_num_pi,def_num_p;
 def_num_p=1.11197e+02;
 def_sig_p=7.27909e-01; 
 def_mean_p=-8.21;
 def_num_pi=6.03914e+01;
 def_mean_pi=2.98314e+00; 
 def_sig_pi=2.57384e-01;
 def_num_k=2.65916e+01;
 def_sig_k=3.74537e-01; 
 def_mean_k=0.0;
 //def_acc=22.7;


   TF1*fk=new TF1("fk","gausn(0)+pol0(3)",min_coin,max_coin);
   TF1*fp=new TF1("fp","gausn(0)+pol0(3)",min_coin,max_coin);
   TF1*fpi=new TF1("fpi","gausn(0)+pol0(3)",min_coin,max_coin);
   TF1*fcoin_acc=new TF1("fcoin_acc","pol0(0)",min_coin,max_coin);
   fk->SetNpx(2000);
   fpi->SetNpx(2000);
   fp->SetNpx(2000);   
   hcoin_ac->Fit("fcoin_acc","Rq","",-20,-15);
   double def_acc=fcoin_acc->GetParameter(0);

   fpi->FixParameter(3,def_acc);
   fp->FixParameter(3,def_acc);
   fk->FixParameter(3,def_acc);




 fpi ->SetParameter(0,def_num_pi);
 fpi ->SetParLimits(0,0.9*def_num_pi,1.5*def_num_pi);
 fpi ->SetParameter(1,def_mean_pi);
 fpi ->SetParLimits(1,def_mean_pi-0.5*def_sig_pi,def_mean_pi+0.5*def_sig_pi);
 fpi ->SetParameter(2,def_sig_pi);
 fpi ->SetParLimits(2,0.8*def_sig_pi,1.2*def_sig_pi);
 fpi ->SetParameter(3,def_acc);

 fp->SetParameter(0,def_num_p);
 fp->SetParameter(1,def_mean_p);
 fp->SetParLimits(1,def_mean_p-0.5*def_sig_p,def_mean_p+0.5*def_sig_p);
 fp->SetParameter(2,def_sig_p);
 fp->SetParLimits(2,0.9*def_sig_p,1.1*def_sig_p);
 fp->SetParameter(3,def_acc);

 fk->SetParameter(0,def_num_k);
 fk->SetParameter(1,def_mean_k);
 fk->SetParLimits(1,def_mean_k-0.5*def_sig_k,def_mean_k+0.5*def_sig_k);
 fk->SetParameter(2,def_sig_k);
 fk->SetParLimits(2,0.8*def_sig_k,1.2*def_sig_k);
 fk->FixParameter(3,def_acc);

 hcoin_ac->Fit("fk","Rq","",def_mean_k-3*def_sig_k,def_mean_k+3*def_sig_k);
 hcoin_ac->Fit("fpi ","Rq","",def_mean_pi-3*def_sig_pi,def_mean_pi+3*def_sig_pi);
 hcoin_ac->Fit("fp","Rq","",def_mean_p-3*def_sig_p,def_mean_p+3*def_sig_p);

 def_num_k=fk->GetParameter(0);
 def_mean_k=fk->GetParameter(1);
 def_sig_k=fk->GetParameter(2);

 def_num_p=fp->GetParameter(0);
 def_mean_p=fp->GetParameter(1);
 def_sig_p=fp->GetParameter(2);

 def_num_pi=fpi ->GetParameter(0);
 def_mean_pi=fpi ->GetParameter(1);
 def_sig_pi=fpi ->GetParameter(2);



 cout<<"===== coincidence Fit infomation ===== "<<endl;
 cout<<"pi number: "<<def_num_pi<<endl;
 cout<<"pi mean: "<<def_mean_pi<<endl;
 cout<<"pi sigma: "<<def_sig_pi<<endl;

 cout<<"k number: "<<def_num_k<<endl;
 cout<<"k mean: "<<def_mean_k<<endl;
 cout<<"k sigma: "<<def_sig_k<<endl;

 cout<<"p number: "<<def_num_p<<endl;
 cout<<"p mean: "<<def_mean_p<<endl;
 cout<<"p sigma: "<<def_sig_p<<endl;

 cout<<"========= coincidence =============== "<<endl;


   

  TF1* fLam=new TF1("fLam","gausn(0)+pol1(3)",1.0,1.3);

  //  TF1* fbg=new TF1("fbg","pol2",1.0,1.3);  
  TF1* fbg=new TF1("fbg","pol1",1.0,1.3);  
  hmm_ac_p->Fit("fbg","","",1.05,1.1);

  double bg[3];
  bg[0]=fbg->GetParameter(0);
  bg[1]=fbg->GetParameter(1);
  //  bg[2]=fbg->GetParameter(2);

  double p_L[3];
  p_L[0]=1.95096;
  p_L[1]=1.11555;
  p_L[2]=4.21431e-03;
  fLam->SetParameter(0,p_L[0]);
  fLam->SetParameter(1,p_L[1]);
  fLam->SetParameter(2,p_L[2]);
  fLam->FixParameter(3,bg[0]);
  fLam->FixParameter(4,bg[1]);
  //  fLam->FixParameter(5,bg[2]);
  //  fLam->SetParLimits(0,0.8*p_L[0],1.5*p_L[0]);
  fLam->SetParLimits(1,0.8*p_L[1],1.2*p_L[1]);
  fLam->SetParLimits(2,0.8*p_L[2],1.2*p_L[2]);
  hmm_ac_p->Fit("fLam","","",1.1,1.2);

  p_L[0]=fLam->GetParameter(0);
  p_L[1]=fLam->GetParameter(1);
  p_L[2]=fLam->GetParameter(2);

  double NL_err=fLam->GetParError(0);


  TF1* fbgS=new TF1("fbgS","pol1(0)",1.15,1.25);  
  double bgS[3];
  hmm_ac2->Fit("fbgS","","",1.1,1.23);
  bgS[0]=fbgS->GetParameter(0);
  bgS[1]=fbgS->GetParameter(1);
  //  bgS[2]=fbgS->GetParameter(2);
  //   bgS[0]=847.014;
  //   bgS[1]=-647.312;
  //  bgS[0]=241.513;
  //  bgS[1]=-186.222;
  TF1* fSig=new TF1("fSig","gausn(0)+pol1(3)",1.15,1.25);
  double p_S[3];
    p_S[0]=2.55776e+00;
    p_S[1]=1.19955e+00;
    p_S[2]=1.06469e-02;

    p_S[0]=9.89687e-01;
    p_S[1]=1.19910e+00;
    p_S[2]=8.89361e-03;

  fSig->SetParameter(0,p_S[0]);
  fSig->SetParameter(1,p_S[1]);
  fSig->SetParameter(2,p_S[2]);
  //fSig->SetParLimits(0,0.5*p_S[2],1.5*p_S[0]);
  fSig->FixParameter(1,p_S[1]);
  // fSig->FixParameter(2,p_S[2]);
  fSig->SetParameter(3,bgS[0]);
  fSig->SetParameter(4,bgS[1]); 
  //fSig->FixParameter(5,bgS[2]);
 //  fSig->FixParameter(2,bg[2]);
  hmm_ac2->Fit("fSig","","",1.15,1.25);

  double NS_err=fSig->GetParError(0);

  p_S[0]=fSig->GetParameter(0);
  p_S[1]=fSig->GetParameter(1);
  p_S[2]=fSig->GetParameter(2);
 
    p_S[0]=2.75048e-01;
    p_S[1]=1.19930e+00;
    p_S[2]=4.34061e-03;

  double NL,NS;
  NL=hmm_ac->Integral(hmm_ac->FindBin(p_L[1]-3*p_L[2]),hmm_ac->FindBin(p_L[1]+3*p_L[2]));
  NS=hmm_ac->Integral(hmm_ac->FindBin(p_S[1]-3*p_S[2]),hmm_ac->FindBin(p_S[1]+3*p_S[2]));

  double LBG=(NL-p_L[0]/0.002);
  double LSN=p_L[0]/0.002/LBG;
  double LSN_err=(1./sqrt(p_L[0]/0.002)+1./sqrt(LBG))*LSN;
 
  double SBG=(NS-p_S[0]/0.002);
  double SSN=p_S[0]/0.002/SBG;
  double SSN_err=(1./sqrt(p_S[0]/0.002)+1./sqrt(SBG))*SSN;


  cout<<"P_L[0]:"<<p_L[0]<<endl;
  cout<<"Lambda: "<<p_L[0]/0.002<<endl;
  cout<<"Lambda Error: "<<NL_err/0.002<<endl;
  cout<<"NL"<<NL<<endl;
  cout<<"Lam S/N"<<p_L[0]/0.002/(NL-p_L[0]/0.002)<<endl;
 cout<<"Lam S/N Error"<<LSN_err<<endl;
  cout<<"P_L[1]:"<<p_L[1]<<endl;
  cout<<"P_L[2]:"<<p_L[2]<<endl;
  cout<<"P_S[0]:"<<p_S[0]<<endl;
  cout<<"Sigma: "<<p_S[0]/0.002<<endl;
  cout<<"Sigma Error: "<<NS_err/0.002<<endl;
  cout<<"NS:"<<NS<<endl;
  cout<<"Sig S/N"<<p_S[0]/0.002/(NS-p_S[0]/0.002)<<endl;
  cout<<"Sig S/N Error"<<SSN_err<<endl;
  cout<<"P_S[1]:"<<p_S[1]<<endl;
  cout<<"P_S[2]:"<<p_S[2]<<endl;
  TCanvas* c4=new TCanvas("c4","c4");
  //  c4->Divide(3,1);
  // c4->cd(1);
  // hmm_ac->Draw();
  
  // c4->cd(2);
  // hmm_ac_acc->Draw();
  c4->cd(1);
  //hcoin_p->Draw();
  hmm_ac_p->Draw();
  fLam->SetLineColor(2);
  fLam->Draw("same");
  fSig->SetLineColor(4);
  fSig->Draw("same");
    fbg->Draw("same");

 TCanvas* c5=new TCanvas("c5","c5");
 c5->cd();
 // hcoin_ac->SetLineColor(1);
 hcoin_ac->Draw();

 fk->SetLineColor(4);
 fk->SetFillColor(4);
 fk->SetFillStyle(3001);
 fp->SetLineColor(3);
 fp->SetFillColor(3);
 fp->SetFillStyle(3001);
 fpi->SetLineColor(2);
 fpi->SetFillColor(2);
 fpi->SetFillStyle(3001);

 fk->Draw("same");
 fpi->Draw("same");
 fp->Draw("same");

TCanvas* cx=new TCanvas("cx","cx");
 cx->cd();
 hmm_ac2->Draw();
 hmm_ac_acc->SetFillStyle(3002);                                  
 hmm_ac_acc->SetFillColor(4);  
 hmm_ac_acc->Draw("same");


}
