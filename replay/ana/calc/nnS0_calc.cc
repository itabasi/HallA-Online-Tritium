

void nnS0_calc(){

  TCanvas*c0=new TCanvas("c0","c0");
  c0->cd();
  double Pe[6]={2.081,2.180,2.278};
  double nnLk[6]={1.975,1.872,1.769};
  TGraph* gnnL=new TGraph();
    for(int i=0;i<3;i++){
  gnnL->SetPoint(i,Pe[i],nnLk[i]);}
  TF1* fnnL=new TF1("fnnL","pol1(0)",1.9,2.4);
  gnnL->Fit("fnnL","","",1.9,2.4);

    
  
  double Ee,Ee_,Ek,pe,pe_,pk,me,mk;
  double ang=13.2*3.14/180.;
  Ee=4.3185;
  pe_=2.1;
  pe_=2.2;  
  pk=1.823;
  me=0.51e-3;
  mk=0.493;
  Ek=sqrt(mk*mk+pk*pk);
  Ee_=sqrt(me*me+pe_*pe_);
  pe=sqrt(Ee*Ee-me*me);
  double mt = 2.808921112; // Tritium mass


  //  TH2D* hmom=new TH2D("hmom","Momentum bite; pe' [GeV/c] ; pK [GeV/c]",1000,pe_ - 0.045*pe_,pe_ + 0.045*pe_, 1000, pk -0.045*pk, pk + 0.045*pk);

  TH2D* hmom=new TH2D("hmom","Momentum bite; pe' [GeV/c] ; pK [GeV/c]",1000, 1.9,2.5, 1000,1.5,2.1);
  
  
  double mm = sqrt((Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
		 (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang)-2*pe*pk*cos(ang)+2*pe_*pk*cos(2*ang)));
  //  double Rp=pk -0.045*pk;
  //  double Lp=pe_-0.045*pe_;
  double Rp=1.6;
  double Lp=2.0;
  while(Rp<2.0){
    while(Lp<2.4){
  double MM = sqrt((Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
		 (pe*pe+Lp*Lp+Rp*Rp-2*pe*Lp*cos(ang)-2*pe*Rp*cos(ang)+2*Lp*Rp*cos(2*ang)));
  Lp=Lp+Lp*0.001;

  if(fabs( MM - 3.05 - 0.1)<0.05)hmom->Fill(Lp,Rp);
  cout<<"MM "<<MM<<endl;
    }
    Rp=Rp+Rp*0.001;
    //Lp=pe_-0.045*pe_;
    Lp=2.0;
  }
  

  hmom->Draw();
  gnnL->Draw("P");
}
