//===================================//
// Aurther Itabashi 2019 June 30th
// nnL experiment mass resolution 
double M_HYP,Mt,M_nnL;
double Ee_,Ek,pe,pe_,pk, th_ee_, th_ek,th_e_k;
double Ee    = 4318.0; // [MeV]
double me = 0.511; // [MeV/c^2]
double mk = 493.7; // [MeV/c^2]




double Theta_ek(double Pe_, double theta_ee_){
  
  double E_e_ = sqrt(Pe_*Pe_ + me *me);
  double Q=sqrt(4*Ee*E_e_*pow(sin(theta_ee_/2),2));
  double Pv=sqrt(pow(Ee-E_e_,2)+pow(Q,2));
  double Ev=Ee - E_e_;
  double Av=acos((Ee-E_e_*(1-pow(Q,2)/(2*Ee*E_e_)))/Pv);//[rad]
  double Avk =0.0;
  return Av + Avk;
}
double Pk(double Pe_, double theta_ee_){
  
  double E_e_ = sqrt(Pe_*Pe_ + me *me);
  double Q=sqrt(4*Ee*E_e_*pow(sin(theta_ee_/2),2));
  double Pv=sqrt(pow(Ee-E_e_,2)+pow(Q,2));
  double Ev=Ee - E_e_;
  double Av=acos((Ee-E_e_*(1-pow(Q,2)/(2*Ee*E_e_)))/Pv);//[rad]
  double Avk =0.0;
  double A,B,C,D;//A,B,C,D is constance
  
  A=pow(M_HYP,2)+pow(Pv,2)-pow(Ev+Mt,2)-pow(mk,2);
  B=(pow(A,2)-4*pow(Ev+Mt,2)*pow(mk,2))/4;
  C=pow(Pv*cos(Avk),2)-pow(Ev+Mt,2);
  D=-A*Pv*cos(Avk);
  double  Pk=(-D-sqrt(pow(D,2)-4*C*B))/(2*C);
  return Pk;
}

void calc_reso(){



  //bool nnL=false;
  bool nnL=true;
  //==== Hydrogen run ======//
  Mt    = 938.27; // [MeV/c^2] proton mass
  M_HYP = 1115.68; //[MeV/c^2] Lambda mass
  pe_   = 2100.0; // [MeV/c] 
  if(nnL){
  Mt    = 2808.921112;
  M_HYP = 2998.;
  pe_   = 2200.0; // [MeV/c]
  }

  // Al //
  //  Mt = 25.126507e3;
  //  M_HYP = 25.312188e3;
  
  double sig_th_ee_, sig_th_ek,sig_th_e_k;  
  double sig_pe_, sig_pk, sig_pe;

  pe = sqrt(Ee*Ee - me*me);
  pk = 1800.;
  th_ee_= 13.2*3.14/180.;
  th_ek = 13.2*3.14/180.;
  th_e_k = th_ee_ + th_ek;
  sig_th_ee_ =  7.0e-3/1.55;
  sig_th_ek  =  7.0e-3/1.55;
  //  sig_th_ee_ =  7.0e-3/2.35;
  //  sig_th_ek  =  7.0e-3/2.35;
  //  sig_th_e_k =  7.0e-3;
  sig_pe  = pe  * 1.2e-4/2.35;
  sig_pe_ = pe_ * 2.5e-4/2.35;
  sig_pk  = pk  * 2.5e-4/2.35;
  
  //==== Gausian function =====//
  TF1* fAee_ = new TF1("fAee_","gausn(0)",-5.0*sig_th_ee_,5.0*sig_th_ee_);
  fAee_  ->SetParameters(1,0.0,sig_th_ee_);
  TF1* fAek = new TF1("fAek","gausn(0)",-5.0*sig_th_ek,5.0*sig_th_ek);
  fAek  ->SetParameters(1,0.0,sig_th_ek);
  TF1* fAe_k = new TF1("fAe_k","gausn(0)",-5.0*sig_th_e_k,5.0*sig_th_e_k);
  fAe_k  ->SetParameters(1,0.0,sig_th_e_k);
  TF1* fpe = new TF1("fpe","gausn(0)",-5.0*sig_pe,5.0*sig_pe);
  fpe  ->SetParameters(1,0.0,sig_pe);
  TF1* fpe_ = new TF1("fpe_","gausn(0)",-5.0*sig_pe_,5.0*sig_pe_);
  fpe_ ->SetParameters(1,0.0,sig_pe_);
  TF1* fpk = new TF1("fpk","gausn(0)",-5.0*sig_pk,5.0*sig_pk);
  fpk  ->SetParameters(1,0.0,sig_pk);

  TF1* fth_ee_ = new TF1("fth_ee_","gausn(0)",-1,1);
  //  fth_ee_  ->SetParameters(1,0.0,0.17);
  fth_ee_  ->SetParameters(1,0.0,1.);
    //====== Fill hist ========//

  TH1D* hmm = new TH1D("hmm","",250,-10,10);
  TH1D* hAee_ = new TH1D("hAee_","",100,10.,15.); 
  TH1D* hpe_ = new TH1D("hpe_","",1000,2000,2400);

  TF1* fmm =new TF1("fmm","gausn(0)",-10,10);
  fmm->SetNpx(2000);
  fmm->SetLineColor(2);
  int nmax=10000;
  double mm ,mm_b;
  double Pe_, PK , Pe, Th_ee_, Th_ek, Th_e_k;
  double dpe_, dpk, dpe, dth_ee_, dth_ek, dth_e_k;
  double Eeb,Ee_b,Ekb;
  double th_ee_cent=th_ee_;
  for (int n =0;n<nmax;n++){

    Pe_  = 0.0;
    dpe_ = 0.0;
    mm   = 0.0;
    mm_b = 0.0;
    th_ee_ = 0.0;

    //    th_ee_ = fth_ee_ ->GetRandom();
    th_ee_ = th_ee_cent + th_ee_;
    
    dpe = fpe ->GetRandom();        Pe = pe + dpe;
    dpe_ = fpe_ ->GetRandom();      Pe_ = pe_ + dpe_;
    dth_ee_ = fAee_ ->GetRandom();  Th_ee_ = th_ee_ + dth_ee_;

    
    Ee  = sqrt(me*me + Pe *Pe );
    Ee_ = sqrt(me*me + Pe_*Pe_);
    Eeb  = sqrt(me*me + pe *pe );
    Ee_b = sqrt(me*me + pe_*pe_);    

    //    pk = Pk(Pe_ , Th_ee_);
    //    th_ek = Theta_ek(Pe_ , Th_ee_);    
    pk = Pk(pe_ , th_ee_);
    th_ek = Theta_ek(pe_ , th_ee_);    
    //    cout<<"th_ek "<<th_ek<<" pe_ "<<pe_<<" th_ee_ "<<th_ee_<<endl;
    dpk = fpk ->GetRandom();      PK = pk + dpk;
    dth_ek = fAek ->GetRandom();  Th_ek = th_ek + dth_ek;
    Ek = sqrt(mk*mk + PK*PK);
    Ekb = sqrt(mk*mk + pk*pk);
    Th_e_k = Th_ek + Th_ee_;
    th_e_k = th_ek + th_ee_;
    mm = sqrt((Ee + Mt - Ee_ - Ek)*(Ee + Mt - Ee_ - Ek)-
    (Pe*Pe+Pe_*Pe_+PK*PK-2*Pe*Pe_*cos(Th_ee_)-2*Pe*PK*cos(Th_ek)+2*Pe_*PK*cos(Th_e_k)));

    mm_b = sqrt((Eeb + Mt - Ee_b - Ekb)*(Eeb + Mt - Ee_b - Ekb)-
    (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(th_ee_)-2*pe*pk*cos(th_ek)+2*pe_*pk*cos(th_e_k)));
    
    hmm->Fill( mm - M_HYP );
    //hmm->Fill( mm_b - M_HYP );
    //    cout<<"mm_b "<<mm_b<<" th_ee_ "<<th_ee_<<" th_ek "<<th_ek<<" th_e_k "<<th_e_k<<endl;
    hAee_ ->Fill( (Th_ee_ )*180./3.14);
    hpe_ ->Fill(Pe_);
  }
 
  TCanvas* c0 =new TCanvas("c0","c0");
  c0->cd();
  hmm->Draw();
  hmm->Fit("fmm","","",-10,10);
  fmm->Draw("same");


  cout<<"=========================="<<endl;
  cout<<"====== COMMENT ==========="<<endl;
  cout<<"=========================="<<endl;
  cout<<"nnL flag "<<nnL<<endl;
  cout<<"mm sigma "<<fmm->GetParameter(2)<<endl;
  cout<<"mm FWHM "<<2.35*(fmm->GetParameter(2))<<" [MeV]"<<endl;
  
}
