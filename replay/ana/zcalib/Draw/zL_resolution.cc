
  int  nfoil=10;

void zL_resolution(){

  char* arm;
  bool rarm=false;
  if(rarm)arm="R";
  else arm="L";
  
  //  TFile*f1=new TFile(Form("../../rootfiles/zcalib/zt_LHRS_132.root"));
  TFile*f1=new TFile(Form("../../rootfiles/zcalib/zt_LHRS_sieve.root"));  
  TTree* t1=(TTree*)f1->Get("T");
  double rvz[100],lvz[100],rvz_c[100],lvz_c[100],cer_c;
  
  if(rarm){
  t1->SetBranchAddress("R.tr.vz", rvz);
  t1->SetBranchAddress("L.tr.vz", lvz);
  t1->SetBranchAddress("R.tr.vz_tuned", rvz_c);

  }else{
  t1->SetBranchAddress("R.tr.vz", rvz);
  t1->SetBranchAddress("L.tr.vz", lvz);
  t1->SetBranchAddress("L.tr.vz_tuned", lvz_c);
  t1->SetBranchAddress("L.cer.asum_c", &cer_c);
  }



  //==== Make Hist =======//
  double min_z=-0.2;
  double max_z=0.2;
  int bin_z=400;
  TH1F* hz=new TH1F("hz",Form("%s-HRS z hist w/o matrix tuning",arm),bin_z,min_z,max_z);
  TH1F* hz_c=new TH1F("hz_c",Form("%s-HRS z hist w/ matrix tunig",arm),bin_z,min_z,max_z);
  
  bool PID=false;
  double zt,zt_c;
  int ENum=t1->GetEntries();
  cout<<"Events : "<<ENum<<endl;

  for(int i=0;i<ENum;i++){
    PID=false;
    t1->GetEntry(i);
    if(rarm){zt=rvz[0]; zt_c=rvz_c[0];}
    else {zt=lvz[0]; zt_c=lvz_c[0];}

    if(rarm==0 && cer_c>3000)PID=true;
    if(PID){
    hz->Fill(zt);
    hz_c->Fill(zt_c);
    }
  }

  

  //====== Fit function ======//

  double z_foil[10]= {-0.125, -0.100, -0.075, -0.050, -0.025,
		       0.00, 0.025, 0.05, 0.10, 0.125};
  //  double zL_foil[10]= { -0.095, -0.0765, -0.0572, -0.0383, -0.0183, 0.00053, 0.0188, 0.0385, 0.079, 0.0975}; //zt_LHRS_132.root
  double zL_foil[10]= { -0.095, -0.071, -0.0525, -0.0383, -0.0183, 0.00053, 0.0188, 0.0385, 0.0767, 0.095};//zt_LHRS_sieve.root
  //  double z_width=0.0125;
    double z_width=0.0125;  
  double z_sig=0.005;
  double zL_sig=0.007;  
  double zL_width=0.01;
  double p0_zc[nfoil],p1_zc[nfoil],p2_zc[nfoil],p0_z[nfoil],p1_z[nfoil],p2_z[nfoil];
  double p0_zc_err[nfoil],p1_zc_err[nfoil],p2_zc_err[nfoil],p0_z_err[nfoil],p1_z_err[nfoil],p2_z_err[nfoil];  
  TF1* fz[nfoil];
  TF1* fz_c[nfoil];  
  for(int i=0;i<nfoil;i++){
    //======= fz ==========//
    //    fz[i]=new TF1(Form("fz[%d]",i),"gausn(0)",zL_foil[i]-zL_width,zL_foil[i]+zL_width);
    fz[i]=new TF1(Form("fz[%d]",i),"gausn(0)",min_z,max_z);    
    fz[i]->SetNpx(2000);
    fz[i]->SetLineColor(i+1);
    fz[i]->SetFillColor(i+1);
    if(i==9){fz[i]->SetLineColor(30);  fz[i]->SetFillColor(30);}
    fz[i]->SetFillStyle(3001);        
    fz[i]->SetParameter(1,zL_foil[i]);
    fz[i]->SetParameter(2,zL_sig);
    hz->Fit(Form("fz[%d]",i),"","RQ",zL_foil[i]-zL_width,zL_foil[i]+zL_width);
    p0_z[i]=fz[i]->GetParameter(0);
    p1_z[i]=fz[i]->GetParameter(1);
    p2_z[i]=fz[i]->GetParameter(2);
    
    //======= fz_c ==========//
    //    fz_c[i]=new TF1(Form("fz_c[%d]",i),"gausn(0)",z_foil[i]-z_width,z_foil[i]+z_width);
    fz_c[i]=new TF1(Form("fz_c[%d]",i),"gausn(0)",min_z,max_z);    
    fz_c[i]->SetNpx(2000);
    fz_c[i]->SetLineColor(i+1);
    fz_c[i]->SetFillColor(i+1);
    if(i==9){fz_c[i]->SetLineColor(30);  fz_c[i]->SetFillColor(30);}
    fz_c[i]->SetFillStyle(3001);        
    fz_c[i]->SetParameter(1,z_foil[i]);
    fz_c[i]->SetParameter(2,z_sig);    
    hz_c->Fit(Form("fz_c[%d]",i),"","RQ",z_foil[i]-z_width,z_foil[i]+z_width);
    p0_zc[i]=fz_c[i]->GetParameter(0);
    p1_zc[i]=fz_c[i]->GetParameter(1);
    p2_zc[i]=fz_c[i]->GetParameter(2);
  }

  TF1* fz_all=new TF1("fz_all","gausn(0)+gausn(3)+gausn(6)+gausn(9)+gausn(12)+gausn(15)+gausn(18)+gausn(21)+gausn(24)+gausn(27)",-0.1,0.1);
  fz_all->SetLineColor(2);
  fz_all->SetNpx(2000);


  
  TF1* fzc_all=new TF1("fzc_all","gausn(0)+gausn(3)+gausn(6)+gausn(9)+gausn(12)+gausn(15)+gausn(18)+gausn(21)+gausn(24)+gausn(27)",-0.15,0.15);
  fzc_all->SetLineColor(2);
  fzc_all->SetNpx(2000);

  for(int i=0;i<nfoil;i++){
    fz_all->SetParameter(3*i,p0_z[i]);
    fz_all->SetParameter(3*i+1,p1_z[i]);
    fz_all->SetParameter(3*i+2,p2_z[i]);
    
    fzc_all->SetParameter(3*i,p0_zc[i]);
    fzc_all->SetParameter(3*i+1,p1_zc[i]);
    fzc_all->SetParameter(3*i+2,p2_zc[i]);
  }
  hz->Fit("fz_all","","RQ",min_z,max_z);  
  hz_c->Fit("fzc_all","","RQ",min_z,max_z);
  for(int i=0;i<nfoil;i++){
    p0_z[i]=fz_all->GetParameter(3*i);
    p1_z[i]=fz_all->GetParameter(3*i+1);
    p2_z[i]=fz_all->GetParameter(3*i+2);
    p0_z_err[i]=fz_all->GetParError(3*i);
    p1_z_err[i]=fz_all->GetParError(3*i+1);
    p2_z_err[i]=fz_all->GetParError(3*i+2);
    fz[i]->SetParameters(p0_z[i],p1_z[i],p2_z[i]);
    
    p0_zc[i]=fzc_all->GetParameter(3*i);
    p1_zc[i]=fzc_all->GetParameter(3*i+1);
    p2_zc[i]=fzc_all->GetParameter(3*i+2);
    p0_zc_err[i]=fzc_all->GetParError(3*i);
    p1_zc_err[i]=fzc_all->GetParError(3*i+1);
    p2_zc_err[i]=fzc_all->GetParError(3*i+2);
    fz_c[i]->SetParameters(p0_zc[i],p1_zc[i],p2_zc[i]);
  }

  

  TCanvas* c10=new TCanvas("c10","c10");
  c10->cd();
  hz->Draw();
  fz_all->Draw("same");
  for(int i=0;i<nfoil;i++)fz[i]->Draw("same");


  
  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  hz_c->Draw();
  fzc_all->Draw("same");
  for(int i=0;i<nfoil;i++)fz_c[i]->Draw("same");


  //======= w/o matrix tuning =======//
  TGraphErrors* gz_sig = new TGraphErrors();
  gz_sig->SetTitle("Multi foil Resolution w/o matrix tuning;#foil ;#sigma_{z} [mm]");
  gz_sig->SetMarkerSize(1.0);
  gz_sig->SetMarkerColor(1);
  gz_sig->SetMarkerStyle(20);    
  TGraphErrors* gz_mean= new TGraphErrors();
  gz_mean->SetTitle("Multi foil mean (Fit - real posi) w/o matrix tuning; #foil ;#mu_{z} [mm]");  
  gz_mean->SetMarkerSize(1.0);
  gz_mean->SetMarkerColor(1);
  gz_mean->SetMarkerStyle(20);    
  for(int i=0;i<nfoil;i++){
    gz_sig->SetPoint(i,i+1,p2_z[i]*1000);
    gz_sig->SetPointError(i,0.0,p2_z_err[i]*1000);
    gz_mean->SetPoint(i,i+1,p1_z[i]*1000-z_foil[i]*1000);
    gz_mean->SetPointError(i,0.0,p1_z_err[i]*1000);}  


    //======= w/o matrix tuning scaling =======//
  TGraphErrors* gz_sig_sc = new TGraphErrors();
  gz_sig_sc->SetTitle("Multi foil Resolution w/o matrix tuning (scaled );#foil ;#sigma_{z} [mm]");
  gz_sig_sc->SetMarkerSize(1.0);
  gz_sig_sc->SetMarkerColor(3);
  gz_sig_sc->SetMarkerStyle(20);
  double scale=(z_foil[nfoil-1]-z_foil[0])/(zL_foil[nfoil-1]-zL_foil[0]);
  for(int i=0;i<nfoil;i++){
    gz_sig_sc->SetPoint(i,i+1,p2_z[i]*1000*scale);
    gz_sig_sc->SetPointError(i,0.0,p2_z_err[i]*1000*scale);}




  

  
  //======= w/ matrix tuning =======//
  TGraphErrors* gz_sig_c = new TGraphErrors();
  gz_sig_c->SetTitle("Multi foil Resolution w/ matrix tuning;#foil ;#sigma_{z} [mm]");
  gz_sig_c->SetMarkerSize(1.0);
  gz_sig_c->SetMarkerColor(2);
  gz_sig_c->SetMarkerStyle(20);    
  TGraphErrors* gz_mean_c= new TGraphErrors();
  gz_mean_c->SetTitle("Multi foil mean (Fit - real posi) w/ matrix tuning ; #foil ;#mu_{z} [mm]");  
  gz_mean_c->SetMarkerSize(1.0);
  gz_mean_c->SetMarkerColor(2);
  gz_mean_c->SetMarkerStyle(20);    
  for(int i=0;i<nfoil;i++){
    gz_sig_c->SetPoint(i,i+1,p2_zc[i]*1000);
    gz_sig_c->SetPointError(i,0.0,p2_zc_err[i]*1000);
    gz_mean_c->SetPoint(i,i+1,p1_zc[i]*1000-z_foil[i]*1000);
    gz_mean_c->SetPointError(i,0.0,p1_zc_err[i]*1000);}


  
  TCanvas* c1=new TCanvas("c1","c1");  
  c1->cd();
  gz_sig->GetXaxis()->CenterTitle();
  gz_sig->GetXaxis()->SetTitleSize(0.04);    
  gz_sig->GetYaxis()->CenterTitle();
  gz_sig->GetYaxis()->SetTitleSize(0.04);      
  gz_sig->GetYaxis()->SetRangeUser(4.0,10.);
  gz_sig->Draw("AP");
  gz_sig_c->Draw("P");
  gz_sig_sc->Draw("P");
  
  TCanvas* c2=new TCanvas("c2","c2");  
  c2->cd();
  gz_mean->GetXaxis()->CenterTitle();  
  gz_mean->GetYaxis()->CenterTitle();
  gz_mean->GetXaxis()->SetTitleSize(0.04);  
  gz_mean->GetYaxis()->SetTitleSize(0.04);  
  gz_mean->GetYaxis()->SetRangeUser(-50.0,50.0);  
  gz_mean->Draw("AP");
  gz_mean_c->Draw("P");  

  //=======================//
  //===== COMMENT =========//
  //=======================//
  for(int i=0;i<nfoil;i++){
    //    cout<<"i : "<<i<<" p0 "<<p0_zc[i]<<" p1 "<<p1_zc[i]<<" p2 "<<p2_zc[i]<<endl;
    //   cout<<"error : "<<" p0 "<<p0_zc[i]<<" p1 "<<p1_zc[i]<<" p2 "<<p2_zc[i]<<endl;
  }
  
  

}
