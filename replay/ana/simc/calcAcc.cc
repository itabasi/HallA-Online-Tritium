
bool rhrs=true;
//bool rhrs=false;

void calcAcc(){


  string ifname ="../rootfiles/simc/single_test.root";
  //ifname ="../rootfiles/simc/single.root";
  //ifname ="../rootfiles/simc/test.root";
    //int Ntot = 17706;
  //    int Ntot = 14422;
  //    Ntot = 1442207;
  int Ntot=0;
  
  if(rhrs){
    ifname ="../rootfiles/simc/single_RHRS.root";
    Ntot = 1790339;
  }else{ 
    ifname ="../rootfiles/simc/single_LHRS.root";
    Ntot = 1816284;
  }
  
  TChain* T =new TChain("SNT");
  T->Add(Form("%s",ifname.c_str()));

	  

  int RevID,LevID;
  float Rp_gen,Rth_gen,Rph_gen,Rx_gen,Ry_gen,Rz_gen;
  float Lp_gen,Lth_gen,Lph_gen,Lx_gen,Ly_gen,Lz_gen;
  
  float Rp_rec,Rth_rec,Rph_rec,Rx_rec,Ry_rec,Rz_rec;
  float Lp_rec,Lth_rec,Lph_rec,Lx_rec,Ly_rec,Lz_rec;    
  float Pm,nu;
  float uq_x,uq_y,uq_z,q;
    T->SetBranchAddress("Pm",&Pm);
    T->SetBranchAddress("uq_x_gen",&uq_x);
    T->SetBranchAddress("uq_y_gen",&uq_y);
    T->SetBranchAddress("uq_z_gen",&uq_z);
    T->SetBranchAddress("q_gen",&q);
    T->SetBranchAddress("nu",&nu);
    T->SetBranchAddress("Rp_gen",&Rp_gen);
    T->SetBranchAddress("Rth_gen",&Rth_gen);
    T->SetBranchAddress("Rph_gen",&Rph_gen);
    //    T->SetBranchAddress("h_xptar",&Rth_gen);
    //    T->SetBranchAddress("h_yptar",&Rph_gen);
    T->SetBranchAddress("Lp_gen",&Lp_gen);
    T->SetBranchAddress("Lth_gen",&Lth_gen);
    T->SetBranchAddress("Lph_gen",&Lph_gen);
    //    T->SetBranchAddress("e_xptar",&Lth_gen);
    //    T->SetBranchAddress("e_yptar",&Lph_gen);

    
  
    //int Ntot = 10000;
    int ENum=T->GetEntries();
    double ratio =(double)ENum/(double)Ntot;
    double   theta_gen_accep = 6.7272080607340551e-2; // rad
    double   phi_gen_accep = 6.9319187228489870e-2; // rad

    // Accepctance calc
    
    double dOmega_gen = 2.0 * theta_gen_accep * 2.0 * phi_gen_accep;
    dOmega_gen = 4.0 * asin(sin(theta_gen_accep)*sin(phi_gen_accep));
    dOmega_gen = 2.0*3.14*(1.0- cos(theta_gen_accep));
    double dOmega = ratio * dOmega_gen;
    double scale  = dOmega_gen/(double)Ntot;
    
    double min_rp,max_rp;
    double min_lp,max_lp;
    double Rp_mean = 1.82;
    double Lp_mean = 2.1;
    double Rp_acc  = 0.045;
    double Lp_acc  = 0.045;
    min_rp = Rp_mean*(1.0-Rp_acc); max_rp =Rp_mean*(1.0 + Rp_acc);
    min_lp = Lp_mean*(1.0-Lp_acc); max_lp =Lp_mean*(1.0 + Lp_acc);
    int nbin =100;
    TH1D* hLp = new TH1D("hLp","",nbin,min_lp,max_lp);
    TH1D* hRp = new TH1D("hRp","",nbin,min_rp,max_rp);
    TH1D* hLp_gen = new TH1D("hLp_gen","",nbin,min_lp,max_lp);
    TH1D* hRp_gen = new TH1D("hRp_gen","",nbin,min_rp,max_rp);
    TH2D* hxpyp = new TH2D("hxpyp","Mom x vs y ; xptar ; yptar",nbin,-theta_gen_accep,theta_gen_accep,nbin,-phi_gen_accep,phi_gen_accep);
    TH2D* hthph = new TH2D("hthph","Theta vs phi ; theta [rad];phi[rad]",nbin,-0.01,theta_gen_accep,nbin,-3.14,3.14);    
    hLp->SetLineColor(2);
    hRp->SetLineColor(4);
    TH1D* hth =new TH1D("hth","",100,0,0.1);
    TH1D* hph =new TH1D("hph","",100,-3.14,3.14);
    TH3D* hLp3 =new TH3D("hLp3","",nbin,-max_lp,max_lp,nbin,-max_lp,max_lp,nbin,-max_lp,max_lp);
    //    TH3D* hRAcc_th =new TH3D("hRAcc_th","",nbin,-max_rp,max_rp,nbin,nbin,0.0,theta_gen_accep,nbin,0,20);
    TH2D* hRp_th =new TH2D("hRp_th","hRp_th;mom [GeV];theta [rad] ; Acceptance [msr]",nbin,min_rp,max_rp,nbin,0.0,theta_gen_accep);
    TH2D* hLp_th =new TH2D("hLp_th","hLp_th;mom [GeV];theta [rad] ; Acceptance [msr]",nbin,min_lp,max_lp,nbin,0.0,theta_gen_accep);
    TH2D* hxy = new TH2D("hxy","",nbin,-0.5,0.5,nbin,-0.5,0.5);
    double Rp_x,Rp_y,Rp_z,Lp_x,Lp_y,Lp_z;
    double theta,phi;
    for(int i=0;i<ENum;i++){
      T->GetEntry(i);
      hLp->Fill(Lp_gen/1000.);
      hRp->Fill(Rp_gen/1000.);

      Rp_z = Rp_gen/sqrt(1.0 + Rth_gen*Rth_gen + Rph_gen*Rph_gen);
      Rp_x = Rp_z * Rth_gen;
      Rp_y = Rp_z * Rph_gen;

      Lp_z = Lp_gen/sqrt(1.0 + Lth_gen*Lth_gen + Lph_gen*Lph_gen);
      Lp_x = Lp_z * Lth_gen;
      Lp_y = Lp_z * Lph_gen;      

      TVector3 L_v,R_v;
      R_v.SetXYZ(Rp_x,Rp_y,Rp_z);
      L_v.SetXYZ(Lp_x,Lp_y,Lp_z);
      if(rhrs){
      theta = R_v.Theta();
      phi = R_v.Phi();
      hxpyp->Fill(Rth_gen,Rph_gen);
      hxy->Fill(Rp_x/1000.,Rp_y/1000.);
      }else{
      theta = L_v.Theta();
      phi = L_v.Phi();
      hxpyp->Fill(Lth_gen,Lph_gen);
      hxy->Fill(Lp_x/1000.,Lp_y/1000.);
      }      

      hthph->Fill(theta,phi);
      hLp3->Fill(Lp_x/1000.,Lp_y/1000.,Lp_z/1000.);

      hRp_th->Fill(Rp_gen/1000.,theta);
      hLp_th->Fill(Lp_gen/1000.,theta);
      
      hth->Fill(theta);
      hph->Fill(phi);
    }


    hRp_th->Scale(dOmega_gen*1000/Ntot*nbin*nbin*ratio);
    hLp_th->Scale(dOmega_gen*1000/Ntot*nbin*nbin*ratio);

    TRandom random;
    int nmax =10000000;
    for(int j=0;j<nmax;j++){
      hRp_gen->Fill(random.Uniform(min_rp,max_rp));
      hLp_gen->Fill(random.Uniform(min_lp,max_lp));
    }
    hRp_gen->Scale(double(Ntot)/double(nmax));
    hLp_gen->Scale(double(Ntot)/double(nmax));
    //  hRp_gen->Scale(double(ENum)/double(nmax));
    //  hLp_gen->Scale(double(ENum)/double(nmax));
    TGraph* gRacc=new TGraph();
    TGraph* gLacc=new TGraph();
    gRacc->SetName("gRacc");
    gLacc->SetName("gLacc");
    gLacc->SetTitle("LHRS Acceptance ; Mom [GeV]; msr");
    gRacc->SetTitle("RHRS Acceptance ; Mom [GeV]; msr");
    gRacc->SetMarkerStyle(20);
    gRacc->SetMarkerColor(4);
    gLacc->SetMarkerStyle(20);
    gLacc->SetMarkerColor(2);
    double xR,xL,yR,yL,xRg,yRg,xLg,yLg;
    for(int i =0; i<nbin;i++){
      
      xR=    hRp->GetXaxis()->GetBinCenter(i);
      yR=    hRp->GetBinContent(i)* dOmega_gen*1000.;
      xL=    hLp->GetXaxis()->GetBinCenter(i);
      yL=    hLp->GetBinContent(i)* dOmega_gen*1000.;    
      xRg=   hRp_gen->GetXaxis()->GetBinCenter(i);
      yRg=   hRp_gen->GetBinContent(i);
      xLg=   hLp_gen->GetXaxis()->GetBinCenter(i);
      yLg=   hLp_gen->GetBinContent(i);    
      
      
      
      if(yRg>0 && xRg) gRacc -> SetPoint(i,xR,yR/yRg);
      else gRacc -> SetPoint(i,xR,0.0);
    //    if(xLg>0 && yLg>0) gLacc -> SetPoint(i,xL,yL/yLg);
      if(xLg>0 && yLg>0) gLacc -> SetPoint(i,xL,yL/yLg);
      else gLacc -> SetPoint(i,xL,0.0);
      //    cout<<"i "<<i<<" xL "<<xL<<" yL "<<yL<<" xLg "<<xLg<<" xRg "<<xRg<<endl;
      //    cout<<"i "<<i<<" yL/yLg "<<yL/yLg<<" yL_gen/Ntot "<<yLg/Ntot*(double)nbin<<endl;
      
  }
    
    cout<<"dOmega_gen "<<dOmega_gen*1000<<" [msr]"<<endl;
    cout<<"Ntot "<<Ntot<<" ENum "<<ENum<<endl;
    TCanvas* c0 =new TCanvas("c0","c0");
    c0->Divide(2,1);
    c0->cd(1);
    gLacc->Draw("AP");
    c0->cd(2);
    gRacc->Draw("AP");  
    gRacc->Draw("P");
    TCanvas* c1 =new TCanvas("c1","c1");
    c1->Divide(2,1);
    c1->cd(1);
    hRp->GetYaxis()->SetRangeUser(0,(double)Ntot/(double)nbin*1.1);
    hRp->Draw("");
    hRp_gen->Draw("same");   
    c1->cd(2);
    hLp->GetYaxis()->SetRangeUser(0,(double)Ntot/(double)nbin*1.1);
    hLp->Draw("");
    hLp_gen->Draw("same");

    TCanvas* c2 =new TCanvas("c2","c2");
    c2->Divide(2,1);
    c2->cd(1);
    hxpyp->Draw("colz");
    c2->cd(2);
    hthph->Draw("colz");

    TCanvas* c3 =new TCanvas("c3","c3");
    c3->Divide(2,1);
    if(rhrs){
    c3->cd(1);
    hRp_th->Draw("LEGO3");
    c3->cd(2);
    hRp_th->Draw("colz");
    }else{
    c3->cd(1);
    hLp_th->Draw("LEGO3");
    c3->cd(2);
    hLp_th->Draw("colz");
    }

    
}
