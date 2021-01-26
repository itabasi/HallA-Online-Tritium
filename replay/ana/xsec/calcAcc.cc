
bool Tkine =true;
//bool Tkine =false;
//bool rhrs = true;
bool rhrs = false;
bool coll = false;
bool test =false;
void calcAcc(){


  string ifname ="../rootfiles/simc/single_test.root";
  //ifname ="../rootfiles/simc/single.root";
  //ifname ="../rootfiles/simc/test.root";
  //  ifname ="../rootfiles/simc/single_gen.root";

    //int Ntot = 17706;
  //    int Ntot = 14422;
  //    Ntot = 1442207;
  int Ntot=0;
  string ofpname; // output file  
  if(rhrs){
    ifname ="../rootfiles/simc/acceptance/single_RHRS_new.root";
    ofpname="./param/accept/RHRS_accept.param"; // output file
    //    Ntot = 1790339;
    Ntot = 2550000;
    
    if(coll){
      ifname ="../rootfiles/simc/acceptance/single_RHRS_coll.root";
      ofpname="./param/accept/RHRS_accept_wcoll.param"; // output file
      Ntot = 3005116;
    }    
  }else{ 
    ifname ="../rootfiles/simc/acceptance/single_LHRS.root";
    ofpname="./param/accept/LHRS_accept.param"; // output file
    //    Ntot = 1816284;
    Ntot = 2580001;

    if(coll){
      ifname ="../rootfiles/simc/acceptance/single_LHRS_coll.root";
      ofpname="./param/accept/LHRS_accept_wcoll.param"; // output file
      Ntot = 3066520;
    }
  }


  if(Tkine && !rhrs){
    ifname ="../rootfiles/simc/acceptance/single_LHRS_Tkine.root";
    ofpname="./param/accept/LHRS_Tkine_accept.param"; // output file
    Ntot = 2616929;
  }

  
  string  ifname_gen ="../rootfiles/simc/acceptance/single_gen.root";

  
  if(test)  ofpname="./param/test.param"; // output file test mode


 
  TChain* T =new TChain("SNT");
  T->Add(Form("%s",ifname.c_str()));
  TChain* Tg =new TChain("SNT");
  Tg->Add(Form("%s",ifname_gen.c_str()));


  // parameter file //
  //  string ofpname="./param/accept/RHRS_accept.param"; // output file

  ofstream ofp(ofpname.c_str());

  ofp <<"######## Acceptance parameters ##########################"<<endl;
  ofp <<"## Auther K. Itabashi "<<endl;
  ofp <<"## Acceptance(mom) =Sum_theta dOmega(mom,theta) "<<endl;
  ofp<<"## mom [GeV]  # theta [rad]  # dOmega(mom,theta) "<<endl;
  //
  //

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

    float Lp_gen_g,Lth_gen_g,Lph_gen_g;
    Tg->SetBranchAddress("Lp_gen",&Lp_gen_g);
    Tg->SetBranchAddress("Lth_gen",&Lth_gen_g);
    Tg->SetBranchAddress("Lph_gen",&Lph_gen_g);    
    //    cout<<"ENum Tg "<<Tg->GetEntries()<<endl;
  
    //int Ntot = 10000;
    int ENum=T->GetEntries();
    double ratio =(double)ENum/(double)Ntot;
    double   theta_gen_accep = 6.7272080607340551e-2; // rad
    double   phi_gen_accep = 6.9319187228489870e-2; // rad

    // Accepctance calc
    
    double dOmega_gen = 2.0 * theta_gen_accep * 2.0 * phi_gen_accep;
    //    dOmega_gen = 4.0 * asin(sin(theta_gen_accep)*sin(phi_gen_accep));
    dOmega_gen = 2.0*3.14*(1.0- cos(theta_gen_accep));
    double dOmega = ratio * dOmega_gen;
    double scale  = dOmega_gen/(double)Ntot;
    
    double min_rp,max_rp;
    double min_lp,max_lp;
    double Rp_mean = 1.82;
    double Lp_mean = 2.1;
    if(Tkine) Lp_mean = 2.21807;
    double Rp_acc  = 0.045;
    double Lp_acc  = 0.045;
    min_rp = Rp_mean*(1.0-Rp_acc); max_rp =Rp_mean*(1.0 + Rp_acc);
    min_lp = Lp_mean*(1.0-Lp_acc); max_lp =Lp_mean*(1.0 + Lp_acc);
    int nbin =100;
    TH1D* hLp = new TH1D("hLp","",nbin,min_lp,max_lp);

    TH1D* hRp = new TH1D("hRp","",nbin,min_rp,max_rp);
    TH1D* hRp2 = new TH1D("hRp2","",nbin,min_rp,max_rp);
    TH1D* hLp_gen = new TH1D("hLp_gen","",nbin,min_lp,max_lp);
    TH1D* hRp_gen = new TH1D("hRp_gen","",nbin,min_rp,max_rp);
    TH2D* hxpyp = new TH2D("hxpyp","Mom x vs y ; xptar ; yptar",nbin,-theta_gen_accep,theta_gen_accep,nbin,-phi_gen_accep,phi_gen_accep);
    TH2D* hthph = new TH2D("hthph","Theta vs phi ; theta [rad];phi[rad]",nbin,-0.01,theta_gen_accep,nbin,-3.14,3.14);    
    hLp->SetLineColor(2);
    hRp->SetLineColor(4);
    TH1D* hth =new TH1D("hth","",100,0,0.1);
    TH1D* hth2 =new TH1D("hth2","",100,0,0.1);
    TH1D* hth_acc =new TH1D("hth_acc","hth_acc;theta [rad];N/Ngen*2PI",100,0,0.1);
    TH1D* hthL =new TH1D("hthL","hthL",100,0,0.1);
    TH1D* hthL_gen =new TH1D("hthL_gen","hthL_gen",100,0,0.1);
    TH1D* hthR_gen =new TH1D("hthR_gen","hthR_gen",100,0,0.1);
    TH1D* hthR =new TH1D("hthR","hthR",100,0,0.1);
    TH1D* hph =new TH1D("hph","",100,-3.14,3.14);
    TH3D* hLp3 =new TH3D("hLp3","",nbin,-max_lp,max_lp,nbin,-max_lp,max_lp,nbin,-max_lp,max_lp);
    
    //    TH3D* hRAcc_th =new TH3D("hRAcc_th","",nbin,-max_rp,max_rp,nbin,nbin,0.0,theta_gen_accep,nbin,0,20);
    TH2D* hRp_th =new TH2D("hRp_th","hRp_th;mom [GeV];theta [rad] ; Acceptance [msr]",nbin,min_rp,max_rp,nbin,0.0,theta_gen_accep);
    TH2D* hRp_th_gen =new TH2D("hRp_th_gen","hRp_th_gen;mom [GeV];theta [rad] ; Acceptance [msr]",nbin,min_rp,max_rp,nbin,0.0,theta_gen_accep);
    TH2D* hRp_th_acc =new TH2D("hRp_th_acc","hRp_th_acc;mom [GeV];theta [rad] ; Acceptance/theta [msr/rad]",nbin,min_rp,max_rp,nbin,0.0,theta_gen_accep);
    TH2D* hLp_th =new TH2D("hLp_th","hLp_th;mom [GeV];theta [rad] ; Acceptance [msr]",nbin,min_lp,max_lp,nbin,0.0,theta_gen_accep);
    TH2D* hLp_th_gen =new TH2D("hLp_th_gen","hLp_th_gen;mom [GeV];theta [rad] ; Acceptance [msr]",nbin,min_lp,max_lp,nbin,0.0,theta_gen_accep);
    TH2D* hLp_th_acc =new TH2D("hLp_th_acc","hLp_th_acc;mom [GeV];theta [rad] ; Acceptance/theta [msr/rad]",nbin,min_lp,max_lp,nbin,0.0,theta_gen_accep);
    TH2D* hxy = new TH2D("hxy","",nbin,-0.5,0.5,nbin,-0.5,0.5);

    

    
    double Rp_x,Rp_y,Rp_z,Lp_x,Lp_y,Lp_z;
    double theta,phi;
    double Rp_xg,Rp_yg,Rp_zg,Lp_xg,Lp_yg,Lp_zg;
    double theta_g,phi_g;
    cout<<" Fill Start "<<endl;
    for(int i=0;i<ENum;i++){
      T->GetEntry(i);
      //      Tg->GetEntry(i);
      hLp->Fill(Lp_gen/1000.);
      hRp->Fill(Rp_gen/1000.);

      Rp_z = Rp_gen/sqrt(1.0 + Rth_gen*Rth_gen + Rph_gen*Rph_gen);
      Rp_x = Rp_z * Rth_gen;
      Rp_y = Rp_z * Rph_gen;

      Lp_z = Lp_gen/sqrt(1.0 + Lth_gen*Lth_gen + Lph_gen*Lph_gen);
      Lp_x = Lp_z * Lth_gen;
      Lp_y = Lp_z * Lph_gen;      


      Lp_zg = Lp_gen_g/sqrt(1.0 + Lth_gen_g*Lth_gen_g + Lph_gen_g*Lph_gen_g);
      Lp_xg = Lp_zg * Lth_gen_g;
      Lp_yg = Lp_zg * Lph_gen_g;      
      
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
      //      hthL->Fill(theta);
      }

      TVector3 L_vg;
      L_vg.SetXYZ(Lp_xg,Lp_yg,Lp_zg);
      theta_g = L_vg.Theta();
      phi_g = L_vg.Phi();
      hthR->Fill(R_v.Theta());
      hthL->Fill(L_v.Theta());
      
      hthph->Fill(theta,phi);
      hLp3->Fill(Lp_x/1000.,Lp_y/1000.,Lp_z/1000.);

      hRp_th->Fill(Rp_gen/1000.,theta);
      hLp_th->Fill(Lp_gen/1000.,theta);
      
      hth->Fill(theta);
      hph->Fill(phi);
      //      cout<<"theta_g "<<theta_g<<" Lp "<<Lp_gen_g<<endl;
      //hLp_th_gen->Fill(Lp_gen_g/1000.,theta_g);
      //      hth2->Fill(theta_g);
    }

    //    hLp->GetXaxis()->FindBin(Lp_mean);
    //    hLp_th_gen->Scale(double(Ntot)/double(ENum));
    //    hth2->Scale(double(Ntot)/double(ENum));
    //    hRp_th->Scale(dOmega_gen*1000/Ntot*nbin*nbin*ratio);
    //    hLp_th->Scale(dOmega_gen*1000/Ntot*nbin*nbin*ratio);

    TRandom random;
    int nmax =3000000;
    double mom_R_gen,mom_L_gen,theta_gen;




    
    for(int k=0;k<nmax;k++){
      mom_R_gen = random.Uniform(min_rp,max_rp);
      mom_L_gen = random.Uniform(min_lp,max_lp);
      theta_gen=-10.0;
      while(theta_gen>theta_gen_accep || theta_gen<0){
	theta_gen = acos(random.Uniform(-1.0,1.0)); }
      hRp_gen->Fill(mom_R_gen);
      hLp_gen->Fill(mom_L_gen);
      hRp_th_gen->Fill(mom_R_gen,theta_gen);
      hLp_th_gen->Fill(mom_L_gen,theta_gen);
      hthL_gen->Fill(theta_gen);
      hthR_gen->Fill(theta_gen);
    }
    
    
    hthL_gen->Scale(double(Ntot)/double(nmax));
    hthR_gen->Scale(double(Ntot)/double(nmax));
    hRp_gen->Scale(double(Ntot)/double(nmax));
    hLp_gen->Scale(double(Ntot)/double(nmax));
    hRp_th_gen->Scale(double(Ntot)/double(nmax));
    hLp_th_gen->Scale(double(Ntot)/double(nmax));


    double x,y,z;
    double x_gen,y_gen,z_gen;
    double sum=0.0;
    double width_y =theta_gen_accep/(double)nbin;
    double accept =0.0;
    
    // x : momentum
    // y : theta
    
    for(int xbin=0;xbin<nbin;xbin++){
      accept =0.0;
      for(int ybin=0;ybin<nbin;ybin++){
	if(rhrs){
	  x = hRp_th_acc->GetXaxis()->GetBinCenter(xbin);
	  y = hRp_th_acc->GetYaxis()->GetBinCenter(ybin);
	  //	  x_gen = hRp_th_gen->GetXaxis()->GetBinCenter(xbin);
	  //	  y_gen = hRp_th_gen->GetYaxis()->GetBinCenter(ybin);
	  z = hRp_th->GetBinContent(hRp_th->GetBin(xbin,ybin));
	  z_gen = hRp_th_gen->GetBinContent(hRp_th_gen->GetBin(xbin,ybin));
	}else{
	  x = hLp_th_acc->GetXaxis()->GetBinCenter(xbin);
	  y = hLp_th_acc->GetYaxis()->GetBinCenter(ybin);
	  //	  x_gen = hLp_th_gen->GetXaxis()->GetBinCenter(xbin);
	  //	  y_gen = hLp_th_gen->GetYaxis()->GetBinCenter(ybin);
	  z = hLp_th->GetBinContent(hLp_th->GetBin(xbin,ybin));
	  z_gen = hLp_th_gen->GetBinContent(hLp_th_gen->GetBin(xbin,ybin));	  
	}
	

	//	if(z_gen>0 && z>0 && z_gen>z)z = z/z_gen*2.0*acos(-1.0);
	//       	else if( z > z_gen && z_gen>0 && z>0)z=  2.0*acos(-1.0);
	//	else z =0.0;


	if(z_gen>0 && z>0 && z_gen>z)z = z/z_gen*2.0*acos(-1.0)*width_y*sin(y);
	else if( z > z_gen && z_gen>0 && z>0)z =  2.0*acos(-1.0)*width_y*sin(y);
	else z =0.0;	

	accept +=z;
	
	if(rhrs){
	  hRp_th_acc->Fill(x,y,z);
	}else{
	  hLp_th_acc->Fill(x,y,z);
	}


	if(z>1.)cout<<"x "<<x<<" y "<<y<<" z "<<z<<" z_gen "<<z_gen<<endl;

	ofp << x <<" " << y <<" "<< z <<endl;	
	if(xbin==nbin/2){
	  hth2->Fill(y,z);
	  sum +=z*sin(y)*width_y;
	}

      }
      
    }

    
    cout<<" integral theta "<<sum<<endl;
    double th1,th2,th0,xx;
    for(int th=0;th<nbin;th++){
      if(rhrs){
      th1 =  hthR_gen->GetBinContent(hthR_gen->GetBin(th));
      th2 =  hthR->GetBinContent(hthL->GetBin(th));
      }else{
      th1 =  hthL_gen->GetBinContent(hthL_gen->GetBin(th));
      th2 =  hthL->GetBinContent(hthL->GetBin(th));
      }
      xx  = hth_acc->GetBinCenter(th);
      if(th2>0 && th1>0)th0=th2/th1;
      else th0=0.0;
      hth_acc->Fill(xx,th0*2.0*acos(-1.0));
    }





    
    

    
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
      if(xLg>0 && yLg>0)gLacc -> SetPoint(i,xL,yL/yLg);
      else gLacc -> SetPoint(i,xL,0.0);
      //      if(xLg>0 && yLg>0) ofp << xL<<" " << yL/yLg/1000. <<endl;

      
      
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
    
    TCanvas* c1 =new TCanvas("c1","c1");
    c1->Divide(2,1);
    c1->cd(1);
    hLp->GetYaxis()->SetRangeUser(0,(double)Ntot/(double)nbin*1.1);
    hLp->Draw("");
    hLp_gen->Draw("same");
    hRp->GetYaxis()->SetRangeUser(0,(double)Ntot/(double)nbin*1.1);
    c1->cd(2);
    hRp->Draw("");
    hRp_gen->Draw("same");   
    TCanvas* c2 =new TCanvas("c2","c2");
    c2->Divide(2,1);
    c2->cd(1);
    hthL_gen->Draw();
    hthL->Draw("same");
    c2->cd(2);
    hthR_gen->Draw();
    hthR->Draw("same");


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

    TCanvas* c4 =new TCanvas("c4","c4");
    c4->Divide(2,1);
    c4->cd(1);
    hth_acc->Draw();
    //    hRp_th->Draw("LEGO3");
    //    c4->cd(2);
    //    hRp_th_gen->Draw("LEGO3");
    //    c4->cd(3);
    //    hLp_th->SetFillColor(4);
    //    hLp_th->Draw("LEGO3");
    c4->cd(2);
    hth2->Draw();
    //    hLp_th_gen->SetFillColor(6);
    //    hLp_th_gen->Draw("LEGO3");

    TCanvas* c5 =new TCanvas("c5","c5");
    c5->Divide(2,1);
    hLp_th_acc->SetFillColor(6);
    hRp_th_acc->SetFillColor(3);
    c5->cd(1);
    if(rhrs)    hRp_th_acc->Draw("SURF2");
    else        hLp_th_acc->Draw("SURF2");
    c5->cd(2);
    if(rhrs) hRp_th_acc->Draw("colz");
    else     hLp_th_acc->Draw("colz");
    
    //    hLp_th_acc->GetYaxis()->SetRangeUser(0.002,theta_gen_accep);
    //    hLp_th_acc->GetZaxis()->SetRangeUser(0.0,20);




    //    TCanvas*c6 =new TCanvas("c6","c6");
    //    c6->Divide(2,1);
    //    c6->cd(1);
    //    hthL_gen->Draw();
    //    c6->cd(2);
    //    hRp2->Draw();
}
