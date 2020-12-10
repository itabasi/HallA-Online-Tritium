

void FittingQF(){



    string rname  = "../rootfiles/mmass/ana_Lambda/best/nnL_small_Ole_all.ro\
ot";
    string rname2 = "../rootfiles/fsi/mmass/nnL_2MeVBin.root";
    string rnameL  = "../rootfiles/fsi/mmass/Lam_2MeVBin.root";

    TFile * f1 =new TFile(rname.c_str());
    TFile * f2 =new TFile(rname2.c_str());
    TFile * fL =new TFile(rnameL.c_str());

    
    TH1D* hexp =(TH1D*)f1->Get("h_peak_L");
    hexp->SetName("hexp");

    TH1F* hmm  =(TH1F*)f2->Get("hmm_L");
    hmm->SetName("hmm");
    TH1F* hmmL  =(TH1F*)fL->Get("hmm_L");
    hmmL->SetName("hmmL");
    TH1D* hexp2 = (TH1D*)hexp->Clone();
    hexp2->SetName("hexp2");


    int nbins = hmm->GetXaxis()->GetNbins();
    int nmin  = hmm->GetXaxis()->GetXmin();
    int nmax  = hmm->GetXaxis()->GetXmax();

    double ymax0 = hexp -> GetBinContent(hexp->GetMaximumBin());
    double ymax1 = hmm  -> GetBinContent(hmm->GetMaximumBin());

    double y0[nbins],y1[nbins];
    int bin_x0=0;
    int bin_x1=nbins;;
    double x0;
    bool max= true;
    int bin_x2 = 0;
    int bin_x3 = 0;
    for(int ibin = 0;ibin<nbins;ibin++){
      x0 = hmm -> GetXaxis()->GetBinCenter(ibin);
      if( x0 < -40.0 )bin_x0 = ibin;
      if( x0 < 50.0)bin_x1 = ibin;
      y0[ibin]  =  hexp  -> GetBinContent(ibin);
      y1[ibin]  =  hmm   -> GetBinContent(ibin);

      if(x0 < -5.0) bin_x2 = ibin;
      if(x0 <  5.0) bin_x3 = ibin;

    }

    
    TGraphErrors* gchi=new TGraphErrors();
    gchi->SetMarkerStyle(3);
    gchi->SetMarkerColor(2);
    double chi2,chi2_min,w,wmin;
    int wmax=100;
    double width = 0.01;
    chi2_min=1e20;
    for(int wi=0; wi<wmax;wi++){
      chi2=0.0;
      w = ymax0/ymax1*(1.0 +(double)(-wmax/2. + wi)*(double)width);
      for(int ibin = bin_x0;ibin<bin_x1;ibin++){
	if(bin_x2< ibin && ibin < bin_x3)continue; // Lambda peak
	if(y0[ibin]!=0)chi2 +=  pow((y0[ibin] -w*y1[ibin])/y0[ibin],2.0)/double(nbins-bin_x0);
      }

      //    cout<<"chi "<<chi2<<endl;
      if(chi2_min > chi2){
	chi2_min = chi2;
	wmin = w;
      }
      gchi->SetPoint(wi,w,chi2);
    } // for w



    hmm->Scale(wmin);
    hmm->SetLineColor(4);
    hexp->SetLineColor(1);
    hexp2->Add(hmm,-1.);
    hexp2->SetLineColor(2);

    
    double nLam = hexp2->Integral(hexp2->FindBin(-5.0),hexp2->FindBin(5.0));
    cout<<"wmin "<<wmin<<" rate "<<ymax0/ymax1<<" chi2 "<<chi2_min<<" nLam "	<<nLam<<endl;



    int binL = hmmL->GetXaxis()->GetNbins();
    double min_Lam = hmmL-> GetXaxis()->GetXmin();
    double max_Lam = hmmL->GetXaxis()->GetXmax();
    double binSize = (max_Lam - min_Lam)/(double)binL;
    double nL = hmmL->GetEntries();
    hmmL->Scale(nLam/nL*binSize);
    //        hmmL->Scale(nLam/nL);

    TCanvas* c0 =new TCanvas("c0","c0");
    c0->cd();
    hexp->Draw();
    //    hexp2->Draw("same");
    hmm->Draw("same");
    hmmL->Draw("same");
    TCanvas* c1 =new TCanvas("c1","c1");
    c1->cd();
    gchi->Draw("AP");
    TCanvas* c2 =new TCanvas("c2","c2");
    c2->cd();
    hexp->Draw();


  ///////////////////////////////////
  /////// Fitting nnL QF ////////////
  //////////////////////////////////



  TH1F* hmm_fsi1  =(TH1F*)f2->Get("hmm_fsi1");
  TH1F* hmm_fsi2  =(TH1F*)f2->Get("hmm_fsi2");
  TH1F* hmm_fsi3  =(TH1F*)f2->Get("hmm_fsi3");
  
  //  TH1F* hmm_Jl_9  =(TH1F*)f2->Get("hmm_Jl_9");
  TH1F* hmm_Jl_9  =(TH1F*)f2->Get("hmm_Jl_1"); //test Verma potental by jost func 
  hmm_fsi1->SetName("hmm_fsi1");
  hmm_fsi2->SetName("hmm_fsi2");
  hmm_fsi3->SetName("hmm_fsi3");
  hmm_Jl_9->SetName("hmm_Jl_9");
  //hmm_Jl_9->SetName("hmm_Jl_1"); //test
  
  TH1D* hexp_nnL      =(TH1D*)f1->Get("h_peak_nnL");
  hexp_nnL->SetName("hexp_nnL");
  TH1D* hL_3H =(TH1D*)fL->Get("hmm_nnL");

  //    hL_3H->Scale(80.5079/1000000.);
  double Nev = (double)hL_3H->GetEntries();
  double scale_Lam = nLam/Nev;
  
  hL_3H->Scale(scale_Lam);
  
  hexp_nnL->Add(hL_3H,-1.);
  //  TH1F* hmm  =(TH1F*)f2->Get("hmm");
  //  hmm->SetName("hmm");
  TH1F* hmm2  =(TH1F*)f2->Get("hmm");
  hmm2->SetName("hmm2");


  //  TH1F* hmm  =(TH1F*)f2->Get("hmm_fsi1");
  //  TH1F* hmm  =(TH1F*)f2->Get("hmm_fsi2");
  //  hmm->SetName("hmm");

  /*
    int nbins = hmm2->GetXaxis()->GetNbins();
    int nmin  = hmm2->GetXaxis()->GetXmin();
    int nmax  = hmm2->GetXaxis()->GetXmax();

    int nbins = hmm2->GetXaxis()->GetNbins();
    int nmin  = hmm2->GetXaxis()->GetXmin();
    int nmax  = hmm2->GetXaxis()->GetXmax();
    double y0[nbins],y1[nbins];
    int bin_x0=0; int bin_x1=nbins;

  */

  nbins = hmm2->GetXaxis()->GetNbins();
  nmin  = hmm2->GetXaxis()->GetXmin();
  nmax  = hmm2->GetXaxis()->GetXmax();  
  ymax0 = hexp_nnL -> GetBinContent(hexp_nnL->GetMaximumBin());
  ymax1 = hmm2  -> GetBinContent(hmm2->GetMaximumBin());

  bin_x0=0;
  bin_x1=nbins;
  max= true;

  
  for(int ibin = 0;ibin<nbins;ibin++){ x0 = hmm2
      -> GetXaxis()->GetBinCenter(ibin); if( x0 < 0.0 )bin_x0 = ibin; //
  if( x0 < 120) bin_x1 = ibin; if( x0 < 100)bin_x1 = ibin; y0[ibin] =
  hexp_nnL -> GetBinContent(ibin); y1[ibin] = hmm2 ->
  GetBinContent(ibin); }

  
  TGraphErrors* gchi2=new TGraphErrors();
  gchi2->SetMarkerStyle(3);
  gchi2->SetMarkerColor(2);
  //  double chi2,chi2_min,w,wmin;
  //  int wmax=100;
  wmax=100;
  //  double width = 0.01;
  width = 0.01;
  chi2_min=1e20;
  for(int wi=0; wi<wmax;wi++){
    chi2=0.0;
    //    w = ymax0/ymax1*(1.0 +(double)(-wmax/2. + wi)*width);
    //    w = ymax0/ymax1*(1.0 +(double)(-wmax/2. + wi)/(double)width);
    w = ymax0/ymax1*(1.0 +(double)(-wmax/2. + wi)*(double)width);
    for(int ibin = bin_x0;ibin<bin_x1;ibin++){
      if(y0[ibin]!=0)chi2 +=  pow((y0[ibin] -w*y1[ibin])/y0[ibin],2.0)/double(nbins-bin_x0);

    }
    //    cout<<"chi "<<chi2<<endl;
    if(chi2_min > chi2){
      chi2_min = chi2;
      wmin = w;
    }
    gchi2->SetPoint(wi,w,chi2);
  } // for w

  cout<<"wmin "<<wmin<<" rate "<<ymax0/ymax1<<" chi2 "<<chi2_min<<endl;
  hmm2->Scale(wmin);
  hmm_fsi1->Scale(wmin);
  hmm_fsi2->Scale(wmin);
  hmm_fsi3->Scale(wmin);
  hmm_Jl_9->Scale(wmin);


  double scale_fsi1 = hmm2 ->GetBinContent(hmm2->GetMaximumBin())/hmm_fsi1->GetBinContent(hmm_fsi1->GetMaximumBin());
  hmm_fsi1->Scale(scale_fsi1);

  double scale_fsi2 = hmm2 ->GetBinContent(hmm2->GetMaximumBin())/hmm_fsi2->GetBinContent(hmm_fsi2->GetMaximumBin());
  hmm_fsi2->Scale(scale_fsi2);

  double scale_fsi3 = hmm2 ->GetBinContent(hmm2->GetMaximumBin())/hmm_fsi3->GetBinContent(hmm_fsi3->GetMaximumBin());
  hmm_fsi3->Scale(scale_fsi3);  
  
  double scale = hmm2 ->GetBinContent(hmm2->GetMaximumBin())/hmm_Jl_9->GetBinContent(hmm_Jl_9->GetMaximumBin());
  hmm_Jl_9->Scale(scale);

  // Peak Significance //

  TH1D* hmm_all   =(TH1D*)f1->Get("h_mm_nnL");
  TH1D* hmm_acc   =(TH1D*)f1->Get("h_acc_nnL");
  hmm_acc->SetName("hmm_acc");
  TH1D* hmm_bg    =(TH1D*)f1->Get("h_acc_nnL");
  hmm_bg->SetName("hmm_bg");
  hmm_bg->Add(hmm2,1);
  hmm_bg->Add(hL_3H,1);
  hmm_Jl_9->Add(hmm_acc,1);
  hmm_Jl_9->Add(hL_3H,1);

  hmm_fsi1->Add(hmm_acc,1);
  hmm_fsi1->Add(hL_3H,1);
  hmm_fsi2->Add(hmm_acc,1);
  hmm_fsi2->Add(hL_3H,1);  
  hmm_fsi3->Add(hmm_acc,1);
  hmm_fsi3->Add(hL_3H,1);
  
  TH1F* hPS =new TH1F("hPS","Peak Significance ; -B_{#Lambda} [MeV] ; Peak Significance",nbins,nmin,nmax);
  TGraphErrors* gPS=new TGraphErrors();
  gPS->SetName("gPS");
  gPS->SetMarkerStyle(22);
  gPS->SetMarkerColor(2);
  
  double tot[nbins],qf[nbins],ps[nbins];
  double ds,dsn,dps,dn,s,n;
  for(int ibin = 0;ibin<nbins;ibin++){

    tot[ibin]   =  hmm_all   -> GetBinContent(ibin);
    qf[ibin]    =  hmm_bg   -> GetBinContent(ibin);
    double x    =  hmm_all   -> GetBinCenter(ibin);
    ps[ibin] = sqrt(pow(tot[ibin]-qf[ibin],2)/(tot[ibin]));

    if(tot[ibin]>1.)hPS->Fill(x,ps[ibin]);
    else hPS->Fill(x,0.0);

    //    ds  = sqrt( fabs(tot[ibin] - qf[ibin] ));
    //    dsn = sqrt(tot[ibin]);
    //    dps = 1./sqrt(tot[ibin])*( ds + 1./2.*fabs(tot[ibin]-qf[ibin])/sqrt(tot[ibin]) );
    s = (tot[ibin] - qf[ibin]);
    n = tot[ibin];
    ds = fabs(1.0 -1./2.*s/(s+n))*sqrt(s/(s+n));
    dn = fabs(1./2.*s/(s+n))*sqrt(n/(s+n));
    dps = ds +dn;
    if(tot[ibin]>1. && s>0){
      gPS->SetPoint(ibin,x,ps[ibin]);
      gPS->SetPointError(ibin,0,dps);
    }else             gPS->SetPoint(ibin,-100,0.0);
    
  }
	
  TCanvas* c10 =new TCanvas("c10","c10");
  c10->cd();
  //  wmin = 1./500.;

  //  hexp_nnL->Draw();
  hmm_all->Draw();
  hmm_bg->Draw("same");
  hmm_acc->Draw("same");
  //  hmm->Draw("same");
  //hmm2->Draw("same");
  hL_3H->Draw("same");
  //  hL_3H->Draw("same");
  //  hmm_fsi1->Draw("same");
  //  hmm_fsi2->Draw("same");
  //  hmm_fsi3->Draw("same");

  TCanvas* c11 =new TCanvas("c11","c11");
  c11->cd();  
  //  gchi2->Draw("AP");
  //  hPS->Draw();
  gPS->Draw("AP");


  

  TCanvas* c12 =new TCanvas("c12","c12");
  c12->cd();
  hmm_all->SetLineColor(1);

  hmm_acc->SetLineColor(4);
  hmm_acc->SetFillColor(4);
  hmm_acc->SetFillStyle(3002);
  hmm_bg->SetLineColor(8);
  hmm_bg->SetFillColor(8);
  hmm_bg->SetFillStyle(3002);
  hmm_Jl_9->SetLineColor(14);
  hmm_Jl_9->SetFillColor(14);
  hmm_Jl_9->SetFillStyle(3002);

  hmm_fsi1->SetLineColor(1);
  hmm_fsi1->SetFillColor(1);
  hmm_fsi1->SetFillStyle(3002);

  hmm_fsi2->SetLineColor(2);
  hmm_fsi2->SetFillColor(2);
  hmm_fsi2->SetFillStyle(3002);    

  hmm_fsi3->SetLineColor(4);
  hmm_fsi3->SetFillColor(4);
  hmm_fsi3->SetFillStyle(3002);    
  
  hmm_all->Draw();
  hmm_all->SetTitle("Missing mass; -B_{#Lambda} [MeV] ; Counts/1 MeV");
  hmm_all->GetXaxis()->SetTitleSize(0.05);
  hmm_all->GetYaxis()->SetTitleSize(0.05);
  hmm_bg->Draw("same");
  hmm_acc->Draw("same");
  hmm_Jl_9->Draw("same");
  hmm_fsi1->Draw("same");
  hmm_fsi2->Draw("same");
  hmm_fsi3->Draw("same");
  
  TFile * ofr =new TFile("./test.root","recreate");

  hmm_all->Write();
  hmm_acc->Write();
  hmm_bg->Write();
  hmm->Write();
  gchi2->SetName("gchi2");
  gchi2->Write();
  hmm_fsi1->Write();
  hmm_fsi2->Write();
  hmm_fsi3->Write();
  
  //  TCanvas* c3 =new TCanvas("c3","c3");
  //  c3->cd();  
  
}
