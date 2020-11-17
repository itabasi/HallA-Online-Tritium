

void Fitting(){

  //  string rname  = "../rootfiles/mmass/ana_Lambda/mmcalib_new/nnL_small_Ole_all.root";
  string rname  = "../rootfiles/mmass/ana_Lambda/2020-09-07_4/nnL_small_Ole_all.root";

  //  
  //   string rname2 = "../rootfiles/fsi/nnL_fsi_wb_l3.root";
  //string rname2 = "../rootfiles/fsi/nnL_fsi_wb_l0.root";

  //   string rname2  = "../rootfiles/fsi/mmass/nnL_Jl0.root";
  //  string rname2  = "../rootfiles/fsi/mmass/test.root"; //test
  string rname2  = "../rootfiles/fsi/mmass/vscale/verma_mmass_test.root"; //Vparam scale
  string rnameL  = "../rootfiles/fsi/mmass/Lam_Jl0.root";
   //    string rnameL  = "../rootfiles/fsi/Lam_fsi.root";



    bool Bin5MeV=false;
    if(Bin5MeV){
      rname  = "../rootfiles/mmass/ana_Lambda/2020-09-11/nnL_small_Ole_all.root";
      rname2 = "../rootfiles/fsi/3HL_5MeVbin.root";
      rnameL  = "../rootfiles/fsi/Lam_5MeVbin.root";
      
    }


  TFile * f1 =new TFile(rname.c_str());
  TFile * f2 =new TFile(rname2.c_str());
  TFile * fL =new TFile(rnameL.c_str());  

  TH1F* hmm_fsi1  =(TH1F*)f2->Get("hmm_fsi1");
  TH1F* hmm_fsi2  =(TH1F*)f2->Get("hmm_fsi2");
  TH1F* hmm_fsi3  =(TH1F*)f2->Get("hmm_fsi3");
  
  //  TH1F* hmm_Jl_9  =(TH1F*)f2->Get("hmm_Jl_9");
  TH1F* hmm_Jl_9  =(TH1F*)f2->Get("hmm_Jl_1"); //test Verma potental by jost func 
  hmm_fsi1->SetName("hmm_fsi1");
  hmm_fsi2->SetName("hmm_fsi2");
  hmm_fsi3->SetName("hmm_fsi3");
  //  hmm_Jl_9->SetName("hmm_Jl_9");
  hmm_Jl_9->SetName("hmm_Jl_1"); //test
  TH1D* hexp      =(TH1D*)f1->Get("h_peak_nnL");
  hexp->SetName("hexp");
  TH1D* hL =(TH1D*)fL->Get("hmm_nnL");
  //  TH1D* hL =(TH1D*)f2->Get("hmm_L");
  //  hL->Scale(43.8*2.0/1000000.);
  //  hL->Scale(36.8*2.0/1000000.);
  //  hL->Scale(36.8*2.0/1000000.*5./2.);
  //  hL->Scale(80.5079/1000000.);
    hL->Scale(80.5079/1000000.);
  hexp->Add(hL,-1.);
  TH1F* hmm  =(TH1F*)f2->Get("hmm");
  hmm->SetName("hmm");
  TH1F* hmm2  =(TH1F*)f2->Get("hmm");
  hmm2->SetName("hmm2");


  //  TH1F* hmm  =(TH1F*)f2->Get("hmm_fsi1");
  //  TH1F* hmm  =(TH1F*)f2->Get("hmm_fsi2");
  //  hmm->SetName("hmm");


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
  for(int ibin = 0;ibin<nbins;ibin++){
    x0 = hmm -> GetXaxis()->GetBinCenter(ibin);
    if( x0 < 0.0 )bin_x0 = ibin;
    //    if( x0 < 120) bin_x1 = ibin;
    if( x0 < 100)bin_x1 = ibin;
    y0[ibin]  =  hexp  -> GetBinContent(ibin);
    y1[ibin]  =  hmm   -> GetBinContent(ibin);
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
    gchi->SetPoint(wi,w,chi2);
  } // for w

  cout<<"wmin "<<wmin<<" rate "<<ymax0/ymax1<<" chi2 "<<chi2_min<<endl;
  hmm->Scale(wmin);
  hmm_fsi1->Scale(wmin);
  hmm_fsi2->Scale(wmin);
  hmm_fsi3->Scale(wmin);
  hmm_Jl_9->Scale(wmin);


  double scale_fsi1 = hmm ->GetBinContent(hmm->GetMaximumBin())/hmm_fsi1->GetBinContent(hmm_fsi1->GetMaximumBin());
  hmm_fsi1->Scale(scale_fsi1);

  double scale_fsi2 = hmm ->GetBinContent(hmm->GetMaximumBin())/hmm_fsi2->GetBinContent(hmm_fsi2->GetMaximumBin());
  hmm_fsi2->Scale(scale_fsi2);

  double scale_fsi3 = hmm ->GetBinContent(hmm->GetMaximumBin())/hmm_fsi3->GetBinContent(hmm_fsi3->GetMaximumBin());
  hmm_fsi3->Scale(scale_fsi3);  
  
  double scale = hmm ->GetBinContent(hmm->GetMaximumBin())/hmm_Jl_9->GetBinContent(hmm_Jl_9->GetMaximumBin());
  hmm_Jl_9->Scale(scale);

  // Peak Significance //

  TH1D* hmm_all   =(TH1D*)f1->Get("h_mm_nnL");
  TH1D* hmm_acc   =(TH1D*)f1->Get("h_acc_nnL");
  hmm_acc->SetName("hmm_acc");
  TH1D* hmm_bg    =(TH1D*)f1->Get("h_acc_nnL");
  hmm_bg->SetName("hmm_bg");
  hmm_bg->Add(hmm2,1);
  hmm_bg->Add(hL,1);
  hmm_Jl_9->Add(hmm_acc,1);
  hmm_Jl_9->Add(hL,1);

  hmm_fsi1->Add(hmm_acc,1);
  hmm_fsi1->Add(hL,1);
  hmm_fsi2->Add(hmm_acc,1);
  hmm_fsi2->Add(hL,1);  
  hmm_fsi3->Add(hmm_acc,1);
  hmm_fsi3->Add(hL,1);
  
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
	
  TCanvas* c0 =new TCanvas("c0","c0");
  c0->cd();
  //  wmin = 1./500.;

  //  hexp->Draw();
  hmm_all->Draw();
  hmm_bg->Draw("same");
  hmm_acc->Draw("same");
  //  hmm->Draw("same");
  //hmm2->Draw("same");
  hL->Draw("same");
  //  hL->Draw("same");
  //  hmm_fsi1->Draw("same");
  //  hmm_fsi2->Draw("same");
  //  hmm_fsi3->Draw("same");

  TCanvas* c1 =new TCanvas("c1","c1");
  c1->cd();  
  //  gchi->Draw("AP");
  //  hPS->Draw();
  gPS->Draw("AP");


  

  TCanvas* c2 =new TCanvas("c2","c2");
  c2->cd();
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
  gchi->SetName("gchip");
  gchi->Write();
  hmm_fsi1->Write();
  hmm_fsi2->Write();
  hmm_fsi3->Write();
  
  //  TCanvas* c3 =new TCanvas("c3","c3");
  //  c3->cd();  
  
}
