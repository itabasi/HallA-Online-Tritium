

void Fitting_Lam(){

  string rname  = "../rootfiles/mmass/ana_Lambda/mmcalib_new/nnL_small_Ole_all.root";
  //string rname2 = "../rootfiles/fsi/nnL_fsi_l3.root";
  string rname2 = "../rootfiles/fsi/nnL_fsi_wb_l3.root";

  string rnameL  = "../rootfiles/fsi/Lam_fsi.root";



  TFile * f1 =new TFile(rname.c_str());
  TFile * f2 =new TFile(rname2.c_str());
  TFile * fL =new TFile(rnameL.c_str());  

  TH1F* hmm_fsi1  =(TH1F*)f2->Get("hmm_fsi1");
  TH1F* hmm_fsi2  =(TH1F*)f2->Get("hmm_fsi2");
  TH1F* hmm_fsi3  =(TH1F*)f2->Get("hmm_fsi3");
  TH1D* hexp =(TH1D*)f1->Get("h_peak_L");
  hexp->SetName("hexp");
  //  TH1D* hL =(TH1D*)fL->Get("hmm_L");
  //  hexp->Add(hL,-80./1000000.);
  TH1F* hmm  =(TH1F*)f2->Get("hmm_L");
  //  TH1F* hmm  =(TH1F*)f2->Get("hmm_fsi1");
  TH1F* hmmL  =(TH1F*)fL->Get("hmm_L");
  //  hmm->SetName("hmm");

  TH1D* hexp2 = (TH1D*)hexp->Clone();
  //  TH1D* hmmL = (TH1D*)hmm->Clone();
  hexp2->SetName("hexp2");
  hmmL->SetName("hmmL");
  
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
    //    w = ymax0/ymax1*(1.0 +(double)(-wmax/2. + wi)*width);
    //    w = ymax0/ymax1*(1.0 +(double)(-wmax/2. + wi)/(double)width);
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
  cout<<"wmin "<<wmin<<" rate "<<ymax0/ymax1<<" chi2 "<<chi2_min<<" nLam "<<nLam<<endl;  

  int binL = hmmL->GetXaxis()->GetNbins();
  //  cout<<" binL "<<binL<<" bin "<<nbins<<endl;
  //  double nL = hmmL->Integral(hmmL->FindBin(-5.0),hmmL->FindBin(5.0))*2.0;
  

  double nL = hmmL->GetEntries();  
  hmmL->Scale(nLam/nL*2.0);
  
  TCanvas* c0 =new TCanvas("c0","c0");
  c0->cd();
  //  wmin = 1./500.;

  
  hexp->Draw();
  hexp2->Draw("same");
  hmm->Draw("same");
  hmmL->Draw("same");
  //  hmm_fsi1->Draw("same");
  //  hmm_fsi2->Draw("same");
  //  hmm_fsi3->Draw("same");

  TCanvas* c1 =new TCanvas("c1","c1");
  c1->cd();  
  gchi->Draw("AP");
}
