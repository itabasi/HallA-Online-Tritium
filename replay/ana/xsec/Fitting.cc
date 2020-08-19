

void Fitting(){

  string rname  = "../rootfiles/mmass/ana_Lambda/mmcalib_new/nnL_small_Ole_all.root";
  //string rname2 = "../rootfiles/fsi/nnL_fsi_l3.root";
  //  string rname2 = "../rootfiles/fsi/nnL_fsi_wb_l3.root";
  string rname2 = "../rootfiles/fsi/nnL_fsi_wb_l0.root";
  string rnameL  = "../rootfiles/fsi/Lam_fsi.root";



  TFile * f1 =new TFile(rname.c_str());
  TFile * f2 =new TFile(rname2.c_str());
  TFile * fL =new TFile(rnameL.c_str());  

  TH1F* hmm_fsi1  =(TH1F*)f2->Get("hmm_fsi1");
  TH1F* hmm_fsi2  =(TH1F*)f2->Get("hmm_fsi2");
  TH1F* hmm_fsi3  =(TH1F*)f2->Get("hmm_fsi3");
  TH1D* hexp =(TH1D*)f1->Get("h_peak_nnL");
  hexp->SetName("hexp");
  TH1D* hL =(TH1D*)fL->Get("hmm_nnL");
  //  TH1D* hL =(TH1D*)f2->Get("hmm_L");
  hL->Scale(43.8*2.0/1000000.);
  hexp->Add(hL,-1.);
  TH1F* hmm  =(TH1F*)f2->Get("hmm");
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
    if( x0 < 80)bin_x1 = ibin;
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

  TCanvas* c0 =new TCanvas("c0","c0");
  c0->cd();
  //  wmin = 1./500.;
  hmm->Scale(wmin);
  hmm_fsi1->Scale(wmin);
  hmm_fsi2->Scale(wmin);
  hmm_fsi3->Scale(wmin);
  hexp->Draw();
  hmm->Draw("same");
  hL->Draw("same");
  //  hL->Draw("");
  hmm_fsi1->Draw("same");
  hmm_fsi2->Draw("same");
  hmm_fsi3->Draw("same");

  TCanvas* c1 =new TCanvas("c1","c1");
  c1->cd();  
  gchi->Draw("AP");
}
