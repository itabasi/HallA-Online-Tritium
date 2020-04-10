
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus1(double x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x,par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x,par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
  //  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus2(double *x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
  //  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  return val;
}

double expgaus3(double *x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
  //  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  val = val + par[4]*x[0]+ par[5];
  return val;
}

double expgaus_mean(double *x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
  //  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  val = val*x[0];
  return val;
}


double expgaus_sigma(double *x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
  //  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  val = val*(x[0] - par[4])*(x[0] - par[4]);
  return val;
}




void Sig_tail(){

  //  string ifname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/mmass/ana_Lambda/mmcalib_new/Lambda_small_OleH_all.root";
  string ifname ="../rootfiles/simc/1Hsimc.root";
  TFile* ofr = new TFile(ifname.c_str());
  //  TH1F* hSig=(TH1F*)ofr->Get("h_peak_L");
  TH1F* hSig=(TH1F*)ofr->Get("hS");
  TF1* fSig=new TF1("fSig","expgaus2",-0.300,0.300,4);
  TF1* fpol1=new TF1("fpol1","[0]*x+[1]",-0.300,0.300);
  TF1* ftot=new TF1("ftot","expgaus3",-0.300,0.300,6);
  fSig->SetNpx(2000);
  fSig->SetLineColor(2);
  fSig->SetParameters(300,0.002,0.0015,0.001);
  ftot->SetNpx(2000);
  ftot->SetLineColor(4);
  ftot->SetParameters(300,0.002,0.0015,0.001);


  hSig->Fit("fSig","","",-0.01,0.03);
  hSig->Fit("ftot","","",-0.01,0.03);

  //  fSig->SetParameter(0,1);
  TF1* fmu=new TF1("fmu","expgaus_mean",-0.300,0.300,4);
  TF1* fsigma=new TF1("fsigma","expgaus_sigma",-0.300,0.300,5);

  double pL[4];
  pL[0]=fSig->GetParameter(0);
  pL[1]=fSig->GetParameter(1);
  pL[2]=fSig->GetParameter(2);
  pL[3]=fSig->GetParameter(3);

  fmu->SetParameters(1.0,pL[1],pL[2],pL[3]);
  double mean=fmu->Integral(-0.05,0.1);

  //  double mean=fSig->Integral(-0.05,0.1);
  cout<<" mean "<<mean<<endl;

  fsigma->SetParameters(1.0,pL[1],pL[2],pL[3],mean);
  double sigma=fsigma->Integral(-0.05,0.1);
  cout<<" sigma "<<sqrt(sigma)<<endl;



  TF1* fmax=new TF1("fmax","expgaus2",-0.300,0.300,4);
  fmax->SetParameters(1.0,pL[1],pL[2],pL[3]);
  double max = fmax->GetMaximum(); 
  cout<<"max "<<max<<endl;
  double zero= fmax->Eval(0);
  cout<<"zero "<<zero<<endl;
  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  //  hSig->GetXaxis()->SetRangeUser(-20,100);
  hSig->Draw();
  fSig->Draw("same");
  ftot->Draw("same");

  double par[4]={1.0,pL[1],pL[2],pL[3]};
  double x[1]={pL[3]};
  double test =expgaus2(x,par);

  cout<<"max "<<expgaus2(x,par);

}


