
double fR_Al(double *x, double *par){

  //  double a=par[0]/1.08628e+02;
  //par[0]=4.37570e-01;
  //   double a=par[0]/1.55801e+02;
  double a=par[0]/9965.;
  double y=a*(152.557 + x[0]*-4.11848 + x[0]*x[0]*-1067.08 + pow(x[0],3.0)*53118.2 +pow(x[0],4.0)*2.58378e+06);
  return y;

}


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
  x[0]=x[0]+par[5];
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-(xx-par[4])/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-(xx-par[4])/par[1]);
  }
  //  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  return val;
}


void Al_macro(){

  //    string fname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/optics/Al_target/tritium_Al_R.root";
  //     string fname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/optics/Al_target/Lambda_small_test2_rasR.root";
  string fname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/optics/Al_target/tritium_Al_Rac.root";
     string fname2="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/optics/Al_target/tritium_Al_R.root";
  //     string fname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/optics/Al_target/tritium_111326_R.root";
  //string fname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/momcalib/Lambda_small_H1_0911_init.root";
  TFile* ofs=new TFile(fname.c_str());
  TFile* ofs2=new TFile(fname2.c_str());
  double par[10];
  //  TH1F* hz=(TH1F*)ofs->Get("hz_rc");
  TH1F* hz=(TH1F*)ofs->Get("hz_ac");
  TF1* f0=new TF1("f0","expgaus2",-0.2,0.2,6);
  TF1* f1=new TF1("f1","gausn(0)",-0.2,0.2);
  TF1* f2=new TF1("f2","gausn(0)",-0.2,0.2);
  TF1* f3=new TF1("f3","[0]*x+[1]",-0.2,0.2);
  TF1* f4=new TF1("f4","expo(0)",-0.2,0.2);
  TF1* f5=new TF1("f5","expo(0)",-0.2,0.2);
  TF1* fpol4=new TF1("fpol4","fR_Al",-0.1,0.1,1);
  TF1* f=new TF1("f","[0]*(152.557 + x*-4.11848 + x*x*-1067.08 + pow(x,3.0)*53118.2 +pow(x,4.0)*2.58378e+06);",-0.2,0.2);
  double foil1 =-0.125;
  double foil2 = 0.125;
  double width = 0.0125;
  double p1[10],p2[10],p3[10];
  
  hz->Fit("f1","RQ","RQ",foil1-width,foil1+width);
  p1[0]= f1->GetParameter(0);
  p1[1]= f1->GetParameter(1);
  p1[2]= f1->GetParameter(2);
  f0->SetParameters(100,0.005,p1[2],foil1+0.041,0.0,0.213);
  //  f0->SetParLimits(0,100,1000);
  //  f0->FixParameter(1,0.01);
  //  f0->FixParameter(2,0.005   );
  hz->Fit("f0","RQ","RQ",-0.15,-0.1);

  //  hz->Fit("fpol4","","",-0.105,0.105);


  hz->Fit("f2","","",foil2-width,foil2+width);
  p2[0]=  f2->GetParameter(0);
  p2[1]=  f2->GetParameter(1);
  p2[2]=  f2->GetParameter(2);
  hz->Fit("f3","RQ","RQ",-0.2,-0.18);
  p3[0]=  f1->GetParameter(0);
  p3[1]=  f1->GetParameter(1);

  double max_peak=hz->GetBinContent(hz->GetMaximumBin());
  cout<<"max_peak "<<max_peak<<" p2[0] "<<p2[0]<<endl;
  //  fpol4->SetParameter(0,p2[0]);
  fpol4->SetParameter(0,max_peak);
  
  TF1*fAl =new TF1("fAl","gausn(0)+gausn(3)+[6]*x+[7]",-0.2,0.2);
  double pAl[10];
  fAl->SetParameters(p1[0],p1[1],p1[2],p2[0],p2[1],p2[2],p3[0],p3[1]);
  fAl->FixParameter(0,p1[0]);
  fAl->FixParameter(1,p1[1]);
  fAl->FixParameter(2,p1[2]);
  fAl->FixParameter(3,p2[0]);
  fAl->FixParameter(4,p2[1]);
  fAl->FixParameter(5,p2[2]);

  //  hz->Fit("fAl","","RQ",-0.2,0.2);

  
  for(int i=0;i<7;i++)pAl[i]= fAl->GetParameter(i);

  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  hz->Draw();
  fAl->SetNpx(2000);
  fAl->SetLineColor(4);
  fAl->SetFillColor(4);
  fAl->SetFillStyle(3002);
  fAl->Draw("same");
  f1->SetLineColor(3);
  f1->Draw("same");
  fpol4->SetLineColor(4);
  fpol4->Draw("same");

    TCanvas* c1=new TCanvas("c1","c1");
    c1->cd();
    hz->Draw();
    fpol4->Draw("same");
    //    f0->Draw("same");
    double num_Al=fpol4->Integral(-0.1,0.1);
    num_Al=num_Al/(0.4/1000.);
    cout<<"Al events "<<num_Al<<endl;
    cout<<"Al all events "<<hz->GetEntries()<<endl;
    cout<<"Al ratio "<<num_Al/hz->GetEntries()<<endl;
}


