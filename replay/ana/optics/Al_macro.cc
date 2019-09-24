
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
  x[0]=x[0]-par[4];
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     //     sum += fland * TMath::Exp(-(xx-par[4])/par[1]);
     sum += fland * TMath::Exp(-(xx)/par[1]);
     //     sum += fland * TMath::Exp((xx-par[4])/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     //     sum += fland * TMath::Exp(-(xx-par[4])/par[1]);
     sum += fland * TMath::Exp(-(xx)/par[1]);
     //     sum += fland * TMath::Exp((xx-par[4])/par[1]);
  }
  //  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  return val;
}


void Al_macro(){

  //    string fname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/optics/Al_target/tritium_Al_R.root";
  //     string fname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/optics/Al_target/Lambda_small_test2_rasR.root";
  string fname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/optics/Al_targetx/tritium_Al_Rac.root";
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
  TF1* fexp=new TF1("fexp","exp(- (x-[1])/[0] )",-0.2,0.2);
  //  TF1* fpol4=new TF1("fpol4","fR_Al",-0.1,0.1,1);
  //  TF1* f=new TF1("f","[0]*(152.557 + x*-4.11848 + x*x*-1067.08 + pow(x,3.0)*53118.2 +pow(x,4.0)*2.58378e+06);",-0.2,0.2);
  TF1* fAl=new TF1("fAl","expgaus2",-0.2,0.2,5);
  double foil1 =-0.125;
  double foil2 = 0.125;
  double width = 0.0125;
  double p1[10],p2[10],p3[10];
  double exp_p[2];

  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  fAl ->SetParameters(5,0.01,0.003,-0.125,-2*0.125-0.005);
  fAl->SetParLimits(0,5.0,100);
  hz->Fit("fAl","","",foil1-width,foil1+width);


  exp_p[0]=fAl->GetParameter(1);
  exp_p[1]=fAl->GetParameter(4);
  fexp->SetParameters(exp_p[0],exp_p[1]);
  fAl->SetLineColor(2);
  hz->Draw();
  fAl->Draw("same");
  //   fexp->Draw("same");
  //   fexp->Draw(); 

  TCanvas* c1=new TCanvas("c1","c1");
  c1->cd();
  //  fAl->Draw();
  //  fexp->Draw("same");
 fexp->Draw(); 
}


