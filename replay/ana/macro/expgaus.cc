//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTF1(TF1 *f, int LColor=2, int LWidth=0, int LStyle=1, int Npx=1000){
  f->SetLineColor(LColor);
  f->SetLineWidth(LWidth);
  f->SetLineStyle(LStyle);
  f->SetNpx(Npx);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus(double *x, double *par) {
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

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus_gaus(double *x, double *par){
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  return expgaus(x,&par[0]) + par[4]*TMath::Gaus(x[0],par[3],par[2],1);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expexpgaus_gaus(double *x, double *par){
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  return expgaus(x,&par[0]) + expgaus(x,&par[4]) + par[8]*TMath::Gaus(x[0],par[10],par[9],1);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
