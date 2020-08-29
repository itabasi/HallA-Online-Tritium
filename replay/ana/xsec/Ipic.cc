#include <complex.h>

void Ipic(){


  TGraphErrors* g0=new TGraphErrors();
  TGraphErrors* g1=new TGraphErrors();
  TGraphErrors* g0R=new TGraphErrors();
  TGraphErrors* g1R=new TGraphErrors();  
  double a0= -2.29;
  double r0 = 3.15;

  double alpha =  (1.+sqrt(1. - 2.*r0/a0) / r0);  
  double beta  = -2.0/r0 +alpha;

  double delta;
  complex<double> Jl;

  double I0,I1;
  double k,p;
  double hbarc=197.;
  complex<double>xi = complex<double>(0.0,1.0);
  for(int i=0; i<100;i++){
  k = ((double)i*10.)/hbarc;
  p = k*hbarc;
  Jl = (k - xi*beta)/(k- xi*alpha);
  I0 = 1./ real( Jl*conj(Jl)  );
  g0->SetPoint(i,p,I0);

  delta = acos(k/(-1./a0 +1./2.*r0*k*k) );
  I1    = pow(sin(delta + r0*k)/sin(r0*k),2.);
  //  I1    = pow(sin(delta + r0*k*hbarc)/sin(r0*k*hbarc),2.);
  //  I1 /= 3.14*4.0;
  if(sin(r0*k)<1.0e-10)I1=0.0;
  g1->SetPoint(i,p,I1);


  }



  for(int i=0;i<100;i++){

    double r =(int)i*0.05;
    k =10./hbarc;
    
    alpha =  (1.+sqrt(1. - 2.*r/a0) / r);  
    beta  = -2.0/r +alpha;

    //test 
    alpha =  (1.+sqrt(1. - 2.*r0/a0) / r0);  
    beta  = -2.0/r0 +alpha;
    
    Jl = (k - xi*beta)/(k- xi*alpha);
    I0 = 1./ real( Jl*conj(Jl)  );
    g0R->SetPoint(i,r,I0);
    
    delta = acos(k/(-1./a0 +1./2.*r0*k*k) );
    I1    = pow(sin(delta + r*k)/sin(r*k),2.);
    //  I1    = pow(sin(delta + r0*k*hbarc)/sin(r0*k*hbarc),2.);
    //  I1 /= 3.14*4.0;
    if(sin(r*k)<1.0e-10)I1=0.0;
    g1R->SetPoint(i,r,I1);    
    
  }

  
  TCanvas*c0=new TCanvas("c0","c0");
  c0->cd();
  g0->SetMarkerColor(2);
  g0->SetMarkerStyle(7);
  g1->SetMarkerColor(4);
  g1->SetMarkerStyle(20);  

  g1->Draw("AP");
  g0->Draw("P");


  
  TCanvas*c1=new TCanvas("c1","c1");
  c1->cd();
  g0R->SetMarkerColor(2);
  g0R->SetMarkerStyle(20);
  g1R->SetMarkerColor(4);
  g1R->SetMarkerStyle(20);  

  g1R->Draw("AP");
  g0R->Draw("P");  
  
}

  

