
double Eloss(double *x ,double* par){

  double phi=x[0]/180.*3.14;
  //  double ang=phi;
  double ang = + phi - 13.2/180.*3.14; //right arm
  double f= -1.31749*sin(-4.61513*ang) +2.0368;

  return f;
};


double Eloss_L(double *x ,double* par){

  double phi=x[0]/180.*3.14;

  double ang = - phi - 13.2/180.*3.14; //right arm
  double f= -1.3576*sin(-4.5957*ang) +2.0909;

  return f;
};


void Eloss_calc(){


  //  TF1* fEloss =new TF1("fEloss","Eloss",-0.5,0,0);
  TF1* fEloss =new TF1("fEloss","Eloss",-30,30,0);
  TF1* fEloss_L =new TF1("fEloss_L","Eloss_L",-30,30,0);
  fEloss_L->SetLineColor(2);
  TCanvas* c0=new TCanvas("c0","c0");
  //  fEloss->Draw();
  //  fEloss_L->Draw("same");
  fEloss_L->Draw("");
}
