
extern double f(double *x, double *par);
extern double df(double *x, double *par);

void mmass_calc(){
  double Ee,Ee_,Ek,pe,pe_,pk,me,mk;
  double ang=13.2*3.14/180.;
  Ee=4.3185;
  pe_=2.1;
  //  pe_=2.2;  
  pk=1.823;
  me=0.51e-3;
  mk=0.493;
  Ek=sqrt(mk*mk+pk*pk);
  Ee_=sqrt(me*me+pe_*pe_);
  pe=sqrt(Ee*Ee-me*me);
  double mt=0.938;
  double mm=sqrt((Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
		 (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang)-2*pe*pk*cos(ang)+2*pe_*pk*cos(2*ang)));
  cout<<"mass (2.1 GeV): "<<mm<<" GeV/c^2"<<endl;

  pe_=2.2;    
  mm=sqrt((Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
	  (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang)-2*pe*pk*cos(ang)+2*pe_*pk*cos(2*ang)));
  cout<<"mass (2.2 GeV): "<<mm<<" GeV/c^2"<<endl;

  
  //  cout<<"P : "<< (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang)-2*pe*pk*cos(ang)+2*pe_*pk*cos(2*ang))<<" E: "<<(Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)<<endl;
  ///  cout<<"pe*pe_ : "<<2*pe*pe_*cos(ang)<<endl;

  //  double test=f(1.0);
  //  cout<<test<<endl;
  TF1* fpe=new TF1("fpe","f",2.0,2.5,1);
  fpe->SetNpx(2000);

  TF1* fdpe=new TF1("fdpe","df",2.0,2.5,1);
  fpe->SetNpx(2000);
  
  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  fdpe->Draw();
  
}


double f(double *x,double *par){

  double Ee,Ee_,Ek,pe,pe_,pk,me,mk;
  double ang=13.2*3.14/180.;
  Ee=4.3185;
  //  pe_=2.1*x[0];
  pe_=x[0];
  par[0]=1;
  //  //  pe_=2.2;
  pk=1.823;
  me=0.51e-3;
  mk=0.493;
  Ek=sqrt(mk*mk+pk*pk);
  Ee_=sqrt(me*me+pe_*pe_);
  pe=sqrt(Ee*Ee-me*me);
  double mt=0.938;


  
  double y=sqrt(fabs(
	      (Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
		 (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang)-2*pe*pk*cos(ang)+2*pe_*pk*cos(2*ang))
		     )
		);

 return y;

}

double df(double *x,double *par){

  double Ee,Ee_,Ek,pe,pe_,pk,me,mk;
  double ang=13.2*3.14/180.;
  Ee=4.3185;
  //  pe_=2.1*x[0];
  pe_=x[0];
  par[0]=1;
  //  //  pe_=2.2;
  pk=1.823;
  me=0.51e-3;
  mk=0.493;
  Ek=sqrt(mk*mk+pk*pk);
  Ee_=sqrt(me*me+pe_*pe_);
  pe=sqrt(Ee*Ee-me*me);
  double mt=0.938;


  
  double y=sqrt(fabs(
	      (Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
		 (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang)-2*pe*pk*cos(ang)+2*pe_*pk*cos(2*ang))
		     )
		);
  pe_=2.1;
  double y0=sqrt(fabs(
	      (Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
		 (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang)-2*pe*pk*cos(ang)+2*pe_*pk*cos(2*ang))
		     )
		);

  
 return y;

}
