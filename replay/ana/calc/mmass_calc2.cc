
extern double f(double *x, double *par);
extern double df(double *x, double *par);

void mmass_calc2(){
  double Ee,Ee_,Ek,pe,pe_,pk,me,mk;
  double ang=13.2*3.14/180.;
  Ee=4.3185;
  // pe_=2.1;
  pe_=2.218;  
  pk=1.823;
  me=0.51e-3;
  mk=0.493;
  Ek=sqrt(mk*mk+pk*pk);
  Ee_=sqrt(me*me+pe_*pe_);
  pe=sqrt(Ee*Ee-me*me);
  //  double mt=0.938;
  //  double mt  = 2.9948138266; // [GeV] nnL threshold mass
  //  double mt2 = 25.3091;      // [GeV] 27MgL threshold mass
  double mt  = 2.808921112; // [GeV] Tritium mass
  double mt2 = 25.13314157;      // [GeV] 27Al threshold mass



  double mm=sqrt((Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
		 (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang)-2*pe*pk*cos(ang)+2*pe_*pk*cos(2*ang)));

  double mm2=sqrt((Ee+mt2-Ee_-Ek)*(Ee+mt2-Ee_-Ek)-
		 (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang)-2*pe*pk*cos(ang)+2*pe_*pk*cos(2*ang)));

  cout<<"mm "<<mm<<" mm2 "<<mm2<<endl;
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
