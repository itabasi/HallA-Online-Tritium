#include <TVector3.h>
#include <TLorentzVector.h> 
void xsec_cm(){


  double pcm,Ecm,theta_cm,beta,gamma;
const  double Me = 0.510998928e-3;
const  double MK = 0.493677;
const  double Mp  = 0.938272046;
const  double ML  =  1.115683;
const  double MT  =  2.808921112;
  TLorentzVector Pv,Pe,Pe_,PK,PL,Pp,PT,Pcm;
  TVector3 B;
  double hrs_ang = 13.2*3.14/180.;
  double Ee = 4.318;
  double pK = 1.823;
  double pe_ = 2.1;
  double pp =0.0;
  double EK  = sqrt(MK*MK + pK*pK  );
  double Ee_ = sqrt(Me*Me + pe_*pe_);
  double Ep  = sqrt(Mp*Mp + pp*pp);
  double pe  = sqrt(Ee*Ee - Me*Me  );


  double acceptance = 6.0e-3; // msr
  double cos = 1.0 - acceptance/2./3.14;  
  Pe.SetPxPyPzE(0.0,0.0,pe,Ee);
  Pp.SetPxPyPzE(0.0,0.0,pp,Ep);
  PT.SetPxPyPzE(0.0,0.0,0.0,MT);
  Pe_.SetPxPyPzE(0.0,0.0,pe_,Ee_);
  Pe_.RotateX( - hrs_ang);
  Pv = Pe - Pe_;
  PK.SetPxPyPzE(0.0,0.0,pK,EK);
  PK.RotateX( + hrs_ang);


  Pcm = Pv + PT;
  //Pcm = Pv + Pp;
  B = Pcm.BoostVector();

  double u = PK.P()/PK.E();    

  //   double u = Pv.P()/Pv.E();
  //  double u = 1.0;

  beta = B.Mag();
  gamma = 1./sqrt(1.0 - beta*beta);
  Pv.Boost(-B);
  PT.Boost(-B);
  PK.Boost(-B);
  
  //  double u = fabs(Pv.P()/Pv.E();  
  //  double u = 1.0;
  double u_ = PK.P()/PK.E();    
  double cos_ = (u*cos- beta)/(1.0 - beta*u*cos)/u_;
  double sin_ = sqrt( 1.0 - cos_*cos_ );
  //   cout<<"u "<<u<<" beta "<<beta<<" cos "<<cos<<" cos_ "<<cos_<<endl;
  cout<<"u "<<u<<" beta "<<beta<<" deg "<<acos(cos)*180./3.14<<" def_ "<<acos(cos_)*180./3.14<<endl;
  cout<<"Pk "<<PK.P()<<" EK "<<PK.E()<<endl;

  
  double accept_cm = 2.0*3.14*(1.0 - cos_);
  double accept_lab =2.0*3.14*(1.0 - cos);
  //double accept_cm = 2.0*3.14*(1.0 - cos);
  double A =1.0;
  double Na = 6.03e23;
  double thickness = 70.8e-3; // [g/cm2]
  double NT = Na* thickness /A;
  double XS = 400.; //[nb/sr]
  double XS_Lab =400.;
  double cm_t_barn =1.0e24;
  double kaon_survival_ratio =0.17;
  double Ngamma =1.13316e14;  // Lambda H kine

  
  //  double Nhyp = 1600.;


  double XS_Lab_to_cm = pow( sin_*sin_ + gamma*gamma*pow(cos_ + beta/u_ ,2.0)  ,3./2.)/gamma/(1.0 + cos_*beta/u_);



  //  XS = XS_Lab*XS_Lab_to_cm;
  double Nhyp_Lab = XS_Lab/cm_to_barn*kaon_survival_ratio*NT*Ngamma/1.0e9*accept_lab;  

  
  double Nhyp = XS/cm_to_barn*kaon_survival_ratio*NT*Ngamma/1.0e9*accept_cm;
  
  cout<<"Nhyp "<<Nhyp<<" XS_cm "<<XS<<" Lab_to_cm "<<XS_Lab_to_cm<<endl;    
  //  cout<<"Nhyp_Lab "<<Nhyp_Lab<<" XS_lab "<<XS_Lab<<endl;  

  
  
  
}
