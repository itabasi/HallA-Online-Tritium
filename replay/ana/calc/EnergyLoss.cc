// Calculation of energy loss of electron beam //
// path through in target material by electron //
// Aurther Itabashi 2019 07 08 //

void EnergyLoss(){


  double Na,me,re,rho,c,beta,tau,I,F,delta,gamma,C,Z,A,alpha,eta,X,C0;
  double N,E0;

  //=====Input paramters =====//
  tau=4300.0; //beam energy [MeV]
  rho=0.144/25.; //hydrogen density [g/cm^3]
  Z=1.;
  A=1.;
  //==========================//


  //==== collison Loss by electron====//
  
  
  Na=6.022e23;//Avogadro's number [1/mol]
  me=0.511; // electron mass [MeV/c^2]
  re=2.817e-13; // classical electron radius [cm]
  c=1.0;
  alpha=1./137.; // fine stracture's number
  I=(12.*Z+7.)*1.0e-6; //[MeV] if Z<13
  //I=(9.76 + 58.8*pow(Z,-1.19))*1.0e-6; // [MeV] if Z>13 

  beta=sqrt(tau*tau-me*me)/tau;
  gamma=1./sqrt(1.-beta*beta);
  eta=beta*gamma;
  X=log10(beta*gamma);
  delta=0.0; // low a few % effects
  C=( 0.42237*pow(eta,-2.) + 0.0304043*pow(eta,-4.) - 0.00038106*pow(eta,-6.))*1.0e-6*I*I
    +( 3.8501908*pow(eta,-2.) - 0.1667989* pow(eta,-4.) + 0.00157955*pow(eta,-6.))*1.0e-9*I*I*I;
  F=1.0-pow(beta,2.)+(pow(tau,2.)/8.-(2*gamma+1)*log(2.0))/pow(tau+1,2);//electon

  double Cb=0.1535; //[MeV cm^2/g]
  double dE_c=2.*3.14*Na*re*re*me*c*c*rho*Z/A/(beta*beta)*
    (log(pow(tau,2.)*(tau+2.)/(2.*(I/(me*pow(c,2))))+F-delta-2.*C/Z));

  
  //====== Bremsstrahlung ===========//

  E0=tau;
  N=rho*Na/A;
  double a=Z/137.;
  double f=a*a*( pow(1+a*a,-1.) + 0.20206 - 0.0396*a*a + 0.0083*pow(a,4.) - 0.002*pow(a,6.));
  double Phi=4*Z*Z*re*re*alpha*(log(183*pow(Z,-1./3.)) + 1./18. -f  ); //if Ee >>137mec^2Z^(-1/3)

  double dE_r=N*E0*Phi;
  double dE=dE_c + dE_r;
  cout<<"total Energy loss [MeV]: "<<dE<<" dE_c: "<<dE_c<<" dE_r: "<<dE_r<<endl; 

  
}
