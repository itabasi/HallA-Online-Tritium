
void photon_flux_cal(){

  double Ee=4300;
  double Ee_=2100;
  double theta_ee=13.2;//degree
  double Q=0;
  double pv=0;
  double mp=948;
  double q=0;
  double Ev=0.0;
  double PI=3.14;
  double ep=0;
  double Flux;
  double alpha=1./137.;
  q=sqrt(Ee*Ee+Ee_*Ee_-2*Ee*Ee_*cos(theta_ee/180*PI));
  Q=sqrt(-(Ee-Ee_)*(Ee-Ee_)+q*q);
  Ev=(Ee-Ee_)+q*q/(2*mp);
  ep=1/(1+2*q*q/(Q*Q)*pow(tan(theta_ee/180*PI/2),2));

  Flux=alpha/(2*PI*PI*Q*Q)*Ev/(1-ep)*Ee_/Ee;

  cout<<"q:"<<q<<endl;
  cout<<"Q"<<Q<<endl;
  cout<<"Ev:"<<Ev<<endl;
  cout<<"epsilon:"<<ep<<endl;
  cout<<"Virtual Photon Flux [/s]:"<<Flux<<endl;



  // Cross Section //

  double NL=1600;
  double NS=396;
  double eff_SR=0.17;
  //double eff_K=
  double Nt=4.26e22;
  double dOme=0.006;//[sr]
  double e=1.602e-19;//[C]
  double C=4.754;//[C]
  double conv=1.0e-24;//[b/cm^2]
  double Nvp=C/e*Flux;

  double CSL=NL/(Nt*Nvp*eff_SR*dOme*conv);
  double CSS=NS/(Nt*Nvp*eff_SR*dOme*conv);  

  cout<<"N virtual photon:"<<Nvp<<endl;
  cout<<"Lambda Cross Sectioon [nb/sr]:"<<CSL*1e9<<endl;
  cout<<"Sigma Cross Sectioon [nb/sr]:"<<CSS*1e9<<endl;  
}
