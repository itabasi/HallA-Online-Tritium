

void virtual_photon_cal(){
  double Gamma;
  double alpha=1./137.;
  double pi=3.14;
  double me=0.510 //[MeV]
  double mp=938.27;
  double Ee=4200 //[MeV]
  double Ee_=2200;
  double deg_ee=12.5; //degree
  double theta_ee=def_ee/360*2*pi;
  double pv;
  double omega=Ee-Ee_;
  pv=sqrt(Ee*E+Ee_*Ee_-2*Ee*Ee_*cos(theta_ee));
  double Ev=omega+pv*pv/(2*mp);


}
