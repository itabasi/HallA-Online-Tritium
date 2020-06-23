

void test(){

  double Rp = 1.8000;
  double Lp = 2.1000;
  double theta,theta_L,phi,phi_L;

  theta = 0.1;
  phi = 0.2;
  theta_L =-0.3;
  phi_L = -0.23;
  

  double Rp_x,Rp_y,Rp_z,Lp_x,Lp_y,Lp_z;
  TVector3 Rv,Lv,Bv;

  Rp_z = Rp/sqrt(1.0 + theta*theta + phi*phi);
  Rp_x = Rp_z * theta;
  Rp_y = Rp_z * phi;

  Lp_z = Lp/sqrt(1.0 + theta_L*theta_L + phi_L*phi_L);
  Lp_x = Lp_z * theta_L;
  Lp_y = Lp_z * phi_L;  

  Rv.SetXYZ(Rp_x,Rp_y,Rp_z);
  Lv.SetXYZ(Lp_x,Lp_y,Lp_z);
  Rv.RotateX(13.2/180.*3.14);
  Lv.RotateX(-13.2/180.*3.14);

  cout<<" Rp "<<Rv.Mag()<<" Rpx "<<Rv.X()<<" Rpy "<<Rv.Y()<<" Rpz "<<Rv.Z()<<endl;
  cout<<" Lp "<<Lv.Mag()<<" Lpx "<<Lv.X()<<" Lpy "<<Lv.Y()<<" Lpz "<<Lv.Z()<<endl;


  double Rp_x2,Rp_y2,Rp_z2,Lp_x2,Lp_y2,Lp_z2;
  TVector3 Rv2,Lv2,Bv2;

  Rp_z2 = Rp/sqrt(1.0 + theta*theta + phi*phi);
  Rp_x2 = Rp_z2 * -theta;
  Rp_y2 = Rp_z2 * -phi;

  Lp_z2 = Lp/sqrt(1.0 + theta_L*theta_L + phi_L*phi_L);
  Lp_x2 = Lp_z2 * -theta_L;
  Lp_y2 = Lp_z2 * -phi_L;  

  Rv2.SetXYZ(Rp_x2,Rp_y2,Rp_z2);
  Lv2.SetXYZ(Lp_x2,Lp_y2,Lp_z2);
  
  Rv2.RotateX(-13.2/180.*3.14);
  Lv2.RotateX(13.2/180.*3.14);

  cout<<" Rp2 "<<Rv2.Mag()<<" Rpx2 "<<Rv2.X()<<" Rpy2 "<<Rv2.Y()<<" Rpz2 "<<Rv2.Z()<<endl;
  cout<<" Lp2 "<<Lv2.Mag()<<" Lpx2 "<<Lv2.X()<<" Lpy2 "<<Lv2.Y()<<" Lpz2 "<<Lv2.Z()<<endl;




  // Eloss //
  double hrs_ang =13.2/180.*3.14;
  double x = -hrs_ang + phi;
  double ELoss_R = -1.31749*sin(-4.61513*x) + 2.03687;

  
}
