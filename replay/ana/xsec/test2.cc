

void test2(){

  TVector3 Rv,Vv,Lv,Pv,Bv,Hv;

  double pk  = 1.823;
  double pe  = 4.32;
  double pe_ = 2.2;
  double mk  = 0.493677;
  double me  = 0.510998928e-3;
  double mp  = 0.938272046;
  double mL  = 1.115683;
  
  double Ek  = sqrt(mk*mk + pk*pk);
  double Ee_ = sqrt(me*me + pe_*pe_);
  double Ee  = sqrt(me*me + pe*pe);
  Bv.SetXYZ(0.0,0.0,pe);
  Rv.SetXYZ(0.0,0.0,pk);
  Lv.SetXYZ(0.0,0.0,pe_);
  Pv.SetXYZ(0.0,0.0,0.0);


  Rv.RotateX( + 13.2*PI/180. );
  //  Hv.SetXYZ(0.0,0.0,0.0);
  Lv.RotateX( - 13.2*PI/180. );
  Vv = Bv -Lv;


  //  cout<<" Ee -Ee_ "<<Ee-Ee_<<" pe - pe_ "<<pe -pe_<<endl;
  //  cout<<" pe -pe_ "<<Bv.Mag() -Lv.Mag()<<" Vv "<<Vv.Mag()<<endl;
  double MH = sqrt( (Ee + mp - Ee_ -Ek)* (Ee + mp - Ee_ -Ek) -   (Bv -Lv - Rv)*(Bv -Lv -Rv)   );
  //  Vv.SetE(-Vv.E());
  //  Hv = Vv + Pv - Rv;
  //  Rv = Vv - Hv + Pv;

  Hv = Vv - Rv -Pv; 
  
  cout<<"Q "<<Vv.Mag() - (Ee -Ee_)<<endl;
  cout<<"ML "<<MH<<" PL "<<Hv.Mag()<<endl;
  cout<<"Magp "<<Pv.Mag()<<endl;
  cout<<"Rp "<<Rv.Mag()<<endl;
  cout<<" theta_vk "<<(Rv.Theta()-Vv.Theta() )*180./PI<<" theta_ek "<<Rv.Theta()*180./PI<<" theta_ee_ "<<Lv.Theta()*180./PI<<endl;



  
  
}
