

void test(){

  TLorentzVector Rv,Vv,Lv,Pv,Bv,Hv;

  double pk  = 1.823;
  double pe  = 4.32;
  double pe_ = 2.1;
  double mk  = 0.493677;
  double me  = 0.510998928e-3;
  double mp  = 0.938272046;
  double mL  = 1.115683;
  
  double Ek  = sqrt(mk*mk + pk*pk);
  double Ee_ = sqrt(me*me + pe_*pe_);
  double Ee =  sqrt(me*me + pe*pe);
  Bv.SetPxPyPzE(0.0,0.0,pe,Ee);
  Rv.SetPxPyPzE(0.0,0.0,pk,Ek);
  Lv.SetPxPyPzE(0.0,0.0,pe_,Ee_);
  Pv.SetPxPyPzE(0.0,0.0,0.0,mp);

  
  Rv.RotateX(  (13.2*PI)/180. );
  Hv.SetPxPyPzE(0.0,0.0,0.0,mL);
  Lv.RotateX( -(13.2*PI)/180. );
  
  Vv = Bv -Lv;
  //  Vv.SetE(-Vv.E());
  Hv = Vv + Pv - Rv;
  //  Rv = Vv - Hv + Pv;

  
  cout<<"Pv "<<Vv.P()<<" Ev "<<Vv.E()<<endl;
  cout<<"PL "<<Hv.P()<<" EL "<<Hv.E()<<endl;
  cout<<"Pp "<<Pv.P()<<" Ep "<<Pv.E()<<endl;
  cout<<"Rp "<<Rv.P()<<" EK "<<Rv.E()<<endl;
  cout<<" theta_vk "<<(Rv.Theta()-Vv.Theta() )*180./PI<<" theta_ek "<<Rv.Theta()*180./PI<<" theta_ee_ "<<Lv.Theta()*180./PI<<endl;
  cout<<"ML "<<sqrt(Hv.E()*Hv.E() -Hv.P()*Hv.P())<<endl;
  cout<<"Q value "<<- Vv.Mag()<<endl;


  
  
}
