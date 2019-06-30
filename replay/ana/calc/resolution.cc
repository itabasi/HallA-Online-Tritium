//===================================//
// Aurther Itabashi 2019 June 30th
// nnL experiment mass resolution 

void resolution(){

  double M_HYP,Mt;
  double Ee,Ee_,Ek,pe,pe_,pk,me,mk, th_ee_, th_ek;
  me = 0.511; // [MeV/c^2]
  mk = 493.7; // [MeV/c^2]
  double c = 1.0; // [/c]
  //==== Hydrogen run ======//
  Mt    = 938.27; // [MeV/c^2] proton mass
  M_HYP = 1115.68; //[MeV/c^2] Lambda mass
  Ee    = 4300.0; // [MeV]
  pe_   = 2100.0; // [MeV/c]
  pk=   = 1800.0; // [MeV/c]
  th_ee_= 13.2*3.14/180.;
  th_ek = 13.2*3.14/180.;
  Ee_= sqrt(me*me + pe_*pe_);
  Ek = sqrt(mk*mk + pk*pk);
  pe = sqrt(Ee*Ee - me*me);
  
  double dm,dEe,dEe_,dEk,dpe,dpe_,dpk,dth_ee_,dth_ek,dth_e_k; //sigma
  doubl sig_m,sig_Ee,sig_Ee_,sig_Ek,sig_pe,sig_pe_,sig_pk,sig_th_ee_,sig_th_ek,sig_th_e_k;


  //====== Resolution of HRS ============//
  sig_Ee  = Ee * 2.0e-4; // [MeV]
  sig_Ee_ = Ee_* 2.0-4;  // [MeV]
  sig_Ek  = Ek * 2.0e-4; // [MeV]

  sig_pe  = pe * 2.0e-4; // [MeV/c]
  sig_pe_ = pe_* 2.0e-4; // [MeV/c]
  sig_pk  = pk * 2.0e-4; // [MeV/c]

  sig_th_ee_ = 6,0e-3; //[rad]
  sig_th_ek  = 6.0e-3; //[rad]
  sig_th_e_k = sqrt(sig_th_ee_*sig_th_ee_ + sig_th_ek*sig_th_ek); //[rad]
  
  double dEe=1./(2.0*M_HYP)*(2Ee)*sig_Ee;
  double dEe_=1./(2.0*M_HYP)*(-2Ee_)*sig_Ee_;
  double dEk=1./(2.0*M_HYP)*(-2Ek)*sig_Ek;  
  double dth_ee_=1./(2.0*M_HYP)*(-2.0*pe*pe_*sin(th_ee_))*sig_th_ee_;
  double dth_ek =1./(2.0*M_HYP)*(-2.0*pe*pk*sin(the_ek))*sig_th_ek;
  double dth_e_k=1./(2.0*M_HYP)*(+2.0*pe_*pk*sin(th_e_k))*sig_th_e_k;
  double dpe =1./(2.0*M_HYP)*(-2.0*pe+2.0*pk*cos(th_ek)+2.0*pe_*cos(th_ee_))*sig_pe;
  double dpe_=1./(2.0*M_HYP)*(-2.0*pe_+2.0*pe*cos(th_ee_)-2.0*pk*cos(th_e_k))*sig_pe_;
  double dpk=1./(2.0*M_HYP)*(-2.0*pk+2.0*pe*cos(th_ek)-2*pe_*cos(th_e_k))*sig_pk;

  //===== Mass resolution no correlation ==============//

  dm=sqrt(dEe*dEe + dEe_*dEe_ + dEk*dEk + dpe*dpe + dpe_*dpe_ + dpk*dpk +
	  dth_ee_*dth_ee_ + dth_ek*dth_ek + dth_e_k*dth_e_k);

  

  
}
