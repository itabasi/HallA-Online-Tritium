//===================================//
// Aurther Itabashi 2019 June 30th
// nnL experiment mass resolution 

void resolution(){

  double M_HYP,Mt,M_nnL;
  double Ee,Ee_,Ek,pe,pe_,pk,me,mk, th_ee_, th_ek,th_e_k;
  me = 0.511; // [MeV/c^2]
  mk = 493.7; // [MeV/c^2]
  double c = 1.0; // [/c]
  bool nnL=false;
  //==== Hydrogen run ======//
  Mt    = 938.27; // [MeV/c^2] proton mass
  M_HYP = 1115.68; //[MeV/c^2] Lambda mass
  if(nnL){
  Mt    = 2808.921112;
  M_HYP = 2998.;
  pe_   = 2200.0; // [MeV/c]
  }
  
  Ee    = 4300.0; // [MeV]
  pe_   = 2100.0; // [MeV/c] 

  pk    = 1800.0; // [MeV/c]
  th_ee_= 13.2*3.14/180.;
  th_ek = 13.2*3.14/180.;


  /*  
  // ====E05-115 Gogami exp. =============
  Ee    = 2344.0; // [MeV]
  pk   = 1200.0; // [MeV/c]
  pe_    = 844.0; // [MeV/c]
  th_ee_= 6.5*3.14/180.;
  th_ek = 5.7*3.14/180.;
  //=====================================
  */
 
  th_e_k = sqrt(th_ek*th_ek + th_ee_*th_ee_);
  Ee_= sqrt(me*me + pe_*pe_);
  Ek = sqrt(mk*mk + pk*pk);
  pe = sqrt(Ee*Ee - me*me);

  
  //  double dm,dEe,dEe_,dEk,dpe,dpe_,dpk,dth_ee_,dth_ek,dth_e_k; //sigma
  double sig_m,sig_Ee,sig_Ee_,sig_Ek,sig_pe,sig_pe_,sig_pk,sig_th_ee_,sig_th_ek,sig_th_e_k;


  //====== Resolution of HRS ============//
  
  sig_pe  = pe *1.0e-4; // [MeV/c]
  sig_pe_ = pe_*2.0e-4; // [MeV/c]
  sig_pk  = pk *2.0e-4; // [MeV/c]

  sig_th_ee_ = sqrt( pow(6.0e-3,2)+pow(2.3e-3,2) ); //[rad]
  sig_th_ek  = sqrt( pow(6.0e-3,2)+pow(2.3e-3,2) ); //[rad]
  sig_th_ek  = 6.67e-3; //[rad]  
  sig_th_e_k = sqrt(sig_th_ee_*sig_th_ee_ + sig_th_ek*sig_th_ek); //[rad]
  
  sig_Ee  = sig_pe * (2.0*pe  /sqrt(me*me + pe *pe )); // [MeV]
  sig_Ee_ = sig_pe_* (2.0*pe_ /sqrt(me*me + pe_*pe_)); // [MeV]
  sig_Ek  = sig_pk * (2.0*pk  /sqrt(mk*mk + pk *pk )); // [MeV]



  /*
  //====== Resolution of Gogami exp ============//
  sig_pe  = pe *1.0e-4; // [MeV/c]
  sig_pe_ = pe_*4.2e-4; // [MeV/c]
  sig_pk  = pk *2.0e-4; // [MeV/c]

  sig_th_ee_ = 4.0e-3/1.54; //[rad] (FWHM)
  sig_th_ek  = 0.4e-3/1.54; //[rad] (FWHM)
  sig_th_e_k = sqrt(sig_th_ee_*sig_th_ee_ + sig_th_ek*sig_th_ek); //[rad]
 
  sig_Ee  = sig_pe  / ( pe  /sqrt(me*me + pe *pe )); // [MeV]
  sig_Ee_ = sig_pe_ / ( pe_ /sqrt(me*me + pe_*pe_)); // [MeV]
  sig_Ek  = sig_pk  / ( pk  /sqrt(mk*mk + pk *pk )); // [MeV]  
  */

  

  
  double dEe  =1./(2.0*M_HYP)*(2.0*Ee  + 2.0*(Mt -Ek -Ee_))*sig_Ee;
  double dEe_ =1./(2.0*M_HYP)*(2.0*Ee_ - 2.0*(Ee +Mt -Ek ))*sig_Ee_;
  double dEk  =1./(2.0*M_HYP)*(2.0*Ek  - 2.0*(Ee +Mt -Ee_))*sig_Ek;  

  double dth_ee_ =1./(2.0*M_HYP)*(-2.0*pe*pe_*sin(th_ee_))*sig_th_ee_;
  double dth_ek  =1./(2.0*M_HYP)*(-2.0*pe*pk*sin(th_ek))*sig_th_ek;
  double dth_e_k =1./(2.0*M_HYP)*(+2.0*pe_*pk*sin(th_e_k))*sig_th_e_k;

  


  double dpe =1./(2.0*M_HYP)*(-2.0 * pe  + 2.0*pk*cos(th_ek)  + 2.0  * pe_* cos(th_ee_)) * sig_pe;
  double dpe_=1./(2.0*M_HYP)*(-2.0 * pe_ + 2.0*pe*cos(th_ee_) - 2.0  * pk * cos(th_e_k)) * sig_pe_;  
  double dpk =1./(2.0*M_HYP)*(-2.0 * pk  + 2.0*pe*cos(th_ek)  - 2.0  * pe_* cos(th_e_k)) * sig_pk;

  //===== Mass resolution no correlation ==============//

  double  dm=sqrt(dEe*dEe + dEe_*dEe_ + dEk*dEk + dpe*dpe + dpe_*dpe_ + dpk*dpk +
		  dth_ee_*dth_ee_ + dth_ek*dth_ek + dth_e_k*dth_e_k);
  

  cout<<"dEe "<<dEe*1000.<<" dEe_ "<<dEe_*1000.<<" dEk "<<dEk*1000.<<" [keV] (FWHM)"<<endl;
  cout<<"dth_ee_ "<<dth_ee_<<" dth_ek "<<dth_ek<<" dth_e_k "<<dth_e_k<<endl;
  cout<<"dpe "<<dpe*1000.<<" dpe_ "<<dpe_*1000.<<" dpk "<<dpk*1000.<<" [keV] (FWHM)"<<endl;
  cout<<"dm_L "<<dm<<" MeV (FWHM)"<<endl;
  cout<<"dm_nnL "<<dm*M_HYP/M_nnL<<" MeV (FWHM)"<<endl;
}
