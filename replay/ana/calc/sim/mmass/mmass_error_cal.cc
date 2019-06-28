// Missing Mass Error Calculation //
// MADE BY ITABASHI              //
// nnL exeriment (E12-18-003)    //


#include <iostream>
using namespace std;

void mmas_error_cal(){

  double Ee,Ee_,Ek,m_H,dPe,dPe_,dPk,dee_,dek,de_k,dm,mt,me,mk,Pe,Pk,Pe_,ee_,ek,e_k,dee_,dek,de_k,DPe,DPe_,DPk,Dee_,Dek,De_k,d_ee_,d_ek,d_e_k,dm,m_H;
  
 
 //======= Given Parameter ======//

  mt=2804; // [MeV]:tritium target mass
  // mt=983.27;//[MeV] :proton target mass
  me=0.510998;// [MeV]:electron mass
  mk=494;// [MeV]:kaion mass



  //========================================//
  //============== Parameter ==============//
  //~======================================//

  Ee=4300;// [MeV]: beam energy
  Pe=sqrt(pow(Ee,2)+pow(me,2));//[MeV] beam momentum
  Pe_=2180;// [MeV]: scattering electrom momentum
  Pk=1800;//  [MeV]: kaion momentum
  ee_=12.5*3.14/180;//[rad] :electron-scattering electron angle
  ek=12.5*3.14/180;//   [rad];electron-kaion angle
  e_k=ee_+ek;// [rad] :scattering electron-kaion angle
  double dpe,dpe_,dpk,dee_,dek,de_k;
  
  /////////resolution rate//////////
  dpe=3e-5;   /////Beam resolution
  dpe_=2.5e-4; /////HRS resolution
  dpk=2.5e-4;  /////HRS resolution

  ////////resolution/////////////
  dPe=Pe*dpe;
  dPe_=Pe_*dpe_;
  dPk=Pk*dpk;

  //dPe=0.9048;//[MeV] :electron momentum difference 
  //dPe_=0.5450;//[MeV]:scattering electron momentum difference
  //dPk=0.3000;//[MeV] : kaion momentum difference

  ////////////////anguler resolution////////////////////////

  dee_=0.003; ///[rad]
  dek=0.003; //[rad]
      de_k=sqrt(pow(dee_,2)+pow(dek,2));//[rad]


  // dee_=0.002;//[rad] :electron-scattering electron angle difference
  //dek=0.002;//[rad] :electron-kaion angle difference
  //de_k=0.0028;//[rad] :scattering electron-kaion angle difference
  Ee_=sqrt(pow(me,2)+pow(Pe_,2));
  Ek=sqrt(pow(mk,2)+pow(Pk,2));
  //cout<<"Ee_="<<Ee_<<endl;
  //cout<<"Ek="<<Ek<<endl;
  m_H=sqrt(pow(mt,2)+pow(me,2)+pow(me,2)+pow(mk,2)+2*mt*(Ee-Ee_-Ek)+2*(-Ee*Ee_-Ee*Ek+Ee_*Ek)+2*(Pe*Pe_*cos(ee_)+Pe*Pk*cos(ek)-Pe_*Pk*cos(e_k)));
  //cout<<"m_H="<<m_H<<endl;
        DPe=dPe/m_H*(Pe_*cos(ee_)+Pk*cos(ek)+mt-Ee_-Ek);
        DPe_=dPe_/m_H*(Pe*cos(ee_)-Pk*cos(e_k)+Pe_/Ee_*(-mt-Ee+Ek));
        DPk=dPk/m_H*(Pe*cos(ek)-Pe_*cos(e_k)+Pk/Ek*(-mt-Ee+Ee_));
        Dee_=-dee_/m_H*(Pe*Pe_*sin(ee_));
	Dek=-dek/m_H*(Pe*Pk*sin(ek));
	De_k=de_k/m_H*(Pe_*Pk*sin(e_k));
	double d_Pe,d_Pe_,d_Pk;
	d_Pe=DPe/dPe;
	d_Pe_=DPe_/dPe_;
	d_Pk=DPk/dPk;
        d_ee_ =-(Pe*Pe_*sin(ee_))/m_H;
	d_ek =-(Pe*Pk*sin(ek))/m_H;
	d_e_k =(Pe_*Pk*sin(e_k))/m_H;

	dm=sqrt(pow(DPe,2)+pow(DPe_,2)+pow(DPk,2)+pow(Dee_,2)+pow(Dek,2)+pow(De_k,2));
	cout<<"dPe="<<DPe<<endl;
	cout<<"dPe'="<<DPe_<<endl;
	cout<<"dPk="<<DPk<<endl;
	cout<<"dee'="<<Dee_<<endl;
	cout<<"dek="<<Dek<<endl;
	cout<<"de'k="<<De_k<<endl;
	cout<<"dm_H= "<<dm<<endl;
	cout<<"m_H="<<m_H<<endl;
	//cout<<"d_ee'="<<d_ee_<<endl;
	cout<<"total error "<<dm<<endl;

	///////////////produce Lambda particle  in center of the target/////////////
	double Pc;//center of targtet by electron
	double dMH;//missing mass error
	
	/*cout<<"MPV [MeV] is "<<endl;
	cin>>Pc;
      
	//assume dP^m_e=dP^m_e_=dP^m_k;


	dMH=-2*Pc*(d_Pe+d_Pe_+d_Pk);
	cout<<"masserror [MeV] "<<dMH<<endl; 
	*/

}

