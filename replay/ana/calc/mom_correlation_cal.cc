// Momentum Setup Calculation //

void mom_correlation_cal(){

  //=== parameters ====//

  double E,Ee,P,Pe,Pv,Ev,Ae,Av,Avk,me,Ib,Ne,edeg;
  double mt,mh,mk,Et,Eh,Ek,Ph,Pt,Pk,Avh,be,W;

  me=0.51;//[MeV] :electron mass
  mk=493.7;//[MeV] kaion mass


  //--- Given parameters----//

  //======= Ca exp. (E12-15-008) SetUp =============================//
  /*
  E=4500;  // [MeV] : beam energy
  Ib=20;//[myuA] :beam flux
  Pe=2180;//[MeV] :scatter electron momentum
  edeg=7.0; // scattered electron degree
  
  //-- mt [MeV] target mass ----------//
   mt=938;// proton mass
  //mt=2832; //3H mass
  //mt=2808;
  Pt=0;//[MeV] target momentum

  //---mh [MeV] Hyperon mass--------//
  //mh=1115;// lambda mass
  // mh=1192; //sigma0 mass
   mh=2995+be; //nnL mass (be is binding energy)
   be=0;//binding energy is nothing
  

 */
  //===== END of Ca exp. Parameters ==================//


  //====== nnL exp. (E12-18-003) SetUp ======================//
  bool nnL=false; // nnL mode
  if(nnL){
  E=4300;  // [MeV] : beam energy
  Ib=20;//[myuA] :beam flux
  //  Pe=2180;//[MeV] :scatter electron momentum
  Pe=2180;//[MeV] :scatter electron momentum
  edeg=13.2; // scattered electron degree

  //-- mt [MeV] target mass ----------//
  //    mt=938;// proton mass
  //  mt=2832; //3H mass
  mt=2808;
  
  Pt=0;//[MeV] target momentum

  //---mh [MeV] Hyperon mass--------//
  // mh=1115;// lambda mass
  //   mh=1192; //sigma0 mass
  // mh=2995+be; //nnL mass (be is binding energy)
  mh= 2994.9 -100.;
   be=0;//binding energy is nothing
  }else{

  E=4300;  // [MeV] : beam energy
  Ib=20;//[myuA] :beam flux
  //  Pe=2180;//[MeV] :scatter electron momentum
  Pe=2100;//[MeV] :scatter electron momentum
  edeg=13.2; // scattered electron degree

  //-- mt [MeV] target mass ----------//
  //    mt=938;// proton mass
  //  mt=2832; //3H mass
  mt=938.;
  
  Pt=0;//[MeV] target momentum

  //---mh [MeV] Hyperon mass--------//
  // mh=1115;// lambda mass
  //   mh=1192; //sigma0 mass
  // mh=2995+be; //nnL mass (be is binding energy)
  mh= 1115.68;


  }

 //==== END of nnL exp. parameters ==============//
  


   

  //---- calcurated parameters---//

  Ne=Ib*6.24e+12;//[1/s] electron current
  // Ee=2725;//[MeV] :scatter electron energy
  P=sqrt(pow(E,2)+pow(me,2)); //[MeV] :beam momentum
  //Ee [MeV] e' energy
  //Ev [MeV] virtualphoton energy
  // Ek [MeV]kaion enrgy
  Ae=edeg*(3.14/180.0); // [rad] e' angle 
  // Av=17.5*(3.14/180.0);//[rad]  virtual photon angle
  Avk=0;//[rad] between virtual and kaion angle
  //cout<<"angle between virtual photon and kaion"<<Avk<<" [rad]"<<endl;


  
  
  //////////////////////////////Scattering plane//////////////////////////////////////

  ////energy consrtvation////
  // E=Ee+Ev
  ////,momentum conservation//////
  //P=Pe+Pk
  //Pe*sin(Ae)=Pv*sin(Av):
  double q;
  // q=sin(Ae)/sin(Av); 
  double pe;
  pe=Pe;

  //==  Coment out ===//
  cout<<endl;
  cout<<"======  Given Parameters ======  "<<endl;
  cout<<"Electron Beam Energy [MeV] is    "<<E<<endl;
  cout<<"Scattered Electron Degree [deg] is       "<<edeg<<endl;
  cout<<endl;
  cout<<"====== Result of Calculation ======  "<<endl;
  //cout<<"Pe is "<<pe<<endl;
  for(int i=-1;i<2;i++)
    {

      cout<<endl;
      cout<<"        Pe +  "<<i<<"*dPe"<<endl;
      double dPe;// Pe' acceptance bite
      //      dPe=0.045*pe;
      dPe=0.45*pe;
      Pe=pe+i*dPe;
      double Q;
      Ee=Pe;
      Q=sqrt(4*E*Ee*pow(sin(Ae/2),2));
	     Pv=sqrt(pow(E-Ee,2)+pow(Q,2));
	     Ev=E-Ee;
      Av=acos((E-Ee*(1-pow(Q,2)/(2*E*Ee)))/Pv);//[rad]

      double Iv,e,alpha;
      alpha=1/137;
      e=1/(1+2*pow(Pv/Q,2)*pow(tan(Ae/2),2));
      Iv=1/(137*2*pow(3.14,2)*pow(Q,2))*(Ev*Ee)/(1-e)/E;
      W=mt*mt +2.*mt*Ev-Ev*Ev;
      cout<<"virtual angule is "<<Av*180/3.14<<endl;
      cout<<"electron momentum is "<<Pe<<endl;
      cout<<"photon momentum is "<<Pv<<endl;
      cout<<"virtual photon flux [1/s]"<<Iv*Ne<<endl;
      cout<<"e is "<<e<<endl;
      cout<<"W is "<<W<<endl;
      cout<<"Q is "<<Q<<endl;

  ////////////////////////////Reaction plane//////////////////////////////////////////
  

  ///////////////////////////////////////////
  // CPk^2+DPk+B=0  (Pk calculation) //////////

  double A,B,C,D;//A,B,C,D is constance
  A=pow(mh,2)+pow(Pv,2)-pow(Ev+mt,2)-pow(mk,2);
  B=(pow(A,2)-4*pow(Ev+mt,2)*pow(mk,2))/4;
  C=pow(Pv*cos(Avk),2)-pow(Ev+mt,2);
  D=-A*Pv*cos(Avk);
  /* cout<<"A is "<<A<<endl;
  cout<<"B is "<<B<<endl;
  cout<<"C is "<<C<<endl;
  cout<<"D is "<<D<<endl;
  cout<<"sqrt is "<<pow(D,2)-4*C*B<<endl;
  */
  Pk=(-D-sqrt(pow(D,2)-4*C*B))/(2*C);
  cout<<"kaion momentum is "<<Pk<<endl;
  cout<<endl;
 

    }



}
