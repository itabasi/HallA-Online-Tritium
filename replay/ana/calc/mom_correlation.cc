const double  me=0.51;//[MeV] :electron mass
const double  mk=493.7;//[MeV] kaion mass

//double Pk(double Ee, double theta_ee_, double pe_, double theta_ek){
//}

double Accep(double* x,double *par){


  double mt = par[0] ; // target mass [MeV]
  double mh = par[1] ; // hypernuclei mass [MeV]
  double Pe = x[0]   ; // [MeV] : scattered electron momentum
  double scale = par[2]; // if scale is 1 scale 2.2 -> 2.1 GeV


  double  E  = 4300.;  // [MeV] : beam energy
  double  Ib = 20;//[myuA] :beam flux
  double edeg = 13.2; // scattered electron degree
  double  Ae = edeg*(3.14/180.0); // [rad] e' angle 
  double  Avk= 0;//[rad] between virtual and kaion angle


  double  Ee=Pe;
  double  Q=sqrt(4*E*Ee*pow(sin(Ae/2),2));
  double  Pv=sqrt(pow(E-Ee,2)+pow(Q,2));
  double  Ev=E-Ee;
  double  Av=acos((E-Ee*(1-pow(Q,2)/(2*E*Ee)))/Pv);//[rad]
  double  P  = sqrt(pow(E,2)+pow(me,2)); //[MeV] :beam momentum


  double alpha=1/137;
  double e=1/(1+2*pow(Pv/Q,2)*pow(tan(Ae/2),2));
  double Iv=1/(137*2*pow(3.14,2)*pow(Q,2))*(Ev*Ee)/(1-e)/E;
  double W=mt*mt +2.*mt*Ev-Ev*Ev;


  // CPk^2+DPk+B=0  (Pk calculation) //////////
  double A,B,C,D;//A,B,C,D is constance
  A=pow(mh,2)+pow(Pv,2)-pow(Ev+mt,2)-pow(mk,2);
  B=(pow(A,2)-4*pow(Ev+mt,2)*pow(mk,2))/4;
  C=pow(Pv*cos(Avk),2)-pow(Ev+mt,2);
  D=-A*Pv*cos(Avk);

  double pk=(-D-sqrt(pow(D,2)-4*C*B))/(2*C);

  return pk;

};



double Pk(double* x,double *par){


  double mt = par[0] ; // target mass [MeV]
  double mh = par[1] ; // hypernuclei mass [MeV]
  double Pe = x[0]   ; // [MeV] : scattered electron momentum

  double  E  = 4300.;  // [MeV] : beam energy
  double  Ib = 20;//[myuA] :beam flux
  double edeg = 13.2; // scattered electron degree
  double  Ae = edeg*(3.14/180.0); // [rad] e' angle 
  double  Avk= 0;//[rad] between virtual and kaion angle


  double  Ee=Pe;
  double  Q=sqrt(4*E*Ee*pow(sin(Ae/2),2));
  double  Pv=sqrt(pow(E-Ee,2)+pow(Q,2));
  double  Ev=E-Ee;
  double  Av=acos((E-Ee*(1-pow(Q,2)/(2*E*Ee)))/Pv);//[rad]

  double  P  = sqrt(pow(E,2)+pow(me,2)); //[MeV] :beam momentum


  double  alpha=1/137;
  double e=1/(1+2*pow(Pv/Q,2)*pow(tan(Ae/2),2));
  double Iv=1/(137*2*pow(3.14,2)*pow(Q,2))*(Ev*Ee)/(1-e)/E;
  double W=mt*mt +2.*mt*Ev-Ev*Ev;


  // CPk^2+DPk+B=0  (Pk calculation) //////////
  double A,B,C,D;//A,B,C,D is constance
  A=pow(mh,2)+pow(Pv,2)-pow(Ev+mt,2)-pow(mk,2);
  B=(pow(A,2)-4*pow(Ev+mt,2)*pow(mk,2))/4;
  C=pow(Pv*cos(Avk),2)-pow(Ev+mt,2);
  D=-A*Pv*cos(Avk);

  double pk=(-D-sqrt(pow(D,2)-4*C*B))/(2*C);

  return pk;

};

double Pk2(double* x,double *par){


  double mt = par[0] ; // target mass [MeV]
  double mh = par[1] ; // hypernuclei mass [MeV]
  double Pe = x[0]/2218.*2100.   ; // [MeV] : scattered electron momentum

  double  E  = 4300.;  // [MeV] : beam energy
  double  Ib = 20;//[myuA] :beam flux
  double edeg = 13.2; // scattered electron degree
  double  Ae = edeg*(3.14/180.0); // [rad] e' angle 
  double  Avk= 0;//[rad] between virtual and kaion angle


  double  Ee=Pe;
  double  Q=sqrt(4*E*Ee*pow(sin(Ae/2),2));
  double  Pv=sqrt(pow(E-Ee,2)+pow(Q,2));
  double  Ev=E-Ee;
  double  Av=acos((E-Ee*(1-pow(Q,2)/(2*E*Ee)))/Pv);//[rad]

  double  P  = sqrt(pow(E,2)+pow(me,2)); //[MeV] :beam momentum


  double  alpha=1/137;
  double e=1/(1+2*pow(Pv/Q,2)*pow(tan(Ae/2),2));
  double Iv=1/(137*2*pow(3.14,2)*pow(Q,2))*(Ev*Ee)/(1-e)/E;
  double W=mt*mt +2.*mt*Ev-Ev*Ev;


  // CPk^2+DPk+B=0  (Pk calculation) //////////
  double A,B,C,D;//A,B,C,D is constance
  A=pow(mh,2)+pow(Pv,2)-pow(Ev+mt,2)-pow(mk,2);
  B=(pow(A,2)-4*pow(Ev+mt,2)*pow(mk,2))/4;
  C=pow(Pv*cos(Avk),2)-pow(Ev+mt,2);
  D=-A*Pv*cos(Avk);

  double pk=(-D-sqrt(pow(D,2)-4*C*B))/(2*C);

  return pk;

};

double mmass(double *x,double *par){


  double mt  = par[0];
  double mh  = par[1];

  double pe_ = par[2];
  double pk=x[0];
  //  double pe_ = x[0];
  //  double pk  = par[2];

  double Ee  = 4300.; 
  double pe  = sqrt(Ee*Ee -me*me);




  double Ee_  = sqrt(pe_*pe_ + me*me);
  double Ek =   sqrt(pk*pk + mk*mk);

  double ang_ee_ = 13.2*3.14/180.;
  double ang_ek  = 13.2*3.14/180.;
  double ang_e_k = 13.2*3.14/180.*2.0;

  double  mm=sqrt((Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
		  (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang_ee_)-2*pe*pk*cos(ang_ek)+2*pe_*pk*cos(ang_e_k)));

  double dm = mm - mh;
  return dm;

};

double mmass2(double *x,double *par){


  double mt  = par[0];
  double mh  = par[1];

  //  double pe_ = par[2];
  //  double pk=x[0];

  double Ee  = 4300.; 
  double pe  = sqrt(Ee*Ee -me*me);
  double pe_ = x[0]/2180.*2100.;
  double pk  = par[2];




  double Ee_  = sqrt(pe_*pe_ + me*me);
  double Ek =   sqrt(pk*pk + mk*mk);

  double ang_ee_ = 13.2*3.14/180.;
  double ang_ek  = 13.2*3.14/180.;
  double ang_e_k = 13.2*3.14/180.*2.0;

  double  mm=sqrt((Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
		  (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang_ee_)-2*pe*pk*cos(ang_ek)+2*pe_*pk*cos(ang_e_k)));

  double dm = mm - mh;
  return dm;

};


double B(double par_0,double par_1, double par_2, double par_3){

  double  mt  = par_0; // target mass
  double  mh  = par_1;
  double  pk  = par_2;
  double  pe_ = par_3;

  // x[0] is missing mass A

  double Ee  = 4300.;
  double pe  = sqrt(Ee*Ee -me*me);
  double Ee_ = sqrt(pe_*pe_ +me*me);//pe_/sqrt(pe_*pe_ + me*me);
  double Ek =  sqrt(pk*pe + mk*mk); //pk/sqrt(pk*pk + mk*mk);

  double ang_ee_ = 13.2*3.14/180.;
  double ang_ek  = 13.2*3.14/180.;
  double ang_e_k = 13.2*3.14/180.*2.0;

  double  mm = sqrt((Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
		    (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang_ee_)-2*pe*pk*cos(ang_ek)+2*pe_*pk*cos(ang_e_k)));
  
  double B = mm - mh;
  //  cout<<"mm "<<mm<<" mh "<<mh<<" mt "<<mt<<" B "<<B<<endl;
  //  cout<<" pe "<<pe<<" pe_ "<<pe_<<" pk "<<pk<<" ang "<<ang_ee_<<" cos "<<cos(ang_ee_)<<endl;
  return  B;

};



void mom_correlation(){


  TF1* fL=new TF1("fL","Pk",2000,2400,2);
  TF1* fL2=new TF1("fL2","Pk2",2000,2400,2);
  TF1* fS=new TF1("fS","Pk",2000,2400,2);
  TF1* fnnL=new TF1("fnnL","Pk",2000,2400,2);
  TF1* fAl=new TF1("fAl","Pk",2000,2400,2);
  TF1* fAl2=new TF1("fAl2","Pk2",2000,2400,2);
  fAl->SetNpx(2000);
  fnnL->SetNpx(2000);
  fL->SetNpx(2000);
  fS->SetNpx(2000);
  fL2->SetNpx(2000);
  fAl2->SetNpx(2000);
  fnnL->SetLineColor(2);
  fL->SetLineColor(3);
  fL2->SetLineColor(5);
  fS->SetLineColor(4);
  fAl2->SetLineColor(7);
  double mp   =  0.938272046e3;
  double mL   =  1.115683e3;
  double mS0  = 1.192642e3;
  double mT   = 2.808921112e3;
  double mnnL = 2.9948138266e3;
  double mAl  = 25.1265e3;
  double mMgL = 25.312188e3;


  fL->SetParameters(mp,mL);
  fL2->SetParameters(mp,mL);
  fS->SetParameters(mp,mS0);
  fnnL->SetParameters(mT,mnnL);
  fAl->SetParameters(mAl,mMgL);
  fAl2->SetParameters(mAl,mMgL+140.);

  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  fL->SetTitle("momentum correlation ; pe' [MeV] ; pK [MeV]");
  fL->Draw();
  fL2->Draw("same");
  fS->Draw("same");
  fnnL->Draw("same");
  fAl->Draw("same");
  fAl2->Draw("same");

  double B_nnL=B(mT,mnnL,1.8e3,2.2e3);
  double B_MgL=B(mAl,mMgL,1.8e3,2.2e3);

  cout<<"B_nnL "<<B_nnL<<" B_MgL "<<B_MgL<<endl;

  TF1* fm_nnL=new TF1("fm_nnL","mmass",1600.,2000.,3);
  TF1* fm_MgL=new TF1("fm_MgL","mmass",1600.,2000.,3);
  TF1* fm_MgL2=new TF1("fm_MgL2","mmass",1600.,2000.,3);
  //  TF1* fm_nnL=new TF1("fm_nnL","mmass",2000.,2400.,3);
  //  TF1* fm_MgL=new TF1("fm_MgL","mmass",2000.,2400.,3);
  //   TF1* fm_MgL2=new TF1("fm_MgL","mmass",2000.,2400.,3);
  fm_nnL->SetParameters(mT ,mnnL,2218.);
  fm_MgL->SetParameters(mAl,mMgL,2218.);
  fm_MgL2->SetParameters(mAl,mMgL,2100.);

  fm_nnL->SetLineColor(2);
  fm_MgL->SetLineColor(4);
  fm_MgL2->SetLineColor(5);
  TCanvas* c1=new TCanvas("c1","c1");
  c1->cd();
  fm_nnL->Draw();  
  fm_MgL->Draw("same");
  fm_MgL2->Draw("same");

}
