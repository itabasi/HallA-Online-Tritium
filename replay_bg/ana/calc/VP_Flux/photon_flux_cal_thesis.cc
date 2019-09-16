double fvp(Double_t *x,Double_t *par);
double Np(Double_t *x,Double_t *par);

  
void photon_flux_cal_thesis(){

  double Ee=4300;
  double Ee_=2100;
  double theta_ee=13.2;//degree
  //  theta_ee=15.2;
  //  Ee_=2000;
  // test //
  //    Ee=4523;
  //   Ee_=2725;
  //  theta_ee=12.5;
  //    double theta_ee=12.5;//degree
  double Q;//Q^2=-q^2 virtual photn 4 vector 
  double pv;//virtual photon momentum
  double mp=938;// proton mass
  double q; //virtual photon 3 vector 
  double Ev;//virtual photon energy(Ee-Ee')
  double PI=3.14;
  double ep;//epsilon 
  double Flux; //virtual photon Flux
  double alpha=0.007297;
  q=sqrt(Ee*Ee+Ee_*Ee_-2*Ee*Ee_*cos(theta_ee/180.*PI));
  Q=sqrt(-(Ee-Ee_)*(Ee-Ee_)+q*q);
  //  Ev=(Ee-Ee_)+q*q/(2.*mp);
  Ev=(Ee-Ee_)-Q*Q/(2.*mp);
  ep=1./(1.+2.*q*q/(Q*Q)*pow(tan(theta_ee/180.*PI/2.),2));
  // ep=0.61;
  Flux=alpha/(2*PI*PI*Q*Q)*Ev/(1-ep)*Ee_/Ee;

  cout<<"theta ee rad :"<<theta_ee/180.*PI<<endl;
  
  // double Flux2=3.8e-6;
  double min_e,max_e;
  min_e=10.0;
  max_e=16.0;

  //rad //
  //  min_e=0.15;
  //  max_e=0.3;
  
  TF1* f_flux=new TF1("f_flux",fvp,min_e,max_e,3);
  f_flux->SetParameter(0,Ee);
  f_flux->SetParameter(1,Ee_);
  f_flux->SetParameter(2,mp);

   TF1* f_flux2=new TF1("f_flux2",fvp,min_e,max_e,3);
  f_flux2->SetParameter(0,Ee);
  f_flux2->SetParameter(1,Ee_+100);
  f_flux2->SetParameter(2,mp);

   TF1* f_flux3=new TF1("f_flux3",fvp,min_e,max_e,3);
  f_flux3->SetParameter(0,Ee);
  f_flux3->SetParameter(1,Ee_-100);
  f_flux3->SetParameter(2,mp);
  
  
  f_flux->SetTitle("Virtul Photon Flux; angle [deg]; /MeV/sr/electron");
  
  TF1* fnp=new TF1("fnp",Np,min_e,max_e,3);
  fnp->SetParameter(0,Ee);
  fnp->SetParameter(1,Ee_);
  fnp->SetParameter(2,mp);

 

 
  
  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  f_flux->Draw();
  f_flux2->Draw("same");
  f_flux3->Draw("same");

  cout<<"Ee: "<<Ee<<endl;
  cout<<"theta_ee [rad]: "<<theta_ee/180.*3.14<<endl;
  cout<<"Ee_: "<<Ee_<<endl;
  cout<<"q: "<<q<<endl;
  cout<<"Q: "<<Q<<endl;
  cout<<"omega (Ee-Ee_):"<<Ee-Ee_<<endl;
  cout<<"Ev:"<<Ev<<endl;
  cout<<"epsilon:"<<ep<<endl;




  // Cross Section //

  //   double NL=1600;
  //   double NS=396;
  //   double NL=355; 
  //   double eff_KL=0.43;
    double NL=440; 
    double eff_KL=0.5;    
   double NS=92;



  //---- Total Events -----//
  //  double  NL=1100.0; //[total events]
  //  double  NS=220.0;
  //  double eff_KL=1.0;  
  //----------------------//

  
  double eff_SR=0.17;
  double eff_KS=0.72;
  double Nt=4.26e22;
  double dOme=0.006;//[sr]
  double e=1.602e-19;//[C]
  //  double C=4.754;//[C]
  double C=2.38;//[C]  half data of hydrogen run
  double conv=1.0e-24;//[b/cm^2]
  double dEe_=Ee_*2*0.045;//HRS Acceptance 4.5 %
  double Flux_sum=Flux*dOme*dEe_;
  //  double Nvp=C/e*Flux_sum;
  double Nvp=C/e*Flux;    
  double CSL=NL/(Nt*Nvp*eff_SR*eff_KL*dOme*conv);
  double CSS=NS/(Nt*Nvp*eff_SR*eff_KS*dOme*conv);  

  cout<<"Virtual Photon Flux [/MeV/sr/electron]:"<<Flux<<endl;
  cout<<"Virtual Photon Flux Integral sum [/e] :"<<Flux_sum<<endl;
  //   cout<<"Virtual Photon Flux [/s]:"<<Flux<<endl;
  
  cout<<"# virtual photon:"<<Nvp<<endl;
  cout<<"# proton in targets: "<<Nt<<endl;
  cout<<"Kaon Survival Rate: "<<eff_SR<<endl;
  cout<<"Kaon AC cut rate : "<<eff_KL<<endl;
  cout<<"R-HRS Accpectance [sr]: "<<dOme<<endl;
  cout<<"# Lambda events: "<<NL<<endl;
  cout<<"Lambda Cross Sectioon [nb/sr]:"<<CSL*1e9<<endl;
  cout<<"Sigma Cross Sectioon [nb/sr]:"<<CSS*1e9<<endl;  
}


double fvp(Double_t *x,Double_t *par){


   
  double q=sqrt(par[0]*par[0]+par[1]*par[1]-2*par[0]*par[1]*cos(x[0]/180*3.14));
  double Q=sqrt(-(par[0]-par[1])*(par[0]-par[1])+pow(q,2));
  double Ev=(par[0]-par[1])-pow(Q,2)/(2*par[2]);
  double ep=1./(1.+2*pow(q/Q,2)*pow(tan(x[0]/180*3.14/2.),2));
  double vp=1./137./(2*pow(3.14*Q,2))*Ev/(1.-ep)*par[1]/par[0];  
  
  /*
  double q=sqrt(par[0]*par[0]+par[1]*par[1]-2*par[0]*par[1]*cos(x[0]));
  double Q=sqrt(-(par[0]-par[1])*(par[0]-par[1])+pow(q,2));
  double Ev=(par[0]-par[1])-pow(Q,2)/(2*par[2]);
   double ep=1./(1.+2*pow(q/Q,2)*pow(tan(x[0]/2.),2));
  */
   //   double ep=1./(1.+2*pow(q/Q,2)*pow(tan(x[0]),2));
  //  double vp=1./137./(2*pow(3.14*Q,2))*Ev/(1.-ep)*par[1]/par[0];  
  
  return(vp);



}



double Np(Double_t *x,Double_t *par){

  
  double q=sqrt(par[0]*par[0]+par[1]*par[1]-2*par[0]*par[1]*cos(x[0]/180*3.14));
  double Q=sqrt(-(par[0]-par[1])*(par[0]-par[1])+pow(q,2));
  double Ev=(par[0]-par[1])-pow(q,2)/(2*Q);
  double ep=1./(1.+2*pow(q/Q,2)*pow(tan(x[0]/180*3.14/2.),2));
  double vp=1./137./(2*pow(3.14*Q,2))*Ev/(1.-ep)*par[1]/par[0];  
  double np=vp*2*3.14*sin(x[0]/180*3.14);
  return(np);



}
