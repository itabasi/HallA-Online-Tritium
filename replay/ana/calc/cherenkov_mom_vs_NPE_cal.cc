// calculation of Cherenkov light momentum vs photon electron noumber
// 2018 Mar. 04 By Itabashi


void cherenkov_mom_vs_NPE_cal(){
  // make TCanvas
  TCanvas *c1=new TCanvas("c1","c1",800,500);


  //=========================//
  // NPE formula             //
  //=========================//

  double pi,alp,z,lam,n,beta,vp,vpi,vk,x,c,gam,ef;

  // Constance value 

  pi=3.14; // ratio of the circle
  alp=1./137; // fine structual constant
  // c=3.0e+8; //photon velocity [m/s] 
  c=1.;

  // Given value



  n=1.055; // Aerogel reflactive index HRS A2 
  //   n=1.015; // Aerogel reflective index HRS A1
  //n=1.05; //Aerogel reflective index HKS AC 
  //  n=1.33;// Water reflactive index HKS WC
  //   n=1.00041; // Gas reflective index HKS Gass

  cout<<"Index is :";
  cin>>n;

    z=1.; // charge 
    x=0.01;  // thickness [m]
    ef=1.; // efficency of AC

    double mp,mpi,mk;   //mass [MeV/c^2]
    
    mp=938.272/pow(c,2);//proton mass [MeV/c^2]
    mpi=139.57/pow(c,2);//pion mass [MeV/c^2]
    mk=493.68/pow(c,2);//kaon mass [MeV/c^2]

   // wave length calculation //
    double wave_min,wave_max;

   wave_min=298e-9; //wave minimum [m]
   wave_max=650e-9; // wave maximum [m]

   lam=1/wave_min-1/wave_max;
   //  cout<<"wave  "<<lam<<endl;
   
   //  double gamma[100];
   //   double N[100];
   double Np[5000],Npi[5000],Nk[5000];
   double Ep[5000],Epi[5000],Ek[5000]; // Energy [MeV]
   double beta_p[5000],beta_pi[5000],beta_k[5000]; //beta
   double p[5000];
   double pmax,pmin;
    pmax=4000; // HKS (Kaon side) momentum max [MeV/c]
    pmin=1; // // min [MeV/c]
    double plot_width;
    plot_width=10; // x axisis wadth interval
  
    //== Photon electron number of the Cosmic ray ==//
    double inf ;
    inf=2*ef*x*pi*alp*pow(z,2)*lam*(1.-1./pow(n,2));
  // Np[i]=2*ef*x*pi*alp*pow(z,2)*lam*(1.-1./pow(beta_p[i]*n,2));
    //==== Coment Out of Parameters ========//
    cout<<endl<<"======= Given Parameters ============="<<
      endl<<endl;
    cout<<"reflective index n = "<<n<<endl;
    cout<<" wave min  "<<wave_min*1e+9<<" [nm]"<<endl;
    cout<<"wave max  "<<wave_max*1e9<<" [nm]"<<endl;
    cout<<"thickness "<<x*100<<" [cm]"<<endl;
    cout<<"PEN of cosmic ray "<<inf<<endl;

    //==============================================================//
    //===== for ====================================================//
    //==============================================================//

 for(int i=pmin/plot_width/c;i<pmax/plot_width/c; i++){
     // gamma[i]=1/sqrt(1-pow(beta[i],2));
     //   for(int j=0;j<pmax; j++){   }

  
   p[i]=i*plot_width; // [MeV/c]

    Ep[i]=sqrt(pow(mp,2)+pow(p[i],2)); //proton energy [MeV]
    Epi[i]=sqrt(pow(mpi,2)+pow(p[i],2)); // pion energy [MeV]
    Ek[i]=sqrt(pow(mk,2)+pow(p[i],2));  //kaon energy [MeV]

 
    beta_p[i]=p[i]/Ep[i]; // proton beta
    beta_k[i]=p[i]/Ek[i]; // kaon beta
    beta_pi[i]=p[i]/Epi[i]; // pio beta
  

  
  

      Np[i]=2*ef*x*pi*alp*pow(z,2)*lam*(1.-1./pow(beta_p[i]*n,2));
     Npi[i]=2*ef*x*pi*alp*pow(z,2)*lam*(1.-1./pow(beta_pi[i]*n,2));
     Nk[i]=2*ef*x*pi*alp*pow(z,2)*lam*(1.-1./pow(beta_k[i]*n,2));
          if(Np[i]<0)
       {Np[i]=0;
}
     if(Npi[i]<0){
       Npi[i]=0;
     }

     if(Nk[i]<0){
       Nk[i]=0;
     }
     


    /*
    cout<<"i  "<<i<<endl;
    cout<<"pi  "<<p[i]<<endl;
    cout<<"beta_pi  "<<beta_pi[i]<<endl;
    cout<<"beta_p   "<<beta_p[i]<<endl;
    cout<<"beta_k    "<<beta_k[i]<<endl;
    cout<<"n  ="<<n<<endl;
    cout<<"Ep   "<<Ep[i]<<endl;
    
    //  cout<<"Calulation result 1 =  "<<(1.-1./pow(beta_pi[i]*n,2))<<endl;
     cout<<"Cal 2  =  "<<2*pi*x*ef*alp*pow(z,2)*lam<<endl;
     cout<<"Npi  "<<Npi[i]<<endl;
     */
}   


 
   TGraph*gp=new TGraph(pmax/plot_width,p,Np);
   TGraph*gpi=new TGraph(pmax/plot_width,p,Npi); 
   TGraph*gk=new TGraph(pmax/plot_width,p,Nk);
    gpi->GetXaxis()->SetRange(0.,4000);
    gpi->GetYaxis()->SetRange(0.,60);
    gpi->GetXaxis()->SetLabelSize(0.05);
    gpi->GetYaxis()->SetLabelSize(0.05);
    gpi->GetYaxis()->SetTitleSize(0.06);
    gpi->GetXaxis()->SetTitleSize(0.06);
    gpi->SetTitle("Aerogel Num. of P.E /cm");
    gpi->GetXaxis()->SetTitle("Momentum  [MeV/c]");
    gpi->GetYaxis()->SetTitle("NPE [ /cm]");
    
    gp->SetLineColor(3);
    gpi->SetLineColor(2);
    gk->SetLineColor(4);
    gpi->SetLineWidth(3);
    gp->SetLineWidth(3);
    gk->SetLineWidth(3);
         gpi->Draw();
	 gp->Draw("same");
	 gk->Draw("same");


    

//=========================================//
//======== fitting =======================//
//=======================================//

	 double d;
	 x=d;
	 /*
	 TH1F *hpi=new TH1F("hpi","hist_pi",pmax-pmin,pmin,pmax);
	 TH1F *hp=new TH1F("hp","hist_p",pmax-pmin,pmin,pmax);
	 TH1F *hk=new TH1F("hk","hist_k",pmax-pmin,pmin,pmax);
	 */

	 //  2*ef*d*pi*alp*pow(z,2)*lam*(1-1/pow(x*n,2))
	 //TF1 *fpi=new TF1("fpi","[0]*(1-1/pow(x*[1],2))",2*ef*d*pi*alp*pow(z,2)*lam,n);
	 //	  fpi->Draw("same");


}
