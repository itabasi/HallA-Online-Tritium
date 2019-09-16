

void check(){


  // string main = "test";
  //string main = "ang_RHRS_sieve_init";
 string main = "angcalib_4th_0915_4_0";
  string end_root = ".root";
  string end_dat = ".dat";
  string ifname="../rootfiles/angcalib/" + main + ".root";
  string ofname="./param/" + main + ".dat";
  TH2D* h2[10];
  TFile* f = new TFile(ifname.c_str() );
  for(int i=0;i<10;i++)h2[i] =(TH2D*)f->Get(Form("hang_%d",i));

  
  const double step = 0.492 * 2.54;
  const int nrow = 11; // the number of row in SS pattern
  const int ncol = 8;  // the number of column in SS pattern
  const int nsshole = nrow*ncol; // the number of holes to consider 
  const int nfoil =10;
  const double hrs_ang = 13.2 * 3.14159 / 180.;
  const double l0 = 1.03;
  bool RHRSTrue=false;
  RHRSTrue=true;

  double fcent_real[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
			    0.000, 0.025, 0.050, 0.100, 0.125};
  TMarker* mark[nsshole][nfoil];
  TMarker* mark_real[nfoil][nsshole];  
  double w[nfoil][nsshole];
  double ssx_off[nfoil][nsshole];
  double ssy_off[nfoil][nsshole];
  bool flag[nfoil][nsshole];
  double refx[nsshole],refy[nsshole],refx_real[nfoil][nsshole],refy_real[nfoil][nsshole];
  double ssy_cent_real[nsshole],ssx_cent_real[nsshole];
  double th[100][10],ph[100][10];
  int nhole=0;
   for(int i=0; i<ncol ; i++){
    for(int j=0; j<nrow; j++){
      ssy_cent_real[i] = -3.0*step + step*i;
      if(j%2==0)ssy_cent_real[i] = ssy_cent_real[i] - step/2.0;
      ssx_cent_real[j] = 5.0*step - step*j;
      refx[nhole] = ssx_cent_real[j];
      refy[nhole] = ssy_cent_real[i];
      nhole++;
      
    }
   }

   double l[10],dth[10],lx;

   for(int i=0;i<nfoil;i++){

    l[i] = 0;
    l[i]=(l0-fcent_real[i]/cos(hrs_ang))*100.;    
    dth[i] = asin(l0*sin(hrs_ang)/(l0*cos(hrs_ang) -fcent_real[i]));
    int  n=0;
    for(int k=0; k<ncol ; k++){
      for(int j=0; j<nrow; j++){
	
      if(RHRSTrue==0) ph[n][i]= -tan( refy[n]*cos(dth[i])/(  refy[n]*sin(dth[i]) + l[i] )  );
      if(RHRSTrue)    ph[n][i]= -tan( refy[n]*cos(dth[i])/(- refy[n]*sin(dth[i]) + l[i] )  );
      if(refy[n]>0) lx=sqrt(pow(l[i],2.0) + pow(refy[n],2.0) + 2.0*l[i]*refy[n]*sin(dth[i]));
      else          lx=sqrt(pow(l[i],2.0) + pow(refy[n],2.0) - 2.0*l[i]*refy[n]*sin(dth[i]));
      th[n][i]=-refx[n]/lx;
      
      mark[n][i] = new TMarker(ph[n][i],th[n][i],28);
      mark[n][i]->SetMarkerColor(1);
      
      //      cout<<"nhole "<<n<<" nfoil "<<i<<" th "<<th[n]<<" ph "<<ph[n]<<endl; 

      n++;
      }
    }

    if(RHRSTrue){
      mark[23][i]->SetMarkerColor(2);
      mark[38][i]->SetMarkerColor(2);
    }else{
      mark[45][i]->SetMarkerColor(2);
      mark[38][i]->SetMarkerColor(2);
    }

    
   }

   



  TCanvas* c1=new TCanvas("c1","c1");  
  c1->Divide(4,3);
  for(int i=0;i<nfoil;i++){
    c1->cd(i+1);
    h2[i]->Draw("colz");
  }


  
    for(int j=0;j<nsshole;j++){
      for(int i=0;i<nfoil;i++){
	c1->cd(i+1);
	mark[j][i]->Draw("same");
      }
    }



}
