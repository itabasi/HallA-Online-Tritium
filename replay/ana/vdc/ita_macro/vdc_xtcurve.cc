////////////////////////////////
// Author Itabshi 2019/10
//Drift velocity check macro
//////////////////////////////

extern double Calc_ddist(double dist, double time, double theta, bool RVDC);

void vdc_xtcurve(){

  string fname="/data1/root/tritium_111500.root";
  //  TFile* ofs=new TFile(fname.c_str());
  TChain* T=new TChain("T");
  T->Add(fname.c_str());
  const int nmax=200;
  double DR_evtypebits,R_vdc_u1_rawtime[nmax],R_vdc_u1_wire[nmax], R_vdc_u1_time[nmax];
  int Ndata_R_vdc_u1_wire;
  double R_tr_th[nmax],L_tr_th[nmax];
  double R_vdc_u1_slope[nmax],R_vdc_u1_ddist[nmax],R_vdc_u1_dist[nmax],R_vdc_u1_trdist[nmax];
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("DR.evtypebits"        ,1);    T->SetBranchAddress("DR.evtypebits"        ,&DR_evtypebits);
  T->SetBranchStatus("R.vdc.u1.rawtime"     ,1);    T->SetBranchAddress("R.vdc.u1.rawtime"     ,R_vdc_u1_rawtime);
  T->SetBranchStatus("R.vdc.u1.time"     ,1);       T->SetBranchAddress("R.vdc.u1.time"     ,R_vdc_u1_time);  
  T->SetBranchStatus("R.vdc.u1.wire"        ,1);    T->SetBranchAddress("R.vdc.u1.wire"        ,R_vdc_u1_wire);
  T->SetBranchStatus("Ndata.R.vdc.u1.wire"  ,1);    T->SetBranchAddress("Ndata.R.vdc.u1.wire"  ,&Ndata_R_vdc_u1_wire);
  T->SetBranchStatus("R.vdc.u1.slope"  ,1);             T->SetBranchAddress("R.vdc.u1.slope"  ,R_vdc_u1_slope);
  T->SetBranchStatus("R.vdc.u1.ddist"  ,1);             T->SetBranchAddress("R.vdc.u1.ddist"  ,R_vdc_u1_ddist);
  T->SetBranchStatus("R.vdc.u1.dist"  ,1);             T->SetBranchAddress("R.vdc.u1.dist"  ,R_vdc_u1_dist);
  T->SetBranchStatus("R.vdc.u1.trdist"  ,1);             T->SetBranchAddress("R.vdc.u1.trdist"  ,R_vdc_u1_trdist);    
  T->SetBranchStatus("R.tr.d_th"  ,1);             T->SetBranchAddress("R.tr.d_th"  ,R_tr_th);
  T->SetBranchStatus("L.tr.d_th"  ,1);             T->SetBranchAddress("L.tr.d_th"  ,L_tr_th);  


  double conv_ns = 0.5; //[ch] -> [ns]
  double v=50.4e-3;// [mm/ns]
  double conv_mm  = conv_ns*v;// [ch]->[mm]
  
  double min_time = 0.0;
  double max_time = 800.0*conv_ns;
  int bin_time =200;
  double min_dist=0.0;
  double max_dist=20.;
  int bin_dist=100;
  
  TH2F*hu1=new TH2F("hu1","Drift velocity [mm/ns];TDC [ns] ; Drift dist [mm]",bin_time,min_time,max_time,bin_dist,min_dist,max_dist);
  TH2F*hu1_c=new TH2F("hu1_c","Drift velocity w/ Drift correction [mm/ns];TDC [ns] ; Drift dist [mm]",bin_time,min_time,max_time,bin_dist,min_dist,max_dist);  
  TH2F* hslope_dist=new TH2F("hslope_dist","Slope vs Dist correlation ; slope [tan] ; distance [mm]",1000,0.0,5.0,1000,0,20);

  TH1F*hres=new TH1F ("hres","Residual ;Residual [mm]",1000,-10,10);
		      //,-1.0e-6,1.0e-6);

  TH2F*hdd=new TH2F("hdd","Drift velocity [mm/ns];Dist [mm] ; Drift dist [mm]",1000,0.0,20,bin_dist,min_dist,max_dist);

  
  int ENum=T->GetEntries();
  cout<<"Events "<<ENum<<endl;

  for(int i=0;i<ENum;i++){

    T->GetEntry(i);
    if(DR_evtypebits==16 && Ndata_R_vdc_u1_wire==5 && (R_vdc_u1_wire[0]==184 || R_vdc_u1_wire[1]==184 ||
						       R_vdc_u1_wire[2]==184 || R_vdc_u1_wire[3]==184 ||R_vdc_u1_wire[4]==184)){
    ///        if(DR_evtypebits==16 && Ndata_R_vdc_u1_wire==5){
      for(int j=0;j<3;j++){
	double jj=(double)j;
	double dist=fabs( (fabs(jj)*( fabs(R_vdc_u1_rawtime[1] - R_vdc_u1_rawtime[0])) + (R_vdc_u1_rawtime[1] - R_vdc_u1_rawtime[3])/2.0)*conv_mm);
	double time=-(R_vdc_u1_rawtime[2+j]-2650.)*conv_ns;
	hu1->Fill(time,dist);
	//	double	ddist=Calc_ddist(dist/1000.,time,R_tr_th[j],true);
	//	double	ddist=Calc_ddist(dist,R_vdc_u1_time[j+2],R_tr_th[j+2],true);
	//	double	ddist=Calc_ddist(dist,time,R_tr_th[0],true);
	double	ddist=Calc_ddist(dist,time,R_vdc_u1_slope[0],true);		
	//	 hu1_c->Fill(time,ddist);
	hu1_c->Fill(time,R_vdc_u1_dist[j+2]*1000.);
	 //	 cout<<"dist "<<dist<<" calc_dist "<<ddist<<" Driftdist "<<R_vdc_u1_dist[j+2]*1000.<<endl;
	//	hu1_c->Fill(time,R_vdc_u1_dist[2+j]*1000.);	
	 hslope_dist->Fill(R_vdc_u1_slope[0],dist);
	 hres->Fill(dist-R_vdc_u1_trdist[j+2]*1000.);
	 hdd->Fill(dist,R_vdc_u1_trdist[j+2]*1000.);
	 //	 hres->Fill(R_vdc_u1_dist[j+2]*1000.-R_vdc_u1_trdist[j+2]*1000.);
	 //	 hdd->Fill(R_vdc_u1_dist[j+2]*1000.,R_vdc_u1_trdist[j+2]*1000.);
      }
    }
  }


  TCanvas* c0=new TCanvas("c0","c0");
  c0->Divide(2,1);
  c0->cd(1);
  hu1->Draw("colz");
  c0->cd(2);
  hu1_c->Draw("colz");  

  TCanvas*c1=new TCanvas("c1","c1");
  c1->cd();
  hslope_dist->Draw("colz");

  TCanvas*c2=new TCanvas("c2","c2");
  c2->cd();
  hres->Draw();

  TCanvas*c3=new TCanvas("c3","c3");
  c3->cd();
  hdd->Draw("colz");
  
}


double Calc_ddist(double dist, double time, double theta, bool RVDC){


  dist=dist/1000.; // mm -> m
  
  Double_t a1 = 0.0, a2 = 0.0;
  theta = 1.0 / theta;


  double fA1tdcCor_R[4]={    2.12e-03,   0.00e+00,   0.00e+00,   0.00e+00};
  double fA2tdcCor_R[4]={   -4.20e-04,   1.30e-03,   1.06e-04,   0.00e+00};
  double fA1tdcCor_L[4]={    2.12e-03,   0.00e+00,   0.00e+00,   0.00e+00};
  double fA2tdcCor_L[4]={   -4.20e-04,   1.30e-03,   1.06e-04,   0.00e+00};
  double fA1tdcCor[4];
  double fA2tdcCor[4];


  
  for(int i=0;i<4;i++){
    if(RVDC){
      fA1tdcCor[i]=fA1tdcCor_R[i];
      fA2tdcCor[i]=fA2tdcCor_R[i];
    }else{
      fA1tdcCor[i]=fA1tdcCor_L[i];
      fA2tdcCor[i]=fA2tdcCor_L[i];
    }
  }


  
  for (Int_t i = 3; i >= 1; i--) {
    a1 = theta * (a1 + fA1tdcCor[i]);
    a2 = theta * (a2 + fA2tdcCor[i]);
  }
  a1 += fA1tdcCor[0];
  a2 += fA2tdcCor[0];


  
  double fDriftVel=50000.;// [m/s]
  double ddist = fDriftVel * time/1.0e9;
  //  cout<<"dist "<<dist<<" ddist"<<ddist<<endl;
  //dist=ddist;

   
  if (dist < 0) {
    // something screwy is going on
  } else if (dist < a1 ) {
    //    distK= fDriftVel * time * (1 + 1 / (a1/a2 + 1));
    dist *= ( 1 + a2 / a1);
  }  else {
    dist +=  a2;
    }

  dist=dist*1000.; //m -> mm

  return dist;


 
}
