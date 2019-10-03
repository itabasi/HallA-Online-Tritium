

void vdc_xtcurve(){

  string fname="/data1/root/tritium_111500.root";
  //  TFile* ofs=new TFile(fname.c_str());
  TChain* T=new TChain("T");
    T->Add(fname.c_str());
  const int nmax=200;
  double DR_evtypebits,R_vdc_u1_rawtime[nmax],R_vdc_u1_wire[nmax];
  int Ndata_R_vdc_u1_wire;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("DR.evtypebits"        ,1);    T->SetBranchAddress("DR.evtypebits"        ,&DR_evtypebits);
  T->SetBranchStatus("R.vdc.u1.rawtime"     ,1);    T->SetBranchAddress("R.vdc.u1.rawtime"     ,R_vdc_u1_rawtime);
  T->SetBranchStatus("R.vdc.u1.wire"        ,1);    T->SetBranchAddress("R.vdc.u1.wire"        ,R_vdc_u1_wire);
  T->SetBranchStatus("Ndata.R.vdc.u1.wire"  ,1);    T->SetBranchAddress("Ndata.R.vdc.u1.wire"  ,&Ndata_R_vdc_u1_wire);



  double conv_ns = 0.5; //[ch] -> [ns]
  //  conv_ns=1.0;
  double v=50.4e-3;// [mm/ns]
  double conv_mm  = conv_ns*v;// [ch]->[mm]
  
  double min_time = 0.0;
  double max_time = 800.0*conv_ns;
  int bin_time =200;
  double min_dist=0.0;
  double max_dist=20.;
  int bin_dist=100;
  
  TH2F*hu1=new TH2F("hu1","Drift velocity [mm/ns];TDC [ns] ; Drift dist [mm]",bin_time,min_time,max_time,bin_dist,min_dist,max_dist);

  int ENum=T->GetEntries();
  cout<<"Events "<<ENum<<endl;

  for(int i=0;i<ENum;i++){

    T->GetEntry(i);
    //    cout<<"i "<<i<<" trig "<<DR_evtypebits<<" Ndata "<<Ndata_R_vdc_u1_wire<<" wire "<<R_vdc_u1_wire[2]<<endl;
    if(DR_evtypebits==16 && Ndata_R_vdc_u1_wire==5 && (R_vdc_u1_wire[0]==184 || R_vdc_u1_wire[1]==184 ||
						       R_vdc_u1_wire[2]==184 || R_vdc_u1_wire[3]==184 ||R_vdc_u1_wire[4]==184)){
    ///        if(DR_evtypebits==16 && Ndata_R_vdc_u1_wire==5){
      for(int j=0;j<3;j++){
	double jj=(double)j;
	double dist=fabs( (fabs(jj)*( fabs(R_vdc_u1_rawtime[1] - R_vdc_u1_rawtime[0])) + (R_vdc_u1_rawtime[1] - R_vdc_u1_rawtime[3])/2.0)*conv_mm);
	double time=-(R_vdc_u1_rawtime[2+j]-2650.)*conv_ns;
	//	cout<<"dist "<<dist<<"time "<<time<<" rawtime "<<R_vdc_u1_rawtime[2+j]-2650.<<endl;
	hu1->Fill(time,dist);

      }
    }
  }


  TCanvas* c0=new TCanvas("c0","c0");
    c0->cd();
    hu1->Draw("colz");
    //  hu1->Draw("");
}
