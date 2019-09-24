// Auther Itabashi 2019/9/20
// Checking Trigger corretion in U1 VDC


extern double vdc_off(int wire);

void vdc_trig_corr(){

  string fname="/data1/root/tritium_111160.root";
  //  TFile* f1=new TFile(fname.c_str());
  TChain* T=new TChain("T");
  T->Add(fname.c_str());
  

  //===============================//
  //======== Set Branch ===========//
  //===============================//
  
  int nmax=100;
  double DR_evtype,R_s2_time[nmax],L_s2_time[nmax],R_vdc_u1_wire[nmax],R_vdc_u1_ritme[nmax],Ndata_R_vdc_u1_rtime,R_vdc_u1_rtime[nmax],R_vdc_u1_time[nmax];
  double R_s2_rt[nmax],R_s2_rt_c[nmax],R_s2_t_pads[nmax],L_s2_rt[nmax],L_s2_rt_c[nmax],L_s2_t_pads[nmax];
  int Ndata_R_vdc_u1_wire,Ndata_R_vdc_rtime;
  T->SetBranchStatus("*",0);  
  T->SetBranchStatus("DR.evtype",1);                 T->SetBranchAddress("DR.evtypebits",&DR_evtype);
  T->SetBranchStatus("R.s2.time"            ,1);     T->SetBranchAddress("R.s2.time"            , R_s2_time           );
  T->SetBranchStatus("L.s2.time"            ,1);  T->SetBranchAddress("L.s2.time"            , L_s2_time           );
  T->SetBranchStatus("R.s2.time"            ,1);     T->SetBranchAddress("R.s2.time"            , R_s2_time           );
  T->SetBranchStatus("L.s2.time"            ,1);  T->SetBranchAddress("L.s2.time"            , L_s2_time           );
  T->SetBranchStatus("R.s2.rt"              ,1);  T->SetBranchAddress("R.s2.rt"              , R_s2_rt             );
  T->SetBranchStatus("R.s2.rt_c"            ,1);  T->SetBranchAddress("R.s2.rt_c"            , R_s2_rt_c           );
  T->SetBranchStatus("R.s2.t_pads"          ,1);  T->SetBranchAddress("R.s2.t_pads"          , R_s2_t_pads         );
  T->SetBranchStatus("L.s2.rt"              ,1);  T->SetBranchAddress("L.s2.rt"              , L_s2_rt             );
  T->SetBranchStatus("L.s2.rt_c"            ,1);  T->SetBranchAddress("L.s2.rt_c"            , L_s2_rt_c           );
  T->SetBranchStatus("L.s2.t_pads"          ,1);  T->SetBranchAddress("L.s2.t_pads"          , L_s2_t_pads         );
  T->SetBranchStatus("R.vdc.u1.wire"          ,1);            T->SetBranchAddress("R.vdc.u1.wire"          ,R_vdc_u1_wire);
  T->SetBranchStatus("Ndata.R.vdc.u1.wire"          ,1);   T->SetBranchAddress("Ndata.R.vdc.u1.wire"          ,&Ndata_R_vdc_u1_wire);         
  // T->SetBranchStatus("R.vdc.u2.wire"          ,1);            T->SetBranchAddress("R.vdc.u2.wire"          ,R_vdc_u2_wire);
  // T->SetBranchStatus("R.vdc.v1.wire"          ,1);            T->SetBranchAddress("R.vdc.v1.wire"          ,R_vdc_v1_wire);
  // T->SetBranchStatus("R.vdc.v2.wire"          ,1);            T->SetBranchAddress("R.vdc.v2.wire"          ,R_vdc_v2_wire);
  T->SetBranchStatus("Ndata.R.vdc.u1.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.u1.rawtime"          ,&Ndata_R_vdc_u1_rtime);         
  // T->SetBranchStatus("Ndata.R.vdc.u2.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.u2.rawtime"          ,&Ndata_R_vdc_u2_rtime);         
  // T->SetBranchStatus("Ndata.R.vdc.v1.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.v1.rawtime"          ,&Ndata_R_vdc_v1_rtime);         
  // T->SetBranchStatus("Ndata.R.vdc.v2.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.v2.rawtime"          ,&Ndata_R_vdc_v2_rtime);
  T->SetBranchStatus("R.vdc.u1.rawtime"          ,1);         T->SetBranchAddress("R.vdc.u1.rawtime"          ,R_vdc_u1_rtime);         
  // T->SetBranchStatus("R.vdc.u2.rawtime"          ,1);         T->SetBranchAddress("R.vdc.u2.rawtime"          ,R_vdc_u2_rtime);         
  // T->SetBranchStatus("R.vdc.v1.rawtime"          ,1);         T->SetBranchAddress("R.vdc.v1.rawtime"          ,R_vdc_v1_rtime);         
  // T->SetBranchStatus("R.vdc.v2.rawtime"          ,1);         T->SetBranchAddress("R.vdc.v2.rawtime"          ,R_vdc_v2_rtime);         
  T->SetBranchStatus("R.vdc.u1.time"          ,1);         T->SetBranchAddress("R.vdc.u1.time"          ,R_vdc_u1_time);         
  // T->SetBranchStatus("R.vdc.u2.time"          ,1);         T->SetBranchAddress("R.vdc.u2.time"          ,R_vdc_u2_time);         
  // T->SetBranchStatus("R.vdc.v1.time"          ,1);         T->SetBranchAddress("R.vdc.v1.time"          ,R_vdc_v1_time);         
  // T->SetBranchStatus("R.vdc.v2.time"          ,1);         T->SetBranchAddress("R.vdc.v2.time"          ,R_vdc_v2_time);         

  double min_tdc=-0.5e-6; //[sec]
  double max_tdc= 0.5e-6; //[sec]
  int    bin_tdc= 1000;//(int)((max_tdc - min_tdc)/5.0e-10);
  TH1F* h_vdc_u1_t=new TH1F("h_vdc_u1_t","VDC U1 Time hist w/o Trigger time correction [sec]; TDC [sec] ; Counts ",bin_tdc,min_tdc,max_tdc);
  TH1F* h_vdc_u1_tc=new TH1F("h_vdc_u1_tc","VDC U1 Time hist [sec] w/ Trigger time correction [sec]; TDC [sec] ; Counts",bin_tdc,min_tdc,max_tdc);

  
  
  //================================//
  //======= Fill Event =============//
  //================================//

  int NEvent=T->GetEntries();
  cout<<"Events "<<NEvent<<endl;
  for(int i=0;i<NEvent;i++){

  //------ Initialization -----//
    for(int j=0;j<nmax;j++){
      R_s2_time[j]     = -1000.0;
      R_s2_rt_c[j]     = -1000.0;
      R_s2_t_pads[j]   = -1000.0;
      L_s2_rt_c[j]     = -1000.0;
      L_s2_time[j]     = -1000.0;
      L_s2_t_pads[j]   = -1000.0;
      R_vdc_u1_rtime[j]= -1000.0;
      R_vdc_u1_time[j] = -1000.0;
      R_vdc_u1_wire[j] = -1000.0;
    }
    Ndata_R_vdc_u1_rtime = -1000.0;
    Ndata_R_vdc_u1_wire  = -1000.0;
    
      
    T->GetEntry(i);

    
    double conv_tdc = 5.0e-10; // [ch] -> [sec]
    double T_trig   = L_s2_rt_c[(int)L_s2_t_pads[0]] - R_s2_rt_c[(int)R_s2_t_pads[0]];

    //    cout<<"DR_evtype "<<DR_evtype<<" Ndata_R_vdc_u1_wire "<<Ndata_R_vdc_u1_wire<<" T_trig "<<T_trig<<endl;
    //    if(DR_evtype==32 && Ndata_R_vdc_u1_wire==5 && R_vdc_u1_wire[0]==305){
    if(DR_evtype==32){
      double T_corr[Ndata_R_vdc_u1_wire];
      for(int j=0;j<Ndata_R_vdc_u1_wire;j++){
	T_corr[j]=0.0; // Initialization
	T_corr[j]=-conv_tdc*(R_vdc_u1_rtime[j]-vdc_off((int)R_vdc_u1_wire[j]))-T_trig;  
	h_vdc_u1_t ->Fill(R_vdc_u1_time[j]);
	h_vdc_u1_tc->Fill(T_corr[j]);
	
      } // end if
    } // end if
  }// End fill  


  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  h_vdc_u1_tc->SetLineColor(4);
  h_vdc_u1_tc->SetFillColor(4);
  h_vdc_u1_tc->SetFillStyle(3002);
  h_vdc_u1_tc->Draw();
  h_vdc_u1_t->Draw("same");

  

  
}




///========================================================//
//================ Function ==============================//
///========================================================//


double vdc_off(int wire){
  // Offset of each wires //
  double offset_u1[368]={
    2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 
    2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 
    2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 
    2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 2644.0, 
    2643.1, 2643.1, 2643.1, 2643.1, 2643.1, 2643.1, 2643.1, 2643.1, 
    2643.1, 2643.1, 2643.1, 2643.1, 2643.1, 2643.1, 2643.1, 2643.1, 
    2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 
    2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 
    2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 
    2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 2643.4, 
    2640.0, 2640.0, 2640.0, 2640.0, 2640.0, 2640.0, 2640.0, 2640.0, 
    2640.0, 2640.0, 2640.0, 2640.0, 2640.0, 2640.0, 2640.0, 2640.0, 
    2643.7, 2643.7, 2643.7, 2643.7, 2643.7, 2643.7, 2643.7, 2643.7, 
    2643.7, 2643.7, 2643.7, 2643.7, 2643.7, 2643.7, 2643.7, 2643.7, 
    2645.6, 2645.6, 2645.6, 2645.6, 2645.6, 2645.6, 2645.6, 2645.6, 
    2645.6, 2645.6, 2645.6, 2645.6, 2645.6, 2645.6, 2645.6, 2645.6, 
    2647.1, 2647.1, 2647.1, 2647.1, 2647.1, 2647.1, 2647.1, 2647.1, 
    2647.1, 2647.1, 2647.1, 2647.1, 2647.1, 2647.1, 2647.1, 2647.1, 
    2645.7, 2645.7, 2645.7, 2645.7, 2645.7, 2645.7, 2645.7, 2645.7, 
    2645.7, 2645.7, 2645.7, 2645.7, 2645.7, 2645.7, 2645.7, 2645.7, 
    2647.3, 2647.3, 2647.3, 2647.3, 2647.3, 2647.3, 2647.3, 2647.3, 
    2647.3, 2647.3, 2647.3, 2647.3, 2647.3, 2647.3, 2647.3, 2647.3, 
    2646.9, 2646.9, 2646.9, 2646.9, 2646.9, 2646.9, 2646.9, 2646.9, 
    2646.9, 2646.9, 2646.9, 2646.9, 2646.9, 2646.9, 2646.9, 2646.9, 
    2649.4, 2649.4, 2649.4, 2649.4, 2649.4, 2649.4, 2649.4, 2649.4, 
    2649.4, 2649.4, 2649.4, 2649.4, 2649.4, 2649.4, 2649.4, 2649.4, 
    2648.0, 2648.0, 2648.0, 2648.0, 2648.0, 2648.0, 2648.0, 2648.0, 
    2648.0, 2648.0, 2648.0, 2648.0, 2648.0, 2648.0, 2648.0, 2648.0, 
    2647.7, 2647.7, 2647.7, 2647.7, 2647.7, 2647.7, 2647.7, 2647.7, 
    2647.7, 2647.7, 2647.7, 2647.7, 2647.7, 2647.7, 2647.7, 2647.7, 
    2648.2, 2648.2, 2648.2, 2648.2, 2648.2, 2648.2, 2648.2, 2648.2, 
    2648.2, 2648.2, 2648.2, 2648.2, 2648.2, 2648.2, 2648.2, 2648.2, 
    2652.6, 2652.6, 2652.6, 2652.6, 2652.6, 2652.6, 2652.6, 2652.6, 
    2652.6, 2652.6, 2652.6, 2652.6, 2652.6, 2652.6, 2652.6, 2652.6, 
    2645.0, 2645.0, 2645.0, 2645.0, 2645.0, 2645.0, 2645.0, 2645.0, 
    2645.0, 2645.0, 2645.0, 2645.0, 2645.0, 2645.0, 2645.0, 2645.0, 
    2655.2, 2655.2, 2655.2, 2655.2, 2655.2, 2655.2, 2655.2, 2655.2, 
    2655.2, 2655.2, 2655.2, 2655.2, 2655.2, 2655.2, 2655.2, 2655.2, 
    2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 
    2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 
    2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 
    2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 
    2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 
    2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 2652.1, 
    2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 
    2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4, 2655.4};
  
  return offset_u1[wire];

}
