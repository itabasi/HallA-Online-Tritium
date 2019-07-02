

double T0time(char* ARM, int wire);

void VDC_tricorr(){

  string ifname="../run_list/nnlambda/Lambda_test.list";
  
  TChain* T=new TChain("T");

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> runname;
    T->Add(runname.c_str());
    //     cout<<buf<<endl;
  }

  int MAX=100;
  int R_Nu1_wire,R_Nu2_wire,R_Nv1_wire,R_Nv2_wire;
  int L_Nu1_wire,L_Nu2_wire,L_Nv1_wire,L_Nv2_wire;
  double R_u1_wire[MAX],R_u2_wire[MAX],R_v1_wire[MAX],R_v2_wire[MAX];
  double L_u1_wire[MAX],L_u2_wire[MAX],L_v1_wire[MAX],L_v2_wire[MAX];

  int R_Nu1_rawtime,R_Nu2_rawtime,R_Nv1_rawtime,R_Nv2_rawtime;
  int L_Nu1_rawtime,L_Nu2_rawtime,L_Nv1_rawtime,L_Nv2_rawtime;  
  double R_u1_rawtime[MAX],R_u2_rawtime[MAX],R_v1_rawtime[MAX],R_v2_rawtime[MAX];
  double L_u1_rawtime[MAX],L_u2_rawtime[MAX],L_v1_rawtime[MAX],L_v2_rawtime[MAX];

  int R_Nu1_time,R_Nu2_time,R_Nv1_time,R_Nv2_time;
  int L_Nu1_time,L_Nu2_time,L_Nv1_time,L_Nv2_time;  
  double R_u1_time[MAX],R_u2_time[MAX],R_v1_time[MAX],R_v2_time[MAX];
  double L_u1_time[MAX],L_u2_time[MAX],L_v1_time[MAX],L_v2_time[MAX];    
  double DR_evtype;
  double R_s2_rt_c[MAX],L_s2_rt_c[MAX];
  double Rs2trpad[MAX],Ls2trpad[MAX];
  double RF1[MAX],LF1[MAX];
  double R_s0_rt_c,L_s0_rt_c;

  
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("DR.evtypebits"          ,1);
  T->SetBranchAddress("DR.evtypebits"          ,&DR_evtype);

  T->SetBranchStatus("R.s0.rt_c"          ,1);
  T->SetBranchAddress("R.s0.rt_c"          , &R_s0_rt_c);  
  T->SetBranchStatus("L.s0.rt_c"          ,1);
  T->SetBranchAddress("L.s0.rt_c"          , &L_s0_rt_c);
  
  //===== RS2 =====//
  T->SetBranchStatus("R.s2.rt_c"          ,1);
  T->SetBranchAddress("R.s2.rt_c"          ,R_s2_rt_c);
  T->SetBranchStatus("R.s2.trpad",1);
  T->SetBranchAddress("R.s2.trpad",Rs2trpad);
  T->SetBranchStatus("RTDC.F1FirstHit",1);
  T->SetBranchAddress("RTDC.F1FirstHit",RF1);
  //==== Ls2 =====//
  T->SetBranchStatus("L.s2.rt_c"          ,1);
  T->SetBranchAddress("L.s2.rt_c"          ,L_s2_rt_c);  
  T->SetBranchStatus("L.s2.trpad",1);
  T->SetBranchAddress("L.s2.trpad",Ls2trpad);
  T->SetBranchStatus("LTDC.F1FirstHit",1);
  T->SetBranchAddress("LTDC.F1FirstHit",LF1);

  
  //========= RVDC ==========//  
  T->SetBranchStatus("Ndata.R.vdc.u1.wire"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.u1.wire"          ,&R_Nu1_wire);  
  T->SetBranchStatus("Ndata.R.vdc.u2.wire"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.u2.wire"          ,&R_Nu2_wire);    
  T->SetBranchStatus("Ndata.R.vdc.v1.wire"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.v1.wire"          ,&R_Nv1_wire);  
  T->SetBranchStatus("Ndata.R.vdc.v2.wire"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.v2.wire"          ,&R_Nv2_wire);  
  T->SetBranchStatus("R.vdc.u1.wire"          ,1);
  T->SetBranchAddress("R.vdc.u1.wire"                   ,R_u1_wire     );
  T->SetBranchStatus("R.vdc.u2.wire"          ,1);
  T->SetBranchAddress("R.vdc.u2.wire"                   ,R_u2_wire     );
  T->SetBranchStatus("R.vdc.v1.wire"          ,1);
  T->SetBranchAddress("R.vdc.v1.wire"                   ,R_v1_wire     );
  T->SetBranchStatus("R.vdc.v2.wire"          ,1);
  T->SetBranchAddress("R.vdc.v2.wire"                   ,R_v2_wire     );

  T->SetBranchStatus("Ndata.R.vdc.u1.rawtime"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.u1.rawtime"          ,&R_Nu1_rawtime);
  T->SetBranchStatus("Ndata.R.vdc.u2.rawtime"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.u2.rawtime"          ,&R_Nu2_rawtime);
  T->SetBranchStatus("Ndata.R.vdc.v1.rawtime"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.v1.rawtime"          ,&R_Nv1_rawtime);
  T->SetBranchStatus("Ndata.R.vdc.v2.rawtime"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.v2.rawtime"          ,&R_Nv2_rawtime);
  T->SetBranchStatus("R.vdc.u1.rawtime"          ,1);
  T->SetBranchAddress("R.vdc.u1.rawtime"                   ,R_u1_rawtime     );
  T->SetBranchStatus("R.vdc.u2.rawtime"          ,1);
  T->SetBranchAddress("R.vdc.u2.rawtime"                   ,R_u2_rawtime     );
  T->SetBranchStatus("R.vdc.v1.rawtime"          ,1);
  T->SetBranchAddress("R.vdc.v1.rawtime"                   ,R_v1_rawtime     );
  T->SetBranchStatus("R.vdc.v2.rawtime"          ,1);
  T->SetBranchAddress("R.vdc.v2.rawtime"                   ,R_v2_rawtime     );

  T->SetBranchStatus("Ndata.R.vdc.u1.time"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.u1.time"          ,&R_Nu1_time);
  T->SetBranchStatus("Ndata.R.vdc.u2.time"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.u2.time"          ,&R_Nu2_time);
  T->SetBranchStatus("Ndata.R.vdc.v1.time"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.v1.time"          ,&R_Nv1_time);
  T->SetBranchStatus("Ndata.R.vdc.v2.time"          ,1);
  T->SetBranchAddress("Ndata.R.vdc.v2.time"          ,&R_Nv2_time);
  T->SetBranchStatus("R.vdc.u1.time"          ,1);
  T->SetBranchAddress("R.vdc.u1.time"                   ,R_u1_time     );
  T->SetBranchStatus("R.vdc.u2.time"          ,1);
  T->SetBranchAddress("R.vdc.u2.time"                   ,R_u2_time     );
  T->SetBranchStatus("R.vdc.v1.time"          ,1);
  T->SetBranchAddress("R.vdc.v1.time"                   ,R_v1_time     );
  T->SetBranchStatus("R.vdc.v2.time"          ,1);
  T->SetBranchAddress("R.vdc.v2.time"                   ,R_v2_time     );
  






  int ENum=  T->GetEntries();
  cout<<"Events : "<<ENum <<endl;
  int k=0;
  int ev=0;
  
  bool trig=false;
  bool wire_u1=false;
  bool wire_u2=false;
  bool wire_v1=false;
  bool wire_v2=false;
  double T0;
  double time_u1;
  double time_u2;
  double time_v1;
  double time_v2;
  double fTDCRes=0.5; //[ns/ch]

  double min_u1=-1000;
  double max_u1=5000;
  int bin_u1=1000;
  TH1F* hu1_time4=new TH1F("hu1_time4","RVDC TDC [ns] wire 305-320 Trigger 4",bin_u1,min_u1,max_u1);
  TH1F* hu1_time4_c=new TH1F("hu1_time4_c","RVDC TDC [ns] wire w/ correction 305-320 Trigger 4",bin_u1,min_u1,max_u1);    
  TH1F* hu1_time=new TH1F("hu1_time","RVDC TDC [ns] wire 305-320",bin_u1,min_u1,max_u1);
  TH1F* hu1_time_c=new TH1F("hu1_time_c","RVDC TDC w/ correction [ns] wire 305-320",bin_u1,min_u1,max_u1);
  for(int i=0;i<ENum;i++){

   if(k==ev*10000){
     cout<<"Fill Event: "<<k<<"/"<<ENum<<endl; 
     ev=ev+1;
   }

    
    trig=false;
    wire_u1=false;
    wire_u2=false;
    wire_v1=false;
    wire_v2=false;

    
    T->GetEntry(i);
    if(DR_evtype==32)trig=true;
    if(R_u1_wire[0]>=305 && R_u1_wire[0]<=320)wire_u1=true;
    if(R_u2_wire[0]>=305 && R_u2_wire[0]<=320)wire_u2=true;
    if(R_v1_wire[0]>=305 && R_v1_wire[0]<=320)wire_v1=true;
    if(R_v2_wire[0]>=305 && R_v2_wire[0]<=320)wire_v2=true;    

    int Rs2pads=(int)Rs2trpad[0];
    int Ls2pads=(int)Ls2trpad[0];


    
     T0=-(R_s2_rt_c[Rs2pads]-L_s2_rt_c[Ls2pads]) *1.0e9;// [ns] S2 offset
     T0=-(R_s0_rt_c-L_s0_rt_c) *1.0e9; // S0 offset   
    if(trig && wire_u1){
      for(int j=0;j<R_Nu1_rawtime;j++){
	time_u1=fTDCRes*( 2528.8 - R_u1_rawtime[j] ) - T0;
	hu1_time_c->Fill(time_u1);
	hu1_time->Fill(R_u1_time[j]*1.0e9);
      } 
    } //end U1 hist
    
    if(wire_u1 && DR_evtype==16){
      for(int j=0;j<R_Nu1_rawtime;j++){
	hu1_time4->Fill(R_u1_time[j]*1.0e9);
	time_u1=fTDCRes*( 2528.8 - R_u1_rawtime[j] ) - T0;
	hu1_time4_c->Fill(time_u1);
      }
    };

  }//END Fill



  
  TCanvas* c0=new TCanvas("c0","c0");

  c0->cd();
  hu1_time->SetLineColor(1);
  hu1_time_c->SetLineColor(4);
  hu1_time_c->SetFillColor(4);
  hu1_time_c->SetLineStyle(3004);  

  hu1_time_c->Draw();
  hu1_time->Draw("same");  
 TCanvas* c1=new TCanvas("c1","c1");
 c1->cd();
 hu1_time4->SetLineColor(1);
 hu1_time4_c->SetLineColor(4);
 hu1_time4_c->SetFillColor(4);
 hu1_time4_c->SetLineStyle(3004);   
 hu1_time4->Draw();
 hu1_time4_c->Draw("same");



}


//========================================================//
//=========== Function ===================================//
//========================================================//


