#include <iostream>
#include <fstream>
//const double MK = 0.493677;            // charged Kaon mass (GeV/c2)
extern double s2f1_off(int i,string ARM, int KINE);

void F1Shift_corr_whist(){


  TChain* oldtree=new TChain("T");
  
  //  string ifname="../run_list/nnlambda/nnL_small_Ole2.list";
  string ifname="../run_list/nnlambda/test_nnL2.list";
  ifname="../run_list/nnlambda/test_nnL4.list";
  ifname="../run_list/nnlambda/nnL_small_Ole2.list";
  string oname;
  string iname;
  //  cout<<"input file name (../run_list/nnlambda/ ):";
  //  cin>> iname;
  cout<<"input file name :"<<ifname<<endl;
  cout<<"output name : ";
  cin >>oname;

  //  ifname="../run_list/nnlambda/" + oname +".list";
  int ENum;
  //  string ifname="../run_list/nnlambda/test_Ole2.list";
  //  string ifname="../run_list/nnlambda/nnL_small_Ole2.list";
  //    string ifname="../run_list/nnlambda/test_nnL2.list";
  //  string ifname="../run_list/nnlambda/Lambda_small_OleH1.list";
  //string ifname="../run_list/nnlambda/test_OleH1.list";
  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  
  if(!ifp){cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  TChain *T = new TChain("T");
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
     if( ifp.eof() ) break;
    stringstream sbuf(buf);
    sbuf >> runname;
    T->Add(runname.c_str());

  }



  TCanvas*c[17];
  for(int i=0;i<17;i++){
    c[i]=new TCanvas(Form("c%d",i),Form("c%d",i));
    c[i]->Divide(4,2);
  }
  c[0]->cd(1);

  
  ENum =T->GetEntries();
  cout<<"Get Entries: "<<ENum<<endl;

  TH1D* hRF1_t[16];
  TH1D* hRF1_b[16];
  TH1D* hLF1_t[16];
  TH1D* hLF1_b[16];
  TH1D* hRF1_ts[16];
  TH1D* hRF1_bs[16];
  TH1D* hLF1_ts[16];
  TH1D* hLF1_bs[16];  
  TH1D* hRF1_tc[16];
  TH1D* hRF1_bc[16];
  TH1D* hLF1_tc[16];
  TH1D* hLF1_bc[16];  
  
  TF1* fRF1_t[16];
  TF1* fRF1_b[16];
  TF1* fLF1_t[16];
  TF1* fLF1_b[16];

  TF1* fRF1_ts[16];
  TF1* fRF1_bs[16];
  TF1* fLF1_ts[16];
  TF1* fLF1_bs[16];
  

  int bin_rtdc,bin_ltdc;
  double min_rtdc,max_rtdc,min_ltdc,max_ltdc;
  int bin_rtdc_s,bin_ltdc_s;
  double min_rtdc_s,max_rtdc_s,min_ltdc_s,max_ltdc_s;  

  min_rtdc=-8200;max_rtdc=-7500;
  min_rtdc_s=57000;max_rtdc_s=57500;
  bin_rtdc=800;bin_rtdc_s=500;
  min_ltdc=-13000.;max_ltdc=-11000.;
  min_ltdc_s=52000.;max_ltdc_s=54000.;
  bin_ltdc=1500;bin_ltdc_s=2000;

  for(int i=0;i<16;i++){
    hRF1_t[i] = new TH1D(Form("hRF1_t_%d",i),Form("hRF1_t_%d",i),bin_rtdc,min_rtdc,max_rtdc);
    hRF1_tc[i] = new TH1D(Form("hRF1_tc_%d",i),Form("hRF1_tc_%d",i),bin_rtdc,min_rtdc,max_rtdc);
    hRF1_ts[i]= new TH1D(Form("hRF1_ts_%d",i),Form("hRF1_ts_%d",i),bin_rtdc_s,min_rtdc_s,max_rtdc_s);
    hRF1_b[i] = new TH1D(Form("hRF1_b_%d",i),Form("hRF1_b_%d",i),bin_rtdc,min_rtdc,max_rtdc);
    hRF1_bc[i] = new TH1D(Form("hRF1_bc_%d",i),Form("hRF1_bc_%d",i),bin_rtdc,min_rtdc,max_rtdc);
    hRF1_bs[i]= new TH1D(Form("hRF1_bs_%d",i),Form("hRF1_bs_%d",i),bin_rtdc_s,min_rtdc_s,max_rtdc_s);
    hLF1_t[i] = new TH1D(Form("hLF1_t_%d",i),Form("hLF1_t_%d",i),bin_ltdc,min_ltdc,max_ltdc);
    hLF1_tc[i] = new TH1D(Form("hLF1_tc_%d",i),Form("hLF1_tc_%d",i),bin_ltdc,min_ltdc,max_ltdc);
    hLF1_ts[i]= new TH1D(Form("hLF1_ts_%d",i),Form("hLF1_ts_%d",i),bin_ltdc_s,min_ltdc_s,max_ltdc_s);
    hLF1_b[i] = new TH1D(Form("hLF1_b_%d",i),Form("hLF1_b_%d",i),bin_ltdc,min_ltdc,max_ltdc);
    hLF1_bc[i] = new TH1D(Form("hLF1_bc_%d",i),Form("hLF1_bc_%d",i),bin_ltdc,min_ltdc,max_ltdc);
    hLF1_bs[i]= new TH1D(Form("hLF1_bs_%d",i),"",bin_ltdc_s,min_ltdc_s,max_ltdc_s);

    fRF1_t[i] = new TF1(Form("fRF1_t_%d",i),"gausn(0)",min_rtdc,max_rtdc);
    fRF1_ts[i] = new TF1(Form("fRF1_ts_%d",i),"gausn(0)",min_rtdc_s,max_rtdc_s);
    fRF1_b[i] = new TF1(Form("fRF1_b_%d",i),"gausn(0)",min_rtdc,max_rtdc);
    fRF1_bs[i] = new TF1(Form("fRF1_bs_%d",i),"gausn(0)",min_rtdc_s,max_rtdc_s);
    fLF1_t[i] = new TF1(Form("fLF1_t_%d",i),"gausn(0)",min_ltdc,max_ltdc);
    fLF1_ts[i] = new TF1(Form("fLF1_ts_%d",i),"gausn(0)",min_ltdc_s,max_ltdc_s);
    fLF1_b[i] = new TF1(Form("fLF1_b_%d",i),"gausn(0)",min_ltdc,max_ltdc);
    fLF1_bs[i] = new TF1(Form("fLF1_bs_%d",i),"gausn(0)",min_ltdc_s,max_ltdc_s);

    fRF1_t[i]->SetLineColor(2);
    fRF1_b[i]->SetLineColor(2);
    fLF1_t[i]->SetLineColor(2);
    fLF1_b[i]->SetLineColor(2);

    fRF1_ts[i]->SetLineColor(2);
    fRF1_bs[i]->SetLineColor(2);
    fLF1_ts[i]->SetLineColor(2);
    fLF1_bs[i]->SetLineColor(2);


    fRF1_t[i]->SetNpx(2000);
    fRF1_b[i]->SetNpx(2000);
    fLF1_t[i]->SetNpx(2000);
    fLF1_b[i]->SetNpx(2000);

    fRF1_ts[i]->SetNpx(2000);
    fRF1_bs[i]->SetNpx(2000);
    fLF1_ts[i]->SetNpx(2000);
    fLF1_bs[i]->SetNpx(2000);        
    
    
  }


  
  
  int nmax=500;
  double RF1[nmax],LF1[nmax];
  double Rs2_pads[100],Ls2_pads[100];
  double R_s2_trpad[100],L_s2_trpad[100];
  double R_tr_pathl[100],L_tr_pathl[100],R_s2_trpath[100],L_s2_trpath[100];
  double ct,ct_c,ct_p,ct_r,rtof,ltof,L_tr_p[100],R_tr_p[100];
  double Rzt[100],Lzt[100];
  double Ra1sum,Ra2sum;
  double Rs2_ra[100],Rs2_la[10],Ls2_ra[100],Ls2_la[100];
  double Rs2_rt[100],Rs2_lt[10],Ls2_rt[100],Ls2_lt[100];
  int runnum;
  double trig;
  
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("DR.evtypebits",1);
  T->SetBranchAddress("DR.evtypebits",&trig);  
  T->SetBranchStatus("RTDC.F1FirstHit",1);
  T->SetBranchAddress("RTDC.F1FirstHit",RF1);
  T->SetBranchStatus("LTDC.F1FirstHit",1);
  T->SetBranchAddress("LTDC.F1FirstHit",LF1); 
  T->SetBranchStatus("R.tr.vz",1);
  T->SetBranchAddress("R.tr.vz",Rzt);
  T->SetBranchStatus("L.tr.vz",1);
  T->SetBranchAddress("L.tr.vz",Lzt);
  T->SetBranchStatus("runnum",1);
  T->SetBranchAddress("runnum",&runnum);
  T->SetBranchStatus("L.s2.t_pads",1);
  T->SetBranchAddress("L.s2.t_pads",Ls2_pads);
  T->SetBranchStatus("L.tr.p"          ,1);
  T->SetBranchAddress("L.tr.p" ,L_tr_p );
  T->SetBranchStatus("L.tr.pathl"      ,1);
  T->SetBranchAddress("L.tr.pathl" ,L_tr_pathl) ; 
  T->SetBranchStatus("L.s2.trpad"      ,1);
  T->SetBranchAddress("L.s2.trpad",L_s2_trpad);
  T->SetBranchStatus("L.s2.trpath"     ,1);
  T->SetBranchAddress("L.s2.trpath",L_s2_trpath);
  
  T->SetBranchStatus("L.s2.ra_p",1);
  T->SetBranchAddress("L.s2.ra_p",Ls2_ra);
  T->SetBranchStatus("L.s2.la_p",1);
  T->SetBranchAddress("L.s2.la_p",Ls2_la);  
  T->SetBranchStatus("L.s2.rt_c",1);
  T->SetBranchAddress("L.s2.rt_c",Ls2_rt);
  T->SetBranchStatus("L.s2.lt_c",1);
  T->SetBranchAddress("L.s2.lt_c",Ls2_lt);  

  T->SetBranchStatus("R.s2.ra_p",1);
  T->SetBranchAddress("R.s2.ra_p",Rs2_ra);
  T->SetBranchStatus("R.s2.la_p",1);
  T->SetBranchAddress("R.s2.la_p",Rs2_la);  
  T->SetBranchStatus("R.s2.rt_c",1);
  T->SetBranchAddress("R.s2.rt_c",Rs2_rt);
  T->SetBranchStatus("R.s2.lt_c",1);
  T->SetBranchAddress("R.s2.lt_c",Rs2_lt);  

  T->SetBranchStatus("R.s2.t_pads",1);
  T->SetBranchAddress("R.s2.t_pads",Rs2_pads);
  T->SetBranchStatus("R.s2.trpad"      ,1);
  T->SetBranchAddress("R.s2.trpad",R_s2_trpad);  
  T->SetBranchStatus("R.tr.p"          ,1);
  T->SetBranchAddress("R.tr.p" ,R_tr_p );
  T->SetBranchStatus("R.s2.trpath"     ,1);
  T->SetBranchAddress("R.s2.trpath",R_s2_trpath);
  T->SetBranchStatus("R.tr.pathl"      ,1);
  T->SetBranchAddress("R.tr.pathl" ,R_tr_pathl) ;
  T->SetBranchStatus("R.a1.asum_p"          ,1);
  T->SetBranchAddress("R.a1.asum_p" ,&Ra1sum );
  T->SetBranchStatus("R.a2.asum_p"          ,1);
  T->SetBranchAddress("R.a2.asum_p" ,&Ra2sum );


  
  //  string ofname="F1shift_corr.root";
  //  ofname="../rootfiles/coin/test_nnL_whist_Ole3.root";
  //  ofname="../rootfiles/coin/nnL_small_Ole2_F1shift.root";
  string ofname ="../rootfiles/coin/" + oname +".root";
  TFile* ofp = new TFile(Form("%s",ofname.c_str()),"recreate");

  TTree* tnew=new TTree("T","F1TDC Calibration");//=T->CloneTree(0);
  double RS2T[16],RS2B[16],LS2T[16],LS2B[16],RS2T_ref,RS2B_ref,LS2T_ref,LS2B_ref;
  double RS2T_c[16],RS2B_c[16],LS2T_c[16],LS2B_c[16];
  int Rs2_pad,Ls2_pad;
  int ac_cut,z_cut,f1_cut;
  int T5,lhit,rhit;
  tnew->Branch("T5",&T5,"T5/I");
  tnew->Branch("RF1",RF1,"RF1[500]");
  tnew->Branch("LF1",LF1,"LF1[500]");
  tnew->Branch("RS2T_ref",&RS2T_ref,"RS2T_ref/D");
  tnew->Branch("RS2B_ref",&RS2B_ref,"RS2B_ref/D");
  tnew->Branch("LS2T_ref",&LS2T_ref,"LS2T_ref/D");
  tnew->Branch("LS2B_ref",&LS2B_ref,"LS2B_ref/D");
  tnew->Branch("RS2T",RS2T,"RS2T[16]/D");  
  tnew->Branch("RS2B",RS2B,"RS2B[16]/D");  
  tnew->Branch("LS2T",LS2T,"LS2T[16]/D");  
  tnew->Branch("LS2B",LS2B,"LS2B[16]/D");
  tnew->Branch("RS2T_c",RS2T_c,"RS2T_c[16]/D");  
  tnew->Branch("RS2B_c",RS2B_c,"RS2B_c[16]/D");  
  tnew->Branch("LS2T_c",LS2T_c,"LS2T_c[16]/D");  
  tnew->Branch("LS2B_c",LS2B_c,"LS2B_c[16]/D");    
  tnew->Branch("Rs2_pad",&Rs2_pad,"Rs2_pad/I");
  tnew->Branch("Ls2_pad",&Ls2_pad,"Ls2_pad/I");
  tnew->Branch("ct",&ct,"ct/D");
  tnew->Branch("ct_c",&ct_c,"ct_c/D");
  tnew->Branch("ct_p",&ct_p,"ct_p/D");
  tnew->Branch("ct_r",&ct_r,"ct_r/D");
  tnew->Branch("rtof",&rtof,"rtof/D");
  tnew->Branch("ltof",&ltof,"ltof/D");
  tnew->Branch("pid_cut",&ac_cut,"ac_cut/i");
  tnew->Branch("z_cut",&z_cut,"z_cut/i");
  tnew->Branch("f1_cut",&f1_cut,"f1_cut/i");
  tnew->Branch("Rzt",Rzt,"Rzt[100]/D");  
  tnew->Branch("Lzt",Lzt,"Lzt[100]/D");  
  tnew->Branch("lhit",&lhit,"lhit/I");
  tnew->Branch("rhit",&rhit,"rhit/I");
  tnew->Branch("Rs2_rt_p",Rs2_rt,"Rs2_rt[100]/D");  
  tnew->Branch("Rs2_lt_p",Rs2_lt,"Rs2_lt[100]/D");
  tnew->Branch("Rs2_ra_p",Rs2_ra,"Rs2_ra[100]/D");  
  tnew->Branch("Rs2_la_p",Rs2_la,"Rs2_la[100]/D");
  tnew->Branch("Ls2_rt_p",Ls2_rt,"Ls2_rt[100]/D");  
  tnew->Branch("Ls2_lt_p",Ls2_lt,"Ls2_lt[100]/D");
  tnew->Branch("Ls2_ra_p",Ls2_ra,"Ls2_ra[100]/D");  
  tnew->Branch("Ls2_la_p",Ls2_la,"Ls2_la[100]/D");  
  
  tnew->Branch("runnum",&runnum);
  
  for(int i=0;i<ENum;i++){
    for(int k=0;k<nmax;k++){
      RF1[k]= 0.0;
      LF1[k]= 0.0;
    }
    
    ct=-1000.0;
    ct_c=-1000.0;
    rtof=-1000;
    ltof=-1000;
    ac_cut=-1;
    z_cut=-1;
    Rs2_pad=-1;
    Ls2_pad=-1;
    f1_cut=-1;
    T5=-1;
    rhit=-1;
    lhit=-1;
    bool Rs2_adc=false;
    bool Ls2_adc=false;


    RS2T_ref = -10000.0;
    RS2B_ref = -10000.0;
    LS2T_ref = -10000.0;
    LS2B_ref = -10000.0;
    
    T->GetEntry(i);
    
    
    RS2T_ref = RF1[9];
    RS2B_ref = RF1[46];
    LS2T_ref = LF1[30];
    LS2B_ref = LF1[37];

    if(trig==32)T5=1;
    
    //    RS2T_ref = RF1[15]; //Reference signal
    //    RS2B_ref = RF1[15]; //Reference signal
    //    LS2T_ref = LF1[47]; //Reference signal
    //    LS2B_ref = LF1[47]; //Reference signal
    double RF_L=LF1[47];
    double RF_R=LF1[15];
    
    //      RS2T_ref = RF1[9];  // L1A_R
    //      RS2B_ref = RF1[9];  // L1A_R
    //      LS2T_ref = LF1[30]; // L1A remote
    //      LS2B_ref = LF1[30]; // L1A remote
  
    Rs2_pad = (int)Rs2_pads[0];
    Ls2_pad = (int)Ls2_pads[0];    


    for(int j=0;j<16;j++){

      RS2T[j]= -10000.0;
      RS2B[j]= -10000.0;
      LS2T[j]= -10000.0;
      LS2B[j]= -10000.0;
      //LS2T[j]= 12315.0;
      //LS2B[j]= 12315.0;
      RS2T_c[j]= -10000.0;
      RS2B_c[j]= -10000.0;
      //      LS2T_c[j]= 12315.0;
      //      LS2B_c[j]= 12315.0;
      LS2T_c[j]= -10000.0;
      LS2B_c[j]= -10000.0;      

      if(j==Rs2_pad){


      RS2T[j]= RF1[j+16];
      RS2B[j]= RF1[j+48];
      RS2T_c[j]= -RS2T[j] + RS2T_ref;
      RS2B_c[j]= -RS2B[j] + RS2B_ref;

      //      RS2T_c[j]= -RS2T[j] + RS2T_ref ;
      //      RS2B_c[j]= -RS2B[j] + RS2B_ref ;

      }else if(j==Ls2_pad){

      LS2T[j]= LF1[j];
      LS2B[j]= LF1[j+48];
      LS2T_c[j]= -LS2T[j] + LS2T_ref;
      LS2B_c[j]= -LS2B[j] + LS2B_ref;


      
      //      LS2T_c[j]= -LS2T[j] + LS2T_ref ;
      //      LS2B_c[j]= -LS2B[j] + LS2B_ref ;      

      }
      
    }
    


    
    const double mk= 0.493677;
    const double cc = 0.299792458;          // speed of light in vacuum (m/ns)
    //    int mode=1;
    int mode=2;
    double tdc_time;
    double coin_offset;
    if(mode==1){
      tdc_time =0.056;// [ns]
      coin_offset = -470.5 + 464.73; // H1 mode
    }else {
      tdc_time=0.058; //[ns]
      coin_offset=470.63 -470.5; // H2 mode
      //      coin_offset-=22.2+1.5;
    }


    
    if(Ra1sum<50. && Ra2sum>2000)ac_cut=1;
    else ac_cut=0;
    if(fabs(Rzt[0]+Lzt[0])/2.0<0.1 && fabs(Rzt[0]-Lzt[0])<0.03)z_cut=1;
    else z_cut=0;
    //    if()f1_cut=1;
    //    else f1_cut=0;    


    if(Rs2_ra[Rs2_pad]>100 && Rs2_la[Rs2_pad]>100)Rs2_adc=true;
    if(Ls2_ra[Ls2_pad]>100 && Ls2_la[Ls2_pad]>100)Ls2_adc=true;

    if(RS2T[Rs2_pad]>0)rhit=1;
    if(LS2T[Ls2_pad]>0)lhit=1;
    
    double Rs2_off = s2f1_off(Rs2_pad,"R",mode);
    double Ls2_off = s2f1_off(Ls2_pad,"L",mode);
    double tof_rb=( ( RS2T_c[Rs2_pad] + RS2B_c[Rs2_pad] + Rs2_off )/2.0 )*tdc_time;
    double tof_lb=( ( LS2T_c[Ls2_pad] + LS2B_c[Ls2_pad] + Ls2_off )/2.0 )*tdc_time;
    double R_betaK=R_tr_p[0]/sqrt(mk*mk + R_tr_p[0]*R_tr_p[0]);      
    double L_tgt_b = tof_lb + ( L_tr_pathl[0] + L_s2_trpath[0] )/cc;
    double R_tgt_b = tof_rb + ( R_tr_pathl[0] + R_s2_trpath[0] )/R_betaK/cc;

    
    //    if(Rs2_adc && Ls2_adc)
      ct = - L_tgt_b + R_tgt_b - coin_offset;
      //    if(Rs2_adc && Ls2_adc)
      ct_r=  R_tgt_b - ( L_tr_pathl[0] + L_s2_trpath[0] )/cc + Ls2_off/2.0*tdc_time -coin_offset ;
    
      if(ct> 3500.) ct_p=ct-3.65049e3;
      if(ct<-3500.) ct_p=ct+3.65063e3;

    
    //===== F1 Correction =======//

      double F1off=350.0-3.0;
      F1off=333;
      double F1off_LS2B;
      double F1off_LS2T;
      double F1off_RS2T;
      double F1off_RS2B;
       if(mode==1){
	F1off_RS2T=0.0;
	//	F1off_RS2B=-2.0;
	F1off_LS2T=0.0;
	//	F1off_LS2B=-3.0;
	F1off_RS2T=0.0;
	F1off_RS2B=0.0;
	
	F1off_LS2T=0.0;
	F1off_LS2B=0.0;
	

      }else if(mode==2){
	//	F1off_RS2T=0.0;
	//	F1off_RS2B=-2.0;
	
	F1off_LS2T=+280.-5.0-14.3;
	F1off_LS2B=+280.-5.0-14.3;
	F1off_RS2T=0.0;
	F1off_RS2B=0.0;	
	F1off_LS2T=0.0;
	F1off_LS2B=0.0;
 
      }
      
      
    for(int j=0;j<16;j++){

      hRF1_t[j]->Fill(RS2T[j]-RS2T_ref);
      hRF1_ts[j]->Fill(RS2T[j]-RS2T_ref);
      hRF1_b[j]->Fill(RS2B[j]-RS2B_ref);
      hRF1_bs[j]->Fill(RS2B[j]-RS2B_ref);
      
      hLF1_t[j]->Fill(LS2T[j]-LS2T_ref);
      hLF1_ts[j]->Fill(LS2T[j]-LS2T_ref);
      hLF1_b[j] ->Fill(LS2B[j]-LS2B_ref);
      hLF1_bs[j]->Fill(LS2B[j]-LS2B_ref);

      
      bool F1shift_flag=true;
      bool uniform=false;
      if(F1shift_flag){

	if(uniform){
	if(RS2T_c[j]>pow(2,15))RS2T_c[j]=RS2T_c[j] -pow(2.0,16) + (F1off + F1off_RS2T);
	else if(RS2T_c[j]<-pow(2,15))RS2T_c[j]=RS2T_c[j] + pow(2.0,16) - (F1off + F1off_RS2T);
	if(RS2B_c[j]>pow(2,15))RS2B_c[j]=RS2B_c[j] -pow(2.0,16) + (F1off + F1off_RS2B);
	else if(RS2B_c[j]<-pow(2,15))RS2B_c[j]=RS2B_c[j] + pow(2.0,16) - (F1off + F1off_RS2B);
	if(LS2T_c[j]>pow(2,15))LS2T_c[j]=LS2T_c[j] -pow(2.0,16) + ( F1off + F1off_LS2T);
	else if(LS2T_c[j]<-pow(2,15))LS2T_c[j]=LS2T_c[j] + pow(2.0,16) - ( F1off +F1off_LS2T);
	if(LS2B_c[j]>pow(2,15))LS2B_c[j]=LS2B_c[j] -pow(2.0,16) + ( F1off + F1off_LS2B);
	else if(LS2B_c[j]<-pow(2,15))LS2B_c[j]=LS2B_c[j] + pow(2.0,16) - ( F1off + F1off_LS2B);
	}else {
	  /*
	double F1shift_RS2T[16]={65619.2,65184.6,65190.3,65190,65189,65189,65190.8,65188.8,65188.9,65188.9,65188.4,65188.9,65188.6,65188,65192.9,65188.2};
	double F1shift_RS2B[16]={65151.2,65198,65184.1,65188.5,65190.9,65192.9,65189.4,65188.9,65187.7,65191.9,65190.6,65186.5,65188,65189.1,65190.2,65193};
	double F1shift_LS2T[16]={65191.4,65189.1,65188.4,65188.4,65186,65189.2,65184.7,65183.5,65191.2,65189.4,65188.5,65187.1,65107.8,65183.7,65117.6,65175.9};
	double F1shift_LS2B[16]={65192.6,65189.2,65189.1,65189.1,65189.2,65189.2,65189.1,65189.1,65189.2,65189.1,65189.2,65189.2,65189.1,65189.1,65189.1,65188.7};
	  */

	  double F1shift_RS2T[16]={65186.8,65188.7,65188.3,65190.1,65189,65188.6,65188.9,65189.8,65189.5,65189.4,65188.8,65188.8,65189.4,65189.2,65188.7,65190.2 };
	  
	  double F1shift_RS2B[16]={65199.2,65187.4,65188.8,65188.8,65189.5,65187.6,65188.2,65188.7,65189.6,65189,65188.3,65188.6,65187.8,65188.8,65189.1,65189.7 };
	  double F1shift_LS2T[16]={64911.1,64910.3,64910.4,64910.7,64910.3,64910.4,64910.7,64910.6,64911.2,64910.6,64910.2,64910.2,64910.7,64910.8,64910.5,64915.8};
	  double F1shift_LS2B[16]={64898.7,64887.9,64911.6,64895.4,64911.6,64906.8,64902.3,64909.3,64916,64911.6,64920.5,64895.6,64911.6,64911.6,64911.6,64912};
	  

	if(RS2T_c[j]!=-10000 && RS2B_c[j]!=-10000 &&
	   LS2T_c[Ls2_pad]==-10000 && LS2B_c[Ls2_pad]==-10000){
	  if(mode==1){
	  LS2T_c[j]=LS2T_c[j]+22315.0;
	  LS2B_c[j]=LS2B_c[j]+22315.0;
	  }else if(mode==2){
	    //	    LS2T_c[j]=LS2T_c[j]+21896.6; 
	    //	    LS2B_c[j]=LS2B_c[j]+21896.6;
	    LS2T_c[j]=LS2T_c[j]+21844.; 
	    LS2B_c[j]=LS2B_c[j]+21844.;
	  }
	  
	  }


	
	if(RS2T_c[j]<-pow(2,15))RS2T_c[j]=RS2T_c[j] + F1shift_RS2T[j];
	//	else if(RS2T_c[j]<-pow(2,15))RS2T_c[j]=RS2T_c[j] + F1shift_RS2T[j];
	if(RS2B_c[j]<-pow(2,15))RS2B_c[j]=RS2B_c[j] + F1shift_RS2B[j];
	//	else if(RS2B_c[j]<-pow(2,15))RS2B_c[j]=RS2B_c[j] + pow(2.0,16) - (F1off + F1off_RS2B);
	if(LS2T_c[j]<-pow(2,15))LS2T_c[j]=LS2T_c[j] + F1shift_LS2T[j];
	//	else if(LS2T_c[j]<-pow(2,15))LS2T_c[j]=LS2T_c[j] + pow(2.0,16) - ( F1off +F1off_LS2T);
	if(LS2B_c[j]<-pow(2,15))LS2B_c[j]=LS2B_c[j] + F1shift_LS2B[j];
	//	else if(LS2B_c[j]<-pow(2,15))LS2B_c[j]=LS2B_c[j] + pow(2.0,16) - ( F1off + F1off_LS2B);


	

	
	}

	
      hRF1_tc[j]->Fill(-RS2T_c[j]);      
      hRF1_bc[j]->Fill(-RS2B_c[j]);
      hLF1_ts[j]->Fill(-LS2T_c[j]);
      hLF1_bc[j]->Fill(-LS2B_c[j]);


      
	
      }//F1shift_flag
    } //end j  


       
    //===== coin time ============//

    

    double tof_r=( ( RS2T_c[Rs2_pad] + RS2B_c[Rs2_pad] + Rs2_off )/2.0 )*tdc_time;
    double tof_l=( ( LS2T_c[Ls2_pad] + LS2B_c[Ls2_pad] + Ls2_off )/2.0 )*tdc_time;
    double L_tgt = tof_l + ( L_tr_pathl[0] + L_s2_trpath[0] )/cc;
    double R_tgt = tof_r + ( R_tr_pathl[0] + R_s2_trpath[0] )/R_betaK/cc;
    
    rtof= R_tgt - coin_offset;
    ltof= L_tgt;
    ct_c = - L_tgt + R_tgt - coin_offset;
    //    if(Rs2_pad==8)ct_c=ct_c-1000.0 +35.3;
    
    tnew->Fill();

   if(i%100000==0)cout<<"Filled Events "<<i<<" / "<<ENum<<endl;
  }


  //========= Fitting =========//
  double mRF1_t[16],mRF1_b[16],mLF1_t[16],mLF1_b[16];
  double mRF1_ts[16],mRF1_bs[16],mLF1_ts[16],mLF1_bs[16];
  double mRF1_t_err[16],mRF1_b_err[16],mLF1_t_err[16],mLF1_b_err[16];
  double mRF1_ts_err[16],mRF1_bs_err[16],mLF1_ts_err[16],mLF1_bs_err[16];

  TGraphErrors* gRF1t=new TGraphErrors();
  TGraphErrors* gRF1b=new TGraphErrors();
  TGraphErrors* gLF1t=new TGraphErrors();
  TGraphErrors* gLF1b=new TGraphErrors();


  gRF1t->SetName("gRF1t");
  gRF1b->SetName("gRF1b");
  gLF1t->SetName("gLF1t");
  gLF1b->SetName("gLF1b");

  gRF1t->SetMarkerStyle(20);
  gRF1b->SetMarkerStyle(20);
  gLF1t->SetMarkerStyle(20);
  gLF1b->SetMarkerStyle(20);  
  
  TGraphErrors* gRF1t_s=new TGraphErrors();
  TGraphErrors* gRF1b_s=new TGraphErrors();
  TGraphErrors* gLF1t_s=new TGraphErrors();
  TGraphErrors* gLF1b_s=new TGraphErrors();

  gRF1t_s->SetName("gRF1t_s");
  gRF1b_s->SetName("gRF1b_s");
  gLF1t_s->SetName("gLF1t_s");
  gLF1b_s->SetName("gLF1b_s");  
  gRF1t_s->SetMarkerStyle(20);
  gRF1b_s->SetMarkerStyle(20);
  gLF1t_s->SetMarkerStyle(20);
  gLF1b_s->SetMarkerStyle(20);  
  
  TGraphErrors* gRF1t_d=new TGraphErrors();
  TGraphErrors* gRF1b_d=new TGraphErrors();
  TGraphErrors* gLF1t_d=new TGraphErrors();
  TGraphErrors* gLF1b_d=new TGraphErrors();
  
  gRF1t_d->SetName("gRF1t_d");
  gRF1b_d->SetName("gRF1b_d");
  gLF1t_d->SetName("gLF1t_d");
  gLF1b_d->SetName("gLF1b_d");  
  gRF1t_d->SetMarkerStyle(20);
  gRF1b_d->SetMarkerStyle(20);
  gLF1t_d->SetMarkerStyle(20);
  gLF1b_d->SetMarkerStyle(20);  

  gRF1t_d->SetMarkerColor(1);
  gRF1b_d->SetMarkerColor(2);
  gLF1t_d->SetMarkerColor(3);
  gLF1b_d->SetMarkerColor(4);  
  

  
  for(int i=0;i<16;i++){

    mRF1_t[i]=hRF1_t[i]->GetBinCenter(hRF1_t[i]->GetMaximumBin());
    mRF1_b[i]=hRF1_b[i]->GetBinCenter(hRF1_b[i]->GetMaximumBin());
    mLF1_t[i]=hLF1_t[i]->GetBinCenter(hLF1_t[i]->GetMaximumBin());
    mLF1_b[i]=hLF1_b[i]->GetBinCenter(hLF1_b[i]->GetMaximumBin());    
    mRF1_ts[i]=hRF1_ts[i]->GetBinCenter(hRF1_ts[i]->GetMaximumBin());
    mRF1_bs[i]=hRF1_bs[i]->GetBinCenter(hRF1_bs[i]->GetMaximumBin());
    mLF1_ts[i]=hLF1_ts[i]->GetBinCenter(hLF1_ts[i]->GetMaximumBin());
    mLF1_bs[i]=hLF1_bs[i]->GetBinCenter(hLF1_bs[i]->GetMaximumBin());    
    
    hRF1_t[i]->Fit(Form("fRF1_t_%d",i),"QR","QR",mRF1_t[i]-50.,mRF1_t[i]+50.);
    hRF1_b[i]->Fit(Form("fRF1_b_%d",i),"QR","QR",mRF1_b[i]-50.,mRF1_b[i]+50.);
    hRF1_ts[i]->Fit(Form("fRF1_ts_%d",i),"QR","QR",mRF1_ts[i]-50.,mRF1_ts[i]+50.);
    hRF1_bs[i]->Fit(Form("fRF1_bs_%d",i),"QR","QR",mRF1_bs[i]-50.,mRF1_bs[i]+50.);    

    if(i<13)hLF1_t[i]->Fit(Form("fLF1_t_%d",i),"QR","QR",mLF1_t[i]-50.,mLF1_t[i]+50.);
    else hLF1_t[i]->Fit(Form("fLF1_t_%d",i),"QR","QR",mLF1_t[i]-200.,mLF1_t[i]+200.);

    hLF1_b[i]->Fit(Form("fLF1_b_%d",i),"QR","QR",mLF1_b[i]-50.,mLF1_b[i]+50.);
    if(i<13)hLF1_ts[i]->Fit(Form("fLF1_ts_%d",i),"QR","QR",mLF1_ts[i]-50.,mLF1_ts[i]+50.);
    else hLF1_ts[i]->Fit(Form("fLF1_ts_%d",i),"QR","QR",mLF1_ts[i]-200.,mLF1_ts[i]+200.);
    hLF1_bs[i]->Fit(Form("fLF1_bs_%d",i),"QR","QR",mLF1_bs[i]-50.,mLF1_bs[i]+50.);        

    mRF1_t[i]=fRF1_t[i]->GetParameter(1);
    mRF1_b[i]=fRF1_b[i]->GetParameter(1);
    mLF1_t[i]=fLF1_t[i]->GetParameter(1);
    mLF1_b[i]=fLF1_b[i]->GetParameter(1);
    mRF1_ts[i]=fRF1_ts[i]->GetParameter(1);
    mRF1_bs[i]=fRF1_bs[i]->GetParameter(1);
    mLF1_ts[i]=fLF1_ts[i]->GetParameter(1);
    mLF1_bs[i]=fLF1_bs[i]->GetParameter(1);

    mRF1_t_err[i]=fRF1_t[i]->GetParError(1);
    mRF1_b_err[i]=fRF1_b[i]->GetParError(1);
    mLF1_t_err[i]=fLF1_t[i]->GetParError(1);
    mLF1_b_err[i]=fLF1_b[i]->GetParError(1);
    mRF1_ts_err[i]=fRF1_ts[i]->GetParError(1);
    mRF1_bs_err[i]=fRF1_bs[i]->GetParError(1);
    mLF1_ts_err[i]=fLF1_ts[i]->GetParError(1);
    mLF1_bs_err[i]=fLF1_bs[i]->GetParError(1);
    

    gRF1t->SetPoint(i,i,mRF1_t[i]);
    gRF1t->SetPointError(i,0,mRF1_t_err[i]);
    gRF1b->SetPoint(i,i,mRF1_b[i]);
    gRF1b->SetPointError(i,0,mRF1_b_err[i]);
    gLF1t->SetPoint(i,i,mLF1_t[i]);
    gLF1t->SetPointError(i,0,mLF1_t_err[i]);
    gLF1b->SetPoint(i,i,mLF1_b[i]);
    gLF1b->SetPointError(i,0,mLF1_b_err[i]);    
    

    gRF1t_s->SetPoint(i,i,mRF1_ts[i]);
    gRF1t_s->SetPointError(i,0,mRF1_ts_err[i]);
    gRF1b_s->SetPoint(i,i,mRF1_bs[i]);
    gRF1b_s->SetPointError(i,0,mRF1_bs_err[i]);
    gLF1t_s->SetPoint(i,i,mLF1_ts[i]);
    gLF1t_s->SetPointError(i,0,mLF1_ts_err[i]);
    gLF1b_s->SetPoint(i,i,mLF1_bs[i]);
    gLF1b_s->SetPointError(i,0,mLF1_bs_err[i]);        

    
    
  }

  for(int i=0;i<16;i++){

    gRF1t_d->SetPoint(i,i,mRF1_ts[i]-mRF1_t[i]-(mRF1_ts[8]-mRF1_t[8]));
    gRF1t_d->SetPointError(i,0,sqrt(mRF1_ts_err[i]*mRF1_ts_err[i] + mRF1_t_err[i]*mRF1_t_err[i]));
    gRF1b_d->SetPoint(i,i,mRF1_bs[i]-mRF1_b[i]-(mRF1_bs[8]-mRF1_b[8]));
    gRF1b_d->SetPointError(i,0,sqrt(mRF1_bs_err[i]*mRF1_bs_err[i] + mRF1_b_err[i]*mRF1_b_err[i]));
    gLF1t_d->SetPoint(i,i,mLF1_ts[i]-mLF1_t[i]-(mLF1_ts[8]-mLF1_t[8]));
    gLF1t_d->SetPointError(i,0,sqrt(mLF1_ts_err[i]*mLF1_ts_err[i] + mLF1_t_err[i]*mLF1_t_err[i]));
    gLF1b_d->SetPoint(i,i,mLF1_bs[i]-mLF1_b[i]-(mLF1_bs[8]-mLF1_b[8]));
    gLF1b_d->SetPointError(i,0,sqrt(mLF1_bs_err[i]*mLF1_bs_err[i] + mLF1_b_err[i]*mLF1_b_err[i]));        


  }


  
  //====== Write ========//
  
  tnew->Write();

  gRF1t->Write();
  gRF1b->Write();
  gLF1t->Write();
  gLF1b->Write();
  gRF1t_s->Write();
  gRF1b_s->Write();
  gLF1t_s->Write();
  gLF1b_s->Write();
  gRF1t_d->Write();
  gRF1b_d->Write();
  gLF1t_d->Write();
  gLF1b_d->Write();

  
  for(int j=0;j<16;j++){
    hRF1_t[j]->Write();
    hRF1_ts[j]->Write();
    hRF1_tc[j]->Write();
    hRF1_b[j]->Write();
    hRF1_bs[j]->Write();
    hRF1_bc[j]->Write();
    hLF1_t[j]->Write();
    hLF1_ts[j]->Write();
    hLF1_tc[j]->Write();
    hLF1_b[j]->Write();
    hLF1_bs[j]->Write();
    hLF1_bc[j]->Write();

    fRF1_t[j]->Write();
    fRF1_ts[j]->Write();
    fRF1_b[j]->Write();
    fRF1_bs[j]->Write();

    
  }


  //======= Print ==============//
  //  string pdf="test.pdf";
  string pdf="../pdf/coin/" +oname +".pdf";

  
  for(int i=0;i<16;i++){
    c[i]->cd(1);
    hRF1_t[i]->Draw();
    c[i]->cd(2);
    hRF1_ts[i]->Draw();
    c[i]->cd(3);
    hRF1_b[i]->Draw();
    c[i]->cd(4);
    hRF1_bs[i]->Draw();
    c[i]->cd(5);
    hLF1_t[i]->Draw();
    c[i]->cd(6);
    hLF1_ts[i]->Draw();
    c[i]->cd(7);
    hLF1_b[i]->Draw();
    c[i]->cd(8);
    hLF1_bs[i]->Draw();        
  }
  c[16]->cd();
  gRF1t_d->GetYaxis()->SetRangeUser(-10,10);
  gRF1t_d->Draw("AP");
  gRF1b_d->Draw("P");
  gLF1t_d->Draw("P");
  gLF1b_d->Draw("P");

  c[0]->Print(Form("%s[",pdf.c_str()));

  for(int i=0;i<17;i++)
   c[i]->Print(Form("%s",pdf.c_str()));
  
  c[16]->Print(Form("%s]",pdf.c_str()));


  
  //======= Parameter Out ==========//

  //  string param="param/nnL_small_Ole2_F1shift.dat";
  //    string param="param/nnL_small_Ole2_F1shift.dat";
  string param="param/" + oname + ".dat";
  ofstream ofparam(param.c_str());
  
  string id[4];
  id[0]="RF1T_shift = ";
  id[1]="RF1B_shift = ";
  id[2]="LF1T_shift = ";
  id[3]="LF1B_shift = ";
  double shift[16];
  for(int j=0;j<4;j++){
    ofparam << id[j];
    for(int i=0;i<16;i++){
      if(j==0)shift[i]= - ( mRF1_ts[i]-mRF1_t[i] );
      else if(j==1)shift[i]= - ( mRF1_bs[i]-mRF1_b[i] );
      else if(j==2)shift[i]=mLF1_ts[i]-mLF1_t[i];
      else if(j==3)shift[i]=mLF1_bs[i]-mLF1_b[i];
      ofparam <<shift[i] <<" ";
      }
    ofparam << endl;
  }
  
  ofparam.close();

  cout<<"===== Comment Out ======="<<endl;
  cout<<"input  file "<<ifname<<endl;
  cout<<"output file "<<ofname<<endl;
  
  ofp->Close();
  
 
 

}





// ####################################################
double s2f1_off(int i,string ARM, int KINE){
// ####################################################

  double RS2_offset[16],LS2_offset[16];
  if(KINE==2){

    double RS2_R_off[16]={-8361.42, -8395.25, -8414.89, -8419.06, -8362.64, -8381.55, -8370.53, -8392.66, -8389.77, -8393.96, -8388.11, -8381.73, -8333.95, -8348.93, -8363.93, -8360.30};
    double RS2_L_off[16]={-8473.92, -8470.25, -8489.89, -8494.06, -8512.64, -8494.05, -8520.53, -8505.16, -8502.27, -8468.96, -8500.61, -8494.23, -8521.45, -8498.93, -8476.43, -8472.80};
    double LS2_R_off[16]={-12441.14, -12490.70, -12579.43, -12601.39, -12471.56, -12471.38, -12658.08, -12656.28, -12690.65, -12489.77, -12701.56, -12675.30, -12696.36, -12491.35, -12709.36, -12539.99 };
    double LS2_L_off[16]={-12141.14, -12190.70, -12091.93, -12076.39, -12209.06, -12208.88, -12058.08, -12056.28, -12015.65, -12227.27, -12026.56, -12000.30, -11983.86, -12191.35, -11996.86, -12239.99 };

    double  RS2_off_H2[16];
    double  LS2_off_H2[16];
    for(int l=0;l<16;l++){
      RS2_off_H2[l]=RS2_R_off[l] + RS2_L_off[l];
      LS2_off_H2[l]=LS2_R_off[l] + LS2_L_off[l];
    }
      
  LS2_offset[i]=LS2_off_H2[i];
  RS2_offset[i]=RS2_off_H2[i];
  }


  if(KINE==1){

double  RS2_off_H1[16]={-16828.7,-16863,-16894,-16893.3,-16870.9,-16867.2,-16900.3,-16876.8,-16895.6,-16861.6,-16895,-16890.7,-16854.6,-16852.9,-16850.5,-16861.9};
double  LS2_off_H1[16]={-25335,-25385.6,-25367,-25392.1,-25391.7,-25386.4,-25422.1,-25428.9,-25414.9,-25424.7,-25436.9, -25381.2,-25390,-25413.4,-25428.7,-26640.8};
  LS2_offset[i]=LS2_off_H1[i];
  RS2_offset[i]=RS2_off_H1[i];
  }

 double s2f1_offset; 
 if(ARM=="R")s2f1_offset=RS2_offset[i];
 else  if(ARM=="L")s2f1_offset=LS2_offset[i];
 else {cout<<"false read out !!"<<endl;}

  return s2f1_offset;

}

