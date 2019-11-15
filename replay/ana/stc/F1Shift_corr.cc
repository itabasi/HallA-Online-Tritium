#include <iostream>
#include <fstream>
//const double MK = 0.493677;            // charged Kaon mass (GeV/c2)
extern double s2f1_off(int i,string ARM, int KINE);

void F1Shift_corr(){


  TChain* oldtree=new TChain("T");

  int ENum;
  //  string ifname="/data2/small/tritium_111160-111220.root";
  //  string ifname="../run_list/nnlambda/Lambda_small_H1.list";
  //    string ifname="../run_list/nnlambda/Lambda_small_H2.list";
  //  string ifname="../run_list/nnlambda/Lambda_small_opt_H2.list";
  //  string ifname="../run_list/nnlambda/nnL_small_opt4.list";
  //   string ifname="../run_list/nnlambda/Lambda_small_optH2_test.list";
  string ifname="../run_list/nnlambda/Lambda_small_optH1test.list";
  //  string ifname="../run_list/nnlambda/test.list";
  //string ifname="../run_list/nnlambda/nnL_small_opt4.list";
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


  

  ENum =T->GetEntries();
  cout<<"Get Entries: "<<ENum<<endl;

  
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
  T->SetBranchStatus("*",0);    
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


  
  string ofname="F1shift_corr.root";
  //  ofname="F1shift_corr_test1104_1.root";
  //  ofname="F1shift_corr_optH2_1.root";
  //  ofname="F1shift_corr_optnnL_4.root";
  //  ofname="F1shift_corr_optH2_test2.root";
  //  ofname="F1shift_corr_optnnL4.root";
  ofname="../rootfiles/tcoin/F1_Lambda_small_optH1.root";  
  TFile* ofp = new TFile(Form("%s",ofname.c_str()),"recreate");

  TTree* tnew=new TTree("T","F1TDC Calibration");//=T->CloneTree(0);
  double RS2T[16],RS2B[16],LS2T[16],LS2B[16],RS2T_ref,RS2B_ref,LS2T_ref,LS2B_ref;
  double RS2T_c[16],RS2B_c[16],LS2T_c[16],LS2B_c[16];
  int Rs2_pad,Ls2_pad;
  int ac_cut,z_cut;
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
  tnew->Branch("Rzt",Rzt,"Rzt[100]/D");  
  tnew->Branch("Lzt",Lzt,"Lzt[100]/D");  

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
    bool Rs2_adc=false;
    bool Ls2_adc=false;

    
    T->GetEntry(i);
    
    RS2T_ref = 0.0;
    RS2B_ref = 0.0;
    LS2T_ref = 0.0;
    LS2B_ref = 0.0;
    
    RS2T_ref = RF1[9];
    RS2B_ref = RF1[46];
    LS2T_ref = LF1[30];
    LS2B_ref = LF1[37];

    
    //    RS2T_ref = RF1[15]; //Reference signal
    //    RS2B_ref = RF1[15]; //Reference signal
    //    LS2T_ref = LF1[47]; //Reference signal
    //    LS2B_ref = LF1[47]; //Reference signal
    double RF_L=LF1[47];
    double RF_R=LF1[15];
    
    //    RS2T_ref = RF1[9];  // L1A_R
    //    RS2B_ref = RF1[9];  // L1A_R
    //    LS2T_ref = LF1[30]; // L1A remote
    //    LS2B_ref = LF1[30]; // L1A remote
  



    for(int j=0;j<16;j++){

      RS2T[j]= 0.0;
      RS2B[j]= 0.0;
      LS2T[j]= 0.0;
      LS2B[j]= 0.0;
      
      RS2T[j]= RF1[j+16];
      RS2B[j]= RF1[j+48];
      LS2T[j]= LF1[j];
      LS2B[j]= LF1[j+48];


      RS2T_c[j]= -RS2T[j] + RS2T_ref;
      RS2B_c[j]= -RS2B[j] + RS2B_ref;
      LS2T_c[j]= -LS2T[j] + LS2T_ref;
      LS2B_c[j]= -LS2B[j] + LS2B_ref;

      RS2T_c[j]= -RS2T[j] + RS2T_ref ;
      RS2B_c[j]= -RS2B[j] + RS2B_ref ;
      LS2T_c[j]= -LS2T[j] + LS2T_ref ;
      LS2B_c[j]= -LS2B[j] + LS2B_ref ;      

    }


    const double mk= 0.493677;
    const double cc = 0.299792458;          // speed of light in vacuum (m/ns)
    int mode=1;
    double tdc_time;
    double coin_offset;
    if(mode==1){
      tdc_time =0.056;// [ns]
      coin_offset = -470.5 + 464.73; // H1 mode
    }else {
      tdc_time=0.058; //[ns]
      coin_offset=470.63 -470.5; // H2 mode 
    }


    Rs2_pad = (int)Rs2_pads[0];
    Ls2_pad = (int)Ls2_pads[0];    

    
    if(Ra1sum<50. && Ra2sum>2000)ac_cut=1;
    else ac_cut=0;
    if(fabs(Rzt[0]+Lzt[0])/2.0<0.1 && fabs(Rzt[0]-Lzt[0])<0.03)z_cut=1;
    else z_cut=0;


    if(Rs2_ra[Rs2_pad]>100 && Rs2_la[Rs2_pad]>100)Rs2_adc=true;
    if(Ls2_ra[Ls2_pad]>100 && Ls2_la[Ls2_pad]>100)Ls2_adc=true;

    
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
    
    //    if(ct> 3500.) ct_p=ct-3.65049e3;
    //    if(ct<-3500.) ct_p=ct+3.65063e3;

    
    //===== F1 Correction =======//

      double F1off=350;
      F1off=0.0;
      double F1off_LS2B;
      double F1off_LS2T;
      double F1off_RS2T;
      double F1off_RS2B;
      if(mode==1){
	//	F1off_RS2T=0.0;
	//	F1off_RS2B=-2.0;
	
	//	F1off_LS2T=0.0;
	//	F1off_LS2B=-3.0;
	F1off_RS2T=0.0;
	F1off_RS2B=0.0;
	
	F1off_LS2T=0.0;
	F1off_LS2B=0.0;
	

      }else if(mode==2){
	//	F1off_RS2T=0.0;
	//	F1off_RS2B=-2.0;
	
	//	F1off_LS2T=+280.-5.0;
	//	F1off_LS2B=+280.-5.0;
	F1off_RS2T=0.0;
	F1off_RS2B=0.0;	
	F1off_LS2T=0.0;
	F1off_LS2B=0.0;
	


      }
      
      
    for(int j=0;j<16;j++){
      bool F1shift_flag=true;

      if(F1shift_flag){


	if(RS2T_c[j]>pow(2,15))RS2T_c[j]=RS2T_c[j] -pow(2.0,16) + (F1off + F1off_RS2T);
	else if(RS2T_c[j]<-pow(2,15))RS2T_c[j]=RS2T_c[j] + pow(2.0,16) - (F1off + F1off_RS2T);
	if(RS2B_c[j]>pow(2,15))RS2B_c[j]=RS2B_c[j] -pow(2.0,16) + (F1off + F1off_RS2B);
      else if(RS2B_c[j]<-pow(2,15))RS2B_c[j]=RS2B_c[j] + pow(2.0,16) - (F1off + F1off_RS2B);

      
      if(LS2T_c[j]>pow(2,15))LS2T_c[j]=LS2T_c[j] -pow(2.0,16) + ( F1off + F1off_LS2T);
      else if(LS2T_c[j]<-pow(2,15))LS2T_c[j]=LS2T_c[j] + pow(2.0,16) - ( F1off +F1off_LS2T);
      if(LS2B_c[j]>pow(2,15))LS2B_c[j]=LS2B_c[j] -pow(2.0,16) + ( F1off + F1off_LS2B);
      else if(LS2B_c[j]<-pow(2,15))LS2B_c[j]=LS2B_c[j] + pow(2.0,16) - ( F1off + F1off_LS2B);

	
      }//F1shift_flag
    }  


       
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

