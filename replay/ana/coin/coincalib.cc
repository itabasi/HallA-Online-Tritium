using namespace std;
#include "coincalib.h"
#include "Param.h"
#include <TMinuit.h>
extern void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
extern double calcf_ct(double* P, double xf, double xpf, double yf, double ypf,double zt);



///////////////////////////////////////////////////////////////////////////

void coincalib::SetRoot(string ifname){
  add_tree(ifname);
  pack_tree();
  readtreeHRSR();
  readtreeHRSL();

}
//////////////////////////////////////////////////////////////////////////

void coincalib::SetRunList(string ifname){

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    add_tree(runname);
  }

  pack_tree();
  readtreeHRSR();
  readtreeHRSL();
}
////////////////////////////////////////////////////////////////////////////

void coincalib::ReadParam(string name){

  param = new ParamMan(name.c_str());
  cout<<"param name : "<<name<<endl;
  if(param -> SetVal())cout<<"F1TDC parameter setted"<<endl; 
  tdc_time=param->F1Res();
  coin_offset=param->GetF1CoinOffset();
  
}

/////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////

void coincalib::MTParam_R(string mtparam){



  //====== RHRS Momentum parameters ========//
    char name_Mpt[500];
    sprintf(name_Mpt, mtparam.c_str()); // optimized
    ifstream Mpt(name_Mpt);
   if (Mpt.fail()){ cerr << "failed open files" <<name_Mpt<<endl; exit(1);}
   for(int i=0;i<nParamTc;i++){
    double par=0.;
    int p=0;
    Mpt >> par >> p >> p >> p >> p >> p;
    Prp[i]=par;
    Pct_R[i]=par;
   }

 
   Mpt.close();

}

////////////////////////////////////////////////////////////////////////////

void coincalib::MTParam_L(string mtparam){

  

   

  //====== LHRS TOF parameters ========//

  char name_Mpt_L[500];
    sprintf(name_Mpt_L, mtparam.c_str()); // optimized
    ifstream Mpt_L(name_Mpt_L);
    if (Mpt_L.fail()){ cerr << "failed open files" <<name_Mpt_L<<endl; exit(1);}
   for(int i=0;i<nParamTc;i++){
    double par=0.;
    int p=0;
    Mpt_L >> par >> p >> p >> p >> p >> p;
    Pct_L[i]=par;
    Pct[i+nParamTc]=par;  // Both momentum paramters
   }

   Mpt_L.close();

}

////////////////////////////////////////////////////////////////////////////

void coincalib::NewRoot(string ofname){

  ofr = new TFile(Form("%s",ofname.c_str()),"recreate");
  tnew=new TTree("T","Coincalib matrix tuning");
  //  tnew=tree->CloneTree(0);
  tnew->Branch("ct",&coint);
  tnew->Branch("ct_c",&coint_c);
  tnew->Branch("RS2T",RS2_F1time);
  tnew->Branch("LS2T",LS2_F1time);
  tnew->Branch("RS2T_ref",&RF1Ref[0]);
  tnew->Branch("LS2T_ref",&LF1Ref[0]);
  tnew->Branch("rtof",&Rtof);
  tnew->Branch("ltof",&Ltof);
  tnew->Branch("rtof_c",&Rtof_c);
  tnew->Branch("ltof_c",&Ltof_c);  
  hcoin_cut=new TH1D("hcoin_cut","Coin with cut ",100,-20.,20.);
}


////////////////////////////////////////////////////////////////////////////

double coincalib::CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit){

  double cointime=0.0;

  convertF1TDCR(param);
  convertF1TDCL(param);

  double Rpathl=R_tr_pathl[rhit]+R_s2_trpath[rhit];
  double Lpathl=L_tr_pathl[lhit]+L_s2_trpath[lhit];
  double Beta_R=R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+MK*MK);
  double Beta_L=L_tr_p[lhit]/sqrt(L_tr_p[lhit]*L_tr_p[lhit]+Me*Me);
  double tof_r=RS2_F1time[RS2_seg] - Rpathl/(Beta_R*LightVelocity);
  double tof_l=LS2_F1time[LS2_seg] - Lpathl/(Beta_L*LightVelocity);
  Rtof=tof_r;
  Ltof=tof_l;

  
  if(RS2_F1time[RS2_seg]!=-9999. &&LS2_F1time[LS2_seg]!=-9999.){
    cointime= - tof_r + tof_l - coin_offset;

  }
  else{
    cointime=-1000;
  }

  //  cout<<"RF1 "<<RS2_F1time[RS2_seg]<<" LF1 "<<LS2_F1time[LS2_seg]<<
  //    " Rseg "<<RS2_seg<<" Lseg "<<LS2_seg<<endl;
  //  cout<<"rtof "<<tof_r<<" ltof "<<tof_l<<" ct "<<cointime<<endl; 
  
  return cointime;
  
}
///////////////////////////////////////////////////////////////////////////

double coincalib::tune(double* pa, int MODE) 
{
  
  double chi = 0.0;
  double arglist[10]; 
  int ierflg = 0;
  int allparam = nParamTc;
  int arm=1;

  if(MODE==0){
    allparam=nParamTc*2;
    arm=2;
  }

  //  for(int i=0;i<allparam;i++)cout<<"i "<<i<<" OptParam "<<pa[i]<<endl;

  cout<<"mode "<<MODE<<" allParam "<<allparam<<endl;


  TMinuit* minuit= new TMinuit(allparam);
  minuit->SetFCN(fcn); // fcn Chi-square function

  
  double start[allparam];
  double step[allparam];

  const int nMatT =nnc;  
  const int nXf   =nnc;
  const int nXpf  =nnc;
  //  const int nYf   =nnc;
  const int nYf   =0;
  //  const int nYpf  =nnc;
  const int nYpf  =0;
  const int nZt   =0;



  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;

  for(int f=0;f<arm;f++){
    for (int n=0;n<nMatT+1;n++){
      for(e=0;e<n+1;e++){
	for (d=0;d<n+1;d++){
	  for (c=0;c<n+1;c++){ 
	    for (b=0;b<n+1;b++){
	      for (a=0;a<n+1;a++){ 
		if (a+b+c+d+e==n){
		  if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		    start[npar] = pa[npar];
		    step[npar] = 1.0e-3;
		    
		  }
		  else{
		    start[npar] = 0.0;
		    step[npar] = 0.0;
		  }
		  npar++;

		}
	      }
	    }
	  }
	}    
      }
    }
  }


  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  minuit -> SetPrintLevel(-1);
  
  double LLim[allparam];// Lower limit for all of the parameter
  double ULim[allparam];// Upper limit for all of the parameter

  char pname[500];
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
 
    //    LLim[i] = pa[i] - pa[i]*0.8;
    //    ULim[i] = pa[i] + pa[i]*0.8;
    

    LLim[i] = pa[i] - 10.0; // temp
    ULim[i] = pa[i] + 10.0; // temp

    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
	  
  }

  
  // ~~~~ Strategy ~~~~
  //    arglist[0] = 2.0; // original
   arglist[0] = 1.0; // test
  //arglist[0] = 0.0;   // test
  minuit->mnexcm("SET STR",arglist,1,ierflg);
 
  // ~~~~ Migrad + Simplex  ~~~~ 
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg); // Chi-square minimization
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double er;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);

  if(amin>0) chi=amin;

  
  for(int i=0 ; i<allparam ; i++){
    if(MODE==0){
      minuit -> GetParameter(i,Pct[i],er);
      if(i<allparam/2){minuit -> GetParameter(i,Pct_R[i],er);}// RHRS momentum parameter
      else{minuit -> GetParameter(i,Pct_L[i-nParamTc],er);}// RHRS momentum parameters
      
      
    }else if(MODE==-1){
      minuit -> GetParameter(i,Pct_R[i],er);// RHRS momentum parameters
      //      cout<<"i "<<i<<" Pct_R "<<Pct_R[i]<<endl;
    }else if(MODE==1){minuit -> GetParameter(i,Pct_L[i],er);// LHRS momentum parameters
    }
  }
  
  return chi;
}

//////////////////////////////////////////////////////////////////////////////

void coincalib::EventSelect(){

  cout<<"==========================="<<endl;
  cout<<"===== Event Selection ======"<<endl;
  cout<<"==========================="<<endl;

  cout<<"Events : "<<ENum<<endl;

  for(int k=0;k<ENum;k++){

    tree->GetEntry(k);

    // ===== Initialization ====== //    

    z_flag=false;
    pid_flag=false;
    ct_flag=false;
    ct_flag=false;
    
 /////////////////////
 //// Coincidence ////
 /////////////////////


    
    if(R_evtype==5){
      int NLtr = (int)L_tr_n;  if(NLtr>MAX) NLtr = MAX;
      int NRtr = (int)R_tr_n;  if(NRtr>MAX) NRtr = MAX;
      
      for(int lt=0;lt<NLtr;lt++){	
        for(int rt=0;rt<NRtr;rt++){

	  if(fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025)z_flag=true;
	  if(R_a1_asum_p<50. && R_a2_asum_p>2000. && L_cer_asum_c>2000.)pid_flag=true;

	  int R_s2pad=(int)R_s2_t_pads[rt];
	  int L_s2pad=(int)L_s2_t_pads[lt]; 
	  
	  if( fabs( CoinCalc(R_s2pad,L_s2pad,rt,lt) )<1.0 )ct_flag=true;
	  if(ntune_event<nmax && z_flag && pid_flag && ct_flag){
	  ct[ntune_event]    = CoinCalc(R_s2pad,L_s2pad,rt,lt);
	  rx_fp[ntune_event] = R_tr_x[rt];
	  ry_fp[ntune_event] = R_tr_y[rt];
	  rth_fp[ntune_event] = R_tr_th[rt];
	  rph_fp[ntune_event] = R_tr_ph[rt];
	  lx_fp[ntune_event] = L_tr_x[lt];
	  ly_fp[ntune_event] = L_tr_y[lt];
	  lth_fp[ntune_event] = L_tr_th[lt];
	  lph_fp[ntune_event] = L_tr_ph[lt];
	  RTOF[ntune_event]=Rtof;
	  RTOF[ntune_event]=Ltof;
	  hcoin_cut->Fill(ct[ntune_event]);
	  ntune_event++;	  
          if(ntune_event%100==0&& ntune_event>10)cout<<"Event Selection : "<<ntune_event<<" / "<<nmax<<endl;


	  }	    
	}//end NRtr
      }//end NLtr
      if(k%100000==0)cout<<"Filled : "<<k<<" / "<<ENum<<endl;
      if(ntune_event==nmax)break;
  }//end loop
  }


  cout<<"Events Selection : "<<ntune_event<<endl;
  
}

////////////////////////////////////////////////////////////////////////////


void coincalib::CoinTuning(string ofname,int MODE){


  if (nite>0){
    cout << "======================================================" <<endl;
    cout << "=================  Tuning started ====================" << endl;
    cout << "======================================================" <<endl;}
    char tempc[500],tempc2[500];
    const  char* new_tempc=ofname.c_str();
    cout<<"new marix file: "<<new_tempc<<endl;
    
    ofstream * ofs1;
    ofstream * ofs2;       
    for(int i=0 ; i<nite+1 ; i++){

      if(i>0){    cout<<"tuning i: "<<i<<" /"<<nite<<endl;
    // --------------------------- //
    // ---- Parameter tuning ----- //
    // --------------------------- //
    
    if(MODE==-1){
    cout<<"------- Rp tuning -----"<<endl;
    chi_sq1[i]=0.0;
    chi_sq1[i] = tune(Pct_R,-1);   // Rp
    }else if(MODE==1){
    cout<<"------- Lp tuning -----"<<endl;    
    chi_sq2[i]=0.0;
    chi_sq2[i] = tune(Pct_L,1);  // Lp
    }else if(MODE==0){
    cout<<"---------pk & pe tuning ------"<<endl;
    chi_sq[i]=0.0;
    chi_sq[i] = tune(Pct,0);
    //    chi_sq[i] = tuning(Pct,i,0); // function
    }

    cout << " Tuning# = " << i << ": chisq = ";
      }
    
    if(MODE==-1 || MODE==0){
    sprintf(tempc,  "%s_Rtof%d.dat",new_tempc,i); 
    cout<<"new matrix Rp: "<<tempc<<endl;
    ofs1 = new ofstream(tempc);}

    if(MODE==1 ||MODE==0){
    sprintf(tempc2, "%s_Ltof%d.dat",new_tempc,i);
    cout<<"new matrix Lp: "<<tempc2<<endl;    
    ofs2 = new ofstream(tempc2);}


    int nppp = 0;
    
    for(int i=0 ; i<nnp+1 ; i++){
      for(int e=0 ; e<nnp+1 ; e++){
	for(int d=0 ; d<nnp+1 ; d++){
	  for(int c=0 ; c<nnp+1 ; c++){
	    for(int b=0 ; b<nnp+1 ; b++){
	      for(int a=0 ; a<nnp+1 ; a++){  
		if(a+b+c+d+e==i){
		  if(MODE==-1){		  
		  *ofs1 << Pct_R[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  }else if(MODE==1){
		  *ofs2 << Pct_L[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		  }else if(MODE==0){

		    *ofs1 << Pct[nppp] 
			  << " " << a 
			  << " " << b
			  << " " << c
			  << " " << d
			  << " " << e << endl;
		    *ofs2 << Pct[nParamTc+nppp] 
			  << " " << a 
			  << " " << b
			  << " " << c
			  << " " << d
			  << " " << e << endl;		    
		    

		    /*
		    *ofs1 << Pct_R[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		    *ofs2 << Pct_L[nppp] 
			<< " " << a 
			<< " " << b
			<< " " << c
			<< " " << d
			<< " " << e << endl;
		    */
		    
		  }

		  nppp++;

		}
	      }
	    }
	  }
	}
      }
    }



      

      
    if(MODE==-1 || MODE==0){
      ofs1->close();
      ofs1->clear();
    }else if(MODE==0 || MODE==1){
      ofs2->close();
      ofs2->clear();}
    

    
}

    
  cout<<"========== Tuning is done ============="<<endl;
  
  

}


///////////////////////////////////////////////////////////////////////////

void coincalib::Fill(){

  cout<<"=================================="<<endl;
  cout<<"======== Fill Events ============="<<endl;
  cout<<"=================================="<<endl;
  double tof_r,tof_l;
  for(int k;k<ENum;k++){
    tof_r= -2222;
    tof_l= +2222;
    tree->GetEntry(k);

    int NLtr = (int)L_tr_n;  //if(NLtr>MAX) NLtr = MAX;
    int NRtr = (int)R_tr_n;  //if(NRtr>MAX) NRtr = MAX;

    for(int lt=0;lt<NLtr;lt++){
      for(int rt=0;rt<NRtr;rt++){

	  int R_s2pad=(int)R_s2_t_pads[rt];
	  int L_s2pad=(int)L_s2_t_pads[lt]; 

      R_tr_x[rt]  = (R_tr_x[rt]-XFPm)/XFPr;
      R_tr_th[rt] = (R_tr_th[rt]-XpFPm)/XpFPr;
      R_tr_y[rt]  = (R_tr_y[rt]-YFPm)/YFPr;
      R_tr_ph[rt] = (R_tr_ph[rt]-YpFPm)/YpFPr;

      L_tr_x[rt]  = (L_tr_x[rt]-XFPm)/XFPr;
      L_tr_th[rt] = (L_tr_th[rt]-XpFPm)/XpFPr;
      L_tr_y[rt]  = (L_tr_y[rt]-YFPm)/YFPr;
      L_tr_ph[rt] = (L_tr_ph[rt]-YpFPm)/YpFPr;	        

      coint   = CoinCalc(R_s2pad,L_s2pad,rt,lt);      
      tof_r=Rtof;
      tof_l=Ltof;
      
      if(MODE<=0)tof_r=+tof_r*calcf_ct(Pct_R,R_tr_x[rt],R_tr_th[rt],R_tr_y[rt],R_tr_ph[rt],1.0);
      if(MODE>=0)tof_l=+tof_l*calcf_ct(Pct_L,L_tr_x[lt],L_tr_th[lt],L_tr_y[lt],L_tr_ph[lt],1.0);

      Rtof_c=tof_r;
      Ltof_c=tof_l;      

      R_tr_x[rt] = R_tr_x[rt] * XFPr + XFPm;
      R_tr_th[rt] = R_tr_th[rt] * XpFPr + XpFPm;
      R_tr_y[rt] = R_tr_y[rt] * YFPr + YFPm;
      R_tr_ph[rt] = R_tr_ph[rt] * YpFPr   + YpFPm;
      
      L_tr_x[rt] = L_tr_x[rt] * XFPr + XFPm;
      L_tr_th[rt] = L_tr_th[rt] * XpFPr + XpFPm;
      L_tr_y[rt] = L_tr_y[rt] * YFPr + YFPm;
      L_tr_ph[rt] = L_tr_ph[rt] * YpFPr   + YpFPm;

	coint_c = tof_r - tof_l;

	tnew->Fill();
      }
    }
          if(k%100000==0)cout<<"Event Fill : "<<k<<" / "<<ENum<<endl;
  }// End Fill
  
  

}
  
//////////////////////////////////////////////////////////////////////////////

void coincalib::Write(){

   tnew->Write();
   hcoin_cut->Write();

}

//////////////////////////////////////////////////////////////////////////////

void coincalib::Close(){
  
  ofr->Close();
  
}

//////////////////////////////////////////////////////////////////////////////

//=======================================================================//
//===============   Main   ==============================================//
//=======================================================================//


int main(int argc, char** argv){



  gStyle->SetOptFit(111111111);
  int ch; char* mode="C";
  //  double tdc_time=58.0e-3;//[ns]
  string ifname;
  string ofname;
  string matrix_name="./matrix /momcalib_matrix.dat";
  string ofMTPname="./matrix/momcalib_test.dat";
  string opt_file="./scale_offset_20190210.dat";  
  string iteration;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = false;
  bool root_flag=false;
  bool matrix_flag=false;
  // bool RHRS_flag=true;
  bool RHRS_flag=false; 
  bool tuning_flag=false;
  string pngname;
  extern char *optarg;
  char* Target="H";
  string F1tdc="1";
  int f1tdc=1;
  //  char *root_init="/w/halla-scifs17exp/triton/itabashi/rootfiles/calib_root/";//ifarm
  string root_init="../rootfiles/";
  string root_end=".root";
  string dat_init="../matrix/";
  string dat_end=".dat";
  string matrix="matrix/zt_RHRS_2.dat";
  string pname;
  string mtparam_R;
  string mtparam_L;
  while((ch=getopt(argc,argv,"h:s:w:t:p:f:r:z:L:R:m:o:O:i:bcoC"))!=-1){
    switch(ch){
      
      
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;

    case 's':
      root_flag = true;
      draw_flag = false;
      single=true;
      ifname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;      

      
    case 'r':
      root_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;

  
    case 'O':
      F1tdc= optarg;
      f1tdc=atoi(F1tdc.c_str());
      //if F1tdc=1 tdc resolution 0.056
      //if F1tdc=2 tdc resolution 0.058
      //if F1tdc=3 tdc resolution 0.058 & Lp scale ON
      cout<<"F1 resolution mode : "<<f1tdc<<endl;      
      break;

    case 'p':
      pname = optarg;
      break;
      
    case 'w':
      ofMTPname  = optarg;
      cout<<"output new parameter filename : "<<ofMTPname<<endl;      
      matrix_flag=true;
      break;

    case 'R':
      MODE=MODE-1;
      mtparam_R =optarg;
      matrix_flag=true;
      break;

    case 'L':
      mtparam_L =optarg;
      matrix_flag=true;
      MODE=MODE+1;

    case 'C':
      MODE=0;      
      
    case 'i':
      iteration =optarg;
      tuning_flag=true;
      nite=atoi(iteration.c_str());
      cout<<"#iteration : "<<nite<<endl;
      break;
      
    case 'o':
      root_flag = true;
      draw_flag = false;
      matrix_flag=true;
      ofname = optarg;
      ofname =root_init+ ofname+ root_end;
      ofMTPname=optarg;
      ofMTPname = dat_init + ofMTPname;//+ dat_end;

      cout<<"output root filename : "<<ofname<<endl;      
      cout<<"output new parameter filename : "<<ofMTPname + ".dat"<<endl;      
      break;
      
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;

    case 'm':
      matrix_name =optarg;
      matrix_flag=true;
      break;
      


    case 'h':
      cout<<"-f : input root  filename"<<endl;
      cout<<"-p : input matrix filename"<<endl;      
      cout<<"-w : output matrix filename"<<endl;
      cout<<"-r : output root filename"<<endl;
      cout<<"-t : target mode H:hydrogen T:trititum"<<endl;
      cout<<"-m : tuning mode C:R&L-HRS, L:L-HRS, R:R-HRS "<<endl;      
      cout<<"-o : output root & tuning file name "<<endl;
      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;

    }
  }




  
  TApplication *theApp =new TApplication("App",&argc,argv);
  gSystem->Load("libMinuit");
  coincalib* Coin=new coincalib();
  if(MODE>=0)Coin->MTParam_L(mtparam_L);
  if(MODE<=0)Coin->MTParam_R(mtparam_R);
  if(single)Coin->SetRoot(ifname);
  else Coin->SetRunList(ifname);
  Coin->NewRoot(ofname);
  Coin->ReadParam(pname);
  Coin->EventSelect();
  if(matrix_flag)Coin->CoinTuning(matrix_name, MODE);
  Coin->Fill();
  Coin->Write();
  Coin->Close();

  cout<<"=====================================" <<endl;
  cout<<"========== Comment Out ==============" <<endl;
  cout<<"=====================================" <<endl;
  cout<<endl;
  cout<<"new root : "<<ofname.c_str()<<endl;
  cout<<"new matrix : "<<ofname.c_str()<<endl;
  
  gSystem->Exit(1);
  theApp->Run();

  return 0;
  
}//end main





//=======================================================================//
//=============== Function ==============================================//
//=======================================================================//

// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
{
  const double sigma=0.2;//[ps]
  double ct_c;
  double rtof_f;
  double ltof_f;
  double Chi2=0.0;
  double chi2;

  double par_R[nParamTp];
  double par_L[nParamTp];

  for(int i=0;i<nParamTc;i++){
    if(mode==0){
      par_R[i]=param[i];
      par_L[i]=param[i+nParamTc];
      
    }else if(mode==-1){
      par_R[i]=param[i];
      par_L[i]=Pct_L[i];
    }else if(mode==1){
      par_R[i]=Pct_R[i];
      par_L[i]=param[i];
    }
  }
  //   for(int i=0;i<nParamTc;i++)cout<<"i "<<i<<" Pct_R "<<par_R[i]<<endl;


  for(int k=0;k<ntune_event;k++){
    ct_c=-1000.;
    rtof_f=-1000.;
    ltof_f=-1000.;
    chi2=10000;


 

    
    
    //==== Calibration ====//

      rx_fp[k]  = (rx_fp[k]-XFPm)/XFPr;
      rth_fp[k] = (rth_fp[k]-XpFPm)/XpFPr;
      ry_fp[k]  = (ry_fp[k]-YFPm)/YFPr;
      rph_fp[k] = (rph_fp[k]-YpFPm)/YpFPr;

      lx_fp[k]  = (lx_fp[k]-XFPm)/XFPr;
      lth_fp[k] = (lth_fp[k]-XpFPm)/XpFPr;
      ly_fp[k]  = (ly_fp[k]-YFPm)/YFPr;
      lph_fp[k] = (lph_fp[k]-YpFPm)/YpFPr;	        

      rtof_f=RTOF[k];
      ltof_f=LTOF[k];
      if(mode<=0)rtof_f+=rtof_f*calcf_ct(par_R,rx_fp[k],rth_fp[k],ry_fp[k],rph_fp[k],1.0);
      if(mode>=0)ltof_f+=ltof_f*calcf_ct(par_L,lx_fp[k],lth_fp[k],ly_fp[k],lph_fp[k],1.0);


    ct_c=rtof_f-ltof_f;

    chi2=pow((ct[k]-ct_c)/sigma,2) / ntune_event;
    //    if(fabs(chi2)<10000)Chi2+=chi2;
    //    else Chi2+=10000;
    Chi2+=chi2;
    //    cout<<"coin "<<ct_c<<" chi2 "<<chi2<<" rotf_f "<<rtof_f<<endl;    
  }
    fval=Chi2;
};


// #################################################
double calcf_ct(double* P, double xf, double xpf, 
		double yf, double ypf, double zt){
// ###############################################



  const int nMatT=nnc; 
  const int nXf=nnc;
  const int nXpf=nnc;
  //  const int nYf=nnc;
  //  const int nYpf=nnc;
  const int nYf=0;
  const int nYpf=0;
  const int nZ=0;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
  for (int n=0;n<nMatT+1;n++){
    for (e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){		
		if (a+b+c+d+e==n){
		  if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf){
		    x = pow(xf,double(a))*pow(xpf,double(b))*
		      pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		  }
		  else{
		    x = 0.;
		  }
		  Y += x*P[npar];
		  npar++;
		}
		
	    }
	  }
	}
      }
    }
  }
  return Y;
  }

// #################################################
double calcf_ct2(double* P, double xf, double xpf, 
		double yf, double ypf, double zt){
// ###############################################


  return P[0]*xf+P1[1];




}
