using namespace std;
#include "coincalib.h"
#include "Param.h"
#include <TMinuit.h>
extern void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
extern double calcf_path(double* P, double xf, double xpf, double yf, double ypf,double zt);
extern double calcf_x(double* P, double xf);

bool DEBUG=true;

  const int nMatT=nnc; 
  const int nXf=nnc;
  const int nXpf=nnc;
  const int nYf=nnc;
  const int nYpf=nnc;
  const int nZt=nnc;


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

   double par=0.;
    int px=0,pth=0,py=0,pph=0,pz=0;
    cout<<"nParam "<<nParamTc<<endl;
   for(int i=0;i<nParamTc;i++){
    par=0.;
    px=0,pth=0,py=0,pph=0,pz=0;
    Mpt >> par >> px >> pth >> py >> pph >> pz;
     Pct[i]=par;
     Pct_R[i]=par;
    //cout<<"i "<<i<<" par "<<par<<endl;
    
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
    int px=0,pth=0,py=0,pph=0,pz=0;
    Mpt_L >> par >> px >> pth >> py >> pph >> pz;
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
  tnew->Branch("rpathl",&R_pathl);
  tnew->Branch("lpathl",&L_pathl);  
  tnew->Branch("rpathl_c",&R_pathl_c);
  tnew->Branch("lpathl_c",&L_pathl_c);
  tnew->Branch("Rp",&Rp);
  tnew->Branch("Beta_R",&Beta_R);
  tnew->Branch("Beta_L",&Beta_L);
  hcoin_cut= new TH1D("hcoin_cut","Coin with cut ",100,-20.,20.);
  hcoin_cut->SetLineColor(2);
  hcoin_cut->SetFillColor(2);
  hcoin_cut->SetFillStyle(3002);
  hcoin    = new TH1D("hcoin","Coin ",100,-20.,20.);
}



////////////////////////////////////////////////////////////////////////////


void coincalib::CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit){



  
  convertF1TDCR(param);
  convertF1TDCL(param);
  
  PathCalc(rhit,lhit);
  
  Beta_R = R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+MK*MK);
  Beta_L = L_tr_p[lhit]/sqrt(L_tr_p[lhit]*L_tr_p[lhit]+Me*Me);
  
  //====== w/o Path Calibration =========//  
  Rtof = RS2_F1time[RS2_seg] - R_pathl/(Beta_R*LightVelocity);
  Ltof = LS2_F1time[LS2_seg] - L_pathl/(Beta_L*LightVelocity);
  
  //====== w/ Path Calibration =========// 
  Rtof_c = RS2_F1time[RS2_seg] - R_pathl/(Beta_R*LightVelocity);
  Ltof_c = LS2_F1time[LS2_seg] - L_pathl/(Beta_L*LightVelocity); 
  
  if(mode<=0) Rtof_c = RS2_F1time[RS2_seg] - R_pathl_c/(Beta_R*LightVelocity);
  if(mode>=0) Ltof_c = LS2_F1time[LS2_seg] - L_pathl_c/(Beta_L*LightVelocity);  
 
  
  Rs2_t=RS2_F1time[RS2_seg];
  Ls2_t=RS2_F1time[LS2_seg];


  if(RS2_F1time[RS2_seg]!=-9999. &&LS2_F1time[LS2_seg]!=-9999.){
    coint = - Rtof + Ltof - coin_offset;
    coint_c = - Rtof_c + Ltof_c - coin_offset;

  }
  else{
    coint=-1000;
    coint_c=-1000;
  }

  //  cout<<" ct "<<coint<<" Rtof "<<Rtof<<" Ltof "<<Ltof<<" Rpath "<<R_pathl<<" Lpath "<<L_pathl<<endl;  
  //  cout<<" Rbeta "<<Beta_R<<" Lbeta "<<Beta_L<<" RS2t "<<Rs2_t<<" Ls2t "<<Ls2_t<<endl;  

  
}
///////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////
/*

double coincalib::CoinCalc_org(int RS2_seg, int LS2_seg, int rhit, int lhit){

  double cointime=0.0;

  convertF1TDCR(param);
  convertF1TDCL(param);

  double Rpathl=R_tr_pathl[rhit]+R_s2_trpath[rhit];
  double Lpathl=L_tr_pathl[lhit]+L_s2_trpath[lhit];
  double Beta_R=R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+MK*MK);
  double Beta_L=L_tr_p[lhit]/sqrt(L_tr_p[lhit]*L_tr_p[lhit]+Me*Me);
  double tof_r=RS2_F1time[RS2_seg] - Rpathl/(Beta_R*LightVelocity);
  double tof_l=LS2_F1time[LS2_seg] - Lpathl/(Beta_L*LightVelocity);
  Rs2_t=RS2_F1time[RS2_seg];
  Ls2_t=RS2_F1time[LS2_seg];
  Rtof=tof_r;
  Ltof=tof_l;
  beta_r=Beta_R;
  beta_l=Beta_L;
  
  if(RS2_F1time[RS2_seg]!=-9999. &&LS2_F1time[LS2_seg]!=-9999.){
    cointime= - tof_r + tof_l - coin_offset;

  }
  else{
    cointime=-1000;
  }
  
  return cointime;
  
}


///////////////////////////////////////////////////////////////////////////

double coincalib::CoinCalc_c(int RS2_seg, int LS2_seg, int rhit, int lhit){

  double cointime=0.0;

  convertF1TDCR(param);
  convertF1TDCL(param);

  
  double Rpathl=R_tr_pathl[rhit]+R_s2_trpath[rhit];
  double Lpathl=L_tr_pathl[lhit]+L_s2_trpath[lhit];

  R_pathl=Rpathl;
  L_pathl=Lpathl;  




      
  
  PathCalc(rhit,lhit);

   
  
  double Beta_R=R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+MK*MK);
  double Beta_L=L_tr_p[lhit]/sqrt(L_tr_p[lhit]*L_tr_p[lhit]+Me*Me);
  double tof_r=RS2_F1time[RS2_seg] - Rpathl/(Beta_R*LightVelocity);
  double tof_l=LS2_F1time[LS2_seg] - Lpathl/(Beta_L*LightVelocity);
  Rs2_t=RS2_F1time[RS2_seg];
  Ls2_t=RS2_F1time[LS2_seg];

  R_pathl_c=Rpathl;
  L_pathl_c=Lpathl;
  
  Rtof_c=tof_r;
  Ltof_c=tof_l;

  
  if(RS2_F1time[RS2_seg]!=-9999. &&LS2_F1time[LS2_seg]!=-9999.){
    cointime= - tof_r + tof_l - coin_offset;

  }
  else{
    cointime=-1000;
  }
  
  return cointime;
  
}

*/


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



  cout<<"mode "<<MODE<<" allParam "<<allparam<<endl;


  TMinuit* minuit= new TMinuit(allparam);
  minuit->SetFCN(fcn); // fcn Chi-square function

  
  double start[allparam];
  double step[allparam];

  /*
  const int nMatT =nnc;  
  const int nXf   =nnc;
  const int nXpf  =nnc;
  const int nYf   =nnc;
  const int nYpf  =nnc;
  const int nZt   =nnc;
  */


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
		    //step[npar] = 0.1;
		    //		    cout<<"n "<<npar<<" start "<<start[npar]<<" step "<<step[npar]<<endl;
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
  //arglist[0] = 2.0; // original
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




//////////////////////////////////////////////////////////////////////////////

void coincalib::EventSelect(){

  cout<<"==========================="<<endl;
  cout<<"===== Event Selection ======"<<endl;
  cout<<"==========================="<<endl;

  cout<<"Events : "<<ENum<<endl;

  for(int k=0;k<ENum;k++){

    // Initialization //
    for(int i=0;i<MAX;i++){
    R_tr_x[i]  = -2000.0;
    R_tr_th[i] = -2000.0;
    R_tr_y[i]  = -2000.0;
    R_tr_ph[i] = -2000.0;
    R_tr_p[i]  = -2000.0;
    L_tr_x[i]  = -2000.0;
    L_tr_th[i] = -2000.0;
    L_tr_y[i]  = -2000.0;
    L_tr_ph[i] = -2000.0;
    L_tr_p[i]  = -2000.0;
    }
    Rtof=-1000.0; Rtof_c=-1000.; Ltof=-1000.; Ltof_c=-1000.;
    coint=-1000.0; coint_c=-1000.0;
    R_pathl=-1000.0; R_pathl_c=-1000; L_pathl=-1000.; L_pathl_c=-1000;
    Rp=-1000.;



    
    tree->GetEntry(k);

    // ===== Initialization ====== //    

    z_flag=false;
    pid_flag=false;
    ct_flag=false;
    ct_flag=false;
    fp_flag=false;
    pion_flag=false;
    kaon_flag=false;

    
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
	  if(fabs(R_tr_x[rt])<10. && fabs(R_tr_th[rt])<10. && fabs(R_tr_y[rt])<10. && fabs(R_tr_ph[rt])<10.
	     && fabs(L_tr_x[rt])<10. && fabs(L_tr_th[rt])<10. && fabs(L_tr_y[rt])<10. && fabs(L_tr_ph[rt])<10.)
	    fp_flag=true;
	  
	  int R_s2pad=(int)R_s2_t_pads[rt];
	  int L_s2pad=(int)L_s2_t_pads[lt]; 
	  CoinCalc(R_s2pad,L_s2pad,rt,lt);
	  //cout<<"k "<<k<<" rt "<<rt<<" lt "<<lt<<" Rseg "<<R_s2pad<<" Lseg "<<L_s2pad<<" coint "<<coint<<endl;
	  if( fabs(coint )<1.0 )kaon_flag=true;
	  //if( fabs( CoinCalc(R_s2pad,L_s2pad,rt,lt) -3)<1.0 )pion_flag=true;

	  if(ntune_event<nmax && z_flag && pid_flag && fp_flag)
	    hcoin->Fill(coint);

	  
	  if(ntune_event<nmax && z_flag && pid_flag && fp_flag && (kaon_flag || pion_flag)){
	    ct[ntune_event]    = coint;
	    rx_fp[ntune_event] = R_tr_x[rt];
	    ry_fp[ntune_event] = R_tr_y[rt];
	    rth_fp[ntune_event] = R_tr_th[rt];
	    rph_fp[ntune_event] = R_tr_ph[rt];
	    rz[ntune_event]=R_tr_vz[rt];
	    lx_fp[ntune_event] = L_tr_x[lt];
	    ly_fp[ntune_event] = L_tr_y[lt];
	    lth_fp[ntune_event] = L_tr_th[lt];
	    lph_fp[ntune_event] = L_tr_ph[lt];
	    lz[ntune_event]=L_tr_vz[lt];
	    rpathl[ntune_event]= R_tr_pathl[rt] + R_s2_trpath[rt];
	    lpathl[ntune_event]= R_tr_pathl[lt] + R_s2_trpath[lt];
	    rs2_t[ntune_event]=Rs2_t;
	    ls2_t[ntune_event]=Ls2_t;
	    rp[ntune_event]=R_tr_p[rt];
	    lp[ntune_event]=L_tr_p[lt];
	    if(kaon_flag){
	      Mass[ntune_event]=MK;
	      Kaon_nev+=1.0;            }
	    if(pion_flag){
	      Mass[ntune_event]=Mpi;
	      Pion_nev+=1.0;            }
	    
	    beta_K[ntune_event]=Beta_R;
	    beta_e[ntune_event]=Beta_L;
	    RTOF[ntune_event]=Rtof;
	    LTOF[ntune_event]=Ltof;
	    hcoin_cut->Fill(ct[ntune_event]);
	    rs2_seg[ntune_event]=R_s2pad;
	    ls2_seg[ntune_event]=L_s2pad;
	    
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
  cout<<"Pion selection : "<<Pion_nev<<" Kaon selection : "<<Kaon_nev<<endl;  
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



////////////////////////////////////////////////////////////////////////////

void coincalib::PathCalc(int rhit, int lhit){


  R_pathl= R_tr_pathl[rhit] + R_s2_trpath[rhit];
  L_pathl= L_tr_pathl[lhit] + L_s2_trpath[lhit];
  
  
  R_tr_x[rhit]  = (R_tr_x[rhit]-XFPm)/XFPr;
  R_tr_th[rhit] = (R_tr_th[rhit]-XpFPm)/XpFPr;
  R_tr_y[rhit]  = (R_tr_y[rhit]-YFPm)/YFPr;
  R_tr_ph[rhit] = (R_tr_ph[rhit]-YpFPm)/YpFPr;
  R_tr_vz[rhit] = (R_tr_vz[rhit] - Ztm)/Ztr;
  
  L_tr_x[lhit]  = (L_tr_x[lhit]-XFPm)/XFPr;
  L_tr_th[lhit] = (L_tr_th[lhit]-XpFPm)/XpFPr;
  L_tr_y[lhit]  = (L_tr_y[lhit]-YFPm)/YFPr;
  L_tr_ph[lhit] = (L_tr_ph[lhit]-YpFPm)/YpFPr;
  L_tr_vz[lhit] = (L_tr_vz[lhit] - Ztm)/Ztr;      
  
  
  
  //==== Calc Path Length =========//
  
      
  
  R_pathl_c = calcf_path(Pct_R,R_tr_x[rhit],R_tr_th[rhit],R_tr_y[rhit],R_tr_ph[rhit],R_tr_vz[rhit]); // ns
  
  L_pathl_c = calcf_path(Pct_L,L_tr_x[lhit],L_tr_th[lhit],L_tr_y[lhit],L_tr_ph[lhit],L_tr_vz[lhit]); // ns
  
      
  R_tr_x[rhit]  = R_tr_x[rhit] * XFPr + XFPm;
  R_tr_th[rhit] = R_tr_th[rhit] * XpFPr + XpFPm;
  R_tr_y[rhit]  = R_tr_y[rhit] * YFPr + YFPm;
  R_tr_ph[rhit] = R_tr_ph[rhit] * YpFPr   + YpFPm;
  R_tr_vz[rhit] = R_tr_vz[rhit] * Ztr + Ztm;
  
  L_tr_x[lhit]  = L_tr_x[lhit] * XFPr + XFPm;
  L_tr_th[lhit] = L_tr_th[lhit] * XpFPr + XpFPm;
  L_tr_y[lhit]  = L_tr_y[lhit] * YFPr + YFPm;
  L_tr_ph[lhit] = L_tr_ph[lhit] * YpFPr   + YpFPm;
  L_tr_vz[lhit] = L_tr_vz[lhit] * Ztr + Ztm;
  
  R_pathl_c = R_pathl_c * PaRr + PaRm;
  L_pathl_c = L_pathl_c * PaLr + PaLm;


  
      
}

///////////////////////////////////////////////////////////////////////////

void coincalib::Fill(){

  cout<<"=================================="<<endl;
  cout<<"======== Fill Events ============="<<endl;
  cout<<"=================================="<<endl;
  double tof_r,tof_l;
  for(int k;k<ENum;k++){

    // Initialization //
    for(int i=0;i<MAX;i++){
    R_tr_x[i]  = -2000.0;
    R_tr_th[i] = -2000.0;
    R_tr_y[i]  = -2000.0;
    R_tr_ph[i] = -2000.0;
    R_tr_p[i]  = -2000.0;
    L_tr_x[i]  = -2000.0;
    L_tr_th[i] = -2000.0;
    L_tr_y[i]  = -2000.0;
    L_tr_ph[i] = -2000.0;
    L_tr_p[i]  = -2000.0;
    }
    Rtof=-1000.0; Rtof_c=-1000.;Ltof=-1000., Ltof_c=-1000.;
    coint=-1000.0; coint_c=-1000.0;
    R_pathl=-1000.0; R_pathl_c=-1000.; L_pathl=-1000.; L_pathl_c=-1000.;
    Rp=-1000.;

    
    tree->GetEntry(k);


    int NLtr = (int)L_tr_n;  if(NLtr>MAX) NLtr = MAX;
    int NRtr = (int)R_tr_n;  if(NRtr>MAX) NRtr = MAX;

    for(int lt=0;lt<NLtr;lt++){
      for(int rt=0;rt<NRtr;rt++){

	int R_s2pad=(int)R_s2_t_pads[rt];
	int L_s2pad=(int)L_s2_t_pads[lt]; 

	Rp=R_tr_p[rt];
	
	CoinCalc(R_s2pad, L_s2pad, rt, lt);	
	tnew->Fill();
	
      }
    }
    if(k%100000==0)cout<<"Event Fill : "<<k<<" / "<<ENum<<endl;
  }// End Fill
  
  
  
}
  
//////////////////////////////////////////////////////////////////////////////

void coincalib::Write(){

   tnew->Write();
   hcoin->Write();
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
  cout<<"new matrix : "<<matrix_name.c_str()<<endl;
  
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
  //  double ct;
  double ct_c;
  double rtof_f;
  double ltof_f;
  double Chi2=0.0;
  double chi2;
  double par_R[nParamTc];
  double par_L[nParamTc];
  double beta_R;
  double beta_L;
  double weight=Pion_nev/Kaon_nev;
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
    chi2=100;
    //    rpathl[k]=0.0;
    //    lpathl[k]=0.0;
    beta_R=0.0;
 

    
    
    //==== Calibration ====//
    //       cout<<"rx_fp "<<rx_fp[k]<<" rth_fp "<<rth_fp[k]<<" ry_fp "<<ry_fp[k]<<" rph_fp "<<rph_fp[k]<<endl;
    //    cout<<"k "<<k<<" lpathl "<<lpathl[k]<<" PaLm "<<PaLm<<" PaLr "<<PaLr<<endl;
            
      rx_fp[k]  = (rx_fp[k]-XFPm)/XFPr;
      rth_fp[k] = (rth_fp[k]-XpFPm)/XpFPr;
      ry_fp[k]  = (ry_fp[k]-YFPm)/YFPr;
      rph_fp[k] = (rph_fp[k]-YpFPm)/YpFPr;
      rz[k]     = (rz[k] - Ztm)/Ztr;
      lx_fp[k]  = (lx_fp[k]-XFPm)/XFPr;
      lth_fp[k] = (lth_fp[k]-XpFPm)/XpFPr;
      ly_fp[k]  = (ly_fp[k]-YFPm)/YFPr;
      lph_fp[k] = (lph_fp[k]-YpFPm)/YpFPr;
      lz[k]     = (lz[k] - Ztm)/Ztr;
      rpathl[k] = (rpathl[k] - PaRm )/PaRr;
      lpathl[k] = (lpathl[k] - PaLm )/PaLr;

      beta_R = rp[k]/sqrt(Mass[k]*Mass[k] + rp[k]*rp[k]);
      beta_L = beta_e[k];

      //      cout<<" betaR "<<beta_R<<" beta_L "<<beta_L<<" Mass "<<Mass[k]<<endl;
 
      if(mode<=0)
 	rpathl[k]= calcf_path(par_R,rx_fp[k],rth_fp[k],ry_fp[k],rph_fp[k],rz[k]); // ns
      if(mode>=0)
	lpathl[k]= calcf_path(par_L,lx_fp[k],lth_fp[k],ly_fp[k],lph_fp[k],lz[k]); // ns



      rx_fp[k]  = rx_fp[k] * XFPr + XFPm;
      rth_fp[k] = rth_fp[k] * XpFPr + XpFPm;
      ry_fp[k]  = ry_fp[k] * YFPr + YFPm;
      rph_fp[k] = rph_fp[k] * YpFPr   + YpFPm;
      rz[k]     = rz[k] * Ztr + Ztm;      
      lx_fp[k]  = lx_fp[k] * XFPr + XFPm;
      lth_fp[k] = lth_fp[k] * XpFPr + XpFPm;
      ly_fp[k]  = ly_fp[k] * YFPr + YFPm;
      lph_fp[k] = lph_fp[k] * YpFPr   + YpFPm;
      lz[k]     = lz[k] * Ztr + Ztm;
      
      rpathl[k] = rpathl[k] * PaRr + PaRm;
      lpathl[k] = lpathl[k] * PaLr + PaLm;


      rtof_f = rs2_t[k] - rpathl[k]/(beta_R*LightVelocity); 
      ltof_f = ls2_t[k] - lpathl[k]/(beta_L*LightVelocity);      
	
      //      ct   = RTOF[k] - LTOF[k];
      ct_c = rtof_f  - ltof_f;
      //      if(fabs(ct_c)>1000.)ct_c=1000.;
      
      if(Mass[k]>0.5 && Pion_nev>10)        chi2= weight * pow((ct_c)/sigma,2) / ntune_event;
      else if(Mass[k]>0.5 && Pion_nev<10) chi2= pow((ct_c)/sigma,2) / ntune_event;
      else  chi2=pow((ct_c-3.0)/sigma,2) / ntune_event;
      //   chi2=fabs(ct-ct_c);
	
    Chi2+=chi2;
    //    cout<<" k "<<k<<" rtof "<<rtof_f<<" ltof"<<ltof_f<<" ct_c "<<ct_c<<" chi2 "<<chi2<<" rpathl "<<rpathl[k]<<endl;
  }


  fval=Chi2;
  //  cout<<"end fcn "<<endl;
};


// #################################################
double calcf_path(double* P, double xf, double xpf, 
		double yf, double ypf, double zt){
// ###############################################


  
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
double calcf_x(double* P, double xf){
// ###############################################

  return P[0]*xf+P[1];




}