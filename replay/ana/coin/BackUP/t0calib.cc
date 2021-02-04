#include "Param.h"
#include "t0calib.h"
#include <TMinuit.h>

int ntuned_event;
int niter=0;
int NMAX =100000;
bool rhrs ;
extern void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);
extern void fcn2(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/);

//////////////////////////////////////////////////////////////


void t0calib::SetRoot(string ifname){

  add_tree(ifname);
  pack_tree();
  readtreeHRSR();
  readtreeHRSL();
  
  tree->SetBranchStatus("LTDC.F1FirstHit",1);
  tree->SetBranchAddress("LTDC.F1FirstHit",LF1);
  tree->SetBranchStatus("RTDC.F1FirstHit",1);
  tree->SetBranchAddress("RTDC.F1FirstHit",RF1);  
  tree->SetBranchStatus("Rs2_pad",1);
  tree->SetBranchAddress("Rs2_pad",&Rs2_pad);
  tree->SetBranchStatus("Ls2_pad",1);
  tree->SetBranchAddress("Ls2_pad",&Ls2_pad);
  tree->SetBranchStatus("Rpathl",1);
  tree->SetBranchAddress("Rpathl",&Rpathl);
  tree->SetBranchStatus("Lpathl",1);
  tree->SetBranchAddress("Lpathl",&Lpathl);
  tree->SetBranchStatus("Rp",1);
  tree->SetBranchAddress("Rp",&momR);
  tree->SetBranchStatus("Lp",1);
  tree->SetBranchAddress("Lp",&momL);
  tree->SetBranchStatus("ct_c",1);
  tree->SetBranchAddress("ct_c",&ct);
}


///////////////////////////////////////////////////////////////



void t0calib::EventSelection(){

  cout<<endl;
  cout<<"==============================================="<<endl;
  cout<<"===========< Event Selection >================="<<endl;
  cout<<"================================================"<<endl;
  
  int ENum = tree->GetEntries();
  if(ENum> NMAX)  ENum = NMAX;
  cout<<"ENum "<<ENum<<endl;
  double LS2_off[16]={2.98,2.6,1.255,2.15,0.282,0.75,0.51,-0.92,-1.41,-0.33,-0.73,3.3,1.916, 1.278,-0.71,0.0};
  int ituned =0;
  double tdc_resolution ;

  
  for(int nev=0;nev<ENum;nev++){

    tree->GetEntry(nev);
    
    ev.Rpathl = Rpathl;
    ev.Lpathl = Lpathl;
    ev.Rp = momR;
    ev.Lp = momL;
    ev.RS2_seg = Rs2_pad;
    ev.LS2_seg = Ls2_pad;
    ev.RF1_ref = RF1[9];
    ev.RF1_t   = RF1[ev.RS2_seg+16];
    ev.RF1_b   = RF1[ev.RS2_seg+48];
    ev.LF1_ref = LF1[40]; // gagami
    //    ev.LF1_ref = LF1[30]; // ana_Lambda
    ev.LF1_ref = 0.0;
    ev.LF1_t   = LF1[ev.LS2_seg];
    ev.LF1_b   = LF1[ev.LS2_seg+48];

    

    bool coin_cut =false;

    if(runnum < 111369) tdc_resolution = 56.23e-3; // nsec
    else                tdc_resolution = 58.e-3; // nsec

    
    double beta_R =  ev.Rp/sqrt(ev.Rp*ev.Rp + Mpi*Mpi);
    double beta_L =  ev.Lp/sqrt(ev.Lp*ev.Lp + Me*Me  );
    double meantime_R =  ev.Rpathl/(beta_R*LightVelocity) - ( (ev.RF1_t +ev.RF1_b )/2.0  -  ev.RF1_ref)* tdc_resolution;
    double meantime_L =  ev.Lpathl/(beta_L*LightVelocity) -                                 ev.LF1_ref * tdc_resolution - LS2_off[ev.LS2_seg];
    double coin_offset = -441.;
    double ctime = meantime_R -meantime_L + coin_offset;

    
    if(fabs(ctime)<2.0) coin_cut = true;
    
    if(coin_cut && ituned<nmax && ev.RS2_seg ){

      rs2_pad[ituned] = ev.RS2_seg;
      ls2_pad[ituned] = ev.LS2_seg;
      rs2_pathl[ituned] = ev.Rpathl;
      ls2_pathl[ituned] = ev.Lpathl;
      beta_r[ituned]    =  ev.Rp/sqrt(ev.Rp*ev.Rp + Mpi*Mpi);
      beta_l[ituned]    =  ev.Lp/sqrt(ev.Lp*ev.Lp + Me*Me  );
      rs2_cor[ituned]   =  ev.Rpathl/(beta_r[ituned]*LightVelocity);
      ls2_cor[ituned]   =  ev.Lpathl/(beta_l[ituned]*LightVelocity);
      RF1_t[ituned]     =  ev.RF1_t;
      RF1_b[ituned]     =  ev.RF1_b;
      RF1_ref[ituned]   =  ev.RF1_ref;
      LF1_t[ituned]     =  ev.LF1_t;
      LF1_b[ituned]     =  ev.LF1_b;
      LF1_ref[ituned]   =  ev.LF1_ref;
      tdc_time[ituned]  =  tdc_resolution;
      event_num[ituned] =  nev;
      ituned ++;
      hct_select->Fill(ctime);
  
    }


    if(ev.RS2_seg==8)hct_LS2[ev.LS2_seg]->Fill(ctime);
    if(ev.LS2_seg==8)hct_RS2[ev.RS2_seg]->Fill(ctime);


    
    hct_ev->Fill(ctime);      






    
    if(nev%10000==0)cout<<"Selected Event "<<nev<<" / "<< ENum<<endl;    
  }// END for


  ntuned_event = ituned;

  cout<<"End Event Selection"<<endl;
  cout<<"ituned : "<<ituned<<" ENum "<<ENum<<endl;
}

///////////////////////////////////////////////////////////////////////////




void t0calib::SetHist(){

  double min_ct =-100.;
  double max_ct = 100.;
  int bin_ct = 1000;

  hct_ev = new TH1D("hct_ev","hist ",bin_ct,min_ct,max_ct);
  hct_ev->SetLineColor(1);
  hct_select = new TH1D("hct_select","hist ",bin_ct,min_ct,max_ct);
  hct_select->SetFillColor(2);
  hct_select->SetLineColor(2);
  hct_select->SetFillStyle(3002);

  hct_select_c = new TH1D("hct_select_c","hist ",bin_ct,min_ct,max_ct);
  hct_select_c->SetFillColor(2);
  hct_select_c->SetLineColor(2);
  hct_select_c->SetFillStyle(3002);  
  
  //  hct_c  = new TH1D("hct_c ","hist ",bin_ct,min_ct,max_ct);
  hct_c  = new TH1D("hct_c ","hist ",bin_ct,min_ct,max_ct);
  hct_c->SetLineColor(4);


  for(int i=0;i<nS2;i++){

    hct_RS2[i] = new TH1D(Form("hct_RS2_%d",i),"hist ",bin_ct,min_ct,max_ct);									       hct_LS2[i] = new TH1D(Form("hct_LS2_%d",i),"hist ",bin_ct,min_ct,max_ct);
    hct_RS2[i] ->SetLineColor(1);
    hct_LS2[i] ->SetLineColor(1);
    fct_RS2[i] =new TF1(Form("fct_RS2_%d",i),"gansn(0)",min_ct,max_ct);
    fct_LS2[i] =new TF1(Form("fct_LS2_%d",i),"gansn(0)",min_ct,max_ct);
    fct_RS2[i] ->SetLineColor(2);
    fct_LS2[i] ->SetLineColor(2);

    
  }
  
}

///////////////////////////////////////////////////////////////////////////

void t0calib::Fill(){

  cout<<endl;
  cout<<"==============================================="<<endl;
  cout<<"===========<   Event Fill    >================="<<endl;
  cout<<"================================================"<<endl;
  
  int ENum = tree->GetEntries();
  int ituned =0;
  if(ENum> NMAX)  ENum = NMAX;
  cout<<"ENum "<<ENum<<endl;
  double tdc_resolution ;

  double ls2_off[16]={2.98,2.6,1.255,2.15,0.282,0.75,0.51,-0.92,-1.41,-0.33,-0.73,3.3,1.916, 1.278,-0.71,0.0};

  double ls2_off[16]={};
  
  for(int nev=0;nev<ENum;nev++){

    tree->GetEntry(nev);
    
    tr.Rpathl = Rpathl;
    tr.Lpathl = Lpathl;
    tr.Rp = momR;
    tr.Lp = momL;
    tr.RS2_seg = Rs2_pad;
    tr.LS2_seg = Ls2_pad;
    tr.RF1_ref = RF1[9];
    tr.RF1_t   = RF1[tr.RS2_seg+16];
    tr.RF1_b   = RF1[tr.RS2_seg+48];
    tr.LF1_ref = LF1[40];
    tr.LF1_ref = 0.0;
    tr.LF1_t   = LF1[tr.LS2_seg];
    tr.LF1_b   = LF1[tr.LS2_seg+48];
    if(runnum < 111369) tdc_resolution = 56.23e-3; // nsec
    else                 tdc_resolution = 58.e-3; // nsec



    
    double time_Rt = (tr.RF1_t - par_off[tr.RS2_seg])*tdc_resolution - tr.RF1_ref * tdc_resolution;
    double time_Rb = (tr.RF1_b - par_off[tr.RS2_seg])*tdc_resolution - tr.RF1_ref * tdc_resolution;

    //    double time_Lt = (tr.LF1_t - par_off[tr.LS2_seg+16])*tdc_resolution - tr.LF1_ref * tdc_resolution;
    //    double time_Lb = (tr.LF1_b - par_off[tr.LS2_seg+16])*tdc_resolution - tr.LF1_ref * tdc_resolution;    

    double time_Lt = (- par_off[tr.LS2_seg+16])*tdc_resolution - tr.LF1_ref * tdc_resolution;
    double time_Lb = (- par_off[tr.LS2_seg+16])*tdc_resolution - tr.LF1_ref * tdc_resolution;    
    

    double Rs2_cor = tr.Rp/sqrt(tr.Rp*tr.Rp + Mpi*Mpi);
    double Ls2_cor = tr.Lp/sqrt(tr.Lp*tr.Lp + Me *Me );
    double coin_off = par_off[32];

    
    double  meantime_R = - (time_Rt + time_Rb)/2.0   + Rs2_cor;
    double  meantime_L = - (time_Lt + time_Lb)/2.0   + Ls2_cor;

    tr.ct_b = meantime_R - meantime_L + coin_off - par_off[tr.LS2_seg+16] + par_off[tr.RS2_seg];
    tr.ct   = meantime_R - meantime_L + coin_off;

    //    if(fabs(ct-3.)<2.)
    //      event_num[ituned] =  nev;


    
    if(event_num[ituned]==nev && ituned<ntuned_event){
      hct_select_c->Fill(tr.ct);
      ituned ++;
    }
      
    hct_c->Fill(tr.ct);
    T ->Fill();

      
    if(nev%10000==0)cout<<"Filled Event "<<nev<<" / "<< ENum<<endl;    
  }// END for


  T->Write();
  hct_c->Write();
  hct_ev->Write();
  hct_select->Write();
  hct_select_c->Write();
}


///////////////////////////////////////////////////////////////////////////

void t0calib::NewRoot(string ofrname){

  ofr = new TFile(ofrname.c_str(),"recreate");
  T = new TTree("T","T");
  T = tree->CloneTree(0);
  T->Branch("cointime",&tr.ct);
  T->Branch("cointime_b",&tr.ct_b);
  

}

///////////////////////////////////////////////////////////////////////////


void t0calib::Tuning(string ofname){

  ofstream * ofparam;
  niter =1;
  if (niter>0){
    cout << "======================================================" <<endl;
    cout << "=================  Tuning started ====================" <<endl;
    cout << "======================================================" <<endl;}
  
  char tempc[500];
  const  char* new_tempc=ofname.c_str();
  cout<<"new marix file: "<<new_tempc<<endl;

  double chi_sq[100];
  niter =2;
  for(int i =0;i<niter;i++){

    sprintf(tempc,"%s_%d.dat",ofname.c_str(),i);
    ofparam = new ofstream(tempc);
    //    rhrs =false;
    //    if(i%2==0)rhrs =true;
    //       else rhrs =false;
    // rhrs =true; // RS2 T0 calibraiton mode
      chi_sq[i] = tune(par_off, i);

    cout<<"i "<<i<<" / "<<niter<<" chi "<<chi_sq[i]<<endl;
    *ofparam <<"# RS2 flag : seg : offset "<<endl;
    for(int j=0;j<nParamTc_off;j++){
      //      cout<<"j "<<j<<" par_off "<<par_off[j]<<endl;
      if(j<16) *ofparam <<1<<" "<<j<<" "<<par_off[j]<<endl;
      else if(16<=j && j <32) *ofparam <<0<<" "<<j-16<<" "<<par_off[j]<<endl;
    }
    *ofparam<<"coin offset "<<par_off[32]<<endl;
;
    ofparam -> close();
    //    ofparam -> clear();
  } // end niter

  
}

////////////////////////////////////////////////////////////////////////////

double t0calib::tune(double *pa, int j){

  double chi =0.0;
  double arglist[10];
  int ierflg =0;
  int allparam = nParamTc_off+1;

  TMinuit* minuit = new TMinuit(allparam);
  minuit -> SetFCN(fcn);

  double start[allparam];
  double step[allparam];

  
  arglist[0] =1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  minuit -> SetPrintLevel(-1);

  double LLim[allparam],ULim[allparam];
  char pname[500];

  for(int i=0;i<allparam;i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i];

    step[i]  = 0.01;
    LLim[i]  = pa[i] - 1000.0;
    ULim[i]  = pa[i] + 1000.0;
    if(i==nParamTc_off){
      start[i]= -441.;
    step[i]  = 0.01;
    LLim[i]  = pa[i] - 1000.0;
    ULim[i]  = pa[i] + 1000.0;
    }
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
   
  }

  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0; // original
  //  arglist[0] = 1.0; // test
  //arglist[0] = 0.0;   // test
  minuit->mnexcm("SET STR",arglist,1,ierflg);

  // ~~~~ Migrad + Simplex  ~~~~
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg); // Chi-square minimizatio\
n

  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double er;

  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);

  if(amin>0) chi=amin;


  for(int i=0 ; i<allparam ; i++){ 
    minuit -> GetParameter(i,pa[i],er);
    minuit -> GetParameter(i,par_off[i],er);
    //    cout<<"i "<<par_off[i]<<endl;
    //    par_off[i]=0.0; //test
  }


  return chi;
  
  
}



////////////////////////////////////////////////////////////////////////////



//===================================================================//
//============================= Main ================================//
//==================================================================//

int main(int argc, char** argv){


  int ch;
  extern char *optarg;
  string ifname;
  string ofrname;
  string ofpname;
  while((ch=getopt(argc,argv,"h:s:w:W:t:p:f:n:r:AlRILbcop"))!=-1){
      switch(ch){
	
      case 'f':
	ifname = optarg;
	cout<<"input filename : "<<ifname<<endl;
	break;

      case 's':
	ifname = optarg;
	cout<<"input root filename : "<<ifname<<endl;
	break;

      case 'w':
	ofpname = optarg;
	cout<<"output param filename : "<<ofpname<<endl;
	break;	
	
      case 'r':
	ofrname = optarg;
	cout<<"output root filename : "<<ofrname<<endl;
	break;	
	

      }

  }

      TApplication *theApp =new TApplication("App",&argc,argv);
      gSystem->Load("libMinuit");

      t0calib* t0 =new t0calib();

      t0 -> SetRoot(ifname);
      t0 -> SetHist();
      t0 -> NewRoot(ofrname);
      //      t0 -> EventSelection();
      //      t0 -> Tuning(ofpname);
      t0 -> Fill();
      
      gSystem->Exit(1);
      theApp->Run();
      

      
  
  }















////////////////////////////////////////////////////////////////////////////

// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param,
	 int /*ilag*/)
// #############################################################
{


  double RS2_off[16],LS2_off[16];
  double Chi_R[16],Chi_L[16];
  int nev_R[16],nev_L[16];

  double ls2_off[16]={2.98,2.6,1.255,2.15,0.282,0.75,0.51,-0.92,-1.41,-0.33,-0.73,3.3,1.916, 1.278,-0.71,0.0};
  
  for(int i=0;i<16; i++){

    /*
    if(rhrs ){
      RS2_off[i] = param[i];
      LS2_off[i] = par_off[i+(int)(nParamTc_off/2)];
    }
    if(!rhrs){
      RS2_off[i] = par_off[i];
      LS2_off[i] = param[i+(int)(nParamTc_off/2)];
      }
    */

    RS2_off[i] = param[i];
    LS2_off[i] = param[i+(int)(nParamTc_off/2)];
    LS2_off[i] = ls2_off[i];
	  nev_R[i] =0;
	  nev_L[i] =0;
	  Chi_R[i] =0.0;
	  Chi_L[i] =0.0;
  }

  double  coin_offset= param[nParamTc_off];
  //double  coin_offset= -441.;
  double time_Rt,time_Rb,time_Lt,time_Lb;
  double meantime_R,meantime_L;
  int RS2_seg, LS2_seg;
  double ctime;
  double chi2 =0.0;
  double ctime_pi=0.0;
  //  double ctime_pi=3.0;
  double sigma = 0.3; // 300 ps

  
    for(int i =0; i< ntuned_event; i++){

      RS2_seg = rs2_pad[i];
      LS2_seg = ls2_pad[i];
      
      time_Rt = (RF1_t[i] - RF1_ref[i] ) * tdc_time[i];
      time_Rb = (RF1_b[i] - RF1_ref[i] ) * tdc_time[i]; 
      time_Lt = ( - LF1_ref[i] ) * tdc_time[i]; 
      time_Lb = ( - LF1_ref[i] ) * tdc_time[i];

      meantime_R = + rs2_cor[i]  - (time_Rt + time_Rb)/2.0 - RS2_off[RS2_seg]* tdc_time[i];
      meantime_L = + ls2_cor[i]  - (time_Lt + time_Lb)/2.0 - LS2_off[LS2_seg]* tdc_time[i];


      ctime = meantime_R - meantime_L + coin_offset;
      

      
      
      Chi_R[RS2_seg] +=  (ctime_pi - ctime)*(ctime_pi - ctime)/sigma/sigma;
      Chi_L[LS2_seg] +=  (ctime_pi - ctime)*(ctime_pi - ctime)/sigma/sigma;

      nev_R[RS2_seg]++;
      nev_L[LS2_seg]++;

      
      // chi2 += (ctime_pi - ctime)*(ctime_pi - ctime)/sigma/sigma/ntuned_event;      
      
    }


    for(int i=1;i<15;i++){
 
      if(nev_R[i]>0 &&  rhrs) chi2 += Chi_R[i]/(double)nev_R[i];
      if(nev_L[i]>0 && !rhrs) chi2 += Chi_L[i]/(double)nev_L[i];
 
      //      if(nev_R[i]>0 &&  rhrs) chi2 += Chi_R[i]/(double)ntuned_event;
      //      if(nev_L[i]>0 && !rhrs) chi2 += Chi_L[i]/(double)ntuned_event;      
 
    }
   
  


  fval = chi2;

}// end fcn





// #############################################################
void fcn2(int &nPar, double* /*grad*/, double &fval, double* param,
	 int /*ilag*/)
// #############################################################
{


  double RS2_off[16],LS2_off[16];
  double Chi_R[16],Chi_L[16];
  int nev_R[16],nev_L[16];



  if(rhrs) RS2_off[i] = param[i];
  else     LS2_off[i] = param[i];

  double  coin_offset= param[nParamTc_off];
  double time_Rt,time_Rb,time_Lt,time_Lb;
  double meantime_R,meantime_L;
  int RS2_seg, LS2_seg;
  double ctime;
  double chi2 =0.0;
  double ctime_pi=0.0;
  //  double ctime_pi=3.0;
  double sigma = 0.3; // 300 ps

  
    for(int i =0; i< ntuned_event; i++){

      RS2_seg = rs2_pad[i];
      LS2_seg = ls2_pad[i];
      
      time_Rt = (RF1_t[i] - RF1_ref[i] ) * tdc_time[i];
      time_Rb = (RF1_b[i] - RF1_ref[i] ) * tdc_time[i]; 
      time_Lt = ( - LF1_ref[i] ) * tdc_time[i]; 
      time_Lb = ( - LF1_ref[i] ) * tdc_time[i];

      meantime_R = + rs2_cor[i]  - (time_Rt + time_Rb)/2.0 - RS2_off[RS2_seg]* tdc_time[i];
      meantime_L = + ls2_cor[i]  - (time_Lt + time_Lb)/2.0 - LS2_off[LS2_seg]* tdc_time[i];


      ctime = meantime_R - meantime_L + coin_offset;

      chi2 += ctime*ctime/sigma/sigma/(double)ntuned_event;

  
    }

    


  fval = chi2;

}// end fcn




