#ifndef zcalib_h
#define zcalib_h 1
using namespace std;
#include "Setting.h"
#include "tree.h"
#include "Param.h"
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <TChain.h>
//#include <fstaream>
#include <iostream>
#include <TMinuit.h>

  TMinuit *minuit;
  double min_chi[100];
  int nite=0;
  double XFP, XpFP;
  double YFP, YpFP;




extern void fcn(int &nPar, double* /*grad*/, 
		double &fval, double* param, int /*iflag*/);
extern double tune(double* pa, int j);
extern double calcf2t_zt(double* P, 
			 double xf, double xpf,
			 double yf, double ypf);




//================================================//
//================== zdata class ================//
//===============================================//

class zcalib : public tree{
 public:
  zcalib();
  ~zcalib();


 public:
  void SetRoot(string ifname);
  void ChainRoot(string ifname);
  //  void SetRunList(string ifname);
  void Mzt(string matrix_name, bool rarm);  
  void Mzt_L(string matrix_name, bool rarm);
  //  void SetBranch(string ifname);
  void ZCorr();
  void NewBranch(string ofname, bool rarm);
  void MakeHist();
  void GetEvent();
  void MTtuning(string ofMTPname, string matrix_name);
  void EventSelection(bool rarm);
  void MakeTree();
  void Fill(bool rarm);
  void Close_tree();
  void Draw();
 public:



  // double xp[nmax], yp[nmax];
  double z_recon[nmax];
  div_t d;

  double ztR[MAX];
  double ztR_opt[MAX];
  double ztR_tuned[MAX];
  TFile* fnew;
  TTree* tnew;
  TCanvas* c1;

  //=======  Mzt ========//
  
  //  double Pzt_opt[nmax];
  // double Pzt_L[nParamTz];


  double refz;
  double residual;

    //===== MakeHist ========//
  TH1F* h1;
  TH1F* h2[nfoil];
  TH1F* hz_tuned; 
  TGraphErrors* gchi=new TGraphErrors();

  TH1F* hz[100];
  //==== New Tree =====//
  //   double Rvz[MAX],Lvz[MAX];
   double XFP_opt[MAX],YFP_opt[MAX],XpFP_opt[MAX],YpFP_opt[MAX];
   char tempc[500];

   //=====Event Selction =======//
     bool rtrig = false;
  bool ltrig = false;

   //===== Tuning ========//

  int evshift = 3000;

 

  //===== Fill_Tuned  =====//
  double ztr_opt;  

};


//========== zcalib class =====================//

zcalib::zcalib(){

 Setting *set = new Setting();
 set->Initialize();
 cout<<"nite: "<<nite<<endl;
};

zcalib::~zcalib(){};


//========== Set Root =============================//
void zcalib::SetRoot(string ifname){
  SetTree(ifname);
  SetBranch();
  
}

//========== Chain Tree =============================//
void zcalib::ChainRoot(string ifname){
  Chain_Tree(ifname);
  SetBranch();
  
}


//======== New Branch ====================//

void zcalib::NewBranch(string ofname, bool rarm){

  // cout<<"Create New Branch "<<endl;


  //  fnew = new TFile("zt_rhrs.root","recreate");
  fnew = new TFile(ofname.c_str(),"recreate");

  if(rarm==true) tnew = new TTree("T","For z calibration (RHRS)");
            else tnew = new TTree("T","For z calibration (LHRS)");

  tnew =t1->CloneTree(0);

  if(rarm==true)  tnew->Branch("R.tr.vz_tuned",ztR_opt, "R.tr.vz_tuned[10]/D" );
  else  tnew->Branch("L.tr.vz_tuned",ztR_opt, "L.tr.vz_tuned[10]/D" );
  /*
  if(rarm==true){
    tnew->Branch("R.tr.vz_opt",ztR, "R.tr.vz_opt[100]/D" );
    tnew->Branch("R.tr.vz_tuned",ztR_opt, "R.tr.vz_tuned[100]/D" );
    tnew->Branch("R.tr.x_opt",XFP_opt, "R.tr.x_opt[100]/D" );
    tnew->Branch("R.tr.y_opt",YFP_opt, "R.tr.y_opt[100]/D" );
    tnew->Branch("R.tr.th_opt",XpFP_opt, "R.tr.th_opt[100]/D" );
    tnew->Branch("R.tr.ph_opt",YpFP_opt, "R.tr.ph_opt[100]/D" );
  }else{
    tnew->Branch("L.tr.x_opt",XFP_opt, "L.tr.x_opt[100]/D" );
    tnew->Branch("L.tr.y_opt",YFP_opt, "L.tr.y_opt[100]/D" );
    tnew->Branch("L.tr.th_opt",XpFP_opt, "L.tr.th_opt[100]/D" );
    tnew->Branch("L.tr.ph_opt",YpFP_opt, "L.tr.ph_opt[100]/D" );
    tnew->Branch("L.tr.vz_opt",ztR, "L.tr.vz_opt[100]/D" );
    tnew->Branch("L.tr.vz_tuned",ztR_opt, "L.tr.vz_tuned[100]/D" );}
  */
  
   
}

//========== Make Hist ================//
void zcalib::MakeHist(){

  h1 = new TH1F("h1","",400,-0.5,0.5);
  h1->GetXaxis()->SetTitle("z-RHRS (m)");
  h1->GetXaxis()->SetRangeUser(-0.2,0.2);
  hz_tuned =new TH1F("hz_tuned","After Matrix tuing z hist",400,-0.5,0.5);
  //  set->SetTH1(hz_tuned,"Afeter Matrix tuing z hist","[m]","Counts");

  for(int i=0;i<nite;i++){
    hz[i]=new TH1F(Form("hz_%d",i),"",400,-0.5,0.5);
    hz[i]->SetLineColor(i+1);}
  
   h2[nfoil];
  for(int i=0 ; i<nfoil ; i++){
    sprintf(tempc,"h2_%d",i);
    h2[i] = new TH1F(tempc,tempc,
		     h1->GetXaxis()->GetNbins(),
		     h1->GetXaxis()->GetXmin(),
		     h1->GetXaxis()->GetXmax());
    h2[i]->GetXaxis()->SetTitle("z-RHRS (m)");
    h2[i]->GetXaxis()->SetRangeUser(-0.2,0.2);
  }

  gchi->SetFillColor(2);
  gchi->SetMarkerStyle(20);
  gchi->SetMarkerColor(2);
  gchi->SetFillStyle(3005);
  gchi->SetMarkerSize(1.0);

  
}




//==================== Mzt =============================//
void zcalib::Mzt(string matrix_name, bool rarm){

  if(rarm==true){
  cout<<"====== R-HRS Matrix Tuning =========="<<endl;
  // cout<<"Mzt analysis"<<endl;
  char name_Mzt[500];
  //sprintf(name_Mzt,"matrices/zt_RHRS_2.dat"); // original
  //  sprintf(name_Mzt,"./matrices/zt_RHRS_opt.dat"); // optimized parameter
  sprintf(name_Mzt,matrix_name.c_str()); // optimized
  ifstream Mzt(name_Mzt);
  cout<<"matrix file: "<<name_Mzt<<endl;
  //double Plenopt[nParamTz];
  for (int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt >> par >> p >> p >> p >> p; 
    Pzt[i]=par;
    Pzt_opt[i]=par;
  }
  Mzt.close();

  }
}

//=========== Mzt_L =======================//
void zcalib::Mzt_L(string matrix_name,bool rarm){
 
  if(rarm==false){
  cout<<"====== L-HRS Matrix Tuning =========="<<endl;
  char name_Mzt_L[500];
  //  sprintf(name_Mzt_L,"./matrices/zt_LHRS_opt.dat"); // optimized
  sprintf(name_Mzt_L,matrix_name.c_str()); // optimized
  ifstream Mzt_L(name_Mzt_L);
    cout<<"matrix file: "<<name_Mzt_L<<endl;
  for (int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt_L >> par >> p >> p >> p >> p;
    // p: matrix elements 4 elements 3order matrix
    // par: value of matrix elements  
    // Pzt_L[i]=par;
    // Pzt_opt[i]=Pzt_L[i];
    Pzt[i]=par;
    Pzt_opt[i]=par;
  
  }
  Mzt_L.close();  
  }

}

//=====================================================//

void zcalib::ZCorr(){

      XFP  = (XFP-XFPm)/XFPr;
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
      //----- Offset tuning ----------//
      ztR[0] = calcf2t_zt(Pzt, XFP, XpFP, YFP, YpFP);
      ztR[0] = ztR[0] * Ztr + Ztm;
      
      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;

     
}

//================ Event Selection Branch ============//
void zcalib::EventSelection(bool rarm){

  cout<<"==========================="<<endl;
  cout<<"===== Event Selection ====="<<endl;
  cout<<"==========================="<<endl;
  
  ntune_event=0;
  
  if(evshift>ENum)evshift=ENum;
   cout<<"Get Entries: "<<ENum<<endl;
   if(ENum<10000)d=div(ENum,1000);
   else   d=div(ENum,10000);

    for(int i=0 ; i<nmax; i++){
    x[i]    = -2222.0;
    y[i]    = -2222.0;
    xp[i]   = -2222.0;
    yp[i]   = -2222.0;
    foil_flag[i] = -1;
    }
  
    for (int i=0 ; i< ENum; i++){
    for(int j=0 ; j<MAX ; j++){
      r_x_fp[j]  = -2222.0;
      r_th_fp[j] = -2222.0;
      r_y_fp[j]  = -2222.0;
      r_ph_fp[j] = -2222.0;
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0;
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;
      ztR[j]=-2222.0;
      ztR_opt[j]=-2222.0;
      Rvz[j]=-2222.0;
      Lvz[j]=-2222.0;
      if(rarm)rvz[j]=-2222.0;
      else lvz[j]=-2222.0;
    }

    trig1 = 0.0;
    trig4 = 0.0;
    trig5 = 0.0;
    rtrig = false;
    ltrig = false;
    

    
    if(i+evshift<ENum) t1->GetEntry(i+evshift);
    else t1->GetEntry(i+evshift-ENum);
    
    if(trig4>1.0) rtrig = true;
    else rtrig = false;
    if(trig1>1.0) ltrig = true;
    else ltrig = false;

    
    if(rarm==true){
    XFP   = r_x_fp[0];
    XpFP  = r_th_fp[0];
    YFP   = r_y_fp[0];
    YpFP  = r_ph_fp[0];
    ztR[0]=rvz[0];
    }
    else{
    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    ztR[0]= lvz[0];
    }

    

    if(
       ( (rtrig==true && rarm==true) //R-HRS tuning
       ||(ltrig==true && rarm==false)) //L-HRS tuning
       && fabs(XFP)  <2.0
       && fabs(XpFP) <0.1
       && fabs(YFP)  <0.5
       && fabs(YpFP) <0.1
       ){


         h1->Fill(ztR[0]); 
      
      XFP  = (XFP-XFPm)/XFPr;
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
      //----- Offset tuning ----------//
      ztR[0] = calcf2t_zt(Pzt, XFP, XpFP, YFP, YpFP);
      ztR[0] = ztR[0] * Ztr + Ztm;
      
      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;

  
      for(int j=0 ; j<nfoil ; j++){
	if(fcent[j]-selection_width<ztR[0] 
	   && ztR[0]< fcent[j]+selection_width){
	   h2[j]->Fill(ztR[0]);
	  if(j+2!=10) h2[j]->SetLineColor(j+2);
	  else h2[j]->SetLineColor(2);
	  h2[j]->SetLineStyle(9);
	  //	  if(ntune_event<nmax){
	  if(ntune_event<nmax){
	    foil_flag[ntune_event] = j;
	    x[ntune_event]  = XFP;
	    y[ntune_event]  = YFP;
	    xp[ntune_event] = XpFP;
	    yp[ntune_event] = YpFP;
	    z_recon[ntune_event] = ztR[0];
	    ntune_event++;
	  }
	}
      }
    }


    //    if(i % 50000 == 0)cout<<i<<" / "<<ent<<endl;
    if(i % (d.quot * 1000) == 0)cout<<i<<" / "<<ENum<<endl;

    
    }
     

}


//=========== Fill ================//

void zcalib::Fill(bool rarm){

  
  cout<<"==========================="<<endl;
  cout<<"====== Fill Events ========"<<endl;
  cout<<"==========================="<<endl;  
  cout<<"rarm: "<<rarm<<endl;  
  
  ENum=t1->GetEntries();
  cout<<"Events : "<<ENum<<endl;


  
  for (int i=0 ; i< ENum; i++){
    
     for(int j=0 ; j<MAX ; j++){
      r_x_fp[j]  = -2222.0;
      r_th_fp[j] = -2222.0;
      r_y_fp[j]  = -2222.0;
      r_ph_fp[j] = -2222.0;
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0;
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;
      ztR[j]=-2222.0;
      ztR_opt[j]=-2222.0;
      //      Rvz[j]=-2222.0;
      //      Lvz[j]=-2222.0;
      
    }

     t1->GetEntry(i);

     if(rarm==true){
    XFP   = r_x_fp[0];
    XpFP  = r_th_fp[0];
    YFP   = r_y_fp[0];
    YpFP  = r_ph_fp[0];
    //   ztR[0]= Rvz[0];
    }else{
    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    //  ztR[0]= Lvz[0];
    }


      XFP  = (XFP-XFPm)/XFPr;
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
       
      //----- Matrix tuning ----------//

      ztR_opt[0] = calcf2t_zt(Pzt_opt, XFP, XpFP, YFP, YpFP);
      ztR_opt[0] = ztR_opt[0] * Ztr + Ztm;
      
      hz_tuned->Fill(ztR_opt[0]);

      if(rarm) Rvz[0] = ztR_opt[0];
      else     Lvz[0] = ztR_opt[0];

     
       
      for(int k=0;k<nite;k++){
	for(int l=0;l<nParamTz;l++)pzt_tuned[l]=Pzt_tuned[l][k];
	ztr_opt=calcf2t_zt(pzt_tuned, XFP, XpFP, YFP, YpFP);	
        ztr_opt = ztr_opt * Ztr + Ztm; 
	hz[k]->Fill(ztr_opt);}

      tnew->Fill();      

    if(i % 100000 == 0)cout<<i<<" / "<<ENum<<endl;
   }
   
 
};

//=========== Make Tree =================//
void zcalib::MakeTree(){
  cout<<"create new tree"<<endl;



  tnew->Write();
  h1->Write();
  hz_tuned->Write();

  for(int i=0;i<nfoil;i++){h2[i]->Write();}
  for(int j=0;j<nite;j++){hz[j]->Write();}
  gchi->SetName("gchi");
  gchi->Write();
}




//========= Matrix Paramter Tuning ==========//

void zcalib::MTtuning(string ofMTPname,string matrix_name){


 const char* temp=ofMTPname.c_str();
 if (nite>0) {
   cout<<"=============================="<<endl;
   cout<<"====== Tuning started: ======" << endl;
   cout<<"============================="<<endl;

   cout<<"  Iteration: "<<nite<<endl;
   cout<<"  Matrix order: "<<nnz<<endl;
   cout<<"  Reading MTtuninf file: "<<temp<<endl;
 }



  for(int i=0 ; i<nite ; i++){

    
    cout<<Form("tuning  %d/%d ",i+1,nite)<<endl;
    sprintf(tempc,"%s_%d.dat",temp,i);
    cout<<"OutPut MTtuning File : "<<tempc<<endl;

        // ---- Parameter tuning ----- //
    min_chi[i]  =  tune(Pzt_opt,i);
    
    for(int j=0 ; j<nParamTz; j++){Pzt_tuned[j][i]=Pzt_opt[j];}

    
    cout<<Form("chi_min %d :",i)<<min_chi[i]<<endl;
    gchi->SetPoint(i,i,min_chi[i]);

    
    // ------ OutPut Parameters ---- //
    ofstream * ofs = new ofstream(tempc);
    //   cout<<"tempc: "<<tempc<<endl;
    //    cout<<"ofs: "<<*ofs<<endl;
    int nppp = 0;
   

    for(int i=0 ; i<nnz+1 ; i++){
      for(int d=0 ; d<nnz+1 ; d++){
	for(int c=0 ; c<nnz+1 ; c++){
	  for(int b=0 ; b<nnz+1 ; b++){
	    for(int a=0 ; a<nnz+1 ; a++){  
	      if(a+b+c+d==i){
	      	       *ofs << Pzt_opt[nppp]
		//	ofs << Pzt_opt[nppp] 
	       	            << " " << a 
		            << " " << b
		            << " " << c
		            << " " << d << endl;
		nppp++;
		// cout << Pzt_opt[nppp]<<endl; 
		//  	      cout << Pzt_opt[nppp] 
		//		   << " " << a 
		//		   << " " << b
		//		   << " " << c
		//		   << " " << d << endl;
	      }
	    }
	  }
	}
      }
    }
    
    ofs->close();
    ofs->clear();
    //    ofs.close();
    //    ofs.clear();

    


  }



  if(nite>0){
  char temp_para[500];
  char temp_para_c[500];
  for(int i=0;i<nParamTz;i++){
    sprintf(temp_para," Pzt_%d :",i);
    sprintf(temp_para_c," Pzt_opt_%d :",i);
    cout<<temp_para<<Pzt[i]<<temp_para_c<<Pzt_opt[i]<<endl;}
  }
  
}



//========== Draw Hist ===============//
void zcalib::Draw(){

  c1=new TCanvas("c1","c1");
  c1->cd();
  h1->Draw();
  for(int i=0 ; i<nfoil ; i++){
    h2[i]->Draw("same");
  }

}


//============ Close_tree  ================//

void zcalib::Close_tree(){

  cout<<"Close Tree"<<endl;
  //f1->Close();
  fnew->Close();
}



//////////////////////////////////////////////////
double calcf2t_zt(double* P, double xf, double xpf, 
                 double yf, double ypf){
  // -----3rd order -----
  //  const int nMatT=3; 
  //  const int nXf=3;
  // const int nXpf=3;
  // const int nYf=3;
  // const int nYpf=3;

  const int nMatT=nnz; 
  const int nXf=nnz;
  const int nXpf=nnz;
  const int nYf=nnz;
  const int nYpf=nnz;

  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0;
  
  for (int n=0;n<nMatT+1;n++){
  	for (d=0;d<n+1;d++){
	  for (c=0;c<n+1;c++){
	    for (b=0;b<n+1;b++){
	      for (a=0;a<n+1;a++){
		
		if (a+b+c+d==n){
		  if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf){
		    x = pow(xf,double(a))*pow(xpf,double(b))*
		      pow(yf,double(c))*pow(ypf,double(d));
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
  
  return Y;
}


// #############################################################
double tune(double* pa, int j){
  double chi2 = 0.0;
  double arglist[10];
  int ierflg = 0;
  int allparam = nParamTz;

  
  minuit=new TMinuit(allparam);
  minuit->SetFCN(fcn);
  
  // ~~~ Chi-square ~~~~
  
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  minuit -> SetPrintLevel(-1);
  //    minuit -> SetPrintLevel(1);
  double start[allparam];
  double step[allparam];
  double LLim[allparam];
  double ULim[allparam];
  char pname[500];

  
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i];
    step[i] = 1.0e-3;
    LLim[i] = pa[i] - pa[i]*0.8;
    ULim[i] = pa[i] + pa[i]*0.8;
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }

  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0;
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  // ~~~~ Migrad + Simplex  ~~~~
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg);
  
  double amin,edm,errdef;
  amin=0.0;
  int nvpar,nparx,icstat;
  double e;
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //  cout<<"amin: "<<amin<<" edm: "<<edm<<" errdef: "<<errdef<<" nvpar: "<<nvpar<<endl;
  minuit -> mnprin(0,amin);
  if(amin>0) {chi2=amin;}

  for(int i=0 ; i<allparam ; i++){
    minuit -> GetParameter(i,Pzt_opt[i],e);



    
  }

  
  
  cout<<"tuning is done "<<endl;
  return chi2;
}

// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/){

 double chi2 = 0.0;
 const double sigma = 0.0045;
  double ztR;
  double refz = 0.0;
  double residual = 0.0;
  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;
    refz = 0.0;
    ztR = 0.0;
    XFP   = x[i];
    XpFP  = xp[i];
    YFP   = y[i];
    YpFP  = yp[i];
    refz  = fcent_real[foil_flag[i]];
    
    XFP   =(XFP -XFPm)/XFPr;
    XpFP  =(XpFP-XpFPm)/XpFPr;
    YFP   =(YFP -YFPm)/YFPr;
    YpFP  =(YpFP-YpFPm)/YpFPr;
    ztR = calcf2t_zt(param, XFP, XpFP, YFP, YpFP);
    
    ztR = ztR * Ztr + Ztm;
    residual = ztR-refz;

    //cout << residual*100 << " cm" << endl;
    if(residual<100 && -100<residual)chi2 = chi2 + pow(residual,2.0);
  
  }

  //  chi2 = sqrt(chi2)/(double)ntune_event/sigma;
  chi2 = sqrt(chi2/(double)ntune_event)/sigma;
  fval = chi2;

}


#endif
