#ifndef Tuning_h
#define Tuning_h 1

#include <iostream>
#include <fstream>
#include "Param.h"
#include "Setting.h"

//=================================================//
//=============     Tuning     ====================//
//=================================================//

class Tuning{

 public:
  Tuning();
  ~Tuning();

  double tune(double * P, int i);
  void Minuit(int param_name, int allparam);
  
  double calcf2t_zt(double * P,
		    double xf, double xpf,
		    double yf, double ypf);

  double calcf2t_ang(double* P,
			    double xf, double xpf, 
			    double yf, double fpf,
			    double zt);  

  double calcRasterCor(double a, double b, double c);

  //--- calcf2t_zt --------//
  const int nMatTz=nnz; 
  const int nXf=nn;
  const int nXpf=nn;
  const int nYf=nn;
  const int nYpf=nn;
  double Y,x; 
  int npar;
  int a,b,c,d,e;
  
  //---- calcf2t_ang -----//
  const int nMatT=nn;
  const int nZt=nn;

  //--- calcRasterCor ----//
  
};

/////////////////////////////////////////

Tuning::Tuning(){};
Tuning::~Tuning(){};

/////////////////////////////////////////
double tune(double *P, int i, int allparam){


  double chi2 = 0.0;
  double arglist[10]; 
  int ierflg = 0;
  //  int allparam = nParamT2;
  //cout << allparam << endl;

  TMinuit* minuit = new TMinuit(allparam);
  minuit->SetFCN(fcn);
  
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1);
  double start[allparam];
  double step[allparam];
  double LLim[allparam];// Lower limit for all of the parameter
  double ULim[allparam];// Upper limit for all of the parameter
  char pname[500];

  int param_name=1;//raster calibration
  Minuit(param_name,allparam);


  
  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0;
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  
  // ~~~~ Migrad + Simplex  ~~~~ 
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg); // Chi-square minimization
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double e;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin);

  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){
    minuit -> GetParameter(i,Opt_Par[i],e);
  }
  
  return chi2;


}

//////////////////////////////////////////////////////////////////

void Tuning::Minuit(int param_name, int allparam){

  //============== param_name =============//
  // 0 : zcalib, 1: rascalib, 2: angcalib  //  
  //======================================//

  if(param_name==1){

    //=================================//
    //==== Raster Calibration =========//
    //================================//
    
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i]; // initial parameters
    //step[i] = 1.0e-5;      
    if(i==0){       // offset parameter for raster x
      step[i] = 1.0e-5;
    }
    else if(i==2){  // gradient parameters for raster x
      step[i] = 1.0e-10;
    }
    else step[i] = 0.0; // raster y is not tuned here
    
    //LLim[i] = pa[i] - pa[i]*0.8;
    //ULim[i] = pa[i] + pa[i]*0.8;
    LLim[i] = pa[i] -1.0; // temporary 
    ULim[i] = pa[i] +1.0; // temporary 
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }


  }else if(param_name==2){}

  
}


#endif
