 #ifndef Param_h
#define Param_h 1


// =================================================== //
// ==== Offset and scaling factors for matrices ====== //
// =================================================== //
const double  XFPm=-0.7,  XpFPm=-0.15; 
const double  YFPm=-0.05, YpFPm=-0.18;
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74; 
const double  XFPr=1.3,   XpFPr=0.27; 
const double  YFPr=0.1,   YpFPr=0.10; 
const double  Xptr=0.15,  Yptr=0.08, Momr=0.18; 
const double  PLm = 25.4, PLr=0.7; 
const double  Ztm = -0.15,Ztr=0.35; 
//const double  Ztm = 0.0,Ztr=1.0; //NO scale  



const int nras=2; // 2nd order matrix

//const int nn = 4; // 4th order matrix using xf, xpf, y, ypf, and zt
//const int nParamT = 126;  // Number of parameters
const int nn = 5; // angle 5th order matrix using xf, xpf, y, ypf, and zt
const int nParamT = 252;  // Number of parameters

const int nnz=5; // 3rd order matrix
const int nParamTz =126;   //Number of z parameters
const int nParamTx=4;
const int nParamTy=4;
const int nParamT_ras=4;  //Nunmber of Raster parameters


  double Pzt_opt[nParamTz];
  double Pzt[nParamTz];
  double Pzt_L[nParamTz];
  double Pzt_tuned[nParamTz][10];
  double pzt_tuned[nParamTz];
  double Pxpt[nParamT];
  double Pypt[nParamT];
  double Pras[nParamT_ras];
  double Pras_L[nParamT_ras];
  double Opt_Par_x[nParamTx],Opt_Par_y[nParamTy];

const int nmax = 500000; // Number of events used for tuning
//const int nmax=100000;// Number of events used for tuning;
const int MAX=10;
//const int nmax = 50000; // Number of events used for tuning
const int nfoil = 10;
double fcent[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
		       0.00, 0.025, 0.05, 0.10, 0.125}; 
double fcent_real[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
			    0.000, 0.025, 0.050, 0.100, 0.125};
//double selection_width = 0.0125; 
//double selection_width = 0.008; // event selection width for z
double selection_width = 0.008; // event selection width for z

const double step = 0.492 * 2.54;
const int nrow = 11; // the number of row in SS pattern
const int ncol = 7;  // the number of column in SS pattern
const int nsshole = nrow*ncol; // the number of holes to consider 
double refx[nsshole];
double refy[nsshole];
double selec_widthx = 0.60; // selection width in x (dispersive plane)
double selec_widthy = 0.45; // selection width in y 
int ntune_event;





double x[nmax], y[nmax];
double xp[nmax], yp[nmax];
double Zrecon[nmax]; // reconstructed z position
int foil_flag[nmax];  
int holegroup[nmax];
double Zt_wRC[nmax];

double l[nfoil];
double projectf[nfoil];
double OptPar1[nParamT];
double OptPar2[nParamT];


const double hrs_ang = 13.2 * 3.14159 / 180.;
const double tdc_time =56e-3;//[ns]

//======= Coincidence Offset Paramters ===========//

const double  RS2_off_H1[16]={-16828.7,-16863,-16894,-16893.3,-16870.9,-16867.2,-16900.3,-16876.8,17554.1,-16861.6,-16895,-16890.7,-16854.6,-16852.9,-16850.5,-16861.9};
const double  LS2_off_H1[16]={-25335,-25385.6,-25367,-25392.1,-25391.7,-25386.4,-25422.1,-25428.9,-25414.9,-25424.7,-25436.9, -25381.2,-25390,-25413.4,-25428.7,-26640.8};

const  double pathl_off=-470.5;
const  double s2_offset=-499.75;
const  double coin_offset=456.65;

#endif
