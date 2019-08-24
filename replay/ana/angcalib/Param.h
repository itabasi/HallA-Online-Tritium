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

#define max 100


//const int nn =  2; // 4th order matrix using xf, xpf, y, ypf, and zt
//const int nnz = 2; // 3th order matrix using zt
//const int nParamT = 21;  // Number of parameters

//const int nn =  3; // 4th order matrix using xf, xpf, y, ypf, and zt
//const int nnz = 3; // 3th order matrix using zt
//const int nParamT = 56;  // Number of parameters

const int nn =  4; // 4th order matrix using xf, xpf, y, ypf, and zt
const int nnz = 4; // 3th order matrix using zt
const int nParamT = 126;  // Number of parameters

//const int nn =  5; // 4th order matrix using xf, xpf, y, ypf, and zt
//const int nnz = 4; // 3th order matrix using zt
//const int nParamT = 252;  // Number of parameters



  double Pxpt[nParamT];
  double Pypt[nParamT];
  double OptPar1[nParamT];
  double OptPar2[nParamT];


//const int nmax = 300; // Number of events used for tuning
const int nmax=10000;// Number of events used for tuning;
//const int nmax=20000;// Number of events used for tuning;
//const int nmax = 50000; // Number of events used for tuning
//const int nite =0;
const int nfoil = 10;
double fcent[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
		       0.00, 0.025, 0.05, 0.10, 0.125}; 
double fcent_real[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
			    0.000, 0.025, 0.050, 0.100, 0.125};
//double selection_width = 0.0125; 
double selection_width = 0.01; // event selection width for z


const double step = 0.492 * 2.54;
const int nrow = 11; // the number of row in SS pattern
const int ncol = 7;  // the number of column in SS pattern
const int nsshole = nrow*ncol; // the number of holes to consider 
double w[nfoil][nsshole];
double refx[nsshole];
double refy[nsshole];
double refx_real[nsshole][nfoil];
double refy_real[nsshole][nfoil];
double ssy_off[nsshole][nfoil];
double ssx_off[nsshole][nfoil];
bool TFlag[nsshole][nfoil];
//double selec_widthx = 0.60; // selection width in x (dispersive plane)
//double selec_widthy = 0.45; // selection width in y
double selec_widthx = 0.60; // selection width in x (dispersive plane) initial
double selec_widthy = 0.45; // selection width in y initial
int ntune_event;





double x[nmax], y[nmax];
double xp[nmax], yp[nmax];
double th[nmax], ph[nmax];
double z_recon[nmax]; // reconstructed z position
int foil_flag[nmax];  
int holegroup[nmax];

double l[nfoil];
double projectf[nfoil];

//const int nParamT2 = 4; 
//double parRaster[nParamT2];
//double Opt_Par[nParamT2];
//double Ras_curx[nmax];
//double Ras_cury[nmax];

const double hrs_ang = 13.2 * 3.14159 / 180.;


#endif
