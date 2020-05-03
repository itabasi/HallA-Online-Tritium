#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
using namespace std;
#include <TApplication.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>
#include <TCut.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>
#include <TText.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TColor.h>
#include <TPaveText.h>
#include <TRandom.h>
#include <TMinuit.h>
#include "recon.h"
#include "Setting.h"
#include "define.h"

extern double calcf2t_zt(double* P, double xf, double xpf, double yf, double ypf);
extern double Calc_ras(double a,double b,double c){return  a *b + c;};   
extern double calcf2t_ang(double* P,double xf, double xpf, double yf, double fpf,double z);
extern double calcf2t_mom(double* P, double xf, double xpf, double yf, double ypf, double zt);



///////////////////////////////////////////////

void recon::SetRoot(string ifname){

  T->Add(ifname.c_str());
  ENum=T->GetEntries();
  ENum=100000.;
  cout<<"Get Entries: "<<ENum<<endl;
}
///////////////////////////////////////////////

void recon::SetBranch(){

 T->SetBranchStatus("*",1);  
 T->SetBranchAddress("h_xfp",&R_tr_x);
 T->SetBranchAddress("h_yfp",&R_tr_y);
 T->SetBranchAddress("h_xpfp",&R_tr_th);
 T->SetBranchAddress("h_ypfp",&R_tr_ph);  
 T->SetBranchAddress("h_xptar",&R_tr_tg_th);
 T->SetBranchAddress("h_yptar",&R_tr_tg_ph);
 
 T->SetBranchAddress("e_xfp",&L_tr_x);
 T->SetBranchAddress("e_yfp",&L_tr_y);
 T->SetBranchAddress("e_xpfp",&L_tr_th);
 T->SetBranchAddress("e_ypfp",&L_tr_ph);   
 T->SetBranchAddress("e_xptar",&L_tr_tg_th);
 T->SetBranchAddress("e_yptar",&L_tr_tg_ph);
 T->SetBranchAddress("zposi",&zposi);
 T->SetBranchAddress("Ein",&Ee);
 T->SetBranchAddress("rasx",&R_Ras_x);
 T->SetBranchAddress("rasx",&L_Ras_x);
 T->SetBranchAddress("mmnuc",&mm); 
 
}


////////////////////////////////////////////////

void recon::SetMatrix(string mtparam){

  cout<<endl;
  cout<<"==============================="<<endl;
  cout<<"=== Input Matrix Parameters ==="<<endl;
  cout<<"==============================="<<endl;

  string buf;
  int s=0;
  ifstream ifp(Form("%s",mtparam.c_str()),ios::in);
  if (ifp.fail()){ cerr << "failed open files" <<mtparam.c_str()<<endl; exit(1);}
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp
.eof() ) break;
    stringstream sbuf(buf);
    sbuf >>param[s];
    cout<<param[s]<<endl;
    s++;
  }


    MTParam_R();cout<<" Input RHRS Matrix parameter "<<endl;
    MTParam_L();cout<<" Input LHRS Matrix parameter "<<endl;
  MTP_mom();

}

/////////////////////////////////////////////////

void recon::MTParam_R(){

  //=================//
  //==== RHRS =======//
  //=================//

  
  //====== RHRS z parameters ======//

    char name_Mzt[500];
    sprintf(name_Mzt, param[0].c_str()); // optimized
    ifstream Mzt(name_Mzt);

   for(int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt >> par >> p >> p >> p >> p;
    Pzt[i]=par;
    //    cout<<"R Mzt : "<<Pzt[i]<<endl;
   }
  Mzt.close();

  
  //====== RHRS raster paramters =======//
    char name_Mras[500];
    sprintf(name_Mras, param[1].c_str()); // optimized
    //    cout<<"RHRS Raster parameters file: "<<name_Mras<<endl;
  ifstream Mras(name_Mras);

  for (int i=0;i<nParamT_ras;i++){

    Mras >> Pras[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) Pras[i] = 0.0;
  }


  Mras.close();    

  
  //===== RHRS theta parameters ======// 
    char name_Mxpt[500];
    sprintf(name_Mxpt, param[2].c_str()); // optimized
  ifstream Mxpt(name_Mxpt);

    for(int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt >> par >> p >> p >> p >> p >> p;
    Pxpt[i]  = par;
  }
  Mxpt.close();  
 //===== RHRS phi parameters ======//
  char name_Mypt[500];
    sprintf(name_Mypt, param[3].c_str()); // optimized  
  ifstream Mypt(name_Mypt);

  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mypt >> par >> p >> p >> p >> p >> p;
    Pypt[i]  = par;
  }
  Mypt.close();    



  
};
//////////////////////////////////////////////////////////////

void recon::MTParam_L(){

  //=================//
  //===== LHRS ======//
  //=================//

  
  //====== LHRS z parameters ======//  
  char name_Mzt_L[500];
  sprintf(name_Mzt_L,param[4].c_str()); // optimized
  ifstream Mzt_L(name_Mzt_L); 
  for (int i=0;i<nParamTz;i++){
    double par=0.;
    int p=0;
    Mzt_L >> par >> p >> p >> p >> p;
    Pzt_L[i]=par;
  }
 Mzt_L.close();

  //====== LHRS raster paramters =======//
    char name_Mras_L[500];
    sprintf(name_Mras_L, param[5].c_str()); // optimized
    //    cout<<"LHRS Raster parameters file: "<<name_Mras_L<<endl;
  ifstream Mras_L(name_Mras_L);

  for (int i=0;i<nParamT_ras;i++){

    Mras_L >> Pras_L[i];
    // ---- raster x is tuned by this macro ---- //
    if(i==1 || i==3) Pras_L[i] = 0.0;
  }
  
  Mras_L.close();    

 
 //===== LHRS theta parameters ======// 
  char name_Mxpt_L[500];
    sprintf(name_Mxpt_L, param[6].c_str()); // optimized
  ifstream Mxpt_L(name_Mxpt_L);
  //  cout<<"LHRS theta parameters file: "<<name_Mxpt_L<<endl;  
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mxpt_L >> par >> p >> p >> p >> p >> p;
    //   cout<<"LHRS theta : "<<par<<endl;
    Pxpt_L[i]  = par;
  }
  Mxpt_L.close();

  
 //===== LHRS phi parameters ===x==//
  char name_Mypt_L[500];
    sprintf(name_Mypt_L, param[7].c_str()); // optimized
  ifstream Mypt_L(name_Mypt_L);

  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mypt_L >> par >> p >> p >> p >> p >> p;
    Pypt_L[i]  = par;    
  }
  Mypt_L.close();    


}

//========================================================//


void recon::MTP_mom(){

  //====== RHRS Momentum parameters ========//
    char name_Mpt[500];
    sprintf(name_Mpt, param[8].c_str()); // optimized
    ifstream Mpt(name_Mpt);
    if (Mpt.fail()){ cerr << "failed open files" <<name_Mpt<<endl; exit(1);}
   for(int i=0;i<nParamTp;i++){
    double par=0.;
    int p=0;
    Mpt >> par >> p >> p >> p >> p >> p;
    Prp[i]=par;
    Opt_par_R[i]=par;
    Opt_par[i]=par;
   }

 
   Mpt.close();

   

  //====== LHRS Momentum parameters ========//
    char name_Mpt_L[500];
    sprintf(name_Mpt_L, param[9].c_str()); // optimized
    ifstream Mpt_L(name_Mpt_L);
    if (Mpt_L.fail()){ cerr << "failed open files" <<name_Mpt_L<<endl; exit(1);}
   for(int i=0;i<nParamTp;i++){
    double par=0.;
    int p=0;
    Mpt_L >> par >> p >> p >> p >> p >> p;
    Plp[i]=par;
    Opt_par_L[i]=par;
    Opt_par[i+nParamTp]=par;  // Both momentum paramters
    //cout<<Form("Opt_par[%d]",i)<<Opt_par[i+nParamTp]<<endl;
   }
  Mpt_L.close();



}

/////////////////////////////////////////////////

void recon::Fill(){

  double cm_to_m=0.01;
  div_t d;
  if(ENum<10000)d=div(ENum,1000);
  else   d=div(ENum,10000);
  
  for(int n =0;n<ENum;n++){
    
    R_tr_x=-100.0;
    R_tr_y=-100.;
    R_tr_th=-100;
    R_tr_ph=-100;
    R_tr_tg_th=-100;
    R_tr_tg_ph=-100;
    L_tr_x=-100.0;
    L_tr_y=-100.;
    L_tr_th=-100;
    L_tr_ph=-100;
    L_tr_tg_th=-100;
    L_tr_tg_ph=-100;
    R_tr_vz=-100.0;
    L_tr_vz=-100.;    
    mm= -100.;
    mm_c=-100;
    zposi=-100;
    R_p=0.0;
    Rp_x=0.0;
    Rp_y=0.0;
    Rp_z=0.0;
    L_p=0.0;
    Lp_x=0.0;
    Lp_y=0.0;
    Lp_z=0.0;
    Ee=-100;
    Ee_=-100;
    Ek=-100;

    
    T->GetEntry(n);
    //====== convet unit cm to m ====//
    R_tr_x *= cm_to_m;
    R_tr_y *= cm_to_m;
    L_tr_x *= cm_to_m;
    L_tr_y *= cm_to_m;

    R_tr_vz = zposi *cm_to_m;
    L_tr_vz = zposi *cm_to_m;
    
    
    Calib();

    
    Ee  *= 1./1000.; // MeV to GeV
    Ee_ = sqrt(L_p*L_p + Me*Me);
    Ek = sqrt(R_p*R_p + MK*MK);
    
    TVector3 B_v(0.0,0.0,sqrt(Ee*Ee-Me*Me));
    TVector3 R_v(Rp_x, Rp_y, Rp_z);
    TVector3 L_v(Lp_x, Lp_y, Lp_z);
    R_v.RotateX(   13.2/180.*3.14);
    L_v.RotateX( - 13.2/180.*3.14);



    if(mm>2.0){
      mm_c = sqrt( (Ee + MTr - Ee_ - Ek)*(Ee + MTr - Ee_ - Ek)
		   - (B_v - L_v - R_v)*(B_v - L_v - R_v) ); // Tritium
      mm_c = (mm_c - MnnL)*1000.;
      mm = (mm - MnnL)*1000.;
    }else{
      mm_c = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
		   - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
      mm_c = (mm_c - ML)*1000.;
      mm   = (mm - ML)*1000.;
    }
    
    


    
    tnew->Fill();
    if(n % (d.quot * 10000) == 0)cout<<n<<" / "<<ENum<<endl;
  }

}



/////////////////////////////////////////////////


void recon::Calib(){

  // ==== Initialization ======//
  R_p = 0.0;
  L_p = 0.0;


  //======= Nomalization ==================//
  R_tr_x    = (R_tr_x-XFPm)/XFPr;
  R_tr_th   = (R_tr_th-XpFPm)/XpFPr;
  R_tr_y    = (R_tr_y-YFPm)/YFPr;
  R_tr_ph   = (R_tr_ph-YpFPm)/YpFPr;
  R_tr_vz   = (R_tr_vz-Ztm)/Ztr;
  R_tr_tg_th= (R_tr_tg_th - Xptm)/Xptr;
  R_tr_tg_ph= (R_tr_tg_ph - Yptm)/Yptr;



  L_tr_x    = (L_tr_x-XFPm)/XFPr; 
  L_tr_th   = (L_tr_th-XpFPm)/XpFPr;
  L_tr_y    = (L_tr_y-YFPm)/YFPr;
  L_tr_vz   = (L_tr_vz-Ztm)/Ztr;
  L_tr_ph   = (L_tr_ph-YpFPm)/YpFPr;
  L_tr_tg_th= (L_tr_tg_th - Xptm)/Xptr;
  L_tr_tg_ph= (L_tr_tg_ph - Yptm)/Yptr;  



  //========================================//
  
  R_tr_vz   = calcf2t_zt(Pzt, R_tr_x, R_tr_th, R_tr_y, R_tr_ph); // nomalized
  L_tr_vz   = calcf2t_zt(Pzt_L, L_tr_x, L_tr_th, L_tr_y, L_tr_ph); //nomalized


    //======== Raster Correction ==========================//    


  // convert units mm to ch //
  R_Ras_x = (R_Ras_x)*(88650 -56600)/0.18 + (88650 -56600)/2.0 ;
  R_Ras_x = (R_Ras_x)*(95000 -60200)/0.18 + (95000 -60200)/2.0 ;
  
    R_tr_vz  = R_tr_vz*Ztr +Ztm; // scaled
    RasterCor= Calc_ras(R_Ras_x, Pras[2], Pras[0]);
    RasterCor = RasterCor/tan(hrs_ang);
    R_tr_vz  = R_tr_vz + RasterCor; // correction
    R_tr_vz  = (R_tr_vz-Ztm)/Ztr;    // nomalization     
    L_tr_vz  = L_tr_vz*Ztr +Ztm;     // scaled
    RasterCor_L = Calc_ras(L_Ras_x, Pras_L[2], Pras_L[0]);
    RasterCor_L = RasterCor_L/tan(hrs_ang);
    L_tr_vz  = L_tr_vz + RasterCor_L;
    L_tr_vz  =  (L_tr_vz  -  Ztm)/Ztr;    // nomalization

    //====================================================//

    


    R_tr_tg_th  = calcf2t_ang(Pxpt,   R_tr_x, R_tr_th, R_tr_y, R_tr_ph,R_tr_vz); // nomalized
    R_tr_tg_ph  = calcf2t_ang(Pypt,   R_tr_x, R_tr_th, R_tr_y, R_tr_ph,R_tr_vz); // nomalized
    L_tr_tg_th  = calcf2t_ang(Pxpt_L, L_tr_x, L_tr_th, L_tr_y, L_tr_ph, L_tr_vz); // nomalized
    L_tr_tg_ph  = calcf2t_ang(Pypt_L, L_tr_x, L_tr_th, L_tr_y, L_tr_ph, L_tr_vz); // nomalized   

    R_p = calcf2t_mom(Opt_par_R, R_tr_x, R_tr_th, R_tr_y, R_tr_ph,R_tr_vz);
    L_p = calcf2t_mom(Opt_par_L, L_tr_x, L_tr_th, L_tr_y, L_tr_ph,L_tr_vz);


    
    //========== Scaled at FP ==================//
    R_tr_x  = R_tr_x  * XFPr + XFPm;
    R_tr_th = R_tr_th * XpFPr + XpFPm;
    R_tr_y  = R_tr_y  * YFPr + YFPm;
    R_tr_ph = R_tr_ph * YpFPr + YpFPm;

    L_tr_x  = L_tr_x  * XFPr + XFPm;
    L_tr_th = L_tr_th * XpFPr + XpFPm;
    L_tr_y  = L_tr_y  * YFPr + YFPm;
    L_tr_ph = L_tr_ph * YpFPr + YpFPm;    

    //=========== Scaled at Taget =============//

    R_tr_vz     = R_tr_vz * Ztr + Ztm; // scaled
    R_tr_tg_th  = R_tr_tg_th * Xptr + Xptm; // scaled
    R_tr_tg_ph  = R_tr_tg_ph * Yptr + Yptm; // scaled
    R_p             = R_p * PRr + PRm; // scaled
    L_tr_vz     = L_tr_vz * Ztr + Ztm; // scaled
    L_tr_tg_th  = L_tr_tg_th * Xptr + Xptm;  // scaled    
    L_tr_tg_ph  = L_tr_tg_ph * Yptr + Yptm;  // scaled    
    L_p             = L_p * PLr + PLm; // scaled    
    



    
    // Lp = 2.2 GeV mode //
    L_p=2.21807/2.1*L_p;


    double dpe = Eloss(0,0,"B");
    double dpk = Eloss(R_tr_tg_ph,R_tr_vz,"R");
    double dpe_= Eloss(R_tr_tg_ph,L_tr_vz,"L");


    Ee += dpe;
    R_p += dpk;
    L_p += dpe_;


    
      //==== Right Hand Coorinate ====//

    Lp_z = L_p/sqrt(1.0*1.0 + tan(L_tr_tg_th)*tan(L_tr_tg_th) + tan(L_tr_tg_ph)*tan(L_tr_tg_ph));
    Rp_z = R_p/sqrt(1.0*1.0 + tan(R_tr_tg_th)*tan(R_tr_tg_th) + tan(R_tr_tg_ph)*tan(R_tr_tg_ph));


   

    Lp_x = Lp_z*    tan( L_tr_tg_th );
    Lp_y = Lp_z*    tan( L_tr_tg_ph );
    Rp_x = Rp_z*    tan( R_tr_tg_th );
    Rp_y = Rp_z*    tan( R_tr_tg_ph ); 


    
}

/////////////////////////////////////////////////////////////////////////////

double recon::Eloss(double yp, double z,char* arm){



  
  double x;
  
  //  if(arm) x = - tan(hrs_ang-yp); //yp : phi [rad] RHRS
  //  else    x = - tan(hrs_ang+yp); //yp : phi [rad] LHRS

  //----- Original coordinate  -------//
  // Definition by K. Suzuki  (fixed Oct. 23rd, 2019)       //
  // R-HRS : right hand coordinate (Unticlockwise rotation) //
  // L-HRS : left  hand coordinate (    Clockwise rotation) //


  
  //  if(arm=="R")        x = - hrs_ang - yp; //yp : phi [rad] RHRS
  //  else if(arm=="L")   x = - hrs_ang + yp; //yp : phi [rad] LHRS
  
  if(arm=="R")        x = - hrs_ang + yp; //yp : phi [rad] RHRS
  else if(arm=="L")   x = - hrs_ang - yp; //yp : phi [rad] LHRS
  else x=0.0;


  double ph[3],pl[2];
  double dEloss;
  bool high;
  
  if(z>0.08)high=false; //low : low  energy Loss (z> 8 cm)  pol1 function
  else high=true;       //high  : high  energy Loss (z< 8 cm) sin function
  
    //==== thickness 0.400 mm ========//

  if(arm=="R"){
    ph[0] = -1.31749;
    ph[1] = -4.61513;
    ph[2] = 2.03687;
    pl[0] = 3.158e-2;
    pl[1] = 4.05819e-1;
    
  }else if(arm=="L"){
    ph[0] = -1.3576;
    ph[1] = -4.5957;
    ph[2] = 2.0909;
    pl[0] = 6.2341e-3;
    pl[1] = 4.0336e-1;
    
  }

  double dEloss_h = ph[0]*sin(ph[1]*x)+ph[2];
  double dEloss_l = pl[0]*x +pl[1];

  
  if(high)dEloss = dEloss_h;
  else dEloss = dEloss_l;
  
  
  
  //==== thickness 0.4 mm in beam energy loss ======//
  if(arm=="B")dEloss=0.1843; //[MeV/c]



  //======== Al Flame Energy Loss ========//
  // Upstream && Downsream //
  // thickness 0.4 mm //
  if(z<-0.1){
    if(arm=="B")dEloss  = 0.1345;
    if(arm=="R")dEloss += 0.803;
    if(arm=="L")dEloss += 0.809;
  }
  else if(z>0.1){
    if(arm=="B")dEloss  = 0.264;
    if(arm=="R")dEloss += 0.313;
    if(arm=="L")dEloss += 0.313;
    
  }
  
  dEloss=dEloss/1000.; // [GeV/c]  
  return dEloss;
  
}

////////////////////////////////////////////////////////////////////////////


void recon::NewRoot(string ofname){

  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  tnew=new TTree("T","Momentum matrix tuning");
  tnew =T->CloneTree(0);
  tnew->Branch("mm", &mm_c,"mm/D");
  tnew->Branch("Rp_rec", &R_p,"Rp_rec/D");
  tnew->Branch("Lp_rec", &L_p,"Lp_rec/D");
  tnew->Branch("Rz_rec", &R_tr_vz,"Rz_rec/f");
  tnew->Branch("Lz_rec", &L_tr_vz,"Lz_rec/f");
  tnew->Branch("Rth_rec", &R_tr_tg_th,"Rth_rec/f");
  tnew->Branch("Lth_rec", &L_tr_tg_th,"Lth_rec/f");
  tnew->Branch("Rph_rec", &R_tr_tg_ph,"Rph_rec/f");
  tnew->Branch("Lph_rec", &L_tr_tg_ph,"Lph_rec/f");
  
}

///////////////////////////////////////////////////////////

void recon::Write(){

  tnew->Write();
  ofp->Close();

}

//////////////////////////////////////////////////////////////

//===================================================================//
//============================= Main ================================//
//==================================================================//

int main(int argc, char** argv){


  gStyle->SetOptFit(111111111);
  int ch; char* mode="C";

  string ifname = "../rootfiles/momcalib/test.root";
  string ofname = "../rootfiles/momcalib/test.root";
  string matrix_name="./matrix/momcalib_matrix.list";
  string ofMTPname="./matrix/test/test.dat";
  string opt_file="./scale_offset_20190210.dat";  
  string iteration;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = false;
  bool root_flag=false;
  bool matrix_flag=false;
  bool RHRS_flag=false; 
  bool tuning_flag=false;
  string pngname;
  extern char *optarg;
  char* Target="H";
  string F1tdc="1";
  int f1tdc=1;
  bool Al1=false;
  bool Al_MODE=false;
  int nmatrix;
  bool nmatrix_flag=false;
  string root_init="../rootfiles/";
  string root_end=".root";
  string dat_init="../matrix/";
  string dat_end=".dat";
  string matrix="matrix/test/test.dat";
  
  while((ch=getopt(argc,argv,"h:s:w:t:p:f:n:r:m:o:O:i:AlRILbcop"))!=-1){
    switch(ch){
      
      
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;



      
    case 'r':
      root_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;


    case 't':
      Target  = optarg;
      if(Target=="H")cout<<"Target : Hydrogen "<<endl;
      if(Target=="T")cout<<"Target : Tritium "<<endl;      
      break;

      
    case 'm':

      matrix_name =optarg;
      break;


      
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;

    case 'p':
      matrix_name =optarg;
      break;
      

    case 'h':
      cout<<"-f : input root  filename"<<endl;
      cout<<"-m : input matrix filename"<<endl;      
      cout<<"-w : output matrix filename"<<endl;
      cout<<"-r : output root filename"<<endl;
      cout<<"-Al: w/ Al events mode    "<<endl;
      cout<<"-n : matrix order "<<endl;
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

  recon* Rec=new recon();

  Rec->SetMatrix(matrix_name);
  Rec->SetRoot(ifname);
  Rec->SetBranch(); 
  if(root_flag)Rec->NewRoot(ofname);
  if(root_flag) Rec->Fill();  
  if(root_flag) Rec->Write();  


  cout<<endl;
  cout<<"============== Input files ==============="<<endl;
  cout<<"run files   : "<<ifname<<endl;
  cout<<"matrix param : "<<matrix_name<<endl;
  cout<<"analysis mode : "<<"Coincidence mode "<<endl;
  cout<<endl;
  cout<<"============== Output files ==============="<<endl;
  cout<<"root_flag "<<root_flag<<endl;
  if(root_flag)cout<<"Root file   : "<<ofname<<endl;
  if(matrix_flag)cout<<"New matrix file : "<<ofMTPname<<endl;
  cout<<"==========================================="<<endl;
  cout<<endl;


  
   gSystem->Exit(1);
   theApp->Run();

 
  return 0;
  
}//end main




///////////////////////////////////////////////////////////

// ###################################################
double calcf2t_mom(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt)
// ###################################################
{
  // ----------------------------------------------------------------- //
  // ------ 4rd order using xf, xpf, yf, ypf, zt, xt, xpt, yt, ytp --- //
  // ----------------------------------------------------------------- //



  
  const int nMatT =nnp;  
  const int nXf   =nnp;
  const int nXpf  =nnp;
  const int nYf   =nnp;
  const int nYpf  =nnp;
  const int nZt   =nnp;

  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
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




// ###################################################
double calcf2t_ang(double* P, double xf, double xpf, 
		     double yf, double ypf, double zt)
// ####################################################
{
  // ------------------------------------------------ //
  // ----- 4rd order using xf, xpf, yf, ypf, zt ----- //
  // ------------------------------------------------ //
  
  const int nMatT=nn;  
  const int nXf=nn;
  const int nXpf=nn;
  const int nYf=nn;
  const int nYpf=nn;
  const int nZt=nn;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
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
double calcf2t_zt(double* P, double xf, double xpf, 
                 double yf, double ypf){
// ###############################################

  int nnz=3;  

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
