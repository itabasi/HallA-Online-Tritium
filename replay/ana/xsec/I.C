#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

Double_t ef_1s( Double_t* x, Double_t* par ) {

  Double_t q = x[0] / 197.3; //  wave length of relative momentum [fm]:  momentum [MeV] / hbar    
  Double_t a = par[0];       //  scattering length [fm]
  Double_t r = par[1];       //  effective range   [fm]

  // Enhancement factor [EPJ A21 (2004) 313-321]
  //
  // | C_FSI |^2 = ( q^2 + beta^2 ) / ( q^2 + alpha^2 )
  //
  //             1                                      1
  //     alpha = - ( 1 - sqrt( 1 - 2r/a ) ),    beta  = - ( 1 + sqrt( 1 - 2r/a ) )
  //             r                                      r
  
  Double_t alpha = ( 1 - TMath::Sqrt( 1 - 2 * r / a ) ) / r;
  Double_t  beta = ( 1 + TMath::Sqrt( 1 - 2 * r / a ) ) / r;

  
  return  ( q*q + beta*beta) / ( q*q + alpha*alpha );

}


Double_t ef_1( Double_t* x, Double_t* par ) {

  Double_t q   = x[0] / 197.3; //  wave length of relative momentum [fm]:  momentum [MeV] / hbar
  Double_t a_s = par[0];       //  scattering length [fm] for single
  Double_t r_s = par[1];       //  effective range   [fm] for single
  Double_t a_t = par[2];       //  scattering length [fm] for triplet
  Double_t r_t = par[3];       //  effective range   [fm] for triplet


  Double_t alpha_s = ( 1 - TMath::Sqrt( 1 - 2 * r_s / a_s ) ) / r_s;
  Double_t  beta_s = ( 1 + TMath::Sqrt( 1 - 2 * r_s / a_s ) ) / r_s;

  Double_t alpha_t = ( 1 - TMath::Sqrt( 1 - 2 * r_t / a_t ) ) / r_t;
  Double_t  beta_t = ( 1 + TMath::Sqrt( 1 - 2 * r_t / a_t ) ) / r_t;


  
  return    0.25 * ( q*q + beta_s*beta_s ) / ( q*q + alpha_s*alpha_s )
          + 0.75 * ( q*q + beta_t*beta_t ) / ( q*q + alpha_t*alpha_t );

}






Double_t ef_2( Double_t*x, Double_t* par ) {

  // Phys. Rev. C76 (2007) 054004, JLab E91-016 collaboration
  // 
  // 
  // I = (  sin (delta0 + k r ) / sin( k r )  )^2
  // 
  //  k / tan(delta0) = - 1/a + r k^2 / 2
  //  1 / tan(delta0) = - 1/(a k) + r k / 2  = ( a r k^2  - 2 ) / ( 2 a k )
  //      tan(delta0) =  2 a k / ( a r k^2 - 2 )
  // 

  Double_t k         = x[0] / 197.3; //  wave length of relative momentum [fm]:  momentum [MeV] / hbar    
  Double_t a_singlet = par[0];       //  scattering length [fm] for single
  Double_t r_singlet = par[1];       //  effective range   [fm] for single
  Double_t a_triplet = par[2];       //  scattering length [fm] for triplet
  Double_t r_triplet = par[3];       //  effective range   [fm] for triplet


  Double_t delta0_s = TMath::ATan(  k / ( -1./a_singlet + 0.5 * r_singlet * k * k ) );
  Double_t delta0_t = TMath::ATan(  k / ( -1./a_triplet + 0.5 * r_triplet * k * k ) );

  return    0.25 * TMath::Power( TMath::Sin( delta0_s + k * r_singlet ) / TMath::Sin( k * r_singlet )  , 2 ) 
          + 0.75 * TMath::Power( TMath::Sin( delta0_t + k * r_triplet ) / TMath::Sin( k * r_triplet )  , 2 );

}



void I()
{

  gROOT -> Reset();
  gROOT -> SetStyle( "Plain" );


  // Enhancement factor I

  // Effective range approximation, a: effective range [fm],   r: scattering length [fm]
  // 
  // 
  //                        _______________ Nijmegen D [PRD15 (1977) 2547]            
  //                       |       ________ Nijmegen F [PRD20 (1979) 1633]            
  //                       |      |       _ NSC89 [PRC40 (1989) 2226]                 
  //                       |      |      |       ____________________________________  NSC97a [PRC59 (1999) 3009]
  //                       |      |      |      |       _____________________________  NSC97b
  //                       |      |      |      |      |        _____________________  NSC97c
  //                       |      |      |      |      |       |       ______________  NSC97d
  //                       |      |      |      |      |       |      |       _______  NSC97e
  //                       |      |      |      |      |       |      |      |       _ NSC97f
  //                       |      |      |      |      |       |      |      |      |            ________________  Juelich A (Lambda-N) [NPA570 (1994) 543]
  //                       |      |      |      |      |       |      |      |      |           /       __________  Juelich A~(Lambda-N) [NPA570 (1994) 543]
  //                       |      |      |      |      |       |      |      |      |          /       /      ____  Juelich B (Lambda-N) [NPA570 (1994) 543]
  //                       |      |      |      |      |       |      |      |      |         /       /      /      Juelich B~(Lambda-N) [NPA570 (1994) 543]
  //                       |      |      |      |      |       |      |      |      |        /       /      /      /      __  Verma [PRC22 (1980)229]  
  //                       |      |      |      |      |       |      |      |      |       /       /      /      /      /     _  Bhaduri (Set I, Lambda-N) [PR 155 (1967) 1671]
  //                       |      |      |      |      |       |      |      |      |      /       /      /      /      /     /
  //                       
  Double_t a_s[15] = { -2.03, -2.40, -2.86, -0.77, -0.97,  -1.28, -1.82, -2.24, -2.68, -1.56,  -2.04, -0.56, -0.40, -2.29, -2.46 };   // S-wave signlet
  Double_t r_s[15] = {  3.66,  3.15,  2.91,  6.09,  5.09,   4.22,  3.52,  3.24,  3.07,  1.43,   0.64,  7.77, 12.28,  3.14,  3.87 };

  Double_t a_t[15] = { -1.84, -1.84, -1.24, -2.15, -2.09,  -2.07, -1.94, -1.83, -1.67, -1.59,  -1.33, -1.91, -2.12, -1.77, -2.07 };    // S-wave triplet
  Double_t r_t[15] = {  3.32,  3.37,  3.33,  2.71,  2.80,   2.86,  3.01,  3.14,  3.34,  3.16,   3.91,  2.43,  2.57,  3.25,  4.50 };

  char* model_name[15] = { "Nijm_D", "Nijm_F", "NSC89", "NSC97a", "NSC97b", "NSC97c", "NSC97d", "NSC97e", "NSC97f", "Jue_A", "Jue_Ac", "Jue_B", "Jue_Bc", "Verma", "Bhaduri" };


  int lcol[15] = { 1, 2, 4, 3, 6,  7, 8, 9, 1, 2,  3, 4, 6, 7, 8 };
  int lsty[15] = { 1, 1, 1, 1, 1,  1, 1, 1, 2, 2,  2, 2, 2, 2, 2 };

  // ___________________________ enhancement factor

  TF1* en1s_fac[15];
  TF1* en1t_fac[15];
  TF1* en1_fac[15];

  TF1* en2_fac[15];

  char fname[20] = "";
  for ( int i=0; i<15; i++ ) {
    sprintf( fname, "en1s_fac_%0d", i );
    en1s_fac[i] = new TF1( fname, ef_1s, 0.0, 1000., 2 );
    en1s_fac[i] -> SetParameter( 0, a_s[i] );
    en1s_fac[i] -> SetParameter( 1, r_s[i] );
    en1s_fac[i] -> SetLineColor( lcol[i] );
    en1s_fac[i] -> SetLineStyle( lsty[i] );
    en1s_fac[i] -> SetLineWidth( 2 );

    sprintf( fname, "en1t_fac_%0d", i );
    en1t_fac[i] = new TF1( fname, ef_1s, 0.0, 1000., 2 );
    en1t_fac[i] -> SetParameter( 0, a_t[i] );
    en1t_fac[i] -> SetParameter( 1, r_t[i] );
    en1t_fac[i] -> SetLineColor( lcol[i] );
    en1t_fac[i] -> SetLineStyle( lsty[i] );
    en1t_fac[i] -> SetLineWidth( 2 );

    sprintf( fname, "en1_fac_%0d", i );
    en1_fac[i] = new TF1( fname, ef_1, 0.0, 1000., 4 );
    en1_fac[i] -> SetParameter( 0, a_s[i] );
    en1_fac[i] -> SetParameter( 1, r_s[i] );
    en1_fac[i] -> SetParameter( 2, a_t[i] );
    en1_fac[i] -> SetParameter( 3, r_t[i] );
    en1_fac[i] -> SetLineColor( lcol[i] );
    en1_fac[i] -> SetLineStyle( lsty[i] );
    en1_fac[i] -> SetLineWidth( 2 );

    sprintf( fname, "en2_fac_%0d", i );
    en2_fac[i] = new TF1( fname, ef_2, 0.0, 1000., 4 );
    en2_fac[i] -> SetParameter( 0, a_s[i] );
    en2_fac[i] -> SetParameter( 1, r_s[i] );
    en2_fac[i] -> SetParameter( 2, a_t[i] );
    en2_fac[i] -> SetParameter( 3, r_t[i] );
    en2_fac[i] -> SetLineColor( lcol[i] );
    en2_fac[i] -> SetLineStyle( lsty[i] );
    en2_fac[i] -> SetLineWidth( 2 );

  }


  TLegend* tlg1 = new TLegend(0.8, 0.3, 0.98, 0.98, "Model" );
  for ( int i=0; i<15; i++ ) {
    tlg1 -> AddEntry( en1s_fac[i], model_name[i], "l" );
  }


  //____________________________________________________________ canvas

  TCanvas* c1 = new TCanvas( "c1", "c1", 600, 800 );
  TPad* p1[3];
  p1[0] = new TPad( "p1_0", "p1_0",  0.0, 2./3., 1.0, 3./3. );
  p1[1] = new TPad( "p1_1", "p1_1",  0.0, 1./3., 1.0, 2./3. );
  p1[2] = new TPad( "p1_2", "p1_2",  0.0, 0./3., 1.0, 1./3. );

  for ( int i=0; i<3; i++ ) {
    p1[i] -> SetTopMargin( 0.10);
    p1[i] -> SetBottomMargin( 0.20);
    p1[i] -> SetLeftMargin( 0.25);
    p1[i] -> SetRightMargin( 0.10);
    p1[i] -> Draw();
  }

  TH1F* hframe[3];


  TPostScript* ps1 = new TPostScript( "../ps/enhancement_factor.ps", -111 );


  // ____________________________________________________________ page 1

  ps1 -> NewPage();

  p1[0] -> cd();
  p1[0] -> SetGridy();

  hframe[0] = p1[0] -> DrawFrame( 0.0,  0.0,  1000., 20.0  );
  hframe[0] -> GetXaxis() -> SetTitle( "Relative momentum of #Lambda-n [MeV/c]" );
  hframe[0] -> GetYaxis() -> SetTitle( "Enhancement factor (^{1}S_{0})" );

  for ( int i=0; i<15; i++ ) {
    en1s_fac[i] -> Draw( "same" );
  }
  tlg1 -> Draw();


  p1[1] -> cd();
  p1[1] -> SetGridy();

  hframe[1] = p1[1] -> DrawFrame( 0.0,  0.0,  1000., 20.0  );
  hframe[1] -> GetXaxis() -> SetTitle( "Relative momentum of #Lambda-n [MeV/c]" );
  hframe[1] -> GetYaxis() -> SetTitle( "Enhancement factor (^{3}S_{1})" );

  //  for ( int i=0; i<15; i++ ) {
  //    en1t_fac[i] -> Draw( "same" );
  //  }

  en1t_fac[13] -> Draw( "same" );
  en1t_fac[9]  -> Draw( "same" );
  en1t_fac[11] -> Draw( "same" );  

  p1[2] -> cd();
  p1[2] -> SetGridy();

  hframe[2] = p1[2] -> DrawFrame( 0.0,  0.0,  1000., 10.0  );
  hframe[2] -> GetXaxis() -> SetTitle( "Relative momentum of #Lambda-n [MeV/c]" );
  hframe[2] -> GetYaxis() -> SetTitle( "Enhancement factor (0.25 * ^{1}S_{0} + 0.75 * ^{3}S_{1})" );



  en1_fac[13] -> Draw( "same" );
  en1_fac[9] -> Draw( "same" );
  en1_fac[11] -> Draw( "same" );
  //  for ( int i=0; i<15; i++ ) {
  //    en1_fac[i] -> Draw( "same" );
  //  }


  c1 -> Update();


  // ____________________________________________________________ page 2

  ps1 -> NewPage();

  p1[0] -> cd();

  hframe[0] = p1[0] -> DrawFrame( 0.0,  0.0,  500., 10.0  );
  hframe[0] -> GetXaxis() -> SetTitle( "Relative momentum of #Lambda-n [MeV/c]" );
  hframe[0] -> GetYaxis() -> SetTitle( "Enhancement factor" );

  en2_fac[13] -> Draw( "same" );
  en2_fac[9] -> Draw( "same" );
  en2_fac[11] -> Draw( "same" );

//  for ( int i=0; i<15; i++ ) {
//    en2_fac[i] -> Draw( "same" );
//  }
  




  c1 -> Update();
  ps1 -> Close();

}
