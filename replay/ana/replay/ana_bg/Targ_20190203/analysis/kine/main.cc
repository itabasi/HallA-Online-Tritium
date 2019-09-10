/*
  main.cc
  
  Toshiyuki Gogami, October 16, 2017
*/

#include "CalcKinematics_eek.hh"
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TObjArray.h>
using namespace std;

int main(int argc, char* argv[]){
  
  //TApplication app("app",&argc, argv);
  
  CalcKinematics_eek* kine = new CalcKinematics_eek();
  ifstream* ifs = new ifstream("input.dat");
  int genflag = 0; // Generation flag
  char input_fname[500];
  char out_fname[500];
  double nmax = 0.0;
  double pcent, pbite;
  double thcent, thbite; // theta_ep
  double phcent, phbite; // phi_ep
  *ifs >> input_fname;
  *ifs >> out_fname;
  *ifs >> genflag;
  *ifs >> nmax;
  *ifs >> pcent >> pbite;
  *ifs >> thcent >> thbite;
  *ifs >> phcent >> phbite;
  //cout << fname << endl;
  
  // 1: p(e,e'K+)Lambda
  kine->SetFlag(genflag); 
  
  TFile* f1 = new TFile(input_fname);
  //TFile* f1 = new TFile("test_Cu.root");
  TTree* t1 = (TTree*)f1->Get("tree");
  Int_t eventid;
  Double_t xBeam, yBeam, zBeam;
  Double_t pVD[20];
  //Double_t pVD;
  ///float pVD[20];
  //Double_t pVD[20];
  Double_t pBeam;
  Double_t pep, theta_ep, phi_ep;
  Int_t eleflag[20];
  Double_t charge[20];
  t1->SetBranchAddress("eventid", &eventid);
  t1->SetBranchAddress("xBeam",   &xBeam);
  t1->SetBranchAddress("yBeam",   &yBeam);
  t1->SetBranchAddress("zBeam",   &zBeam);
  t1->SetBranchAddress("pVD",     &pVD);
  t1->SetBranchAddress("pBeam",   &pBeam);
  t1->SetBranchAddress("eleflag", &eleflag);
  t1->SetBranchAddress("charge",  &charge);
  double ent = t1->GetEntries();
  double nloop = ent;
  
  TRandom3* rnd1 = new TRandom3();
  rnd1->SetSeed();
  
  TRandom3* rnd2 = new TRandom3();
  rnd2->SetSeed();
  
  TRandom3* rnd3 = new TRandom3();
  rnd3->SetSeed();
  
  double thmin, thmax;
  thmin = cos(thcent-thbite);
  thmax = cos(thcent+thbite);
  
  double xp_k, yp_k;
  double xp_ep, yp_ep;
  double pk;
  double xt, yt, zt;
  double mhyp, mm, mm2;
  double Emom0;
  double Qsq;
  double omega;
  double th_gk;
  
  TH1F* hmm = new TH1F("hmm","",1000,-5.,5.);
  TH1F* hmm2 = (TH1F*)hmm->Clone("hmm2");
  TH1F* h_theta_ep = new TH1F("h_theta_ep","",100,0.0,0.5);
  //TH1F* h_theta_ep = new TH1F("h_theta_ep","",100,0.9,1.);
  
  char newROOTfname[600];
  sprintf(newROOTfname,"%s_seedcheck.root",out_fname);
  TFile* fnew = new TFile(newROOTfname,"recreate");
  TTree* tnew = new TTree("tree","ROOT file for check");
  tnew->Branch("evID", &eventid,"evID/I");
  tnew->Branch("BMom0",&pBeam,  "BMom0/D");
  tnew->Branch("BMom", &pVD[0], "BMom/D");
  tnew->Branch("EMom", &pep,    "EMom/D");
  tnew->Branch("EXpt", &xp_ep,  "EXpt/D");
  tnew->Branch("EYpt", &yp_ep,  "EYpt/D");
  tnew->Branch("KMom", &pk,     "KMom/D");
  tnew->Branch("KXpt", &xp_k,   "KXpt/D");
  tnew->Branch("KYpt", &yp_k,   "KYpt/D");
  tnew->Branch("mm",   &mm,     "mm/D");
  tnew->Branch("mm2",   &mm2,     "mm2/D");
  tnew->Branch("mhyp", &mhyp,   "mhyp/D");
  tnew->Branch("Qsq",  &Qsq,    "Qsq/D");
  tnew->Branch("omega",&omega,  "omega/D");
  tnew->Branch("th_gk",&th_gk,  "th_gk/D");
  
  char Seedfname[600];
  sprintf(Seedfname,"%s.seed",out_fname);
  ofstream* dragon_seed = new ofstream(Seedfname);
  int totn = 0;
  
  if(nloop>nmax) nloop = nmax;
  for(int i=0 ; i<nloop ; i++){
    mm    = -2222.0;
    Qsq   = -2222.0;
    omega = -2222.0;
    th_gk = -2222.0;
    pk = -2222.0;
    xp_k = -2222.0;
    yp_k = -2222.0;
    for(int j=0 ; j<20 ; j++){
      pVD[j] = -2222.0;
      charge[j] = -10.0;
    }
    t1->GetEntry(i);
    //cout << pVD[0] << endl;
    ///cout << eleflag << endl;
    //if(eleflag==1){
    //cout << type(pVD[0]) << endl;
    //cout << charge[0] << " " << pVD[0]  << endl;
    //cout << charge[0] << " " << pVD << endl;
    //cout << eleflag[0] << endl;
    //if(charge[0]==-1 && pVD[0]>300.){
    if(charge[0]==-1 && pVD[0]>300.0){
    //if(charge[0]==-1 && eleflag[0]==1 && pVD[0]>300.0){
      kine->SetEventID(eventid);
      kine->SetGenPos(xBeam, yBeam, zBeam);
      pep = rnd1->Uniform(pcent*(1.-pbite/100.),
			  pcent*(1.+pbite/100.));
      
      theta_ep = acos( rnd2->Uniform(thmin,thmax) );
      h_theta_ep -> Fill(theta_ep);
      //h_theta_ep -> Fill(cos(theta_ep));
      
      phi_ep = rnd3->Uniform(phcent-phbite,
			     phcent+phbite);
      //      cout << xBeam << " " << yBeam << " " << zBeam << " "
      //	   << pBeam << " " << pVD << " " 
      //	   << pep << " " << theta_ep << " " 
      //	   << phi_ep<< endl;
      //      
      
      // ==================================== //
      // ===== Kinematics calculation ======= //
      // ==================================== //
      //kine->Calc(pVD[0], pep, theta_ep, phi_ep);
      //kine->Calc(pVD, pep, theta_ep, phi_ep);
      int ok = 0;
      ok = kine->Calc(pVD[0], pep, theta_ep, phi_ep);	
      // ==================================== //
      
      TVector3 veck = kine->GetKMom();
      pk   = veck.X();
      xp_k = veck.Y();
      yp_k = veck.Z();
      TVector3 vecep = kine->GetEpMom();
      //pep  = vecep.X();
      xp_ep = vecep.Y();
      yp_ep = vecep.Z();
      //cout << pep << " " << vecep.X() << endl;
      TVector3 vecgenpos = kine->GetGenPos();
      xt = vecgenpos.X();
      yt = vecgenpos.Y();
      zt = vecgenpos.Z();

      if(ok==1){
	cout << eventid << " " 
	     << xt << " " << yt << " " << zt << " "
	     << pBeam << " " << pVD[0] << " "
	     << pep<< " " << xp_ep<< " " << yp_ep << " "
	     << pk << " " << xp_k << " " << yp_k  << endl;
	
	//dragon_seed->precision(10);
	dragon_seed->precision(6);
	*dragon_seed << eventid << " "
		     << xt << " " << yt << " " << zt << " ";
	dragon_seed->precision(10);
	*dragon_seed << pBeam << " " << pVD[0] << " "
		     << pep<< " " << xp_ep<< " " << yp_ep << " "
		     << pk << " " << xp_k << " " << yp_k  << endl;
	
	mm    = kine->CheckMM();
	mm2   = kine->CheckMM2();
	//cout << " ------- " << mm << " " << mm2 << endl;
	mhyp = kine->GetAssumedMhyp();
	omega = kine->GetOmega();
	Qsq = kine->GetQsquare();
	th_gk = kine->GetThetaGammaK();
	
	hmm->Fill(mm-mhyp);
	hmm2->Fill(mm2-mhyp);
	tnew->Fill();
	totn++;
      }
      
    }
  }
  dragon_seed->close();
  
  tnew->Write();
  TObjArray h(1);
  h.Add(hmm);
  h.Add(hmm2);
  h.Add(h_theta_ep);
  h.Write();
  fnew->Close();

  cout << endl;
  cout << " input:  " << input_fname << endl;
  cout << " flag:   " << kine->GetFlag() << endl;
  cout << "         "  << "[1:p(e,e'K+)L, 300:3H(e,e'K+)nnL, 4000:40Ca(e,e'K+)40LK, 4800:48Ca(e,e'K+)48LK ]" << endl;
  cout << " output: " << out_fname << ".seed (.root)" << endl;
  cout << " total:  " << totn << " events out of "
       <<  ent << " events " << endl;
  
  
  //TCanvas* c1 = new TCanvas("c1","c1");
  //hmm->Draw();
  //c1->WaitPrimitive();
      
  return 0;
}
