//------------------//
// VP Flux Integral //
//------------------//
// Okuyama code : 2021/01/28


#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <time.h>
#include <string>

#include "TStyle.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TCut.h"
#include "TLine.h"
#include "TRandom.h"
#include "TLorentzVector.h"
using namespace std;
void SetTitle( TH2D*, const char*, const char*, const char*);
void SetTitle( TH1D*, const char*, const char*, const char*);

static const double PI = 4.*atan(1.0);

double vpflux_lab(double *par){
	double Einc  = 4.3;//[GeV]
	//double Escat = 2.1;//[GeV]
	double Escat = par[0];//[GeV]
	double theta = par[1];	

	double Me=511.*pow(10.,-6.);//[GeV/c^2]
	double Mp=0.9382720;//[GeV/c^2]
	double Qsq=2*Einc*Escat*(1-cos(theta));
	double omega = Einc - Escat;
	double q2=Qsq+omega*omega;
	double kg=omega-Qsq/(2*Mp);
	double eps=1/(1+2*(q2/Qsq)*tan(theta/2)*tan(theta/2));
//if(theta>0.230&&theta<0.231){
//cout<<"Qsq="<<Qsq<<endl;
//cout<<"q2="<<q2<<endl;
//cout<<"kg="<<kg<<endl;
//cout<<"eps="<<eps<<endl;
//}

	double vpflux=Escat*kg/(137*2*PI*PI*Einc*Qsq*(1-eps));
	return vpflux;
}
//TF1* func3 = new TF1("func3",vpflux_lab, 0.001, 0.3,1);
//func3->SetNpx(600);
//func3->SetParameter(0,2.1);
//func3->SetLineColor(kGreen);
//func3->SetLineWidth(4);
//func3->Draw("same");

int main(int argc, char** argv){
	//////////////////////
	// command argument //
	//////////////////////
	int option;
	//string filename = "1H_kaon";
	//label
	//string filename = "LHRS_new";//tehta_width=0.067 rad
	//string filename = "LHRS_big";//thetawidth=0.1 rad
	string filename = "LHRS_Lambda_grnd";//thetawidth=0.1 rad
	//string filename = "LHRS_datafit";//thetawidth=0.1 rad
	//string filename = "LHRS_SigmaZ_grnd";//thetawidth=0.1 rad
	//string filename = "BOTH";
	//string filename = "RHRS";
	bool BatchFlag = false;
	bool PDFFlag = true;
	while((option=getopt(argc, argv, "f:bp"))!=-1){
		switch(option){
			case 'f':
				filename = optarg;
				break;
			case 'b':
				cout << "batch mode" << endl;
				BatchFlag = true;
				break;
			case 'p':
				PDFFlag = true;
				break;
			case 'h':
				cout<<"-f : input root filename"<<endl;
				cout<<"-b : execute in batch mode"<<endl;
				cout<<"-p : print pdf file"<<endl;
				return 0;
				break;
			case '?':
				cout << "unknown option: " << option << endl;
				return 0;
				break;
			default:
				cout<<"type -h to see help!!"<<endl;
				return 0;
				break;
		}
	}
  TApplication theApp("App", &argc, argv);
	ostringstream RootFile;
	ostringstream PDFFile;
	RootFile << "../" << filename << ".root";	
	PDFFile << "./vpflux.pdf";	

	///////////////////////
	// General condition //
	///////////////////////
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleSize(0.05, "X");
  gStyle->SetTitleSize(0.05, "Y");
  gStyle->SetTitleOffset(0.9, "X");
  gStyle->SetTitleOffset(0.9, "Y");
	gStyle->SetOptStat(0);
	gStyle->SetPadGridX(0);
	gStyle->SetPadGridY(0);

	FILE *fp;
	const int SizeOfBuffer = 32;
    char str[SizeOfBuffer];
	double centralmom, centraltheta, thetawidth;
	double centralphi, phiwidth;


	
	centraltheta=13.2*PI/180.;
	//thetawidth=0.067;
	thetawidth=0.1;
	centralphi=0.;
	phiwidth=2*PI;
	cout<< "central theta = " << centraltheta*180/PI << " [deg]" << endl;
	cout<< "theta width= " << thetawidth*180/PI << " [deg]" << endl;
	centralphi = 0.;

	TFile *f = new TFile(RootFile.str().c_str());
	TTree *t = (TTree*)f->Get("SNT");
	TTree *tree = (TTree*)f->Get("SNT");

	double thetamin = 1.*(1.*centraltheta - thetawidth);
	double thetamax = 1.*(1.*centraltheta + thetawidth);
	double phimin = 1.*(centralphi - phiwidth);
	double phimax = 1.*(centralphi + phiwidth);
	double Domega = 2.*PI* (1-cos(thetawidth))*1000.; // [msr]
	double momwidth = 0.5;
	double MAX = 0.01;
	double volume = Domega*MAX/1000.; 
	cout<< "thetamin = " << thetamin*180/PI << " [deg]" <<endl;
	cout<< "thetamax = " << thetamax*180/PI << " [deg]" <<endl;
	cout<< "phimin = " << phimin*180/PI << " [deg]" <<endl;
	cout<< "phimax = " << phimax*180/PI << " [deg]" <<endl;
	cout<< "Domega = " << Domega << " [msr]" <<endl;
	cout<< "volume = " << volume <<endl;


	int bin_mom = 150;
	double min_mom = 1.9;//nnL
	double max_mom = 2.3;//nnL
	int bin_th = 150;
	//double min_th = 0.10;
	//double max_th = 0.35;
	double min_th = 0.95;//cos
	double max_th = 1.0;//cos
	double min_ph = -0.4292;//default=0.57
	double max_ph = 0.4292;
	int bin_2D_mom = 50;
	int bin_2D_th = 50;
	int bin_2D_ph = 50;

	TH1D *h_mom_gen = new TH1D( "h_mom_gen", "", bin_mom, min_mom, max_mom);
	TH1D *h_mom_result = new TH1D( "h_mom_result", "", bin_mom, min_mom, max_mom);
	TH1D *h_sa_mom_result = new TH1D("h_sa_mom_result", "", bin_mom, min_mom, max_mom);
	//t->Project("h_mom_result","Lp_orig");
	
	//TH1D *h_mom_gen_cm = new TH1D( "h_mom_gen_cm", "", bin_mom, min_mom_cm, max_mom_cm);
	//TH1D *h_mom_result_cm = new TH1D( "h_mom_result_cm", "", bin_mom, min_mom_cm, max_mom_cm);
	//TH1D *h_sa_mom_result_cm = new TH1D("h_sa_mom_result_cm", "", bin_mom, min_mom_cm, max_mom_cm);

	TH1D *h_vp_mom = new TH1D( "h_vp_mom", "", bin_mom, min_mom, max_mom);
	TH1D *h_vp_mom2 = new TH1D( "h_vp_mom2", "w/ HRS-R Acceptance", bin_mom, min_mom, max_mom);
	TH1D *h_vp_mom_result = new TH1D("h_vp_mom_result", "", bin_mom, min_mom, max_mom);
	TH1D *h_vp_mom_result2 = new TH1D("h_vp_mom_result2", "", bin_mom, min_mom, max_mom);//w/ RHRS Acceptance
	TH1D *h_vp_mom_result3 = new TH1D("h_vp_mom_result3", "", bin_mom, min_mom*1000., max_mom*1000.);//w/ RHRS Acceptance
	SetTitle(h_vp_mom_result3, "", "Momentum [MeV/c]", "Integrated VP Flux [#gamma#lower[-0.2]{#scale[0.5]{#bullet}}MeV^{-1}#lower[-0.2]{#scale[0.5]{#bullet}}electron^{-1}]");
	//h_vp_mom_result3->GetXaxis()->SetNdivisions(505, kFALSE);
	//h_vp_mom_result3->SetMarkerColor(kBlack);
	h_vp_mom_result3->SetLineColor(kBlack);
	h_vp_mom_result3->SetLineWidth(2);
	TH2D *h_th_mom = new TH2D( "h_th_mom", "", bin_2D_th, 0.,20.,bin_2D_mom, 1.95, 2.25);
	TH2D *h_momL_momR = new TH2D( "h_momL_momR", "", 50, 1.73, 1.93 ,50., 1.95, 2.25);
	TH1D *h_momL = new TH1D( "h_momL", "", 200,1.95,2.25);
	TH1D *h_th_cm = new TH1D( "h_th_cm", "", 100.,0.,20.);
	TH1D *h_th_cm2 = new TH1D( "h_th_cm2", "", 100.,0.,20.);

	int n1, n2;
	double val;
	double err;
	const double min_mom_gen = h_mom_gen->GetBinCenter(h_mom_gen->FindFirstBinAbove()) - 0.05;
	const double max_mom_gen = h_mom_gen->GetBinCenter(h_mom_gen->FindLastBinAbove()) + 0.05;
	const double min_mom_bin = h_mom_gen->FindBin(min_mom_gen);
	const double max_mom_bin = h_mom_gen->FindBin(max_mom_gen);
	cout <<"mom range: "<< min_mom_gen << " " << max_mom_gen << endl;
	cout <<"mom_bin: "<< min_mom_bin << " " << max_mom_bin << endl;
//----------------------------//
//------------Fill------------//
//----------------------------//
	int ENum=0;
	float R_mom=0.;
	float R_th=0.;
	float R_ph=0.;
	float L_mom=0.;
	float L_th=0.;
	float L_ph=0.;
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Rp_gen",1);tree->SetBranchAddress("Rp_gen",&R_mom);
	tree->SetBranchStatus("Rth_gen",1);tree->SetBranchAddress("Rth_gen",&R_th);
    tree->SetBranchStatus("Rph_gen"    ,1);tree->SetBranchAddress("Rph_gen"    ,&R_ph     );
	tree->SetBranchStatus("Lp_gen",1);tree->SetBranchAddress("Lp_gen",&L_mom);
	tree->SetBranchStatus("Lth_gen",1);tree->SetBranchAddress("Lth_gen",&L_th);
    tree->SetBranchStatus("Lph_gen"    ,1);tree->SetBranchAddress("Lph_gen"    ,&L_ph     );
    ENum = tree->GetEntries();
  for(int i=0;i<ENum;i++){
    tree->GetEntry(i);
    if(i%100000==0)cout<<i<<" / "<<ENum<<endl;
	
	L_mom /= 1000.;//MeV-->GeV
	R_mom /= 1000.;//MeV-->GeV
	//R_mom = 1.85;//GeV/c
	//R_th=0.;
	//R_ph=0.;
	double Einc  = 4.318;//[GeV]
	//double Escat = 2.1;//[GeV]
	double Escat = L_mom;//[GeV]
	double theta = L_th+centraltheta;	

	double Me=511.*pow(10.,-6.);//[GeV/c^2]
	double Mp=0.9382720;//[GeV/c^2]
	double MK=0.494677;//[GeV/c^2]
    double ML = 1.115683;//[GeV/c2]
    double mh = ML;//hypernuclei
    double mt = Mp;//target mass
	double Qsq=2*Einc*Escat*(1-cos(theta));
	double deltaE = Einc - Escat;
	double q2=Qsq+deltaE*deltaE;
	double kg=deltaE-Qsq/(2*Mp);
	double eps=1/(1+2*(q2/Qsq)*tan(theta/2)*tan(theta/2));
		//h_mom_gen->SetBinContent(h_mom_gen->FindBin(mom),deltaE);
		//h_mom_result->SetBinContent(h_mom_gen->FindBin(mom),mom);
		for(int i=0;i<bin_mom;i++){
		//h_mom_gen->SetBinContent(i,26147*(1./189.)*1000.*(max_mom-min_mom)/bin_mom);//first try
		//h_mom_gen->SetBinContent(i,2583235*(1./189.)*1000.*(max_mom-min_mom)/bin_mom);//LHRS 1,000,000 (2020/10/4)
		//h_mom_gen->SetBinContent(i,2582007*(1./189.)*1000.*(max_mom-min_mom)/bin_mom);//LHRS_new 1,000,000 (2020/10/17)// true density
		//h_mom_gen->SetBinContent(i,5688035*(1./189.)*1000.*(max_mom-min_mom)/bin_mom);//LHRS_big 1,000,000 (2020/10/27)
		h_mom_gen->SetBinContent(i,5695419*(1./189.)*1000.*(max_mom-min_mom)/bin_mom);//LHRS_Lambda_grnd 1,000,000 (2020/11/22)
		//h_mom_gen->SetBinContent(i,5990331*(1./189.)*1000.*(max_mom-min_mom)/bin_mom);//LHRS_datafit 1,000,000 (2021/1/6)
		//h_mom_gen->SetBinContent(i,5696881*(1./189.)*1000.*(max_mom-min_mom)/bin_mom);//LHRS_SigmaZ_grnd 1,000,000 (2020/11/22)
		//h_mom_gen->SetBinContent(i,2549979*(1./189.)*1000.*(max_mom-min_mom)/bin_mom);//RHRS 1,000,000 (2020/10/4)
		}
	double vpflux=Escat*kg/(137*2*PI*PI*Einc*Qsq*(1-eps));
	double k = MAX*gRandom->Uniform();
		//h_mom_result->SetBinContent(h_mom_gen->FindBin(mom),mom);
//cout<<"Qsq="<<Qsq<<endl;
//cout<<"q2="<<q2<<endl;
//cout<<"kg="<<kg<<endl;
//cout<<"eps="<<eps<<endl;

	//cout<<"vpflux vs k = "<<vpflux<<" : "<<k<<endl;
	
	    //===== Right Hand Coordinate ====//
	    //th and phi are originally meant tan(theta) and tan(phi),
	    //so, they should not be treated like tan(R_tr_tr_th) //2020.6.30 Okuyama

		double B_mom = Einc;
		
	    double R_pz = R_mom/sqrt(1.0*1.0 + pow(R_th, 2.0) + pow( R_ph,2.0));
	    double R_px = R_pz * ( R_th );
	    double R_py = R_pz * ( R_ph );

	    double L_pz = L_mom/sqrt(1.0*1.0 + pow(L_th, 2.0) + pow(L_ph,2.0));
	    double L_px = L_pz * ( L_th );
	    double L_py = L_pz * ( L_ph );


	    double B_E =sqrt(B_mom*B_mom + Me*Me);
	    double R_E =sqrt(R_mom*R_mom + MK*MK);
	    double L_E =sqrt(L_mom*L_mom + Me*Me);


		TLorentzVector L_4vec;//Left
		TLorentzVector R_4vec;//Right
		TLorentzVector B_4vec;//Beam
		TLorentzVector T_4vec;//Target
		TLorentzVector G_4vec;//Gamma (Virtual Photon)
		L_4vec.SetPxPyPzE(L_px, L_py, L_pz, L_E);
        R_4vec.SetPxPyPzE(R_px, R_py, R_pz, R_E);
        B_4vec.SetPxPyPzE(0.0 ,  0.0,B_mom, B_E);
        T_4vec.SetPxPyPzE(0.0 ,  0.0,  0.0,  mt);




	    double pL    = L_mom;//GeV
	    double pR    = R_mom;//GeV
		theta = L_th;
		double theta_R = R_th;
		double phi = L_ph;
		double phi_R = R_ph;
		double phi0=13.2*PI/180;//rad
		double phi_L = L_4vec.Phi();//LHRS frame
		double phi_RHRS = R_4vec.Phi();//RHRS frame
	    L_4vec.RotateX( -13.2/180.*PI );
	    R_4vec.RotateX(  13.2/180.*PI );
        double mass,mm;
		TLorentzVector Missing;
		Missing = B_4vec + T_4vec - L_4vec - R_4vec;
		mass = Missing.M();
        //mass = sqrt( (Ee + mt - L_E - R_E)*(Ee + mt - L_E - R_E)-(B_v - L_v - R_v)*(B_v - L_v - R_v) );
	    mm=mass - mh;//shift by ML

		double theta_ee = L_4vec.Theta();
		//double theta_ee = acos((-phi*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		double theta_ek = R_4vec.Theta();
		//double theta_ek = acos((phi_R*sin(phi0)+cos(phi0))/(sqrt(1+theta*theta+phi*phi)));//original frame
		double phi_ee = L_4vec.Phi();//original frame
//cout<<"phi_ee="<<phi_ee<<endl;
		double phi_ek = R_4vec.Phi()+2*PI;//original frame
//cout<<"phi_ek="<<phi_ek<<endl;

		G_4vec = B_4vec - L_4vec;
		double mom_g=sqrt(G_4vec.Px()*G_4vec.Px()+G_4vec.Py()*G_4vec.Py()+G_4vec.Pz()*G_4vec.Pz());
		Qsq = G_4vec.M()*G_4vec.M();
		double phi_g = G_4vec.Phi()+2*PI;
		double theta_g = G_4vec.Theta();
//cout<<"theta_gk="<<(theta_ek-theta_g)*180./PI<<endl;
		double theta_gk_lab = G_4vec.Angle(R_4vec.Vect());
//cout<<"theta_gk_lab(TLorentz)="<<theta_gk_lab*180./PI<<endl;
		double omega=G_4vec.E();
		double pY = sqrt((omega+Mp-sqrt(pR*pR+MK*MK))*(omega+Mp-sqrt(pR*pR+MK*MK))-ML*ML);
		double W = sqrt((omega+Mp)*(omega+Mp)-mom_g*mom_g);
		//double theta_gk_lab_test = acos((mom_g*mom_g+pR*pR-pY*pY)/(2.*pR*mom_g));
		double theta_eg_lab = acos((B_E-L_E*(1.-Qsq/2./B_E/L_E))/mom_g);
//cout<<"theta_eg_lab="<<theta_eg_lab*180./PI<<endl;
		double theta_gk_lab_test = 13.2*PI/180.-theta_eg_lab;
//cout<<"theta_gk_lab="<<theta_gk_lab_test*180./PI<<endl;
		double beta=mom_g/(omega+Mp);
		double ER=sqrt(pR*pR+MK*MK);
		double gamma=1./sqrt(1-beta*beta);
		double theta_gk_cm_test = atan((pR*sin(theta_gk_lab_test))/(-1.*gamma*beta*ER+gamma*pR*cos(theta_gk_lab_test)));
//cout<<"theta_gk_cm="<<theta_gk_cm_test*180./PI<<endl;
	
		TVector3 boost;
		TLorentzVector GT_4vec;
		GT_4vec=G_4vec+T_4vec;
		boost=GT_4vec.BoostVector();
		R_4vec.Boost(-boost);
		L_4vec.Boost(-boost);
		B_4vec.Boost(-boost);
		double theta_gk_cm = G_4vec.Angle(R_4vec.Vect());
		double pR_cm=sqrt(R_4vec.Px()*R_4vec.Px()+R_4vec.Py()*R_4vec.Py()+R_4vec.Pz()*R_4vec.Pz());
		double pL_cm=sqrt(L_4vec.Px()*L_4vec.Px()+L_4vec.Py()*L_4vec.Py()+L_4vec.Pz()*L_4vec.Pz());
		double pB_cm=sqrt(B_4vec.Px()*B_4vec.Px()+B_4vec.Py()*B_4vec.Py()+B_4vec.Pz()*B_4vec.Pz());

//cout<<"theta_gk_cm(TLorentz)="<<theta_gk_cm*180./PI<<endl;
		double p_cm=sqrt(GT_4vec.Px()*GT_4vec.Px()+GT_4vec.Py()*GT_4vec.Py()+GT_4vec.Pz()*GT_4vec.Pz());
		double E_cm = GT_4vec.E();
		double ER_cm=sqrt(pR_cm*pR_cm+MK*MK);
//cout<<"beta="<<beta<<endl;
//cout<<"gamma="<<gamma<<endl;

		double labtocm = (gamma*pR_cm*pR_cm*(pR_cm*cos(theta_gk_cm)+beta*ER_cm))/(pow(sqrt(pR_cm*pR_cm*sin(theta_gk_cm)*sin(theta_gk_cm)+gamma*gamma*(pR_cm*cos(theta_gk_cm)+beta*ER_cm)*(pR_cm*cos(theta_gk_cm)+beta*ER_cm)),3.));
//cout<<"labtocm="<<labtocm<<endl;
		double tan_lab1 = sin(theta_gk_cm)/(gamma*(cos(theta_gk_cm)+beta*sqrt(MK*MK+pR_cm*pR_cm)/pR_cm));
		double tan_lab2 = sin(theta_gk_cm)/(gamma*(cos(theta_gk_cm)+(omega*Mp-Qsq*Qsq)/(omega*Mp+Mp*Mp)));

//===Angle partition===//
		//h_mom_result->Fill(L_mom);
		//if(theta_gk_cm*180./PI>=6.&&theta_gk_cm*180./PI<10.)h_mom_result->Fill(L_mom);
		h_mom_result->Fill(L_mom);
		//h_mom_result->Fill(L_mom);
		h_th_mom->Fill(theta_gk_cm*180./PI,L_mom);
		h_momL_momR->Fill(R_mom,L_mom);
		if(R_mom>1.74&&R_mom<1.91)h_momL->Fill(L_mom);
		h_th_cm->Fill(theta_gk_cm*180./PI);
		if(L_mom<2.12)h_th_cm2->Fill(theta_gk_cm*180./PI);
	 if(vpflux>k){//Full
		h_vp_mom->Fill(L_mom);
		//if(L_mom>2.08&&theta_gk_cm*180./PI>=6.&&theta_gk_cm*180./PI<10.){
		//if(L_mom<2.12&&theta_gk_cm*180./PI>=6.&&theta_gk_cm*180./PI<10.){//partition
//Partition templete
//L_mom>2.08&&Qsq<0.5 (Lambda, Qsq)
//L_mom<2.12&&Qsq<0.5 (Sigma0, Qsq)
//L_mom>2.08&&theta_gk_cm*180./PI<8 (Lambda, theta_gk_cm*180./PI)
//L_mom<2.12&&theta_gk_cm*180./PI<8 (Sigma0, theta_gk_cm*180./PI)
		//if(L_mom>2.010&&L_mom<2.160){//Full Mom. cut
		//if(L_mom>2.092&&L_mom<2.160){//L_mom>2.08(Lambda), L_mom<2.12(Simga0)
		//if(L_mom<2.108&&L_mom>2.010){//L_mom>2.08(Lambda), L_mom<2.12(Simga0)
		//change
		//if(L_mom<2.160&&L_mom>2.092&&theta_gk_cm*180./PI>=8.){//L_mom>2.082(Lambda), L_mom<2.108(Simga0)
		//if(L_mom<2.160&&L_mom>2.092&&Qsq>=0.5){//L_mom>2.082(Lambda), L_mom<2.108(Simga0)
		//if(L_mom<2.108&&L_mom>2.010&&theta_gk_cm*180./PI>=8.){//L_mom>2.082(Lambda), L_mom<2.108(Simga0)
		//if(L_mom<2.108&&L_mom>2.010&&Qsq>=0.5){//L_mom>2.082(Lambda), L_mom<2.108(Simga0)
		h_vp_mom2->Fill(L_mom);
		//if(L_mom>2.1)h_vp_mom2->Fill(L_mom);
		//}
	}
	}

	double vpflux_tot_mom=0.;
	double vpflux_tot_mom_err=0.;
	for(int i=0; i<bin_mom; i++){
		n1 = 0;
		n1 = (int)h_mom_gen->GetBinContent(i+1);

		// === all cuts ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_mom_result->GetBinContent(i+1);
		//cout<<"nbin_mom:"<<i<<endl;
		//cout<<"n1="<<n1<<endl;
		//cout<<"n2="<<n2<<endl;
		if(n1!=0 && n2!=0)val = Domega*(1.0*n2/n1);
		h_sa_mom_result->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_sa_mom_result->SetBinError(i+1, err);

		// === VP Flux ===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_vp_mom->GetBinContent(i+1);
		//cout<<"nbin_mom:"<<i<<endl;
		//cout<<"n1="<<n1<<endl;
		//cout<<"n2="<<n2<<endl;
		if(n1!=0 && n2!=0)val = volume*(1.0*n2/n1);//*1000.;//MeV->GeV
		h_vp_mom_result->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_vp_mom_result->SetBinError(i+1, err);

		// === VP Flux (w/ RHRS Acceptance)===
		n2 = 0;
		val = 0.;
		err = 0.;
		n2 = (int)h_vp_mom2->GetBinContent(i+1);
		//cout<<"nbin_mom:"<<i<<endl;
		//cout<<"n1="<<n1<<endl;
		//cout<<"n2="<<n2<<endl;
		if(n1!=0 && n2!=0)val = volume*(1.0*n2/n1);//*1000;//MeV->GeV
		h_vp_mom_result2->SetBinContent(i+1, val);
		if(n1!=0 && n2!=0)err = val * sqrt(1./n2 + 1./n1 - 2./sqrt(1.*n1*n2));
		h_vp_mom_result2->SetBinError(i+1, err);

		double bin_width=(max_mom-min_mom)/bin_mom;//GeV/c
		vpflux_tot_mom+=val*bin_width;
		vpflux_tot_mom_err+=sqrt(err*bin_width*err*bin_width);
	}
	cout<<"vpflux_tot_mom="<<vpflux_tot_mom<<endl;
	cout<<"vpflux_tot_mom_err="<<vpflux_tot_mom_err<<endl;

////////////////////
// Draw histgrams //
////////////////////
	TH1D *hframe;
	TLine *lmom = new TLine( centralmom, 0., centralmom, Domega);
	TLine *lth = new TLine( centraltheta, 0., centraltheta, Domega);
	lmom->SetLineColor(4);
	lth->SetLineColor(4);

	TCanvas *c1 = new TCanvas("c1", "Momentum Acceptance (Lab)");
	TCanvas *c2 = new TCanvas("c2", "Virtual Photon Flux");


	//======= Momentum Acceptance (Lab) ======
	c1->Divide(2,2);
	c1->cd(1);
//	hframe = (TH1D*)gPad->DrawFrame( min_mom_gen, 0., max_mom_gen, Domega + 0.5);
	//SetTitle(h_sa_mom_result, "Solid Angle vs. Momentum at Reference Plane", "Momentum [GeV/c]", "Solid Angle [msr]");
	SetTitle(h_sa_mom_result, "", "Momentum [GeV/c]", "Solid Angle [msr]");
//	h_sa_mom_result->GetXaxis()->SetNdivisions(506, kFALSE);
	h_sa_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_sa_mom_result->SetMarkerColor(kRed);
	h_sa_mom_result->SetLineColor(kRed);
	h_sa_mom_result->SetLineWidth(2);
//	h_sa_mom_result->Draw("same");
	h_sa_mom_result->Draw("e2");
	h_sa_mom_result->Draw("p same");
//	lmom->Draw("same");
	c1->cd(2);
	SetTitle(h_mom_gen, "", "Momentum [GeV/c]", "");
	h_mom_gen->GetXaxis()->SetNdivisions(505, kFALSE);
	h_mom_gen->Draw();
	c1->cd(3);
	SetTitle(h_mom_result, "", "Momentum [GeV/c]", "");
	h_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_mom_result->Draw();
	c1->cd(4);
	h_mom_gen->Draw();
	h_mom_result->SetLineColor(kAzure);
	h_mom_result->Draw("same");

	//======= VP Flux  ======
	c2->Divide(2,2);
	c2->cd(1);
	//SetTitle(h_vp_mom_result, "Integrated VP Flux vs. Momentum (w/ all Cuts)", "Momentum [GeV/c]", "Integrated VP Flux [/GeV]");
		//h_vp_mom_result->Scale(bin_mom);
	h_vp_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_vp_mom_result2->Draw();
	c2->cd(2);
	h_vp_mom_result->GetXaxis()->SetNdivisions(505, kFALSE);
	h_vp_mom_result->Draw();
		//h_vp_mom_result2->Scale(bin_mom);
	h_vp_mom_result2->GetXaxis()->SetNdivisions(505, kFALSE);
	h_vp_mom_result2->SetLineColor(kAzure);
	h_vp_mom_result2->Draw("");
	double ymax = (h_vp_mom_result->GetBinContent(h_vp_mom_result->GetMaximumBin()));
	TLine *RHRS_min = new TLine( 2.08, 0., 2.08, 1.1*ymax);
	TLine *RHRS_max = new TLine( 2.12, 0., 2.12, 1.1*ymax);
	RHRS_min->SetLineColor(kRed);
	RHRS_max->SetLineColor(kRed);
	RHRS_min->Draw("same");//Lambda
	//RHRS_max->Draw("same");//Sigma0
	c2->cd(3);
	h_vp_mom->Draw();

	TCanvas *c3 = new TCanvas("c3", "test");
	c3->Divide(2,2);
	c3->cd(1);
	h_th_mom->Draw("colz");
	c3->cd(2);
	h_th_cm->Draw("");
	h_th_cm2->SetLineColor(kRed);
	h_th_cm2->Draw("same");
	c3->cd(3);
	h_momL_momR->Draw("colz");
	c3->cd(4);
	h_momL->Draw("");

	TCanvas *c4 = new TCanvas("c4", "c4");
	h_sa_mom_result->Draw("e2");
	h_sa_mom_result->Draw("p same");

	TCanvas *c5 = new TCanvas("c5", "c5");
	h_mom_gen->Draw();
	h_mom_result->SetLineColor(kAzure);
	h_mom_result->Draw("same");

	//label
	ofstream fout("vpflux_SIMC.dat");//LHRS
	//ofstream fout("vpflux_dummy.dat");//RHRS
		fout<<"#VP Flux Calculation with SIMC Acceptance"<<endl;
		fout<<"#1.8<pe[GeV/c]<2.4, 150 partition --> 1bin=4MeV/c"<<endl;
		double vpflux_temp;
		double vpflux_total=0.;
	for(int i=0; i<bin_mom; i++){
		vpflux_temp = h_vp_mom_result2->GetBinContent(i+1);
		double vpflux_tempe = h_vp_mom_result2->GetBinError(i+1);
		h_vp_mom_result3->SetBinContent(i+1,vpflux_temp*0.001);
		h_vp_mom_result3->SetBinError(i+1,vpflux_tempe*0.001);
		double vpflux_temp_ = vpflux_temp;// * 0.004/bin_mom;
		fout<<i+1<<" "<<vpflux_temp_<<endl;
		vpflux_total+=vpflux_temp_;
	}
	
	cout<<"vpflux_total="<<vpflux_total<<endl;
	TCanvas *c6 = new TCanvas("c6", "c6");
	h_vp_mom_result3->Draw("E");
c6->Print("vpflux_full_MeV.pdf");

//c4->Print("LHRS_acceptance.pdf");
//c5->Print("LHRS_mc_gen.pdf");

	if(!BatchFlag){
		theApp.Run();
	}

	if(PDFFlag){
		cout << "Creating a PDF file ... " << endl;
		c1 ->Print(Form("%s[",PDFFile.str().c_str()) );
		c1 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c2 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c3 ->Print(Form("%s" ,PDFFile.str().c_str()) );
		c3 ->Print(Form("%s]" ,PDFFile.str().c_str()) );
	}
//c1->Close();
//c2->Close();
//c3->Close();
	

	//label
	ofstream fout2("LHRS_SIMC.dat");//LHRS
	//ofstream fout2("RHRS_SIMC.dat");//RHRS
		fout2<<"#Acceptance by SIMC"<<endl;
		fout2<<"#1.8<pe[GeV/c]<2.4, 150 partition --> 1bin=4MeV/c"<<endl;
		double Acceptance_temp;
		double Acceptance_total=0.;
		int Acceptance_bin=0;
	for(int i=0; i<bin_mom; i++){
		Acceptance_temp = h_sa_mom_result->GetBinContent(i+1);
		Acceptance_temp = Acceptance_temp;
		fout2<<i+1<<" "<<Acceptance_temp<<endl;
		Acceptance_total+=Acceptance_temp;
		if(Acceptance_temp!=0)Acceptance_bin++;
	}
	cout<<"Acceptance (average)="<<Acceptance_total/(double)Acceptance_bin<<endl;


	cout<< "Finish!" << endl;
	return 0;
}

void SetTitle( TH2D *h, const char *title, const char *xaxis, const char *yaxis){
	h->GetXaxis()->SetTitle(xaxis);
	h->GetYaxis()->SetTitle(yaxis);
	h->SetTitle(title);
}

void SetTitle( TH1D *h, const char *title, const char *xaxis, const char *yaxis){
	h->GetXaxis()->SetTitle(xaxis);
	h->GetYaxis()->SetTitle(yaxis);
	h->SetTitle(title);
}
