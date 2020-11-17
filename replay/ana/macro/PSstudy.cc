#include <TROOT.h>
#include <TMath.h>
#include "Math/Math.h"
const int nmax=5;
  // nnL peak significance //
  double fit_range=1.0 ; //+- 1 Mev
  int imax=300;

void PSstudy(){

  string ifname;
  string file;
  TFile* ifp;

  TH1D* hmm_nnL;
  TH1D* hacc_nnL;
  TH1D* hpeak_nnL;

  string dir ="../rootfiles/mmass/ana_Lambda/2020-09-02_2/";
  file = "nnL_small_Ole_all.root";
  

  double mm,mm_mixed[10];
  
    ifname =  dir+file;
    ifp =new TFile(ifname.c_str());

    TChain* T=new TChain("T");
    T->Add(ifname.c_str());


    T->SetBranchStatus("*",0);
    T->SetBranchStatus("mm_nnL",1);
    T->SetBranchAddress("mm_nnL",&mm);
    T->SetBranchStatus("mm_mix",1);
    T->SetBranchAddress("mm_mix",mm_mixed);
    




    

    int ENum = T->GetEntries();
    cout<<"ENum "<<ENum<<endl;


    TH1D* hmm =new TH1D("hmm","nnL missing mass ;-B_{#Lambda} [MeV] ; Counts",10000,-150,150);
    TH1D* hacc =new TH1D("hacc","nnL missing mass (Accidental B.G.) ;-B_{#Lambda} [MeV] ; Counts",10000,-150,150);

    TH1D* hmm_p =new TH1D("hmm_p","nnL missing mass Peak;-B_{#Lambda} [MeV] ; Counts",10000,-150,150);

    for(int k=0;k<ENum;k++){

      T->GetEntry(k);

      hmm->Fill(mm);
      hmm_p->Fill(mm);
      for(int i=0;i<10;i++)
	hacc->Fill(mm_mixed[i],1./200.);

    }
    


    //    hmm_p->Add(hacc,-1);



    TCanvas* c0=new TCanvas("c0","c0");
    c0->cd();

    hmm_p->Draw();


    // Peak Significance Study //


    TGraphAsymmErrors*gPS=new TGraphAsymmErrors();
    gPS->SetName("gPS");
    gPS->SetMarkerColor(2);
    gPS->SetMarkerStyle(22);
    //
    double sigma =1.0; // MeV
    double side = 2.0*sigma;
    int nbin = hmm_p->GetXaxis()->GetNbins();
    int min_bin = hmm_p->GetXaxis()->FindBin(-100.);
    int max_bin = hmm_p->GetXaxis()->FindBin(100.);
    double mean;
    int bin_side_min,bin_side_max;
    double ymin,ymax,y;
    double Nbg,Np,Nps;
    double dNbg_u,dNbg_l,dNp_u,dNp_l,dPS_u,dPS_l;
    int i=0;
    for(int ibin=min_bin;ibin<max_bin;ibin++){
      
      mean =    hmm_p->GetBinCenter(ibin);
      y    =    hmm_p->GetBinContent(ibin); 
      bin_side_min = hmm_p->FindBin(mean - side);
      bin_side_max = hmm_p->FindBin(mean + side);

      ymin = hmm_p ->GetBinContent(bin_side_min);
      ymax = hmm_p ->GetBinContent(bin_side_max);

      double a = (ymax - ymin)/(2.0*side);
      double b = ymax - a * (mean + side);

      Nbg  = a*mean +b;
      Np   = y - Nbg;

      Nps =Np/sqrt(Np+Nbg);
	
      dNbg_l =sqrt(Nbg);
      dNbg_u =sqrt(Nbg);
      dNp_l =sqrt(Np);
      dNp_u =sqrt(Np);


      dPS_u = (fabs(1.0 - Np/(2.0*(Np+Nbg)))*dNp_u + Nbg/(2.*(Np+Nbg)*dNbg_u))/sqrt(Nbg +Np);
            dPS_l = (fabs(1.0 - Np/(2.0*(Np+Nbg)))*dNp_l + Nbg/(2.*(Np+Nbg)*dNbg_l))/sqrt(Nbg +Np);
	    //	    cout<<"mean "<<mean<<endl;
	    if(Np+Nbg>1 && Np-Nbg>0){
	      //	      cout<<"Nps "<<Nps<<" Np "<<Np<<" Nbg "<<Nbg<<endl;
	    gPS->SetPoint(i, mean, Nps);
	    gPS->SetPointError(i,0.0,0.0,dPS_l,dPS_u);
	    }
      //      dNbg = TMath::Poisson(1,sigma);
      //      dNbg = ROOT::Math::poisson_cdf(,y);

	    i++;
      
    }


    TCanvas* c1 =new TCanvas("c1","c1");
    c1->cd();
    gPS->Draw("AP");
    
    
    string ofname ="hist.root";
    TFile* ofp =new TFile(ofname.c_str(),"recreate");
    hmm_p->Write();
    
    
}
