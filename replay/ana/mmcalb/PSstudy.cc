
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
  

    ifname =  dir+file;
    ifp =new TFile(ifname.c_str());

    hmm_nnL    = (TH1D*)ifp->Get("h_mm_nnL");
    hacc_nnL   = (TH1D*)ifp->Get("h_acc_nnL");
    hpeak_nnL  = (TH1D*)ifp->Get("h_peak_nnL")->Clone();
    
    hmm_nnL->SetName("hmm_nnL");
    hacc_nnL->SetName("hacc_nnL");
    hpeak_nnL->SetName("hpeak_nnL");

 
  
  TF1* fnnL;
  TF1* fbg;
  //  TF1
  double Ntot[imax],Npeak[imax],Nbg[imax];
  double sig_nnL[imax], mean_nnL[imax];
  TGraphErrors* gnnL;
  TGraphErrors* gPS_max =new TGraphErrors();
  gPS_max->SetName("gPS_max");
  gPS_max->SetMarkerStyle(20);

  TGraphErrors* gPS_n3 =new TGraphErrors();
  gPS_n3->SetName("gPS_n3");
  gPS_n3->SetMarkerStyle(20);
  gPS_n3->SetMarkerColor(4);
  gPS_n3->SetTitle("BL =-3 MeV PS; weight ; Peak Significane");
  TGraphErrors* gPS_p3 =new TGraphErrors();
  gPS_p3->SetName("gPS_p3");
  gPS_p3->SetTitle("BL =+3 MeV PS; weight ; Peak Significane");
  gPS_p3->SetMarkerStyle(20);
  gPS_p3->SetMarkerColor(2);


  TGraphErrors* gPS[imax];

  TGraphErrors* gmean =new TGraphErrors();
  gmean->SetName("gmean");
  gmean->SetMarkerStyle(20);  
  double PS[imax];
  double min_fit,max_fit;
  double max_PS;
  double max_mean;
  int n=0;

    fnnL=new TF1("fnnL","gausn(0)",-300,300);
    gnnL=new TGraphErrors();
    gnnL->SetName("gnnL");
    gnnL->SetMarkerStyle(20);
    gnnL->SetMarkerColor(1);
    
    n=0;
    
    for(int i=1;i<imax;i++){
      
      min_fit = -300. + (double)(2*i)* fit_range;
      max_fit = -300. + (double)(2*i +2)* fit_range;
      double min_fit_l = -300. + (double)(2*(i-1))* fit_range;
      double max_fit_l = -300. + (double)(2*(i-1) +2)* fit_range;
      double min_fit_h = -300. + (double)(2*(i+1))* fit_range;
      double max_fit_h = -300. + (double)(2*(i+1) +2)* fit_range;      


	gPS[i]=new TGraphErrors();
	gPS[i]->SetName(Form("gPS_%d",i));
	gPS[i]->SetTitle(Form("Peak Significance at %f; weight ; PS ",(min_fit+max_fit)/2.0));
	gPS[i]->SetMarkerStyle(20);  

      
      int bin_min_nnL,bin_max_nnL;      
      bin_min_nnL    = hmm_nnL->GetXaxis()->FindBin(min_fit);
      bin_max_nnL    = hmm_nnL->GetXaxis()->FindBin(max_fit);
      Ntot[i]     = hmm_nnL->Integral(bin_min_nnL,bin_max_nnL);

      
      int bin_min_nnL_acc_l,bin_max_nnL_acc_l;      
      bin_min_nnL_acc_l    = hmm_nnL->GetXaxis()->FindBin(min_fit_l);
      bin_max_nnL_acc_l    = hmm_nnL->GetXaxis()->FindBin(max_fit_l);

      int bin_min_nnL_acc_h,bin_max_nnL_acc_h;      
      bin_min_nnL_acc_h    = hmm_nnL->GetXaxis()->FindBin(min_fit_h);
      bin_max_nnL_acc_h    = hmm_nnL->GetXaxis()->FindBin(max_fit_h);   
      Nbg[i]      = (hmm_nnL->Integral(bin_min_nnL_acc_l,bin_max_nnL_acc_l) + hmm_nnL->Integral(bin_min_nnL_acc_h,bin_max_nnL_acc_h) )/2.0;

      
      Npeak[i] = Ntot[i]-Nbg[i];

      if(Ntot[i]>0 && Npeak[i]>0)PS[i]    = Npeak[i]/sqrt(Ntot[i]);
      else PS[i]=0.0;
      gnnL->SetPoint(n,(max_fit+min_fit)/2.,PS[i]);

      
      if(max_PS< PS[i]){
	max_PS=PS[i];
	max_mean = (max_fit+min_fit)/2.0;
      }

      
      n++;  
    }


  
  //=====================================================//
  
  string ofname=dir +"hist.root";
  TFile* ofp=new TFile(ofname.c_str(),"recreate");  


  for(int i=1;i<imax;i++)gPS[i]->Write();


    hmm_nnL->Write();
    hacc_nnL->Write();
    hpeak_nnL->Write();
    gnnL->Write();

    

    

    cout<<"output root "<<ofname<<endl;
    gSystem->Exit(1);
}
