
const int nmax=5;
  // nnL peak significance //
  double fit_range=1.0 ; //+- 1 Mev
  int imax=300;

void Alstudy(){

  string ifname[nmax];
  string file;
  TFile* ifp[nmax];
  TH1D* hmm_c[nmax];
  TH1D* hmm_Al_c[nmax];
  TH1D* hmm_nnL_c[nmax];
  TH1D* hacc_nnL_c[nmax];
  TH1D* hpeak_nnL_c[nmax];
  TH1D* hmm_Al_zp[nmax];
  TH1D* hmm_Al_zn[nmax];
  TH1D* hmm_Al_zp_c[nmax];
  TH1D* hmm_Al_zn_c[nmax];  
  TH2D* hmm_Al_nnL[nmax];
  TH2D* hmm_Al_nnL_zp[nmax];
  TH2D* hmm_Al_nnL_zn[nmax];
  //  string dir ="../rootfiles/momcalib/wAl_p18_w2/";
  //  string dir ="../rootfiles/momcalib/wAl_n17/";
  string dir ="../rootfiles/momcalib/Altuning/";
  for(int w=0;w<5;w++){
    //    file = Form("momcalib_5th_wAl_p18_w%d.root",w*10);
    //  file = Form("momcalib_5th_wAl_n17_w%d.root",w*10);
    if(w!=0)    file = Form("momcalib_5th_wAl_w1_m-%d_r2.root",w);
    else if(w==0)file = Form("momcalib_5th_wAl_w1_m%d_r2.root",w);
    ifname[w] =  dir+file;
    ifp[w] =new TFile(ifname[w].c_str());
    hmm_c[w]        = (TH1D*)ifp[w]->Get("hmm_cut");
    hmm_nnL_c[w]    = (TH1D*)ifp[w]->Get("hmm_nnL_cut");
    hacc_nnL_c[w]   = (TH1D*)ifp[w]->Get("hmm_nnL_acc");
    hpeak_nnL_c[w]  = (TH1D*)ifp[w]->Get("hmm_nnL_cut")->Clone();
    hpeak_nnL_c[w]->Add(hacc_nnL_c[w],-1);
    hmm_Al_c[w]    =(TH1D*)ifp[w]->Get("hmm_Al_cut");
    hmm_Al_zp[w]   =(TH1D*)ifp[w]->Get("hmm_Al_zp");
    hmm_Al_zp_c[w] =(TH1D*)ifp[w]->Get("hmm_Al_zp_c");
    hmm_Al_zn[w]   =(TH1D*)ifp[w]->Get("hmm_Al_zn");
    hmm_Al_zn_c[w] =(TH1D*)ifp[w]->Get("hmm_Al_zn_c");    
    hmm_c[w]->SetName(Form("hmm_c_%d",w));
    hmm_Al_nnL[w]   =(TH2D*)ifp[w]->Get("hmm_Al_nnL");
    hmm_Al_nnL_zp[w]   =(TH2D*)ifp[w]->Get("hmm_Al_nnL_zp");
    hmm_Al_nnL_zn[w]   =(TH2D*)ifp[w]->Get("hmm_Al_nnL_zn");
    hmm_nnL_c[w]->SetName(Form("hmm_nnL_c_%d",w));
    hacc_nnL_c[w]->SetName(Form("hacc_nnL_c_%d",w));
    hpeak_nnL_c[w]->SetName(Form("hpeak_nnL_c_%d",w));
    hmm_Al_c[w]->SetName(Form("hmm_Al_c_%d",w));
    hmm_Al_zp[w]->SetName(Form("hmm_Al_zp_%d",w));
    hmm_Al_zp_c[w]->SetName(Form("hmm_Al_zp_c_%d",w));
    hmm_Al_zn[w]->SetName(Form("hmm_Al_zn_%d",w));
    hmm_Al_zn_c[w]->SetName(Form("hmm_Al_zn_c_%d",w));
    hmm_Al_nnL[w]->SetName(Form("hmm_Al_nnL_%d",w));
    hmm_Al_nnL_zp[w]->SetName(Form("hmm_Al_nnL_zp_%d",w));
    hmm_Al_nnL_zn[w]->SetName(Form("hmm_Al_nnL_zn_%d",w));
    hmm_Al_zp[w]->SetLineColor(6);
    hmm_Al_zp_c[w]->SetLineColor(2);
    hmm_Al_zn[w]->SetLineColor(9);
    hmm_Al_zn_c[w]->SetLineColor(4);    
    
  } 




  ///===== Analysis ========//
  TF1* fL[nmax];
  TF1* fS[nmax];
  TGraphErrors* gL_mean=new TGraphErrors();
  TGraphErrors* gS_mean=new TGraphErrors();
  TGraphErrors* gL_sig=new TGraphErrors();
  TGraphErrors* gS_sig=new TGraphErrors();

  gL_mean->SetName("gL_mean");
  gL_mean->SetMarkerColor(2);
  gL_mean->SetMarkerStyle(20);
  gS_mean->SetName("gS_mean");
  gS_mean->SetMarkerColor(4);
  gS_mean->SetMarkerStyle(20);
  gL_sig->SetName("gL_sig");
  gL_sig->SetMarkerColor(2);
  gL_sig->SetMarkerStyle(20);
  gS_sig->SetName("gS_sig");
  gS_sig->SetMarkerColor(4);
  gS_sig->SetMarkerStyle(20);  
  double mean_sigma=76.959;  //MeV
  double Lam_mean[nmax],Lam_sig[nmax],Sig_mean[nmax],Sig_sig[nmax];
  for(int w=0;w<nmax;w++){
    fL[w]=new TF1(Form("fL_%d",w),"gausn(0)",-50,50);
    fS[w]=new TF1(Form("fS_%d",w),"gausn(0)",50,100);

    hmm_c[w]->Fit(Form("fL_%d",w),"QR","QR",-5,5);
    Lam_mean[w] = fL[w]->GetParameter(1);
    Lam_sig[w]  = fL[w]->GetParameter(2);
    hmm_c[w]->Fit(Form("fS_%d",w),"QR","QR",-5+mean_sigma,5+mean_sigma);
    Sig_mean[w] = fS[w]->GetParameter(1) - mean_sigma;
    Sig_sig[w]  = fS[w]->GetParameter(2);

    gL_mean->SetPoint(w,w*10.,Lam_mean[w]);
    gL_sig ->SetPoint(w,w*10.,Lam_sig[w]);
    gS_mean->SetPoint(w,w*10.,Sig_mean[w]);
    gS_sig ->SetPoint(w,w*10.,Sig_sig[w]);
    
 
  }


  
  TF1* fnnL[nmax];
  TF1* fbg[nmax];
  //  TF1
  double Ntot[nmax][imax],Npeak[nmax][imax],Nbg[nmax][imax];
  double sig_nnL[nmax][imax], mean_nnL[nmax][imax];
  TGraphErrors* gnnL[nmax];
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
  double PS[nmax][imax];
  double min_fit,max_fit;
  double max_PS[nmax];
  double max_mean[nmax];
  int n=0;
  for(int w=0;w<nmax;w++){
    fnnL[w]=new TF1(Form("fnnL_%d",w),"gausn(0)",-300,300);
    gnnL[w]=new TGraphErrors();
    gnnL[w]->SetName(Form("gnnL_%d",w));
    gnnL[w]->SetMarkerStyle(20);
    gnnL[w]->SetMarkerColor(w+1);
    
    n=0;
    for(int i=1;i<imax;i++){
      
      min_fit = -300. + (double)(2*i)* fit_range;
      max_fit = -300. + (double)(2*i +2)* fit_range;
      double min_fit_l = -300. + (double)(2*(i-1))* fit_range;
      double max_fit_l = -300. + (double)(2*(i-1) +2)* fit_range;
      double min_fit_h = -300. + (double)(2*(i+1))* fit_range;
      double max_fit_h = -300. + (double)(2*(i+1) +2)* fit_range;      

      if(w==0){
	gPS[i]=new TGraphErrors();
	gPS[i]->SetName(Form("gPS_%d",i));
	gPS[i]->SetTitle(Form("Peak Significance at %f; weight ; PS ",(min_fit+max_fit)/2.0));
	gPS[i]->SetMarkerStyle(20);  
      }

      //      cout<<"w "<<w<<" i "<<i<<" min_fit"<<min_fit<<" max_fit "<<max_fit<<endl;
      //      hpeak_nnL->Fit(Form("fnnL_%d",w),min_fit,max_fit);
      //      Npeak[w][i]    = fnnL[w]->Integral(min_fit,max_fit);
      //      mean_nnL[w][i] = fnnL[w] ->GetParameter(1);
      //      sig_nnL[w][i]  = fnnL[w] ->GetParameter(2);
      
      int bin_min_nnL,bin_max_nnL;      
      bin_min_nnL    = hmm_nnL_c[w]->GetXaxis()->FindBin(min_fit);
      bin_max_nnL    = hmm_nnL_c[w]->GetXaxis()->FindBin(max_fit);
      Ntot[w][i]     = hmm_nnL_c[w]->Integral(bin_min_nnL,bin_max_nnL);

      
      int bin_min_nnL_acc_l,bin_max_nnL_acc_l;      
      bin_min_nnL_acc_l    = hmm_nnL_c[w]->GetXaxis()->FindBin(min_fit_l);
      bin_max_nnL_acc_l    = hmm_nnL_c[w]->GetXaxis()->FindBin(max_fit_l);

      int bin_min_nnL_acc_h,bin_max_nnL_acc_h;      
      bin_min_nnL_acc_h    = hmm_nnL_c[w]->GetXaxis()->FindBin(min_fit_h);
      bin_max_nnL_acc_h    = hmm_nnL_c[w]->GetXaxis()->FindBin(max_fit_h);   
      Nbg[w][i]      = (hmm_nnL_c[w]->Integral(bin_min_nnL_acc_l,bin_max_nnL_acc_l) + hmm_nnL_c[w]->Integral(bin_min_nnL_acc_h,bin_max_nnL_acc_h) )/2.0;

      //      cout<<"w "<<w<<" i "<<i<<" bin_min_acc "<<bin_min_nnL_acc<<" bin_max_acc "<<bin_max_nnL_acc<<" nbg "<<Nbg[w][i]<<endl;
      
      Npeak[w][i] = Ntot[w][i]-Nbg[w][i];

      if(Ntot[w][i]>0 && Npeak[w][i]>0)PS[w][i]    = Npeak[w][i]/sqrt(Ntot[w][i]);
      else PS[w][i]=0.0;
      gnnL[w]->SetPoint(n,(max_fit+min_fit)/2.,PS[w][i]);
      //      cout<<"w "<<w<<"i "<<i<<"Npeak "<<Npeak[w][i]<<" Nbg "<<Nbg[w][i]<<" PS "<<PS[w][i]<<" fit_min "<<min_fit<<" fit_max "<<max_fit<<endl;
      gPS[i]->SetPoint(w,w,PS[w][i]);
      
      if(max_PS[w]< PS[w][i]){
	max_PS[w]=PS[w][i];
	max_mean[w] = (max_fit+min_fit)/2.0;
      }

      if(fabs((max_fit+min_fit)/2.-3.0 )<fit_range)
	gPS_p3->SetPoint(w,w*10.,PS[w][i]);
      if(fabs((max_fit+min_fit)/2.+3.0 )<fit_range)
	gPS_n3->SetPoint(w,w*10.,PS[w][i]);      
      
      
      n++;  
    }

    gPS_max->SetPoint(w,w*10.,max_PS[w]);
    gmean->SetPoint(w,w*10.,max_mean[w]);
    
  }

  
  //=====================================================//
  
  string ofname=dir +"hist.root";
  TFile* ofp=new TFile(ofname.c_str(),"recreate");  


  for(int i=1;i<imax;i++)gPS[i]->Write();
    for(int w=0;w<nmax;w++){
    hmm_c[w]->Write();
    hmm_nnL_c[w]->Write();
    hacc_nnL_c[w]->Write();
    hpeak_nnL_c[w]->Write();
    hmm_Al_c[w]->Write();
    hmm_Al_zp[w]->Write();
    hmm_Al_zn[w]->Write();
    hmm_Al_zp_c[w]->Write();
    hmm_Al_zn_c[w]->Write();
    fL[w]->Write();
    fS[w]->Write();
    hmm_Al_nnL[w]->GetXaxis()->SetRangeUser(-50,50);
    hmm_Al_nnL_zn[w]->GetXaxis()->SetRangeUser(-50,50);
    hmm_Al_nnL_zp[w]->GetXaxis()->SetRangeUser(-50,50);
    hmm_Al_nnL[w]->GetYaxis()->SetRangeUser(-50,50);
    hmm_Al_nnL_zn[w]->GetYaxis()->SetRangeUser(-50,50);
    hmm_Al_nnL_zp[w]->GetYaxis()->SetRangeUser(-50,50);    
    hmm_Al_nnL[w]->Write();
    hmm_Al_nnL_zn[w]->Write();
    hmm_Al_nnL_zp[w]->Write();
    gnnL[w]->Write();
    }
    

    
    gL_mean->Write();
    gS_mean->Write();
    gL_sig->Write();
    gS_sig->Write();
    gPS_max->Write();
    gPS_n3->Write();
    gPS_p3->Write();
    gmean->Write();
    cout<<"output root "<<ofname<<endl;
    gSystem->Exit(1);
}
