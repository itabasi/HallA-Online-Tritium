// Sieve hole Resolution caluclation macro
// Aurther Itabshi 2019 07 /01 

extern double sqgaus(double *x, double *par);
extern double square(double x, double r);

void sieve_resolution(){

  bool rarm=false;
  TFile*f1 ;
  if(rarm)f1=new TFile(Form("../rootfiles/angcalib/ang_RHRS_woRas.root"));
  else {
    
    //    f1=new TFile(Form("../rootfiles/angcalib/ang_LHRS_woRas.root"));}
    //    f1=new TFile(Form("../rootfiles/angcalib/ang_LHRS_sieve.root"));}
    f1=new TFile(Form("../rootfiles/angcalib/angcalib_LHRS_4th_0915_0.root"));}
    
  TTree* t1=(TTree*)f1->Get("T");


  Double_t trig5;
  Double_t trig4;
  Double_t trig1;
  double ps_asum,sh_asum;
  double gs_asum;
  double Zt[100];
  double Rvz[100];
  double Lvz[100];
  double Zt_tuned[100];
  double ztR_opt[100];
  double a1,a2;
  double ssx,ssy;
  double cer;
  t1->SetBranchAddress("DR.T1", &trig1);
  t1->SetBranchAddress("DR.T4", &trig4);
  t1->SetBranchAddress("R.ps.asum_c", &ps_asum);
  t1->SetBranchAddress("R.sh.asum_c", &sh_asum);
  t1->SetBranchAddress("L.cer.asum_c", &cer);  
  t1->SetBranchAddress("DR.T5", &trig5);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);    
  t1->SetBranchAddress("ss_x", &ssx);
  t1->SetBranchAddress("ss_y", &ssy);      
  if(rarm==true){
   t1->SetBranchAddress("R.tr.vz_opt",ztR_opt);
   t1->SetBranchAddress("R.tr.vz_tuned",Zt);
  } else{
   t1->SetBranchAddress("L.tr.vz_opt",ztR_opt);
   t1->SetBranchAddress("L.tr.vz_tuned",Zt);
  }


  //======= Make Hist ======//
  double min_x=-10.0; //[cm]
  double max_x=10.0;
  double min_y=-10.0; //[cm]
  double max_y=10.0;  
  int bin_x=500;
  int bin_y=500;
  int j=0;
  TH1F* hth=new TH1F("hth","theta",bin_x,min_x,max_x);
  hth->SetTitle("#theta at Target hist; Sieve slit [cm] ; Counts");
  hth->SetTitleSize(0.06,"x");
  hth->SetTitleSize(0.06,"y");
  hth->GetXaxis()->SetTitleOffset(0.8);  
  hth->GetYaxis()->SetTitleOffset(0.8);    
  TH1F* hph=new TH1F("hph","phi",bin_y,min_y,max_y);
  hph->SetTitle("#phi at Target hist; Sieve slit [cm] ; Counts");  
  hph->SetTitleSize(0.06,"x");
  hph->SetTitleSize(0.06,"y");
  hph->GetXaxis()->SetTitleOffset(0.8);  
  hph->GetYaxis()->SetTitleOffset(0.8);    
  TH2F* hth_ph=new TH2F("hth_ph","theta vs phi",bin_x,min_x,max_x,bin_y,min_y,max_y);
  
  
  int ENum=t1->GetEntries();
  //  ENum=10;
  cout<<"Entries : "<<ENum<<endl;
  bool PID;
  bool ph_cut;
  bool th_cut;
  bool z_cut;
  for(int k=0;k<ENum;k++){

    PID=false;
    z_cut=false;
    th_cut=false;
    ph_cut=false;
    
    t1->GetEntry(k);
   if(cer>2000 && rarm==0)PID=true;
    if(Zt[0]<0.01 && -0.015< Zt[0])z_cut=true;
    if(-0.6 < ssx && ssx < 0.6)th_cut=true;
    if(-0.4 < ssy && ssy < 0.4)ph_cut=true;    
   
   if( PID ){
     if(z_cut && ph_cut) hth->Fill(ssx);
     if(z_cut && th_cut) hph->Fill(ssy);
       hth_ph->Fill(ssx,ssy);
   }
   
   if(k%(ENum/10)==0){
     cout<<"FIlled : "<<j*10<<" %"<<endl;
     j++;}
  }// Filled

	

  

  
  
   double hole1_size=0.4;// hole radius [cm]
   double hole2_size=0.2;// hole radius [cm]

  double min, max;

  double ss_x[5]={-5.0, -2.5, 0.0, 2.5, 5.0};
  double ss_y[5]={-2.5, -1.25, 0.0, 1.25, 2.5};

  int  nth = 5;
  int  nph = 5;  
  TF1* fth[nth];
  TF1* fph[nth];
  TGraphErrors* gth = new TGraphErrors();
  gth->SetMarkerColor(2);
  gth->SetMarkerStyle(20);
  gth->SetMarkerSize(1.0);

  gth->SetTitle("#theta resolution [cm]; #peak ; sigma [cm]");
  TGraphErrors* gph = new TGraphErrors();
  //  gth->SetTitleSize(0.06,"x");
  //  gth->SetTitleSize(0.06,"y");
  gth->GetXaxis()->SetTitleOffset(0.8);  
  gth->GetYaxis()->SetTitleOffset(0.8);    
  gph->SetMarkerColor(2);
  gph->SetMarkerStyle(20);
  gph->SetMarkerSize(1.0);
  gph->SetTitle("#phi resolution [cm]; #peak ; sigma [cm]");
  //  gph->SetTitleSize(0.06,"x");
  //  gph->SetTitleSize(0.06,"y");
  gph->GetXaxis()->SetTitleOffset(0.8);  
  gph->GetYaxis()->SetTitleOffset(0.8);
  
  double sig_th[nth] , sig_th_err[nth];
  double sig_ph[nph] , sig_ph_err[nph];  
  for(int i=0; i<nth;i++){

  fth[i] = new TF1(Form("fth_%d",i),"sqgaus",-10.,10.,5);    
  fth[i]->SetNpx(2000);
  fth[i]->SetParameter(0,300);
  //  fth[i]->SetParLimits(0,100,1000);    
  if(i==2)fth[i]->FixParameter(1,hole1_size);
  else fth[i]->FixParameter(1,hole2_size);
  fth[i]->SetParameter(2,0.2);
  fth[i]->SetParLimits(2,0.1,0.5);    
  fth[i]->SetParameter(3,ss_x[i]);
  //  fth[i]->SetParLimits(3,ss_x[i]-hole1_size,ss_x[i]+hole1_size);  
  fth[i]->SetLineColor(i+2);
  fth[i]->SetFillColor(i+2);
  fth[i]->SetFillStyle(3001);
  hth->Fit(Form("fth_%d",i),"","",ss_x[i]-0.3,ss_x[i]+0.3);
  sig_th[i]=fth[i]->GetParameter(2);
  sig_th_err[i]=fth[i]->GetParError(2);  
  gth->SetPoint(i,i,sig_th[i]);
  //  if(sig_th[i]!=hole2_size)
    gth->SetPointError(i,0,sig_th_err[i]);
  }


  cout<<endl;
  cout<<endl;
  cout<<"======== phi fitting ======"<<endl;
  cout<<endl;
  cout<<endl;  

  //====== phi =================//

  for(int i=0; i<nph;i++){
  fph[i] = new TF1(Form("fph_%d",i),"sqgaus",-10.,10.,4);
  fph[i]->SetNpx(2000);
  fph[i]->SetParameter(0,1000);
  if(i==2)fph[i]->FixParameter(1,hole1_size);
  else fph[i]->FixParameter(1,hole2_size);
  fph[i]->SetParameter(2,0.1);
  fph[i]->SetParameter(3,ss_y[i]);
  fph[i]->SetLineColor(i+2);
  fph[i]->SetFillColor(i+2);
  fph[i]->SetFillStyle(3001);
  hph->Fit(Form("fph_%d",i),"","",ss_y[i]-0.3,ss_y[i]+0.3);
  sig_ph[i]=fph[i]->GetParameter(2);
  gph->SetPoint(i,i,sig_ph[i]);
  sig_ph_err[i]=fph[i]->GetParError(2);
  if(sig_ph[i]!=hole2_size)gph->SetPointError(i,0,sig_ph_err[i]);    
  }
 



  //================================//
  //=========== Draw ===============//
  //===============================//
  
  

  
  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  hth->Draw();
  for(int i=0;i<nth;i++)fth[i]->Draw("same");

  TCanvas* c1=new TCanvas("c1","c1");
  c1->cd();
  hph->Draw();
  for(int i=0;i<nph;i++)fph[i]->Draw("same");  


  
  TCanvas* c2=new TCanvas("c2","c2"); 
  c2->cd();
  gth->SetMinimum(0.1);
  gth->SetMaximum(0.4);			 
  gth->Draw("AP");
  

  TCanvas* c3=new TCanvas("c3","c3"); 
  c3->cd();
  gph->SetMinimum(0.0);
  gph->SetMaximum(0.2);			   
  gph->Draw("AP");  
 
  
}//end main


//====== Function ========//

double sqgaus(double *x, double *par) {
  //par[0]=Total area
  //par[1]=hole size of sieve slits
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;
  double sq;
// Range of convolution integral
//  xlow = x[0] - sc * par[2];
  xlow =  - sc * par[2];
  xupp =  + sc * par[2];
  //  if(par[3]==-5.0 || par[3]==5.0)  cout<<"par[3] : "<<par[3]<<" xlow : "<<xlow<<" xupp : "<<xupp<<endl;
  //  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step;
     fland = TMath::Gaus(xx,x[0]-par[3],par[2]);
     sq=square(xx, par[1] ); 
     sum += fland *sq;
     xx = xupp - (i-0.5) * step;
     fland = TMath::Gaus(xx,x[0]-par[3],par[2]);
     sq=square(xx , par[1] ); 
     sum += fland *sq;
  }
     val = par[0] * step * sum * invsq2pi / par[2];
  return val;
  }


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


double square(double x, double r){
  double y;
  if(-r < x  && x < r ) y=1.0;
  else y=0.0;

  return y;

}
