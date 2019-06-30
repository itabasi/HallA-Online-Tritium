// Angle at Target resolution //
// Draw macro                 //

double mean_th(int i);
double mean_ph(int i);

const   int nmax=100;

void vdc_ang_resolution(){

  
  bool rarm=false;
  string arm;
  if(rarm)arm="R";
  else arm="L";

  TFile*f2=new TFile(Form("../../rootfiles/angcalib/ang_%sHRS_721.root",arm.c_str()));
  //  TFile*f2=new TFile(Form("../../rootfiles/angcalib/ang_%sHRS_sieve.root",arm.c_str()));
  TTree* t2=(TTree*)f2->Get("T");
  
  double runnum;
  double hallap;
  Double_t trig5;
  Double_t trig4;
  Double_t trig1; 
  int Nrvz,Nlvz;
  int NRvz,NLvz;
  double a1, a2;
  double r_s2_t_pads[nmax];
  double l_s2_t_pads[nmax];
  double r_s2_nthit;
  double l_s2_nthit;
  double r_th_fp[nmax];
  double l_th_fp[nmax];
  double r_ph_fp[nmax];
  double l_ph_fp[nmax];
  double l_x_fp[nmax];
  double r_x_fp[nmax];
  double l_y_fp[nmax];
  double r_y_fp[nmax];
  double ztR_opt[nmax];
  double lx[nmax];
  double ly[nmax];
  double rx[nmax];
  double ry[nmax];
  double ps_asum,sh_asum;
  double gs_asum;
  double Zt[nmax];
  double Rvz[nmax];
  double Lvz[nmax];
  double Rth[nmax], Rph[nmax];
  double Lth[nmax], Lph[nmax];  
    double Xpt[nmax], Ypt[nmax];
  double Xpt_tuned[nmax], Ypt_tuned[nmax];
  double Xpt_init,Ypt_init;
  double Xt, Yt;
  double Zt_tuned[nmax];
  double XFP[nmax],  YFP[nmax];
  double XpFP[nmax], YpFP[nmax];


  
  //======= Get Branch ===============//
  t2->SetBranchAddress("fEvtHdr.fRun", &runnum);
  t2->SetBranchAddress("HALLA_p", &hallap);
  t2->SetBranchAddress("DR.T1", &trig1);
  t2->SetBranchAddress("DR.T4", &trig4);
  t2->SetBranchAddress("R.ps.asum_c", &ps_asum);
  t2->SetBranchAddress("R.sh.asum_c", &sh_asum);  
  t2->SetBranchAddress("DR.T5", &trig5);
  t2->SetBranchAddress("R.a1.asum_c", &a1);
  t2->SetBranchAddress("R.a2.asum_c", &a2);  
  t2->SetBranchAddress("Ndata.R.tr.vz", &NRvz);
  t2->SetBranchAddress("Ndata.L.tr.vz", &NLvz);
  t2->SetBranchAddress("L.cer.asum_c",&gs_asum);  
  t2->SetBranchAddress("R.tr.x",   &r_x_fp);
  t2->SetBranchAddress("L.tr.x",   &l_x_fp);
  t2->SetBranchAddress("R.tr.y",   &r_y_fp);
  t2->SetBranchAddress("L.tr.y",   &l_y_fp);
  t2->SetBranchAddress("R.tr.th",  &r_th_fp);
  t2->SetBranchAddress("L.tr.th",  &l_th_fp);
  t2->SetBranchAddress("R.tr.ph",  &r_ph_fp);
  t2->SetBranchAddress("L.tr.ph",  &l_ph_fp);

  //------ At Target ------------//
  t2->SetBranchAddress("R.tr.vx",   &rx);
  t2->SetBranchAddress("L.tr.vx",   &lx);
  t2->SetBranchAddress("R.tr.vy",   &ry);
  t2->SetBranchAddress("L.tr.vy",   &ly);
  t2->SetBranchAddress("R.tr.vz", &Rvz);
  t2->SetBranchAddress("L.tr.vz", &Lvz);
  t2->SetBranchAddress("R.tr.tg_th", &Rth);
  t2->SetBranchAddress("R.tr.tg_ph", &Rph);
  t2->SetBranchAddress("L.tr.tg_th", &Lth);
  t2->SetBranchAddress("L.tr.tg_ph", &Lph);
  
  if(rarm==true){
   t2->SetBranchAddress("R.tr.tg_th_tuned",Xpt);
   t2->SetBranchAddress("R.tr.tg_ph_tuned",Ypt);   
   t2->SetBranchAddress("R.tr.vz_tuned",Zt);
  } else{
   t2->SetBranchAddress("L.tr.tg_th_tuned",Xpt);
   t2->SetBranchAddress("L.tr.tg_ph_tuned",Ypt);   
   t2->SetBranchAddress("L.tr.vz_tuned",Zt);
  }
   

  //====== Make Hist =====//
  double min_ph=-0.06;
  double max_ph=0.06;
  int bin_ph=300;
  double min_th=-0.1;
  double max_th=0.1;
  int bin_th=300;


  
  TH2F* hLth_ph=new TH2F("hLth_ph","LHRS Theta vs Phi Sieve Slit 2D hist",bin_th,min_th,max_th,bin_ph,min_ph,max_ph);
  hLth_ph->SetTitle("LHRS Theta vs Phi Sieve Slit 2D hist; theta [rad] ;phi [rad]");
  hLth_ph->SetTitleSize(0.06,"x");
  hLth_ph->SetTitleSize(0.06,"y");
  hLth_ph->GetXaxis()->SetTitleOffset(0.8);  
  hLth_ph->GetYaxis()->SetTitleOffset(0.8);    
  TH1F* hLth_c =new TH1F("hLth_c","LHRS theta hist w/ cut",bin_th,min_th,max_th);
  hLth_c->SetTitle("LHRS theta hist w/ cut; theta [rad] ;Counts");
  hLth_c->SetTitleSize(0.06,"x");
  hLth_c->SetTitleSize(0.06,"y");
  hLth_c->GetXaxis()->SetTitleOffset(0.8);  
  hLth_c->GetYaxis()->SetTitleOffset(0.8);      
  TH1F* hLph_c =new TH1F("hLph_c","LHRS phi hist w/ cut",bin_ph,min_ph,max_ph);
  hLph_c->SetTitle("LHRS phi hist w/ cut; phi [rad] ;Counts");  
  hLph_c->SetTitleSize(0.06,"x");
  hLph_c->SetTitleSize(0.06,"y");
  hLph_c->GetXaxis()->SetTitleOffset(0.8);  
  hLph_c->GetYaxis()->SetTitleOffset(0.8);        
  int ENum=t2->GetEntries();
  cout<<"Entries : "<<ENum<<endl;
  bool test=false;
  if(test)ENum=100;
  //===== cut flag ======//
  bool z_cut=false;
  bool th_cut=false;
  bool ph_cut=false;
  bool foil_cut=false;
  bool pid_cut=false;

  //====== Initialization ====//

  for(int j=0;j<nmax;j++){
    Xpt[j]=-2222.0;
    Ypt[j]=-2222.0;    
    Zt[j]=-2222.0;
  }
  gs_asum=-2222.0;
  
  for(int i=0;i<ENum;i++){
  z_cut=false;
  th_cut=false;
  ph_cut=false;
  foil_cut=false;
  pid_cut=false;

    t2->GetEntry(i);

    //===== LHRS Cut definition =====//
    if(rarm == 0){
    if(gs_asum>2000)pid_cut=true;
    if(Zt[0]<0.01 && -0.015< Zt[0])z_cut=true;
    if(-0.005<Xpt[0] && Xpt[0]<0.005)th_cut=true;
    if(-0.005<Ypt[0] && Ypt[0]<0.005)ph_cut=true;    
      }

    //====== FIll hist ===========//
    
    //---- Theta vs phi hist ----//
    if(pid_cut && z_cut && rarm==0)hLth_ph->Fill(Xpt[0],Ypt[0]);
    //---- theta hist ----------//
    if(pid_cut && z_cut && rarm==0 && ph_cut)hLth_c->Fill(Xpt[0]);        
    //---- phi hist ------------//
    if(pid_cut && z_cut && rarm==0 && th_cut)hLph_c->Fill(Ypt[0]);

    if(i%(ENum/10)==0)cout<<"i : "<<i<<" / "<<ENum<<endl;
    
  }//end for


  


  //===== Fit hist ======//
  int nth=11;
  int nph=6;
  TF1* fth[nth];
  TF1* fph[nph];
  TF1* fth_all=new TF1("fth_all","gausn(0)+gausn(3)+gausn(6)+gausn(9)+gausn(12)+gausn(15)+gausn(18)+gausn(21)+gausn(24)+gausn(27)+gausn(30)",min_th,max_th);
  fth_all->SetNpx(2000);
  fth_all->SetLineColor(2);  
  TF1* fph_all=new TF1("fph_all","gausn(0)+gausn(3)+gausn(6)+gausn(9)+gausn(12)+gausn(15)",min_ph,max_ph);
  fph_all->SetNpx(2000);
  fph_all->SetLineColor(2);

  
  for(int i=0;i<nth;i++){
    fth[i]=new TF1(Form("fth_%d",i),"gausn(0)",min_th,max_th);
    fth[i]->SetNpx(2000);
    fth[i]->SetLineColor(i+3);
    fth[i]->SetFillColor(i+3);
    fth[i]->SetFillStyle(3001);
  }
    fth[7]->SetLineColor(46);
    fth[7]->SetFillColor(46);

    
  for(int i=0;i<nph;i++){
    fph[i]=new TF1(Form("fph_%d",i),"gausn(0)",min_ph,max_ph);
    fph[i]->SetNpx(2000);
    fph[i]->SetLineColor(i+3);
    fph[i]->SetFillColor(i+3);
    fph[i]->SetFillStyle(3001);      
  }

  
  double Mean_th;
  double sig_th=0.006;
  double th0[nth],th1[nth],th2[nth];
  for(int i=0;i<nth;i++){
    Mean_th=mean_th(i);
    fth[i]->SetParameter(2,Mean_th);
    fth[i]->SetParameter(3,sig_th);
    if(rarm==0){
      hLth_c->Fit(Form("fth_%d",i),"QR","",Mean_th-sig_th,Mean_th+sig_th);
      th0[i]=fth[i]->GetParameter(0);
      th1[i]=fth[i]->GetParameter(1);
      th2[i]=fth[i]->GetParameter(2);
      fth_all->SetParameter(3*i,   th0[i]);
      fth_all->SetParameter(3*i+1, th1[i]);
      fth_all->SetParameter(3*i+2, th2[i]);
    }
  }

  hLth_c->Fit("fth_all","RQ","",min_th,max_th);
  for(int i=0;i<nth;i++){
    th0[i]=fth_all->GetParameter(3*i  );
    th1[i]=fth_all->GetParameter(3*i+1);
    th2[i]=fth_all->GetParameter(3*i+2);
    fth[i]->SetParameter(3*i,   th0[i]);
    fth[i]->SetParameter(3*i+1, th1[i]);
    fth[i]->SetParameter(3*i+2, th2[i]);
  }
  
  double Mean_ph;
  double sig_ph=0.003;
  double ph0[nph],ph1[nph],ph2[nph];
  for(int i=0;i<nph;i++){
    Mean_ph=mean_ph(i);
    fph[i]->SetParameter(2,Mean_ph);
    fph[i]->SetParameter(3,sig_ph);
    if(rarm==0){
      hLph_c->Fit(Form("fph_%d",i),"QR","",Mean_ph-sig_ph,Mean_ph+sig_ph);
      ph0[i]=fph[i]->GetParameter(0);
      ph1[i]=fph[i]->GetParameter(1);
      ph2[i]=fph[i]->GetParameter(2);
      fph_all->SetParameter(3*i,   ph0[i]);
      fph_all->SetParameter(3*i+1, ph1[i]);
      fph_all->SetParameter(3*i+2, ph2[i]);    }

  }

  hLph_c->Fit("fph_all","RQ","",min_ph,max_ph);
  for(int i=0;i<nph;i++){
    ph0[i]=fph_all->GetParameter(3*i  );
    ph1[i]=fph_all->GetParameter(3*i+1);
    ph2[i]=fph_all->GetParameter(3*i+2);
    fph[i]->SetParameter(3*i,   ph0[i]);
    fph[i]->SetParameter(3*i+1, ph1[i]);
    fph[i]->SetParameter(3*i+2, ph2[i]);
  }
  
  

  TCanvas* c0=new TCanvas("c0","c0");
  TCanvas* cth=new TCanvas("cth","cth");
  TCanvas* cph=new TCanvas("cph","cph");
  c0->cd();
  hLth_ph->Draw("colz");
  cth->cd();
  cth->SetGridx();
  
  hLth_c->Draw();
  for(int i=0;i<nth;i++)fth[i]->Draw("same");
  fth_all->Draw("same");
  
  cph->cd();
  cph->SetGridx();
  hLph_c->Draw();
  for(int i=0;i<nph;i++)fph[i]->Draw("same");
  fph_all->Draw("same");

  
  //============= COMMENT OUT ===============//
  cout<<"================================================"<<endl;
  cout<<"================================================"<<endl;
  cout<<"================================================"<<endl;
  for(int i=0;i<nth;i++){
    cout<<Form("======== theta peak[%d] =========",i)<<endl;
    cout<<"th0 : "<<th0[i]<<" th1: "<<th1[i]<<" th2: "<<th2[i]<<endl;
  }
    cout<<"================================================"<<endl;
    cout<<"================================================"<<endl;
    cout<<"================================================"<<endl;
  for(int i=0;i<nph;i++){
    cout<<Form("======== phi peak[%d] =========",i)<<endl;
    cout<<"ph0 : "<<ph0[i]<<" ph1: "<<ph1[i]<<" ph2: "<<ph2[i]<<endl;
  }
  
  
}//End main


//===================================================//
//============ Function =============================//
//===================================================//

double mean_th(int i){

  double mean;

  double mean_th[11]={-6.29076e-02, -4.97639e-02, -3.74836e-02,-2.47726e-02, -1.22645e-02
		      ,1.68988e-04, 1.23193e-02, 2.47212e-02, 3.65164e-02,  4.92898e-02, 6.34132e-02};

  mean=mean_th[i];
  return mean;
}

double mean_ph(int i){

  double mean;

  double mean_ph[11]={-3.38453e-02, -2.43430e-02, -1.16963e-02,  5.38166e-04,  1.28959e-02, 2.52151e-02};
  mean=mean_ph[i];
  return mean;
}
