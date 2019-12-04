////////////////////////////////////////
//  Coin time Run dependece correction
//  Auther Itabashi Nov. 25th, 2019
////////////////////////////////////////

void coin_run(){

  string ifname="../rootfiles/mmass/ana_Lambda/Lambda_small_optH_1119.root";
  TFile* ifs=new TFile(ifname.c_str());
  TChain* T=new TChain("T");
  T->Add(ifname.c_str());

  int ENum=T->GetEntries();
  cout<<"Events : "<<ENum<<endl;
  TCanvas* c0=new TCanvas("c0","c0");

  double ct,Rx_fp,ac1,ac2,trig;
  int z_cut,runnum;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("ct_b",1);
  T->SetBranchAddress("ct_b",&ct);
  T->SetBranchStatus("runnum",1);
  T->SetBranchAddress("runnum",&runnum);  
  T->SetBranchStatus("Rx_fp",1);
  T->SetBranchAddress("Rx_fp",&Rx_fp);  
  T->SetBranchStatus("Rx_fp",1);
  T->SetBranchAddress("Rx_fp",&Rx_fp);
  T->SetBranchStatus("trig",1);
  T->SetBranchAddress("trig",&trig);
  T->SetBranchStatus("z_cut",1);
  T->SetBranchAddress("z_cut",&z_cut);    
  T->SetBranchStatus("ac1_npe_sum",1);
  T->SetBranchAddress("ac1_npe_sum",&ac1);
  T->SetBranchStatus("ac2_npe_sum",1);
  T->SetBranchAddress("ac2_npe_sum",&ac2);


  int bin_run,min_run,max_run;

  double min_ct=-5.0;
  double max_ct=5.0;

  int bin_ct=(int)(max_ct-min_ct)*20;
  min_run=111160;
  max_run=111220;
  bin_run=max_run-min_run;

  


  
  TH2D*hcoin_run=new TH2D("hcoin_run","Coin time Run dependece ; Runnum ; coin time [ns]",bin_run,min_run,max_run,bin_ct,min_ct,max_ct);


  bool pi_cut=false;
  for(int k=0; k<ENum;k++){

    ct=-1000.;
    Rx_fp=-10.;
    z_cut=-1;
    trig=0;
    pi_cut=false;
    T->GetEntry(k);
    if(ac1>1.0 && ac2>5.0)pi_cut=true;
    if(pi_cut && z_cut>0 && trig==5)hcoin_run->Fill(runnum,ct);
  }


  TH1D* hcoin_s[bin_run];
  TF1* fpi[bin_run];
  double mean[bin_run],mean_err[bin_run];
  double mean_pi=3.0;
  TGraphErrors* gcoin=new TGraphErrors();
  int nrun[bin_run];
  for(int i=0;i<bin_run;i++){
    
    hcoin_s[i]=hcoin_run->ProjectionY(Form("hcoin_s_%d",i),i,i+1);
    fpi[i]=new TF1(Form("fpi_%d",i),"gausn(0)",min_ct,max_ct);
    fpi[i]->SetLineColor(2);

    hcoin_s[i]->Fit(Form("fpi_%d",i),"RQ","RQ",mean_pi-0.2,mean_pi+0.2);
    mean[i]=fpi[i]->GetParameter(1);
    mean_err[i]=fpi[i]->GetParError(1);
    gcoin->SetPoint(i,i,mean[i]);
    gcoin->SetPointError(i,0,mean_err[i]);
    nrun[i]=min_run + i;
  }  


  string ofname="./test.param";
  ofstream ofs(ofname.c_str());

  ofs << " #run  pi-offset "<< endl; 
  for(int i=0;i<bin_run;i++){

    ofs << min_run + i << " " << mean[i] -3.0 << endl;
    

  }

  ofs.close();

  c0->cd();
  gcoin->SetMarkerStyle(20);
  gcoin->SetMarkerColor(4);
  gcoin->Draw("AP");
  
}
