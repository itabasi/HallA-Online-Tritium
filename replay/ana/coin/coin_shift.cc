

void coin_shift(){

  string ifname="../rootfiles/mmass/ana_Lambda/Lambda_small_optH_1119.root";

  TFile* ifs=new TFile(ifname.c_str());
  TChain* T=new TChain("T");
  T->Add(ifname.c_str());


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


  int bin_run;
  int min_run,max_run;
  min_run=111160;
  max_run=111220;
  bin_run=(int)(max_run-min_run);
  int bin_ct;
  double min_ct,max_ct;
  min_ct=-5000.;
  max_ct=5000;
  bin_ct=(int)(max_ct-min_ct)*50;
  
  TH2D* hcoin=new TH2D("hcoin","Coin-time rundependece ; runnum ; coin time [ns]",bin_run,min_run,max_run,bin_ct,min_ct,max_ct);

  int ENum=T->GetEntries();
  cout<<"Events : "<<ENum<<endl;

  for(int k=0;k<ENum;k++){
    ct=-1000.;
    Rx_fp=-10.;
    z_cut=-1;
    trig=0;
    
    T->GetEntry(k);
    if(trig==5)hcoin->Fill(runnum,ct);
    }

  //=====================================//
  //=========== Project =================//
  //=====================================//

  TCanvas* c0=new TCanvas("c0","c0");
  TCanvas* c1=new TCanvas("c1","c1");
  TCanvas* c2=new TCanvas("c2","c2");
  TCanvas* c3=new TCanvas("c3","c3");
  TCanvas* c4=new TCanvas("c4","c4");
  
  TH1D* hcoin_s[bin_run];
  TF1* fpi[bin_run];
  TF1* fpi2[bin_run];
  double mean_pi=3.0;
  double mean_pi2= 3.66864e+03;
  double mean[bin_run],mean2[bin_run],mean_err[bin_run],mean2_err[bin_run];
  TGraphErrors* gcoin=new TGraphErrors();
  TGraphErrors* gcoin2=new TGraphErrors();
  TGraphErrors* gcoin3=new TGraphErrors();
  gcoin->SetName("gcoin");
  gcoin2->SetName("gcoin2");
  gcoin2->SetName("gcoin3");
  gcoin->SetMarkerColor(2);
  gcoin->SetMarkerStyle(20);
  gcoin->SetFillStyle(3005);  
  gcoin->SetFillColor(2);
  gcoin2->SetMarkerColor(3);
  gcoin2->SetMarkerStyle(20);
  gcoin2->SetFillStyle(3005);  
  gcoin2->SetFillColor(3);
  gcoin3->SetMarkerColor(4);
  gcoin3->SetMarkerStyle(20);
  gcoin3->SetFillStyle(3005);  
  gcoin3->SetFillColor(4);  
  
  double offset=mean_pi2;
  
  for(int i=0;i<bin_run;i++){
    
    hcoin_s[i]=hcoin->ProjectionY(Form("hcoin_s_%d",i),i,i+1);
    
    fpi[i]=new TF1(Form("fpi_%d",i),"gausn(0)",min_ct,max_ct);
    fpi2[i]=new TF1(Form("fpi2_%d",i),"gausn(0)",min_ct,max_ct); 
    fpi[i]->SetLineColor(2);
    fpi2[i]->SetLineColor(2);
    hcoin_s[i]->Fit(Form("fpi_%d",i),"RQ","RQ",mean_pi-0.2,mean_pi+0.2);
    hcoin_s[i]->Fit(Form("fpi2_%d",i),"RQ","RQ",mean_pi2-0.5,mean_pi2+0.5);
    mean[i]=fpi[i]->GetParameter(1);
    mean2[i]=fpi2[i]->GetParameter(1);
    mean_err[i]=fpi[i]->GetParError(1);
    mean2_err[i]=fpi2[i]->GetParError(1);
    gcoin->SetPoint(i,i,mean[i]-3.1);
    gcoin->SetPointError(i,0,mean_err[i]);
    if(mean2_err[i]<3.0){
    gcoin2->SetPoint(i,i,-mean2[i]+mean[i]+offset-3.11);
    gcoin2->SetPointError(i,0,sqrt( pow(mean2_err[i],2)+ pow(mean_err[i],2) ));
    gcoin3->SetPoint(i,i,mean2[i]-offset);
    gcoin3->SetPointError(i,0, mean2_err[i]);
    }else {
    gcoin2->SetPoint(i,i,0);
    gcoin2->SetPointError(i,0,0);
    gcoin3->SetPoint(i,i,0);
    gcoin3->SetPointError(i,0,0);     }
  }  


  c0->cd();
  gcoin2->Draw("AP");
  c1->cd();
  gcoin->Draw("AP");
  //  gcoin2->Draw("P");
  gcoin3->Draw("P");



  c2->Divide(4,4);
  c3->Divide(4,4);
  c4->Divide(4,4);

  
  for(int i=0;i<16;i++){
    c2->cd(i+1);
    hcoin_s[i]->Draw();
    fpi[i]->Draw("same");
    //hcoin_s[i]->GetXaxis()->SetRangeUser(-5,5);
    hcoin_s[i]->GetXaxis()->SetRangeUser(3665,3675);
    c3->cd(i+1);
    hcoin_s[i+16]->Draw();
    fpi[i+16]->Draw("same");
    hcoin_s[i+16]->GetXaxis()->SetRangeUser(3665,3675);    
    //hcoin_s[i+16]->GetXaxis()->SetRangeUser(-5,5);
    c4->cd(i+1);
    hcoin_s[i+32]->Draw();
    fpi[i+32]->Draw("same");
    hcoin_s[i+32]->GetXaxis()->SetRangeUser(3665,3675);
    //hcoin_s[i+32]->GetXaxis()->SetRangeUser(-5,5);
  }
  

  
}
