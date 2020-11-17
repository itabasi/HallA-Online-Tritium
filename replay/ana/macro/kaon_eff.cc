// This macro estimate kaon efficiency //


void kaon_eff(){


  string ifr = "../rootfiles/mmass/ana_Lambda/2020-09-06/nnL_small_Ole_all.root";


  TFile* ifp =new TFile(ifr.c_str());

  TChain* T = new TChain("T");
  TH1D* hct =new TH1D("hct","",2000,-50,50);
  TH1D* hct_acc =new TH1D("hct_acc","",2000,-50,50);
  TH1D* hct_acc2 =new TH1D("hct_acc2","",2000,-50,50);
  TH1D* hct_k =new TH1D("hct_k","",2000,-50,50);
  TH1D* hct_acc_k =new TH1D("hct_acc_k","",2000,-50,50);
  
  T->Add(ifr.c_str());
  int nrun;
  int pid_cut,z_cut;
  double ct;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("pid_cut",1);
  T->SetBranchAddress("pid_cut",&pid_cut);
  T->SetBranchStatus("z_cut",1);
  T->SetBranchAddress("z_cut",&z_cut);
  T->SetBranchStatus("ct_c",1);
  T->SetBranchAddress("ct_c",&ct);


  int ENum = T->GetEntries();
  
  cout<<"ENum : "<<ENum<<endl;
  //  ENum=100000;  //test
  for(int n=0;n<ENum;n++){

    T->GetEntry(n);

    //cout<<"ct "<<ct<<"z_cut "<<z_cut<<" pid_cut "<<pid_cut<<endl;
    if(pid_cut>0 && z_cut>0)hct_k->Fill(ct);
    if(z_cut>0)hct->Fill(ct);    
    double ct_acc =-100.;
    ct_acc= ct;

    //    if(((-50.<ct_acc && ct_acc<-30.) || (30.<ct_acc && ct_acc<50.)) && z_cut>0){
    if((-50.<ct_acc && ct_acc<-30.) && z_cut>0){
      if(pid_cut>0 && z_cut>0)hct_acc_k->Fill(ct);
      hct_acc->Fill(ct);
      while(1){
	if(-1.0<ct_acc && ct_acc<1.0 )break;
	else if(ct_acc<-1.0)ct_acc += 2.0;
	else if(ct_acc> 1.0)ct_acc +=-2.0;
      }
      //      cout<<"ct "<<ct<<" ct_acc "<<ct_acc<<endl;
      if(pid_cut>0 && z_cut>0)hct_acc_k->Fill(ct_acc,0.1);
      hct_acc->Fill(ct_acc,0.1);
    }
    

    
  }


  //  cout<<"Integral hct_acc -1 to 1"<<hct_acc->Integral(hct_acc->FindBin(-1.),hct_acc->FindBin(1.))<<endl;
  //  cout<<"Integral acc "<<hct_acc->Integral(hct_acc->FindBin(-50),hct_acc->FindBin(-30.))+hct_acc->Integral(hct_acc->FindBin(30),hct_acc->FindBin(50.))<<endl;
  //    cout<<"Integral acc "<<hct_acc->Integral(hct_acc->FindBin(-32),hct_acc->FindBin(-30.))<<endl;
  //  cout<<"Integral acc "<<hct_acc->Integral(hct_acc->FindBin(-50),hct_acc->FindBin(-30.))+hct_acc->Integral(hct_acc->FindBin(30),hct_acc->FindBin(50.))<<endl;


  //  hct_acc->Scale(0.1);

    hct->Add(hct_acc,-1);
    hct_k->Add(hct_acc_k,-1);
    TCanvas* c0=new TCanvas("c0","c0");
    c0->cd();
    hct->SetLineColor(1);
    hct_acc->SetLineColor(2);
    hct_acc->SetFillColor(2);
    hct_acc->SetFillStyle(3002);
    hct->Draw();
    hct_acc->Draw("same");
    

    TCanvas* c1=new TCanvas("c1","c1");
    c1->cd();
    hct_k->SetLineColor(1);
  hct_acc_k->SetLineColor(2);
  hct_acc_k->SetFillColor(2);
  hct_acc_k->SetFillStyle(3002);
  hct_k->Draw();
  hct_acc_k->Draw("same");






  
}
