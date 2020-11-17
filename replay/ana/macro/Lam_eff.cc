// This macro estimate kaon efficiency //


void Lam_eff(){


  string ifr = "../rootfiles/mmass/ana_Lambda/2020-09-06/Lambda_small_OleH_all.root";


  TFile* ifp =new TFile(ifr.c_str());

  TChain* T = new TChain("T");
  TH1D* hmm =new TH1D("hmm","",300,-100,200);
  TH1D* hmm_acc =new TH1D("hmm_acc","",300,-100,200);
  TH1D* hmm_k =new TH1D("hmm_k","",300,-100,200);
  TH1D* hmm_acc_k =new TH1D("hmm_acc_k","",300,-100,200);
  
  T->Add(ifr.c_str());
  int nrun;
  int pid_cut,z_cut,ct_cut;
  double ct,mm;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("pid_cut",1);
  T->SetBranchAddress("pid_cut",&pid_cut);
  T->SetBranchStatus("z_cut",1);
  T->SetBranchAddress("z_cut",&z_cut);
  T->SetBranchStatus("ct_cut",1);
  T->SetBranchAddress("ct_cut",&ct_cut);  
  T->SetBranchStatus("mm",1);
  T->SetBranchAddress("mm",&mm);
  T->SetBranchStatus("ct_c",1);
  T->SetBranchAddress("ct_c",&ct);


  int ENum = T->GetEntries();
  
  cout<<"ENum : "<<ENum<<endl;
  //  ENum=100000;  //test
  for(int n=0;n<ENum;n++){

    T->GetEntry(n);

    //cout<<"ct "<<ct<<"z_cut "<<z_cut<<" pid_cut "<<pid_cut<<endl;
    if(pid_cut>0 && z_cut>0 && ct_cut)hmm_k->Fill(mm);
    if(z_cut>0 && ct_cut)hmm->Fill(mm);    
    double ct_acc =-100.;
    ct_acc= ct;

    if(((-50.<ct_acc && ct_acc<-30.) || (30.<ct_acc && ct_acc<50.)) && z_cut>0){
      if(pid_cut>0 && z_cut>0)hmm_acc_k->Fill(mm,0.1);
      hmm_acc->Fill(mm,0.1);
      
    }
    

    
  }


  hmm->Add(hmm_acc,-1);
  hmm_k->Add(hmm_acc_k,-1);


  
    TCanvas* c0=new TCanvas("c0","c0");
    c0->cd();
    hmm->SetLineColor(1);
    hmm_acc->SetLineColor(2);
    hmm_acc->SetFillColor(2);
    hmm_acc->SetFillStyle(3002);
    hmm->Draw();
    //hmm_acc->Draw("same");
    

    TCanvas* c1=new TCanvas("c1","c1");
    c1->cd();
    hmm_k->SetLineColor(1);
    hmm_acc_k->SetLineColor(2);
    hmm_acc_k->SetFillColor(2);
    hmm_acc_k->SetFillStyle(3002);
    hmm_k->Draw();
    //    hmm_acc_k->Draw("same");
    

    //    double Nmm = hmm->Integral(hmm->FindBin(-100.),hmm->FindBin(100.));
    //    double Nmm_k = hmm_k->Integral(hmm_k->FindBin(-100.),hmm_k->FindBin(100.));

    double Nmm = hmm->Integral(hmm->FindBin(-5.),hmm->FindBin(5.));
    double Nmm_k = hmm_k->Integral(hmm_k->FindBin(-5.),hmm_k->FindBin(5.));
    cout<<"Integral hmm "<<Nmm<<endl;
    cout<<"Integral hmm_k "<<Nmm_k<<endl;



  
}
