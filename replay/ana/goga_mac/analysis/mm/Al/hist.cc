/*
  hist.cc
  
  T. Gogami, Aug 9, 2019
*/


void hist(){
  gStyle->SetOptStat(0);
  //gROOT->SetSetStyle("");
  
  // ============================== //
  // ===== Open ROOT file ======== //
  // ============================== //
  //TFile* f1 = new TFile("T2.root");
  TFile* f1 = new TFile("Al_20190813.root");
  TTree* t1 = (TTree*)f1->Get("tree");
  
  
  // ============================ //
  // ===== 
  // ============================ //
  const int max = 100;
  double a1, a2;
  double mm[max];
  double ctime[max];
  double vz[max];
  
  t1->SetBranchAddress("ctime", &ctime);
  t1->SetBranchAddress("mm", &mm);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);
  t1->SetBranchAddress("vz_mean", &vz);

  TH1F* hist1  = new TH1F("hist1","",200,-100.0,100.0);
  hist1->GetXaxis()->SetTitle("-B_{#Lambda} (MeV)");
  hist1->GetYaxis()->SetTitle("Counts /  MeV");
  //TH1F* hist1  = new TH1F("hist1","",200,-100.0,100.0);
  TH1F* hist1_ = (TH1F*)hist1->Clone("hist1_");
  TH1F* hist2  = new TH1F("hist2","",300,-15.0,15.0);
  hist2->GetXaxis()->SetTitle("coin time (ns)");
  hist2->GetYaxis()->SetTitle("Counts / 100 ps");
  TH1F* hist2_ = (TH1F*)hist2->Clone("hist2_");
  TH1F* hist2__= (TH1F*)hist2->Clone("hist2__");
  TH1F* hist2_acc= (TH1F*)hist2->Clone("hist2_acc");

  hist2->GetXaxis()->SetRangeUser(-1.0,1.0);
  
  hist1_->SetFillColor(9);
  hist1_->SetFillStyle(1001);
  hist2->SetLineColor(2);
  hist2_->SetLineColor(1);
  hist2__->SetLineColor(1);
  hist2_->SetFillColor(9);
  hist2_->SetFillStyle(1001);
  hist2_acc->SetFillColor(8);
  hist2_acc->SetFillStyle(3001);
  
  double step = 2.0;
  
  double ent = t1->GetEntries();
  
  for(int i=0 ; i<ent ; i++){
    for(int j=0 ; j<max ; j++){
      ctime[j] = -2222.0;
      mm[j]    = -2222.0;
      vz[j]    = -2222.0;
    }
    a1 = -22222.0;
    a2 = -22222.0;
    
    
    t1->GetEntry(i);
    
    if(a1<1.0
       && a2>3.0
       && a2<13.0
       && fabs(vz[0])>0.05 
       && fabs(vz[0])<0.25 ){
      hist2__->Fill(ctime[0]);
      
      if(fabs(ctime[0]+step*4.0)<1.0
	 || fabs(ctime[0]+step*5.0)<1.0
	 || fabs(ctime[0]+step*6.0)<1.0
	 || fabs(ctime[0]+step*7.0)<1.0
	 || fabs(ctime[0]-step*2.0)<1.0
	 || fabs(ctime[0]-step*7.0)<1.0
	 ){
	hist1_->Fill(mm[0]);
	hist2_->Fill(ctime[0]);
	if(fabs(ctime[0]+step*4.0)<1.0) hist2_acc->Fill(ctime[0]+step*4.0);
	else if(fabs(ctime[0]+step*6.0)<1.0) hist2_acc->Fill(ctime[0]+step*6.0);
	else if(fabs(ctime[0]+step*7.0)<1.0) hist2_acc->Fill(ctime[0]+step*7.0);
	else if(fabs(ctime[0]-step*2.0)<1.0) hist2_acc->Fill(ctime[0]-step*2.0);
	else if(fabs(ctime[0]-step*7.0)<1.0) hist2_acc->Fill(ctime[0]-step*7.0);

      }
      
      if(abs(ctime[0])<1.0){
	hist1 ->Fill(mm[0]);
	hist2->Fill(ctime[0]);
      }
    }
    
  }

  hist1_->Scale(1./6.);
  hist2_acc->Scale(1./6.);

  hist1 ->SetLineColor(1);
  hist1_->SetMarkerStyle(25);
  hist1_->SetMarkerSize(0.5);
  hist1_->SetLineColor(1);
  
  
  TCanvas* c1 = new TCanvas("c1","c1");
  hist1 ->Draw();
  hist1_->Draw("same");

  TCanvas* c2 = new TCanvas("c2","c2");
  hist2__->Draw();
  //hist2_acc->Draw("same");
  hist2_->Draw("same");
  hist2->Draw("same");
  
}
  

  
