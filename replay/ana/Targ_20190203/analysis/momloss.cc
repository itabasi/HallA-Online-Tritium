/*
  momloss.cc
  
  Toshiyuki Gogami, November 23, 2018
*/

double mk = 493.677;
double me = 0.511;

void momloss(){
  // ====================================== //
  // ======= General conditions =========== //
  // ====================================== //
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);


  // =========================== //
  // ===== Open ROOT file ====== //
  // =========================== //
  double mass = 0.0;
  int flag = 0;
  TFile* f1 = new TFile("./root/kaon1.82GeVc_right_noRaster_highstat.root");mass = mk; flag = 2;
  //TFile* f1 = new TFile("./root/scate2.218GeVc_left_noRaster_highstat.root");mass = me; flag = 1;
  //TFile* f1 = new TFile("./root/beam4.32GeVc_noRaster.root");mass = me; flag = 3;
  TTree* t1 = (TTree*)f1->Get("tree");

  
  // ================================ //
  // ====== Create histograms  ====== //
  // ================================ //
  double xmin = 0.0, xmax = 30.0;
  int xbin = 1000;
  //int xbin = 2000;
  TH1F* h1 = new TH1F("h1","",xbin,xmin,xmax);
  h1->GetXaxis()->SetTitle("dp (MeV)");
  //h1->GetYaxis()->SetTitle("Counts / 15 keV");
  h1->GetYaxis()->SetTitle("Counts / 30 keV");
  h1->GetXaxis()->SetRangeUser(0.0,3.0);
  h1->SetFillStyle(1001);
  h1->SetFillColor(9);

  TH1F* h1_ = (TH1F*)h1->Clone("h1_");
  
  TH1F* h2 = (TH1F*)h1->Clone("h2");
  h2->GetXaxis()->SetTitle("dE (MeV)");

  TH2F* h3 = new TH2F("h3","",100,-20.0,20.0,xbin,xmin,xmax);
  h3->GetYaxis()->SetTitle("dp (MeV)");
  h3->GetXaxis()->SetTitle("z (cm)");
  h3->GetYaxis()->SetRangeUser(0.0,3.0);
  h3->GetXaxis()->SetRangeUser(-15.,15.0);

  const int max = 20;
  double pid[max];
  double pmom[max], pene[max];
  double bmom, bene;
  double z;
  double dE;
  double ent = t1->GetEntries();
  //ent = 5000;

  if(flag==2){
    t1->SetBranchAddress("rpid",   &pid);
    t1->SetBranchAddress("rp",     &pmom);
  }
  else if (flag==1){
    t1->SetBranchAddress("lpid",   &pid);
    t1->SetBranchAddress("lp",     &pmom);
  }
  else{
    t1->SetBranchAddress("bpid",   &pid);
    t1->SetBranchAddress("bp",     &pmom);
  }
  t1->SetBranchAddress("pBeam",  &bmom);
  t1->SetBranchAddress("zBeam",  &z);
  bool ok = false;
  double gopid = 0;
  if(flag==1) gopid = 0.0;
  else if(flag==2) gopid = 2.0;
  else gopid = 0.0;
  
  for(int i=0 ; i<ent ; i++){
    for(int j=0 ; j<max ; j++){
      pid[j]  = -1;
      pmom[j] = -2222.0;
      pene[j] = -2222.0;
    }
    z = -2222.0;
    ok = false;
    dE = -1.0;

    t1->GetEntry(i);
    //cout << bmom-pmom[0] << endl;
    
    if(pid[0]==gopid) ok = true;
    else ok = false;
    
    if(ok==true){
      if(flag==3) bmom = -bmom;
      h1 ->Fill(bmom-pmom[0]);
      h3 ->Fill(z,bmom-pmom[0]);
      //if(bmom-pmom[0]>0.6){
	if(z<9.0){
	h1_->Fill(bmom-pmom[0]);
      }

      bene    = sqrt(mass*mass + bmom*bmom);
      for(int j=0 ; j<max ; j++){
	pene[j] = sqrt(mass*mass+pmom[j]*pmom[j]);
      }
      dE = bene - pene[0];
      //cout << dE << endl;
      h2->Fill(dE);
    }
  }

  TCanvas* c1 = new TCanvas("c1","c1");
  h1->Draw();

  TCanvas* c2 = new TCanvas("c2","c2");
  h2->Draw();

  c1->cd();
  TF1* func1 = new TF1("func1","landau",0.0,10.0);
  func1->SetNpx(2000);
  func1->SetParameters(10000.,1.1,9.8e-2);
  //h1->Fit("func1","","",0.7,1.5);

  TCanvas* c3 = new TCanvas("c3","c3");
  h3->Draw("col");

  TCanvas* c4 = new TCanvas("c4","c4");
  h1_->Draw();

  // ===== Print ====== //
  //c1->Print("momloss_1.png");
  //c2->Print("momloss_2.png");
  //c3->Print("momloss_3.png");
  //c4->Print("momloss_4.png");
  
   
}

