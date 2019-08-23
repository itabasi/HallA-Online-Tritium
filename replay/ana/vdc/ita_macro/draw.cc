

void draw(){

  string ifname  ="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/VDC/vdct0_2run.root";
  string ifname2 ="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/VDC/vdct0_10run.root";
  TFile * fin =new TFile(ifname.c_str(),"read");  
  TFile * fin2 =new TFile(ifname2.c_str(),"read");  


  string ifname3  ="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/VDC/initial_10run/111180-111189.root";
  string ifname4 ="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/VDC/initial_10run/111190-111199.root";
  TFile * fin3 =new TFile(ifname3.c_str(),"read");  
  TFile * fin4 =new TFile(ifname4.c_str(),"read");  

  string ifname5  ="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/VDC/initial_2run/11116X_t0tuned.root";
  string ifname6 ="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/VDC/initial_10run/111160-111169.root";
  TFile * fin5 =new TFile(ifname5.c_str(),"read");  
  TFile * fin6 =new TFile(ifname6.c_str(),"read");  


  
  TGraphErrors* g1=(TGraphErrors*)fin->Get("gLu1_ac_200");
  TGraphErrors* g2=(TGraphErrors*)fin2->Get("gLu1_ac_200");
  TH1D* h1 =(TH1D*)fin3->Get("hLu1_rtime_50");
  TF1* f1=(TF1*)fin3->Get("fLu1_rt0_50");
  h1->SetLineColor(2);
  h1->SetFillColor(2);
  h1->SetFillStyle(3002);
  f1->SetLineColor(2);
  TH1D* h2 =(TH1D*)fin4->Get("hLu1_rtime_50");
  TF1* f2=(TF1*)fin4->Get("fLu1_rt0_50");
  h2->SetLineColor(4);
  h2->SetFillColor(4);
  h2->SetFillStyle(3002);  
  f2->SetLineColor(4);  

  g1->SetMarkerSize(1.);
  g1->SetMarkerColor(2);
  g1->SetMarkerColor(1);
  g2->SetMarkerSize(1.);
  g2->SetMarkerColor(4);



  TChain * T5=new TChain("T");
  //  TChain * T6=new TChain("T");

  double Lu1rt5[500],Lu1rt6[500];
  T5->Add(ifname5.c_str());
  T5->SetBranchStatus("*",0);
  T5->SetBranchStatus("Lu1_rt_p"          ,1);
  T5->SetBranchAddress("Lu1_rt_p"          ,Lu1rt5);

  //  T6->Add(ifname6.c_str());
  //  T6->SetBranchStatus("*",0);
  //  T6->SetBranchStatus("Lu1_rt_p"          ,1);
  //  T6->SetBranchAddress("Lu1_rt_p"          ,Lu1rt6);

  
  //  int  ENum5=T5->GetEntries();
  //  int  ENum6=T6->GetEntries();
  //  TH1F* h5=new TH1F("h5","",900,2100,3000);
  //  TH1F* h6=new TH1F("h6","",1000,-900,100);
  
  //  for(int i=0;i<ENum5;i++){T5->GetEntry(i);h5->Fill(Lu1rt5[50]+2936.91);}
  //  for(int i=0;i<ENum6;i++){T6->GetEntry(i);h6->Fill(Lu1rt6[200]);}
  //  TH1D* h6 =(TH1D*)fin6->Get("hLu1_rtime_50");
  //  ENum6=T6->GetEntries();



  
  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  g1->Draw("AP");
  g2->Draw("P");

  TCanvas* c1=new TCanvas("c1","c1");
  c1->cd();
  h1->Draw();
  h2->Draw("same");
  f1->Draw("same");
  f2->Draw("same");

  /*
  TCanvas* c2=new TCanvas("c2","c2");
  c2->cd();
  h5->SetLineColor(2);
  h6->SetLineColor(6);
  h5->Draw();
  h6->Draw("same");
  */


  
}

