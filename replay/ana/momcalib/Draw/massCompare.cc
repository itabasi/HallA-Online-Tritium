

void massCompare(){

  //  TFile*f1=new TFile("../root/Lambda_small_H1_woRangcalib.root");
  //  TFile*f2=new TFile("../root/Lambda_small_H1_woMomcalib2.root" );


  //  TFile*f1=new TFile("../root/Lambda_small_H_woMomcalib_new.root");
  TFile* f1=new TFile(Form("/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/mmass/ana_Lambda/Lambda_small_H_0th.root"));
  //  TFile*f1=new TFile("../root/Lambda_small_H_woMomcalib_new.root");
  TFile*f2=new TFile("../root/momcalib_4th_0909_0.root" );


  TH1F* hist_a=(TH1F*)f1->Get("h_mm_L");
  TH1F* hist_b=(TH1F*)f2->Get("hmm_cut");    


  TCanvas*c0=new TCanvas("c0","c0");
  c0->cd();
  hist_a->SetLineColor(4);
  hist_a->SetFillColor(4);
  //  hist_a->SetFillStyle(3002);
  hist_b->SetLineColor(2);
  hist_b->SetFillColor(2);
  //  hist_b->SetFillStyle(3002);
  hist_b->Draw();
  hist_a->Draw("same");
  //  hist_b->Draw("same");
  
}
