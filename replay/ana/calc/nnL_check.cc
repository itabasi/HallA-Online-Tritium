


void nnL_check(){


  string ifr_10_name = "../rootfiles/momcalib/momcalib_5th_wAl_037_w10.root";
  string ifr_20_name = "../rootfiles/momcalib/momcalib_5th_wAl_037_w20.root";
  string ifr_30_name = "../rootfiles/momcalib/momcalib_5th_wAl_037.root";


  TFile *ifr_10 =new TFile(ifr_10_name.c_str());
  TFile *ifr_20 =new TFile(ifr_20_name.c_str());
  TFile *ifr_30 =new TFile(ifr_30_name.c_str());
  TH1F* h00 = (TH1F*)ifr_10->Get("hmm_nnL");
  TH1F* h10 = (TH1F*)ifr_10->Get("hmm_nnL_cut");
  TH1F* h20 = (TH1F*)ifr_20->Get("hmm_nnL_cut");
  TH1F* h30 = (TH1F*)ifr_30->Get("hmm_nnL_cut");

  h00->SetLineColor(1);
  h00->SetFillStyle(0);
  h10->SetLineColor(2);
  h10->SetFillStyle(0);
  h20->SetLineColor(3);
  h20->SetFillStyle(0);
  h30->SetLineColor(4);
  h30->SetFillStyle(0);
  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  h00->Draw();
  h10->Draw("same");
  h20->Draw("same");
  h30->Draw("same");

  TGraph * g23=new TGraph();
  TGraph * g03=new TGraph();

  double y23[4], y03[4];

  y03[0]= h00->Integral(h00->GetXaxis()->FindBin(1),h00->GetXaxis()->FindBin(5));
  g03->SetPoint(0,0,y03[0]);
  y03[1]= h10->Integral(h10->GetXaxis()->FindBin(1),h10->GetXaxis()->FindBin(5));
  g03->SetPoint(1,10,y03[1]);
  y03[2]= h20->Integral(h20->GetXaxis()->FindBin(1),h20->GetXaxis()->FindBin(5));
  g03->SetPoint(2,20,y03[2]);
  y03[3]= h30->Integral(h30->GetXaxis()->FindBin(1),h30->GetXaxis()->FindBin(5));
  g03->SetPoint(3,30,y03[3]);


  y23[0]= h00->Integral(h00->GetXaxis()->FindBin(21),h00->GetXaxis()->FindBin(25));
  g23->SetPoint(0,0,y23[0]);
  y23[1]= h10->Integral(h10->GetXaxis()->FindBin(21),h10->GetXaxis()->FindBin(25));
  g23->SetPoint(1,10,y23[1]);
  y23[2]= h20->Integral(h20->GetXaxis()->FindBin(21),h20->GetXaxis()->FindBin(25));
  g23->SetPoint(2,20,y23[2]);
  y23[3]= h30->Integral(h30->GetXaxis()->FindBin(21),h30->GetXaxis()->FindBin(25));
  g23->SetPoint(3,30,y23[3]);



  TCanvas* c1=new TCanvas("c1","c1");
  c1->cd();
  g03->SetMarkerStyle(20);
  g03->Draw("AP");
  g23->SetMarkerStyle(20);
  g23->SetMarkerColor(2);
  g23->Draw("AP");

}
