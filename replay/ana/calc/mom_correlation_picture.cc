

void mom_correlation_picture(){


  double pe[6]={1200,2081,2180,2278};
  double Lk[6]={2842,1907,1798,1688};
  double Sk[6]={2736,1786,1675,1561};
  double nnLk[6]={2885,1975,1872,1769};
  TGraph* gL=new TGraph();
  TGraph* gS=new TGraph();
  TGraph* gnnL=new TGraph();

  for(int i=0;i<3;i++){
    gL->SetPoint(i,pe[i],Lk[i]);
  gS->SetPoint(i,pe[i],Sk[i]);
  gnnL->SetPoint(i,pe[i],nnLk[i]);}

  TF1* fL=new TF1("fL","pol1(0)",1900,2400);
  TF1* fS=new TF1("fS","pol1(0)",1900,2400);
  TF1* fnnL=new TF1("fnnL","pol1(0)",1900,2400);
  gS->Fit("fS","","",1900,2400);
  gL->Fit("fL","","",1900,2400);
  gnnL->Fit("fnnL","","",1900,2400);
  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  //  gL->GetXaxis()->SetRangeUser(1600,2500);
  // gL->GetYaxis()->SetRangeUser(1600,2000);
  //  gnnL->SetMarkerSize(1.0);
  
  // gnnL->SetMarkerStyle(21);
  // gnnL->SetFillColor(4);
  //gnnL->SetMarkerColor(4);
  //gnnL->SetFillStyle(1.0);
  gS->SetLineColor(2);
  gS->SetLineWidth(2);
  gL->SetLineColor(4);
  gnnL->SetLineColor(1);

  fL->Draw();
  fS->Draw("same");
  fnnL->Draw("same");
  // gL->Draw("APL");
  // gS->Draw("PL");
  // gnnL->Draw("PL");
}
