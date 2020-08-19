

void fsi_param(){

  TGraphErrors* g1 =new TGraphErrors();
  g1->SetMarkerStyle(3);
  int nmax =9;
  
  double  p1[nmax];
  p1[0] = 6.0;   // 0 MeV
  p1[1] = 5.0;   // 25 MeV
  p1[2] = 3.4;   // 50 MeV
  p1[3] = 2.4;   // 75 MeV
  p1[4] = 1.8;   // 100 MeV
  p1[5] = 1.5;   // 125 MeV
  p1[6] = 1.4;   // 150 MeV
  p1[7] = 1.2;   // 150 MeV
  p1[8] = 1.2;   // 175 MeV

  double p[nmax];

  p[0] = 0.0;
  p[1] = 25.0;
  p[2] = 50.;
  p[3] = 75.;
  p[4] = 100.;
  p[5] = 125.;
  p[6] = 150.;
  p[7] = 175.;
  p[8] = 200;

  for(int i=0; i<nmax;i++)
    g1->SetPoint(i,p[i],p1[i]);


  TCanvas* c0 = new TCanvas("c0","c0");

  c0->cd();
  g1->Draw("AP");

  //  ofstream ofr("./test.root","recreate");`
  //  g1->Write();
  
}
