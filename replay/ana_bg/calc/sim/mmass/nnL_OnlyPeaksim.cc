/*
Made by K.Itabashi Apr. 30th 2018
nnL peak simulation

*/
void nnL_OnlyPeaksim(){

  int n; // number of nnL signal
  n=120;
  double wid1,wid2,nwid;
 
  nwid=1.;//natural width [MeV] 
  wid1=1.5;
    //sqrt(pow(nwid,2)+pow(1.5,2));
    wid2=20;
    //sqrt(pow(nwid,2)+pow(20,2));
 
  double sig1,sig2,a_sig;
  a_sig=2*sqrt(2*log(2));

  // cout<<"a_sig is "<<a_sig<<endl;
  sig1=wid1/a_sig;
  sig2=wid2/a_sig;

  TF1* fun1=new TF1("fun1","gaus",-50,50);
  TF1* fun2=new TF1("fun2","gaus",-50,50);
  fun1->SetParameters(100,0,sig1);
  fun2->SetParameters(100,0,sig2);

  TH1F* hist1 = new TH1F("hist1","gaus",200,-50,50);
  TH1F* hist2 = new TH1F("hist2","gaus",200,-50,50);


  hist1->FillRandom("fun1",n);
  hist1->SetFillColor(0);

  hist2->FillRandom("fun2",n);
  hist2->SetFillColor(2);

  TCanvas* c0=new TCanvas("c0","c0");
 hist1->SetTitle("nnL Signal");
 hist1->SetTitle(0);
  hist1->SetXTitle("Missing Mass [MeV/c^2]");
  //hist1->GetXaxis()->CenterTitle();
  //hist1->GetYaxis()->CenterTitle();

  hist1->SetYTitle("Counts");


  hist1->Draw();
  hist2->Draw("same");
 }
