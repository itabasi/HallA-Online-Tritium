/*
  ofseteval.cc
  Toshi Gogami, Aug 29, 2019
*/

void offseteval(){
  ifstream* ifs = new ifstream("scale_offset_20190829.dat");
  const int nfoil = 10;
  double y1[nfoil];
  double y2[nfoil];
  double x[nfoil];
  double temp;
  for(int i=0 ; i<nfoil ; i++){
    *ifs >> temp >> y1[i] >> temp >> y2[i] >> temp;
    cout << y1[i] << " " << y2[i] << endl;
    x[i] = i;
  }

  TGraph* gr1 = new TGraph(nfoil,x,y1);
  gr1->SetName("gr1");
  TGraph* gr2 = new TGraph(nfoil,x,y2);
  gr2->SetName("gr2");

  TCanvas* c1 = new TCanvas("c1","c1");
  gr1->Draw("a*");
  
  TCanvas* c2 = new TCanvas("c2","c2");
  gr2->Draw("a*");
  

  TF1* func1 = new TF1("func1","[0]",0.0,10.0);
  gr1->Fit("func1","","",0.5,8.5);
  TF1* func2 = new TF1("func2","[0]+[1]*x",0.0,10.0);
  gr2->Fit("func2","","",0.5,8.5);
  
  cout << "1-0: " << func1->Eval(0) << ", "
       << "1-9: " << func1->Eval(9) << endl;
  cout << "2-0: " << func2->Eval(0) << ", "
       << "2-9: " << func2->Eval(9) << endl;
  
  
  


  
}
