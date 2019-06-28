#include <iostream>
#include <string>
using namespace std;
#include "TApplication.h"
#include <TMinuit.h>
//#include <TFitter.h>
void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  f = par[0]*par[0] + par[0] + 1;
}

 int main(int argc, char** argv){
    TApplication *theApp =new TApplication("App",&argc,argv);


    
  TMinuit *min = new TMinuit(1);
  cout<<"TMinuit is OK"<<endl;
  int p_level= min->SetPrintLevel(1);
  cout<<"SetPrintLevel"<<endl;
  delete min;
  /*
  double par, parErr;
  string parName = "test_para";
  double stepSize = 0.1, minVal = -100.0, maxVal = 100.0;

  min->DefineParameter(0, parName.c_str(), par, stepSize, minVal, maxVal);

  min->SetFCN(chi2);

  int migrad_stats = min->Migrad();

  min->GetParameter(0, par, parErr);

  cout << "Result: " <<par << " +/- " << parErr << endl;
  cout << "Status of Migrad: " << migrad_stats << endl;

  delete min;
  */
   theApp->Run();


  return 0;
}
