#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
using namespace std;

#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGaxis.h"
#include <TMinuit.h>
#include <TFitter.h>
//////////////////////////// main //////////////////////////////////////////
int main(int argc, char** argv){
  TApplication *theApp = new TApplication("App", &argc, argv);

  //TMinuit *min = new TMinuit();
  TMinuit *min = new TMinuit(1);
  //TFitter *min = new TFitter();
  int p_level = min->SetPrintLevel(-1);
  delete min;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}
