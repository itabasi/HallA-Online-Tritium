using namespace std;
#include "Setting.h"
#include "Math/ProbFunc.h"
#include "Math/DistFunc.h"

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
int main(int argc, char** argv){

//  gErrorIgnoreLevel = kError;

  int ch;
  extern char *optarg;
  string ifname, ofname;

  while((ch=getopt(argc,argv,"f:"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      break;
    case 'h':
      cout<<"-f (inputfile): input ROOT file name"<<endl;
      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;
    }
  }

  Setting *set = new Setting();
  set->Initialize();

  TApplication theApp("App", &argc, argv);

  TFile *ifp = new TFile("hist.root","read");

  TH1D *h_mm = (TH1D*)gROOT->FindObject("hmm_p");

  const int Nbin = h_mm->GetXaxis()->GetNbins();
  double S[Nbin], N[Nbin], Pval[Nbin], MM[Nbin], PS[Nbin];
  double obin = (h_mm->GetBinCenter(2) - h_mm->GetBinCenter(1));
  double width = 1.00;
  double window = width / obin;
  double bgr = 2; // background range with sigma
  cout << "One bin size: " << obin   << " MeV" << endl;
  cout << "peak resolution: " << window <<endl;
  TF1 *f;

  for(int i = 0;i<Nbin;i++){
    Pval[i] = PS[i] = 0;
    double bg1 = h_mm->Integral(i-int((2+bgr)*window),i-int(2*window));
    double bg2 = h_mm->Integral(i+int(2*window),i+int((2+bgr)*window));
    double bg = (bg1 + bg2) / (bgr);
    double pk = h_mm->Integral(i-int(window),i+int(window));
    double mm = h_mm->GetBinCenter(i);
    f = new TF1("f_pois","TMath::Poisson(x,[0])",0,1000);
    f->SetParameter(0,bg);
    double pval = f->Integral(pk,1000.,1E-12);
    if( pval == 0. ){ pval = 1.; }
    if( i-int(6*window) <= 0    ){ pval = 1.; }
    if( i+int(6*window) >= Nbin ){ pval = 1.; }
    double ps = ROOT::Math::gaussian_quantile_c(pval,1.);
    if( pval>=1 || pval<=0 ){ ps = 0; }
    if( ps < 0 ){ ps = 0; }
//    cout<<i<<"  "<<mm<<"  "<<bg<<"  "<<pk<<"  "<<pval<<"  "<<ps<<endl;
    S[i]    = pk;
    N[i]    = bg;
    Pval[i] = pval;
    PS[i]   = ps;
    MM[i]   = mm;
    f->Clear();
  }

  TGraph *g_s   = new TGraph(Nbin,MM,S);
  TGraph *g_n   = new TGraph(Nbin,MM,N);
  TGraph *g_val = new TGraph(Nbin,MM,Pval);
  TGraph *g_ps  = new TGraph(Nbin,MM,PS);
  set->SetGr(g_s  ,"","#DeltaB_{#Lambda} (MeV/c^{2})^{2}","");
  set->SetGr(g_n  ,"","#DeltaB_{#Lambda} (MeV/c^{2})^{2}","",2);
  set->SetGr(g_val,"p-vale","#DeltaB_{#Lambda} (MeV/c^{2})^{2}","p-value");
  set->SetGr(g_ps ,"P.S","#DeltaB_{#Lambda} (MeV/c^{2})^{2}","peak significance");
  g_val->SetMinimum(1E-10);
  g_val->GetXaxis()->SetRangeUser(-100,150);
  g_ps ->GetXaxis()->SetRangeUser(-100,150);
  g_s  ->GetXaxis()->SetRangeUser(-100,150);
  g_n  ->GetXaxis()->SetRangeUser(-100,150);

  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(1,2);
  c1->cd(1); gPad->SetLogy(); g_val->Draw("AL");
  c1->cd(2); g_ps->Draw("AL");

  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c2->Divide(1,1);
  c2->cd(1);
  g_s->Draw("AL");
  g_n->Draw("sameL");

  cout<<"Done!"<<endl;
  ifp->Close();

//  gSystem->Exit(1);
  theApp.Run();
  return 0;

}

