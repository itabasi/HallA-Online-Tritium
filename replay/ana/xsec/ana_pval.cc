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

  //  TFile *ifp = new TFile("hist.root","read");
  //  TH1D *h_mm = (TH1D*)gROOT->FindObject("hmm_p");

  TFile *ifp = new TFile(ifname.c_str(),"read");
  TH1D *h_nnL = (TH1D*)gROOT->FindObject("h_mm_nnL");

  double min_mm = h_nnL->GetXaxis()->GetXmin();
  double max_mm = h_nnL->GetXaxis()->GetXmax();
  const int Nbin = (int)(max_mm - min_mm)*30;
  TH1D *h_mm =new TH1D("h_mm","",Nbin,min_mm,max_mm);

  double S[Nbin], N[Nbin], Pval[Nbin], MM[Nbin], PS[Nbin];
  double obin = (h_mm->GetBinCenter(2) - h_mm->GetBinCenter(1));
  double width = 1.00;
  double window = width / obin;
  double bgr = 2; // background range with sigma
  cout << "One bin size: " << obin   << " MeV" << endl;
  cout << "peak resolution: " << window <<" "<<endl;
  
  //===< Get Tree >=========//

  TChain* T =new TChain("T");
  T->Add(ifname.c_str());
  int ENum = T->GetEntries();
  cout<<"ENum "<<ENum<<endl;
  int z_cut,pid_cut,ct_cut;
  double mm_nnL;
  T->SetBranchAddress("pid_cut",&pid_cut);
  T->SetBranchAddress("z_cut"  ,&z_cut);
  T->SetBranchAddress("ct_cut" ,&ct_cut);
  T->SetBranchAddress("mm_nnL",&mm_nnL);


  
  for(int iev=0;iev<ENum;iev++){

    z_cut=-1;
    ct_cut=-1;
    pid_cut=-1;
    mm_nnL=-100.;
    
    T->GetEntry(iev);
    if(z_cut>0 && ct_cut>0 && pid_cut>0){
      h_mm->Fill(mm_nnL);
    }

  }

  cout<<"End Fill "<<endl;
  

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

    S[i]    = pk;
    N[i]    = bg;
    Pval[i] = pval;
    PS[i]   = ps;
    MM[i]   = mm;
    f->Clear();
  }

  TGraph *g_s   = new TGraph(Nbin,MM,S);
  g_s->SetName("g_s");
  TGraph *g_n   = new TGraph(Nbin,MM,N);
  g_s->SetName("g_n");
  TGraph *g_val = new TGraph(Nbin,MM,Pval);
  g_s->SetName("g_val");
  TGraph *g_ps  = new TGraph(Nbin,MM,PS);
  g_s->SetName("g_ps");
  set->SetGr(g_s  ,"g_s","#DeltaB_{#Lambda} (MeV/c^{2})^{2}","");
  set->SetGr(g_n  ,"g_n","#DeltaB_{#Lambda} (MeV/c^{2})^{2}","",2);
  set->SetGr(g_val,"g_val","#DeltaB_{#Lambda} (MeV/c^{2})^{2}","p-value");
  set->SetGr(g_ps ,"g_ps","#DeltaB_{#Lambda} (MeV/c^{2})^{2}","peak significance");
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

  string orname="./test.root";
  TFile* ofr = new TFile(orname.c_str(),"recreate");

  // TGraph Point Setting //
  g_ps->SetMarkerStyle(7);
  g_ps->SetMarkerColor(2);

  //  h_mm->Write();
  ofr->Write();
  g_s->Write();
  g_n->Write();
  g_ps->Write();
  g_val->Write();
  ofr->Close();    

  
  gSystem->Exit(1);
  theApp.Run();
  return 0;

}

