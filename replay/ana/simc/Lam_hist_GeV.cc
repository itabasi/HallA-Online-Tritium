

void Lam_hist_GeV(){

  string ifr_L = "../rootfiles/simc/1H_Lam.root";
  string ifr_S = "../rootfiles/simc/1H_Sig.root";

  TFile *ifL =new TFile(ifr_L.c_str());
  TFile *ifS =new TFile(ifr_S.c_str());

  TTree* TL= (TTree*)ifL->Get("SNT");
  int ENumL=  TL->GetEntries();
  cout<<"ENumL "<<ENumL<<endl;
  TH1D* hL=new TH1D("hL","",250,-0.3,0.2);
  TH1D* hL_MeV=new TH1D("hL_MeV","",250,-300,200);

  double ML = 1115.683e-3; // GeV
  float weight,normfac,mm;

  TL->SetBranchAddress("mmnuc",&mm);
  TL->SetBranchAddress("Weight",&weight);
  TL->SetBranchAddress("Normfac",&normfac);

  for(int nev=0;nev<ENumL;nev++){
    mm=-1000.0;
    weight=0.0;
    normfac=0.0;
    TL->GetEntry(nev);
    hL->Fill(mm-ML, 1./weight/normfac);
    hL_MeV->Fill((mm-ML)*1000., 1./weight/normfac);

  }


  TTree* TS= (TTree*)ifS->Get("SNT");
  int ENumS=  TS->GetEntries();
  TH1D* hS_L=new TH1D("hS_L","",250,-0.300,0.200);
  TH1D* hS=new TH1D("hS","",250,-0.300,0.200);
  TH1D* hS_MeV=new TH1D("hS_MeV","",250,-300,200);
  TH1D* hS_L_MeV=new TH1D("hS_L_MeV","",250,-300,200);
  cout<<"ENumS "<<ENumS<<endl;
  double MS0 = 1192.642e-3; //GeV

  float weight_s,normfac_s,mm_s;

  TS->SetBranchAddress("mmnuc",&mm_s);
  TS->SetBranchAddress("Weight",&weight_s);
  TS->SetBranchAddress("Normfac",&normfac_s);

  for(int nev=0;nev<ENumL;nev++){
    mm_s=-1000.;
    weight_s=0.0;
    normfac_s=0.0;
    TS->GetEntry(nev);
    hS_L->Fill(mm_s-ML, 1./weight_s/normfac_s);
    hS->Fill(mm_s-MS0, 1./weight_s/normfac_s);
    hS_L_MeV->Fill((mm_s-ML)*1000., 1./weight_s/normfac_s);
    hS_MeV->Fill((mm_s-MS0)*1000., 1./weight_s/normfac_s);
  }

  hL->Add(hS_L);
  hL->Scale(1./100.);
  hL_MeV->Add(hS_L_MeV);
  hL_MeV->Scale(1./100.);
  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  hL->Draw();
  //  hS->Draw("same");

  TFile* ofr =new TFile("../rootfiles/simc/1Hsimc.root","recreate");
  hL->Write();
  hS->Write();
  hL_MeV->Write();
  hS_MeV->Write();
}
