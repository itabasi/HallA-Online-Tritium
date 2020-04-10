

void Lam_hist(){

  string ifr_L = "../rootfiles/simc/1H_Lam.root";
  string ifr_S = "../rootfiles/simc/1H_Sig.root";

  TFile *ifL =new TFile(ifr_L.c_str());
  TFile *ifS =new TFile(ifr_S.c_str());

  TTree* TL= (TTree*)ifL->Get("SNT");
  int ENumL=  TL->GetEntries();
  cout<<"ENumL "<<ENumL<<endl;
  TH1D* hL=new TH1D("hL","",250,-300,200);
  
  double ML = 1115.683; // MeV

  float weight,normfac,mm;

  TL->SetBranchAddress("mmnuc",&mm);
  TL->SetBranchAddress("Weight",&weight);
  TL->SetBranchAddress("Normfac",&normfac);

  for(int nev=0;nev<ENumL;nev++){
    mm=-1000.0;
    weight=0.0;
    normfac=0.0;
    TL->GetEntry(nev);
    hL->Fill(mm*1000.-ML, 1./weight/normfac);
  }


  TTree* TS= (TTree*)ifS->Get("SNT");
  int ENumS=  TS->GetEntries();
  TH1D* hS=new TH1D("hS","",250,-300,200);
  cout<<"ENumS "<<ENumS<<endl;


  float weight_s,normfac_s,mm_s;

  TS->SetBranchAddress("mmnuc",&mm_s);
  TS->SetBranchAddress("Weight",&weight_s);
  TS->SetBranchAddress("Normfac",&normfac_s);

  for(int nev=0;nev<ENumL;nev++){
    mm_s=-1000.;
    weight_s=0.0;
    normfac_s=0.0;
    TS->GetEntry(nev);
    hS->Fill(mm_s*1000.-ML, 1./weight_s/normfac_s);
  }

  hL->Add(hS);
  hL->Scale(1./100.);

  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  hL->Draw();
  //  hS->Draw("same");

  TFile* ofr =new TFile("../rootfiles/simc/1Hsimc.root","recreate");
  hL->Write();

}
