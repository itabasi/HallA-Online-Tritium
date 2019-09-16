void QFhist(){
  TCanvas *QFhist=new TCanvas("QFhist","QFhist");
  TFile*f1 = new TFile("test_kaon_he3_goga3.root");
  TTree*t1 = (TTree*)f1->Get("h666");
  float mm;
  t1->SetBranchAddress("mmnuc",&mm);
  ////////////////////////////hist////////////////////////////////
  TH1F*h1=new TH1F("h1","bounding energy",200,-0.01,0.09); 
  h1->SetFillColor(3);
  ///////////////////// //QF hist////////////////////////////////////////
  int ent;
  ent=t1->GetEntries();
  for(int i=0;i<2000;i++){
    t1->GetEntry(i);
    h1->Fill(mm-0.938272-0.939565-1.115683);
  }
  h1->GetXaxis()->SetTitle("Binding Energy [GeV]");
  h1->GetYaxis()->SetTitle("Counts /0.5 MeV ");
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();
  h1->Draw();
}
