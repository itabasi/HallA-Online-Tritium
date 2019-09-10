
void multi_foil(){

  TFile* f1=new TFile("/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/zcalib/zt_LHRS_132.root");
  TTree* t1=(TTree*)f1->Get("T");
  double Lz[100],Lz_c[100];
  t1->SetBranchStatus("*",0);
  t1->SetBranchStatus("L.tr.vz",1);
  t1->SetBranchAddress("L.tr.vz",Lz);
  t1->SetBranchStatus("L.tr.vz_tuned",1);
  t1->SetBranchAddress("L.tr.vz_tuned",Lz_c);

  //  int ENum=t1->GetEntires();
  TH1F* hLz=(TH1F*)f1->Get("h1");
  TH1F* hLz_c=(TH1F*)f1->Get("h2");  



  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  hLz->GetYaxis()->SetRangeUser(0,350);
  hLz->Draw();  

  //---- Multi Foil Line ----//
  TLine* l[10];
  for(int i=0;i<10;i++){
    l[i]=new TLine();
    l[i]->SetLineWidth(2);
    l[i]->SetLineColor(2);
    l[i]->SetLineStyle(2);
    if(i<8) l[i]->DrawLine(-0.125+i*0.025,0,-0.125+i*0.025,350);
    else    l[i]->DrawLine(-0.125+i*0.025+0.025,0,-0.125+i*0.025+0.025,350);
  }
   
 
  

  
  

}
