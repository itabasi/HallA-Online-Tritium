//////////////////////////////////////////////////////
// In Coin time, Correction of Rx at FP dependence
// Auther Itabashi Nov. 21st, 2019
//////////////////////////////////////////////////////

void Rx_corr(){
  
  string ifname="../rootfiles/mmass/ana_Lambda/Lambda_small_optH_1119.root";

  TFile* ifs=new TFile(ifname.c_str());
  TChain* T=new TChain("T");
  T->Add(ifname.c_str());

  double ct,Rx_fp,ac1,ac2,trig;
  int z_cut;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("ct",1);
  T->SetBranchAddress("ct",&ct);
  T->SetBranchStatus("Rx_fp",1);
  T->SetBranchAddress("Rx_fp",&Rx_fp);  
  T->SetBranchStatus("Rx_fp",1);
  T->SetBranchAddress("Rx_fp",&Rx_fp);
  T->SetBranchStatus("trig",1);
  T->SetBranchAddress("trig",&trig);
  T->SetBranchStatus("z_cut",1);
  T->SetBranchAddress("z_cut",&z_cut);    
  T->SetBranchStatus("ac1_npe_sum",1);
  T->SetBranchAddress("ac1_npe_sum",&ac1);
  T->SetBranchStatus("ac2_npe_sum",1);
  T->SetBranchAddress("ac2_npe_sum",&ac2);    


  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();  

  TH2D*hcoin=new TH2D("hcoin","Coin vs Rx_FP;coin time [ns] ; Rx_fp [m]",20,-1,1,20,-0.8,0.8);
  
  int ENum=T->GetEntries();
  cout<<"Events : "<<ENum<<endl;

  for(int k=0;k<ENum;k++){
    ct=-1000.;
    Rx_fp=-10.;
    z_cut=-1;
    trig=0;
    
    T->GetEntry(k);
    //    cout<<"trig "<<trig<<" z_cut "<<z_cut<<" ac1 "<<ac1<<" ac2 "<<ac2<<" Rx_fp "<<Rx_fp<<" ct "<<ct<<endl;
    if(trig==5 && z_cut>0 && ac1<1.0 && ac2>1.25 && fabs(ct-Rx_fp)<0.5)hcoin->Fill(ct,Rx_fp);
    }




  TObjArray aSlices;
  TF1*g=new TF1("g","gausn(0)",-1,1);
  hcoin->FitSlicesY(g,0,-1,0,"QRG2");
  TH1D*  hslice=(TH1D*)gROOT->FindObject("hcoin_1");

  TF1* fit=new TF1("fit","[0]*x+[1]",-0.6,0.6);
  fit->SetLineColor(6);
  fit->SetParameter(0,1.0);
  hslice->Fit("fit","","QR");

  double p0=fit->GetParameter(0);
  double p1=fit->GetParameter(1);
  cout<<"p0 "<<p0<<" p1 "<<p1<<endl;


  hcoin->Draw("colz");
  hslice->Draw("same");
  fit->Draw("same");

  
  //====================================//
  //====== Correction =================//
  //===================================//

  TH2D*hcoin_c=new TH2D("hcoin_c","Coin vs Rx_FP;coin time [ns] ; Rx_fp [m]",1000,-20,20,20,-1,1);
  TH2D*hcoin_b=new TH2D("hcoin_b","Coin vs Rx_FP;coin time [ns] ; Rx_fp [m]",1000,-20,20,20,-1,1);
  for(int k=0;k<ENum;k++){
    ct=-1000.;
    Rx_fp=-10.;
    z_cut=-1;
    trig=0;
    T->GetEntry(k);
    //    Rx_fp=Rx_fp-p0*ct + p1;
   
    if(trig==5 && z_cut>0 && ac1<1.0 && ac2>1.25)hcoin_b->Fill(ct,Rx_fp);
    ct=(ct-Rx_fp-p1)/p0;
    if(trig==5 && z_cut>0 && ac1<1.0 && ac2>1.25)hcoin_c->Fill(ct,Rx_fp);
  }



  TCanvas* c1=new TCanvas("c1","c1");
  c1->Divide(1,2);
  c1->cd(1);    
  hcoin_b->Draw("colz");  
  c1->cd(2);    
  hcoin_c->Draw("colz");

}
