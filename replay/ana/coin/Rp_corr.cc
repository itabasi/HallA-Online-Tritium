/////////////////////////////////////////////////////
// In Coin time, Correction of Rp at FP dependence
// Auther Itabashi Nov. 21st, 2019
//////////////////////////////////////////////////////

void Rp_corr(){
  
  string ifname="../rootfiles/mmass/ana_Lambda/Lambda_small_optH_1119.root";

  TFile* ifs=new TFile(ifname.c_str());
  TChain* T=new TChain("T");
  T->Add(ifname.c_str());

  double ct,Rp,ac1,ac2,trig;
  int z_cut;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("ct",1);
  T->SetBranchAddress("ct",&ct);
  T->SetBranchStatus("Rp",1);
  T->SetBranchAddress("Rp",&Rp);  
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

  TH2D*hcoin=new TH2D("hcoin","Coin vs Rx_FP;coin time [ns] ; Rx_fp [m]",20,-1,1,20,1.6,2.0);
  
  int ENum=T->GetEntries();
  cout<<"Events : "<<ENum<<endl;

  for(int k=0;k<ENum;k++){
    ct=-1000.;
    Rp=-10.;
    z_cut=-1;
    trig=0;
    
    T->GetEntry(k);
    //        cout<<"trig "<<trig<<" z_cut "<<z_cut<<" ac1 "<<ac1<<" ac2 "<<ac2<<" Rp "<<Rp<<" ct "<<ct<<endl;
    if(trig==5 && z_cut>0 && ac1<1.0 && ac2>1.25 && fabs(Rp-0.16*ct-1.85)<0.05)hcoin->Fill(ct,Rp);
    }




  TObjArray aSlices;
  TF1*g=new TF1("g","gausn(0)",1.7,1.9);
  hcoin->FitSlicesY(g,0,-1,0,"QRG2");
  TH1D*  hslice=(TH1D*)gROOT->FindObject("hcoin_1");

  TF1* fit=new TF1("fit","[0]*x+[1]",1.75,1.95);
  fit->SetLineColor(6);
  fit->SetParameter(0,1.0);
  hslice->Fit("fit","","QR");

  double p0=fit->GetParameter(0);
  double p1=fit->GetParameter(1);
  cout<<"p0 "<<p0<<" p1 "<<p1<<endl;
  p0=0.25;

  hcoin->Draw("colz");
  hslice->Draw("same");
  fit->Draw("same");

  
  //====================================//
  //====== Correction =================//
  //===================================//

  TH2D*hcoin_c=new TH2D("hcoin_c","Coin vs Rx_FP;coin time [ns] ; Rp [m]",1000,-20,20,20,1.75,1.95);
  TH2D*hcoin_b=new TH2D("hcoin_b","Coin vs Rx_FP;coin time [ns] ; Rp [m]",1000,-20,20,20,1.75,1.95);
  for(int k=0;k<ENum;k++){
    ct=-1000.;
    Rp=-10.;
    z_cut=-1;
    trig=0;
    T->GetEntry(k);
    //    Rp=Rp-p0*ct + p1;
   
    if(trig==5 && z_cut>0 && ac1<1.0 && ac2>1.25)hcoin_b->Fill(ct,Rp);
    //    ct=(ct-Rp-p1)/p0;
    ct=ct-Rp/p0;
    if(trig==5 && z_cut>0 && ac1<1.0 && ac2>1.25)hcoin_c->Fill(ct,Rp);
  }



  TCanvas* c1=new TCanvas("c1","c1");
  c1->Divide(1,2);
  c1->cd(1);    
  hcoin_b->Draw("colz");  
  c1->cd(2);    
  hcoin_c->Draw("colz");

}
