/////////////////////////////////////////////////////
// In Coin time, Correction of Rp at FP dependence
// Auther Itabashi Nov. 21st, 2019
//////////////////////////////////////////////////////

void Rp_corr2(){
  
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
  int  nn=100;
  TH1D* hCoin[nn];//=new TH1D("hCoin","",1000,-20,20);
  TH2D* hCoin2[nn];
  for(int i=0;i<nn;i++){
    hCoin2[i]=new TH2D(Form("hCoin2_%d",i),"",1000,-20,30,20,1.6,2.0);  
    hCoin[i]=new TH1D(Form("hCoin_%d",i),"",1000,-20,30);  }
  int ENum=T->GetEntries();
  cout<<"Events : "<<ENum<<endl;

  for(int k=0;k<ENum;k++){
    ct=-1000.;
    Rp=-10.;
    z_cut=-1;
    trig=0;
    
    T->GetEntry(k);

    if(trig==5 && z_cut>0 && ac1<1.0 && ac2>1.25){
      for(int i=0;i<nn;i++){
	double ct_c=ct-Rp/( (double)(-nn + i*2.0)*0.05+0.3 ) +7.5;
	hCoin[i]->Fill(ct_c);
	hCoin2[i]->Fill(ct_c,Rp);

      }

      
      hcoin->Fill(ct,Rp);}
    }

  double pi[nn],mean[nn],sig[nn],sig_err[nn];
  TGraphErrors* gcoin=new TGraphErrors();
  TF1* fit[nn];
  for(int i=0;i<nn;i++){
    
    mean[i]=hCoin[i]->GetXaxis()->GetBinCenter(hCoin[i]->GetMaximumBin()) -3.0;    
    fit[i]=new TF1(Form("fit_%d",i),"gausn(0)",-5,5);
    fit[i]->SetLineColor(2);
    hCoin[i]->Fit(Form("fit_%d",i),"","QR",-1.0+mean[i], +1.0 +mean[i]);
    mean[i]= fit[i]->GetParameter(1);
    sig[i] = fit[i]->GetParameter(2);
    sig_err[i] = fit[i]->GetParError(2);
    gcoin->SetPoint(i,i,sig[i]);
    gcoin->SetPointError(i,0,sig_err[i]);
  }



  
  

  c0->cd();
  gcoin->SetMarkerColor(2);
  gcoin->SetMarkerStyle(20);
  
  gcoin->Draw("AP");
  //  hCoin[9]->Draw();

  TCanvas* c1=new TCanvas("c1","c1");
  c1->Divide(5,2);
    TCanvas* c2=new TCanvas("c2","c2");
  c2->Divide(5,2);
  for(int i=0;i<nn;i++){
    c1->cd(i+1);
    hCoin[i]->Draw();
    c2->cd(i+1);
    hCoin2[i]->Draw("colz");    
  }

  TCanvas* c3=new TCanvas("c3","c3");
  double mean_min=100;
  int min=0;
  for(int i=0;i<nn;i++){
    if(mean[i]<mean_min)mean_min=mean[i];
    min=i;
      }
  
  c3->cd();
  hCoin[min]->Draw();
}//end main
