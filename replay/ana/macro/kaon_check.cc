// check pion tail
// itabashi 2020/06/04

void kaon_check(){

  string ifname = "../rootfiles/mmass/ana_Lambda/mmcalib_new/nnL_small_Ole_all.root";

  TFile* ifp=new TFile(ifname.c_str());

  TChain* T=new TChain("T");
  T->Add(ifname.c_str());

  double ct;
  int pid_cut,z_cut,nrun;

  T->SetBranchStatus("*",0); 
  T->SetBranchStatus("pid_cut",1);
  T->SetBranchAddress("pid_cut",&pid_cut);
  T->SetBranchStatus("z_cut",1);
  T->SetBranchAddress("z_cut",&z_cut); 
  T->SetBranchStatus("ct",1);
  T->SetBranchAddress("ct",&ct);
  T->SetBranchStatus("runnum",1);
  T->SetBranchAddress("runnum",&nrun);

  int ENum = T->GetEntries();

  TH1F* hct=new TH1F("hct","",2000,-100,100);
  TH1F* hct_acc=new TH1F("hct_acc","",2000,-100,100);

  for(int i=0;i<ENum;i++){
    T->GetEntry(i);
    if(pid_cut>0 && z_cut>0){
      hct->Fill(ct);
      double ct_acc=ct;
      if(fabs( ct_acc )<80. && fabs(ct_acc)>20.){
	while(1){
	  if(ct_acc>20.)
	    ct_acc = ct_acc - 20.;
	  if(ct_acc<-20.)
	    ct_acc = ct_acc +20.;
	  if(fabs(ct_acc)<20.){
	    hct_acc ->Fill(ct_acc);
	    break;
	  }
	    
	}// while
	
      } // if ct_cut
    }// if cut

  }
    hct_acc->Scale(1./3.);
    TH1F* hct_peak = (TH1F*)hct->Clone();
    hct_peak->SetName("hct_peak");
    hct_peak->Add(hct_acc,-1.);
    TCanvas* c0=new TCanvas("c0","c0");
    c0->cd();
    
    //====== Pion Fitting ======//
    TF1* fpi_1= new TF1("fpi_1","gausn(0)",-5,10);
    TF1* fpi_2= new TF1("fpi_2","gausn(0)",-5,10);
    TF1* fpi= new TF1("fpi","gausn(0)+ gausn(3)",-5,10);
    hct_peak->Fit("fpi_1","","",2,4);
    double Ppi_1[3],Ppi_2[3],Ppi[3];

    Ppi_1[0] = fpi_1->GetParameter(0);
    Ppi_1[1] = fpi_1->GetParameter(1);
    Ppi_1[2] = fpi_1->GetParameter(2);

    fpi->FixParameter(0,Ppi_1[0]);
    fpi->FixParameter(1,Ppi_1[1]);
    fpi->FixParameter(2,Ppi_1[2]);


    fpi->SetParameter(3,500.);
    fpi->SetParameter(4,2.5);
    fpi->SetParameter(5,2.);
    hct_peak->Fit("fpi","","",-5,2);
    Ppi_1[0] = fpi->GetParameter(0);
    Ppi_1[1] = fpi->GetParameter(1);
    Ppi_1[2] = fpi->GetParameter(2);
    Ppi_2[0] = fpi->GetParameter(3);
    Ppi_2[1] = fpi->GetParameter(4);
    Ppi_2[2] = fpi->GetParameter(5);        


    fpi_1->SetParameter(0,Ppi_1[0]);
    fpi_1->SetParameter(1,Ppi_1[1]);
    fpi_1->SetParameter(2,Ppi_1[2]);
    fpi_2->SetParameter(0,Ppi_2[0]);
    fpi_2->SetParameter(1,Ppi_2[1]);
    fpi_2->SetParameter(2,Ppi_2[2]);        

    hct->GetXaxis()->SetRangeUser(-20,20);
    hct->SetLineColor(1);
    hct_peak->SetLineColor(4);
    hct->Draw("");
    hct_peak->GetXaxis()->SetRangeUser(-20,20);
    hct_peak->Draw("same");
    
    //    hct_acc->Draw("same");
    fpi->SetLineColor(2);
    fpi_1->SetLineColor(3);
    fpi_2->SetLineColor(4);
    fpi->Draw("same"); 
    fpi_1->Draw("same");
    fpi_2->Draw("same");


    TFile* ofp=new TFile("./pion.root","recreate");
    hct->Write();
    hct_acc->Write();
    hct_peak->Write();
  } //main
  

