////////////////////////
// Empty Target study 
// Auther K. Itabashi
/////////////////////////

void empty(){

  string ifname="../rootfiles/optics/Al_target/empty_small.root";

  double Rz,Lz;
    TChain *T = new TChain("T");
    T->Add(ifname.c_str());
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("Rz",1);
    T->SetBranchAddress("Rz",&Rz);
    T->SetBranchStatus("Lz",1);
    T->SetBranchAddress("Lz",&Lz);
    
  string ofname="../rootfiles/optics/Al_target/empty_fit.root";
  TFile* ofp= new TFile(ofname.c_str(),"recreate");
  TTree* tnew=new TTree("T","Al events estimation");
  tnew=T->CloneTree(0);
  double Zm;
  tnew->Branch("Zm",&Zm,"Zm/D");
    
    double min_al1=-0.2;
    double max_al1=-0.0;
    double min_al2=0.0;
    double max_al2=0.2;
    int bin_al=400;
    double width=((double)bin_al)/(max_al2-min_al1);
    //    TH1D* hAl=new TH1D("hAl1","(ZR+RL)/2.0 hist",bin_al,min_al1,max_al1);
    TH1D* hAl=new TH1D("hAl","(ZR+RL)/2.0 hist",400,-0.2,0.2);
    TGraphErrors* gAl=new TGraphErrors();
    
    int ENum=T->GetEntries();
    cout<<"Events "<<ENum<<endl;
    int nn=1;
    int nENum=ENum/10.;
        cout<<"nEvents "<<nENum<<endl;
    for(int k=0;k<ENum;k++){
      bool z_cut=false;
      Zm=-1;
      T->GetEntry(k);

      if(fabs(Rz-Lz)<0.025)Zm=(Lz+Rz)/2.0;
      hAl->Fill(Zm);
      tnew->Fill();
      if(k>nn*nENum){
	double bg=hAl->Integral(hAl->GetXaxis()->FindBin(-0.1),hAl->GetXaxis()->FindBin(0.1)); 
	//	double nAl=hAl->Integral(hAl->GetXaxis()->FindBin(-0.2),hAl->GetXaxis()->FindBin(0.2)); 
	double nAl=hAl->Integral(hAl->GetXaxis()->FindBin(-0.2),hAl->GetXaxis()->FindBin(-0.125))
	  + hAl->Integral(hAl->GetXaxis()->FindBin(0.125),hAl->GetXaxis()->FindBin(0.2)); 
	gAl->SetPoint(nn-1,k,bg);
    cout<<"Backgrond events "<<nAl<<" ev "<<bg<<endl;    
	//	cout<<"Al ratio "<<nAl<<endl;
	nn++;
      }


      
    }

    cout<<"nn "<<nn<<endl;
    //======= Fitting Hist ============//

    //    TF1* fAl1=new TF1("fAl1","gausn(0)",min_al1,max_al1);
    //    TF1* fAl2=new TF1("fAl2","gausn(0)",min_al2,max_al2);
    //======= Exp func fit ============//
    min_al1=-0.12;
    max_al2= 0.12;
    TF1* fAl1=new TF1("fAl1","expo(0)",min_al1,max_al1);
    TF1* fAl2=new TF1("fAl2","expo(0)",min_al2,max_al2);
    fAl1->SetName("fAl1");
    fAl2->SetName("fAl2");
    fAl1->SetParameters(-2.14091e+01 , -2.16077e+02);
    fAl2->SetParameters( 2.50582e+01 ,  2.49322e+02);
    //=================================//

    hAl->Fit("fAl1","","QR",min_al1,max_al1);
    hAl->Fit("fAl2","","QR",min_al2,max_al2);    
    fAl1->SetLineColor(2);
    fAl1->SetFillColor(2);
    fAl1->SetFillStyle(3002);
    fAl2->SetLineColor(6);
    fAl2->SetFillColor(6);
    fAl2->SetFillStyle(3002);




    TCanvas* c0=new TCanvas("c0","c0");
    c0->cd();
    hAl->Draw();
    fAl1->Draw("same");
    fAl2->Draw("same");

    TCanvas* c1=new TCanvas("c1","c1");
    c1->cd();
    gAl->SetName("gAl");
    gAl->SetTitle("Al BG contaminate ; Total Events; # Al B.G.");
    gAl->SetMarkerStyle(3);
    gAl->SetMarkerColor(2);
    gAl->SetMarkerSize(1);
    gAl->Draw("AP");

    
    cout<<"=================="<<endl;
    cout<<"==== COMMENT======"<<endl;
    cout<<"=================="<<endl;

    double int1=fAl1->Integral(-0.1,0);
    double int2=fAl2->Integral(0.0,0.1);
    double int0 =hAl->Integral(hAl->GetXaxis()->FindBin(-0.1),hAl->GetXaxis()->FindBin(0.1));

    cout<<"int 0 "<<int0<<endl;
    cout<<"Int 1 "<<int1*width<<endl;
    cout<<"Int 2 "<<int2*width<<endl;
			       
    //    tnew->Write();
    //    hAl->Write();
    //    fAl1->Write();
    //    fAl2->Write();
    //    ofp->Close();


    
}
