// S2 Offset tuning 

//#include "Setting.h"

double offset_fadc(int  RHRS, int RPMT, int Seg);


void offset(){


  //  Setting* set=new Setting;
  
  string ifname="/home/itabashi/ana/E12-17-003/HallA-Online-Tritium/replay/ana/rootfiles/tcoin/Lambda_small1.root";
  TFile* ft= new TFile(ifname.c_str());
  TTree* tree=(TTree*)ft->Get("T");

 double f1_s2r;//in real this is left PMT
 double f1_s2l;
 double f1_Rs2r[16],f1_Rs2l[16],f1_Ls2r[16],f1_Ls2l[16];
 double fbus_Rs2r[16],fbus_Rs2l[16],fbus_Ls2r[16],fbus_Ls2l[16];
 double fadc_Rs2r[16],fadc_Rs2l[16],fadc_Ls2r[16],fadc_Ls2l[16];
 double coin_f1,coin_fadc,coin_fbus;
 double coin_f1_c,coin_fadc_c,coin_fbus_c;
 int event,trig;
 int Ls2pads,Rs2pads;
   tree->SetBranchAddress("nev",&event);
   tree->SetBranchAddress("trig",&trig);
   tree->SetBranchAddress("Rs2_pads",&Rs2pads);
   tree->SetBranchAddress("Ls2_pads",&Ls2pads);
   
 //-------- F1TDC -----------//
  tree->SetBranchAddress("Ls2_f1_rt", f1_Ls2r);
  tree->SetBranchAddress("Ls2_f1_lt", f1_Ls2l);
  tree->SetBranchAddress("Rs2_f1_rt", f1_Rs2r);
  tree->SetBranchAddress("Rs2_f1_lt", f1_Rs2l);
  tree->SetBranchAddress("coin_f1", &coin_f1);
  tree->SetBranchAddress("coin_f1_c", &coin_f1_c);  
  //------ FBUS -------------//
  tree->SetBranchAddress("Ls2_fbus_rt", fbus_Ls2r);
  tree->SetBranchAddress("Ls2_fbus_lt", fbus_Ls2l);
  tree->SetBranchAddress("Rs2_fbus_rt", fbus_Rs2r);
  tree->SetBranchAddress("Rs2_fbus_lt", fbus_Rs2l);
  tree->SetBranchAddress("coin_fbus", &coin_fbus);
  tree->SetBranchAddress("coin_fbus_c", &coin_fbus_c);   
  //------- FADC -------------//
  tree->SetBranchAddress("Ls2_fadc_rt", fadc_Ls2r);
  tree->SetBranchAddress("Ls2_fadc_lt", fadc_Ls2l);
  tree->SetBranchAddress("Rs2_fadc_rt", fadc_Rs2r);
  tree->SetBranchAddress("Rs2_fadc_lt", fadc_Rs2l);
  tree->SetBranchAddress("coin_fadc", &coin_fadc);
  tree->SetBranchAddress("coin_fadc_c", &coin_fadc_c);


  
  TH1F* hRs2r_fadc[16];
  TH1F* hRs2l_fadc[16];
  TH1F* hLs2r_fadc[16];
  TH1F* hLs2l_fadc[16];
  TH1F* hRs2r_fadc_c[16];
  TH1F* hRs2l_fadc_c[16];
  TH1F* hLs2r_fadc_c[16];
  TH1F* hLs2l_fadc_c[16];  

  double min_fadc_r=240.0;
  double max_fadc_r=260.0;
  double min_fadc_l=130.0;
  double max_fadc_l=150.0;
  int bin_fadc_r=(int)(max_fadc_r-min_fadc_r)/0.0625;
  int bin_fadc_l=(int)(max_fadc_l-min_fadc_l)/0.0625;
  char temp0[500],temp1[500],temp2[500],temp3[500];
  char temp0_c[500],temp1_c[500],temp2_c[500],temp3_c[500];
  for(int i=0;i<16;i++){
    sprintf(temp0,"hRs2r_fadc_%d",i);
    hRs2r_fadc[i]=new TH1F(temp0,"",bin_fadc_r,min_fadc_r,max_fadc_r);
    hRs2r_fadc[i]->SetTitle(Form("RHRS S2 Seg %d R-PMT FADC TDC Hist; TDC [ns];Counts/62.5 ps",i));
    sprintf(temp1,"hRs2l_fadc_%d",i);
    hRs2l_fadc[i]=new TH1F(temp1,"",bin_fadc_r,min_fadc_r,max_fadc_r);    
    hRs2l_fadc[i]->SetTitle(Form("RHRS S2 Seg %d L-PMT FADC TDC Hist; TDC [ns];Counts/62.5 ps",i));
    sprintf(temp2,"hLs2r_fadc_%d",i);
    hLs2r_fadc[i]=new TH1F(temp2,"",bin_fadc_l,min_fadc_l,max_fadc_l);
    hLs2r_fadc[i]->SetTitle(Form("LHRS S2 Seg %d R-PMT FADC TDC Hist;TDC [ns]; Counts/62.5ps",i));
    sprintf(temp3,"hLs2l_fadc_%d",i);
    hLs2l_fadc[i]=new TH1F(temp3,"",bin_fadc_l,min_fadc_l,max_fadc_l);    
    hLs2l_fadc[i]->SetTitle(Form("LHRS S2 Seg %d L-PMT FADC TDC Hist;TDC [ns];Counts/62.5 ps",i));


    sprintf(temp0_c,"hRs2r_fadc_c%d",i);
    hRs2r_fadc_c[i]=new TH1F(temp0_c,"",bin_fadc_r,min_fadc_r,max_fadc_r);
    hRs2r_fadc_c[i]->SetTitle(Form("RHRS S2 Seg %d R-PMT FADC TDC Hist; TDC [ns];Counts/62.5 ps",i));
    sprintf(temp1_c,"hRs2l_fadc_c%d",i);
    hRs2l_fadc_c[i]=new TH1F(temp1_c,"",bin_fadc_r,min_fadc_r,max_fadc_r);    
    hRs2l_fadc_c[i]->SetTitle(Form("RHRS S2 Seg %d L-PMT FADC TDC Hist; TDC [ns];Counts/62.5 ps",i));
    sprintf(temp2_c,"hLs2r_fadc_c%d",i);
    hLs2r_fadc_c[i]=new TH1F(temp2_c,"",bin_fadc_l,min_fadc_l,max_fadc_l);
    hLs2r_fadc_c[i]->SetTitle(Form("LHRS S2 Seg %d R-PMT FADC TDC Hist;TDC [ns]; Counts/62.5ps",i));
    sprintf(temp3_c,"hLs2l_fadc_c%d",i);
    hLs2l_fadc_c[i]=new TH1F(temp3_c,"",bin_fadc_l,min_fadc_l,max_fadc_l);    
    hLs2l_fadc_c[i]->SetTitle(Form("LHRS S2 Seg %d L-PMT FADC TDC Hist;TDC [ns];Counts/62.5 ps",i));    
  }


    bool test =true;
			     //    test=false;
    int nev=tree->GetEntries();
    if(test)nev=100000;
    cout<<"Events: "<<nev<<endl;

    for(int k=0;k<nev;k++){
      tree->GetEntry(k);
      if(trig==5){
      for(int i=0;i<16;i++){
	hRs2r_fadc[i]->Fill(fadc_Rs2r[i]);
	hRs2l_fadc[i]->Fill(fadc_Rs2l[i]);	
	hLs2r_fadc[i]->Fill(fadc_Ls2r[i]);
	hLs2l_fadc[i]->Fill(fadc_Ls2l[i]);

	hRs2r_fadc_c[i]->Fill(fadc_Rs2r[i]-offset_fadc(1,1,i));
	hRs2l_fadc_c[i]->Fill(fadc_Rs2l[i]-offset_fadc(1,0,i));	
	hLs2r_fadc_c[i]->Fill(fadc_Ls2r[i]-offset_fadc(0,1,i));
	hLs2l_fadc_c[i]->Fill(fadc_Ls2l[i]-offset_fadc(0,0,i));
	
      }
      }
      if(k % 100000 ==0)cout<<k<<" / "<<nev<<endl;
    }

    double Rs2r_p[16],Rs2l_p[16],Ls2r_p[16],Ls2l_p[16];
          for(int i=0;i<16;i++){
	Rs2r_p[i]= hRs2r_fadc[i]-> GetXaxis()-> GetBinCenter(hRs2r_fadc[i]->GetMaximumBin());
	Rs2l_p[i]= hRs2l_fadc[i]-> GetXaxis()-> GetBinCenter(hRs2l_fadc[i]->GetMaximumBin());
	Ls2r_p[i]= hLs2r_fadc[i]-> GetXaxis()-> GetBinCenter(hLs2r_fadc[i]->GetMaximumBin());
	Ls2l_p[i]= hLs2l_fadc[i]-> GetXaxis()-> GetBinCenter(hLs2l_fadc[i]->GetMaximumBin());
	  }
  





	  //======================================//
	  //========= COMMENT OUT ================//
	  //======================================//
	  
	  cout<<"======= Rs2r_offset ========= "<<endl;
	  for(int i=0;i<16;i++){
	  if(Rs2r_p[i]-Rs2r_p[6]>1.0) Rs2r_p[i]=Rs2r_p[i]-4.0;
          else if(Rs2r_p[i]-Rs2r_p[6]<-1.0) Rs2r_p[i]=Rs2r_p[i]+4.0;
	  cout<<Rs2r_p[i]-Rs2r_p[6]<<", ";}
	  cout<<endl;
	  cout<<"======= Rs2l_offset ========= "<<endl;
	  for(int i=0;i<16;i++){
	if(Rs2l_p[i]-Rs2l_p[5]>1.0)Rs2l_p[i]=Rs2l_p[i]-4.0;
	else if(Rs2l_p[i]-Rs2l_p[5]<-1.0)Rs2l_p[i]=Rs2l_p[i]+4.0;
	  cout<<Rs2l_p[i]-Rs2l_p[5]<<", ";	  }
	  cout<<endl;
	  cout<<"======= Ls2r_offset ========= "<<endl;
	  for(int i=0;i<16;i++)cout<<Form("Ls2r_p%d : ",i)<<Ls2r_p[i]<<endl;
	  for(int i=0;i<16;i++){
	  if(Ls2r_p[i]-Ls2r_p[6]>1.0) Ls2r_p[i]=Ls2r_p[i]-4.0;
          else if(Ls2r_p[i]-Ls2r_p[6]<-1.0) Ls2r_p[i]=Ls2r_p[i]+4.0;
	  cout<<Ls2r_p[i]-Ls2r_p[6]<<", ";}
	  cout<<endl;
	  cout<<"======= Ls2l_offset ========= "<<endl;
	  for(int i=0;i<16;i++)cout<<Form("Ls2l_p%d : ",i)<<Ls2l_p[i]<<endl;
	  for(int i=0;i<16;i++){
	  if(Ls2l_p[i]-Ls2l_p[6]>1.0)Ls2l_p[i]=Ls2l_p[i]-4.0;
	  else if(Ls2l_p[i]-Ls2l_p[6]<-1.0)Ls2l_p[i]=Ls2l_p[i]+4.0;
	    cout<<Ls2l_p[i]-Ls2l_p[6]<<", ";	  }
	  cout<<endl;






         //=======================//
         //======== Fit ==========//
         //=======================//


        TF1* fit_Rs2r[16];
        TF1* fit_Rs2l[16];
        TF1* fit_Ls2r[16];
	TF1* fit_Ls2l[16];
	TGraph* gRs2r_off=new TGraph(); gRs2r_off->SetTitle("RHRS S2 RPMT offset");
	gRs2r_off->SetMarkerStyle(21);
	gRs2r_off->SetMarkerColor(kRed);
	gRs2r_off->SetMarkerSize(0.5);
	gRs2r_off->SetMinimum(Rs2r_p[6]-1.0); 	gRs2r_off->SetMaximum(Rs2r_p[6]+1.0);
	TGraph* gRs2l_off=new TGraph(); gRs2l_off->SetTitle("RHRS S2 LPMT offset");
	gRs2l_off->SetMarkerStyle(21);
	gRs2l_off->SetMarkerColor(kRed);
	gRs2l_off->SetMarkerSize(0.5);
 	gRs2l_off->SetMinimum(Rs2l_p[5]-1.0); 	gRs2l_off->SetMaximum(Rs2l_p[5]+1.0);
	TGraph* gLs2r_off=new TGraph(); gLs2r_off->SetTitle("LHRS S2 RPMT offset");
	gLs2r_off->SetMarkerStyle(21);
	gLs2r_off->SetMarkerColor(kRed);
	gLs2r_off->SetMarkerSize(0.5);
	gLs2r_off->SetMinimum(Ls2r_p[6]-1.0); 	gLs2r_off->SetMaximum(Ls2r_p[6]+1.0);
	TGraph* gLs2l_off=new TGraph(); gLs2l_off->SetTitle("LHRS S2 LPMT offset");
	gLs2l_off->SetMarkerStyle(21);
	gLs2l_off->SetMarkerColor(kRed);
	gLs2l_off->SetMarkerSize(0.5);
	gLs2l_off->SetMinimum(Ls2l_p[6]-1.0); 	gLs2l_off->SetMaximum(Ls2l_p[6]+1.0);

	
	double Rs2r_p1[16],Rs2l_p1[16],Ls2r_p1[16],Ls2l_p1[16];
	for(int i=0;i<16;i++){
	  fit_Rs2r[i]=new TF1(Form("fit_Rs2r[%d]",i),"gausn(0)",min_fadc_r,max_fadc_r);
	  fit_Rs2r[i]->SetLineColor(2);
	  hRs2r_fadc_c[i]->Fit(fit_Rs2r[i],"Rq","",Rs2r_p[6]-0.5,Rs2r_p[6]+0.5);
	  Rs2r_p1[i]=fit_Rs2r[i]->GetParameter(1);
	  gRs2r_off->SetPoint(i,i,Rs2r_p1[i]);
	  fit_Rs2l[i]=new TF1(Form("fit_Rs2l[%d]",i),"gausn(0)",min_fadc_r,max_fadc_r);
	  fit_Rs2l[i]->SetLineColor(2);
	  hRs2l_fadc_c[i]->Fit(fit_Rs2l[i],"Rq","",Rs2l_p[5]-0.5,Rs2l_p[5]+0.5);
	  Rs2l_p1[i]=fit_Rs2l[i]->GetParameter(1);
	  gRs2l_off->SetPoint(i,i,Rs2l_p1[i]);
	  fit_Ls2r[i]=new TF1(Form("fit_Ls2r[%d]",i),"gausn(0)",min_fadc_l,max_fadc_l);
	  fit_Ls2r[i]->SetLineColor(2);
	  hLs2r_fadc_c[i]->Fit(fit_Ls2r[i],"Rq","",Ls2r_p[6]-0.5,Ls2r_p[6]+0.5);	  
	  Ls2r_p1[i]=fit_Ls2r[i]->GetParameter(1);
	  gLs2r_off->SetPoint(i,i,Ls2r_p1[i]);
	  fit_Ls2l[i]=new TF1(Form("fit_Ls2l[%d]",i),"gausn(0)",min_fadc_l,max_fadc_l);
     	  fit_Ls2l[i]->SetLineColor(2);
	  hLs2l_fadc_c[i]->Fit(fit_Ls2l[i],"Rq","",Ls2l_p[6]-0.5,Ls2l_p[6]+0.5);	  
	  Ls2l_p1[i]=fit_Ls2l[i]->GetParameter(1);
	  gLs2l_off->SetPoint(i,i,Ls2l_p1[i]);
	}
	

	
	
        //=============================//
        //========== Draw =============//
        //=============================//


	  gROOT->SetBatch();
	  TCanvas* c0=new TCanvas("c0","RHRS S2 RPMT ");
	  c0->Divide(4,4);
	  TCanvas* c1=new TCanvas("c1","RHRS S2 LPMT ");
	  c1->Divide(4,4);
	  TCanvas* c2=new TCanvas("c2","LHRS S2 RPMT ");
	  c2->Divide(4,4);
	  TCanvas* c3=new TCanvas("c3","LHRS S2 LPMT ");
	  c3->Divide(4,4);
	  TCanvas* c4=new TCanvas("c4","Offset position Graph");
	  c4->Divide(2,2);
	  
	  for(int i=0;i<16;i++){
	    c0->cd(i+1);
	    hRs2r_fadc[i]->SetLineColor(1);
	    hRs2r_fadc[i]->Draw();
	    hRs2r_fadc_c[i]->SetLineColor(3);
	    hRs2r_fadc_c[i]->Draw("same"); 
	    c1->cd(i+1);
	    hRs2l_fadc[i]->SetLineColor(1);	    
	    hRs2l_fadc[i]->Draw();
	    hRs2l_fadc_c[i]->SetLineColor(3);
	    hRs2l_fadc_c[i]->Draw("same"); 	    
            c2->cd(i+1);
	    hLs2r_fadc[i]->SetLineColor(1);	    
	    hLs2r_fadc[i]->Draw();
	    hLs2r_fadc_c[i]->SetLineColor(3);
	    hLs2r_fadc_c[i]->Draw("same"); 	    
	    c3->cd(i+1);
	    hLs2l_fadc[i]->SetLineColor(1);	    
	    hLs2l_fadc[i]->Draw();
	    hLs2l_fadc_c[i]->SetLineColor(3);
	    hLs2l_fadc_c[i]->Draw("same"); 	    
	  }

	  c4->cd(1);
	  gRs2r_off->Draw("AP");
	  c4->cd(2);
	  gRs2l_off->Draw("AP");	  
	  c4->cd(3);
	  gLs2r_off->Draw("AP");	  
	  c4->cd(4);
	  gLs2l_off->Draw("AP");

	  
	  string ofname="./offset.pdf";
	  TString name;
	  name.Form(ofname.c_str());
	  c0->Print(name + "[","pdf");
	  c0->Print(name,"pdf");	  
	  c1->Print(name,"pdf");
	  c2->Print(name,"pdf");
	  c3->Print(name,"pdf");	  
	  c4->Print(name,"pdf");	  
	  c4->Print(name +"]","pdf");	  


}


double offset_fadc(int RHRS,int  RPMT, int Seg){

  double Rs2r_offset[16]={-0.5625, -0.125, -0.0625, 0, -0.0625, -0.0625, 0, -0.0625, -0.0625, 0, 0, 0.0625, 0.0625, 0, -0.0625, 0};
  //   {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double Rs2l_offset[16]={0.0625, 0.0625, 0, 0, 0, 0, 4.0625, 0.0625, 0, 0.0625, -0.0625, 4, 3.9375, 4, -0.0625, -0.0625};
   //       {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  
   double Ls2r_offset[16]={3.9375, 4, 0, 0, 0, 0.0625, 0, 4.0625, 0.0625, 0.0625, 0, 0, 0, 0, 0.0625, 0.0625};
   //     {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
   double Ls2l_offset[16]={0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0.0625, 0, 0.0625, 4, 0};
   //     {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};


  if(RHRS==1 && RPMT==1)          return Rs2r_offset[Seg];
  else if (RHRS==1 && RPMT==0)    return Rs2l_offset[Seg];
  else if (RHRS==0 && RPMT==1)    return Ls2r_offset[Seg];
  else if (RHRS==0 && RPMT==0)    return Ls2l_offset[Seg];
 
}
