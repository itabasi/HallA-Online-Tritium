//Auther Itabashi 2019/10/02
// VDC wire efficiency study 

void vdc_eff(){

  string fname="/data1/root/tritium_111500.root";
  TChain* T=new TChain("T");
  //    T->Add(fname.c_str());
  for(int i=111500;i<111510;i++){  
    string name = Form("/data1/root/tritium_%d.root",i);
    T->Add(name.c_str());
  }
  int nmax=500;
  double Ru1_wire[nmax],Ru2_wire[nmax],Rv1_wire[nmax],Rv2_wire[nmax],Lu1_wire[nmax],Lu2_wire[nmax],Lv1_wire[nmax],Lv2_wire[nmax];
  double Ru1_nclust,Ru2_nclust,Rv1_nclust,Rv2_nclust,Lu1_nclust,Lu2_nclust,Lv1_nclust,Lv2_nclust;
  int NRu1_wire,NRu2_wire,NRv1_wire,NRv2_wire, NLu1_wire,NLu2_wire,NLv1_wire,NLv2_wire;
  double evtype;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("DR.evtypebits"          ,1);            T->SetBranchAddress("DR.evtypebits"          ,&evtype);
  T->SetBranchStatus("R.vdc.u1.wire"          ,1);            T->SetBranchAddress("R.vdc.u1.wire"          ,Ru1_wire);
  T->SetBranchStatus("R.vdc.u2.wire"          ,1);            T->SetBranchAddress("R.vdc.u2.wire"          ,Ru2_wire);
  T->SetBranchStatus("R.vdc.v1.wire"          ,1);            T->SetBranchAddress("R.vdc.v1.wire"          ,Rv1_wire);
  T->SetBranchStatus("R.vdc.v2.wire"          ,1);            T->SetBranchAddress("R.vdc.v2.wire"          ,Rv2_wire);
  T->SetBranchStatus("Ndata.R.vdc.u1.wire"    ,1);            T->SetBranchAddress("Ndata.R.vdc.u1.wire"    ,&NRu1_wire);
  T->SetBranchStatus("Ndata.R.vdc.u2.wire"    ,1);            T->SetBranchAddress("Ndata.R.vdc.u2.wire"    ,&NRu2_wire);
  T->SetBranchStatus("Ndata.R.vdc.v1.wire"    ,1);            T->SetBranchAddress("Ndata.R.vdc.v1.wire"    ,&NRv1_wire);
  T->SetBranchStatus("Ndata.R.vdc.v2.wire"    ,1);            T->SetBranchAddress("Ndata.R.vdc.v2.wire"    ,&NRv2_wire);
  T->SetBranchStatus("R.vdc.u1.nclust"        ,1);            T->SetBranchAddress("R.vdc.u1.nclust"        ,&Ru1_nclust);
  T->SetBranchStatus("R.vdc.u2.nclust"        ,1);            T->SetBranchAddress("R.vdc.u2.nclust"        ,&Ru2_nclust);
  T->SetBranchStatus("R.vdc.v1.nclust"        ,1);            T->SetBranchAddress("R.vdc.v1.nclust"        ,&Rv1_nclust);
  T->SetBranchStatus("R.vdc.v2.nclust"        ,1);            T->SetBranchAddress("R.vdc.v2.nclust"        ,&Rv2_nclust);
  T->SetBranchStatus("L.vdc.u1.wire"          ,1);            T->SetBranchAddress("L.vdc.u1.wire"          ,Lu1_wire);
  T->SetBranchStatus("L.vdc.u2.wire"          ,1);            T->SetBranchAddress("L.vdc.u2.wire"          ,Lu2_wire);
  T->SetBranchStatus("L.vdc.v1.wire"          ,1);            T->SetBranchAddress("L.vdc.v1.wire"          ,Lv1_wire);
  T->SetBranchStatus("L.vdc.v2.wire"          ,1);            T->SetBranchAddress("L.vdc.v2.wire"          ,Lv2_wire);
  T->SetBranchStatus("L.vdc.u1.nclust"        ,1);            T->SetBranchAddress("L.vdc.u1.nclust"        ,&Lu1_nclust);
  T->SetBranchStatus("L.vdc.u2.nclust"        ,1);            T->SetBranchAddress("L.vdc.u2.nclust"        ,&Lu2_nclust);
  T->SetBranchStatus("L.vdc.v1.nclust"        ,1);            T->SetBranchAddress("L.vdc.v1.nclust"        ,&Lv1_nclust);
  T->SetBranchStatus("L.vdc.v2.nclust"        ,1);            T->SetBranchAddress("L.vdc.v2.nclust"        ,&Lv2_nclust);
  T->SetBranchStatus("Ndata.L.vdc.u1.wire"    ,1);            T->SetBranchAddress("Ndata.L.vdc.u1.wire"    ,&NLu1_wire);
  T->SetBranchStatus("Ndata.L.vdc.u2.wire"    ,1);            T->SetBranchAddress("Ndata.L.vdc.u2.wire"    ,&NLu2_wire);
  T->SetBranchStatus("Ndata.L.vdc.v1.wire"    ,1);            T->SetBranchAddress("Ndata.L.vdc.v1.wire"    ,&NLv1_wire);
  T->SetBranchStatus("Ndata.L.vdc.v2.wire"    ,1);            T->SetBranchAddress("Ndata.L.vdc.v2.wire"    ,&NLv2_wire);

  int nwire=400;
  int bin_wire=nwire;
  int min_wire=0;
  int max_wire=nwire;
  TH1F* hRu1_eff=new TH1F("hRu1_eff","RHRS VDC U1 wire efficienecy ; wires ; Efficiency [%]",bin_wire,min_wire,max_wire);
  TH1F* hRv1_eff=new TH1F("hRv1_eff","RHRS VDC V1 wire efficienecy ; wires ; Efficiency [%]",bin_wire,min_wire,max_wire);
  TH1F* hRu2_eff=new TH1F("hRu2_eff","RHRS VDC U2 wire efficienecy ; wires ; Efficiency [%]",bin_wire,min_wire,max_wire);
  TH1F* hRv2_eff=new TH1F("hRv2_eff","RHRS VDC V2 wire efficienecy ; wires ; Efficiency [%]",bin_wire,min_wire,max_wire);
  TH1F* hLu1_eff=new TH1F("hLu1_eff","LHRS VDC U1 wire efficienecy ; wires ; Efficiency [%]",bin_wire,min_wire,max_wire);
  TH1F* hLv1_eff=new TH1F("hLv1_eff","LHRS VDC V1 wire efficienecy ; wires ; Efficiency [%]",bin_wire,min_wire,max_wire);
  TH1F* hLu2_eff=new TH1F("hLu2_eff","LHRS VDC U2 wire efficienecy ; wires ; Efficiency [%]",bin_wire,min_wire,max_wire);
  TH1F* hLv2_eff=new TH1F("hLv2_eff","LHRS VDC V2 wire efficienecy ; wires ; Efficiency [%]",bin_wire,min_wire,max_wire);

  TGraphErrors* gRu1_eff=new TGraphErrors();
  //  gRu1_eff->SetTitle("RHRS VDC U1 Efficiency");
  gRu1_eff->SetTitle("RHRS VDC U1 Efficiency ; wires ; Efficiency");
  TGraphErrors* gRv1_eff=new TGraphErrors();
  gRv1_eff->SetTitle("RHRS VDC V1 Efficiency ; wires ; Efficiency");
  TGraphErrors* gRu2_eff=new TGraphErrors();
  gRu2_eff->SetTitle("RHRS VDC U2 Efficiency ; wires ; Efficiency");  
  TGraphErrors* gRv2_eff=new TGraphErrors();
  gRv2_eff->SetTitle("RHRS VDC V2 Efficiency ; wires ; Efficiency");    
  TGraphErrors* gLu1_eff=new TGraphErrors();
  gLu1_eff->SetTitle("LHRS VDC U1 Efficiency ; wires ; Efficiency");  
  TGraphErrors* gLv1_eff=new TGraphErrors();
  gLv1_eff->SetTitle("LHRS VDC V1 Efficiency ; wires ; Efficiency");    
  TGraphErrors* gLu2_eff=new TGraphErrors();
  gLu2_eff->SetTitle("LHRS VDC U2 Efficiency ; wires ; Efficiency");  
  TGraphErrors* gLv2_eff=new TGraphErrors();
  gLv2_eff->SetTitle("LHRS VDC V2 Efficiency ; wires ; Efficiency");  
  
  
  double Ru1[nwire],Ru2[nwire],Rv1[nwire],Rv2[nwire],Lu1[nwire],Lu2[nwire],Lv1[nwire],Lv2[nwire];
  double Ru1_err[nwire],Ru2_err[nwire],Rv1_err[nwire],Rv2_err[nwire],Lu1_err[nwire],Lu2_err[nwire],Lv1_err[nwire],Lv2_err[nwire];

  //==== Initialization ===//
  for(int i=0;i<nwire;i++){
    Ru1[i]=0.0; Ru1_err[i]=0.0;
    Ru2[i]=0.0; Ru2_err[i]=0.0;
    Rv1[i]=0.0; Rv1_err[i]=0.0;
    Rv2[i]=0.0; Rv2_err[i]=0.0;
    Lu1[i]=0.0; Lu1_err[i]=0.0;
    Lu2[i]=0.0; Lu2_err[i]=0.0;
    Lv1[i]=0.0; Lv1_err[i]=0.0;
    Lv2[i]=0.0; Lv2_err[i]=0.0;    
  }
  
  int ENum=T->GetEntries();
  cout<<"Events : "<<ENum<<endl;

  for(int k=0;k<ENum;k++){

    for(int i=0;i<nmax;i++){
      Ru1_wire[i]=-1.0;
      Ru2_wire[i]=-1.0;
      Rv1_wire[i]=-1.0;
      Rv2_wire[i]=-1.0;
      Lu1_wire[i]=-1.0;
      Lu2_wire[i]=-1.0;
      Lv1_wire[i]=-1.0;
      Lv2_wire[i]=-1.0;
      
    }
    
    T->GetEntry(k);

    bool trig=false;
    if(evtype==32)trig=true;
    trig=    true;
    //======= RU1 ======================//

    if(Ru1_nclust==0 && trig){
	if(Ru1_wire[1]-Ru1_wire[0]!=1.0 && Ru1_wire[1]-Ru1_wire[0]<=3.0)
	  for(int j=(int)Ru1_wire[0]+1;j<Ru1_wire[1];j++){
	    Ru1_err[j] +=1.0;
	  }

      for(int i=0;i<NRu1_wire;i++){
	Ru1[(int)Ru1_wire[i]] += 1.0;

      }
    }

    /*
    if(Ru1_nclust==0 && trig){
      for(int i=0;i<NRu1_wire;i++){
	Ru1[(int)Ru1_wire[i]] += 1.0;
	if(Ru1_wire[i+1]-Ru1_wire[i]!=1.0 && Ru1_wire[i+1]-Ru1_wire[i]<=3.0 && i<NRu1_wire-1)
	  for(int j=(int)Ru1_wire[i]+1;j<Ru1_wire[i+1];j++){
	    Ru1_err[j] +=1.0;
	    //	    cout<<"NRu1_wire "<<NRu1_wire<<" Ru1_wire "<<Ru1_wire[i]<<" Ru1 "<<Ru1[(int)Ru1_wire[i]]<<" j "<<j<<" Ru1_err "<<Ru1_err[j]<<endl;
	  }
      }
    }
*/
    //======= RU2 ======================//
    if(Ru2_nclust==0 && trig){
      for(int i=0;i<NRu2_wire;i++){
	Ru2[(int)Ru2_wire[i]] += 1.0;
	if(Ru2_wire[i+1]-Ru2_wire[i]!=1.0 && Ru2_wire[i+1]-Ru2_wire[i]<=3.0 && Ru2_nclust==0 && i<NRu2_wire-1)
	  for(int j=(int)Ru2_wire[i]+1;j<Ru2_wire[i+1];j++)Ru2_err[j] +=1.0;
      }    
    }
    //======= Rv1 ======================//
    if(Rv1_nclust==0 && trig){
      for(int i=0;i<NRv1_wire;i++){
	Rv1[(int)Rv1_wire[i]] += 1.0;
	if(Rv1_wire[i+1]-Rv1_wire[i]!=1.0 && Rv1_wire[i+1]-Rv1_wire[i]<=3.0 && Rv1_nclust==0 && i<NRv1_wire-1)
	  for(int j=(int)Rv1_wire[i]+1;j<Rv1_wire[i+1];j++)Rv1_err[j] +=1.0;
      }    
    }
    //======= Rv2 ======================//
    if(Rv2_nclust==0 && trig){
      for(int i=0;i<NRv2_wire;i++){
	Rv2[(int)Rv2_wire[i]] += 1.0;
	if(Rv2_wire[i+1]-Rv2_wire[i]!=1.0 && Rv2_wire[i+1]-Rv2_wire[i]<=3.0 && Rv2_nclust==0 && i<NRv2_wire-1)
	  for(int j=(int)Rv2_wire[i]+1;j<Rv2_wire[i+1];j++)Rv2_err[j] +=1.0;
      }    
    }
    //======= LU1 ======================//
    if(Lu1_nclust==0 && trig){
      for(int i=0;i<NLu1_wire;i++){
	Lu1[(int)Lu1_wire[i]] += 1.0;
	if(Lu1_wire[i+1]-Lu1_wire[i]!=1.0 && Lu1_wire[i+1]-Lu1_wire[i]<=3.0 && Lu1_nclust==0 && i<NLu1_wire-1)
	  for(int j=(int)Lu1_wire[i]+1;j<Lu1_wire[i+1];j++)Lu1_err[j] +=1.0;
      }
    }
    //======= LU2 ======================//
    if(Lu2_nclust==0 && trig){
      for(int i=0;i<NLu2_wire;i++){
	Lu2[(int)Lu2_wire[i]] += 1.0;
	if(Lu2_wire[i+1]-Lu2_wire[i]!=1.0 && Lu2_wire[i+1]-Lu2_wire[i]<=3.0 && Lu2_nclust==0  && i<NLu2_wire-1)
	  for(int j=(int)Lu2_wire[i]+1;j<Lu2_wire[i+1];j++)Lu2_err[j] +=1.0;
      }    
    }
    //======= Lv1 ======================//
    if(Lv1_nclust==0 && trig){
      for(int i=0;i<NLv1_wire;i++){
	Lv1[(int)Lv1_wire[i]] += 1.0;
	if(Lv1_wire[i+1]-Lv1_wire[i]!=1.0 && Lv1_wire[i+1]-Lv1_wire[i]<=3.0 && Lv1_nclust==0  && i<NLv1_wire-1 )
	  for(int j=(int)Lv1_wire[i]+1;j<Lv1_wire[i+1];j++)Lv1_err[j] +=1.0;
      }    
    }
    //======= Lv2 ======================//
    if(Lv2_nclust==0 && trig){
      for(int i=0;i<NLv2_wire;i++){
	Lv2[(int)Lv2_wire[i]] += 1.0;
	if(Lv2_wire[i+1]-Lv2_wire[i]!=1.0 && Lv2_wire[i+1]-Lv2_wire[i]<=3.0 && Lv2_nclust==0 && i<NLv2_wire-1){
	  for(int j=(int)Lv2_wire[i]+1;j<Lv2_wire[i+1];j++){Lv2_err[j] +=1.0;
	    //	    cout<<"Lv2_wire[i] "<<Lv2_wire[i]<<" Lv2_wire[i+1] "<<Lv2_wire[i+1]<<" j "<<j<<endl;	  
	    //	    for(int l=0;l<NLv2_wire;l++)cout<<" Lv2_wire "<<l<<" "<<Lv2_wire[l];
	    //	    cout<<endl;
	    //	    cout<<endl;
	  }
	}
      }
      
    }     
 
  }//End Fill

  for(int i=0;i<nwire;i++){
    hRu1_eff->Fill(i,Ru1[i]/(Ru1[i]+Ru1_err[i]));
    hRu2_eff->Fill(i,Ru2[i]/(Ru2[i]+Ru2_err[i]));
    hRv1_eff->Fill(i,Rv1[i]/(Rv1[i]+Rv1_err[i]));
    hRv2_eff->Fill(i,Rv2[i]/(Rv2[i]+Rv2_err[i]));
    hLu1_eff->Fill(i,Lu1[i]/(Lu1[i]+Lu1_err[i]));
    hLu2_eff->Fill(i,Lu2[i]/(Lu2[i]+Lu2_err[i]));
    hLv1_eff->Fill(i,Lv1[i]/(Lv1[i]+Lv1_err[i]));
    hLv2_eff->Fill(i,Lv2[i]/(Lv2[i]+Lv2_err[i]));

    if(Ru1[i]>0)gRu1_eff->SetPoint(i,i,Ru1[i]/(Ru1[i]+Ru1_err[i]));
    else gRu1_eff->SetPoint(i,i,0.0);
    if(Ru2[i]>0)gRu2_eff->SetPoint(i,i,Ru2[i]/(Ru2[i]+Ru2_err[i]));
    else gRu2_eff->SetPoint(i,i,0.0);
    if(Rv1[i]>0)gRv1_eff->SetPoint(i,i,Rv1[i]/(Rv1[i]+Rv1_err[i]));
    else gRv1_eff->SetPoint(i,i,0.0);
    if(Rv2[i]>0)gRv2_eff->SetPoint(i,i,Rv2[i]/(Rv2[i]+Rv2_err[i]));
    else gRv2_eff->SetPoint(i,i,0.0);


    if(Lu1[i]>0)gLu1_eff->SetPoint(i,i,Lu1[i]/(Lu1[i]+Lu1_err[i]));
    else gLu1_eff->SetPoint(i,i,0.0);
    if(Lu2[i]>0)gLu2_eff->SetPoint(i,i,Lu2[i]/(Lu2[i]+Lu2_err[i]));
    else gLu2_eff->SetPoint(i,i,0.0);
    if(Lv1[i]>0)gLv1_eff->SetPoint(i,i,Lv1[i]/(Lv1[i]+Lv1_err[i]));
    else gLv1_eff->SetPoint(i,i,0.0);
    if(Lv2[i]>0)gLv2_eff->SetPoint(i,i,Lv2[i]/(Lv2[i]+Lv2_err[i]));
    else gLv2_eff->SetPoint(i,i,0.0);

  }
  
  TCanvas* c0=new TCanvas("c0","c0");  
  c0->Divide(2,2);
  c0->cd(1);
  gRu1_eff->Draw("AP");
  c0->cd(2);
  gRv1_eff->Draw("AP");  
  c0->cd(3);
  gRu2_eff->Draw("AP");
  c0->cd(4);
  gRv2_eff->Draw("AP");

  
  TCanvas* c1=new TCanvas("c1","c1");  
  c1->Divide(2,2);
  c1->cd(1);
  gLu1_eff->Draw("AP");
  c1->cd(2);
  gLv1_eff->Draw("AP");  
  c1->cd(3);
  gLu2_eff->Draw("AP");
  c1->cd(4);
  gLv2_eff->Draw("AP");

  
}
